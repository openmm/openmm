/* -------------------------------------------------------------------------- *
 * OpenMM — Metal Platform                                                    *
 * Portions copyright (c) 2026 Chun-Chi Hung. Authors: Chun-Chi Hung.           *
 * This is part of the OpenMM molecular simulation toolkit.                    *
 * See https://openmm.org/development.                                         *
 * This program is free software: you can redistribute it and/or modify        *
 * it under the terms of the GNU Lesser General Public License as published    *
 * by the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version.                                         *
 * This program is distributed WITHOUT ANY WARRANTY; without even the          *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   *
 * See the GNU Lesser General Public License for more details.                  *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.         *
 * -------------------------------------------------------------------------- */

#include "MetalPlatform.h"
#include "MetalBuildConfiguration.h"
#include "openmm/CompoundIntegrator.h"
#include "openmm/Context.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include "SHA1.h"
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <sys/utsname.h>

using namespace OpenMM;
using namespace std;

namespace {
/** @brief Escape JSON metadata, including externally supplied filenames. */
string quoted(const string& value) {
    ostringstream output;
    output << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') output << '\\' << c;
        else if (c < 32) output << "\\u" << hex << setw(4) << setfill('0') << int(c);
        else output << c;
    }
    return output.str()+'"';
}

/** @brief Read a caller-supplied OpenMM XML object without requiring Python bindings. */
template <class T> unique_ptr<T> readXML(const string& path) {
    ifstream input(path.c_str());
    if (!input) throw OpenMMException("Cannot open benchmark input: "+path);
    return unique_ptr<T>(XmlSerializer::deserialize<T>(input));
}

/** @brief Content identity for reproducibility, not a cryptographic authenticity check. */
template <class T> string fingerprint(const T& object) {
    ostringstream xml;
    XmlSerializer::serialize(&object, "BenchmarkInput", xml);
    const string text = xml.str();
    CSHA1 hash;
    for (size_t offset = 0; offset < text.size(); offset += 1048576)
        hash.Update(reinterpret_cast<const unsigned char*>(text.data()+offset),
                static_cast<unsigned int>(min(size_t(1048576), text.size()-offset)));
    hash.Final();
    unsigned char digest[20];
    hash.GetHash(digest);
    ostringstream result;
    for (unsigned char byte : digest) result << hex << setw(2) << setfill('0') << int(byte);
    return result.str();
}

/** @brief Reject auto-generated stochastic seeds, including nested Force/CompoundIntegrator nodes. */
void requireExplicitSeeds(const SerializationNode& node) {
    if (node.hasProperty("randomSeed") && node.getIntProperty("randomSeed") == 0)
        throw OpenMMException("Benchmark input uses an automatic randomSeed=0. Set explicit nonzero seeds before XML export.");
    for (const SerializationNode& child : node.getChildren()) requireExplicitSeeds(child);
}
template <class T> void requireExplicitSeeds(const T& object) {
    SerializationNode node;
    SerializationProxy::getProxy(typeid(object)).serialize(&object, node);
    requireExplicitSeeds(node);
}

/** @brief Match reserved random variables, not longer user-defined identifiers. */
bool usesStochasticVariables(const string& expression) {
    for (size_t begin = 0; begin < expression.size();) {
        const unsigned char first = expression[begin];
        if (!isalpha(first) && first != '_') {
            ++begin;
            continue;
        }
        size_t end = begin+1;
        while (end < expression.size() &&
                (isalnum(static_cast<unsigned char>(expression[end])) || expression[end] == '_'))
            ++end;
        const string variable = expression.substr(begin, end-begin);
        if (variable == "uniform" || variable == "gaussian") return true;
        begin = end;
    }
    return false;
}

/**
 * @brief Limit benchmark repeats to validated CustomIntegrator checkpoint state.
 * Common's separate per-DOF uniform RNG is not included in its checkpoint.  This
 * harness conservatively excludes random expressions; the backend still supports
 * them. Check every CompoundIntegrator child, including currently inactive ones.
 */
void requireCheckpointSafeIntegrator(const Integrator& integrator) {
    const CompoundIntegrator* compound = dynamic_cast<const CompoundIntegrator*>(&integrator);
    if (compound != nullptr) {
        for (int i = 0; i < compound->getNumIntegrators(); ++i)
            requireCheckpointSafeIntegrator(compound->getIntegrator(i));
    }
    const CustomIntegrator* custom = dynamic_cast<const CustomIntegrator*>(&integrator);
    if (custom == nullptr) return;
    bool stochastic = usesStochasticVariables(custom->getKineticEnergyExpression());
    for (int i = 0; i < custom->getNumComputations() && !stochastic; ++i) {
        CustomIntegrator::ComputationType type;
        string variable, expression;
        custom->getComputationStep(i, type, variable, expression);
        stochastic = usesStochasticVariables(expression);
    }
    if (stochastic)
        throw OpenMMException("Benchmark checkpoint repeats do not support stochastic CustomIntegrator "
                "expressions (uniform/gaussian), including kinetic-energy expressions. "
                "Use a deterministic CustomIntegrator or a supported seeded integrator; "
                "this is a benchmark limitation, not a Metal backend limitation.");
}

/** @brief Create a neutral, deterministic synthetic workload, not a molecular benchmark. */
void createSynthetic(const string& name, int count, System& system, vector<Vec3>& positions) {
    if (name != "cutoff" && name != "pme" && name != "ljpme" && name != "gbsa")
        throw OpenMMException("Synthetic case must be cutoff, pme, ljpme, or gbsa");
    const int side = int(ceil(cbrt(double(count))));
    const double spacing = 0.4, edge = max(3.0, side*spacing);
    system.setDefaultPeriodicBoxVectors(Vec3(edge, 0, 0), Vec3(0, edge, 0), Vec3(0, 0, edge));
    NonbondedForce* nonbonded = new NonbondedForce();
    system.addForce(nonbonded);
    nonbonded->setNonbondedMethod(name == "pme" ? NonbondedForce::PME :
        name == "ljpme" ? NonbondedForce::LJPME : name == "gbsa" ?
        NonbondedForce::CutoffNonPeriodic : NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    nonbonded->setEwaldErrorTolerance(1e-4);
    GBSAOBCForce* gbsa = nullptr;
    if (name == "gbsa") {
        gbsa = new GBSAOBCForce();
        gbsa->setNonbondedMethod(GBSAOBCForce::CutoffNonPeriodic);
        gbsa->setCutoffDistance(1.0);
        system.addForce(gbsa);
    }
    for (int i = 0; i < count; i++) {
        system.addParticle(39.9);
        double charge = count%2 && i == count-1 ? 0 : (i%2 ? -0.1 : 0.1);
        nonbonded->addParticle(charge, 0.3, 0.1);
        if (gbsa != nullptr) gbsa->addParticle(charge, 0.15, 0.8);
        positions.push_back(Vec3(i%side, (i/side)%side, i/(side*side))*spacing);
    }
}
} // namespace

/**
 * @brief Measure synchronized wall time after compilation/warmup, restoring the
 * pre-warmup checkpoint for every repeat. No kernel-only speedup is claimed.
 * Serialized System/State/Integrator inputs support the existing OpenMM benchmark
 * systems and additional core force workloads without cloning their setup code.
 */
int main(int argc, char* argv[]) {
    try {
        map<string, string> options = {{"case", "pme"}, {"particles", "4096"},
            {"steps", "100"}, {"warmup", "20"}, {"repeats", "5"}};
        for (int i = 1; i < argc; i += 2) {
            const string key = argv[i];
            if (key == "--help") {
                cout << "BenchmarkMetal [--case cutoff|pme|ljpme|gbsa] [--particles N]\n"
                     << "  [--system system.xml --state state.xml [--integrator integrator.xml]]\n"
                     << "  [--steps N --warmup N --repeats N]\n"
                     << "Emits JSON lines; compare identical inputs with one build switch changed.\n";
                return 0;
            }
            if (i+1 >= argc || key.size() < 3 || key.substr(0, 2) != "--")
                throw OpenMMException("Expected --option value; use --help");
            const string name = key.substr(2);
            if (options.count(name) == 0 && name != "system" && name != "state" && name != "integrator")
                throw OpenMMException("Unknown option: "+key);
            options[name] = argv[i+1];
        }
        const int steps = stoi(options["steps"]), warmup = stoi(options["warmup"]), repeats = stoi(options["repeats"]);
        const int count = stoi(options["particles"]);
        if (steps <= 0 || warmup < 1 || repeats < 1 || count < 1)
            throw OpenMMException("Steps, warmup, repeats, and particle count must be positive");
        unique_ptr<System> system;
        unique_ptr<State> state;
        unique_ptr<Integrator> integrator;
        vector<Vec3> positions;
        if (options.count("system") || options.count("state")) {
            if (!options.count("system") || !options.count("state"))
                throw OpenMMException("Provide both --system and --state XML inputs");
            system = readXML<System>(options["system"]);
            state = readXML<State>(options["state"]);
        }
        else {
            system.reset(new System());
            createSynthetic(options["case"], count, *system, positions);
        }
        if (options.count("integrator")) integrator = readXML<Integrator>(options["integrator"]);
        else integrator.reset(new VerletIntegrator(0.001));
        requireExplicitSeeds(*system);
        requireExplicitSeeds(*integrator);
        requireCheckpointSafeIntegrator(*integrator);
        const string systemHash = fingerprint(*system);
        const Integrator& integratorObject = *integrator;
        const string integratorType = SerializationProxy::getProxy(typeid(integratorObject)).getTypeName();
        MetalPlatform platform;
        Context context(*system, *integrator, platform);
        if (state) context.setState(*state);
        else context.setPositions(positions);
        // State restoration may change CustomIntegrator parameters or the active
        // CompoundIntegrator, so identify the actual initial state after setState().
        const string integratorHash = fingerprint(*integrator);
        const double initialStepSize = integrator->getStepSize();
        const string stateHash = fingerprint(context.getState(State::Positions|State::Velocities|
                State::Parameters|State::IntegratorParameters));
        context.getState(State::Energy);
        stringstream checkpoint(ios_base::in | ios_base::out | ios_base::binary);
        context.createCheckpoint(checkpoint);
        integrator->step(warmup);
        context.getState(State::Energy); // Finish compilation, then discard warmup's state in every repeat.
        struct utsname os;
        if (uname(&os) != 0) throw OpenMMException("Cannot read OS version");
        cout << setprecision(12) << "{\"type\":\"configuration\",\"case\":"
             << quoted(state ? options["system"] : options["case"])
             << ",\"device\":" << quoted(platform.getPropertyValue(context, "DeviceName"))
             << ",\"os\":" << quoted(string(os.sysname)+" "+os.release)
             << ",\"particles\":" << system->getNumParticles()
             << ",\"system_file\":" << quoted(options.count("system") ? options["system"] : "")
             << ",\"state_file\":" << quoted(options.count("state") ? options["state"] : "")
             << ",\"integrator_file\":" << quoted(options.count("integrator") ? options["integrator"] : "")
             << ",\"system_xml_sha1\":" << quoted(systemHash)
             << ",\"initial_state_xml_sha1\":" << quoted(stateHash)
             << ",\"integrator_xml_sha1\":" << quoted(integratorHash)
             << ",\"integrator_type\":" << quoted(integratorType)
             << ",\"initial_step_size_ps\":" << initialStepSize
             << ",\"warmup_steps\":" << warmup << ",\"flags\":" << METAL_BUILD_CONFIGURATION << "}" << endl;
        for (int repeat = 0; repeat < repeats; repeat++) {
            checkpoint.clear();
            checkpoint.seekg(0);
            context.loadCheckpoint(checkpoint);
            context.getState(State::Energy); // Restore and invalidate/rebuild outside timing.
            const auto start = chrono::steady_clock::now();
            integrator->step(steps);
            const State finalState = context.getState(State::Energy);
            const double elapsed = chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if (!isfinite(finalState.getPotentialEnergy()) || !isfinite(finalState.getKineticEnergy()))
                throw OpenMMException("Non-finite benchmark energy; timing is invalid");
            cout << "{\"type\":\"measurement\",\"repeat\":" << repeat
                 << ",\"steps\":" << steps << ",\"wall_seconds\":" << elapsed
                 << ",\"steps_per_second\":" << steps/elapsed
                 << ",\"potential_energy\":" << finalState.getPotentialEnergy()
                 << ",\"kinetic_energy\":" << finalState.getKineticEnergy() << "}" << endl;
        }
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}
