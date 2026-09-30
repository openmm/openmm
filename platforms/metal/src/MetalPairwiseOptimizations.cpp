/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Adapted from Common pairwise kernels used by CUDA, HIP, and OpenCL.        *
 * Original OpenMM code:                                                      *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.       *
 * Authors: Peter Eastman                                                     *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 * This program is distributed WITHOUT ANY WARRANTY; without even the        *
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. *
 * See the GNU Lesser General Public License for more details.                *
 * You should have received a copy of the GNU Lesser General Public License  *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#include "MetalPairwiseOptimizations.h"
#include "CommonKernelSources.h"
#include "openmm/OpenMMException.h"
#include <cctype>
#include <regex>
#include <set>
#include <vector>

using namespace OpenMM;
using namespace std;

namespace {

/** Match all literal template text, allowing only Common's existing holes. */
bool matchesTemplate(const string& source, const string& pattern, const set<string>& holes = {}) {
    vector<string> literals;
    size_t literalStart = 0;
    for (size_t i = 0; i < pattern.size();) {
        if (!isalpha(static_cast<unsigned char>(pattern[i])) && pattern[i] != '_') {
            i++;
            continue;
        }
        size_t end = i+1;
        while (end < pattern.size() && (isalnum(static_cast<unsigned char>(pattern[end])) || pattern[end] == '_')) end++;
        if (holes.count(pattern.substr(i, end-i))) {
            literals.push_back(pattern.substr(literalStart, i-literalStart));
            literalStart = end;
        }
        i = end;
    }
    literals.push_back(pattern.substr(literalStart));
    if (literals.size() == 1)
        return source == pattern;
    size_t cursor = 0;
    for (size_t i = 0; i < literals.size(); i++) {
        size_t next = source.find(literals[i], cursor);
        if (next == string::npos || (i == 0 && next != 0)) return false;
        cursor = next+literals[i].size();
    }
    return cursor == source.size();
}

/** Replace a reviewed literal and optionally check its template occurrence count. */
void replace(string& source, const string& oldText, const string& newText, int expected = -1) {
    size_t pos = 0;
    int count = 0;
    while ((pos = source.find(oldText, pos)) != string::npos) {
        source.replace(pos, oldText.size(), newText);
        pos += newText.size();
        count++;
    }
    if (expected >= 0 && count != expected)
        throw OpenMMException("Common template changed: review Metal pairwise transformation of "+oldText);
}

/** Read-only lane data or a lane-owned accumulator that travels with its particle. */
struct Register {
    string type, name;
    bool readOnly;
};

string shuffleAssignment(const string& destination, const string& value, const string& lane) {
    return destination+" = simdShuffle("+value+", "+lane+");\n";
}

/** MetalKernel separately verifies the hardware SIMD width and padded launch. */
string checkedSource(const string& source) {
    return "#if defined(TILE_SIZE) && TILE_SIZE != 32\n"
            "#error Metal pairwise register paths require TILE_SIZE=32\n#endif\n"+source;
}

/** Replace a known array's audited subscripts after diagonal broadcasts are made. */
void scalarize(string& source, const string& name) {
    for (const string& index : {"LOCAL_ID", "localAtomIndex", "atom2", "tbx+tj"})
        replace(source, name+"["+index+"]", name);
}

/** Keep diagonal broadcast shuffles before any cutoff or particle validity branch. */
string customGB(string source) {
    const regex declaration("LOCAL ([A-Za-z0-9_]+) (local_[A-Za-z0-9_]+)\\[LOCAL_BUFFER_SIZE\\];");
    vector<Register> registers;
    for (sregex_iterator i(source.begin(), source.end(), declaration), end; i != end; ++i) {
        const string name = (*i)[2];
        const bool readOnly = name == "local_pos" || name.find("local_params") == 0 || name.find("local_values") == 0;
        const bool accumulator = name == "local_value" || name == "local_force" ||
                name.find("local_deriv") == 0 || name.find("local_dValue0dParam") == 0;
        if (!readOnly && !accumulator)
            throw OpenMMException("Unrecognized Common CustomGB lane field: "+name);
        registers.push_back({(*i)[1], name, readOnly});
    }
    if (registers.empty()) throw OpenMMException("Missing Common CustomGB lane data");
    source = regex_replace(source, declaration, "$1 $2 = $1(0);");
    const string diagonalEnd = "\n        else {\n            // This is an off-diagonal tile.";
    const size_t end = source.find(diagonalEnd);
    if (end == string::npos) throw OpenMMException("Missing Common CustomGB diagonal tile");
    string diagonal = source.substr(0, end);
    string rest = source.substr(end);
    string broadcast;
    for (const Register& reg : registers) {
        if (reg.readOnly) {
            const string value = "_metal_broadcast_"+reg.name;
            broadcast += reg.type+" "+shuffleAssignment(value, reg.name, "j");
            replace(diagonal, reg.name+"[atom2]", value);
        }
    }
    const string diagonalLoop = "for (unsigned int j = 0; j < TILE_SIZE; j++) {";
    replace(diagonal, diagonalLoop, diagonalLoop+"\n"+broadcast, 1);
    // The diagonal no longer communicates through local memory.
    replace(diagonal, "SYNC_WARPS;", "");
    source = diagonal+rest;
    for (const Register& reg : registers) scalarize(source, reg.name);
    // Common emits this token-pasting macro before the energy kernel.
    replace(source, "local_deriv##INDEX[LOCAL_ID]", "local_deriv##INDEX");
    replace(source, "LOCAL int atomIndices[LOCAL_BUFFER_SIZE];", "", 1);
    replace(source, "const unsigned int tbx = LOCAL_ID - tgx;", "const unsigned int tbx = LOCAL_ID - tgx;\nint atomIndices = 0;", 1);
    scalarize(source, "atomIndices");
    string rotate;
    for (const Register& reg : registers)
        rotate += shuffleAssignment(reg.name, reg.name, "(tgx+1)&31");
    rotate += shuffleAssignment("atomIndices", "atomIndices", "(tgx+1)&31");
    const regex next("tj = \\(tj \\+ 1\\) & \\(TILE_SIZE - 1\\);\\s*SYNC_WARPS;");
    if (distance(sregex_iterator(source.begin(), source.end(), next), sregex_iterator()) != 3)
        throw OpenMMException("Missing Common CustomGB ring steps");
    source = regex_replace(source, next, "tj = (tj + 1) & (TILE_SIZE - 1);\n"+rotate);
    return source;
}

/** DPD preserves i=0..31 visitation, hence the original per-lane RNG consumption. */
string dpdParticles(string source) {
    replace(source, "LOCAL mixed3 localPos[WORK_GROUP_SIZE];", "mixed3 localPos = make_mixed3(0);", 1);
    replace(source, "LOCAL mixed4 localVel[WORK_GROUP_SIZE];", "mixed4 localVel = make_mixed4(0);", 1);
    replace(source, "LOCAL volatile int localType[WORK_GROUP_SIZE];", "int localType = 0;", 1);
    replace(source, "LOCAL int atomIndices[WORK_GROUP_SIZE];", "int atomIndices = 0;", 1);
    for (const string& name : {"localPos", "localVel", "localType", "atomIndices"}) {
        replace(source, name+"[LOCAL_ID]", name);
        replace(source, name+"[tbx+i]", "_metal_broadcast_"+name);
    }
    // Padded atom1 lanes still provide their particle to other active lanes.
    replace(source, "if (atom1 < numAtoms) {", "{", 3);
    replace(source, "if ((x != y || atom1 < atom2) && atom2 < numAtoms)",
            "if (atom1 < numAtoms && (x != y || atom1 < atom2) && atom2 < numAtoms)", 1);
    replace(source, "if (atom2 < numAtoms)", "if (atom1 < numAtoms && atom2 < numAtoms)", 2);
    const string loop = "for (int i = 0; i < TILE_SIZE; i++) {";
    const string broadcast = "\nmixed3 _metal_broadcast_localPos = simdShuffle(localPos, i);\n"
            "mixed4 _metal_broadcast_localVel = simdShuffle(localVel, i);\n"
            "int _metal_broadcast_localType = simdShuffle(localType, i);\n";
    const size_t split = source.find("// Second loop: process tiles from the neighbor list.");
    string first = source.substr(0, split), second = source.substr(split);
    replace(first, loop, loop+broadcast, 1);
    replace(second, loop, loop+broadcast+"int _metal_broadcast_atomIndices = simdShuffle(atomIndices, i);\n", 2);
    return first+second;
}

const set<string> customGBHoles = {"PARAMETER_ARGUMENTS", "ATOM_PARAMETER_DATA", "LOAD_ATOM1_PARAMETERS",
    "LOAD_ATOM2_PARAMETERS", "LOAD_LOCAL_PARAMETERS_FROM_1", "LOAD_LOCAL_PARAMETERS_FROM_GLOBAL",
    "COMPUTE_VALUE", "ADD_TEMP_DERIVS1", "ADD_TEMP_DERIVS2", "STORE_PARAM_DERIVS1", "STORE_PARAM_DERIVS2",
    "INIT_PARAM_DERIVS", "DECLARE_ATOM1_DERIVATIVES", "CLEAR_LOCAL_DERIVATIVES", "COMPUTE_INTERACTION",
    "RECORD_DERIVATIVE_2", "STORE_DERIVATIVES_1", "STORE_DERIVATIVES_2", "SAVE_PARAM_DERIVS"};

} // namespace

MetalPairwiseOptimizations::Settings MetalPairwiseOptimizations::getBuildSettings() {
    Settings settings;
#if OPENMM_METAL_FAST_CUSTOM_GB_VALUE_SHUFFLE
    settings.customGBValue = true;
#endif
#if OPENMM_METAL_FAST_CUSTOM_GB_ENERGY_SHUFFLE
    settings.customGBEnergy = true;
#endif
#if OPENMM_METAL_FAST_DPD_PARTICLE_SHUFFLE
    settings.dpdParticles = true;
#endif
    return settings;
}

string MetalPairwiseOptimizations::apply(const string& source) {
    return apply(source, getBuildSettings());
}

string MetalPairwiseOptimizations::apply(const string& source, const Settings& settings) {
    if (settings.customGBValue && matchesTemplate(source, CommonKernelSources::customGBValueN2, customGBHoles))
        return checkedSource(customGB(source));
    if (settings.customGBEnergy && matchesTemplate(source, CommonKernelSources::customGBEnergyN2, customGBHoles))
        return checkedSource(customGB(source));
    if (settings.dpdParticles && source == CommonKernelSources::dpd)
        return checkedSource(dpdParticles(source));
    return source;
}
