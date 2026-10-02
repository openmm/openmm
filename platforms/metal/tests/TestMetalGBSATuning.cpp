/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * Metal Platform tests: copyright (c) 2026 Chun-Chi Hung.                      *
 * Author: Chun-Chi Hung.                                                       *
 * This program is free software under the GNU Lesser General Public License, *
 * version 3 or (at your option) any later version, without any warranty.       *
 * See <http://www.gnu.org/licenses/> for the license.                          *
 * -------------------------------------------------------------------------- */

#include "MetalContext.h"
#include "MetalNonbondedUtilities.h"
#include "CommonKernelSources.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

using namespace OpenMM;
using namespace std;

/** @brief Inspect the actual generated source before normal Metal compilation. */
class RecordingContext : public MetalContext {
public:
    RecordingContext(const System& system, bool floating, bool reference) :
            MetalContext(system, nullptr, nullptr, floating), reference(reference) {}
    ComputeProgram compileProgram(const string source, const map<string, string>& defines={}) override {
        if (source.find("global_obc0_bornForce") != string::npos)
            nonbondedSource = source;
        map<string, string> selected = defines;
        if (reference) {
            selected["USE_GBSA_CHAIN_RULE_GUARD"] = "0";
        }
        return MetalContext::compileProgram(source, selected);
    }
    string nonbondedSource;
    bool reference;
};

/**
 * @brief Compare guarded pair arithmetic on tail tiles and both buffer ABIs.
 * The reference explicitly disables the shader-side guard.
 * Global fixed-point input is unchanged, including signed conversion boundaries.
 */
vector<double> evaluateChainRule(bool reference, bool floating, int atoms) {
    System system;
    for (int i = 0; i < atoms; i++) system.addParticle(1);
    RecordingContext context(system, floating, reference);
    MetalNonbondedUtilities& nb = static_cast<MetalNonbondedUtilities&>(context.getNonbondedUtilities());
    const int padded = context.getPaddedNumAtoms();
    ComputeArray params, born;
    params.initialize<mm_float2>(context, padded, "obc0_obcParams");
    born.initialize(context, padded, floating ? sizeof(float) : sizeof(int64_t), "obc0_bornForce");
    map<string, string> names = {{"OBC_PARAMS1", "obc0_obcParams1"}, {"OBC_PARAMS2", "obc0_obcParams2"},
        {"BORN_FORCE1", "obc0_bornForce1"}, {"BORN_FORCE2", "obc0_bornForce2"}};
    string source = context.replaceStrings(CommonKernelSources::gbsaObc2, names);
    nb.addInteraction(false, false, false, 1.0, vector<vector<int> >(), source, 0);
    nb.addParameter(ComputeParameterInfo(params, "obc0_obcParams", "float", 2));
    nb.addParameter(ComputeParameterInfo(born, "obc0_bornForce", "mm_long", 1));
    context.initialize();
    vector<mm_float4> positions(padded, mm_float4(0, 0, 0, 0));
    vector<mm_float2> radii(padded, mm_float2(0.1f, 0.08f));
    vector<int64_t> raw(padded, 0);
    vector<float> floats(padded, 0);
    const int64_t cases[] = {0, 1, -1, (int64_t(1)<<24)-1, (int64_t(1)<<24)+1,
        -(int64_t(1)<<24)-1, (int64_t(1)<<32)+257, -(int64_t(1)<<32)-257,
        (int64_t(1)<<45)+12345, -(int64_t(1)<<45)-12345};
    for (int i = 0; i < atoms; i++) {
        positions[i] = mm_float4(0.4f*(i%7), 0.43f*((i/7)%5), 0.47f*(i/35), 0);
        radii[i] = mm_float2(0.1f+0.003f*(i%5), 0.08f+0.002f*(i%7));
        raw[i] = cases[i%(sizeof(cases)/sizeof(cases[0]))];
        floats[i] = float(raw[i])*(1.0f/4294967296.0f);
    }
    context.getPosq().upload(positions);
    params.upload(radii);
    if (floating) born.upload(floats);
    else born.upload(raw);
    context.clearBuffer(context.getLongForceBuffer());
    context.clearBuffer(context.getEnergyBuffer());
    nb.prepareInteractions(1);
    nb.computeInteractions(1, true, false);
    const string& generated = context.nonbondedSource;
    ASSERT(!generated.empty());
    ASSERT(generated.find("mm_long* restrict global_obc0_bornForce") != string::npos);
    vector<double> forces;
    context.downloadFixedPointBuffer(context.getLongForceBuffer(), forces);
    for (double value : forces) ASSERT(isfinite(value));
    // Kernels only read this input; fixed-point storage and all low bits survive.
    if (!floating) {
        vector<int64_t> after;
        born.download(after);
        ASSERT_EQUAL(raw.size(), after.size());
        for (int i = 0; i < padded; i++) ASSERT_EQUAL(raw[i], after[i]);
    }
    return forces;
}

/**
 * @brief Compare formulas without requiring deterministic float atomic order.
 * The large signed inputs deliberately create cancellation. Repeated runs of
 * the unchanged float reference differ across tiles, but the fixed-point and
 * single-tile cases do not need this allowance. The RMS-scaled floor is a
 * conservative fixture allowance, not a bound on arbitrary summation error.
 */
void compareForces(const vector<double>& reference, const vector<double>& actual, bool floating, int atoms) {
    ASSERT_EQUAL(reference.size(), actual.size());
    double referenceSquared = 0, errorSquared = 0;
    for (int i = 0; i < actual.size(); i++) {
        referenceSquared += reference[i]*reference[i];
        const double error = actual[i]-reference[i];
        errorSquared += error*error;
    }
    const double referenceRms = sqrt(referenceSquared/reference.size());
    const int atomBlocks = (atoms+31)/32;
    const double roundoff = floating && atomBlocks > 1 ?
            2*atomBlocks*numeric_limits<float>::epsilon()*referenceRms : 0;
    ASSERT_EQUAL_TOL(0.0, sqrt(errorSquared/reference.size()), 2e-5*max(1.0, referenceRms));
    for (int i = 0; i < actual.size(); i++)
        ASSERT_EQUAL_TOL(0.0, actual[i]-reference[i], 2e-5*max(1.0, abs(reference[i]))+roundoff);
}

int main() {
    try {
        for (bool floating : {false, true}) {
            for (int atoms : {31, 97}) {
                vector<double> reference = evaluateChainRule(true, floating, atoms);
                vector<double> actual = evaluateChainRule(false, floating, atoms);
                compareForces(reference, actual, floating, atoms);
            }
        }
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
    cout << "Metal GBSA chain-rule guard tests passed" << endl;
    return 0;
}
