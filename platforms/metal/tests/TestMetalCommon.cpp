/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from platforms/opencl/tests/OpenCLTests.h and its test wrappers.
 * Original OpenCL Platform test code:
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.
 * Authors: Peter Eastman, Evan Pretti
 *
 * Metal Platform tests:
 * Portions copyright (c) 2026 Chun-Chi Hung.
 * Authors: Chun-Chi Hung
 *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "MetalPlatform.h"
#include "openmm/OpenMMException.h"
#include <string>

#ifndef OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
#define OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS 1
#endif

#ifdef METAL_TEST_FLOAT_ACCUMULATORS
OpenMM::MetalPlatform platform(true);
#else
OpenMM::MetalPlatform platform;
#endif

/**
 * @brief Initialize an existing Common test for the scoped Metal baseline.
 * @param argc Number of arguments; optional arguments are single and device index 0.
 * @param argv Program name followed by optional precision and device index.
 * @throws OpenMMException If a test requests an unsupported execution mode.
 *
 * The Common test assertions and tolerances are reused unchanged. Passing them
 * is provisional regression evidence, not a new numerical-accuracy contract.
 */
void initializeTests(int argc, char* argv[]) {
    if (argc > 3)
        throw OpenMM::OpenMMException("Metal tests accept only optional precision and device-index arguments");
    if (argc > 1 && std::string(argv[1]) != "single")
        throw OpenMM::OpenMMException("Metal tests support only single precision");
    if (argc > 2 && std::string(argv[2]) != "0")
        throw OpenMM::OpenMMException("Metal tests support only the default GPU (DeviceIndex 0)");
    platform.setPropertyDefaultValue("Precision", "single");
    platform.setPropertyDefaultValue("DeviceIndex", "0");
    platform.setPropertyDefaultValue("UseCpuPme", "false");
}

#ifndef METAL_COMMON_TEST_HEADER
#error "Define METAL_COMMON_TEST_HEADER to an existing Common test header"
#endif
#if defined(METAL_TEST_CHECKPOINTS)
// The generic checkpoint main unconditionally requests two devices. Keep its
// single-device test functions, but provide a runner matching this platform.
#define main runUnscopedCheckpointSuite
#elif defined(METAL_TEST_MINIMIZER) && !OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
// The fixed-only configuration cannot use Common's CPU extreme-force fallback.
// Keep all ordinary assertions and test explicit rejection of that case below.
#define main runFloatingAccumulatorMinimizerSuite
#endif
#include METAL_COMMON_TEST_HEADER
#if defined(METAL_TEST_MINIMIZER) && !OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
#undef main
#endif
#if defined(METAL_TEST_CHECKPOINTS)
#undef main

/** @brief Preserve checkpoint and cross-Context continuation coverage on one GPU. */
void testMetalSingleDeviceCheckpoints() {
    System system;
    NonbondedForce* force = new NonbondedForce();
    system.addForce(force);
    system.addParticle(1.0);
    system.addParticle(1.0);
    force->addParticle(0.1, 0.2, 0.1);
    force->addParticle(-0.1, 0.2, 0.1);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(0, 0, 0), Vec3(1, 0, 0)});
    context.setVelocities({Vec3(0.1, 0, 0), Vec3(-0.1, 0, 0)});
    integrator.step(5);
    const int stateTypes = State::Positions | State::Velocities | State::Parameters;
    State saved = context.getState(stateTypes);
    stringstream checkpoint(ios_base::out | ios_base::in | ios_base::binary);
    context.createCheckpoint(checkpoint);
    integrator.step(3);
    State expected = context.getState(stateTypes);

    context.loadCheckpoint(checkpoint);
    State restored = context.getState(stateTypes);
    compareStates(saved, restored);
    integrator.step(3);
    State continued = context.getState(stateTypes);
    compareStates(expected, continued);

    VerletIntegrator secondIntegrator(0.001);
    Context second(system, secondIntegrator, platform);
    checkpoint.seekg(0, checkpoint.beg);
    second.loadCheckpoint(checkpoint);
    State transferred = second.getState(stateTypes);
    compareStates(saved, transferred);
    secondIntegrator.step(3);
    State secondContinued = second.getState(stateTypes);
    compareStates(expected, secondContinued);
}

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testSetState();
        testMetalSingleDeviceCheckpoints();
        testLangevin();
        runPlatformTests();
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
#endif

#if defined(METAL_TEST_NONBONDED)
#include "openmm/CustomNonbondedForce.h"

/**
 * @brief Exercise neighbor-list overflow recovery and the full energy/derivative buffers.
 *
 * All 4096 particles fit inside one cutoff sphere. The initial 20 tiles per
 * atom block cannot hold the dense pair list, so a valid result requires the
 * GPU tile counter, array growth, argument rebinding, and repeated evaluation.
 */
void testMetalDenseNeighborList() {
    const int count = 4096;
    System system;
    CustomNonbondedForce* force = new CustomNonbondedForce("a+b");
    force->addGlobalParameter("a", 1.0);
    force->addGlobalParameter("b", 0.0);
    force->addEnergyParameterDerivative("a");
    force->addEnergyParameterDerivative("b");
    force->setNonbondedMethod(CustomNonbondedForce::CutoffNonPeriodic);
    force->setCutoffDistance(1.0);
    vector<Vec3> positions(count);
    for (int i = 0; i < count; i++) {
        system.addParticle(1.0);
        force->addParticle(vector<double>());
        positions[i] = Vec3(0.005*(i%16), 0.005*((i/16)%16), 0.005*(i/256));
    }
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    double pairs = 0.5*count*(count-1);
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            context.setParameter("a", 2.0);
            context.setParameter("b", 1.0);
        }
        State state = context.getState(State::Energy | State::Forces | State::ParameterDerivatives);
        ASSERT_EQUAL_TOL((pass == 0 ? 1.0 : 3.0)*pairs, state.getPotentialEnergy(), 1e-5);
        ASSERT_EQUAL_TOL(pairs, state.getEnergyParameterDerivatives().at("a"), 1e-5);
        ASSERT_EQUAL_TOL(pairs, state.getEnergyParameterDerivatives().at("b"), 1e-5);
        const Vec3 zero(0, 0, 0);
        for (const Vec3& value : state.getForces()) {
            ASSERT_EQUAL_VEC(zero, value, 1e-6);
        }
    }
}
#endif

#if defined(METAL_TEST_MINIMIZER)
#include "openmm/ATMForce.h"
#include "openmm/CustomCVForce.h"
#include "openmm/CustomIntegrator.h"

/** @brief Recover large bonded forces without substituting a special nonbonded solver. */
void testMetalLargeBondMinimization() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.0, 1e22);
    system.addForce(bonds);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(0, 0, 0), Vec3(2, 0, 0)});
    LocalEnergyMinimizer::minimize(context, 1.0, 1000);
    vector<Vec3> positions = context.getState(State::Positions).getPositions();
    Vec3 delta = positions[1]-positions[0];
    ASSERT_EQUAL_TOL(1.0, sqrt(delta.dot(delta)), 1e-5);
}

/** @brief Preserve current parameters, selected groups, and non-position Context state. */
void testMetalLargeCustomMinimizationState() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x-target)^2");
    force->addGlobalParameter("k", 1e22);
    force->addGlobalParameter("target", -2.0);
    force->addParticle(0);
    force->setForceGroup(2);
    system.addForce(force);
    CustomExternalForce* excluded = new CustomExternalForce("1e22*(x+4)^2");
    excluded->addParticle(0);
    excluded->setForceGroup(1);
    system.addForce(excluded);
    system.setDefaultPeriodicBoxVectors(Vec3(2, 0, 0), Vec3(0.1, 3, 0), Vec3(0.2, 0.1, 4));
    VerletIntegrator integrator(0.001);
    integrator.setIntegrationForceGroups(1<<2);
    integrator.setConstraintTolerance(2e-5);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(2, 0.125, -0.25)});
    context.setVelocities({Vec3(0.25, 0.5, 0.75)});
    context.setParameter("target", 0.25);
    context.setPeriodicBoxVectors(Vec3(3, 0, 0), Vec3(0.1, 4, 0), Vec3(0.2, 0.3, 5));
    context.setTime(3.25);
    context.setStepCount(17);
    State before = context.getState(State::Positions | State::Velocities | State::Parameters);
    LocalEnergyMinimizer::minimize(context, 1.0, 1000);
    State after = context.getState(State::Positions | State::Velocities | State::Parameters);
    ASSERT_EQUAL_TOL(0.25, after.getPositions()[0][0], 1e-5);
    ASSERT_EQUAL(before.getPositions()[0][1], after.getPositions()[0][1]);
    ASSERT_EQUAL(before.getPositions()[0][2], after.getPositions()[0][2]);
    ASSERT_EQUAL_VEC(before.getVelocities()[0], after.getVelocities()[0], 0.0);
    ASSERT_EQUAL(before.getTime(), after.getTime());
    ASSERT_EQUAL(17, context.getStepCount());
    ASSERT_EQUAL(0.25, context.getParameter("target"));
    ASSERT_EQUAL(1<<2, integrator.getIntegrationForceGroups());
    ASSERT_EQUAL(2e-5, integrator.getConstraintTolerance());
    Vec3 beforeA, beforeB, beforeC, afterA, afterB, afterC;
    before.getPeriodicBoxVectors(beforeA, beforeB, beforeC);
    after.getPeriodicBoxVectors(afterA, afterB, afterC);
    ASSERT_EQUAL_VEC(beforeA, afterA, 0.0);
    ASSERT_EQUAL_VEC(beforeB, afterB, 0.0);
    ASSERT_EQUAL_VEC(beforeC, afterC, 0.0);
}

/** @brief Keep Common virtual-site distribution and constraint restraints in the recovery path. */
void testMetalLargeVirtualSiteMinimization() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.addConstraint(0, 1, 1.0);
    system.setVirtualSite(2, new TwoParticleAverageSite(0, 1, 0.5, 0.5));
    CustomExternalForce* force = new CustomExternalForce("0.5*k*x*x");
    force->addGlobalParameter("k", 1e20);
    force->addParticle(2);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    integrator.setConstraintTolerance(1e-5);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(1.5, 0, 0), Vec3(2.5, 0, 0), Vec3(2, 0, 0)});
    LocalEnergyMinimizer::minimize(context, 1.0, 1000);
    vector<Vec3> positions = context.getState(State::Positions).getPositions();
    Vec3 delta = positions[1]-positions[0];
    ASSERT_EQUAL_TOL(1.0, sqrt(delta.dot(delta)), 2e-5);
    const Vec3 midpoint = 0.5*(positions[0]+positions[1]);
    ASSERT_EQUAL_VEC(midpoint, positions[2], 1e-6);
    ASSERT_EQUAL_TOL(0.0, positions[2][0], 1e-5);
}

/** @brief Expose each reported position through the caller's Context and retain early stops. */
void testMetalMinimizerReporterContext() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("(x-1)^4");
    force->addParticle(0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(5, 0, 0)});
    class Reporter : public MinimizationReporter {
    public:
        Reporter(Context& context) : context(context), calls(0) {
        }
        bool report(int iteration, const vector<double>& x, const vector<double>& gradient,
                map<string, double>& statistics) override {
            calls++;
            reported = Vec3(x[0], x[1], x[2]);
            vector<Vec3> current = context.getState(State::Positions).getPositions();
            ASSERT_EQUAL_VEC(reported, current[0], 1e-6);
            return true;
        }
        Context& context;
        int calls;
        Vec3 reported;
    } reporter(context);
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000, &reporter);
    ASSERT_EQUAL(1, reporter.calls);
    vector<Vec3> finalPositions = context.getState(State::Positions).getPositions();
    ASSERT_EQUAL_VEC(reporter.reported, finalPositions[0], 1e-6);
}

/** @brief Minimize the applied force parameters, not later un-applied owner edits. */
void testMetalMinimizerAppliedParameters() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.0, 100.0);
    system.addForce(bonds);
    VerletIntegrator integrator(0.001);
    integrator.setIntegrationForceGroups(1);
    Context context(system, integrator, platform);
    bonds->setBondParameters(0, 0, 1, 2.0, 100.0);
    bonds->setForceGroup(1);
    context.setPositions({Vec3(0, 0, 0), Vec3(3, 0, 0)});
    LocalEnergyMinimizer::minimize(context, 1e-4, 1000);
    vector<Vec3> positions = context.getState(State::Positions).getPositions();
    Vec3 delta = positions[1]-positions[0];
    ASSERT_EQUAL_TOL(1.0, sqrt(delta.dot(delta)), 1e-5);

    bonds->updateParametersInContext(context);
    bonds->setBondParameters(0, 0, 1, 3.0, 100.0);
    context.setPositions({Vec3(0, 0, 0), Vec3(4, 0, 0)});
    LocalEnergyMinimizer::minimize(context, 1e-4, 1000);
    positions = context.getState(State::Positions).getPositions();
    delta = positions[1]-positions[0];
    ASSERT_EQUAL_TOL(2.0, sqrt(delta.dot(delta)), 1e-5);
}

/** @brief Reporter updates affect subsequent gradients in the same Context. */
void testMetalMinimizerReporterUpdates() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("k*(x-target)^4");
    force->addGlobalParameter("k", 1.0);
    force->addPerParticleParameter("target");
    force->addParticle(0, {0.0});
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(5, 0, 0)});
    class Reporter : public MinimizationReporter {
    public:
        Reporter(Context& context, CustomExternalForce& force) : context(context), force(force), calls(0) {
        }
        bool report(int iteration, const vector<double>& x, const vector<double>& gradient,
                map<string, double>& statistics) override {
            calls++;
            if (calls == 1) {
                context.setParameter("k", 2.0);
                force.setParticleParameters(0, 0, {1.0});
                force.updateParametersInContext(context);
                return false;
            }
            double delta = x[0]-1.0;
            ASSERT_EQUAL_TOL(8.0*delta*delta*delta, gradient[0], 2e-4);
            return true;
        }
        Context& context;
        CustomExternalForce& force;
        int calls;
    } reporter(context, *force);
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000, &reporter);
    ASSERT_EQUAL(2, reporter.calls);
    ASSERT_EQUAL(2.0, context.getParameter("k"));
}

/** @brief A reporter exception must not prevent subsequent force evaluations or minimization. */
void testMetalMinimizerReporterException() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("(x-1)^4");
    force->addParticle(0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(5, 0, 0)});
    class Reporter : public MinimizationReporter {
    public:
        bool report(int iteration, const vector<double>& x, const vector<double>& gradient,
                map<string, double>& statistics) override {
            throw OpenMMException("intentional Metal minimizer reporter failure");
        }
    } reporter;
    bool caught = false;
    try {
        LocalEnergyMinimizer::minimize(context, 1e-5, 1000, &reporter);
    }
    catch (const OpenMMException& error) {
        ASSERT_EQUAL(string("intentional Metal minimizer reporter failure"), string(error.what()));
        caught = true;
    }
    ASSERT(caught);
    State state = context.getState(State::Positions | State::Forces);
    double delta = state.getPositions()[0][0]-1.0;
    ASSERT_EQUAL_TOL(-4.0*delta*delta*delta, state.getForces()[0][0], 1e-5);
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000);
    state = context.getState(State::Forces);
    ASSERT(fabs(state.getForces()[0][0]) < 2e-5);
}

/** @brief Switch cached parent and linked CustomCV force pipelines in both directions. */
void testMetalMinimizerLinkedCustomCV() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.0, 100.0);
    CustomCVForce* cv = new CustomCVForce("scale*bond");
    cv->addGlobalParameter("scale", 2.0);
    cv->addCollectiveVariable("bond", bonds);
    system.addForce(cv);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setParameter("scale", 3.0);
    for (double distance : {2.0, 1.5}) {
        context.setPositions({Vec3(0, 0, 0), Vec3(distance, 0, 0)});
        State initial = context.getState(State::Energy | State::Forces);
        const double delta = distance-1.0;
        ASSERT_EQUAL_TOL(150.0*delta*delta, initial.getPotentialEnergy(), 1e-5);
        ASSERT_EQUAL_VEC(Vec3(300.0*delta, 0, 0), initial.getForces()[0], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(-300.0*delta, 0, 0), initial.getForces()[1], 1e-5);
        LocalEnergyMinimizer::minimize(context, 1e-4, 1000);
        State minimized = context.getState(State::Positions | State::Forces);
        Vec3 separation = minimized.getPositions()[1]-minimized.getPositions()[0];
        ASSERT_EQUAL_TOL(1.0, sqrt(separation.dot(separation)), 1e-5);
        ASSERT(sqrt(minimized.getForces()[0].dot(minimized.getForces()[0])) < 2e-4);
        ASSERT_EQUAL(3.0, context.getParameter("scale"));
    }
}

/** @brief Preserve the two displaced ATM linked contexts across repeated minimization. */
void testMetalMinimizerLinkedATM() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->addBond(0, 1, 1.0, 100.0);
    ATMForce* atm = new ATMForce("0.5*(u0+u1)");
    atm->addParticle();
    atm->addParticle(new ATMForce::FixedDisplacement(Vec3(0.5, 0, 0)));
    atm->addForce(bonds);
    system.addForce(atm);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    for (double distance : {2.0, 1.5}) {
        context.setPositions({Vec3(0, 0, 0), Vec3(distance, 0, 0)});
        State initial = context.getState(State::Energy | State::Forces);
        const double delta0 = distance-1.0;
        const double delta1 = distance-0.5;
        ASSERT_EQUAL_TOL(25.0*(delta0*delta0+delta1*delta1), initial.getPotentialEnergy(), 1e-5);
        ASSERT_EQUAL_VEC(Vec3(50.0*(delta0+delta1), 0, 0), initial.getForces()[0], 1e-5);
        ASSERT_EQUAL_VEC(Vec3(-50.0*(delta0+delta1), 0, 0), initial.getForces()[1], 1e-5);
        LocalEnergyMinimizer::minimize(context, 1e-4, 1000);
        State minimized = context.getState(State::Positions | State::Energy | State::Forces);
        Vec3 separation = minimized.getPositions()[1]-minimized.getPositions()[0];
        ASSERT_EQUAL_TOL(0.75, sqrt(separation.dot(separation)), 1e-5);
        ASSERT_EQUAL_TOL(3.125, minimized.getPotentialEnergy(), 1e-5);
        ASSERT(sqrt(minimized.getForces()[0].dot(minimized.getForces()[0])) < 2e-4);
    }
}

/** @brief Invalidate a CustomIntegrator's cached forces before resuming ordinary integration. */
void testMetalMinimizerCustomIntegratorCache() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("(x-0.75)^4");
    force->addParticle(0);
    system.addForce(force);
    CustomIntegrator integrator(0.1);
    integrator.addPerDofVariable("savedForce", 0.0);
    integrator.addComputePerDof("savedForce", "f");
    integrator.addComputePerDof("v", "v+dt*f/m");
    Context context(system, integrator, platform);
    context.setPositions({Vec3(2, 0, 0)});
    context.setVelocities({Vec3(0.5, 0, 0)});
    integrator.step(2);
    vector<Vec3> savedForces;
    integrator.getPerDofVariableByName("savedForce", savedForces);
    ASSERT_EQUAL_VEC(Vec3(-7.8125, 0, 0), savedForces[0], 1e-6);
    Vec3 velocityBefore = context.getState(State::Velocities).getVelocities()[0];
    class Reporter : public MinimizationReporter {
    public:
        bool report(int iteration, const vector<double>& x, const vector<double>& gradient,
                map<string, double>& statistics) override {
            return true;
        }
    } reporter;
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000, &reporter);
    const Vec3 position = context.getState(State::Positions).getPositions()[0];
    const double delta = position[0]-0.75;
    const Vec3 expectedForce(-4.0*delta*delta*delta, 0, 0);
    // A nonzero residual distinguishes a recomputed force from a cleared buffer.
    ASSERT(fabs(expectedForce[0]) > 1e-4);
    // Do not request State::Forces here: the next integrator step must refresh them.
    integrator.step(1);
    integrator.getPerDofVariableByName("savedForce", savedForces);
    ASSERT_EQUAL_VEC(expectedForce, savedForces[0], 1e-5);
    State resumed = context.getState(State::Positions | State::Velocities);
    ASSERT_EQUAL_VEC(position, resumed.getPositions()[0], 1e-6);
    ASSERT_EQUAL_VEC(velocityBefore+0.1*expectedForce, resumed.getVelocities()[0], 1e-5);
    ASSERT_EQUAL(3, context.getStepCount());
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000);
    integrator.step(1);
    integrator.getPerDofVariableByName("savedForce", savedForces);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), savedForces[0], 2e-5);
}

#if !OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
/** @brief Reject unrepresentable forces instead of silently converging or using a CPU force path. */
template <class Test>
void testMetalMinimizerOverflow(Test test) {
    bool caught = false;
    try {
        test();
    }
    catch (const OpenMMException& error) {
        const string message = error.what();
        ASSERT(message.find("Minimization exceeded the GPU ") == 0);
        ASSERT(message.find(" range") != string::npos);
        ASSERT(message.find("CPU fallback is disabled") != string::npos);
        caught = true;
    }
    ASSERT(caught);
}

/** @brief Reset sticky range diagnostics before reusing a Context after an overflow. */
void testMetalMinimizerRangeRecovery() {
    System system;
    system.addParticle(1.0);
    CustomExternalForce* force = new CustomExternalForce("0.5*k*(x-0.5)^2");
    force->addGlobalParameter("k", 1e22);
    force->addParticle(0);
    system.addForce(force);
    VerletIntegrator integrator(0.001);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(2, 0, 0)});
    testMetalMinimizerOverflow([&]() { LocalEnergyMinimizer::minimize(context, 1.0, 100); });
    context.setParameter("k", 4.0);
    context.setPositions({Vec3(2, 0, 0)});
    LocalEnergyMinimizer::minimize(context, 1e-5, 1000);
    State recovered = context.getState(State::Positions | State::Forces);
    ASSERT_EQUAL_VEC(Vec3(0.5, 0, 0), recovered.getPositions()[0], 1e-5);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), recovered.getForces()[0], 2e-5);
}

/** @brief Minimize bounded nonbonded forces with a constraint and a virtual site. */
void testMetalFixedMinimizerConstrainedVirtualSite() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.addConstraint(0, 1, 1.0);
    system.setVirtualSite(2, new TwoParticleAverageSite(0, 1, 0.5, 0.5));
    NonbondedForce* nonbonded = new NonbondedForce();
    // Put the LJ minimum at the constrained distance so the positive fixture
    // does not require cancellation against a stiff constraint restraint.
    const double sigma = pow(0.5, 1.0/6.0);
    nonbonded->addParticle(0.0, sigma, 0.1);
    nonbonded->addParticle(0.0, sigma, 0.1);
    nonbonded->addParticle(0.0, 0.2, 0.0);
    system.addForce(nonbonded);
    CustomExternalForce* external = new CustomExternalForce("5*x*x");
    external->addParticle(2);
    system.addForce(external);
    VerletIntegrator integrator(0.001);
    integrator.setConstraintTolerance(1e-5);
    Context context(system, integrator, platform);
    context.setPositions({Vec3(1.5, 0, 0), Vec3(2.5, 0, 0), Vec3(2, 0, 0)});
    double initialEnergy = context.getState(State::Energy).getPotentialEnergy();
    LocalEnergyMinimizer::minimize(context, 1e-4, 1000);
    State minimized = context.getState(State::Positions | State::Forces | State::Energy);
    const vector<Vec3>& positions = minimized.getPositions();
    const Vec3 separation = positions[1]-positions[0];
    const Vec3 midpoint = 0.5*(positions[0]+positions[1]);
    ASSERT_EQUAL_TOL(1.0, sqrt(separation.dot(separation)), 2e-5);
    ASSERT_EQUAL_VEC(midpoint, positions[2], 1e-6);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), midpoint, 2e-5);
    ASSERT(minimized.getPotentialEnergy() < initialEnergy);
    ASSERT_EQUAL_TOL(-0.1, minimized.getPotentialEnergy(), 1e-5);
    ASSERT_EQUAL_VEC(-minimized.getForces()[0], minimized.getForces()[1], 2e-4);
    ASSERT(sqrt(minimized.getForces()[0].dot(minimized.getForces()[0])) < 10.0);
}

/** @brief Run the Common minimizer checks with explicit fixed-point overflow expectations. */
int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        int failures = 0;
        auto run = [&](const char* name, void (*test)()) {
            cout << "Fixed-point " << name << endl;
            try {
                test();
            }
            catch (const exception& error) {
                cerr << name << ": " << error.what() << endl;
                failures++;
            }
        };
        run("testHarmonicBonds", testHarmonicBonds);
        // These unchanged seeded random-cloud cases encounter trial forces
        // beyond Q32.32 or Common's stricter force safety guard. The ON build
        // retains their original success assertions; this limited-range build
        // must reject them instead of taking Common's CPU recovery path.
        run("testLargeSystem (expected overflow)", []() { testMetalMinimizerOverflow(testLargeSystem); });
        run("testVirtualSites (expected overflow)", []() { testMetalMinimizerOverflow(testVirtualSites); });
        run("testLargeForces (expected overflow)", []() { testMetalMinimizerOverflow(testLargeForces); });
        run("bounded nonbonded/constraint/virtual-site case", testMetalFixedMinimizerConstrainedVirtualSite);
        run("testForceGroups", testForceGroups);
        run("testMasslessParticles", testMasslessParticles);
        run("testReporter", testReporter);
        run("Metal regression cases", runPlatformTests);
        if (failures != 0)
            return 1;
    }
    catch (const exception& error) {
        cerr << "exception: " << error.what() << endl;
        return 1;
    }
    cout << "Done (fixed-point minimization; extreme-force cases rejected)" << endl;
    return 0;
}
#endif
#endif

#if defined(METAL_TEST_CONSTANT_POTENTIAL)

/** @brief No platform-specific setup is needed beyond the common test initialization. */
void platformInitialize() {
}

/** @brief Retain the OpenCL single-device ConstantPotentialForce regression coverage. */
void runPlatformTests(ConstantPotentialForce::ConstantPotentialMethod method, bool usePreconditioner) {
    testEnergyConservation(method, usePreconditioner, 10);
    testCompareToReferencePlatform(method, usePreconditioner);
    testLargeNeighborList(method, usePreconditioner);
}

#else

/** @brief Run the useful OpenCL extra checks without its multi-device cases. */
void runPlatformTests() {
#if defined(METAL_TEST_NONBONDED)
    testReordering();
    testMetalDenseNeighborList();
#elif defined(METAL_TEST_MONTE_CARLO_BAROSTAT)
    testWater();
    testLJPressure();
#elif defined(METAL_TEST_MONTE_CARLO_ANISOTROPIC_BAROSTAT)
    testLJPressure();
#elif defined(METAL_TEST_MINIMIZER)
#if OPENMM_METAL_MINIMIZE_FLOAT_ACCUMULATORS
    testMetalLargeBondMinimization();
    testMetalLargeCustomMinimizationState();
    testMetalLargeVirtualSiteMinimization();
#else
    testMetalMinimizerOverflow(testMetalLargeBondMinimization);
    testMetalMinimizerOverflow(testMetalLargeCustomMinimizationState);
    testMetalMinimizerOverflow(testMetalLargeVirtualSiteMinimization);
    testMetalMinimizerRangeRecovery();
#endif
    testMetalMinimizerReporterContext();
    testMetalMinimizerAppliedParameters();
    testMetalMinimizerReporterUpdates();
    testMetalMinimizerReporterException();
    testMetalMinimizerLinkedCustomCV();
    testMetalMinimizerLinkedATM();
    testMetalMinimizerCustomIntegratorCache();
#endif
}

#endif
