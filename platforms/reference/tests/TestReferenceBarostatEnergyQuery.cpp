#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/NonbondedForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/LangevinMiddleIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using namespace OpenMM;
class DerivedMiddle : public LangevinMiddleIntegrator {
public: DerivedMiddle() : LangevinMiddleIntegrator(300, 1, 0.001) {}
};
static int changedBoxes = 0;
static std::vector<double> simulate(int kind, bool enabled, bool rigid, int groups, double pressure) {
#ifdef _WIN32
    _putenv_s("OPENMM_EXPERIMENT_BAROSTAT_POTENTIAL_ONLY", enabled ? "1" : "0");
#else
    setenv("OPENMM_EXPERIMENT_BAROSTAT_POTENTIAL_ONLY", enabled ? "1" : "0", 1);
#endif
    System system;
    auto* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.0);
    nonbonded->setForceGroup(1);
    auto* bonds = new HarmonicBondForce();
    bonds->setForceGroup(2);
    std::vector<Vec3> positions, velocities;
    for (int i=0; i<12; ++i) {
        system.addParticle(18.0);
        nonbonded->addParticle(0, 0.3, 0.05);
        positions.push_back(Vec3(0.7*(i%3), 0.7*((i/3)%2), 0.7*(i/6)));
        velocities.push_back(Vec3(0.001*i, -0.002*i, 0.003*i));
        if (i%3 == 1) bonds->addBond(i-1, i, 0.7, 10.0);
    }
    system.addForce(nonbonded);
    system.addForce(bonds);
    system.setDefaultPeriodicBoxVectors(Vec3(3.5,0,0), Vec3(0,3.5,0), Vec3(0,0,3.5));
    auto* barostat = new MonteCarloBarostat(pressure, 300, 1);
    barostat->setRandomNumberSeed(44221);
    barostat->setScaleMoleculesAsRigid(rigid);
    system.addForce(barostat);
    std::unique_ptr<Integrator> integrator;
    if (kind == 0) integrator.reset(new LangevinMiddleIntegrator(300,1,0.001));
    else if (kind == 1) integrator.reset(new DerivedMiddle());
    else integrator.reset(new VerletIntegrator(0.001));
    if (kind != 2) dynamic_cast<LangevinMiddleIntegrator&>(*integrator).setRandomNumberSeed(99231);
    integrator->setIntegrationForceGroups(groups);
    Context context(system,*integrator,Platform::getPlatformByName("Reference"));
    context.setPositions(positions);
    context.setVelocities(velocities);
    std::vector<double> record;
    double previousBox = 3.5;
    for (int step=0; step<40; ++step) {
        integrator->step(1);
        auto state=context.getState(State::Positions|State::Velocities|State::Energy,false,groups);
        for (const auto& p:state.getPositions()) for (int j=0;j<3;++j) record.push_back(p[j]);
        for (const auto& v:state.getVelocities()) for (int j=0;j<3;++j) record.push_back(v[j]);
        Vec3 a,b,c;state.getPeriodicBoxVectors(a,b,c);
        for (auto v : {a,b,c}) for (int j=0;j<3;++j) record.push_back(v[j]);
        record.push_back(state.getPotentialEnergy());
        record.push_back(state.getKineticEnergy());
        if (a[0]!=previousBox) ++changedBoxes;
        previousBox=a[0];
    }
    return record;
}
int main() {
    try {
        int cases=0;
        for (int kind=0;kind<3;++kind) for (bool rigid : {false,true})
            for (int groups : {2,6}) for (double pressure : {2.0,10000.0}) {
                auto off=simulate(kind,false,rigid,groups,pressure);
                auto on=simulate(kind,true,rigid,groups,pressure);
                if (off.size()!=on.size() || std::memcmp(off.data(),on.data(),off.size()*sizeof(double))!=0)
                    throw std::runtime_error("Enabled and disabled trajectories differ");
                ++cases;
            }
        if (!changedBoxes) throw std::runtime_error("No accepted volume changes were exercised");
        std::cout<<"{\"status\":\"PASS\",\"paired_cases\":"<<cases<<",\"total_steps\":1920,\"trajectory_bytes_equal\":true,\"box_changes\":"<<changedBoxes<<",\"platform\":\"Reference\",\"gpu_executed\":false}"<<std::endl;
    } catch (const std::exception& error) { std::cerr<<error.what()<<std::endl;return 1; }
}
