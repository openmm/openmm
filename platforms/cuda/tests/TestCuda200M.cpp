#include <fstream>
#include <iostream>
#include "CudaTests.h"
#include "OpenMM.h"
using namespace OpenMM;
using namespace std;

int main(int argc, char** argv) {
    initializeTests(argc, argv);

    auto integrator = VerletIntegrator(0.001);
    auto system = System();
    auto* force = new NonbondedForce();

    force->setNonbondedMethod(NonbondedForce::NonbondedMethod::CutoffPeriodic);
    force->setCutoffDistance(1.3);

    system.setDefaultPeriodicBoxVectors(
        Vec3(2150.4 / 10,0,0),
        Vec3(0,1075.2 / 10,0),
        Vec3(0,0,1075.2 / 10));

    vector<Vec3> pos, vel;

    ifstream in("/home/pavlov/work/lj_test/dump.10061824.xyz.100");
    string s;
    getline(in, s);
    int N = stoi(s);
    getline(in, s);
    cout << N << s << endl;
    for(int i = 0; i < N; i++) {
        ///*
        int t;
        double x, y, z;
        in >> t >> x >> y >> z;
        if(in.fail()) {
            cout << "Input failed!" << endl;
            return 0;
        }
        pos.emplace_back(x / 10, y / 10, z / 10);
        //*/
        system.addParticle(39.948);
        force->addParticle(0, 0.3405, 0.2381 * 4.184);
    }
    in.close();
    system.addForce(force);
    cout << "File imported!" << endl;
    Context ctx(system, integrator, platform, {{"Precision", "mixed"}});
    cout << "Context created!" << endl;

    //ifstream checkp("1.checkpoint");
    //ctx.loadCheckpoint(checkp);
    //cout << "Loaded checkpoint!" << endl;

    ctx.setPositions(pos);
    ctx.setVelocitiesToTemperature(1500);

    cout << "Simulation set up!" << endl;
    auto state = ctx.getState(State::Energy);
    cout << ctx.getStepCount() << ' ' << state.getKineticEnergy() << ' ' << state.getPotentialEnergy() << endl;
    for(int i = ctx.getStepCount(), delta = 100; i < 50000; i = ctx.getStepCount()) {
        integrator.step(delta);
        int step = ctx.getStepCount();
        //ofstream out(to_string(step) + ".checkpoint");
        //ctx.createCheckpoint(out);
        cout << "Step taken!" << endl;
        auto state = ctx.getState(State::Energy);
        cout << step << ' ' << state.getKineticEnergy() << ' ' << state.getPotentialEnergy() << endl;
    }
    cout << "All done! Destructing!" << endl;
}
