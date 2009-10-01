/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "OpenCLKernels.h"
#include "OpenCLForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLNonbondedUtilities.h"
#include <cmath>

using namespace OpenMM;
using namespace std;

static const double KILO = 1e3;                      // Thousand
static const double BOLTZMANN = 1.380658e-23;            // (J/K)
static const double AVOGADRO = 6.0221367e23;            // ()
static const double RGAS = BOLTZMANN*AVOGADRO;     // (J/(mol K))
static const double BOLTZ = (RGAS/KILO);            // (kJ/(mol K))

void OpenCLCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void OpenCLCalcForcesAndEnergyKernel::beginForceComputation(ContextImpl& context) {
    cl.clearBuffer(cl.getForceBuffers());
    cl.getNonbondedUtilties().prepareInteractions();
}

void OpenCLCalcForcesAndEnergyKernel::finishForceComputation(ContextImpl& context) {
    cl.reduceBuffer(cl.getForceBuffers(), cl.getNumForceBuffers());
    cl.getNonbondedUtilties().prepareInteractions();
}

void OpenCLCalcForcesAndEnergyKernel::beginEnergyComputation(ContextImpl& context) {
    cl.clearBuffer(cl.getEnergyBuffer());
}

double OpenCLCalcForcesAndEnergyKernel::finishEnergyComputation(ContextImpl& context) {
    OpenCLArray<cl_float>& energy = cl.getEnergyBuffer();
    energy.download();
    double sum = 0.0f;
    for (int i = 0; i < energy.getSize(); i++)
        sum += energy[i];
    return sum;
}

void OpenCLUpdateStateDataKernel::initialize(const System& system) {
}

double OpenCLUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cl.getTime();
}

void OpenCLUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    cl.setTime(time);
}

void OpenCLUpdateStateDataKernel::getPositions(ContextImpl& context, std::vector<Vec3>& positions) {
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    posq.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 pos = posq[i];
//        int3 offset = gpu->posCellOffsets[i];
//        positions[order[i]] = Vec3(pos.x-offset.x*gpu->sim.periodicBoxSizeX, pos.y-offset.y*gpu->sim.periodicBoxSizeY, pos.z-offset.z*gpu->sim.periodicBoxSizeZ);
        positions[order[i]] = Vec3(pos.x, pos.y, pos.z);
    }
}

void OpenCLUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4& pos = posq[i];
        const Vec3& p = positions[order[i]];
        pos.x = p[0];
        pos.y = p[1];
        pos.z = p[2];
    }
    posq.upload();
//    for (int i = 0; i < gpu->posCellOffsets.size(); i++)
//        gpu->posCellOffsets[i] = make_int3(0, 0, 0);
}

void OpenCLUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    OpenCLArray<mm_float4>& velm = cl.getVelm();
    velm.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 vel = velm[i];
        velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
    }
}

void OpenCLUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    OpenCLArray<mm_float4>& velm = cl.getVelm();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4& vel = velm[i];
        const Vec3& p = velocities[order[i]];
        vel.x = p[0];
        vel.y = p[1];
        vel.z = p[2];
    }
    velm.upload();
}

void OpenCLUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    OpenCLArray<mm_float4>& force = cl.getForce();
    force.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 f = force[i];
        forces[order[i]] = Vec3(f.x, f.y, f.z);
    }
}

class OpenCLBondForceInfo : public OpenCLForceInfo {
public:
    OpenCLBondForceInfo(int requiredBuffers, const HarmonicBondForce& force) : OpenCLForceInfo(requiredBuffers, false, 0.0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getBondParameters(group1, particle1, particle2, length1, k1);
        force.getBondParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const HarmonicBondForce& force;
};

OpenCLCalcHarmonicBondForceKernel::~OpenCLCalcHarmonicBondForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    numBonds = force.getNumBonds();
    params = new OpenCLArray<mm_float2>(cl, numBonds, "bondParams");
    indices = new OpenCLArray<mm_int4>(cl, numBonds, "bondIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float2> paramVector(numBonds);
    vector<mm_int4> indicesVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        paramVector[i] = (mm_float2) {length, k};
        indicesVector[i] = (mm_int4) {particle1, particle2, forceBufferCounter[particle1]++, forceBufferCounter[particle2]++};

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < forceBufferCounter.size(); i++) {
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLBondForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("harmonicBondForce.cl"));
    kernel = cl::Kernel(program, "calcHarmonicBondForce");
}

void OpenCLCalcHarmonicBondForceKernel::executeForces(ContextImpl& context) {
    kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
    kernel.setArg<cl_int>(1, numBonds);
    kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    cl.executeKernel(kernel, numBonds);
}

double OpenCLCalcHarmonicBondForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLAngleForceInfo : public OpenCLForceInfo {
public:
    OpenCLAngleForceInfo(int requiredBuffers, const HarmonicAngleForce& force) : OpenCLForceInfo(requiredBuffers, false, 0.0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(index, particle1, particle2, particle3, angle, k);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        double angle1, angle2, k1, k2;
        force.getAngleParameters(group1, particle1, particle2, particle3, angle1, k1);
        force.getAngleParameters(group2, particle1, particle2, particle3, angle2, k2);
        return (angle1 == angle2 && k1 == k2);
    }
private:
    const HarmonicAngleForce& force;
};

OpenCLCalcHarmonicAngleForceKernel::~OpenCLCalcHarmonicAngleForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    params = new OpenCLArray<mm_float2>(cl, numAngles, "angleParams");
    indices = new OpenCLArray<mm_int8>(cl, numAngles, "angleIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float2> paramVector(numAngles);
    vector<mm_int8> indicesVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        paramVector[i] = (mm_float2) {angle, k};
        indicesVector[i] = (mm_int8) {particle1, particle2, particle3,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, 0, 0};

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < forceBufferCounter.size(); i++) {
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLAngleForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("harmonicAngleForce.cl"));
    kernel = cl::Kernel(program, "calcHarmonicAngleForce");
}

void OpenCLCalcHarmonicAngleForceKernel::executeForces(ContextImpl& context) {
    kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
    kernel.setArg<cl_int>(1, numAngles);
    kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    cl.executeKernel(kernel, numAngles);
}

double OpenCLCalcHarmonicAngleForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLPeriodicTorsionForceInfo : public OpenCLForceInfo {
public:
    OpenCLPeriodicTorsionForceInfo(int requiredBuffers, const PeriodicTorsionForce& force) : OpenCLForceInfo(requiredBuffers, false, 0.0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, periodicity, phase, k);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, periodicity1, periodicity2;
        double phase1, phase2, k1, k2;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, periodicity1, phase1, k1);
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, periodicity2, phase2, k2);
        return (periodicity1 == periodicity2 && phase1 == phase2 && k1 == k2);
    }
private:
    const PeriodicTorsionForce& force;
};

OpenCLCalcPeriodicTorsionForceKernel::~OpenCLCalcPeriodicTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    params = new OpenCLArray<mm_float4>(cl, numTorsions, "periodicTorsionParams");
    indices = new OpenCLArray<mm_int8>(cl, numTorsions, "periodicTorsionIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float4> paramVector(numTorsions);
    vector<mm_int8> indicesVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        paramVector[i] = (mm_float4) {k, phase, (float) periodicity};
        indicesVector[i] = (mm_int8) {particle1, particle2, particle3, particle4,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, forceBufferCounter[particle4]++};

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < forceBufferCounter.size(); i++) {
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLPeriodicTorsionForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("periodicTorsionForce.cl"));
    kernel = cl::Kernel(program, "calcPeriodicTorsionForce");
}

void OpenCLCalcPeriodicTorsionForceKernel::executeForces(ContextImpl& context) {
    kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
    kernel.setArg<cl_int>(1, numTorsions);
    kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    cl.executeKernel(kernel, numTorsions);
}

double OpenCLCalcPeriodicTorsionForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLRBTorsionForceInfo : public OpenCLForceInfo {
public:
    OpenCLRBTorsionForceInfo(int requiredBuffers, const RBTorsionForce& force) : OpenCLForceInfo(requiredBuffers, false, 0.0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        double c0a, c0b, c1a, c1b, c2a, c2b, c3a, c3b, c4a, c4b, c5a, c5b;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, c0a, c1a, c2a, c3a, c4a, c5a);
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, c0b, c1b, c2b, c3b, c4b, c5b);
        return (c0a == c0b && c1a == c1b && c2a == c2b && c3a == c3b && c4a == c4b && c5a == c5b);
    }
private:
    const RBTorsionForce& force;
};

OpenCLCalcRBTorsionForceKernel::~OpenCLCalcRBTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    params = new OpenCLArray<mm_float8>(cl, numTorsions, "rbTorsionParams");
    indices = new OpenCLArray<mm_int8>(cl, numTorsions, "rbTorsionIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float8> paramVector(numTorsions);
    vector<mm_int8> indicesVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        paramVector[i] = (mm_float8) {c0, c1, c2, c3, c4, c5};
        indicesVector[i] = (mm_int8) {particle1, particle2, particle3, particle4,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, forceBufferCounter[particle4]++};

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < forceBufferCounter.size(); i++) {
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLRBTorsionForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("rbTorsionForce.cl"));
    kernel = cl::Kernel(program, "calcRBTorsionForce");
}

void OpenCLCalcRBTorsionForceKernel::executeForces(ContextImpl& context) {
    kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
    kernel.setArg<cl_int>(1, numTorsions);
    kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
    kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    cl.executeKernel(kernel, numTorsions);
}

double OpenCLCalcRBTorsionForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

OpenCLCalcNonbondedForceKernel::~OpenCLCalcNonbondedForceKernel() {
    if (sigmaEpsilon != NULL)
        delete sigmaEpsilon;
}

void OpenCLCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    vector<pair<int, int> > exclusions;
    vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0 || epsilon != 0.0)
            exceptions.push_back(i);
    }

    // Initialize nonbonded interactions.

    int numParticles = force.getNumParticles();
    sigmaEpsilon = new OpenCLArray<mm_float2>(cl, numParticles, "sigmaEpsilon");
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    vector<mm_float2> sigmaEpsilonVector(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        posq[i].w = (float) charge;
        sigmaEpsilonVector[i] = (mm_float2) {(float) (0.5*sigma), (float) (2.0*sqrt(epsilon))};
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < (int) exclusions.size(); i++) {
        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
    }
    posq.upload();
    sigmaEpsilon->upload(sigmaEpsilonVector);
    bool useCutoff = (force.getNonbondedMethod() != NonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != NonbondedForce::NoCutoff && force.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic);
//    if (force.getNonbondedMethod() != NonbondedForce::NoCutoff) {
//        method = CUTOFF;
//    }
//    if (force.getNonbondedMethod() == NonbondedForce::CutoffPeriodic) {
//        method = PERIODIC;
//    }
//    if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
//        double ewaldErrorTol = force.getEwaldErrorTolerance();
//        double alpha = (1.0/force.getCutoffDistance())*std::sqrt(-std::log(ewaldErrorTol));
//        double mx = boxVectors[0][0]/force.getCutoffDistance();
//        double my = boxVectors[1][1]/force.getCutoffDistance();
//        double mz = boxVectors[2][2]/force.getCutoffDistance();
//        double pi = 3.1415926535897932385;
//        int kmaxx = (int)std::ceil(-(mx/pi)*std::log(ewaldErrorTol));
//        int kmaxy = (int)std::ceil(-(my/pi)*std::log(ewaldErrorTol));
//        int kmaxz = (int)std::ceil(-(mz/pi)*std::log(ewaldErrorTol));
//        if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
//            if (kmaxx%2 == 0)
//                kmaxx++;
//            if (kmaxy%2 == 0)
//                kmaxy++;
//            if (kmaxz%2 == 0)
//                kmaxz++;
//            gpuSetEwaldParameters(gpu, (float) alpha, kmaxx, kmaxy, kmaxz);
//            method = EWALD;
//        }
//        else {
//            int gridSizeX = -0.5*kmaxx*std::log(ewaldErrorTol);
//            int gridSizeY = -0.5*kmaxy*std::log(ewaldErrorTol);
//            int gridSizeZ = -0.5*kmaxz*std::log(ewaldErrorTol);
//            gpuSetPMEParameters(gpu, (float) alpha, gridSizeX, gridSizeY, gridSizeZ);
//            method = PARTICLE_MESH_EWALD;
//        }
//    }
//    data.nonbondedMethod = method;
//    gpuSetCoulombParameters(gpu, 138.935485f, particle, c6, c12, q, symbol, exclusionList, method);
    cl.getNonbondedUtilties().addInteraction(useCutoff, usePeriodic, force.getCutoffDistance(), exclusionList);
    cl.getNonbondedUtilties().addParameter("sigmaEpsilon", "float2", 8, sigmaEpsilon->getDeviceBuffer());

    // Compute the Ewald self energy.

    ewaldSelfEnergy = 0.0;
    if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
//        double selfEnergyScale = gpu->sim.epsfac*gpu->sim.alphaEwald/std::sqrt(PI);
//            for (int i = 0; i < numParticles; i++)
//                ewaldSelfEnergy -= selfEnergyScale*q[i]*q[i];
    }

    // Initialize 1-4 nonbonded interactions.

    {
        int numExceptions = exceptions.size();
        vector<int> particle1(numExceptions);
        vector<int> particle2(numExceptions);
        vector<float> c6(numExceptions);
        vector<float> c12(numExceptions);
        vector<float> q1(numExceptions);
        vector<float> q2(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double charge, sig, eps;
            force.getExceptionParameters(exceptions[i], particle1[i], particle2[i], charge, sig, eps);
            c6[i] = (float) (4*eps*pow(sig, 6.0));
            c12[i] = (float) (4*eps*pow(sig, 12.0));
            q1[i] = (float) charge;
            q2[i] = 1.0f;
        }
//        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, particle1, particle2, c6, c12, q1, q2);
    }
}

void OpenCLCalcNonbondedForceKernel::executeForces(ContextImpl& context) {
    cl.getNonbondedUtilties().computeInteractions();
}

double OpenCLCalcNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return ewaldSelfEnergy;
}

//OpenCLCalcCustomNonbondedForceKernel::~OpenCLCalcCustomNonbondedForceKernel() {
//}
//
//void OpenCLCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
//    data.primaryKernel = this; // This must always be the primary kernel so it can update the global parameters
//    data.hasCustomNonbonded = true;
//    numParticles = force.getNumParticles();
//    _gpuContext* gpu = data.gpu;
//
//    // Identify which exceptions are actual interactions.
//
//    vector<pair<int, int> > exclusions;
//    vector<int> exceptions;
//    {
//        vector<double> parameters;
//        for (int i = 0; i < force.getNumExceptions(); i++) {
//            int particle1, particle2;
//            force.getExceptionParameters(i, particle1, particle2, parameters);
//            exclusions.push_back(pair<int, int>(particle1, particle2));
//            if (parameters.size() > 0)
//                exceptions.push_back(i);
//        }
//    }
//
//    // Initialize nonbonded interactions.
//
//    vector<int> particle(numParticles);
//    vector<vector<double> > parameters(numParticles);
//    vector<vector<int> > exclusionList(numParticles);
//    for (int i = 0; i < numParticles; i++) {
//        force.getParticleParameters(i, parameters[i]);
//        particle[i] = i;
//        exclusionList[i].push_back(i);
//    }
//    for (int i = 0; i < (int)exclusions.size(); i++) {
//        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
//        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
//    }
//    Vec3 boxVectors[3];
//    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
//    gpuSetPeriodicBoxSize(gpu, (float)boxVectors[0][0], (float)boxVectors[1][1], (float)boxVectors[2][2]);
//    OpenCLNonbondedMethod method = NO_CUTOFF;
//    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
//        method = CUTOFF;
//    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic) {
//        method = PERIODIC;
//    }
//    data.customNonbondedMethod = method;
//
//    // Initialize exceptions.
//
//    int numExceptions = exceptions.size();
//    vector<int> exceptionParticle1(numExceptions);
//    vector<int> exceptionParticle2(numExceptions);
//    vector<vector<double> > exceptionParams(numExceptions);
//    for (int i = 0; i < numExceptions; i++)
//        force.getExceptionParameters(exceptions[i], exceptionParticle1[i], exceptionParticle2[i], exceptionParams[i]);
//
//    // Record the tabulated functions.
//
//    for (int i = 0; i < force.getNumFunctions(); i++) {
//        string name;
//        vector<double> values;
//        double min, max;
//        bool interpolating;
//        force.getFunctionParameters(i, name, values, min, max, interpolating);
//        gpuSetTabulatedFunction(gpu, i, name, values, min, max, interpolating);
//    }
//
//    // Record information for the expressions.
//
//    vector<string> paramNames;
//    vector<string> combiningRules;
//    for (int i = 0; i < force.getNumParameters(); i++) {
//        paramNames.push_back(force.getParameterName(i));
//        combiningRules.push_back(force.getParameterCombiningRule(i));
//    }
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = force.getGlobalParameterDefaultValue(i);
//    }
//    gpuSetCustomNonbondedParameters(gpu, parameters, exclusionList, exceptionParticle1, exceptionParticle2, exceptionParams, method,
//            (float)force.getCutoffDistance(), force.getEnergyFunction(), combiningRules, paramNames, globalParamNames);
//    if (globalParamValues.size() > 0)
//        SetCustomNonbondedGlobalParams(&globalParamValues[0]);
//}
//
//void OpenCLCalcCustomNonbondedForceKernel::executeForces(ContextImpl& context) {
//    if (data.primaryKernel == this) {
//        updateGlobalParams(context);
//        calcForces(context, data);
//    }
//}
//
//double OpenCLCalcCustomNonbondedForceKernel::executeEnergy(ContextImpl& context) {
//    if (data.primaryKernel == this) {
//        updateGlobalParams(context);
//        return calcEnergy(context, data, system);
//    }
//    return 0.0;
//}
//
//void OpenCLCalcCustomNonbondedForceKernel::updateGlobalParams(ContextImpl& context) {
//    bool changed = false;
//    for (int i = 0; i < globalParamNames.size(); i++) {
//        float value = (float) context.getParameter(globalParamNames[i]);
//        if (value != globalParamValues[i])
//            changed = true;
//        globalParamValues[i] = value;
//    }
//    if (changed)
//        SetCustomNonbondedGlobalParams(&globalParamValues[0]);
//}
//
//OpenCLCalcGBSAOBCForceKernel::~OpenCLCalcGBSAOBCForceKernel() {
//}
//
//void OpenCLCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
//
//    int numParticles = system.getNumParticles();
//    _gpuContext* gpu = data.gpu;
//    vector<float> radius(numParticles);
//    vector<float> scale(numParticles);
//    vector<float> charge(numParticles);
//    for (int i = 0; i < numParticles; i++) {
//        double particleCharge, particleRadius, scalingFactor;
//        force.getParticleParameters(i, particleCharge, particleRadius, scalingFactor);
//        radius[i] = (float) particleRadius;
//        scale[i] = (float) scalingFactor;
//        charge[i] = (float) particleCharge;
//    }
//    gpuSetObcParameters(gpu, (float) force.getSoluteDielectric(), (float) force.getSolventDielectric(), radius, scale, charge);
//}
//
//void OpenCLCalcGBSAOBCForceKernel::executeForces(ContextImpl& context) {
//}
//
//static void initializeIntegration(const System& system, OpenCLPlatform::PlatformData& data, const Integrator& integrator) {
//
//    // Initialize any terms that haven't already been handled by a Force.
//
//    _gpuContext* gpu = data.gpu;
//    if (!data.hasBonds)
//        gpuSetBondParameters(gpu, vector<int>(), vector<int>(), vector<float>(), vector<float>());
//    if (!data.hasAngles)
//        gpuSetBondAngleParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>());
//    if (!data.hasPeriodicTorsions)
//        gpuSetDihedralParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<int>());
//    if (!data.hasRB)
//        gpuSetRbDihedralParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>(),
//                vector<float>(), vector<float>(), vector<float>(), vector<float>());
//    if (!data.hasNonbonded) {
//        gpuSetCoulombParameters(gpu, 138.935485f, vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<char>(), vector<vector<int> >(), NO_CUTOFF);
//        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<float>());
//    }
//
//    // Set masses.
//
//    int numParticles = system.getNumParticles();
//    vector<float> mass(numParticles);
//    for (int i = 0; i < numParticles; i++)
//        mass[i] = (float) system.getParticleMass(i);
//    gpuSetMass(gpu, mass);
//
//    // Set constraints.
//
//    int numConstraints = system.getNumConstraints();
//    vector<int> particle1(numConstraints);
//    vector<int> particle2(numConstraints);
//    vector<float> distance(numConstraints);
//    vector<float> invMass1(numConstraints);
//    vector<float> invMass2(numConstraints);
//    for (int i = 0; i < numConstraints; i++) {
//        int particle1Index, particle2Index;
//        double constraintDistance;
//        system.getConstraintParameters(i, particle1Index, particle2Index, constraintDistance);
//        particle1[i] = particle1Index;
//        particle2[i] = particle2Index;
//        distance[i] = (float) constraintDistance;
//        invMass1[i] = 1.0f/mass[particle1Index];
//        invMass2[i] = 1.0f/mass[particle2Index];
//    }
//    gpuSetConstraintParameters(gpu, particle1, particle2, distance, invMass1, invMass2, (float)integrator.getConstraintTolerance());
//
//    // Finish initialization.
//
//    gpuBuildThreadBlockWorkList(gpu);
//    gpuBuildExclusionList(gpu);
//    gpuBuildOutputBuffers(gpu);
//    gpuSetConstants(gpu);
//    kClearBornForces(gpu);
//    kClearForces(gpu);
//    cudaThreadSynchronize();
//}
//
//double OpenCLCalcGBSAOBCForceKernel::executeEnergy(ContextImpl& context) {
//	return 0.0;
//}

OpenCLIntegrateVerletStepKernel::~OpenCLIntegrateVerletStepKernel() {
}

void OpenCLIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cl.initialize(system);
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("verlet.cl"));
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
}

void OpenCLIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilties();
    int numAtoms = cl.getNumAtoms();
    double dt = integrator.getStepSize();

    // Call the first integration kernel.

    kernel1.setArg<cl_int>(0, numAtoms);
    kernel1.setArg<cl_float>(1, dt);
    kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(4, cl.getForce().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2.setArg<cl_int>(0, numAtoms);
    kernel2.setArg<cl_float>(1, dt);
    kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(4, integration.getPosDelta().getDeviceBuffer());
    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLIntegrateLangevinStepKernel::~OpenCLIntegrateLangevinStepKernel() {
    if (params != NULL)
        delete params;
    if (xVector != NULL)
        delete xVector;
    if (vVector != NULL)
        delete vVector;
}

void OpenCLIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cl.initialize(system);
    cl.getIntegrationUtilties().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("langevin.cl"));
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    kernel3 = cl::Kernel(program, "integrateLangevinPart3");
    params = new OpenCLArray<cl_float>(cl, 11, "langevinParams");
    xVector = new OpenCLArray<mm_float4>(cl, cl.getPaddedNumAtoms(), "xVector");
    vVector = new OpenCLArray<mm_float4>(cl, cl.getPaddedNumAtoms(), "vVector");
    vector<mm_float4> initialXVector(xVector->getSize(), (mm_float4) {0.0f, 0.0f, 0.0f, 0.0f});
    xVector->upload(initialXVector);
    prevStepSize = -1.0;
}

void OpenCLIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilties();
    int numAtoms = cl.getNumAtoms();
    int numThreads = cl.getNumThreadBlocks()*cl.ThreadBlockSize;
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        double kT = BOLTZ*temperature;
        double GDT = stepSize/tau;
        double EPH = exp(0.5*GDT);
        double EMH = exp(-0.5*GDT);
        double EP = exp(GDT);
        double EM = exp(-GDT);
        double B, C, D;
        if (GDT >= 0.1)
        {
            double term1 = EPH - 1.0;
            term1 *= term1;
            B = GDT*(EP - 1.0) - 4.0*term1;
            C = GDT - 3.0 + 4.0*EMH - EM;
            D = 2.0 - EPH - EMH;
        }
        else
        {
            double term1 = 0.5*GDT;
            double term2 = term1*term1;
            double term4 = term2*term2;

            double third = 1.0/3.0;
            double o7_9 = 7.0/9.0;
            double o1_12 = 1.0/12.0;
            double o17_90 = 17.0/90.0;
            double o7_30 = 7.0/30.0;
            double o31_1260 = 31.0/1260.0;
            double o_360 = 1.0/360.0;
            B = term4*(third + term1*(third + term1*(o17_90 + term1*o7_9)));
            C = term2*term1*(2.0*third + term1*(-0.5 + term1*(o7_30 + term1*(-o1_12 + term1*o31_1260))));
            D = term2*(-1.0 + term2*(-o1_12 - term2*o_360));
        }
        double DOverTauC = D/(tau*C);
        double TauOneMinusEM = tau*(1.0-EM);
        double TauDOverEMMinusOne = tau*D/(EM - 1.0);
        double fix1 = tau*(EPH - EMH);
        if (fix1 == 0.0)
            fix1 = stepSize;
        double oneOverFix1 = 1.0/fix1;
        double V = sqrt(kT*(1.0 - EM));
        double X = tau*sqrt(kT*C);
        double Yv = sqrt(kT*B/C);
        double Yx = tau*sqrt(kT*B/(1.0 - EM));
        vector<cl_float> p(params->getSize());
        p[0] = EM;
        p[1] = EM;
        p[2] = DOverTauC;
        p[3] = TauOneMinusEM;
        p[4] = TauDOverEMMinusOne;
        p[5] = V;
        p[6] = X;
        p[7] = Yv;
        p[8] = Yx;
        p[9] = fix1;
        p[10] = oneOverFix1;
        params->upload(p);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    kernel1.setArg<cl_int>(0, numAtoms);
    kernel1.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(2, cl.getForce().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(3, integration.getPosDelta().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(4, params->getDeviceBuffer());
    kernel1.setArg(5, params->getSize()*sizeof(cl_float), NULL);
    kernel1.setArg<cl::Buffer>(6, xVector->getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(7, vVector->getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(8,integration.getRandom().getDeviceBuffer());
    kernel1.setArg<cl_uint>(9, integration.prepareRandomNumbers(2*numThreads));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2.setArg<cl_int>(0, numAtoms);
    kernel2.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(3, params->getDeviceBuffer());
    kernel2.setArg(4, params->getSize()*sizeof(cl_float), NULL);
    kernel2.setArg<cl::Buffer>(5, xVector->getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(6, vVector->getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(7,integration.getRandom().getDeviceBuffer());
    kernel2.setArg<cl_uint>(8, integration.prepareRandomNumbers(2*numThreads));
    cl.executeKernel(kernel2, numAtoms);

    // Reapply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the third integration kernel.

    kernel3.setArg<cl_int>(0, numAtoms);
    kernel3.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    kernel3.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
    cl.executeKernel(kernel3, numAtoms);

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
}
//
//OpenCLIntegrateBrownianStepKernel::~OpenCLIntegrateBrownianStepKernel() {
//}
//
//void OpenCLIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
//    initializeIntegration(system, data, integrator);
//    _gpuContext* gpu = data.gpu;
//    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
//    gpuInitializeRandoms(gpu);
//    prevStepSize = -1.0;
//}
//
//void OpenCLIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
//    _gpuContext* gpu = data.gpu;
//    double temperature = integrator.getTemperature();
//    double friction = integrator.getFriction();
//    double stepSize = integrator.getStepSize();
//    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
//        // Initialize the GPU parameters.
//
//        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
//        gpuSetBrownianIntegrationParameters(gpu, (float) tau, (float) stepSize, (float) temperature);
//        gpuSetConstants(gpu);
//        kGenerateRandoms(gpu);
//        prevTemp = temperature;
//        prevFriction = friction;
//        prevStepSize = stepSize;
//    }
//    kBrownianUpdatePart1(gpu);
//    kApplyFirstShake(gpu);
//    kApplyFirstSettle(gpu);
//    kApplyFirstCCMA(gpu);
//    if (data.removeCM)
//        if (data.stepCount%data.cmMotionFrequency == 0)
//            gpu->bCalculateCM = true;
//    kBrownianUpdatePart2(gpu);
//    data.time += stepSize;
//    data.stepCount++;
//}
//
//OpenCLIntegrateVariableVerletStepKernel::~OpenCLIntegrateVariableVerletStepKernel() {
//}
//
//void OpenCLIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
//    initializeIntegration(system, data, integrator);
//    prevErrorTol = -1.0;
//}
//
//void OpenCLIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
//    _gpuContext* gpu = data.gpu;
//    double errorTol = integrator.getErrorTolerance();
//    if (errorTol != prevErrorTol) {
//        // Initialize the GPU parameters.
//
//        gpuSetVerletIntegrationParameters(gpu, 0.0f, (float) errorTol);
//        gpuSetConstants(gpu);
//        prevErrorTol = errorTol;
//    }
//    float maxStepSize = (float)(maxTime-data.time);
//    kSelectVerletStepSize(gpu, maxStepSize);
//    kVerletUpdatePart1(gpu);
//    kApplyFirstShake(gpu);
//    kApplyFirstSettle(gpu);
//    kApplyFirstCCMA(gpu);
//    if (data.removeCM)
//        if (data.stepCount%data.cmMotionFrequency == 0)
//            gpu->bCalculateCM = true;
//    kVerletUpdatePart2(gpu);
//    gpu->psStepSize->Download();
//    data.time += (*gpu->psStepSize)[0].y;
//    if ((*gpu->psStepSize)[0].y == maxStepSize)
//        data.time = maxTime; // Avoid round-off error
//    data.stepCount++;
//}
//
//OpenCLIntegrateVariableLangevinStepKernel::~OpenCLIntegrateVariableLangevinStepKernel() {
//}
//
//void OpenCLIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
//    initializeIntegration(system, data, integrator);
//    _gpuContext* gpu = data.gpu;
//    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
//    gpuInitializeRandoms(gpu);
//    prevErrorTol = -1.0;
//}
//
//void OpenCLIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
//    _gpuContext* gpu = data.gpu;
//    double temperature = integrator.getTemperature();
//    double friction = integrator.getFriction();
//    double errorTol = integrator.getErrorTolerance();
//    if (temperature != prevTemp || friction != prevFriction || errorTol != prevErrorTol) {
//        // Initialize the GPU parameters.
//
//        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
//        gpuSetLangevinIntegrationParameters(gpu, (float) tau, 0.0f, (float) temperature, errorTol);
//        gpuSetConstants(gpu);
//        kGenerateRandoms(gpu);
//        prevTemp = temperature;
//        prevFriction = friction;
//        prevErrorTol = errorTol;
//    }
//    float maxStepSize = (float)(maxTime-data.time);
//    kSelectLangevinStepSize(gpu, maxStepSize);
//    kLangevinUpdatePart1(gpu);
//    kApplyFirstShake(gpu);
//    kApplyFirstSettle(gpu);
//    kApplyFirstCCMA(gpu);
//    if (data.removeCM)
//        if (data.stepCount%data.cmMotionFrequency == 0)
//            gpu->bCalculateCM = true;
//    kLangevinUpdatePart2(gpu);
//    kApplySecondShake(gpu);
//    kApplySecondSettle(gpu);
//    kApplySecondCCMA(gpu);
//    gpu->psStepSize->Download();
//    data.time += (*gpu->psStepSize)[0].y;
//    if ((*gpu->psStepSize)[0].y == maxStepSize)
//        data.time = maxTime; // Avoid round-off error
//    data.stepCount++;
//}
//
//OpenCLApplyAndersenThermostatKernel::~OpenCLApplyAndersenThermostatKernel() {
//}
//
//void OpenCLApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
//    _gpuContext* gpu = data.gpu;
//    gpu->seed = (unsigned long) thermostat.getRandomNumberSeed();
//    gpuInitializeRandoms(gpu);
//    prevStepSize = -1.0;
//}
//
//void OpenCLApplyAndersenThermostatKernel::execute(ContextImpl& context) {
//    _gpuContext* gpu = data.gpu;
//    double temperature = context.getParameter(AndersenThermostat::Temperature());
//    double frequency = context.getParameter(AndersenThermostat::CollisionFrequency());
//    double stepSize = context.getIntegrator().getStepSize();
//    if (temperature != prevTemp || frequency != prevFrequency || stepSize != prevStepSize) {
//        // Initialize the GPU parameters.
//
//        gpuSetAndersenThermostatParameters(gpu, (float) temperature, frequency);
//        gpuSetConstants(gpu);
//        kGenerateRandoms(gpu);
//        prevTemp = temperature;
//        prevFrequency = frequency;
//        prevStepSize = stepSize;
//    }
//    kCalculateAndersenThermostat(gpu);
//}
//
void OpenCLCalcKineticEnergyKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double OpenCLCalcKineticEnergyKernel::execute(ContextImpl& context) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.

    OpenCLArray<mm_float4>& velm = cl.getVelm();
    velm.download();
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i) {
        mm_float4 v = velm[i];
        energy += masses[i]*(v.x*v.x+v.y*v.y+v.z*v.z);
    }
    return 0.5*energy;
}
//
//void OpenCLRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
//    data.removeCM = true;
//    data.cmMotionFrequency = force.getFrequency();
//}
//
//void OpenCLRemoveCMMotionKernel::execute(ContextImpl& context) {
//}
