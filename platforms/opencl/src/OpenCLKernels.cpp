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
#include "openmm/internal/NonbondedForceImpl.h"
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

static string doubleToString(double value) {
    stringstream s;
    s.precision(8);
    s << scientific << value << "f";
    return s.str();
}

static string intToString(int value) {
    stringstream s;
    s << value;
    return s.str();
}

void OpenCLCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void OpenCLCalcForcesAndEnergyKernel::beginForceComputation(ContextImpl& context) {
    if (cl.getNonbondedUtilities().getUseCutoff() && cl.getComputeForceCount()%100 == 0)
        cl.reorderAtoms();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    cl.clearBuffer(cl.getForceBuffers());
    cl.getNonbondedUtilities().prepareInteractions();
}

void OpenCLCalcForcesAndEnergyKernel::finishForceComputation(ContextImpl& context) {
    cl.getNonbondedUtilities().computeInteractions();
    cl.reduceBuffer(cl.getForceBuffers(), cl.getNumForceBuffers());
}

void OpenCLCalcForcesAndEnergyKernel::beginEnergyComputation(ContextImpl& context) {
    if (cl.getNonbondedUtilities().getUseCutoff() && cl.getComputeForceCount()%100 == 0)
        cl.reorderAtoms();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    cl.clearBuffer(cl.getEnergyBuffer());
    cl.getNonbondedUtilities().prepareInteractions();
}

double OpenCLCalcForcesAndEnergyKernel::finishEnergyComputation(ContextImpl& context) {
    cl.getNonbondedUtilities().computeInteractions();
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
    mm_float4 periodicBoxSize = cl.getNonbondedUtilities().getPeriodicBoxSize();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 pos = posq[i];
        mm_int4 offset = cl.getPosCellOffsets()[i];
        positions[order[i]] = Vec3(pos.x-offset.x*periodicBoxSize.x, pos.y-offset.y*periodicBoxSize.y, pos.z-offset.z*periodicBoxSize.z);
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
    for (int i = 0; i < cl.getPosCellOffsets().size(); i++)
        cl.getPosCellOffsets()[i] = (mm_int4) {0, 0, 0, 0};
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
    OpenCLBondForceInfo(int requiredBuffers, const HarmonicBondForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    for (int i = 0; i < forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
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
    OpenCLAngleForceInfo(int requiredBuffers, const HarmonicAngleForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    for (int i = 0; i < forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
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
    OpenCLPeriodicTorsionForceInfo(int requiredBuffers, const PeriodicTorsionForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    for (int i = 0; i < forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
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
    OpenCLRBTorsionForceInfo(int requiredBuffers, const RBTorsionForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    for (int i = 0; i < forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
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

class OpenCLNonbondedForceInfo : public OpenCLForceInfo {
public:
    OpenCLNonbondedForceInfo(int requiredBuffers, const NonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, sigma1, sigma2, epsilon1, epsilon2;
        force.getParticleParameters(particle1, charge1, sigma1, epsilon1);
        force.getParticleParameters(particle2, charge2, sigma2, epsilon2);
        return (charge1 == charge2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double chargeProd1, chargeProd2, sigma1, sigma2, epsilon1, epsilon2;
        force.getExceptionParameters(group1, particle1, particle2, chargeProd1, sigma1, epsilon1);
        force.getExceptionParameters(group2, particle1, particle2, chargeProd2, sigma2, epsilon2);
        return (chargeProd1 == chargeProd2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
private:
    const NonbondedForce& force;
};

OpenCLCalcNonbondedForceKernel::~OpenCLCalcNonbondedForceKernel() {
    if (sigmaEpsilon != NULL)
        delete sigmaEpsilon;
    if (exceptionParams != NULL)
        delete exceptionParams;
    if (exceptionIndices != NULL)
        delete exceptionIndices;
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
    double sumSquaredCharges = 0.0;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        posq[i].w = (float) charge;
        sigmaEpsilonVector[i] = (mm_float2) {(float) (0.5*sigma), (float) (2.0*sqrt(epsilon))};
        exclusionList[i].push_back(i);
        sumSquaredCharges += charge*charge;
    }
    for (int i = 0; i < (int) exclusions.size(); i++) {
        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
    }
    posq.upload();
    sigmaEpsilon->upload(sigmaEpsilonVector);
    bool useCutoff = (force.getNonbondedMethod() != NonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != NonbondedForce::NoCutoff && force.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic);
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    map<string, string> defines;
    if (useCutoff) {
        // Compute the reaction field constants.

        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
        defines["REACTION_FIELD_K"] = doubleToString(reactionFieldK);
        defines["REACTION_FIELD_C"] = doubleToString(reactionFieldC);
    }
    if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
        // Compute the Ewald parameters.

        double alpha;
        int kmaxx, kmaxy, kmaxz;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
        defines["EWALD_ALPHA"] = doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        double selfEnergyScale = 138.935485*alpha/std::sqrt(M_PI);
        ewaldSelfEnergy = - 138.935485*alpha*sumSquaredCharges/std::sqrt(M_PI);

        // Create the reciprocal space kernels.

        map<string, string> replacements;
        replacements["NUM_ATOMS"] = intToString(numParticles);
        replacements["KMAX_X"] = intToString(kmaxx);
        replacements["KMAX_Y"] = intToString(kmaxy);
        replacements["KMAX_Z"] = intToString(kmaxz);
        replacements["RECIPROCAL_BOX_SIZE_X"] = doubleToString(2.0*M_PI/boxVectors[0][0]);
        replacements["RECIPROCAL_BOX_SIZE_Y"] = doubleToString(2.0*M_PI/boxVectors[1][1]);
        replacements["RECIPROCAL_BOX_SIZE_Z"] = doubleToString(2.0*M_PI/boxVectors[2][2]);
        replacements["RECIPROCAL_COEFFICIENT"] = doubleToString(138.935485*4*M_PI/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]));
        replacements["EXP_COEFFICIENT"] = doubleToString(-1.0/(4.0*alpha*alpha));
        cl::Program program = cl.createProgram(cl.loadSourceFromFile("ewald.cl"), replacements);
        ewaldSumsKernel = cl::Kernel(program, "calculateEwaldCosSinSums");
        ewaldForcesKernel = cl::Kernel(program, "calculateEwaldForces");
        cosSinSums = new OpenCLArray<mm_float2>(cl, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), "cosSinSums");
    }
    else
        ewaldSelfEnergy = 0.0;

    // Add the interaction to the default nonbonded kernel.
    
    string source = cl.loadSourceFromFile("coulombLennardJones.cl", defines);
    cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source);
    cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo("sigmaEpsilon", "float2", sizeof(cl_float2), sigmaEpsilon->getDeviceBuffer()));
    cutoffSquared = force.getCutoffDistance()*force.getCutoffDistance();

    // Initialize the exceptions.

    int numExceptions = exceptions.size();
    int maxBuffers = cl.getNonbondedUtilities().getNumForceBuffers();
    if (numExceptions > 0) {
        exceptionParams = new OpenCLArray<mm_float4>(cl, numExceptions, "exceptionParams");
        exceptionIndices = new OpenCLArray<mm_int4>(cl, numExceptions, "exceptionIndices");
        vector<mm_float4> exceptionParamsVector(numExceptions);
        vector<mm_int4> exceptionIndicesVector(numExceptions);
        vector<int> forceBufferCounter(system.getNumParticles(), 0);
        for (int i = 0; i < numExceptions; i++) {
            int particle1, particle2;
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[i], particle1, particle2, chargeProd, sigma, epsilon);
            exceptionParamsVector[i] = (mm_float4) {(float) (138.935485*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f};
            exceptionIndicesVector[i] = (mm_int4) {particle1, particle2, forceBufferCounter[particle1]++, forceBufferCounter[particle2]++};
        }
        exceptionParams->upload(exceptionParamsVector);
        exceptionIndices->upload(exceptionIndicesVector);
        for (int i = 0; i < forceBufferCounter.size(); i++)
            maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLNonbondedForceInfo(maxBuffers, force));
    if (useCutoff) {
        defines["USE_CUTOFF"] = "1";
    }
    if (usePeriodic)
        defines["USE_PERIODIC"] = "1";
    cl::Program program = cl.createProgram(cl.loadSourceFromFile("nonbondedExceptions.cl"), defines);
    exceptionsKernel = cl::Kernel(program, "computeNonbondedExceptions");
}

void OpenCLCalcNonbondedForceKernel::executeForces(ContextImpl& context) {
    if (exceptionIndices != NULL) {
        int numExceptions = exceptionIndices->getSize();
        exceptionsKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        exceptionsKernel.setArg<cl_int>(1, numExceptions);
        exceptionsKernel.setArg<cl_float>(2, cutoffSquared);
        exceptionsKernel.setArg<mm_float4>(3, cl.getNonbondedUtilities().getPeriodicBoxSize());
        exceptionsKernel.setArg<cl::Buffer>(4, cl.getForceBuffers().getDeviceBuffer());
        exceptionsKernel.setArg<cl::Buffer>(5, cl.getEnergyBuffer().getDeviceBuffer());
        exceptionsKernel.setArg<cl::Buffer>(6, cl.getPosq().getDeviceBuffer());
        exceptionsKernel.setArg<cl::Buffer>(7, exceptionParams->getDeviceBuffer());
        exceptionsKernel.setArg<cl::Buffer>(8, exceptionIndices->getDeviceBuffer());
        cl.executeKernel(exceptionsKernel, numExceptions);
    }
    if (cosSinSums != NULL) {
        ewaldSumsKernel.setArg<cl::Buffer>(0, cl.getEnergyBuffer().getDeviceBuffer());
        ewaldSumsKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
        ewaldSumsKernel.setArg<cl::Buffer>(2, cosSinSums->getDeviceBuffer());
        cl.executeKernel(ewaldSumsKernel, cosSinSums->getSize());
        ewaldForcesKernel.setArg<cl::Buffer>(0, cl.getForceBuffers().getDeviceBuffer());
        ewaldForcesKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
        ewaldForcesKernel.setArg<cl::Buffer>(2, cosSinSums->getDeviceBuffer());
        cl.executeKernel(ewaldForcesKernel, cl.getNumAtoms());
    }
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

class OpenCLGBSAOBCForceInfo : public OpenCLForceInfo {
public:
    OpenCLGBSAOBCForceInfo(int requiredBuffers, const GBSAOBCForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, scale1, scale2;
        force.getParticleParameters(particle1, charge1, radius1, scale1);
        force.getParticleParameters(particle2, charge2, radius2, scale2);
        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
    }
private:
    const GBSAOBCForce& force;
};

OpenCLCalcGBSAOBCForceKernel::~OpenCLCalcGBSAOBCForceKernel() {
    if (params != NULL)
        delete params;
    if (bornSum != NULL)
        delete bornSum;
    if (bornRadii != NULL)
        delete bornRadii;
    if (bornForce != NULL)
        delete bornForce;
    if (obcChain != NULL)
        delete obcChain;
}

void OpenCLCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    params = new OpenCLArray<mm_float2>(cl, cl.getPaddedNumAtoms(), "gbsaObcParams");
    bornRadii = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms(), "bornRadii");
    obcChain = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms(), "obcChain");
    bornSum = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornSum");
    bornForce = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornForce");
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    int numParticles = force.getNumParticles();
    vector<mm_float2> paramsVector(numParticles);
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        paramsVector[i] = (mm_float2) {(float) radius, (float) (scalingFactor*radius)};
        posq[i].w = (float) charge;
    }
    posq.upload();
    params->upload(paramsVector);
    prefactor = 2.0*-166.02691*0.4184*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    string source = cl.loadSourceFromFile("gbsaObc2.cl");
    nb.addInteraction(useCutoff, usePeriodic, false, force.getCutoffDistance(), vector<vector<int> >(), source);
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo("obcParams", "float2", sizeof(cl_float2), params->getDeviceBuffer()));;
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo("bornForce", "float", sizeof(cl_float), bornForce->getDeviceBuffer()));;
    cl.addForce(new OpenCLGBSAOBCForceInfo(nb.getNumForceBuffers(), force));
}

void OpenCLCalcGBSAOBCForceKernel::executeForces(ContextImpl& context) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    if (!hasCreatedKernels) {
        // These Kernels cannot be created in initialize(), because the OpenCLNonbondedUtilities has not been initialized yet then.

        hasCreatedKernels = true;
        map<string, string> defines;
        if (nb.getForceBufferPerAtomBlock())
            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
        if (nb.getUseCutoff())
            defines["USE_CUTOFF"] = "1";
        if (nb.getUsePeriodic())
            defines["USE_PERIODIC"] = "1";
        defines["PERIODIC_BOX_SIZE_X"] = doubleToString(nb.getPeriodicBoxSize().x);
        defines["PERIODIC_BOX_SIZE_Y"] = doubleToString(nb.getPeriodicBoxSize().y);
        defines["PERIODIC_BOX_SIZE_Z"] = doubleToString(nb.getPeriodicBoxSize().z);
        defines["CUTOFF_SQUARED"] = doubleToString(nb.getCutoffDistance()*nb.getCutoffDistance());
        defines["PREFACTOR"] = doubleToString(prefactor);
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
        string filename = (cl.getSIMDWidth() == 32 ? "gbsaObc_nvidia.cl" : "gbsaObc_default.cl");
        cl::Program program = cl.createProgram(cl.loadSourceFromFile(filename), defines);
        computeBornSumKernel = cl::Kernel(program, "computeBornSum");
        computeBornSumKernel.setArg<cl::Buffer>(0, bornSum->getDeviceBuffer());
        computeBornSumKernel.setArg(1, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        computeBornSumKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        computeBornSumKernel.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        computeBornSumKernel.setArg<cl::Buffer>(4, params->getDeviceBuffer());
        computeBornSumKernel.setArg(5, OpenCLContext::ThreadBlockSize*sizeof(cl_float2), NULL);
        computeBornSumKernel.setArg(6, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        if (nb.getUseCutoff()) {
            computeBornSumKernel.setArg<cl::Buffer>(7, nb.getInteractingTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(8, nb.getInteractionFlags().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(9, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            computeBornSumKernel.setArg<cl::Buffer>(7, nb.getTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl_uint>(8, nb.getTiles().getSize());
        }
        force1Kernel = cl::Kernel(program, "computeGBSAForce1");
        force1Kernel.setArg<cl::Buffer>(0, cl.getForceBuffers().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        force1Kernel.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        force1Kernel.setArg(4, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        force1Kernel.setArg<cl::Buffer>(5, bornRadii->getDeviceBuffer());
        force1Kernel.setArg(6, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        force1Kernel.setArg<cl::Buffer>(7, bornForce->getDeviceBuffer());
        force1Kernel.setArg(8, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        force1Kernel.setArg(9, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        if (nb.getUseCutoff()) {
            force1Kernel.setArg<cl::Buffer>(10, nb.getInteractingTiles().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(11, nb.getInteractionFlags().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(12, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            force1Kernel.setArg<cl::Buffer>(10, nb.getTiles().getDeviceBuffer());
            force1Kernel.setArg<cl_uint>(11, nb.getTiles().getSize());
        }
        program = cl.createProgram(cl.loadSourceFromFile("gbsaObcReductions.cl"), defines);
        reduceBornSumKernel = cl::Kernel(program, "reduceBornSum");
        reduceBornSumKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        reduceBornSumKernel.setArg<cl_int>(1, cl.getNumForceBuffers());
        reduceBornSumKernel.setArg<cl_float>(2, 1.0f);
        reduceBornSumKernel.setArg<cl_float>(3, 0.8f);
        reduceBornSumKernel.setArg<cl_float>(4, 4.85f);
        reduceBornSumKernel.setArg<cl::Buffer>(5, bornSum->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(6, params->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(7, bornRadii->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(8, obcChain->getDeviceBuffer());
        reduceBornForceKernel = cl::Kernel(program, "reduceBornForce");
        reduceBornForceKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        reduceBornForceKernel.setArg<cl_int>(1, cl.getNumForceBuffers());
        reduceBornForceKernel.setArg<cl::Buffer>(2, bornForce->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(4, params->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(5, bornRadii->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(6, obcChain->getDeviceBuffer());
    }
    cl.clearBuffer(*bornSum);
    cl.clearBuffer(*bornForce);
    cl.executeKernel(computeBornSumKernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
    cl.executeKernel(reduceBornSumKernel, cl.getPaddedNumAtoms());
    cl.executeKernel(force1Kernel, cl.getPaddedNumAtoms());
    cl.executeKernel(reduceBornForceKernel, cl.getPaddedNumAtoms());
}

double OpenCLCalcGBSAOBCForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

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
