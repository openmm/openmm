/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "AmoebaCudaKernels.h"
#include "CudaAmoebaKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/AmoebaTorsionTorsionForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaForceInfo.h"
#include "CudaKernelSources.h"
#include "CudaNonbondedUtilities.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<errorMessage<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

/* -------------------------------------------------------------------------- *
 *                           AmoebaHarmonicBond                               *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaHarmonicBondForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaHarmonicBondForce& force) : force(force) {
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
    const AmoebaHarmonicBondForce& force;
};

CudaCalcAmoebaHarmonicBondForceKernel::CudaCalcAmoebaHarmonicBondForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
                CalcAmoebaHarmonicBondForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcAmoebaHarmonicBondForceKernel::~CudaCalcAmoebaHarmonicBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaHarmonicBondForceKernel::initialize(const System& system, const AmoebaHarmonicBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = CudaArray::create<float2>(cu, numBonds, "bondParams");
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaAmoebaKernelSources::amoebaBondForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    replacements["CUBIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicBondCubic());
    replacements["QUARTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicBondQuartic());
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}
//
///* -------------------------------------------------------------------------- *
// *                           AmoebaUreyBradley                               *
// * -------------------------------------------------------------------------- */
//
//class CudaCalcAmoebaUreyBradleyForceKernel::ForceInfo : public CudaForceInfo {
//public:
//    ForceInfo(const AmoebaUreyBradleyForce& force) : force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumInteractions();
//    }
//    void getParticlesInGroup(int index, std::vector<int>& particles) {
//        int particle1, particle2;
//        double length, k;
//        force.getUreyBradleyParameters(index, particle1, particle2, length, k);
//        particles.resize(2);
//        particles[0] = particle1;
//        particles[1] = particle2;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2;
//        double length1, length2, k1, k2;
//        force.getUreyBradleyParameters(group1, particle1, particle2, length1, k1);
//        force.getUreyBradleyParameters(group2, particle1, particle2, length2, k2);
//        return (length1 == length2 && k1 == k2);
//    }
//private:
//    const AmoebaUreyBradleyForce& force;
//};
//
//CudaCalcAmoebaUreyBradleyForceKernel::CudaCalcAmoebaUreyBradleyForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
//                CalcAmoebaUreyBradleyForceKernel(name, platform), cu(cu), system(system) {
//    data.incrementKernelCount();
//}
//
//CudaCalcAmoebaUreyBradleyForceKernel::~CudaCalcAmoebaUreyBradleyForceKernel() {
//    data.decrementKernelCount();
//}
//
//void CudaCalcAmoebaUreyBradleyForceKernel::initialize(const System& system, const AmoebaUreyBradleyForce& force) {
//
//    data.setAmoebaLocalForcesKernel( this );
//
//    numInteractions = force.getNumInteractions();
//    std::vector<int>   particle1(numInteractions);
//    std::vector<int>   particle2(numInteractions);
//    std::vector<float> length(numInteractions);
//    std::vector<float> quadratic(numInteractions);
//
//    for (int i = 0; i < numInteractions; i++) {
//
//        int particle1Index, particle2Index;
//        double lengthValue, kValue;
//        force.getUreyBradleyParameters(i, particle1Index, particle2Index, lengthValue, kValue );
//
//        particle1[i]     = particle1Index; 
//        particle2[i]     = particle2Index; 
//        length[i]        = static_cast<float>( lengthValue );
//        quadratic[i]     = static_cast<float>( kValue );
//    } 
//    gpuSetAmoebaUreyBradleyParameters( data.getAmoebaGpu(), particle1, particle2, length, quadratic, 
//                                       static_cast<float>(force.getAmoebaGlobalUreyBradleyCubic()),
//                                       static_cast<float>(force.getAmoebaGlobalUreyBradleyQuartic()) );
//    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
//}
//
//double CudaCalcAmoebaUreyBradleyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if( data.getAmoebaLocalForcesKernel() == this ){
//        computeAmoebaLocalForces( data );
//    }
//    return 0.0;
//}

/* -------------------------------------------------------------------------- *
 *                           AmoebaHarmonicAngle                              *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaHarmonicAngleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaHarmonicAngleForce& force) : force(force) {
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
    const AmoebaHarmonicAngleForce& force;
};

CudaCalcAmoebaHarmonicAngleForceKernel::CudaCalcAmoebaHarmonicAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
            CalcAmoebaHarmonicAngleForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcAmoebaHarmonicAngleForceKernel::~CudaCalcAmoebaHarmonicAngleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaHarmonicAngleForceKernel::initialize(const System& system, const AmoebaHarmonicAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = CudaArray::create<float2>(cu, numAngles, "angleParams");
    vector<float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = make_float2((float) angle, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaAmoebaKernelSources::amoebaAngleForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    replacements["CUBIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicAngleCubic());
    replacements["QUARTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicAngleQuartic());
    replacements["PENTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicAnglePentic());
    replacements["SEXTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicAngleSextic());
    replacements["RAD_TO_DEG"] = cu.doubleToString(180/M_PI);
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaHarmonicInPlaneAngle                       *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaHarmonicInPlaneAngleForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        double angle, k;
        force.getAngleParameters(index, particle1, particle2, particle3, particle4, angle, k);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        double angle1, angle2, k1, k2;
        force.getAngleParameters(group1, particle1, particle2, particle3, particle4, angle1, k1);
        force.getAngleParameters(group2, particle1, particle2, particle3, particle4, angle2, k2);
        return (angle1 == angle2 && k1 == k2);
    }
private:
    const AmoebaHarmonicInPlaneAngleForce& force;
};

CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::CudaCalcAmoebaHarmonicInPlaneAngleForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
          CalcAmoebaHarmonicInPlaneAngleForceKernel(name, platform), cu(cu), system(system) {
}

CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::~CudaCalcAmoebaHarmonicInPlaneAngleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::initialize(const System& system, const AmoebaHarmonicInPlaneAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(4));
    params = CudaArray::create<float2>(cu, numAngles, "angleParams");
    vector<float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], angle, k);
        paramVector[i] = make_float2((float) angle, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    replacements["CUBIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicInPlaneAngleCubic());
    replacements["QUARTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicInPlaneAngleQuartic());
    replacements["PENTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicInPlaneAnglePentic());
    replacements["SEXTIC_K"] = cu.doubleToString(force.getAmoebaGlobalHarmonicInPlaneAngleSextic());
    replacements["RAD_TO_DEG"] = cu.doubleToString(180/M_PI);
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaAmoebaKernelSources::amoebaInPlaneForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

///* -------------------------------------------------------------------------- *
// *                               AmoebaTorsion                                *
// * -------------------------------------------------------------------------- */
//
//class CudaCalcAmoebaTorsionForceKernel::ForceInfo : public CudaForceInfo {
//public:
//    ForceInfo(const AmoebaTorsionForce& force) : force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumTorsions();
//    }
//    void getParticlesInGroup(int index, std::vector<int>& particles) {
//        int particle1, particle2, particle3, particle4;
//        vector<double> torsion1, torsion2, torsion3;
//        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3);
//        particles.resize(4);
//        particles[0] = particle1;
//        particles[1] = particle2;
//        particles[2] = particle3;
//        particles[3] = particle4;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2, particle3, particle4;
//        vector<double> torsion11, torsion21, torsion31;
//        vector<double> torsion12, torsion22, torsion32;
//        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, torsion11, torsion21, torsion31);
//        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, torsion12, torsion22, torsion32);
//        for (int i = 0; i < (int) torsion11.size(); ++i)
//            if (torsion11[i] != torsion12[i])
//                return false;
//        for (int i = 0; i < (int) torsion21.size(); ++i)
//            if (torsion21[i] != torsion22[i])
//                return false;
//        for (int i = 0; i < (int) torsion31.size(); ++i)
//            if (torsion31[i] != torsion32[i])
//                return false;
//        return true;
//    }
//private:
//    const AmoebaTorsionForce& force;
//};
//
//CudaCalcAmoebaTorsionForceKernel::CudaCalcAmoebaTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
//             CalcAmoebaTorsionForceKernel(name, platform), cu(cu), system(system) {
//    data.incrementKernelCount();
//}
//
//CudaCalcAmoebaTorsionForceKernel::~CudaCalcAmoebaTorsionForceKernel() {
//    data.decrementKernelCount();
//}
//
//void CudaCalcAmoebaTorsionForceKernel::initialize(const System& system, const AmoebaTorsionForce& force) {
//
//    data.setAmoebaLocalForcesKernel( this );
//    numTorsions                     = force.getNumTorsions();
//    std::vector<int> particle1(numTorsions);
//    std::vector<int> particle2(numTorsions);
//    std::vector<int> particle3(numTorsions);
//    std::vector<int> particle4(numTorsions);
//
//    std::vector< std::vector<float> > torsionParameters1(numTorsions);
//    std::vector< std::vector<float> > torsionParameters2(numTorsions);
//    std::vector< std::vector<float> > torsionParameters3(numTorsions);
//
//    for (int i = 0; i < numTorsions; i++) {
//
//        std::vector<double> torsionParameter1;
//        std::vector<double> torsionParameter2;
//        std::vector<double> torsionParameter3;
//
//        std::vector<float> torsionParameters1F(3);
//        std::vector<float> torsionParameters2F(3);
//        std::vector<float> torsionParameters3F(3);
//
//        force.getTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], torsionParameter1, torsionParameter2, torsionParameter3 );
//        for ( unsigned int jj = 0; jj < torsionParameter1.size(); jj++) {
//            torsionParameters1F[jj] = static_cast<float>(torsionParameter1[jj]);
//            torsionParameters2F[jj] = static_cast<float>(torsionParameter2[jj]);
//            torsionParameters3F[jj] = static_cast<float>(torsionParameter3[jj]);
//        }
//        torsionParameters1[i] = torsionParameters1F;
//        torsionParameters2[i] = torsionParameters2F;
//        torsionParameters3[i] = torsionParameters3F;
//    }
//    gpuSetAmoebaTorsionParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, torsionParameters1, torsionParameters2, torsionParameters3 );
//    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
//}
//
//double CudaCalcAmoebaTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if( data.getAmoebaLocalForcesKernel() == this ){
//        computeAmoebaLocalForces( data );
//    }
//    return 0.0;
//}

/* -------------------------------------------------------------------------- *
  *                              AmoebaPiTorsion                              *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaPiTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaPiTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumPiTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4, particle5, particle6;
        double k;
        force.getPiTorsionParameters(index, particle1, particle2, particle3, particle4, particle5, particle6, k);
        particles.resize(6);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
        particles[4] = particle5;
        particles[5] = particle6;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, particle5, particle6;
        double k1, k2;
        force.getPiTorsionParameters(group1, particle1, particle2, particle3, particle4, particle5, particle6, k1);
        force.getPiTorsionParameters(group2, particle1, particle2, particle3, particle4, particle5, particle6, k2);
        return (k1 == k2);
    }
private:
    const AmoebaPiTorsionForce& force;
};

CudaCalcAmoebaPiTorsionForceKernel::CudaCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
         CalcAmoebaPiTorsionForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcAmoebaPiTorsionForceKernel::~CudaCalcAmoebaPiTorsionForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaPiTorsionForceKernel::initialize(const System& system, const AmoebaPiTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumPiTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumPiTorsions()/numContexts;
    numPiTorsions = endIndex-startIndex;
    if (numPiTorsions == 0)
        return;
    vector<vector<int> > atoms(numPiTorsions, vector<int>(6));
    params = CudaArray::create<float>(cu, numPiTorsions, "piTorsionParams");
    vector<float> paramVector(numPiTorsions);
    for (int i = 0; i < numPiTorsions; i++) {
        double k;
        force.getPiTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], k);
        paramVector[i] = (float) k;
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaAmoebaKernelSources::amoebaPiTorsionForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaPiTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaStretchBend                                *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaStretchBendForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaStretchBendForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumStretchBends();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3;
        double lengthAB, lengthCB, angle, k;
        force.getStretchBendParameters(index, particle1, particle2, particle3, lengthAB, lengthCB, angle, k);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        double lengthAB1, lengthAB2, lengthCB1, lengthCB2, angle1, angle2, k1, k2;
        force.getStretchBendParameters(group1, particle1, particle2, particle3, lengthAB1, lengthCB1, angle1, k1);
        force.getStretchBendParameters(group2, particle1, particle2, particle3, lengthAB2, lengthCB2, angle2, k2);
        return (lengthAB1 == lengthAB2 && lengthCB1 == lengthCB2 && angle1 == angle2 && k1 == k2);
    }
private:
    const AmoebaStretchBendForce& force;
};

CudaCalcAmoebaStretchBendForceKernel::CudaCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
                   CalcAmoebaStretchBendForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcAmoebaStretchBendForceKernel::~CudaCalcAmoebaStretchBendForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaStretchBendForceKernel::initialize(const System& system, const AmoebaStretchBendForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumStretchBends()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumStretchBends()/numContexts;
    numStretchBends = endIndex-startIndex;
    if (numStretchBends == 0)
        return;
    vector<vector<int> > atoms(numStretchBends, vector<int>(3));
    params = CudaArray::create<float4>(cu, numStretchBends, "stretchBendParams");
    vector<float4> paramVector(numStretchBends);
    for (int i = 0; i < numStretchBends; i++) {
        double lengthAB, lengthCB, angle, k;
        force.getStretchBendParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], lengthAB, lengthCB, angle, k);
        paramVector[i] = make_float4((float) lengthAB, (float) lengthCB, (float) angle, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float4");
    replacements["RAD_TO_DEG"] = cu.doubleToString(180/M_PI);
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaAmoebaKernelSources::amoebaStretchBendForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaStretchBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaOutOfPlaneBend                             *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaOutOfPlaneBendForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaOutOfPlaneBendForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumOutOfPlaneBends();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        double k;
        force.getOutOfPlaneBendParameters(index, particle1, particle2, particle3, particle4, k);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        double k1, k2;
        force.getOutOfPlaneBendParameters(group1, particle1, particle2, particle3, particle4, k1);
        force.getOutOfPlaneBendParameters(group2, particle1, particle2, particle3, particle4, k2);
        return (k1 == k2);
    }
private:
    const AmoebaOutOfPlaneBendForce& force;
};

CudaCalcAmoebaOutOfPlaneBendForceKernel::CudaCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
          CalcAmoebaOutOfPlaneBendForceKernel(name, platform), cu(cu), system(system), params(NULL) {
}

CudaCalcAmoebaOutOfPlaneBendForceKernel::~CudaCalcAmoebaOutOfPlaneBendForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcAmoebaOutOfPlaneBendForceKernel::initialize(const System& system, const AmoebaOutOfPlaneBendForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumOutOfPlaneBends()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumOutOfPlaneBends()/numContexts;
    numOutOfPlaneBends = endIndex-startIndex;
    if (numOutOfPlaneBends == 0)
        return;
    vector<vector<int> > atoms(numOutOfPlaneBends, vector<int>(4));
    params = CudaArray::create<float>(cu, numOutOfPlaneBends, "outOfPlaneParams");
    vector<float> paramVector(numOutOfPlaneBends);
    for (int i = 0; i < numOutOfPlaneBends; i++) {
        double k;
        force.getOutOfPlaneBendParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], k);
        paramVector[i] = (float) k;
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float");
    replacements["CUBIC_K"] = cu.doubleToString(force.getAmoebaGlobalOutOfPlaneBendCubic());
    replacements["QUARTIC_K"] = cu.doubleToString(force.getAmoebaGlobalOutOfPlaneBendQuartic());
    replacements["PENTIC_K"] = cu.doubleToString(force.getAmoebaGlobalOutOfPlaneBendPentic());
    replacements["SEXTIC_K"] = cu.doubleToString(force.getAmoebaGlobalOutOfPlaneBendSextic());
    replacements["RAD_TO_DEG"] = cu.doubleToString(180/M_PI);
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaAmoebaKernelSources::amoebaOutOfPlaneBendForce, replacements), force.getForceGroup());
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaOutOfPlaneBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaTorsionTorsion                             *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaTorsionTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaTorsionTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsionTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex;
        force.getTorsionTorsionParameters(index, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
        particles.resize(5);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
        particles[4] = particle5;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, particle5;
        int chiral1, chiral2, grid1, grid2;
        force.getTorsionTorsionParameters(group1, particle1, particle2, particle3, particle4, particle5, chiral1, grid1);
        force.getTorsionTorsionParameters(group2, particle1, particle2, particle3, particle4, particle5, chiral2, grid2);
        return (grid1 == grid2);
    }
private:
    const AmoebaTorsionTorsionForce& force;
};

CudaCalcAmoebaTorsionTorsionForceKernel::CudaCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
                CalcAmoebaTorsionTorsionForceKernel(name, platform), cu(cu), system(system), gridValues(NULL), gridParams(NULL), torsionParams(NULL) {
}

CudaCalcAmoebaTorsionTorsionForceKernel::~CudaCalcAmoebaTorsionTorsionForceKernel() {
    cu.setAsCurrent();
    if (gridValues != NULL)
        delete gridValues;
    if (gridParams != NULL)
        delete gridParams;
    if (torsionParams != NULL)
        delete torsionParams;
}

void CudaCalcAmoebaTorsionTorsionForceKernel::initialize(const System& system, const AmoebaTorsionTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsionTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsionTorsions()/numContexts;
    numTorsionTorsions = endIndex-startIndex;
    if (numTorsionTorsions == 0)
        return;
    
    // Record torsion parameters.
    
    vector<vector<int> > atoms(numTorsionTorsions, vector<int>(5));
    vector<int2> torsionParamsVec(numTorsionTorsions);
    torsionParams = CudaArray::create<int2>(cu, numTorsionTorsions, "torsionTorsionParams");
    for (int i = 0; i < numTorsionTorsions; i++)
        force.getTorsionTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], torsionParamsVec[i].x, torsionParamsVec[i].y);
    torsionParams->upload(torsionParamsVec);
    
    // Record the grids.
    
    vector<float4> gridValuesVec;
    vector<float4> gridParamsVec;
    for (int i = 0; i < force.getNumTorsionTorsionGrids(); i++) {
        const TorsionTorsionGrid& initialGrid = force.getTorsionTorsionGrid(i);

        // check if grid needs to be reordered: x-angle should be 'slow' index

        bool reordered = false;
        TorsionTorsionGrid reorderedGrid;
        if (initialGrid[0][0][0] != initialGrid[0][1][0]) {
            AmoebaTorsionTorsionForceImpl::reorderGrid(initialGrid, reorderedGrid);
            reordered = true;
        }
        const TorsionTorsionGrid& grid = (reordered ? reorderedGrid : initialGrid);
        float range = grid[0][grid[0].size()-1][1] - grid[0][0][1];
        gridParamsVec.push_back(make_float4(gridValuesVec.size(), grid[0][0][0], range/(grid.size()-1), grid.size()));
        for (int j = 0; j < grid.size(); j++)
            for (int k = 0; k < grid[j].size(); k++)
                gridValuesVec.push_back(make_float4((float) grid[j][k][2], (float) grid[j][k][3], (float) grid[j][k][4], (float) grid[j][k][5]));
    }
    gridValues = CudaArray::create<float4>(cu, gridValuesVec.size(), "torsionTorsionGridValues");
    gridParams = CudaArray::create<float4>(cu, gridParamsVec.size(), "torsionTorsionGridParams");
    gridValues->upload(gridValuesVec);
    gridParams->upload(gridParamsVec);
    map<string, string> replacements;
    replacements["GRID_VALUES"] = cu.getBondedUtilities().addArgument(gridValues->getDevicePointer(), "float4");
    replacements["GRID_PARAMS"] = cu.getBondedUtilities().addArgument(gridParams->getDevicePointer(), "float4");
    replacements["TORSION_PARAMS"] = cu.getBondedUtilities().addArgument(torsionParams->getDevicePointer(), "int2");
    replacements["RAD_TO_DEG"] = cu.doubleToString(180/M_PI);
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaAmoebaKernelSources::amoebaTorsionTorsionForce, replacements), force.getForceGroup());
    cu.getBondedUtilities().addPrefixCode(CudaAmoebaKernelSources::bicubic);
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

///* -------------------------------------------------------------------------- *
// *                             AmoebaMultipole                                *
// * -------------------------------------------------------------------------- */
//
//static void computeAmoebaMultipoleForce( CudaContext& cu ) {
//
//    amoebaGpuContext gpu = data.getAmoebaGpu();
//    data.incrementMultipoleForceCount();
//
//    if( 0 && data.getLog() ){
//        (void) fprintf( data.getLog(), "In computeAmoebaMultipoleForce hasAmoebaGeneralizedKirkwood=%d\n", data.getHasAmoebaGeneralizedKirkwood() );
//        (void) fflush( data.getLog());
//    }
//
//    data.initializeGpu();
//
//    // calculate Born radii using either the Grycuk or OBC algorithm if GK is active
//
//    if( data.getHasAmoebaGeneralizedKirkwood() ){
//        kClearBornSum( gpu->gpuContext );
//        if( data.getUseGrycuk() ){
//            kCalculateAmoebaGrycukBornRadii( gpu );
//            kReduceGrycukGbsaBornSum( gpu );
//        } else {
//            kCalculateObcGbsaBornSum(gpu->gpuContext);
//            kReduceObcGbsaBornSum(gpu->gpuContext);
//       }
//    }   
//
//    // multipoles
//
//    kCalculateAmoebaMultipoleForces(gpu, data.getHasAmoebaGeneralizedKirkwood() );
//
//    // GK
//
//    if( data.getHasAmoebaGeneralizedKirkwood() ){
//        kCalculateAmoebaKirkwood(gpu);
//    }
//
//    if( 0 && data.getLog() ){
//        (void) fprintf( data.getLog(), "completed computeAmoebaMultipoleForce\n" );
//        (void) fflush( data.getLog());
//    }
//}
//
//static void computeAmoebaMultipolePotential( CudaContext& cu, const std::vector< Vec3 >& inputGrid,
//                                             std::vector< double >& outputElectrostaticPotential) {
//
//    amoebaGpuContext gpu = data.getAmoebaGpu();
//
//    // load grid to board and allocate board memory for potential buffers
//    // calculate potential
//    // load potential into return vector
//    // deallocate board memory
//
//    gpuSetupElectrostaticPotentialCalculation( gpu, inputGrid );
//    data.setGpuInitialized( false );
//    data.initializeGpu();
// 
//    kCalculateAmoebaMultipolePotential( gpu );
//    gpuLoadElectrostaticPotential( gpu, inputGrid.size(), outputElectrostaticPotential );
//    gpuCleanupElectrostaticPotentialCalculation( gpu );
//
//    if( 0 && data.getLog() ){
//        (void) fprintf( data.getLog(), "completed computeAmoebaMultipolePotential\n" );
//        (void) fflush( data.getLog());
//    }
//}
//
//static void computeAmoebaSystemMultipoleMoments( CudaContext& cu, const Vec3& origin,
//                                                 std::vector< double >& outputMultipoleMonents) {
//
//    amoebaGpuContext gpu = data.getAmoebaGpu();
//
//    data.setGpuInitialized( false );
//    data.initializeGpu();
//    kCalculateAmoebaSystemMultipoleMoments( gpu, origin, outputMultipoleMonents );
//
//}
//
//class CudaCalcAmoebaMultipoleForceKernel::ForceInfo : public CudaForceInfo {
//public:
//    ForceInfo(const AmoebaMultipoleForce& force) : force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        double charge1, charge2, thole1, thole2, damping1, damping2, polarity1, polarity2;
//        int axis1, axis2, multipole11, multipole12, multipole21, multipole22, multipole31, multipole32;
//        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
//        force.getMultipoleParameters(particle1, charge1, dipole1, quadrupole1, axis1, multipole11, multipole21, multipole31, thole1, damping1, polarity1);
//        force.getMultipoleParameters(particle2, charge2, dipole2, quadrupole2, axis2, multipole12, multipole22, multipole32, thole2, damping2, polarity2);
//        if (charge1 != charge2 || thole1 != thole2 || damping1 != damping2 || polarity1 != polarity2 || axis1 != axis2){
//            return false;
//        }
//        for (int i = 0; i < (int) dipole1.size(); ++i){
//            if (dipole1[i] != dipole2[i]){
//                return false;
//            }
//        }
//        for (int i = 0; i < (int) quadrupole1.size(); ++i){
//            if (quadrupole1[i] != quadrupole2[i]){
//                return false;
//            }
//        }
//        return true;
//    }
//private:
//    const AmoebaMultipoleForce& force;
//};
//
//CudaCalcAmoebaMultipoleForceKernel::CudaCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
//         CalcAmoebaMultipoleForceKernel(name, platform), cu(cu), system(system) {
//    data.incrementKernelCount();
//}
//
//CudaCalcAmoebaMultipoleForceKernel::~CudaCalcAmoebaMultipoleForceKernel() {
//    data.decrementKernelCount();
//}
//
//void CudaCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {
//
//    numMultipoles   = force.getNumMultipoles();
//
//    data.setHasAmoebaMultipole( true );
//
//    std::vector<float> charges(numMultipoles);
//    std::vector<float> dipoles(3*numMultipoles);
//    std::vector<float> quadrupoles(9*numMultipoles);
//    std::vector<float> tholes(numMultipoles);
//    std::vector<float> dampingFactors(numMultipoles);
//    std::vector<float> polarity(numMultipoles);
//    std::vector<int>   axisTypes(numMultipoles);
//    std::vector<int>   multipoleAtomZs(numMultipoles);
//    std::vector<int>   multipoleAtomXs(numMultipoles);
//    std::vector<int>   multipoleAtomYs(numMultipoles);
//    std::vector< std::vector< std::vector<int> > > multipoleAtomCovalentInfo(numMultipoles);
//    std::vector<int> minCovalentIndices(numMultipoles);
//    std::vector<int> minCovalentPolarizationIndices(numMultipoles);
//
//    //float scalingDistanceCutoff = static_cast<float>(force.getScalingDistanceCutoff());
//    float scalingDistanceCutoff = 50.0f;
//
//    std::vector<AmoebaMultipoleForce::CovalentType> covalentList;
//    covalentList.push_back( AmoebaMultipoleForce::Covalent12 );
//    covalentList.push_back( AmoebaMultipoleForce::Covalent13 );
//    covalentList.push_back( AmoebaMultipoleForce::Covalent14 );
//    covalentList.push_back( AmoebaMultipoleForce::Covalent15 );
//
//    std::vector<AmoebaMultipoleForce::CovalentType> polarizationCovalentList;
//    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent11 );
//    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent12 );
//    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent13 );
//    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent14 );
//
//    std::vector<int> covalentDegree;
//    AmoebaMultipoleForceImpl::getCovalentDegree( force, covalentDegree );
//    int dipoleIndex      = 0;
//    int quadrupoleIndex  = 0;
//    int maxCovalentRange = 0;
//    double totalCharge   = 0.0;
//    for (int i = 0; i < numMultipoles; i++) {
//
//        // multipoles
//
//        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
//        double charge, tholeD, dampingFactorD, polarityD;
//        std::vector<double> dipolesD;
//        std::vector<double> quadrupolesD;
//        force.getMultipoleParameters(i, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
//                                     tholeD, dampingFactorD, polarityD );
//
//        totalCharge                       += charge;
//        axisTypes[i]                       = axisType;
//        multipoleAtomZs[i]                 = multipoleAtomZ;
//        multipoleAtomXs[i]                 = multipoleAtomX;
//        multipoleAtomYs[i]                 = multipoleAtomY;
//
//        charges[i]                         = static_cast<float>(charge);
//        tholes[i]                          = static_cast<float>(tholeD);
//        dampingFactors[i]                  = static_cast<float>(dampingFactorD);
//        polarity[i]                        = static_cast<float>(polarityD);
//
//        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[0]);
//        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[1]);
//        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[2]);
//        
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[0]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[1]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[2]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[3]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[4]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[5]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[6]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[7]);
//        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[8]);
//
//        // covalent info
//
//        std::vector< std::vector<int> > covalentLists;
//        force.getCovalentMaps(i, covalentLists );
//        multipoleAtomCovalentInfo[i] = covalentLists;
//
//        int minCovalentIndex, maxCovalentIndex;
//        AmoebaMultipoleForceImpl::getCovalentRange( force, i, covalentList, &minCovalentIndex, &maxCovalentIndex );
//        minCovalentIndices[i] = minCovalentIndex;
//        if( maxCovalentRange < (maxCovalentIndex - minCovalentIndex) ){
//            maxCovalentRange = maxCovalentIndex - minCovalentIndex;
//        }
//
//        AmoebaMultipoleForceImpl::getCovalentRange( force, i, polarizationCovalentList, &minCovalentIndex, &maxCovalentIndex );
//        minCovalentPolarizationIndices[i] = minCovalentIndex;
//        if( maxCovalentRange < (maxCovalentIndex - minCovalentIndex) ){
//            maxCovalentRange = maxCovalentIndex - minCovalentIndex;
//        }
//    }
//
//    int polarizationType = static_cast<int>(force.getPolarizationType());
//    int nonbondedMethod = static_cast<int>(force.getNonbondedMethod());
//    if( nonbondedMethod != 0 && nonbondedMethod != 1 ){
//         throw OpenMMException("AmoebaMultipoleForce nonbonded method not recognized.\n");
//    }
//
//    if( polarizationType != 0 && polarizationType != 1 ){
//         throw OpenMMException("AmoebaMultipoleForce polarization type not recognized.\n");
//    }
//
//    gpuSetAmoebaMultipoleParameters(data.getAmoebaGpu(), charges, dipoles, quadrupoles, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
//                                    tholes, scalingDistanceCutoff, dampingFactors, polarity,
//                                    multipoleAtomCovalentInfo, covalentDegree, minCovalentIndices, minCovalentPolarizationIndices, (maxCovalentRange+2),
//                                    0, force.getMutualInducedMaxIterations(),
//                                    static_cast<float>( force.getMutualInducedTargetEpsilon()),
//                                    nonbondedMethod, polarizationType,
//                                    static_cast<float>( force.getCutoffDistance()),
//                                    static_cast<float>( force.getAEwald()) );
//    if (nonbondedMethod == AmoebaMultipoleForce::PME) {
//        double alpha = force.getAEwald();
//        int xsize, ysize, zsize;
//        NonbondedForce nb;
//        nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
//        nb.setCutoffDistance(force.getCutoffDistance());
//        std::vector<int> pmeGridDimension;
//        force.getPmeGridDimensions( pmeGridDimension );
//        int pmeParametersSetBasedOnEwaldErrorTolerance;
//        if( pmeGridDimension[0] == 0 || alpha == 0.0 ){
//            NonbondedForceImpl::calcPMEParameters(system, nb, alpha, xsize, ysize, zsize);
//            pmeParametersSetBasedOnEwaldErrorTolerance = 1;
//        } else {
//            alpha = force.getAEwald();
//            xsize = pmeGridDimension[0];
//            ysize = pmeGridDimension[1];
//            zsize = pmeGridDimension[2];
//            pmeParametersSetBasedOnEwaldErrorTolerance = 0;
//        }
//
//        gpuSetAmoebaPMEParameters(data.getAmoebaGpu(), (float) alpha, xsize, ysize, zsize);
//
//        if( data.getLog() ){
//            (void) fprintf( data.getLog(), "AmoebaMultipoleForce: PME parameters tol=%12.3e cutoff=%12.3f alpha=%12.3f [%d %d %d]\n",
//                            force.getEwaldErrorTolerance(), force.getCutoffDistance(),  alpha, xsize, ysize, zsize );
//            if( pmeParametersSetBasedOnEwaldErrorTolerance  ){
//                 (void) fprintf( data.getLog(), "Parameters based on error tolerance and OpenMM algorithm.\n" );
//            } else {
//                 double alphaT;
//                 int xsizeT, ysizeT, zsizeT;
//                 NonbondedForceImpl::calcPMEParameters(system, nb, alphaT, xsizeT, ysizeT, zsizeT);
//                 double impliedTolerance  = alpha*force.getCutoffDistance();
//                        impliedTolerance  = 0.5*exp( -(impliedTolerance*impliedTolerance) );
//                 (void) fprintf( data.getLog(), "Using input parameters implied tolerance=%12.3e;", impliedTolerance );
//                 (void) fprintf( data.getLog(), "OpenMM param: aEwald=%12.3f [%6d %6d %6d]\n", alphaT, xsizeT, ysizeT, zsizeT);
//            }
//            (void) fprintf( data.getLog(), "\n" );
//            (void) fflush( data.getLog() );
//        }
//
//        data.setApplyMultipoleCutoff( 1 );
//
//        data.cudaPlatformData.nonbondedMethod = PARTICLE_MESH_EWALD;
//        amoebaGpuContext amoebaGpu            = data.getAmoebaGpu();
//        gpuContext gpu                        = amoebaGpu->gpuContext;
//        gpu->sim.nonbondedCutoffSqr           = static_cast<float>(force.getCutoffDistance()*force.getCutoffDistance());
//        gpu->sim.nonbondedMethod              = PARTICLE_MESH_EWALD;
//    }
//    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
//}
//
//double CudaCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    computeAmoebaMultipoleForce( data );
//    return 0.0;
//}
//
//void CudaCalcAmoebaMultipoleForceKernel::getElectrostaticPotential(ContextImpl& context,  const std::vector< Vec3 >& inputGrid,
//                                                                   std::vector< double >& outputElectrostaticPotential) {
//    computeAmoebaMultipolePotential( data, inputGrid, outputElectrostaticPotential );
//    return;
//}
//
//void CudaCalcAmoebaMultipoleForceKernel::getSystemMultipoleMoments(ContextImpl& context,  const Vec3& origin,
//                                                                   std::vector< double >& outputMultipoleMonents) {
//    computeAmoebaSystemMultipoleMoments( data, origin, outputMultipoleMonents);
//    return;
//}
//
///* -------------------------------------------------------------------------- *
// *                       AmoebaGeneralizedKirkwood                            *
// * -------------------------------------------------------------------------- */
//
//class CudaCalcAmoebaGeneralizedKirkwoodForceKernel::ForceInfo : public CudaForceInfo {
//public:
//    ForceInfo(const AmoebaGeneralizedKirkwoodForce& force) : force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        double charge1, charge2, radius1, radius2, scale1, scale2;
//        force.getParticleParameters(particle1, charge1, radius1, scale1);
//        force.getParticleParameters(particle2, charge2, radius2, scale2);
//        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
//    }
//private:
//    const AmoebaGeneralizedKirkwoodForce& force;
//};
//
//CudaCalcAmoebaGeneralizedKirkwoodForceKernel::CudaCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
//           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), cu(cu), system(system) {
//    data.incrementKernelCount();
//}
//
//CudaCalcAmoebaGeneralizedKirkwoodForceKernel::~CudaCalcAmoebaGeneralizedKirkwoodForceKernel() {
//    data.decrementKernelCount();
//}
//
//void CudaCalcAmoebaGeneralizedKirkwoodForceKernel::initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) {
//
//    data.setHasAmoebaGeneralizedKirkwood( true );
//
//    int numParticles = system.getNumParticles();
//
//    std::vector<float> radius(numParticles);
//    std::vector<float> scale(numParticles);
//    std::vector<float> charge(numParticles);
//
//    for( int ii = 0; ii < numParticles; ii++ ){
//        double particleCharge, particleRadius, scalingFactor;
//        force.getParticleParameters(ii, particleCharge, particleRadius, scalingFactor);
//        radius[ii]  = static_cast<float>( particleRadius );
//        scale[ii]   = static_cast<float>( scalingFactor );
//        charge[ii]  = static_cast<float>( particleCharge );
//    }   
//    if( data.getUseGrycuk() ){
//
//        gpuSetAmoebaGrycukParameters( data.getAmoebaGpu(), static_cast<float>(force.getSoluteDielectric() ), 
//                                      static_cast<float>( force.getSolventDielectric() ), 
//                                      radius, scale, charge,
//                                      force.getIncludeCavityTerm(),
//                                      static_cast<float>( force.getProbeRadius() ), 
//                                      static_cast<float>( force.getSurfaceAreaFactor() ) ); 
//        
//    } else {
//
//        gpuSetAmoebaObcParameters( data.getAmoebaGpu(), static_cast<float>(force.getSoluteDielectric() ), 
//                                   static_cast<float>( force.getSolventDielectric() ), 
//                                   radius, scale, charge,
//                                   force.getIncludeCavityTerm(),
//                                   static_cast<float>( force.getProbeRadius() ), 
//                                   static_cast<float>( force.getSurfaceAreaFactor() ) ); 
//
//    }
//    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
//}
//
//double CudaCalcAmoebaGeneralizedKirkwoodForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    // handled in computeAmoebaMultipoleForce()
//    return 0.0;
//}
//
//static void computeAmoebaVdwForce( CudaContext& cu ) {
//
//    amoebaGpuContext gpu = data.getAmoebaGpu();
//    data.initializeGpu();
//
//    // Vdw14_7F
//    kCalculateAmoebaVdw14_7Forces(gpu, data.getUseVdwNeighborList());
//}

/* -------------------------------------------------------------------------- *
 *                           AmoebaVdw                                        *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaVdwForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaVdwForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        int iv1, iv2, class1, class2;
        double sigma1, sigma2, epsilon1, epsilon2, reduction1, reduction2;
        force.getParticleParameters(particle1, iv1, class1, sigma1, epsilon1, reduction1);
        force.getParticleParameters(particle2, iv2, class2, sigma2, epsilon2, reduction2);
        return (class1 == class2 && sigma1 == sigma2 && epsilon1 == epsilon2 && reduction1 == reduction2);
    }
private:
    const AmoebaVdwForce& force;
};

CudaCalcAmoebaVdwForceKernel::CudaCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) :
        CalcAmoebaVdwForceKernel(name, platform), cu(cu), system(system), hasInitializedNonbonded(false), sigmaEpsilon(NULL),
        bondReductionAtoms(NULL), bondReductionFactors(NULL), tempPosq(NULL), tempForces(NULL), nonbonded(NULL) {
}

CudaCalcAmoebaVdwForceKernel::~CudaCalcAmoebaVdwForceKernel() {
    cu.setAsCurrent();
    if (sigmaEpsilon != NULL)
        delete sigmaEpsilon;
    if (bondReductionAtoms != NULL)
        delete bondReductionAtoms;
    if (bondReductionFactors != NULL)
        delete bondReductionFactors;
    if (tempPosq != NULL)
        delete tempPosq;
    if (tempForces != NULL)
        delete tempForces;
    if (nonbonded != NULL)
        delete nonbonded;
}

void CudaCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {
    cu.setAsCurrent();
    sigmaEpsilon = CudaArray::create<float2>(cu, cu.getPaddedNumAtoms(), "sigmaEpsilon");
    bondReductionAtoms = CudaArray::create<int>(cu, cu.getPaddedNumAtoms(), "bondReductionAtoms");
    bondReductionFactors = CudaArray::create<float>(cu, cu.getPaddedNumAtoms(), "sigmaEpsilon");
    tempPosq = new CudaArray(cu, cu.getPaddedNumAtoms(), cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4), "tempPosq");
    tempForces = CudaArray::create<long long>(cu, 3*cu.getPaddedNumAtoms(), "tempForces");
    
    // Record atom parameters.
    
    vector<float2> sigmaEpsilonVec(cu.getPaddedNumAtoms(), make_float2(0, 1));
    vector<int> bondReductionAtomsVec(cu.getPaddedNumAtoms(), 0);
    vector<float> bondReductionFactorsVec(cu.getPaddedNumAtoms(), 0);
    vector<vector<int> > exclusions(cu.getNumAtoms());
    for (int i = 0; i < force.getNumParticles(); i++) {
        int ivIndex, classIndex;
        double sigma, epsilon, reductionFactor;
        force.getParticleParameters(i, ivIndex, classIndex, sigma, epsilon, reductionFactor);
        sigmaEpsilonVec[i] = make_float2((float) sigma, (float) epsilon);
        bondReductionAtomsVec[i] = ivIndex;
        bondReductionFactorsVec[i] = (float) reductionFactor;
        force.getParticleExclusions(i, exclusions[i]);
        exclusions[i].push_back(i);
    }
    sigmaEpsilon->upload(sigmaEpsilonVec);
    bondReductionAtoms->upload(bondReductionAtomsVec);
    bondReductionFactors->upload(bondReductionFactorsVec);
 
    // This force is applied based on modified atom positions, where hydrogens have been moved slightly
    // closer to their parent atoms.  We therefore create a separate CudaNonbondedUtilities just for
    // this force, so it will have its own neighbor list and interaction kernel.
    
    nonbonded = new CudaNonbondedUtilities(cu);
    nonbonded->addParameter(CudaNonbondedUtilities::ParameterInfo("sigmaEpsilon", "float", 2, sizeof(float2), sigmaEpsilon->getDevicePointer()));
    
    // Create the interaction kernel.
    
    map<string, string> replacements;
    string sigmaCombiningRule = force.getSigmaCombiningRule();
    if (sigmaCombiningRule == "ARITHMETIC")
        replacements["SIGMA_COMBINING_RULE"] = "1";
    else if (sigmaCombiningRule == "GEOMETRIC")
        replacements["SIGMA_COMBINING_RULE"] = "2";
    else if (sigmaCombiningRule == "CUBIC-MEAN")
        replacements["SIGMA_COMBINING_RULE"] = "3";
    else
        throw OpenMMException("Illegal combining rule for sigma: "+sigmaCombiningRule);
    string epsilonCombiningRule = force.getEpsilonCombiningRule();
    if (epsilonCombiningRule == "ARITHMETIC")
        replacements["EPILON_COMBINING_RULE"] = "1";
    else if (epsilonCombiningRule == "GEOMETRIC")
        replacements["EPILON_COMBINING_RULE"] = "2";
    else if (epsilonCombiningRule == "HARMONIC")
        replacements["EPILON_COMBINING_RULE"] = "3";
    else if (epsilonCombiningRule == "HHG")
        replacements["EPILON_COMBINING_RULE"] = "4";
    else
        throw OpenMMException("Illegal combining rule for sigma: "+sigmaCombiningRule);
    double cutoff = force.getCutoff();
    double taperCutoff = cutoff*0.9;
    replacements["CUTOFF_DISTANCE"] = cu.doubleToString(force.getCutoff());
    replacements["TAPER_CUTOFF"] = cu.doubleToString(taperCutoff);
    double cutoff2 = cutoff*cutoff;
    double taperCutoff2 = taperCutoff*taperCutoff;
    double denom = pow(cutoff-taperCutoff, -5.0);
    replacements["TAPER_C0"] = cu.doubleToString(cutoff*cutoff2 * (cutoff2-5.0*cutoff*taperCutoff+10.0*taperCutoff2)*denom);
    replacements["TAPER_C1"] = cu.doubleToString(-30.0*cutoff2*taperCutoff2*denom);
    replacements["TAPER_C2"] = cu.doubleToString(30.0*(cutoff2*taperCutoff+cutoff*taperCutoff2)*denom);
    replacements["TAPER_C3"] = cu.doubleToString(-10.0*(cutoff2+4.0*cutoff*taperCutoff+taperCutoff2)*denom);
    replacements["TAPER_C4"] = cu.doubleToString(15.0*(cutoff+taperCutoff)*denom);
    replacements["TAPER_C5"] = cu.doubleToString(-6.0*denom);
    nonbonded->addInteraction(force.getUseNeighborList(), force.getPBC(), true, force.getCutoff(), exclusions,
        cu.replaceStrings(CudaAmoebaKernelSources::amoebaVdwForce2, replacements), force.getForceGroup());
    
    // Create the other kernels.
    
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaAmoebaKernelSources::amoebaVdwForce1, defines);
    prepareKernel = cu.getKernel(module, "prepareToComputeForce");
    spreadKernel = cu.getKernel(module, "spreadForces");
    cu.addForce(new ForceInfo(force));
}

double CudaCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedNonbonded) {
        hasInitializedNonbonded = true;
        nonbonded->initialize(system);
    }
    const char* errorMessage = "Error copying array";
    CHECK_RESULT(cuMemcpyDtoDAsync(tempPosq->getDevicePointer(), cu.getPosq().getDevicePointer(), tempPosq->getSize()*tempPosq->getElementSize(), 0));
    CHECK_RESULT(cuMemcpyDtoDAsync(tempForces->getDevicePointer(), cu.getForce().getDevicePointer(), tempForces->getSize()*tempForces->getElementSize(), 0));
    void* prepareArgs[] = {&cu.getForce().getDevicePointer(), &cu.getPosq().getDevicePointer(), &tempPosq->getDevicePointer(),
        &bondReductionAtoms->getDevicePointer(), &bondReductionFactors->getDevicePointer()};
    cu.executeKernel(prepareKernel, prepareArgs, cu.getPaddedNumAtoms());
    nonbonded->prepareInteractions();
    nonbonded->computeInteractions();
    void* spreadArgs[] = {&cu.getForce().getDevicePointer(), &tempForces->getDevicePointer(), &bondReductionAtoms->getDevicePointer(), &bondReductionFactors->getDevicePointer()};
    cu.executeKernel(spreadKernel, spreadArgs, cu.getPaddedNumAtoms());
    CHECK_RESULT(cuMemcpyDtoDAsync(cu.getPosq().getDevicePointer(), tempPosq->getDevicePointer(), tempPosq->getSize()*tempPosq->getElementSize(), 0));
    CHECK_RESULT(cuMemcpyDtoDAsync(cu.getForce().getDevicePointer(), tempForces->getDevicePointer(), tempForces->getSize()*tempForces->getElementSize(), 0));
    return 0.0;
}

///* -------------------------------------------------------------------------- *
// *                           AmoebaWcaDispersion                              *
// * -------------------------------------------------------------------------- */
//
//static void computeAmoebaWcaDispersionForce( CudaContext& cu ) {
//
//    data.initializeGpu();
//    if( 0 && data.getLog() ){
//        (void) fprintf( data.getLog(), "Calling computeAmoebaWcaDispersionForce  " ); (void) fflush( data.getLog() );
//    }
//
//    kCalculateAmoebaWcaDispersionForces( data.getAmoebaGpu() );
//
//    if( 0 && data.getLog() ){
//        (void) fprintf( data.getLog(), " -- completed\n" ); (void) fflush( data.getLog() );
//    }
//}
//
//class CudaCalcAmoebaWcaDispersionForceKernel::ForceInfo : public CudaForceInfo {
//public:
//    ForceInfo(const AmoebaWcaDispersionForce& force) : force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        double radius1, radius2, epsilon1, epsilon2;
//        force.getParticleParameters(particle1, radius1, epsilon1);
//        force.getParticleParameters(particle2, radius2, epsilon2);
//        return (radius1 == radius2 && epsilon1 == epsilon2);
//    }
//private:
//    const AmoebaWcaDispersionForce& force;
//};
//
//CudaCalcAmoebaWcaDispersionForceKernel::CudaCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, CudaContext& cu, System& system) : 
//           CalcAmoebaWcaDispersionForceKernel(name, platform), cu(cu), system(system) {
//    data.incrementKernelCount();
//}
//
//CudaCalcAmoebaWcaDispersionForceKernel::~CudaCalcAmoebaWcaDispersionForceKernel() {
//    data.decrementKernelCount();
//}
//
//void CudaCalcAmoebaWcaDispersionForceKernel::initialize(const System& system, const AmoebaWcaDispersionForce& force) {
//
//    // per-particle parameters
//
//    int numParticles = system.getNumParticles();
//    std::vector<float> radii(numParticles);
//    std::vector<float> epsilons(numParticles);
//    for( int ii = 0; ii < numParticles; ii++ ){
//
//        double radius, epsilon;
//        force.getParticleParameters( ii, radius, epsilon );
//
//        radii[ii]         = static_cast<float>( radius );
//        epsilons[ii]      = static_cast<float>( epsilon );
//    }   
//    float totalMaximumDispersionEnergy =  static_cast<float>( AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy( force ) );
//    gpuSetAmoebaWcaDispersionParameters( data.getAmoebaGpu(), radii, epsilons, totalMaximumDispersionEnergy,
//                                          static_cast<float>( force.getEpso( )),
//                                          static_cast<float>( force.getEpsh( )),
//                                          static_cast<float>( force.getRmino( )),
//                                          static_cast<float>( force.getRminh( )),
//                                          static_cast<float>( force.getAwater( )),
//                                          static_cast<float>( force.getShctd( )),
//                                          static_cast<float>( force.getDispoff( ) ) );
//    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
//}
//
//double CudaCalcAmoebaWcaDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    computeAmoebaWcaDispersionForce( data );
//    return 0.0;
//}
