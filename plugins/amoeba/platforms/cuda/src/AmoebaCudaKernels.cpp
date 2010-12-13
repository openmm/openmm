/* -------------------------------------------------------------------------- *
 *                               AmoebaOpenMM                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors:                                                                   *
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
#include "openmm/internal/ContextImpl.h"
#include "kernels/amoebaGpuTypes.h"
#include "kernels/cudaKernels.h"
#include "kernels/amoebaCudaKernels.h"
#include "internal/AmoebaMultipoleForceImpl.h"
#include "internal/AmoebaWcaDispersionForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaForceInfo.h"

#include <cmath>
#ifdef _MSC_VER
#include <windows.h>
#endif

extern "C" int gpuSetConstants( gpuContext gpu );

using namespace OpenMM;
using namespace std;

/* -------------------------------------------------------------------------- *
 *                           Calculates bonded forces                         *
 * -------------------------------------------------------------------------- */

static void computeAmoebaLocalForces( AmoebaCudaData& data ) {
    amoebaGpuContext gpu = data.getAmoebaGpu();

    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "computeAmoebaLocalForces\n" ); (void) fflush( data.getLog() );
    }

    data.initializeGpu();
    kCalculateAmoebaLocalForces(gpu);

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

CudaCalcAmoebaHarmonicBondForceKernel::CudaCalcAmoebaHarmonicBondForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
                CalcAmoebaHarmonicBondForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaHarmonicBondForceKernel::~CudaCalcAmoebaHarmonicBondForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaHarmonicBondForceKernel::initialize(const System& system, const AmoebaHarmonicBondForce& force) {

    data.setAmoebaLocalForcesKernel( this );

    numBonds = force.getNumBonds();
    std::vector<int>   particle1(numBonds);
    std::vector<int>   particle2(numBonds);
    std::vector<float> length(numBonds);
    std::vector<float> quadratic(numBonds);

    for (int i = 0; i < numBonds; i++) {

        int particle1Index, particle2Index;
        double lengthValue, kValue;
        force.getBondParameters(i, particle1Index, particle2Index, lengthValue, kValue );

        particle1[i]     = particle1Index; 
        particle2[i]     = particle2Index; 
        length[i]        = static_cast<float>( lengthValue );
        quadratic[i]     = static_cast<float>( kValue );
    } 
    gpuSetAmoebaBondParameters( data.getAmoebaGpu(), particle1, particle2, length, quadratic, 
                                static_cast<float>(force.getAmoebaGlobalHarmonicBondCubic()),
                                static_cast<float>(force.getAmoebaGlobalHarmonicBondQuartic()) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaUreyBradley                               *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaUreyBradleyForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaUreyBradleyForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumInteractions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getUreyBradleyParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getUreyBradleyParameters(group1, particle1, particle2, length1, k1);
        force.getUreyBradleyParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const AmoebaUreyBradleyForce& force;
};

CudaCalcAmoebaUreyBradleyForceKernel::CudaCalcAmoebaUreyBradleyForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
                CalcAmoebaUreyBradleyForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaUreyBradleyForceKernel::~CudaCalcAmoebaUreyBradleyForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaUreyBradleyForceKernel::initialize(const System& system, const AmoebaUreyBradleyForce& force) {

    data.setAmoebaLocalForcesKernel( this );

    numInteractions = force.getNumInteractions();
    std::vector<int>   particle1(numInteractions);
    std::vector<int>   particle2(numInteractions);
    std::vector<float> length(numInteractions);
    std::vector<float> quadratic(numInteractions);

    for (int i = 0; i < numInteractions; i++) {

        int particle1Index, particle2Index;
        double lengthValue, kValue;
        force.getUreyBradleyParameters(i, particle1Index, particle2Index, lengthValue, kValue );

        particle1[i]     = particle1Index; 
        particle2[i]     = particle2Index; 
        length[i]        = static_cast<float>( lengthValue );
        quadratic[i]     = static_cast<float>( kValue );
    } 
    gpuSetAmoebaUreyBradleyParameters( data.getAmoebaGpu(), particle1, particle2, length, quadratic, 
                                       static_cast<float>(force.getAmoebaGlobalUreyBradleyCubic()),
                                       static_cast<float>(force.getAmoebaGlobalUreyBradleyQuartic()) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaUreyBradleyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
}

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

CudaCalcAmoebaHarmonicAngleForceKernel::CudaCalcAmoebaHarmonicAngleForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
            CalcAmoebaHarmonicAngleForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaHarmonicAngleForceKernel::~CudaCalcAmoebaHarmonicAngleForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaHarmonicAngleForceKernel::initialize(const System& system, const AmoebaHarmonicAngleForce& force) {

    data.setAmoebaLocalForcesKernel( this );

    numAngles                     = force.getNumAngles();
    std::vector<int> particle1(numAngles);
    std::vector<int> particle2(numAngles);
    std::vector<int> particle3(numAngles);
    std::vector<float> angle(numAngles);
    std::vector<float> k(numAngles);

    for (int i = 0; i < numAngles; i++) {
        double angleValue, kQuadratic;
        force.getAngleParameters(i, particle1[i], particle2[i], particle3[i], angleValue, kQuadratic);
        angle[i]            = static_cast<float>( angleValue );
        k[i]                = static_cast<float>( kQuadratic );
    }
    gpuSetAmoebaAngleParameters(data.getAmoebaGpu(), particle1, particle2, particle3, angle, k,
                                static_cast<float>(force.getAmoebaGlobalHarmonicAngleCubic()),
                                static_cast<float>(force.getAmoebaGlobalHarmonicAngleQuartic()),
                                static_cast<float>(force.getAmoebaGlobalHarmonicAnglePentic()),
                                static_cast<float>(force.getAmoebaGlobalHarmonicAngleSextic()) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
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

CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::CudaCalcAmoebaHarmonicInPlaneAngleForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
          CalcAmoebaHarmonicInPlaneAngleForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::~CudaCalcAmoebaHarmonicInPlaneAngleForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::initialize(const System& system, const AmoebaHarmonicInPlaneAngleForce& force) {

    data.setAmoebaLocalForcesKernel( this );

    numAngles = force.getNumAngles();

    std::vector<int> particle1(numAngles);
    std::vector<int> particle2(numAngles);
    std::vector<int> particle3(numAngles);
    std::vector<int> particle4(numAngles);
    std::vector<float> angle(numAngles);
    std::vector<float> k(numAngles);

    for (int i = 0; i < numAngles; i++) {
        double angleValue, kQuadratic;
        force.getAngleParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], angleValue, kQuadratic);
        //angle[i]            = static_cast<float>( (angleValue*RadiansToDegrees) );
        angle[i]            = static_cast<float>( angleValue );
        k[i]                = static_cast<float>( kQuadratic );
    }
    gpuSetAmoebaInPlaneAngleParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, angle, k,
                                       static_cast<float>( force.getAmoebaGlobalHarmonicInPlaneAngleCubic()),
                                       static_cast<float>( force.getAmoebaGlobalHarmonicInPlaneAngleQuartic()),
                                       static_cast<float>( force.getAmoebaGlobalHarmonicInPlaneAnglePentic()),
                                       static_cast<float>( force.getAmoebaGlobalHarmonicInPlaneAngleSextic() ) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaHarmonicInPlaneAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaHarmonicTorsion                            *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        vector<double> torsion1, torsion2, torsion3;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, torsion1, torsion2, torsion3);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        vector<double> torsion11, torsion21, torsion31;
        vector<double> torsion12, torsion22, torsion32;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, torsion11, torsion21, torsion31);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, torsion12, torsion22, torsion32);
        for (int i = 0; i < (int) torsion11.size(); ++i)
            if (torsion11[i] != torsion12[i])
                return false;
        for (int i = 0; i < (int) torsion21.size(); ++i)
            if (torsion21[i] != torsion22[i])
                return false;
        for (int i = 0; i < (int) torsion31.size(); ++i)
            if (torsion31[i] != torsion32[i])
                return false;
        return true;
    }
private:
    const AmoebaTorsionForce& force;
};

CudaCalcAmoebaTorsionForceKernel::CudaCalcAmoebaTorsionForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
             CalcAmoebaTorsionForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaTorsionForceKernel::~CudaCalcAmoebaTorsionForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaTorsionForceKernel::initialize(const System& system, const AmoebaTorsionForce& force) {

    data.setAmoebaLocalForcesKernel( this );
    numTorsions                     = force.getNumTorsions();
    std::vector<int> particle1(numTorsions);
    std::vector<int> particle2(numTorsions);
    std::vector<int> particle3(numTorsions);
    std::vector<int> particle4(numTorsions);

    std::vector< std::vector<float> > torsionParameters1(numTorsions);
    std::vector< std::vector<float> > torsionParameters2(numTorsions);
    std::vector< std::vector<float> > torsionParameters3(numTorsions);

    for (int i = 0; i < numTorsions; i++) {

        std::vector<double> torsionParameter1;
        std::vector<double> torsionParameter2;
        std::vector<double> torsionParameter3;

        std::vector<float> torsionParameters1F(3);
        std::vector<float> torsionParameters2F(3);
        std::vector<float> torsionParameters3F(3);

        force.getTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], torsionParameter1, torsionParameter2, torsionParameter3 );
        for ( unsigned int jj = 0; jj < torsionParameter1.size(); jj++) {
            torsionParameters1F[jj] = static_cast<float>(torsionParameter1[jj]);
            torsionParameters2F[jj] = static_cast<float>(torsionParameter2[jj]);
            torsionParameters3F[jj] = static_cast<float>(torsionParameter3[jj]);
        }
        torsionParameters1[i] = torsionParameters1F;
        torsionParameters2[i] = torsionParameters2F;
        torsionParameters3[i] = torsionParameters3F;
    }
    gpuSetAmoebaTorsionParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, torsionParameters1, torsionParameters2, torsionParameters3 );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaHarmonicPiTorsion                          *
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

CudaCalcAmoebaPiTorsionForceKernel::CudaCalcAmoebaPiTorsionForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
         CalcAmoebaPiTorsionForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaPiTorsionForceKernel::~CudaCalcAmoebaPiTorsionForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaPiTorsionForceKernel::initialize(const System& system, const AmoebaPiTorsionForce& force) {

    data.setAmoebaLocalForcesKernel( this );
    numPiTorsions                     = force.getNumPiTorsions();

    std::vector<int> particle1(numPiTorsions);
    std::vector<int> particle2(numPiTorsions);
    std::vector<int> particle3(numPiTorsions);
    std::vector<int> particle4(numPiTorsions);
    std::vector<int> particle5(numPiTorsions);
    std::vector<int> particle6(numPiTorsions);

    std::vector<float> torsionKParameters(numPiTorsions);

    for (int i = 0; i < numPiTorsions; i++) {

        double torsionKParameter;

        force.getPiTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], particle5[i], particle6[i], torsionKParameter);
        torsionKParameters[i] = static_cast<float>(torsionKParameter);
    }
    gpuSetAmoebaPiTorsionParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, particle5, particle6, torsionKParameters);
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaPiTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
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

CudaCalcAmoebaStretchBendForceKernel::CudaCalcAmoebaStretchBendForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
                   CalcAmoebaStretchBendForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaStretchBendForceKernel::~CudaCalcAmoebaStretchBendForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaStretchBendForceKernel::initialize(const System& system, const AmoebaStretchBendForce& force) {

    data.setAmoebaLocalForcesKernel( this );
    numStretchBends                     = force.getNumStretchBends();

    std::vector<int>   particle1(numStretchBends);
    std::vector<int>   particle2(numStretchBends);
    std::vector<int>   particle3(numStretchBends);
    std::vector<float> lengthABParameters(numStretchBends);
    std::vector<float> lengthCBParameters(numStretchBends);
    std::vector<float> angleParameters(numStretchBends);
    std::vector<float> kParameters(numStretchBends);

    for (int i = 0; i < numStretchBends; i++) {

        double lengthAB, lengthCB, angle, k;

        force.getStretchBendParameters(i, particle1[i], particle2[i], particle3[i], lengthAB, lengthCB, angle, k);
        lengthABParameters[i] = static_cast<float>(lengthAB);
        lengthCBParameters[i] = static_cast<float>(lengthCB);
        angleParameters[i]    = static_cast<float>(angle);
        kParameters[i]        = static_cast<float>(k);
    }
    gpuSetAmoebaStretchBendParameters(data.getAmoebaGpu(), particle1, particle2, particle3, lengthABParameters, lengthCBParameters, angleParameters, kParameters);
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaStretchBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
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

CudaCalcAmoebaOutOfPlaneBendForceKernel::CudaCalcAmoebaOutOfPlaneBendForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
          CalcAmoebaOutOfPlaneBendForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaOutOfPlaneBendForceKernel::~CudaCalcAmoebaOutOfPlaneBendForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaOutOfPlaneBendForceKernel::initialize(const System& system, const AmoebaOutOfPlaneBendForce& force) {

    data.setAmoebaLocalForcesKernel( this );
    numOutOfPlaneBends                     = force.getNumOutOfPlaneBends();

    std::vector<int>   particle1(numOutOfPlaneBends);
    std::vector<int>   particle2(numOutOfPlaneBends);
    std::vector<int>   particle3(numOutOfPlaneBends);
    std::vector<int>   particle4(numOutOfPlaneBends);
    std::vector<float> kParameters(numOutOfPlaneBends);

    for (int i = 0; i < numOutOfPlaneBends; i++) {

        double k;

        force.getOutOfPlaneBendParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], k);
        kParameters[i] = static_cast<float>(k);
    }
    gpuSetAmoebaOutOfPlaneBendParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, kParameters,
                                         static_cast<float>( force.getAmoebaGlobalOutOfPlaneBendCubic()),
                                         static_cast<float>( force.getAmoebaGlobalOutOfPlaneBendQuartic()),
                                         static_cast<float>( force.getAmoebaGlobalOutOfPlaneBendPentic()),
                                         static_cast<float>( force.getAmoebaGlobalOutOfPlaneBendSextic() ) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaOutOfPlaneBendForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
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

CudaCalcAmoebaTorsionTorsionForceKernel::CudaCalcAmoebaTorsionTorsionForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
                CalcAmoebaTorsionTorsionForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaTorsionTorsionForceKernel::~CudaCalcAmoebaTorsionTorsionForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaTorsionTorsionForceKernel::initialize(const System& system, const AmoebaTorsionTorsionForce& force) {

    data.setAmoebaLocalForcesKernel( this );
    numTorsionTorsions = force.getNumTorsionTorsions();

    // torsion-torsion parameters

    std::vector<int>   particle1(numTorsionTorsions);
    std::vector<int>   particle2(numTorsionTorsions);
    std::vector<int>   particle3(numTorsionTorsions);
    std::vector<int>   particle4(numTorsionTorsions);
    std::vector<int>   particle5(numTorsionTorsions);
    std::vector<int>   chiralCheckAtomIndex(numTorsionTorsions);
    std::vector<int>   gridIndices(numTorsionTorsions);

    for (int i = 0; i < numTorsionTorsions; i++) {
        force.getTorsionTorsionParameters(i, particle1[i], particle2[i], particle3[i],
                                             particle4[i], particle5[i],
                                             chiralCheckAtomIndex[i], gridIndices[i]);
    }
    gpuSetAmoebaTorsionTorsionParameters(data.getAmoebaGpu(), particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndices );

    // torsion-torsion grids

    numTorsionTorsionGrids = force.getNumTorsionTorsionGrids();
    std::vector< std::vector< std::vector< std::vector<float> > > > floatGrids;

    floatGrids.resize(numTorsionTorsionGrids);
    for (int gridIndex = 0; gridIndex < numTorsionTorsionGrids; gridIndex++) {

        const TorsionTorsionGrid& grid = force.getTorsionTorsionGrid( gridIndex );

        floatGrids[gridIndex].resize( grid.size() );
        for (unsigned int ii = 0; ii < grid.size(); ii++) {

            floatGrids[gridIndex][ii].resize( grid[ii].size() );
            for (unsigned int jj = 0; jj < grid[ii].size(); jj++) {

                floatGrids[gridIndex][ii][jj].resize( grid[ii][jj].size() );
                for (unsigned int kk = 0; kk < grid[ii][kk].size(); kk++) {
                    floatGrids[gridIndex][ii][jj][kk] = static_cast<float>(grid[ii][jj][kk]);
                }
            }
        }
    }
    gpuSetAmoebaTorsionTorsionGrids(data.getAmoebaGpu(), floatGrids );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaTorsionTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( data.getAmoebaLocalForcesKernel() == this ){
        computeAmoebaLocalForces( data );
    }
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                             AmoebaMultipole                                *
 * -------------------------------------------------------------------------- */

static void computeAmoebaMultipoleForce( AmoebaCudaData& data ) {

    amoebaGpuContext gpu = data.getAmoebaGpu();
    if( data.getMultipoleForceCount() == 0 ){
        gpuCopyWorkUnit( gpu );
    }
    //if( data.getApplyCutoff() && (data.getMultipoleForceCount() % 100) == 0 ){
        //gpuReorderAtoms(gpu->gpuContext);
    //}
    data.incrementMultipoleForceCount();
    data.initializeGpu();

    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "computeAmoebaMultipoleForce\n" );
        (void) fflush( data.getLog());
    }

    // calculate Born radii

    if( data.getHasAmoebaGeneralizedKirkwood() ){
        kCalculateObcGbsaBornSum(gpu->gpuContext);
        kReduceObcGbsaBornSum(gpu->gpuContext);
    }   

    // multipoles

    kCalculateAmoebaMultipoleForces(gpu, data.getHasAmoebaGeneralizedKirkwood() );

//kClearForces(gpu->gpuContext);
//kClearEnergy(gpu->gpuContext);
//(void) fprintf( data.getLog(), "computeAmoebaMultipoleForce clearing forces/energy after kCalculateAmoebaMultipoleForces()\n" );

    // GK

    if( data.getHasAmoebaGeneralizedKirkwood() ){
        kCalculateAmoebaKirkwood(gpu);
    }

    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "completed computeAmoebaMultipoleForce\n" );
        (void) fflush( data.getLog());
    }
}

class CudaCalcAmoebaMultipoleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaMultipoleForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, thole1, thole2, damping1, damping2, polarity1, polarity2;
        int axis1, axis2, multipole11, multipole12, multipole21, multipole22, multipole31, multipole32;
        vector<double> dipole1, dipole2, quadrupole1, quadrupole2;
        force.getMultipoleParameters(particle1, charge1, dipole1, quadrupole1, axis1, multipole11, multipole21, multipole31, thole1, damping1, polarity1);
        force.getMultipoleParameters(particle2, charge2, dipole2, quadrupole2, axis2, multipole12, multipole22, multipole32, thole2, damping2, polarity2);
        if (charge1 != charge2 || thole1 != thole2 || damping1 != damping2 || polarity1 != polarity2 || axis1 != axis2)
            return false;
        for (int i = 0; i < (int) dipole1.size(); ++i)
            if (dipole1[i] != dipole2[i])
                return false;
        for (int i = 0; i < (int) quadrupole1.size(); ++i)
            if (quadrupole1[i] != quadrupole2[i])
                return false;
        return true;
    }
private:
    const AmoebaMultipoleForce& force;
};

CudaCalcAmoebaMultipoleForceKernel::CudaCalcAmoebaMultipoleForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
         CalcAmoebaMultipoleForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaMultipoleForceKernel::~CudaCalcAmoebaMultipoleForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaMultipoleForceKernel::initialize(const System& system, const AmoebaMultipoleForce& force) {

    numMultipoles   = force.getNumMultipoles();

    data.setHasAmoebaMultipole( true );

    std::vector<float> charges(numMultipoles);
    std::vector<float> dipoles(3*numMultipoles);
    std::vector<float> quadrupoles(9*numMultipoles);
    std::vector<float> tholes(numMultipoles);
    std::vector<float> dampingFactors(numMultipoles);
    std::vector<float> polarity(numMultipoles);
    std::vector<int>   axisTypes(numMultipoles);
    std::vector<int>   multipoleAtomZs(numMultipoles);
    std::vector<int>   multipoleAtomXs(numMultipoles);
    std::vector<int>   multipoleAtomYs(numMultipoles);
    std::vector< std::vector< std::vector<int> > > multipoleAtomCovalentInfo(numMultipoles);
    std::vector<int> minCovalentIndices(numMultipoles);
    std::vector<int> minCovalentPolarizationIndices(numMultipoles);

    float scalingDistanceCutoff = static_cast<float>(force.getScalingDistanceCutoff());

    std::vector<AmoebaMultipoleForce::CovalentType> covalentList;
    covalentList.push_back( AmoebaMultipoleForce::Covalent12 );
    covalentList.push_back( AmoebaMultipoleForce::Covalent13 );
    covalentList.push_back( AmoebaMultipoleForce::Covalent14 );
    covalentList.push_back( AmoebaMultipoleForce::Covalent15 );

    std::vector<AmoebaMultipoleForce::CovalentType> polarizationCovalentList;
    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent11 );
    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent12 );
    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent13 );
    polarizationCovalentList.push_back( AmoebaMultipoleForce::PolarizationCovalent14 );

    std::vector<int> covalentDegree;
    AmoebaMultipoleForceImpl::getCovalentDegree( force, covalentDegree );
    int dipoleIndex      = 0;
    int quadrupoleIndex  = 0;
    int maxCovalentRange = 0;
    double totalCharge   = 0.0;
    for (int i = 0; i < numMultipoles; i++) {

        // multipoles

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, tholeD, dampingFactorD, polarityD;
        std::vector<double> dipolesD;
        std::vector<double> quadrupolesD;
        force.getMultipoleParameters(i, charge, dipolesD, quadrupolesD, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                     tholeD, dampingFactorD, polarityD );

        totalCharge                       += charge;
        axisTypes[i]                       = axisType;
        multipoleAtomZs[i]                 = multipoleAtomZ;
        multipoleAtomXs[i]                 = multipoleAtomX;
        multipoleAtomYs[i]                 = multipoleAtomY;

        charges[i]                         = static_cast<float>(charge);
        tholes[i]                          = static_cast<float>(tholeD);
        dampingFactors[i]                  = static_cast<float>(dampingFactorD);
        polarity[i]                        = static_cast<float>(polarityD);

        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[0]);
        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[1]);
        dipoles[dipoleIndex++]             = static_cast<float>(dipolesD[2]);
        
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[0]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[1]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[2]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[3]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[4]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[5]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[6]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[7]);
        quadrupoles[quadrupoleIndex++]     = static_cast<float>(quadrupolesD[8]);

        // covalent info

        std::vector< std::vector<int> > covalentLists;
        force.getCovalentMaps(i, covalentLists );
        multipoleAtomCovalentInfo[i] = covalentLists;

        int minCovalentIndex, maxCovalentIndex;
        AmoebaMultipoleForceImpl::getCovalentRange( force, i, covalentList, &minCovalentIndex, &maxCovalentIndex );
        minCovalentIndices[i] = minCovalentIndex;
        if( maxCovalentRange < (maxCovalentIndex - minCovalentIndex) ){
            maxCovalentRange = maxCovalentIndex - minCovalentIndex;
        }

        AmoebaMultipoleForceImpl::getCovalentRange( force, i, polarizationCovalentList, &minCovalentIndex, &maxCovalentIndex );
        minCovalentPolarizationIndices[i] = minCovalentIndex;
        if( maxCovalentRange < (maxCovalentIndex - minCovalentIndex) ){
            maxCovalentRange = maxCovalentIndex - minCovalentIndex;
        }
    }

    int iterativeMethod = static_cast<int>(force.getMutualInducedIterationMethod());
    if( iterativeMethod != 0 ){
         throw OpenMMException("Iterative method for mutual induced dipoles not recognized.\n");
    }

    int nonbondedMethod = static_cast<int>(force.getNonbondedMethod());
    if( nonbondedMethod != 0 && nonbondedMethod != 1 ){
         throw OpenMMException("AmoebaMultipoleForce nonbonded method not recognized.\n");
    }

    gpuSetAmoebaMultipoleParameters(data.getAmoebaGpu(), charges, dipoles, quadrupoles, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                    tholes, scalingDistanceCutoff, dampingFactors, polarity,
                                    multipoleAtomCovalentInfo, covalentDegree, minCovalentIndices, minCovalentPolarizationIndices, (maxCovalentRange+2),
                                    static_cast<int>(force.getMutualInducedIterationMethod()),
                                    force.getMutualInducedMaxIterations(),
                                    static_cast<float>( force.getMutualInducedTargetEpsilon()),
                                    nonbondedMethod,
                                    static_cast<float>( force.getCutoffDistance()),
                                    static_cast<float>( force.getAEwald()),
                                    static_cast<float>( force.getElectricConstant()) );
    if (nonbondedMethod == AmoebaMultipoleForce::PME) {
        double alpha;
        int xsize, ysize, zsize;
        NonbondedForce nb;
        nb.setEwaldErrorTolerance(force.getEwaldErrorTolerance());
        nb.setCutoffDistance(force.getCutoffDistance());
        std::vector<int> pmeGridDimension;
        force.getPmeGridDimensions( pmeGridDimension );
        if( 1 || pmeGridDimension[0] == 0 ){
            NonbondedForceImpl::calcPMEParameters(system, nb, alpha, xsize, ysize, zsize);
/*
alpha = 5.446;
xsize = 60;
ysize = 48;
zsize = 48;
*/
        } else {
            alpha = force.getAEwald();
            xsize = pmeGridDimension[0];
            ysize = pmeGridDimension[1];
            zsize = pmeGridDimension[2];
        }
        if( data.getLog() ){
            (void) fprintf( data.getLog(), "AmoebaMultipoleForce: PME parameters tol=%12.3e cutoff=%12.3f alpha=%12.3f [%d %d %d]\n",
                            force.getEwaldErrorTolerance(), force.getCutoffDistance(),  alpha, xsize, ysize, zsize );
            (void) fflush( data.getLog() );
        }
        gpuSetAmoebaPMEParameters(data.getAmoebaGpu(), (float) alpha, xsize, ysize, zsize);
        data.setApplyMultipoleCutoff( 1 );

        data.cudaPlatformData.nonbondedMethod = PARTICLE_MESH_EWALD;
        amoebaGpuContext amoebaGpu            = data.getAmoebaGpu();
        gpuContext gpu                        = amoebaGpu->gpuContext;
        gpu->sim.nonbondedCutoffSqr           = force.getCutoffDistance()*force.getCutoffDistance();
        gpu->sim.nonbondedMethod              = PARTICLE_MESH_EWALD;
    }
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaMultipoleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    computeAmoebaMultipoleForce( data );
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                       AmoebaGeneralizedKirkwood                            *
 * -------------------------------------------------------------------------- */

class CudaCalcAmoebaGeneralizedKirkwoodForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaGeneralizedKirkwoodForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, scale1, scale2;
        force.getParticleParameters(particle1, charge1, radius1, scale1);
        force.getParticleParameters(particle2, charge2, radius2, scale2);
        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
    }
private:
    const AmoebaGeneralizedKirkwoodForce& force;
};

CudaCalcAmoebaGeneralizedKirkwoodForceKernel::CudaCalcAmoebaGeneralizedKirkwoodForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
           CalcAmoebaGeneralizedKirkwoodForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaGeneralizedKirkwoodForceKernel::~CudaCalcAmoebaGeneralizedKirkwoodForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaGeneralizedKirkwoodForceKernel::initialize(const System& system, const AmoebaGeneralizedKirkwoodForce& force) {

    data.setHasAmoebaGeneralizedKirkwood( true );

    int numParticles = system.getNumParticles();

    std::vector<float> radius(numParticles);
    std::vector<float> scale(numParticles);
    std::vector<float> charge(numParticles);

    for( int ii = 0; ii < numParticles; ii++ ){
        double particleCharge, particleRadius, scalingFactor;
        force.getParticleParameters(ii, particleCharge, particleRadius, scalingFactor);
        radius[ii]  = static_cast<float>( particleRadius );
        scale[ii]   = static_cast<float>( scalingFactor );
        charge[ii]  = static_cast<float>( particleCharge );
    }   
    gpuSetAmoebaObcParameters( data.getAmoebaGpu(), static_cast<float>(force.getSoluteDielectric() ), 
                               static_cast<float>( force.getSolventDielectric() ), 
                               static_cast<float>( force.getDielectricOffset() ), radius, scale, charge,
                               force.getIncludeCavityTerm(),
                               static_cast<float>( force.getProbeRadius() ), 
                               static_cast<float>( force.getSurfaceAreaFactor() ) ); 
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaGeneralizedKirkwoodForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // handled in computeAmoebaMultipoleForce()
    return 0.0;
}

static void computeAmoebaVdwForce( AmoebaCudaData& data ) {

    amoebaGpuContext gpu = data.getAmoebaGpu();
    data.initializeGpu();

    // Vdw14_7F

    kCalculateAmoebaVdw14_7Forces(gpu, data.getUseVdwNeighborList());
}

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
        return (iv1 == iv2 && class1 == class2 && sigma1 == sigma2 && epsilon1 == epsilon2 && reduction1 == reduction2);
    }
private:
    const AmoebaVdwForce& force;
};

CudaCalcAmoebaVdwForceKernel::CudaCalcAmoebaVdwForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) :
       CalcAmoebaVdwForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaVdwForceKernel::~CudaCalcAmoebaVdwForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaVdwForceKernel::initialize(const System& system, const AmoebaVdwForce& force) {

    // per-particle parameters

    int numParticles = system.getNumParticles();

    std::vector<int> indexIVs(numParticles);
    std::vector<int> indexClasses(numParticles);
    std::vector< std::vector<int> > allExclusions(numParticles);
    std::vector<float> sigmas(numParticles);
    std::vector<float> epsilons(numParticles);
    std::vector<float> reductions(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        int indexIV, indexClass;
        double sigma, epsilon, reduction;
        std::vector<int> exclusions;

        force.getParticleParameters( ii, indexIV, indexClass, sigma, epsilon, reduction );
        force.getParticleExclusions( ii, exclusions );
        for( unsigned int jj = 0; jj < exclusions.size(); jj++ ){
           allExclusions[ii].push_back( exclusions[jj] );
        }

        indexIVs[ii]      = indexIV;
        indexClasses[ii]  = indexClass;
        sigmas[ii]        = static_cast<float>( sigma );
        epsilons[ii]      = static_cast<float>( epsilon );
        reductions[ii]    = static_cast<float>( reduction );
    }   

    gpuSetAmoebaVdwParameters( data.getAmoebaGpu(), indexIVs, indexClasses, sigmas, epsilons, reductions,
                               force.getSigmaCombiningRule(), force.getEpsilonCombiningRule(),
                               allExclusions, force.getPBC(), static_cast<float>(force.getCutoff()) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
    if( data.getLog() ){
        (void) fprintf( data.getLog(), "CudaCalcAmoebaVdwForceKernel PBC=%d getUseNeighborList=%d\n",
                        force.getPBC(), force.getUseNeighborList() );
    }
    data.setUseVdwNeighborList( force.getUseNeighborList() );
}

double CudaCalcAmoebaVdwForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    computeAmoebaVdwForce( data );
    return 0.0;
}

/* -------------------------------------------------------------------------- *
 *                           AmoebaWcaDispersion                              *
 * -------------------------------------------------------------------------- */

static void computeAmoebaWcaDispersionForce( AmoebaCudaData& data ) {

    data.initializeGpu();
    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "Calling computeAmoebaWcaDispersionForce  " ); (void) fflush( data.getLog() );
    }

    kCalculateAmoebaWcaDispersionForces( data.getAmoebaGpu() );

    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), " -- completed\n" ); (void) fflush( data.getLog() );
    }
}

class CudaCalcAmoebaWcaDispersionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const AmoebaWcaDispersionForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double radius1, radius2, epsilon1, epsilon2;
        force.getParticleParameters(particle1, radius1, epsilon1);
        force.getParticleParameters(particle2, radius2, epsilon2);
        return (radius1 == radius2 && epsilon1 == epsilon2);
    }
private:
    const AmoebaWcaDispersionForce& force;
};

CudaCalcAmoebaWcaDispersionForceKernel::CudaCalcAmoebaWcaDispersionForceKernel(std::string name, const Platform& platform, AmoebaCudaData& data, System& system) : 
           CalcAmoebaWcaDispersionForceKernel(name, platform), data(data), system(system) {
    data.incrementKernelCount();
}

CudaCalcAmoebaWcaDispersionForceKernel::~CudaCalcAmoebaWcaDispersionForceKernel() {
    data.decrementKernelCount();
}

void CudaCalcAmoebaWcaDispersionForceKernel::initialize(const System& system, const AmoebaWcaDispersionForce& force) {

    // per-particle parameters

    int numParticles = system.getNumParticles();
    std::vector<float> radii(numParticles);
    std::vector<float> epsilons(numParticles);
    for( int ii = 0; ii < numParticles; ii++ ){

        double radius, epsilon;
        force.getParticleParameters( ii, radius, epsilon );

        radii[ii]         = static_cast<float>( radius );
        epsilons[ii]      = static_cast<float>( epsilon );
    }   
    float totalMaximumDispersionEnergy =  static_cast<float>( AmoebaWcaDispersionForceImpl::getTotalMaximumDispersionEnergy( force ) );
    gpuSetAmoebaWcaDispersionParameters( data.getAmoebaGpu(), radii, epsilons, totalMaximumDispersionEnergy,
                                          static_cast<float>( force.getEpso( )),
                                          static_cast<float>( force.getEpsh( )),
                                          static_cast<float>( force.getRmino( )),
                                          static_cast<float>( force.getRminh( )),
                                          static_cast<float>( force.getAwater( )),
                                          static_cast<float>( force.getShctd( )),
                                          static_cast<float>( force.getDispoff( ) ) );
    data.getAmoebaGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaCalcAmoebaWcaDispersionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    computeAmoebaWcaDispersionForce( data );
    return 0.0;
}
