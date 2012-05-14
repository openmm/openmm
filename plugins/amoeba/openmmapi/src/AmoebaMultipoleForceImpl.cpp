/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
 * Contributors:                                                              *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/AmoebaMultipoleForceImpl.h"
#include "openmm/amoebaKernels.h"

using namespace OpenMM;

using std::vector;

bool AmoebaMultipoleForceImpl::initializedCovalentDegrees = false;
int AmoebaMultipoleForceImpl::CovalentDegrees[]           = { 1,2,3,4,0,1,2,3};

AmoebaMultipoleForceImpl::AmoebaMultipoleForceImpl(AmoebaMultipoleForce& owner) : owner(owner) {
}

AmoebaMultipoleForceImpl::~AmoebaMultipoleForceImpl() {
}

void AmoebaMultipoleForceImpl::initialize(ContextImpl& context) {

    System& system = context.getSystem();
    if (owner.getNumMultipoles() != system.getNumParticles())
        throw OpenMMException("AmoebaMultipoleForce must have exactly as many particles as the System it belongs to.");

    for( int ii = 0; ii < system.getNumParticles(); ii++ ){

        int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
        double charge, thole, dampingFactor, polarity ;
        std::vector<double> molecularDipole;
        std::vector<double> molecularQuadrupole;

        owner.getMultipoleParameters( ii, charge, molecularDipole, molecularQuadrupole, axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY,
                                      thole, dampingFactor, polarity );

       // only 'Z-then-X', 'Bisector', Z-Bisect, ThreeFold  currently handled

        if( axisType != AmoebaMultipoleForce::ZThenX     && axisType != AmoebaMultipoleForce::Bisector &&
            axisType != AmoebaMultipoleForce::ZBisect    && axisType != AmoebaMultipoleForce::ThreeFold &&
            axisType != AmoebaMultipoleForce::ZOnly      && axisType != AmoebaMultipoleForce::NoAxisType ) {
             std::stringstream buffer;
             buffer << "AmoebaMultipoleForce: axis type=" << axisType;
             buffer << " not currently handled - only axisTypes[ ";
             buffer << AmoebaMultipoleForce::ZThenX   << ", " << AmoebaMultipoleForce::Bisector  << ", ";
             buffer << AmoebaMultipoleForce::ZBisect  << ", " << AmoebaMultipoleForce::ThreeFold << ", ";
             buffer << AmoebaMultipoleForce::NoAxisType;
             buffer << "] (ZThenX, Bisector, Z-Bisect, ThreeFold, NoAxisType) currently handled .";
             throw OpenMMException(buffer.str());
        }
    }
    kernel = context.getPlatform().createKernel(CalcAmoebaMultipoleForceKernel::Name(), context);
    kernel.getAs<CalcAmoebaMultipoleForceKernel>().initialize(context.getSystem(), owner);
}

double AmoebaMultipoleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcAmoebaMultipoleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> AmoebaMultipoleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcAmoebaMultipoleForceKernel::Name());
    return names;
}

const int* AmoebaMultipoleForceImpl::getCovalentDegrees( void ) {
    if( !initializedCovalentDegrees ){
        initializedCovalentDegrees                                      = true;
        CovalentDegrees[AmoebaMultipoleForce::Covalent12]               = 1;
        CovalentDegrees[AmoebaMultipoleForce::Covalent13]               = 2;
        CovalentDegrees[AmoebaMultipoleForce::Covalent14]               = 3;
        CovalentDegrees[AmoebaMultipoleForce::Covalent15]               = 4;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent11]   = 0;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent12]   = 1;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent13]   = 2;
        CovalentDegrees[AmoebaMultipoleForce::PolarizationCovalent14]   = 3;
    }
    return CovalentDegrees;
}

void AmoebaMultipoleForceImpl::getCovalentRange( const AmoebaMultipoleForce& force, int atomIndex, const std::vector<AmoebaMultipoleForce::CovalentType>& lists,
                                                 int* minCovalentIndex, int* maxCovalentIndex ){

    *minCovalentIndex =  999999999;
    *maxCovalentIndex = -999999999;
    for( unsigned int kk = 0; kk < lists.size(); kk++ ){
        AmoebaMultipoleForce::CovalentType jj = lists[kk];
        std::vector<int> covalentList;
        force.getCovalentMap( atomIndex, jj, covalentList );
        for( unsigned int ii = 0; ii < covalentList.size(); ii++ ){
            if( *minCovalentIndex > covalentList[ii] ){
               *minCovalentIndex = covalentList[ii];
            }
            if( *maxCovalentIndex < covalentList[ii] ){
               *maxCovalentIndex = covalentList[ii];
            }
        }
    }   
    return;
}

void AmoebaMultipoleForceImpl::getCovalentDegree( const AmoebaMultipoleForce& force, std::vector<int>& covalentDegree ){
    covalentDegree.resize( AmoebaMultipoleForce::CovalentEnd );
    const int* CovalentDegrees = AmoebaMultipoleForceImpl::getCovalentDegrees();
    for( unsigned int kk = 0; kk < AmoebaMultipoleForce::CovalentEnd; kk++ ){
        covalentDegree[kk] = CovalentDegrees[kk];
    }   
    return;
}
