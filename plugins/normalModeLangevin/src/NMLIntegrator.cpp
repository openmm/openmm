/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Chris Sweet                                                       *
 * Contributors: Christopher Bruns                                            *
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

#include "NMLIntegrator.h"
#include "IntegrateNMLStepKernel.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

NMLIntegrator::NMLIntegrator(double temperature, double frictionCoeff, double stepSize, ProtoMol::EigenvectorInfo* projectionVectorInfo ) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setProjectionVectorInfo(projectionVectorInfo);
    setConstraintTolerance(1e-4);
}

void NMLIntegrator::initialize(ContextImpl& contextRef) {
    context = &contextRef;
    kernel = context->getPlatform().createKernel(IntegrateNMLStepKernel::Name(), contextRef);
    dynamic_cast<IntegrateNMLStepKernel&>(kernel.getImpl()).initialize(contextRef.getSystem(), *this);
}

vector<string> NMLIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateNMLStepKernel::Name());
    return names;
}

void NMLIntegrator::step(int steps) {
    

    std::cout << "One loop " << std::endl;
    
    //updates per loop. Need it here? ####Moved within loop. Single v multiple steps error if not!
    //context->updateContextState();
    //context->calcForces();
    
    //loop
    for (int i = 0; i < steps; ++i) {

        //updates per loop. Need it here?
        context->updateContextState();
        context->calcForces();

        //do half kick and full drift
        dynamic_cast<IntegrateNMLStepKernel&>(kernel.getImpl()).execute(*context, *this, 0.0, 1 );    //stepType 1 is halfKick/drift LangevinLeapfrog
        //in case projection vectors changed, clear flag
        setProjVecChanged(false);

        //updates per loop, ####done in minimizer
        //context->updateContextState();
        //context->calcForces();
        
        //minimize compliment space
        minimize(50);

        //do half kick
        dynamic_cast<IntegrateNMLStepKernel&>(kernel.getImpl()).execute(*context, *this, 0.0, 2 );    //stepType 1 is halfKick LangevinLeapfrog
              
    }
    
    //update time
    context->setTime(context->getTime()+getStepSize() * steps);
    
}

void NMLIntegrator::minimize(int maxsteps) {
    
    //minimum limit
    const double minlim = getMinimumLimit();
    
    //loop
    for (int i = 0; i < maxsteps; ++i) {
        
        //updates per loop
        context->updateContextState();
        context->calcForces();
        
        //get initial PE
        const double initialPE = context->calcPotentialEnergy();
        
        //minimize
        dynamic_cast<IntegrateNMLStepKernel&>(kernel.getImpl()).execute(*context, *this, initialPE, 3 );    //stepType 3 is simple minimizer
        setProjVecChanged(false);
        
        //get PE difference
        const double currentPE = context->calcPotentialEnergy();
        const double peDifference =  initialPE - currentPE;
        
        //std::cout << "Loop " << i << ", Current " << currentPE << ", Old " << initialPE << ", lim " << minlim << std::endl;
        //end condition met?
        if(peDifference < minlim && peDifference >= 0.0){
            //std::cout << "Breaking at " << i << ", diff " << peDifference << std::endl;
            break;
        }
        
        //Require quadratic solution
        if(peDifference < 0.0){
            //std::cout << "Quadratic at " << i << std::endl;
            
            //update for new positions
            context->updateContextState();
            context->calcForces();
            
            //minimize, uses quadratic soultion as 'stepsize' not forced
            dynamic_cast<IntegrateNMLStepKernel&>(kernel.getImpl()).execute(*context, *this, currentPE, 4 ); //stepType 4 is quadratic minimizer   
            setProjVecChanged(false);
            
            std::cout << "Quadratic " << i << ", Current " << context->calcPotentialEnergy() << ", Old " << initialPE << ", lim " << minlim << std::endl;
            
            const double quadraticPE = context->calcPotentialEnergy();
            
            if(quadraticPE > initialPE){
                break;
            }
            
        }
        
        
    }
    
}
