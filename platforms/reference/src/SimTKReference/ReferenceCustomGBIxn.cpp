
/* Portions copyright (c) 2009 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string.h>
#include <sstream>

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomGBIxn.h"

using std::map;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

/**---------------------------------------------------------------------------------------

   ReferenceCustomGBIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomGBIxn::ReferenceCustomGBIxn(const vector<Lepton::ExpressionProgram>& valueExpressions,
                     const vector<vector<Lepton::ExpressionProgram> >& valueDerivExpressions,
                     const vector<string>& valueNames,
                     const vector<OpenMM::CustomGBForce::ComputationType>& valueTypes,
                     const vector<Lepton::ExpressionProgram>& energyExpressions,
                     const vector<vector<Lepton::ExpressionProgram> > energyDerivExpressions,
                     const vector<OpenMM::CustomGBForce::ComputationType>& energyTypes,
                     const vector<string>& parameterNames) :
            cutoff(false), periodic(false), valueExpressions(valueExpressions), valueDerivExpressions(valueDerivExpressions), valueNames(valueNames),
            valueTypes(valueTypes), energyExpressions(energyExpressions), energyDerivExpressions(energyDerivExpressions),
            energyTypes(energyTypes), paramNames(parameterNames) {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomGBIxn::ReferenceCustomGBIxn";

   // ---------------------------------------------------------------------------------------

    for (int i = 0; i < (int) paramNames.size(); i++) {
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << paramNames[i] << j;
            particleParamNames.push_back(name.str());
        }
    }
    for (int i = 0; i < (int) valueNames.size(); i++) {
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << valueNames[i] << j;
            particleValueNames.push_back(name.str());
        }
    }
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomGBIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomGBIxn::~ReferenceCustomGBIxn( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nReferenceCustomGBIxn::~ReferenceCustomGBIxn";

   // ---------------------------------------------------------------------------------------

}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

  int ReferenceCustomGBIxn::setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors ) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;

    return ReferenceForce::DefaultReturn;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param boxSize             the X, Y, and Z widths of the periodic box

     @return ReferenceForce::DefaultReturn

     --------------------------------------------------------------------------------------- */

  int ReferenceCustomGBIxn::setPeriodic( RealOpenMM* boxSize ) {

    assert(cutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    periodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
    return ReferenceForce::DefaultReturn;

  }


/**---------------------------------------------------------------------------------------

   Calculate the custom pair ixn

   @param numberOfAtoms    number of atoms
   @param atomCoordinates  atom coordinates
   @param atomParameters   atom parameters                             atomParameters[atomIndex][paramterIndex]
   @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                           exclusions[atomIndex][0] = number of exclusions
                           exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                           interacting w/ atom atomIndex
   @param fixedParameters  non atom parameters (not currently used)
   @param globalParameters the values of global parameters
   @param forces           force array (forces added)
   @param energyByAtom     atom energy
   @param totalEnergy      total energy

   @return ReferenceForce::DefaultReturn

   --------------------------------------------------------------------------------------- */

int ReferenceCustomGBIxn::calculateIxn(int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
                                           const vector<set<int> >& exclusions, map<string, double>& globalParameters, RealOpenMM** forces,
                                           RealOpenMM* energyByAtom, RealOpenMM* totalEnergy) const {
    int numValues = valueTypes.size();
    vector<vector<ComputedValue> > values(numValues);
    for (int valueIndex = 0; valueIndex < numValues; valueIndex++) {
        if (valueTypes[valueIndex] == OpenMM::CustomGBForce::SingleParticle)
            calculateSingleParticleValue(valueIndex, numberOfAtoms, values, globalParameters, atomParameters);
        else if (valueTypes[valueIndex] == OpenMM::CustomGBForce::ParticlePair)
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, values, globalParameters, exclusions, true);
        else
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, values, globalParameters, exclusions, false);
    }
    for (int termIndex = 0; termIndex < (int) energyExpressions.size(); termIndex++) {
        if (energyTypes[termIndex] == OpenMM::CustomGBForce::SingleParticle)
            calculateSingleParticleEnergyTerm(termIndex, numberOfAtoms, values, globalParameters, atomParameters, forces, totalEnergy);
        else if (energyTypes[termIndex] == OpenMM::CustomGBForce::ParticlePair)
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, values, globalParameters, exclusions, true, forces, totalEnergy);
        else
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, values, globalParameters, exclusions, false, forces, totalEnergy);
    }
   return ReferenceForce::DefaultReturn;
}

void ReferenceCustomGBIxn::calculateSingleParticleValue(int index, int numAtoms, vector<vector<ComputedValue> >& values,
        const map<string, double>& globalParameters, RealOpenMM** atomParameters) const {
    values[index].resize(numAtoms);
    map<string, double> variables = globalParameters;
    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < (int) paramNames.size(); j++)
            variables[paramNames[j]] = atomParameters[i][j];
        for (int j = 0; j < index; j++)
            variables[valueNames[j]] = values[j][i].value;
        values[index][i].value = (RealOpenMM) valueExpressions[index].evaluate(variables);
        for (int j = 0; j < index; j++) {
            RealOpenMM scale = (RealOpenMM) valueDerivExpressions[index][j].evaluate(variables);
            values[index][i].gradient[0] += scale*values[j][i].gradient[0];
            values[index][i].gradient[1] += scale*values[j][i].gradient[1];
            values[index][i].gradient[2] += scale*values[j][i].gradient[2];
        }
    }
}

void ReferenceCustomGBIxn::calculateParticlePairValue(int index, int numAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
        vector<vector<ComputedValue> >& values, const map<string, double>& globalParameters, const vector<set<int> >& exclusions, bool useExclusions) const {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        values[index][i].value = (RealOpenMM) 0.0;
    if (cutoff) {
        for (int i = 0; i < (int) neighborList->size(); i++) {
            OpenMM::AtomPair pair = (*neighborList)[i];
            if (useExclusions && exclusions[pair.first].find(pair.second) != exclusions[pair.first].end())
                continue;
            calculateOnePairValue(index, pair.first, pair.second, atomCoordinates, atomParameters, globalParameters, values);
            calculateOnePairValue(index, pair.second, pair.first, atomCoordinates, atomParameters, globalParameters, values);
        }
    }
    else {
        for (int i = 0; i < numAtoms; i++){
            for (int j = i+1; j < numAtoms; j++ ){
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairValue(index, i, j, atomCoordinates, atomParameters, globalParameters, values);
                calculateOnePairValue(index, j, i, atomCoordinates, atomParameters, globalParameters, values);
           }
        }
    }
}

void ReferenceCustomGBIxn::calculateOnePairValue(int index, int atom1, int atom2, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
        const map<string, double>& globalParameters, vector<vector<ComputedValue> >& values) const {
    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxSize, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    RealOpenMM r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;
    map<string, double> variables = globalParameters;
    for (int i = 0; i < paramNames.size(); i++) {
        variables[particleParamNames[i*2]] = atomParameters[atom1][i];
        variables[particleParamNames[i*2+1]] = atomParameters[atom2][i];
    }
    variables["r"] = r;
    for (int i = 0; i < index; i++) {
        variables[particleValueNames[i*2]] = values[i][atom1].value;
        variables[particleValueNames[i*2+1]] = values[i][atom2].value;
    }
    values[index][atom1].value += (RealOpenMM) valueExpressions[index].evaluate(variables);
    RealOpenMM dVdR = (RealOpenMM) valueDerivExpressions[index][0].evaluate(variables);
    RealOpenMM rinv = 1/r;
    for (int i = 0; i < 3; i++)
       values[index][atom1].gradient[i] += dVdR*deltaR[i]*rinv;
    for (int i = 0; i < index; i++) {
        RealOpenMM scale = (RealOpenMM) valueDerivExpressions[index][i+1].evaluate(variables);
        values[index][atom1].gradient[0] += scale*values[i][atom1].gradient[0];
        values[index][atom1].gradient[1] += scale*values[i][atom1].gradient[1];
        values[index][atom1].gradient[2] += scale*values[i][atom1].gradient[2];
    }
}

void ReferenceCustomGBIxn::calculateSingleParticleEnergyTerm(int index, int numAtoms, const vector<vector<ComputedValue> >& values,
        const map<string, double>& globalParameters, RealOpenMM** atomParameters, RealOpenMM** forces, RealOpenMM* totalEnergy) const {
    map<string, double> variables = globalParameters;
    for (int i = 0; i < numAtoms; i++) {
        for (int j = 0; j < (int) paramNames.size(); j++)
            variables[paramNames[j]] = atomParameters[i][j];
        for (int j = 0; j < (int) valueNames.size(); j++)
            variables[valueNames[j]] = values[j][i].value;
        if (totalEnergy != NULL)
            *totalEnergy += (RealOpenMM) energyExpressions[index].evaluate(variables);
        for (int j = 0; j < (int) valueNames.size(); j++) {
            RealOpenMM scale = (RealOpenMM) energyDerivExpressions[index][j].evaluate(variables);
            forces[i][0] -= scale*values[j][i].gradient[0];
            forces[i][1] -= scale*values[j][i].gradient[1];
            forces[i][2] -= scale*values[j][i].gradient[2];
        }
    }
}

void ReferenceCustomGBIxn::calculateParticlePairEnergyTerm(int index, int numAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
        const vector<vector<ComputedValue> >& values, const map<string, double>& globalParameters, const vector<set<int> >& exclusions, bool useExclusions,
        RealOpenMM** forces, RealOpenMM* totalEnergy) const {
    if (cutoff) {
        for (int i = 0; i < (int) neighborList->size(); i++) {
            OpenMM::AtomPair pair = (*neighborList)[i];
            if (useExclusions && exclusions[pair.first].find(pair.second) != exclusions[pair.first].end())
                continue;
            calculateOnePairEnergyTerm(index, pair.first, pair.second, atomCoordinates, atomParameters, globalParameters, values, forces, totalEnergy);
        }
    }
    else {
        for (int i = 0; i < numAtoms; i++){
            for (int j = i+1; j < numAtoms; j++ ){
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairEnergyTerm(index, i, j, atomCoordinates, atomParameters, globalParameters, values, forces, totalEnergy);
           }
        }
    }
}

void ReferenceCustomGBIxn::calculateOnePairEnergyTerm(int index, int atom1, int atom2, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
        const map<string, double>& globalParameters, const vector<vector<ComputedValue> >& values, RealOpenMM** forces, RealOpenMM* totalEnergy) const {
    RealOpenMM deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxSize, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    RealOpenMM r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;
    map<string, double> variables = globalParameters;
    for (int i = 0; i < paramNames.size(); i++) {
        variables[particleParamNames[i*2]] = atomParameters[atom1][i];
        variables[particleParamNames[i*2+1]] = atomParameters[atom2][i];
    }
    variables["r"] = r;
    for (int i = 0; i < (int) valueNames.size(); i++) {
        variables[particleValueNames[i*2]] = values[i][atom1].value;
        variables[particleValueNames[i*2+1]] = values[i][atom2].value;
    }
    if (totalEnergy != NULL)
        *totalEnergy += (RealOpenMM) energyExpressions[index].evaluate(variables);
    RealOpenMM dEdR = (RealOpenMM) energyDerivExpressions[index][0].evaluate(variables);
    dEdR *= 1/r;
    for (int i = 0; i < 3; i++) {
       forces[atom1][i] -= dEdR*deltaR[i];
       forces[atom2][i] += dEdR*deltaR[i];
    }
    for (int i = 0; i < (int) valueNames.size(); i++) {
        RealOpenMM scale1 = (RealOpenMM) energyDerivExpressions[index][2*i+1].evaluate(variables);
        RealOpenMM scale2 = (RealOpenMM) energyDerivExpressions[index][2*i+2].evaluate(variables);
        forces[atom1][0] -= scale1*values[i][atom1].gradient[0];
        forces[atom1][1] -= scale1*values[i][atom1].gradient[1];
        forces[atom1][2] -= scale1*values[i][atom1].gradient[2];
        forces[atom2][0] += scale2*values[i][atom2].gradient[0];
        forces[atom2][1] += scale2*values[i][atom2].gradient[1];
        forces[atom2][2] += scale2*values[i][atom2].gradient[2];
    }
}
