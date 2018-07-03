
/* Portions copyright (c) 2009-2018 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "ReferenceCustomGBIxn.h"

using std::map;
using std::set;
using std::string;
using std::stringstream;
using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceCustomGBIxn constructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomGBIxn::ReferenceCustomGBIxn(const vector<Lepton::CompiledExpression>& valueExpressions,
                     const vector<vector<Lepton::CompiledExpression> > valueDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> > valueGradientExpressions,
                     const vector<vector<Lepton::CompiledExpression> > valueParamDerivExpressions,
                     const vector<string>& valueNames,
                     const vector<OpenMM::CustomGBForce::ComputationType>& valueTypes,
                     const vector<Lepton::CompiledExpression>& energyExpressions,
                     const vector<vector<Lepton::CompiledExpression> > energyDerivExpressions,
                     const vector<vector<Lepton::CompiledExpression> > energyGradientExpressions,
                     const vector<vector<Lepton::CompiledExpression> > energyParamDerivExpressions,
                     const vector<OpenMM::CustomGBForce::ComputationType>& energyTypes,
                     const vector<string>& parameterNames) :
            cutoff(false), periodic(false), valueExpressions(valueExpressions), valueDerivExpressions(valueDerivExpressions), valueGradientExpressions(valueGradientExpressions), valueParamDerivExpressions(valueParamDerivExpressions),
            valueTypes(valueTypes), energyExpressions(energyExpressions), energyDerivExpressions(energyDerivExpressions), energyGradientExpressions(energyGradientExpressions), energyParamDerivExpressions(energyParamDerivExpressions),
            energyTypes(energyTypes) {

    for (int i = 0; i < this->valueExpressions.size(); i++)
        expressionSet.registerExpression(this->valueExpressions[i]);
    for (int i = 0; i < this->valueDerivExpressions.size(); i++)
        for (int j = 0; j < this->valueDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->valueDerivExpressions[i][j]);
    for (int i = 0; i < this->valueGradientExpressions.size(); i++)
        for (int j = 0; j < this->valueGradientExpressions[i].size(); j++)
            expressionSet.registerExpression(this->valueGradientExpressions[i][j]);
    for (int i = 0; i < this->valueParamDerivExpressions.size(); i++)
        for (int j = 0; j < this->valueParamDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->valueParamDerivExpressions[i][j]);
    for (int i = 0; i < this->energyExpressions.size(); i++)
        expressionSet.registerExpression(this->energyExpressions[i]);
    for (int i = 0; i < this->energyDerivExpressions.size(); i++)
        for (int j = 0; j < this->energyDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->energyDerivExpressions[i][j]);
    for (int i = 0; i < this->energyGradientExpressions.size(); i++)
        for (int j = 0; j < this->energyGradientExpressions[i].size(); j++)
            expressionSet.registerExpression(this->energyGradientExpressions[i][j]);
    for (int i = 0; i < this->energyParamDerivExpressions.size(); i++)
        for (int j = 0; j < this->energyParamDerivExpressions[i].size(); j++)
            expressionSet.registerExpression(this->energyParamDerivExpressions[i][j]);
    rIndex = expressionSet.getVariableIndex("r");
    xIndex = expressionSet.getVariableIndex("x");
    yIndex = expressionSet.getVariableIndex("y");
    zIndex = expressionSet.getVariableIndex("z");
    for (auto& param : parameterNames) {
        paramIndex.push_back(expressionSet.getVariableIndex(param));
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << param << j;
            particleParamIndex.push_back(expressionSet.getVariableIndex(name.str()));
        }
    }
    for (auto& value : valueNames) {
        valueIndex.push_back(expressionSet.getVariableIndex(value));
        for (int j = 1; j < 3; j++) {
            stringstream name;
            name << value << j;
            particleValueIndex.push_back(expressionSet.getVariableIndex(name.str()));
        }
    }
}

/**---------------------------------------------------------------------------------------

   ReferenceCustomGBIxn destructor

   --------------------------------------------------------------------------------------- */

ReferenceCustomGBIxn::~ReferenceCustomGBIxn() {
}

  /**---------------------------------------------------------------------------------------

     Set the force to use a cutoff.

     @param distance            the cutoff distance
     @param neighbors           the neighbor list to use

     --------------------------------------------------------------------------------------- */

  void ReferenceCustomGBIxn::setUseCutoff(double distance, const OpenMM::NeighborList& neighbors) {

    cutoff = true;
    cutoffDistance = distance;
    neighborList = &neighbors;
  }

  /**---------------------------------------------------------------------------------------

     Set the force to use periodic boundary conditions.  This requires that a cutoff has
     also been set, and the smallest side of the periodic box is at least twice the cutoff
     distance.

     @param vectors    the vectors defining the periodic box

     --------------------------------------------------------------------------------------- */

  void ReferenceCustomGBIxn::setPeriodic(Vec3* vectors) {

    if (cutoff) {
        assert(vectors[0][0] >= 2.0*cutoffDistance);
        assert(vectors[1][1] >= 2.0*cutoffDistance);
        assert(vectors[2][2] >= 2.0*cutoffDistance);
    }
    periodic = true;
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];
  }

void ReferenceCustomGBIxn::calculateIxn(int numberOfAtoms, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
                                           const vector<set<int> >& exclusions, map<string, double>& globalParameters, vector<Vec3>& forces,
                                           double* totalEnergy, double* energyParamDerivs) {
    for (auto& param : globalParameters)
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first), param.second);
    
    // Initialize arrays for storing values.
    
    int numValues = valueTypes.size();
    int numDerivs = valueParamDerivExpressions[0].size();
    values.resize(numValues);
    dEdV.resize(numValues, vector<double>(numberOfAtoms, 0.0));
    dValuedParam.resize(numValues);
    for (int i = 0; i < numValues; i++)
        dValuedParam[i].resize(numDerivs, vector<double>(numberOfAtoms, 0.0));

    // First calculate the computed values.

    for (int valueIndex = 0; valueIndex < numValues; valueIndex++) {
        if (valueTypes[valueIndex] == OpenMM::CustomGBForce::SingleParticle)
            calculateSingleParticleValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters);
        else if (valueTypes[valueIndex] == OpenMM::CustomGBForce::ParticlePair)
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, exclusions, true);
        else
            calculateParticlePairValue(valueIndex, numberOfAtoms, atomCoordinates, atomParameters, exclusions, false);
    }

    // Now calculate the energy and its derivates.

    for (int termIndex = 0; termIndex < (int) energyExpressions.size(); termIndex++) {
        if (energyTypes[termIndex] == OpenMM::CustomGBForce::SingleParticle)
            calculateSingleParticleEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, forces, totalEnergy, energyParamDerivs);
        else if (energyTypes[termIndex] == OpenMM::CustomGBForce::ParticlePair)
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, exclusions, true, forces, totalEnergy, energyParamDerivs);
        else
            calculateParticlePairEnergyTerm(termIndex, numberOfAtoms, atomCoordinates, atomParameters, exclusions, false, forces, totalEnergy, energyParamDerivs);
    }

    // Apply the chain rule to evaluate forces.

    calculateChainRuleForces(numberOfAtoms, atomCoordinates, atomParameters, exclusions, forces, energyParamDerivs);
}

void ReferenceCustomGBIxn::calculateSingleParticleValue(int index, int numAtoms, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xIndex, atomCoordinates[i][0]);
        expressionSet.setVariable(yIndex, atomCoordinates[i][1]);
        expressionSet.setVariable(zIndex, atomCoordinates[i][2]);
        for (int j = 0; j < (int) paramIndex.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < index; j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        values[index][i] = valueExpressions[index].evaluate();

        // Calculate derivatives with respect to parameters.

        for (int j = 0; j < valueParamDerivExpressions[index].size(); j++)
            dValuedParam[index][j][i] += valueParamDerivExpressions[index][j].evaluate();
        for (int j = 0; j < index; j++) {
            double dVdV = valueDerivExpressions[index][j].evaluate();
            for (int k = 0; k < valueParamDerivExpressions[index].size(); k++)
                dValuedParam[index][k][i] += dVdV*dValuedParam[j][k][i];
        }
    }
}

void ReferenceCustomGBIxn::calculateParticlePairValue(int index, int numAtoms, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
        const vector<set<int> >& exclusions, bool useExclusions) {
    values[index].resize(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        values[index][i] = 0.0;
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (auto& pair : *neighborList) {
            if (useExclusions && exclusions[pair.first].find(pair.second) != exclusions[pair.first].end())
                continue;
            calculateOnePairValue(index, pair.first, pair.second, atomCoordinates, atomParameters);
            calculateOnePairValue(index, pair.second, pair.first, atomCoordinates, atomParameters);
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairValue(index, i, j, atomCoordinates, atomParameters);
                calculateOnePairValue(index, j, i, atomCoordinates, atomParameters);
           }
        }
    }
}

void ReferenceCustomGBIxn::calculateOnePairValue(int index, int atom1, int atom2, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters) {
    double deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    double r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;
    for (int i = 0; i < (int) paramIndex.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rIndex, r);
    for (int i = 0; i < index; i++) {
        expressionSet.setVariable(particleValueIndex[i*2], values[i][atom1]);
        expressionSet.setVariable(particleValueIndex[i*2+1], values[i][atom2]);
    }
    values[index][atom1] += valueExpressions[index].evaluate();
    
    // Calculate derivatives with respect to parameters.
    
    for (int i = 0; i < valueParamDerivExpressions[index].size(); i++)
        dValuedParam[index][i][atom1] += valueParamDerivExpressions[index][i].evaluate();
}

void ReferenceCustomGBIxn::calculateSingleParticleEnergyTerm(int index, int numAtoms, vector<Vec3>& atomCoordinates,
        vector<vector<double> >& atomParameters, vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xIndex, atomCoordinates[i][0]);
        expressionSet.setVariable(yIndex, atomCoordinates[i][1]);
        expressionSet.setVariable(zIndex, atomCoordinates[i][2]);
        for (int j = 0; j < (int) paramIndex.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 0; j < valueIndex.size(); j++)
            expressionSet.setVariable(valueIndex[j], values[j][i]);
        
        // Compute energy and force.
        
        if (totalEnergy != NULL)
            *totalEnergy += energyExpressions[index].evaluate();
        for (int j = 0; j < (int) valueIndex.size(); j++)
            dEdV[j][i] += energyDerivExpressions[index][j].evaluate();
        forces[i][0] -= energyGradientExpressions[index][0].evaluate();
        forces[i][1] -= energyGradientExpressions[index][1].evaluate();
        forces[i][2] -= energyGradientExpressions[index][2].evaluate();
        
        // Compute derivatives with respect to parameters.
        
        for (int k = 0; k < energyParamDerivExpressions[index].size(); k++)
            energyParamDerivs[k] += energyParamDerivExpressions[index][k].evaluate();
    }
}

void ReferenceCustomGBIxn::calculateParticlePairEnergyTerm(int index, int numAtoms, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
        const vector<set<int> >& exclusions, bool useExclusions, vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (auto& pair : *neighborList) {
            if (useExclusions && exclusions[pair.first].find(pair.second) != exclusions[pair.first].end())
                continue;
            calculateOnePairEnergyTerm(index, pair.first, pair.second, atomCoordinates, atomParameters, forces, totalEnergy, energyParamDerivs);
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                if (useExclusions && exclusions[i].find(j) != exclusions[i].end())
                    continue;
                calculateOnePairEnergyTerm(index, i, j, atomCoordinates, atomParameters, forces, totalEnergy, energyParamDerivs);
           }
        }
    }
}

void ReferenceCustomGBIxn::calculateOnePairEnergyTerm(int index, int atom1, int atom2, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
        vector<Vec3>& forces, double* totalEnergy, double* energyParamDerivs) {
    // Compute the displacement.

    double deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    double r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramIndex.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rIndex, r);
    for (int i = 0; i < (int) valueIndex.size(); i++) {
        expressionSet.setVariable(particleValueIndex[i*2], values[i][atom1]);
        expressionSet.setVariable(particleValueIndex[i*2+1], values[i][atom2]);
    }

    // Evaluate the energy and its derivatives.

    if (totalEnergy != NULL)
        *totalEnergy += energyExpressions[index].evaluate();
    double dEdR = energyDerivExpressions[index][0].evaluate();
    dEdR *= 1/r;
    for (int i = 0; i < 3; i++) {
       forces[atom1][i] -= dEdR*deltaR[i];
       forces[atom2][i] += dEdR*deltaR[i];
    }
    for (int i = 0; i < (int) valueIndex.size(); i++) {
        dEdV[i][atom1] += energyDerivExpressions[index][2*i+1].evaluate();
        dEdV[i][atom2] += energyDerivExpressions[index][2*i+2].evaluate();
    }
        
    // Compute derivatives with respect to parameters.

    for (int i = 0; i < energyParamDerivExpressions[index].size(); i++)
        energyParamDerivs[i] += energyParamDerivExpressions[index][i].evaluate();
}

void ReferenceCustomGBIxn::calculateChainRuleForces(int numAtoms, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
        const vector<set<int> >& exclusions, vector<Vec3>& forces, double* energyParamDerivs) {
    if (cutoff) {
        // Loop over all pairs in the neighbor list.

        for (auto& pair : *neighborList) {
            bool isExcluded = (exclusions[pair.first].find(pair.second) != exclusions[pair.first].end());
            calculateOnePairChainRule(pair.first, pair.second, atomCoordinates, atomParameters, forces, isExcluded);
            calculateOnePairChainRule(pair.second, pair.first, atomCoordinates, atomParameters, forces, isExcluded);
        }
    }
    else {
        // Perform an O(N^2) loop over all atom pairs.

        for (int i = 0; i < numAtoms; i++) {
            for (int j = i+1; j < numAtoms; j++) {
                bool isExcluded = (exclusions[i].find(j) != exclusions[i].end());
                calculateOnePairChainRule(i, j, atomCoordinates, atomParameters, forces, isExcluded);
                calculateOnePairChainRule(j, i, atomCoordinates, atomParameters, forces, isExcluded);
           }
        }
    }

    // Compute chain rule terms for computed values that depend explicitly on particle coordinates.

    for (int i = 0; i < numAtoms; i++) {
        expressionSet.setVariable(xIndex, atomCoordinates[i][0]);
        expressionSet.setVariable(yIndex, atomCoordinates[i][1]);
        expressionSet.setVariable(zIndex, atomCoordinates[i][2]);
        vector<double> dVdX(valueDerivExpressions.size(), 0.0);
        vector<double> dVdY(valueDerivExpressions.size(), 0.0);
        vector<double> dVdZ(valueDerivExpressions.size(), 0.0);
        for (int j = 0; j < (int) paramIndex.size(); j++)
            expressionSet.setVariable(paramIndex[j], atomParameters[i][j]);
        for (int j = 1; j < (int) valueIndex.size(); j++) {
            expressionSet.setVariable(valueIndex[j-1], values[j-1][i]);
            for (int k = 1; k < j; k++) {
                double dVdV = valueDerivExpressions[j][k].evaluate();
                dVdX[j] += dVdV*dVdX[k];
                dVdY[j] += dVdV*dVdY[k];
                dVdZ[j] += dVdV*dVdZ[k];
            }
            dVdX[j] += valueGradientExpressions[j][0].evaluate();
            dVdY[j] += valueGradientExpressions[j][1].evaluate();
            dVdZ[j] += valueGradientExpressions[j][2].evaluate();
            forces[i][0] -= dEdV[j][i]*dVdX[j];
            forces[i][1] -= dEdV[j][i]*dVdY[j];
            forces[i][2] -= dEdV[j][i]*dVdZ[j];
        }
    }
        
    // Compute chain rule terms for derivatives with respect to parameters.

    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < (int) valueIndex.size(); j++)
            for (int k = 0; k < dValuedParam[j].size(); k++)
                energyParamDerivs[k] += dEdV[j][i]*dValuedParam[j][k][i];
}

void ReferenceCustomGBIxn::calculateOnePairChainRule(int atom1, int atom2, vector<Vec3>& atomCoordinates, vector<vector<double> >& atomParameters,
        vector<Vec3>& forces, bool isExcluded) {
    // Compute the displacement.

    double deltaR[ReferenceForce::LastDeltaRIndex];
    if (periodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom2], atomCoordinates[atom1], periodicBoxVectors, deltaR);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom2], atomCoordinates[atom1], deltaR);
    double r = deltaR[ReferenceForce::RIndex];
    if (cutoff && r >= cutoffDistance)
        return;

    // Record variables for evaluating expressions.

    for (int i = 0; i < (int) paramIndex.size(); i++) {
        expressionSet.setVariable(particleParamIndex[i*2], atomParameters[atom1][i]);
        expressionSet.setVariable(particleParamIndex[i*2+1], atomParameters[atom2][i]);
    }
    expressionSet.setVariable(rIndex, r);
    expressionSet.setVariable(particleValueIndex[0], values[0][atom1]);
    expressionSet.setVariable(particleValueIndex[1], values[0][atom2]);

    // Evaluate the derivative of each parameter with respect to position and apply forces.

    double rinv = 1/r;
    deltaR[0] *= rinv;
    deltaR[1] *= rinv;
    deltaR[2] *= rinv;
    vector<double> dVdR1(valueDerivExpressions.size(), 0.0);
    vector<double> dVdR2(valueDerivExpressions.size(), 0.0);
    if (!isExcluded || valueTypes[0] != OpenMM::CustomGBForce::ParticlePair) {
        dVdR1[0] = valueDerivExpressions[0][0].evaluate();
        dVdR2[0] = -dVdR1[0];
        for (int i = 0; i < 3; i++) {
            forces[atom1][i] -= dEdV[0][atom1]*dVdR1[0]*deltaR[i];
            forces[atom2][i] -= dEdV[0][atom1]*dVdR2[0]*deltaR[i];
        }
    }
    for (int i = 0; i < (int) paramIndex.size(); i++)
        expressionSet.setVariable(paramIndex[i], atomParameters[atom1][i]);
    expressionSet.setVariable(valueIndex[0], values[0][atom1]);
    for (int i = 1; i < (int) valueIndex.size(); i++) {
        expressionSet.setVariable(valueIndex[i], values[i][atom1]);
        expressionSet.setVariable(xIndex, atomCoordinates[atom1][0]);
        expressionSet.setVariable(yIndex, atomCoordinates[atom1][1]);
        expressionSet.setVariable(zIndex, atomCoordinates[atom1][2]);
        for (int j = 0; j < i; j++) {
            double dVdV = valueDerivExpressions[i][j].evaluate();
            dVdR1[i] += dVdV*dVdR1[j];
            dVdR2[i] += dVdV*dVdR2[j];
        }
        for (int k = 0; k < 3; k++) {
            forces[atom1][k] -= dEdV[i][atom1]*dVdR1[i]*deltaR[k];
            forces[atom2][k] -= dEdV[i][atom1]*dVdR2[i]*deltaR[k];
        }
    }
}
