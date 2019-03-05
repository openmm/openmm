
/* Portions copyright (c) 2006-2019 Stanford University and Simbios.
 * Contributors: Peter Eastman, Pande Group
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
#include "ReferenceCCMAAlgorithm.h"
#include "ReferenceDynamics.h"
#include "quern.h"
#include "openmm/OpenMMException.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ThreadPool.h"
#include <map>
#include <set>
#include <utility>

using namespace OpenMM;
using namespace std;

ReferenceCCMAAlgorithm::ReferenceCCMAAlgorithm(int numberOfAtoms,
                                               int numberOfConstraints,
                                               const vector<pair<int, int> >& atomIndices,
                                               const vector<double>& distance,
                                               vector<double>& masses,
                                               vector<AngleInfo>& angles,
                                               double elementCutoff) {
    _numberOfConstraints = numberOfConstraints;
    _elementCutoff = elementCutoff;
    _atomIndices = atomIndices;
    _distance = distance;

    _maximumNumberOfIterations = 150;
    _hasInitializedMasses = false;

    // work arrays

    if (_numberOfConstraints > 0) {
        _r_ij.resize(numberOfConstraints);
        _d_ij2 = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray(numberOfConstraints, NULL, 1, 0.0, "dij_2");
        _distanceTolerance = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray(numberOfConstraints, NULL, 1, 0.0, "distanceTolerance");
        _reducedMasses = SimTKOpenMMUtilities::allocateOneDRealOpenMMArray(numberOfConstraints, NULL, 1, 0.0, "reducedMasses");
    }
    if (numberOfConstraints > 0)
    {
        // Compute the constraint coupling matrix

        vector<vector<int> > atomAngles(numberOfAtoms);
        for (int i = 0; i < (int) angles.size(); i++)
            atomAngles[angles[i].atom2].push_back(i);
        vector<vector<pair<int, double> > > matrix(numberOfConstraints);
        vector<set<int> > atomConstraints(numberOfAtoms);
        for (int j = 0; j < numberOfConstraints; j++) {
            atomConstraints[_atomIndices[j].first].insert(j);
            atomConstraints[_atomIndices[j].second].insert(j);
        }
        for (int j = 0; j < numberOfConstraints; j++) {
            int atomj0 = _atomIndices[j].first;
            int atomj1 = _atomIndices[j].second;
            double invMass0 = 1/masses[atomj0];
            double invMass1 = 1/masses[atomj1];
            set<int> constraints = atomConstraints[atomj0];
            constraints.insert(atomConstraints[atomj1].begin(), atomConstraints[atomj1].end());
            for (int k : constraints) {
                if (j == k) {
                    matrix[j].push_back(pair<int, double>(j, 1.0));
                    continue;
                }
                double scale;
                int atomk0 = _atomIndices[k].first;
                int atomk1 = _atomIndices[k].second;
                int atoma, atomb, atomc;
                if (atomj0 == atomk0) {
                    atoma = atomj1;
                    atomb = atomj0;
                    atomc = atomk1;
                    scale = invMass0/(invMass0+invMass1);
                }
                else if (atomj1 == atomk1) {
                    atoma = atomj0;
                    atomb = atomj1;
                    atomc = atomk0;
                    scale = invMass1/(invMass0+invMass1);
                }
                else if (atomj0 == atomk1) {
                    atoma = atomj1;
                    atomb = atomj0;
                    atomc = atomk0;
                    scale = invMass0/(invMass0+invMass1);
                }
                else if (atomj1 == atomk0) {
                    atoma = atomj0;
                    atomb = atomj1;
                    atomc = atomk1;
                    scale = invMass1/(invMass0+invMass1);
                }
                else
                    continue; // These constraints are not connected.

                // Look for a third constraint forming a triangle with these two.

                bool foundConstraint = false;
                for (int other : atomConstraints[atoma]) {
                    if (_atomIndices[other].first == atomc || _atomIndices[other].second == atomc) {
                        double d1 = _distance[j];
                        double d2 = _distance[k];
                        double d3 = _distance[other];
                        matrix[j].push_back(pair<int, double>(k, scale*(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2)));
                        foundConstraint = true;
                        break;
                    }
                }
                if (!foundConstraint) {
                    // We didn't find one, so look for an angle force field term.

                    const vector<int>& angleCandidates = atomAngles[atomb];
                    for (int candidate : angleCandidates) {
                        const AngleInfo& angle = angles[candidate];
                        if ((angle.atom1 == atoma && angle.atom3 == atomc) || (angle.atom3 == atoma && angle.atom1 == atomc)) {
                            matrix[j].push_back(pair<int, double>(k, scale*cos(angle.angle)));
                            break;
                        }
                    }
                }
            }
        }

        // Invert it using QR.

        vector<int> matrixRowStart;
        vector<int> matrixColIndex;
        vector<double> matrixValue;
        for (int i = 0; i < numberOfConstraints; i++) {
            matrixRowStart.push_back(matrixValue.size());
            for (auto& element : matrix[i]) {
                matrixColIndex.push_back(element.first);
                matrixValue.push_back(element.second);
            }
        }
        matrixRowStart.push_back(matrixValue.size());
        int *qRowStart, *qColIndex, *rRowStart, *rColIndex;
        double *qValue, *rValue;
        QUERN_compute_qr(numberOfConstraints, numberOfConstraints, &matrixRowStart[0], &matrixColIndex[0], &matrixValue[0], NULL,
                &qRowStart, &qColIndex, &qValue, &rRowStart, &rColIndex, &rValue);
        vector<vector<pair<int, double> > > transposedMatrix(numberOfConstraints);
        _matrix.resize(numberOfConstraints);

        // Extract columns from the inverse matrix one at a time.  It is done in parallel,
        // since this can be very slow.
        
        ThreadPool threads;
        threads.execute([&] (ThreadPool& pool, int threadIndex) {
            vector<double> rhs(numberOfConstraints);
            for (int i = threadIndex; i < numberOfConstraints; i += pool.getNumThreads()) {
                // Extract column i of the inverse matrix.

                for (int j = 0; j < numberOfConstraints; j++)
                    rhs[j] = (i == j ? 1.0 : 0.0);
                QUERN_multiply_with_q_transpose(numberOfConstraints, qRowStart, qColIndex, qValue, &rhs[0]);
                QUERN_solve_with_r(numberOfConstraints, rRowStart, rColIndex, rValue, &rhs[0], &rhs[0]);
                for (int j = 0; j < numberOfConstraints; j++) {
                    double value = rhs[j]*distance[i]/distance[j];
                    if (fabs(value) > elementCutoff)
                        transposedMatrix[i].push_back(pair<int, double>(j, value));
                }
            }
        });
        threads.waitForThreads();

        // For purposes of thread safety we extracted the matrix in transposed form, so we need to transpose it again.

        for (int i = 0; i < numberOfConstraints; i++) {
            for (int j = 0; j < transposedMatrix[i].size(); j++) {
                pair<int, double> value = transposedMatrix[i][j];
                _matrix[value.first].push_back(make_pair(i, value.second));
            }
        }
        QUERN_free_result(qRowStart, qColIndex, qValue);
        QUERN_free_result(rRowStart, rColIndex, rValue);
    }
}

ReferenceCCMAAlgorithm::~ReferenceCCMAAlgorithm() {
    if (_numberOfConstraints > 0) {
        SimTKOpenMMUtilities::freeOneDRealOpenMMArray(_d_ij2, "d_ij2");
        SimTKOpenMMUtilities::freeOneDRealOpenMMArray(_distanceTolerance, "distanceTolerance");
        SimTKOpenMMUtilities::freeOneDRealOpenMMArray(_reducedMasses, "reducedMasses");
    }
}

int ReferenceCCMAAlgorithm::getNumberOfConstraints() const {
    return _numberOfConstraints;
}

int ReferenceCCMAAlgorithm::getMaximumNumberOfIterations() const {
    return _maximumNumberOfIterations;
}

void ReferenceCCMAAlgorithm::setMaximumNumberOfIterations(int maximumNumberOfIterations) {
    _maximumNumberOfIterations = maximumNumberOfIterations;
}

void ReferenceCCMAAlgorithm::apply(vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& atomCoordinatesP,
                                         vector<double>& inverseMasses, double tolerance) {
    applyConstraints(atomCoordinates, atomCoordinatesP, inverseMasses, false, tolerance);
}

void ReferenceCCMAAlgorithm::applyToVelocities(std::vector<OpenMM::Vec3>& atomCoordinates,
               std::vector<OpenMM::Vec3>& velocities, std::vector<double>& inverseMasses, double tolerance) {
    applyConstraints(atomCoordinates, velocities, inverseMasses, true, tolerance);
}

void ReferenceCCMAAlgorithm::applyConstraints(vector<Vec3>& atomCoordinates,
                                         vector<Vec3>& atomCoordinatesP,
                                         vector<double>& inverseMasses, bool constrainingVelocities, double tolerance) {
    // temp arrays

    vector<Vec3>& r_ij = _r_ij;
    double* d_ij2 = _d_ij2;
    double* reducedMasses = _reducedMasses;

    // calculate reduced masses on 1st pass

    if (!_hasInitializedMasses) {
        _hasInitializedMasses = true;
        for (int ii = 0; ii < _numberOfConstraints; ii++) {
           int atomI = _atomIndices[ii].first;
           int atomJ = _atomIndices[ii].second;
           reducedMasses[ii] = 0.5/(inverseMasses[atomI] + inverseMasses[atomJ]);
        }
    }

    // setup: r_ij for each (i,j) constraint

    for (int ii = 0; ii < _numberOfConstraints; ii++) {
        int atomI = _atomIndices[ii].first;
        int atomJ = _atomIndices[ii].second;
        r_ij[ii] = atomCoordinates[atomI] - atomCoordinates[atomJ];
        d_ij2[ii] = r_ij[ii].dot(r_ij[ii]);
    }
    double lowerTol = 1-2*tolerance+tolerance*tolerance;
    double upperTol = 1+2*tolerance+tolerance*tolerance;

    // main loop

    int iterations = 0;
    int numberConverged = 0;
    vector<double> constraintDelta(_numberOfConstraints);
    vector<double> tempDelta(_numberOfConstraints);
    while (iterations < getMaximumNumberOfIterations()) {
        numberConverged  = 0;
        for (int ii = 0; ii < _numberOfConstraints; ii++) {
            int atomI = _atomIndices[ii].first;
            int atomJ = _atomIndices[ii].second;
            Vec3 rp_ij = atomCoordinatesP[atomI] - atomCoordinatesP[atomJ];
            if (constrainingVelocities) {
                double rrpr = rp_ij.dot(r_ij[ii]);
                constraintDelta[ii] = -2*reducedMasses[ii]*rrpr/d_ij2[ii];
                if (fabs(constraintDelta[ii]) <= tolerance)
                    numberConverged++;
            }
            else {
                double rp2  = rp_ij.dot(rp_ij);
                double dist2 = _distance[ii]*_distance[ii];
                double diff = dist2 - rp2;
                constraintDelta[ii] = 0;
                double rrpr = DOT3(rp_ij, r_ij[ii]);
                constraintDelta[ii] = reducedMasses[ii]*diff/rrpr;
                if (rp2 >= lowerTol*dist2 && rp2 <= upperTol*dist2)
                    numberConverged++;
            }
        }
        if (numberConverged == _numberOfConstraints)
            break;
        iterations++;

        if (_matrix.size() > 0) {
            for (int i = 0; i < _numberOfConstraints; i++) {
                double sum = 0.0;
                for (auto& element : _matrix[i])
                    sum += element.second*constraintDelta[element.first];
                tempDelta[i] = sum;
            }
            constraintDelta = tempDelta;
        }
        for (int ii = 0; ii < _numberOfConstraints; ii++) {
            int atomI = _atomIndices[ii].first;
            int atomJ = _atomIndices[ii].second;
            Vec3 dr = r_ij[ii]*constraintDelta[ii];
            atomCoordinatesP[atomI] += dr*inverseMasses[atomI];
            atomCoordinatesP[atomJ] -= dr*inverseMasses[atomJ];
        }
    }
}

const vector<vector<pair<int, double> > >& ReferenceCCMAAlgorithm::getMatrix() const {
    return _matrix;
}
