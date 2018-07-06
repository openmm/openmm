/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2018 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "CpuBondForce.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace std;

CpuBondForce::CpuBondForce() {
}

void CpuBondForce::initialize(int numAtoms, int numBonds, int numAtomsPerBond, vector<vector<int> >& bondAtoms, ThreadPool& threads) {
    this->numBonds = numBonds;
    this->numAtomsPerBond = numAtomsPerBond;
    this->bondAtoms = &bondAtoms[0];
    this->threads = &threads;
    int numThreads = threads.getNumThreads();
    int targetBondsPerThread = numBonds/numThreads;
    
    // Record the bonds that include each atom.
    
    vector<set<int> > atomBonds(numAtoms);
    for (int bond = 0; bond < numBonds; bond++) {
        for (int i = 0; i < numAtomsPerBond; i++)
            atomBonds[bondAtoms[bond][i]].insert(bond);
    }
    
    // Divide bonds into groups.
    
    vector<int> atomThread(numAtoms, -1);
    vector<int> bondThread(numBonds, -1);
    threadBonds.resize(numThreads);
    int numProcessed = 0;
    int thread = 0;
    list<int> candidateBonds;
    while (thread < numThreads) {
        // Find the next unassigned bond.
        
        while (numProcessed < numBonds && bondThread[numProcessed] != -1)
            numProcessed++;
        if (numProcessed == numBonds)
            break; // We've gone through the whole list of bonds.
        
        // See if this bond can be assigned to this thread.
        
        if (!canAssignBond(numProcessed, thread, atomThread)) {
            numProcessed++;
            continue;
        }
        
        // Assign this bond to the thread.
        
        assignBond(numProcessed++, thread, atomThread, bondThread, atomBonds, candidateBonds);
        
        // Assign additional bonds that have been identified as involving atoms assigned to this thread.
        
        while (!candidateBonds.empty() && threadBonds[thread].size() < targetBondsPerThread) {
            int bond = *candidateBonds.begin();
            if (bondThread[bond] == -1 && canAssignBond(bond, thread, atomThread))
                assignBond(bond, thread, atomThread, bondThread, atomBonds, candidateBonds);
            candidateBonds.pop_front();
        }
        
        // If we have assigned enough bonds to this thread, move on to the next one.
        
        if (threadBonds[thread].size() >= targetBondsPerThread) {
            candidateBonds.clear();
            thread++;
        }
    }
    
    // Look through the remaining bonds and see whether any of them can be assigned.
    
    candidateBonds.clear();
    for (int bond = 0; bond < numBonds; bond++) {
        if (bondThread[bond] == -1) {
            // See whether this bond can be assigned to a thread.
            
            bool canAssign = true;
            int assignment = -1;
            for (int i = 0; i < numAtomsPerBond; i++) {
                int thread = atomThread[bondAtoms[bond][i]];
                if (thread == -1 || thread == assignment)
                    continue;
                if (assignment == -1)
                    assignment = thread;
                else {
                    canAssign = false;
                    break;
                }
            }
            if (canAssign) {
                // Assign this bond to a thread.
                
                if (assignment == -1)
                    assignment = numThreads-1;
                assignBond(bond, assignment, atomThread, bondThread, atomBonds, candidateBonds);
            }
            else {
                // Add it to the list of "extra" bonds.
                
                extraBonds.push_back(bond);
            }
        }
    }
}

bool CpuBondForce::canAssignBond(int bond, int thread, vector<int>& atomThread) {
    for (int i = 0; i < numAtomsPerBond; i++) {
        int atom = bondAtoms[bond][i];
        if (atomThread[atom] != -1 && atomThread[atom] != thread)
            return false;
    }
    return true;
}

void CpuBondForce::assignBond(int bond, int thread, vector<int>& atomThread, vector<int>& bondThread, vector<set<int> >& atomBonds, list<int>& candidateBonds) {
    // Assign the bond to a thread.
    
    bondThread[bond] = thread;
    threadBonds[thread].push_back(bond);
    
    // Mark every atom in this bond as also belonging to the thread, and add all of their
    // bonds to the list of candidates.
    
    for (int i = 0; i < numAtomsPerBond; i++) {
        int& atom = atomThread[bondAtoms[bond][i]];
        if (atom == thread)
            continue;
        if (atom != -1)
            throw OpenMMException("CpuBondForce: Internal error: atoms assigned to threads incorrectly");
        atom = thread;
        for (int bond : atomBonds[atom])
            candidateBonds.push_back(bond);
    }
}

void CpuBondForce::calculateForce(vector<Vec3>& atomCoordinates, vector<vector<double> >& parameters, vector<Vec3>& forces, 
        double* totalEnergy, ReferenceBondIxn& referenceBondIxn) {
    // Have the worker threads compute their forces.
    
    vector<double> threadEnergy(threads->getNumThreads(), 0);
    threads->execute([&] (ThreadPool& threads, int threadIndex) {
        double* energy = (totalEnergy == NULL ? NULL : &threadEnergy[threadIndex]);
        threadComputeForce(threads, threadIndex, atomCoordinates, parameters, forces, energy, referenceBondIxn);
    });
    threads->waitForThreads();
    
    // Compute any "extra" bonds.
    
    for (int i = 0; i < extraBonds.size(); i++) {
        int bond = extraBonds[i];
        referenceBondIxn.calculateBondIxn(bondAtoms[bond], atomCoordinates, parameters[bond], forces, totalEnergy, NULL);
    }

    // Compute the total energy.
    
    if (totalEnergy != NULL)
        for (int i = 0; i < threads->getNumThreads(); i++)
            *totalEnergy += threadEnergy[i];
}

void CpuBondForce::threadComputeForce(ThreadPool& threads, int threadIndex, vector<Vec3>& atomCoordinates, vector<vector<double> >& parameters, vector<Vec3>& forces, 
            double* totalEnergy, ReferenceBondIxn& referenceBondIxn) {
    vector<int>& bonds = threadBonds[threadIndex];
    int numBonds = bonds.size();
    for (int i = 0; i < numBonds; i++) {
        int bond = bonds[i];
        referenceBondIxn.calculateBondIxn(bondAtoms[bond], atomCoordinates, parameters[bond], forces, totalEnergy, NULL);
    }
}