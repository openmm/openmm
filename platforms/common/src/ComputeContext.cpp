/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2024 Stanford University and the Authors.      *
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

#include "openmm/common/ComputeContext.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "hilbert.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;

const int ComputeContext::ThreadBlockSize = 64;
const int ComputeContext::TileSize = 32;

ComputeContext::ComputeContext(const System& system) : system(system), time(0.0), stepCount(0), computeForceCount(0), stepsSinceReorder(99999),
        forceNextReorder(false), atomsWereReordered(false), forcesValid(false), thread(NULL) {
    thread = new WorkThread();
}

ComputeContext::~ComputeContext() {
    if (thread != NULL)
        delete thread;
}

void ComputeContext::addForce(ComputeForceInfo* force) {
    forces.push_back(force);
}

void ComputeContext::setAtomIndex(std::vector<int>& index){
    atomIndex = index;
    getAtomIndexArray().upload(atomIndex);
    for (auto listener : reorderListeners)
        listener->execute();
}

string ComputeContext::replaceStrings(const string& input, const std::map<std::string, std::string>& replacements) const {
    static set<char> symbolChars;
    if (symbolChars.size() == 0) {
        symbolChars.insert('_');
        for (char c = 'a'; c <= 'z'; c++)
            symbolChars.insert(c);
        for (char c = 'A'; c <= 'Z'; c++)
            symbolChars.insert(c);
        for (char c = '0'; c <= '9'; c++)
            symbolChars.insert(c);
    }
    string result = input;
    for (auto& pair : replacements) {
        int index = 0;
        int size = pair.first.size();
        do {
            index = result.find(pair.first, index);
            if (index != result.npos) {
                if ((index == 0 || symbolChars.find(result[index-1]) == symbolChars.end()) && (index == result.size()-size || symbolChars.find(result[index+size]) == symbolChars.end())) {
                    // We have found a complete symbol, not part of a longer symbol.

                    // Do not allow to replace a symbol contained in single-line comments with a multi-line content
                    // because only the first line will be commented
                    // (the check is used to prevent incorrect commenting during development).
                    if (pair.second.find('\n') != pair.second.npos) {
                        int prevIndex = index;
                        while (prevIndex > 1 && result[prevIndex] != '\n') {
                            if (result[prevIndex] == '/' && result[prevIndex - 1] == '/') {
                                throw OpenMMException("Symbol " + pair.first + " is contained in a single-line comment");
                            }
                            prevIndex--;
                        }
                    }

                    result.replace(index, size, pair.second);
                    index += pair.second.size();
                }
                else
                    index++;
            }
        } while (index != result.npos);
    }
    return result;
}

string ComputeContext::doubleToString(double value, bool mixedIsDouble) const {
    stringstream s;
    bool useDouble = (getUseDoublePrecision() || (mixedIsDouble && getUseMixedPrecision()));
    s.precision(useDouble ? 16 : 8);
    s << scientific << value;
    if (!useDouble)
        s << "f";
    return s.str();
}

string ComputeContext::intToString(int value) const {
    stringstream s;
    s << value;
    return s.str();
}

/**
 * This class ensures that atom reordering doesn't break virtual sites.
 */
class ComputeContext::VirtualSiteInfo : public ComputeForceInfo {
public:
    VirtualSiteInfo(const System& system) {
        for (int i = 0; i < system.getNumParticles(); i++) {
            if (system.isVirtualSite(i)) {
                const VirtualSite& vsite = system.getVirtualSite(i);
                siteTypes.push_back(&typeid(vsite));
                vector<int> particles;
                particles.push_back(i);
                for (int j = 0; j < vsite.getNumParticles(); j++)
                    particles.push_back(vsite.getParticle(j));
                siteParticles.push_back(particles);
                vector<double> weights;
                if (dynamic_cast<const TwoParticleAverageSite*>(&vsite) != NULL) {
                    // A two particle average.

                    const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(vsite);
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                }
                else if (dynamic_cast<const ThreeParticleAverageSite*>(&vsite) != NULL) {
                    // A three particle average.

                    const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(vsite);
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                    weights.push_back(site.getWeight(2));
                }
                else if (dynamic_cast<const OutOfPlaneSite*>(&vsite) != NULL) {
                    // An out of plane site.

                    const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(vsite);
                    weights.push_back(site.getWeight12());
                    weights.push_back(site.getWeight13());
                    weights.push_back(site.getWeightCross());
                }
                siteWeights.push_back(weights);
            }
        }
    }
    int getNumParticleGroups() {
        return siteTypes.size();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        particles = siteParticles[index];
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (siteTypes[group1] != siteTypes[group2])
            return false;
        int numParticles = siteWeights[group1].size();
        if (siteWeights[group2].size() != numParticles)
            return false;
        for (int i = 0; i < numParticles; i++)
            if (siteWeights[group1][i] != siteWeights[group2][i])
                return false;
        return true;
    }
private:
    vector<const type_info*> siteTypes;
    vector<vector<int> > siteParticles;
    vector<vector<double> > siteWeights;
};

void ComputeContext::findMoleculeGroups() {
    // The first time this is called, we need to identify all the molecules in the system.

    if (moleculeGroups.size() == 0) {
        // Add a ForceInfo that makes sure reordering doesn't break virtual sites.

        addForce(new VirtualSiteInfo(system));

        // First make a list of every other atom to which each atom is connect by a constraint or force group.

        vector<vector<int> > atomBonds(system.getNumParticles());
        for (int i = 0; i < system.getNumConstraints(); i++) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            atomBonds[particle1].push_back(particle2);
            atomBonds[particle2].push_back(particle1);
        }
        for (auto force : forces) {
            vector<int> particles;
            for (int j = 0; j < force->getNumParticleGroups(); j++) {
                force->getParticlesInGroup(j, particles);
                for (int k = 1; k < (int) particles.size(); k++)
                    for (int m = 0; m < k; m++) {
                        atomBonds[particles[k]].push_back(particles[m]);
                        atomBonds[particles[m]].push_back(particles[k]);
                    }
            }
        }

        // Now identify atoms by which molecule they belong to.

        vector<vector<int> > atomIndices = ContextImpl::findMolecules(numAtoms, atomBonds);
        int numMolecules = atomIndices.size();
        vector<int> atomMolecule(numAtoms);
        for (int i = 0; i < (int) atomIndices.size(); i++)
            for (int j = 0; j < (int) atomIndices[i].size(); j++)
                atomMolecule[atomIndices[i][j]] = i;

        // Construct a description of each molecule.

        molecules.resize(numMolecules);
        for (int i = 0; i < numMolecules; i++) {
            molecules[i].atoms = atomIndices[i];
            molecules[i].groups.resize(forces.size());
        }
        for (int i = 0; i < system.getNumConstraints(); i++) {
            int particle1, particle2;
            double distance;
            system.getConstraintParameters(i, particle1, particle2, distance);
            molecules[atomMolecule[particle1]].constraints.push_back(i);
        }
        for (int i = 0; i < (int) forces.size(); i++) {
            vector<int> particles;
            for (int j = 0; j < forces[i]->getNumParticleGroups(); j++) {
                forces[i]->getParticlesInGroup(j, particles);
                if (particles.size() > 0)
                    molecules[atomMolecule[particles[0]]].groups[i].push_back(j);
            }
        }
    }

    // Sort them into groups of identical molecules.

    vector<Molecule> uniqueMolecules;
    vector<vector<int> > moleculeInstances;
    vector<vector<int> > moleculeOffsets;
    for (int molIndex = 0; molIndex < (int) molecules.size(); molIndex++) {
        Molecule& mol = molecules[molIndex];

        // See if it is identical to another molecule.

        bool isNew = true;
        for (int j = 0; j < (int) uniqueMolecules.size() && isNew; j++) {
            Molecule& mol2 = uniqueMolecules[j];
            bool identical = (mol.atoms.size() == mol2.atoms.size() && mol.constraints.size() == mol2.constraints.size());

            // See if the atoms are identical.

            int atomOffset = mol2.atoms[0]-mol.atoms[0];
            for (int i = 0; i < (int) mol.atoms.size() && identical; i++) {
                if (mol.atoms[i] != mol2.atoms[i]-atomOffset || system.getParticleMass(mol.atoms[i]) != system.getParticleMass(mol2.atoms[i]))
                    identical = false;
                for (int k = 0; k < (int) forces.size(); k++)
                    if (!forces[k]->areParticlesIdentical(mol.atoms[i], mol2.atoms[i]))
                        identical = false;
            }

            // See if the constraints are identical.

            for (int i = 0; i < (int) mol.constraints.size() && identical; i++) {
                int c1particle1, c1particle2, c2particle1, c2particle2;
                double distance1, distance2;
                system.getConstraintParameters(mol.constraints[i], c1particle1, c1particle2, distance1);
                system.getConstraintParameters(mol2.constraints[i], c2particle1, c2particle2, distance2);
                if (c1particle1 != c2particle1-atomOffset || c1particle2 != c2particle2-atomOffset || distance1 != distance2)
                    identical = false;
            }

            // See if the force groups are identical.

            for (int i = 0; i < (int) forces.size() && identical; i++) {
                if (mol.groups[i].size() != mol2.groups[i].size())
                    identical = false;
                vector<int> p1, p2;
                for (int k = 0; k < (int) mol.groups[i].size() && identical; k++) {
                    if (!forces[i]->areGroupsIdentical(mol.groups[i][k], mol2.groups[i][k]))
                        identical = false;
                    forces[i]->getParticlesInGroup(mol.groups[i][k], p1);
                    forces[i]->getParticlesInGroup(mol2.groups[i][k], p2);
                    for (int m = 0; m < p1.size(); m++)
                        if (p1[m] != p2[m]-atomOffset)
                            identical = false;
                }
            }
            if (identical) {
                moleculeInstances[j].push_back(molIndex);
                moleculeOffsets[j].push_back(mol.atoms[0]);
                isNew = false;
            }
        }
        if (isNew) {
            uniqueMolecules.push_back(mol);
            moleculeInstances.push_back(vector<int>());
            moleculeInstances[moleculeInstances.size()-1].push_back(molIndex);
            moleculeOffsets.push_back(vector<int>());
            moleculeOffsets[moleculeOffsets.size()-1].push_back(mol.atoms[0]);
        }
    }
    moleculeGroups.resize(moleculeInstances.size());
    for (int i = 0; i < (int) moleculeInstances.size(); i++)
    {
        moleculeGroups[i].instances = moleculeInstances[i];
        moleculeGroups[i].offsets = moleculeOffsets[i];
        vector<int>& atoms = uniqueMolecules[i].atoms;
        moleculeGroups[i].atoms.resize(atoms.size());
        for (int j = 0; j < (int) atoms.size(); j++)
            moleculeGroups[i].atoms[j] = atoms[j]-atoms[0];
    }
}

void ComputeContext::invalidateMolecules() {
    for (int i = 0; i < forces.size(); i++)
        if (invalidateMolecules(forces[i]))
            return;
}

bool ComputeContext::invalidateMolecules(ComputeForceInfo* force, bool checkAtoms, bool checkGroups) {
    if (numAtoms == 0 || !getNonbondedUtilities().getUseCutoff())
        return false;
    bool valid = true;
    int forceIndex = -1;
    for (int i = 0; i < forces.size(); i++)
        if (forces[i] == force)
            forceIndex = i;
    getThreadPool().execute([&] (ThreadPool& threads, int threadIndex) {
        for (int group = 0; valid && group < (int) moleculeGroups.size(); group++) {
            MoleculeGroup& mol = moleculeGroups[group];
            vector<int>& instances = mol.instances;
            vector<int>& offsets = mol.offsets;
            vector<int>& atoms = mol.atoms;
            int numMolecules = instances.size();
            Molecule& m1 = molecules[instances[0]];
            int offset1 = offsets[0];
            int numThreads = threads.getNumThreads();
            int start = max(1, threadIndex*numMolecules/numThreads);
            int end = (threadIndex+1)*numMolecules/numThreads;
            for (int j = start; j < end; j++) {
                // See if the atoms are identical.

                Molecule& m2 = molecules[instances[j]];
                if (checkAtoms) {
                    int offset2 = offsets[j];
                    for (int i = 0; i < (int) atoms.size() && valid; i++) {
                        if (!force->areParticlesIdentical(atoms[i]+offset1, atoms[i]+offset2))
                            valid = false;
                    }
                }

                // See if the force groups are identical.

                if (valid && forceIndex > -1 && checkGroups) {
                    for (int k = 0; k < (int) m1.groups[forceIndex].size() && valid; k++)
                        if (!force->areGroupsIdentical(m1.groups[forceIndex][k], m2.groups[forceIndex][k]))
                            valid = false;
                }
            }
        }
    });
    getThreadPool().waitForThreads();
    if (valid)
        return false;

    // The list of which molecules are identical is no longer valid.  We need to restore the
    // atoms to their original order, rebuild the list of identical molecules, and sort them
    // again.

    resetAtomOrder();
    findMoleculeGroups();
    reorderAtoms();
    return true;
}

void ComputeContext::resetAtomOrder() {
    ContextSelector selector(*this);
    vector<mm_int4> newCellOffsets(numAtoms);
    if (getUseDoublePrecision()) {
        vector<mm_double4> oldPosq(paddedNumAtoms);
        vector<mm_double4> newPosq(paddedNumAtoms, mm_double4(0,0,0,0));
        vector<mm_double4> oldVelm(paddedNumAtoms);
        vector<mm_double4> newVelm(paddedNumAtoms, mm_double4(0,0,0,0));
        getPosq().download(oldPosq);
        getVelm().download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        getPosq().upload(newPosq);
        getVelm().upload(newVelm);
    }
    else if (getUseMixedPrecision()) {
        vector<mm_float4> oldPosq(paddedNumAtoms);
        vector<mm_float4> newPosq(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_float4> oldPosqCorrection(paddedNumAtoms);
        vector<mm_float4> newPosqCorrection(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_double4> oldVelm(paddedNumAtoms);
        vector<mm_double4> newVelm(paddedNumAtoms, mm_double4(0,0,0,0));
        getPosq().download(oldPosq);
        getPosqCorrection().download(oldPosqCorrection);
        getVelm().download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newPosqCorrection[index] = oldPosqCorrection[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        getPosq().upload(newPosq);
        getPosqCorrection().upload(newPosqCorrection);
        getVelm().upload(newVelm);
    }
    else {
        vector<mm_float4> oldPosq(paddedNumAtoms);
        vector<mm_float4> newPosq(paddedNumAtoms, mm_float4(0,0,0,0));
        vector<mm_float4> oldVelm(paddedNumAtoms);
        vector<mm_float4> newVelm(paddedNumAtoms, mm_float4(0,0,0,0));
        getPosq().download(oldPosq);
        getVelm().download(oldVelm);
        for (int i = 0; i < numAtoms; i++) {
            int index = atomIndex[i];
            newPosq[index] = oldPosq[i];
            newVelm[index] = oldVelm[i];
            newCellOffsets[index] = posCellOffsets[i];
        }
        getPosq().upload(newPosq);
        getVelm().upload(newVelm);
    }
    for (int i = 0; i < numAtoms; i++) {
        atomIndex[i] = i;
        posCellOffsets[i] = newCellOffsets[i];
    }
    getAtomIndexArray().upload(atomIndex);
    for (auto listener : reorderListeners)
        listener->execute();
    forceNextReorder = true;
}

void ComputeContext::validateAtomOrder() {
    for (auto& mol : moleculeGroups) {
        for (int atom : mol.atoms) {
            set<int> identical;
            for (int offset : mol.offsets)
                identical.insert(atom+offset);
            for (int i : identical)
                if (identical.find(atomIndex[i]) == identical.end()) {
                    resetAtomOrder();
                    reorderAtoms();
                    return;
                }
        }
    }
}

void ComputeContext::forceReorder() {
    forceNextReorder = true;
}

void ComputeContext::reorderAtoms() {
    atomsWereReordered = false;
    if (numAtoms == 0 || !getNonbondedUtilities().getUseCutoff() || (stepsSinceReorder < 250 && !forceNextReorder)) {
        stepsSinceReorder++;
        return;
    }
    forceNextReorder = false;
    atomsWereReordered = true;
    stepsSinceReorder = 0;
    if (getUseDoublePrecision())
        reorderAtomsImpl<double, mm_double4, double, mm_double4>();
    else if (getUseMixedPrecision())
        reorderAtomsImpl<float, mm_float4, double, mm_double4>();
    else
        reorderAtomsImpl<float, mm_float4, float, mm_float4>();
}

template <class Real, class Real4, class Mixed, class Mixed4>
void ComputeContext::reorderAtomsImpl() {

    // Find the range of positions and the number of bins along each axis.

    vector<Real4> oldPosq(paddedNumAtoms);
    vector<Real4> oldPosqCorrection(paddedNumAtoms);
    vector<Mixed4> oldVelm(paddedNumAtoms);
    getPosq().download(oldPosq);
    getVelm().download(oldVelm);
    if (getUseMixedPrecision())
        getPosqCorrection().download(oldPosqCorrection);
    Real minx = oldPosq[0].x, maxx = oldPosq[0].x;
    Real miny = oldPosq[0].y, maxy = oldPosq[0].y;
    Real minz = oldPosq[0].z, maxz = oldPosq[0].z;
    Vec3 periodicBoxX, periodicBoxY, periodicBoxZ;
    getPeriodicBoxVectors(periodicBoxX, periodicBoxY, periodicBoxZ);
    Vec3 invPeriodicBoxSize(1.0/periodicBoxX[0], 1.0/periodicBoxY[1], 1.0/periodicBoxZ[2]);
    if (getNonbondedUtilities().getUsePeriodic()) {
        minx = miny = minz = 0.0;
        maxx = periodicBoxX[0];
        maxy = periodicBoxY[1];
        maxz = periodicBoxZ[2];
    }
    else {
        for (int i = 1; i < numAtoms; i++) {
            const Real4& pos = oldPosq[i];
            minx = min(minx, pos.x);
            maxx = max(maxx, pos.x);
            miny = min(miny, pos.y);
            maxy = max(maxy, pos.y);
            minz = min(minz, pos.z);
            maxz = max(maxz, pos.z);
        }
    }

    // Loop over each group of identical molecules and reorder them.

    
    vector<int> originalIndex(numAtoms);
    vector<Real4> newPosq(paddedNumAtoms, Real4(0,0,0,0));
    vector<Real4> newPosqCorrection(paddedNumAtoms, Real4(0,0,0,0));
    vector<Mixed4> newVelm(paddedNumAtoms, Mixed4(0,0,0,0));
    vector<mm_int4> newCellOffsets(numAtoms);
    for (auto& mol : moleculeGroups) {
        // Find the center of each molecule.

        int numMolecules = mol.offsets.size();
        vector<int>& atoms = mol.atoms;
        vector<Real4> molPos(numMolecules);
        Real invNumAtoms = (Real) (1.0/atoms.size());
        for (int i = 0; i < numMolecules; i++) {
            molPos[i].x = 0.0f;
            molPos[i].y = 0.0f;
            molPos[i].z = 0.0f;
            for (int j = 0; j < (int)atoms.size(); j++) {
                int atom = atoms[j]+mol.offsets[i];
                const Real4& pos = oldPosq[atom];
                molPos[i].x += pos.x;
                molPos[i].y += pos.y;
                molPos[i].z += pos.z;
            }
            molPos[i].x *= invNumAtoms;
            molPos[i].y *= invNumAtoms;
            molPos[i].z *= invNumAtoms;
            if (molPos[i].x != molPos[i].x)
                throw OpenMMException("Particle coordinate is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan");
        }
        if (getNonbondedUtilities().getUsePeriodic()) {
            // Move each molecule position into the same box.

            for (int i = 0; i < numMolecules; i++) {
                Real4 center = molPos[i];
                int zcell = (int) floor(center.z*invPeriodicBoxSize[2]);
                center.x -= zcell*periodicBoxZ[0];
                center.y -= zcell*periodicBoxZ[1];
                center.z -= zcell*periodicBoxZ[2];
                int ycell = (int) floor(center.y*invPeriodicBoxSize[1]);
                center.x -= ycell*periodicBoxY[0];
                center.y -= ycell*periodicBoxY[1];
                int xcell = (int) floor(center.x*invPeriodicBoxSize[0]);
                center.x -= xcell*periodicBoxX[0];
                if (xcell != 0 || ycell != 0 || zcell != 0) {
                    Real dx = molPos[i].x-center.x;
                    Real dy = molPos[i].y-center.y;
                    Real dz = molPos[i].z-center.z;
                    molPos[i] = center;
                    for (int j = 0; j < (int) atoms.size(); j++) {
                        int atom = atoms[j]+mol.offsets[i];
                        Real4 p = oldPosq[atom];
                        p.x -= dx;
                        p.y -= dy;
                        p.z -= dz;
                        oldPosq[atom] = p;
                        posCellOffsets[atom].x -= xcell;
                        posCellOffsets[atom].y -= ycell;
                        posCellOffsets[atom].z -= zcell;
                    }
                }
            }
        }

        // Select a bin for each molecule, then sort them by bin.

        bool useHilbert = (numMolecules > 5000 || atoms.size() > 8); // For small systems, a simple zigzag curve works better than a Hilbert curve.
        Real binWidth;
        if (useHilbert)
            binWidth = (Real) (max(max(maxx-minx, maxy-miny), maxz-minz)/255.0);
        else
            binWidth = (Real) (0.2*getNonbondedUtilities().getMaxCutoffDistance());
        Real invBinWidth = (Real) (1.0/binWidth);
        int xbins = 1 + (int) ((maxx-minx)*invBinWidth);
        int ybins = 1 + (int) ((maxy-miny)*invBinWidth);
        vector<pair<int, int> > molBins(numMolecules);
        bitmask_t coords[3];
        for (int i = 0; i < numMolecules; i++) {
            int x = (int) ((molPos[i].x-minx)*invBinWidth);
            int y = (int) ((molPos[i].y-miny)*invBinWidth);
            int z = (int) ((molPos[i].z-minz)*invBinWidth);
            int bin;
            if (useHilbert) {
                coords[0] = x;
                coords[1] = y;
                coords[2] = z;
                bin = (int) hilbert_c2i(3, 8, coords);
            }
            else {
                int yodd = y&1;
                int zodd = z&1;
                bin = z*xbins*ybins;
                bin += (zodd ? ybins-y : y)*xbins;
                bin += (yodd ? xbins-x : x);
            }
            molBins[i] = pair<int, int>(bin, i);
        }
        sort(molBins.begin(), molBins.end());

        // Reorder the atoms.

        for (int i = 0; i < numMolecules; i++) {
            for (int atom : atoms) {
                int oldIndex = mol.offsets[molBins[i].second]+atom;
                int newIndex = mol.offsets[i]+atom;
                originalIndex[newIndex] = atomIndex[oldIndex];
                newPosq[newIndex] = oldPosq[oldIndex];
                if (getUseMixedPrecision())
                    newPosqCorrection[newIndex] = oldPosqCorrection[oldIndex];
                newVelm[newIndex] = oldVelm[oldIndex];
                newCellOffsets[newIndex] = posCellOffsets[oldIndex];
            }
        }
    }

    // Update the arrays.

    ContextSelector selector(*this);
    for (int i = 0; i < numAtoms; i++) {
        atomIndex[i] = originalIndex[i];
        posCellOffsets[i] = newCellOffsets[i];
    }
    getPosq().upload(newPosq);
    if (getUseMixedPrecision())
        getPosqCorrection().upload(newPosqCorrection);
    getVelm().upload(newVelm);
    getAtomIndexArray().upload(atomIndex);
    for (auto listener : reorderListeners)
        listener->execute();
}

void ComputeContext::addReorderListener(ReorderListener* listener) {
    reorderListeners.push_back(listener);
}

void ComputeContext::addPreComputation(ForcePreComputation* computation) {
    preComputations.push_back(computation);
}

void ComputeContext::addPostComputation(ForcePostComputation* computation) {
    postComputations.push_back(computation);
}

int ComputeContext::findLegalFFTDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

struct ComputeContext::WorkThread::ThreadData {
    ThreadData(std::queue<ComputeContext::WorkTask*>& tasks, bool& waiting,  bool& finished, bool& threwException, OpenMMException& stashedException,
            pthread_mutex_t& queueLock, pthread_cond_t& waitForTaskCondition, pthread_cond_t& queueEmptyCondition) :
        tasks(tasks), waiting(waiting), finished(finished), threwException(threwException), stashedException(stashedException),
        queueLock(queueLock), waitForTaskCondition(waitForTaskCondition), queueEmptyCondition(queueEmptyCondition) {
    }
    std::queue<ComputeContext::WorkTask*>& tasks;
    bool& waiting;
    bool& finished;
    bool& threwException;
    OpenMMException& stashedException;
    pthread_mutex_t& queueLock;
    pthread_cond_t& waitForTaskCondition;
    pthread_cond_t& queueEmptyCondition;
};

static void* threadBody(void* args) {
    ComputeContext::WorkThread::ThreadData& data = *reinterpret_cast<ComputeContext::WorkThread::ThreadData*>(args);
    while (!data.finished || data.tasks.size() > 0) {
        pthread_mutex_lock(&data.queueLock);
        while (data.tasks.empty() && !data.finished) {
            data.waiting = true;
            pthread_cond_signal(&data.queueEmptyCondition);
            pthread_cond_wait(&data.waitForTaskCondition, &data.queueLock);
        }
        // If we keep going after having caught an exception once, next tasks will likely throw too and we don't want the initial exception overshadowed.
        while (data.threwException && !data.tasks.empty()) {
            delete data.tasks.front();
            data.tasks.pop();
        }
        ComputeContext::WorkTask* task = NULL;
        if (!data.tasks.empty()) {
            data.waiting = false;
            task = data.tasks.front();
            data.tasks.pop();
        }
        pthread_mutex_unlock(&data.queueLock);
        if (task != NULL) {
            try {
                task->execute();
            }
            catch (const OpenMMException& e) {
                data.threwException = true;
                data.stashedException = e;
            }
            delete task;
        }
    }
    data.waiting = true;
    pthread_cond_signal(&data.queueEmptyCondition);
    delete &data;
    return 0;
}

ComputeContext::WorkThread::WorkThread() : waiting(true), finished(false), threwException(false), stashedException("Default WorkThread exception. This should never be thrown.") {
    pthread_mutex_init(&queueLock, NULL);
    pthread_cond_init(&waitForTaskCondition, NULL);
    pthread_cond_init(&queueEmptyCondition, NULL);
    ThreadData* data = new ThreadData(tasks, waiting, finished, threwException, stashedException, queueLock, waitForTaskCondition, queueEmptyCondition);
    pthread_create(&thread, NULL, threadBody, data);
}

ComputeContext::WorkThread::~WorkThread() {
    pthread_mutex_lock(&queueLock);
    finished = true;
    pthread_cond_broadcast(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
    pthread_join(thread, NULL);
    pthread_mutex_destroy(&queueLock);
    pthread_cond_destroy(&waitForTaskCondition);
    pthread_cond_destroy(&queueEmptyCondition);
}

void ComputeContext::WorkThread::addTask(ComputeContext::WorkTask* task) {
    pthread_mutex_lock(&queueLock);
    tasks.push(task);
    waiting = false;
    pthread_cond_signal(&waitForTaskCondition);
    pthread_mutex_unlock(&queueLock);
}

bool ComputeContext::WorkThread::isWaiting() {
    return waiting;
}

bool ComputeContext::WorkThread::isFinished() {
    return finished;
}

bool ComputeContext::WorkThread::isCurrentThread() {
    return (pthread_self() == thread);
}

void ComputeContext::WorkThread::flush() {
    pthread_mutex_lock(&queueLock);
    while (!waiting)
       pthread_cond_wait(&queueEmptyCondition, &queueLock);
    pthread_mutex_unlock(&queueLock);
    if (threwException) {
        threwException = false;
        throw stashedException;
    }
}
