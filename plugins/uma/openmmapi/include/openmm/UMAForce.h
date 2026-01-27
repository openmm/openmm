#ifndef OPENMM_UMAFORCE_H_
#define OPENMM_UMAFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                             *
 * Portions copyright (c) 2025 Stanford University and the Authors.            *
 * Authors: Muhammad Hasyim                                                    *
 * Contributors:                                                               *
 *                                                                             *
 * Permission is hereby granted, free of charge, to any person obtaining a     *
 * copy of this software and associated documentation files (the "Software"),  *
 * to deal in the Software without restriction, including without limitation   *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,    *
 * and/or sell copies of the Software, and to permit persons to whom the       *
 * Software is furnished to do so, subject to the following conditions:        *
 *                                                                             *
 * The above copyright notice and this permission notice shall be included in  *
 * all copies or substantial portions of the Software.                         *
 *                                                                             *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL     *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR       *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE   *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                      *
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "internal/windowsExportUMA.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class implements the Universal Molecular Atomistic (UMA) potential
 * from FAIRChem using PyTorch models integrated via libtorch C++ API.
 * 
 * UMA is a family of pre-trained machine learning potentials that can predict
 * energies and forces for molecular systems, materials, and catalytic systems.
 * This implementation provides native C++ integration for maximum performance,
 * including zero-copy GPU memory sharing between OpenMM and PyTorch.
 * 
 * To use this force, specify a model name and task type:
 * 
 * \code
 * UMAForce* force = new UMAForce("uma-s-1p1", UMAForce::OMol);
 * force->setCharge(0);
 * force->setSpin(1);
 * system.addForce(force);
 * \endcode
 * 
 * The force supports multiple task types for different applications:
 * - OMol: Molecular systems (requires charge and spin)
 * - OMat: Materials and crystals
 * - OC20: Open Catalyst 2020 dataset
 * - ODAC: Direct Air Capture reactions
 * - OMC: Multi-component systems
 * 
 * Models are automatically loaded from the HuggingFace cache or can be
 * specified with a custom path using setModelPath().
 */
class OPENMM_EXPORT_UMA UMAForce : public Force {
public:
    /**
     * Task types supported by UMA models.
     */
    enum TaskType {
        OMol = 0,   ///< Molecular systems
        OMat = 1,   ///< Materials
        OC20 = 2,   ///< Open Catalyst 2020
        ODAC = 3,   ///< Catalyst reactions
        OMC = 4     ///< Multi-component
    };
    
    /**
     * Create a UMAForce.
     *
     * @param modelName   the name of the UMA model (e.g., "uma-s-1p1", "uma-m-1p1")
     * @param taskType    the task type for this force
     */
    UMAForce(const std::string& modelName, TaskType taskType);
    
    /**
     * Get the model name.
     */
    const std::string& getModelName() const {
        return modelName;
    }
    
    /**
     * Get the task type.
     */
    TaskType getTaskType() const {
        return taskType;
    }
    
    /**
     * Set a custom model path instead of using the default HuggingFace cache.
     *
     * @param path   the path to the model .pt file
     */
    void setModelPath(const std::string& path);
    
    /**
     * Get the custom model path, or empty string if using default.
     */
    const std::string& getModelPath() const {
        return modelPath;
    }
    
    /**
     * Set inference settings for the model ("default" or "turbo").
     *
     * @param settings   the inference settings string
     */
    void setInferenceSettings(const std::string& settings);
    
    /**
     * Get the inference settings.
     */
    const std::string& getInferenceSettings() const {
        return inferenceSettings;
    }
    
    /**
     * Set the CUDA device index to use for inference.
     *
     * @param deviceIndex   the CUDA device index (default: 0)
     */
    void setDevice(int deviceIndex);
    
    /**
     * Get the device index.
     */
    int getDevice() const {
        return deviceIndex;
    }
    
    /**
     * Set the total charge for OMol tasks.
     *
     * @param charge   the total system charge
     */
    void setCharge(int charge);
    
    /**
     * Get the charge.
     */
    int getCharge() const {
        return charge;
    }
    
    /**
     * Set the spin multiplicity for OMol tasks.
     *
     * @param spin   the spin multiplicity
     */
    void setSpin(int spin);
    
    /**
     * Get the spin multiplicity.
     */
    int getSpin() const {
        return spin;
    }
    
    /**
     * Set a subset of atoms to use for the force calculation.
     * If not set, all atoms in the system are used.
     *
     * @param atoms   vector of atom indices
     */
    void setAtomSubset(const std::vector<int>& atoms);
    
    /**
     * Get the atom subset, or empty vector if all atoms are used.
     */
    const std::vector<int>& getAtomSubset() const {
        return atomSubset;
    }
    
    /**
     * Set atomic reference energies for energy corrections.
     *
     * @param refs   map from atomic number to reference energy in eV
     */
    void setAtomicReferences(const std::map<int, double>& refs);
    
    /**
     * Get the atomic reference energies.
     */
    const std::map<int, double>& getAtomicReferences() const {
        return atomicRefs;
    }
    
    /**
     * Set the cutoff radius for neighbor list construction.
     *
     * @param cutoff   the cutoff radius in nm
     */
    void setCutoffRadius(double cutoff);
    
    /**
     * Get the cutoff radius in nm.
     */
    double getCutoffRadius() const {
        return cutoffRadius;
    }
    
    /**
     * Set the maximum number of neighbors per atom.
     *
     * @param maxNeighbors   the maximum number of neighbors
     */
    void setMaxNeighbors(int maxNeighbors);
    
    /**
     * Get the maximum number of neighbors.
     */
    int getMaxNeighbors() const {
        return maxNeighbors;
    }
    
    /**
     * Returns whether this force uses periodic boundary conditions.
     */
    bool usesPeriodicBoundaryConditions() const override;
    
    /**
     * Set whether this force uses periodic boundary conditions.
     *
     * @param periodic   whether to use periodic boundaries
     */
    void setUsesPeriodicBoundaryConditions(bool periodic);
    
protected:
    ForceImpl* createImpl() const override;
    
private:
    std::string modelName;
    std::string modelPath;
    std::string inferenceSettings;
    TaskType taskType;
    int charge;
    int spin;
    int deviceIndex;
    bool usePeriodic;
    double cutoffRadius;
    int maxNeighbors;
    std::vector<int> atomSubset;
    std::map<int, double> atomicRefs;
};

} // namespace OpenMM

#endif /*OPENMM_UMAFORCE_H_*/
