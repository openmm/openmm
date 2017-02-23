/* 
 * Reference implementation of PME reciprocal space interactions.
 * 
 * Copyright (c) 2009, Erik Lindahl, Rossen Apostolov, Szilard Pall
 * All rights reserved.
 * Contact: lindahl@cbr.su.se Stockholm University, Sweden.
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this 
 * list of conditions and the following disclaimer. Redistributions in binary 
 * form must reproduce the above copyright notice, this list of conditions and 
 * the following disclaimer in the documentation and/or other materials provided 
 * with the distribution.
 * Neither the name of the author/university nor the names of its contributors may 
 * be used to endorse or promote products derived from this software without 
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __ReferencePME_H__
#define __ReferencePME_H__

#include "openmm/Vec3.h"
#include "openmm/internal/windowsExport.h"
#include <vector>

namespace OpenMM {

typedef double rvec[3];


typedef struct pme *
pme_t;

/*
 * Initialize a PME calculation and set up data structures
 *
 * Arguments:
 * 
 * ppme        Pointer to an opaque pme_t object
 * ewaldcoeff  Coefficient derived from the beta factor to participate
 *             direct/reciprocal space. See gromacs code for documentation!
 *             We assume that you are using nm units...
 * natoms      Number of atoms to set up data structure sof
 * ngrid       Size of the full pme grid
 * pme_order   Interpolation order, almost always 4
 * epsilon_r   Dielectric coefficient, typically 1.0.
 */
int OPENMM_EXPORT
pme_init(pme_t* ppme,
         double ewaldcoeff,
         int natoms,
         const int ngrid[3],
         int pme_order,
         double epsilon_r);

/*
 * Evaluate reciprocal space PME energy and forces.
 *
 * Args:
 *
 * pme         Opaque pme_t object, must have been initialized with pme_init()
 * x           Pointer to coordinate data array (nm)
 * f           Pointer to force data array (will be written as kJ/mol/nm)
 * charge      Array of charges (units of e)
 * box         Simulation cell dimensions (nm)
 * energy      Total energy (will be written in units of kJ/mol)
 */
int OPENMM_EXPORT
pme_exec(pme_t pme,
         const std::vector<OpenMM::Vec3>& atomCoordinates,
         std::vector<OpenMM::Vec3>& forces,
         const std::vector<double>& charges,
         const OpenMM::Vec3 periodicBoxVectors[3],
         double* energy);


/**
 * Evaluate reciprocal space PME dispersion energy and forces.
 *
 * Args:
 *
 * pme         Opaque pme_t object, must have been initialized with pme_init()
 * x           Pointer to coordinate data array (nm)
 * f           Pointer to force data array (will be written as kJ/mol/nm)
 * c6s         Array of c6 coefficients (units of sqrt(kJ/mol).nm^3 )
 * box         Simulation cell dimensions (nm)
 * energy      Total energy (will be written in units of kJ/mol)
 */
int OPENMM_EXPORT
pme_exec_dpme(pme_t pme,
              const std::vector<OpenMM::Vec3>& atomCoordinates,
              std::vector<OpenMM::Vec3>& forces,
              const std::vector<double>& c6s,
              const OpenMM::Vec3 periodicBoxVectors[3],
              double* energy);




/* Release all memory in pme structure */
int OPENMM_EXPORT
pme_destroy(pme_t    pme);

} // namespace OpenMM

#endif // __ReferencePME_H__
