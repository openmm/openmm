/* 
 * Reference implementation of PME reciprocal space interactions.
 * 
 * Copyright (c) 2009-2025 Erik Lindahl, Rossen Apostolov, Szilard Pall, Peter Eastman
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
#include <array>
#include <complex>
#include <vector>

namespace OpenMM {

class OPENMM_EXPORT ReferencePME {
public:
    /*
     * Initialize a PME calculation and set up data structures
     *
     * Arguments:
     * 
     * ewaldcoeff  Coefficient derived from the beta factor to participate
     *             direct/reciprocal space. See gromacs code for documentation!
     *             We assume that you are using nm units...
     * natoms      Number of atoms to set up data structure sof
     * ngrid       Size of the full pme grid
     * pme_order   Interpolation order, almost always 4
     * epsilon_r   Dielectric coefficient, typically 1.0.
     */
    ReferencePME(double ewaldcoeff, int natoms, const int ngrid[3], int pme_order, double epsilon_r);

    /*
     * Evaluate reciprocal space PME energy and forces.
     *
     * Args:
     *
     * atomCoordinates     Pointer to coordinate data array (nm)
     * forces              Pointer to force data array (will be written as kJ/mol/nm)
     * charges             Array of charges (units of e)
     * periodicBoxVectors  Simulation cell dimensions (nm)
     * energy              Total energy (will be written in units of kJ/mol)
     */
    void exec(const std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces, const std::vector<double>& charges,
            const OpenMM::Vec3 periodicBoxVectors[3], double& energy);

    /*
     * Evaluate reciprocal space PME energy and charge derivatives.
     *
     * Args:
     *
     * atomCoordinates     Pointer to coordinate data array (nm)
     * chargeDerivatives   Pointer to charge derivative data array (will be written as kJ/mol/e)
     * chargeIndices       Pointer to array of indices of particles to compute charge derivatives for
     * charges             Array of charges (units of e)
     * periodicBoxVectors  Simulation cell dimensions (nm)
     */
    void exec_charge_derivatives(const std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<double>& chargeDerivatives,
            const std::vector<int>& chargeIndices, const std::vector<double>& charges, const OpenMM::Vec3 periodicBoxVectors[3]);

    /**
     * Evaluate reciprocal space PME dispersion energy and forces.
     *
     * Args:
     *
     * atomCoordinates     Pointer to coordinate data array (nm)
     * forces              Pointer to force data array (will be written as kJ/mol/nm)
     * c6s                 Array of c6 coefficients (units of sqrt(kJ/mol).nm^3 )
     * periodicBoxVectors  Simulation cell dimensions (nm)
     * energy              Total energy (will be written in units of kJ/mol)
     */
    void exec_dpme(const std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces, const std::vector<double>& c6s,
            const OpenMM::Vec3 periodicBoxVectors[3], double& energy);
private:
    void calculate_bsplines_moduli();
    void update_grid_index_and_fraction(const std::vector<Vec3>& atomCoordinates, const Vec3 recipBoxVectors[3]);
    void update_bsplines();
    void grid_spread_charge(const std::vector<double>& charges);
    void pme_reciprocal_convolution(const Vec3 periodicBoxVectors[3], const Vec3 recipBoxVectors[3], double& energy);
    void dpme_reciprocal_convolution(const Vec3 periodicBoxVectors[3], const Vec3 recipBoxVectors[3], double& energy);
    void grid_interpolate_force(const Vec3 recipBoxVectors[3], const std::vector<double>& charges, std::vector<Vec3>& forces);
    void grid_interpolate_charge_derivatives(const Vec3 recipBoxVectors[3], const std::vector<double>& charges,
            std::vector<double>& chargeDerivatives, const std::vector<int>& chargeIndices);
    int natoms;
    double ewaldcoeff;
    std::vector<std::complex<double> > grid; /* Memory for the grid we spread charges on.
                                              * Element (i,j,k) is accessed as:
                                              * grid[i*ngrid[1]*ngrid[2] + j*ngrid[2] + k] */
    int ngrid[3];             /* Total grid dimensions (all data is complex!) */
    int order;                /* PME interpolation order. Almost always 4 */

    /* Data for bspline interpolation, see the Essman PME paper */
    std::vector<double> bsplines_moduli[3];   /* 3 pointers, to x/y/z bspline moduli, each of length ngrid[x/y/z]   */
    std::vector<double> bsplines_theta[3];    /* each of x/y/z has length order*natoms */
    std::vector<double> bsplines_dtheta[3];   /* each of x/y/z has length order*natoms */

    std::vector<std::array<int, 3> > particleindex; /* Array of length natoms. Each element is
                                                     * three ints that specify the grid
                                                     * indices for that particular atom. Updated every step! */
    std::vector<Vec3> particlefraction; /* Array of length natoms. Fractional offset in the grid for
                                         * each atom in all three dimensions. */

    /* Further explanation of index/fraction:
     *
     * Assume we have a cell of size 10*10*10nm, and a total grid dimension of 100*100*100 cells.
     * In other words, each cell is 0.1*0.1*0.1 nm.
     *
     * If particle i has coordinates { 0.543 , 6.235 , -0.73 }, we will get:
     *
     * particleindex[i]    = { 5 , 62 , 92 }         (-0.73 + 10 = 9.27, we always apply PBC for grid calculations!)
     * particlefraction[i] = { 0.43 , 0.35 , 0.7 }   (this is the fraction of the cell length where the atom is)
     *
     * (The reason for precaculating / storing these is that it gets a bit more complex for triclinic cells :-)
     *
     * In the current code version we might assume that a particle is not more than a whole box length away from
     * the central cell, i.e., in this case we would assume all coordinates fall in -10 nm < x,y,z < 20 nm.
     */

    double epsilon_r; /* Dielectric coefficient to use, typically 1.0 */
};

} // namespace OpenMM

#endif // __ReferencePME_H__
