/*
 * Reference implementation of PME reciprocal space interactions.
 *
 * Copyright (c) 2009-2025, Erik Lindahl, Rossen Apostolov, Szilard Pall, Peter Eastman, Evan Pretti
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

#include <cmath>
#include <vector>

#include "ReferencePME.h"
#include "ReferenceForce.h"
#include "SimTKOpenMMRealType.h"

#ifdef _MSC_VER
  #define POCKETFFT_NO_VECTORS
#endif
#include "pocketfft_hdronly.h"

using namespace std;

namespace OpenMM {

/* Only called once from constructor, performance does not matter! */
void ReferencePME::calculate_bsplines_moduli() {
    int nmax = 0;
    for (int d = 0; d < 3; d++) {
        nmax = (ngrid[d] > nmax) ? ngrid[d] : nmax;
        bsplines_moduli[d].resize(ngrid[d]);
    }

    /* temp storage in this routine */
    vector<double> data(order);
    vector<double> ddata(order);
    vector<double> bsplines_data(nmax);

    data[order-1]=0;
    data[1]=0;
    data[0]=1;

    for (int k = 3; k < order; k++) {
        double div=1.0/(k-1.0);
        data[k-1]=0;
        for (int l = 1; l < (k-1); l++)
            data[k-l-1]=div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
        data[0]=div*data[0];
    }

    /* differentiate */
    ddata[0]=-data[0];
    for (int k = 1 ; k < order; k++)
        ddata[k]=data[k-1]-data[k];

    double div=1.0/(order-1);
    data[order-1]=0;

    for (int l = 1; l < (order-1); l++)
        data[order-l-1]=div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
    data[0]=div*data[0];

    for (int i = 0; i < nmax; i++)
        bsplines_data[i]=0;
    for (int i = 1; i <= order; i++)
        bsplines_data[i]=data[i-1];

    /* Evaluate the actual bspline moduli for X/Y/Z */
    for (int d = 0; d < 3; d++) {
        int ndata = ngrid[d];
        for (int i = 0; i < ndata; i++) {
            double sc = 0;
            double ss = 0;
            for (int j = 0; j < ndata; j++) {
                double arg=(2.0*M_PI*i*j)/ndata;
                sc+=bsplines_data[j]*cos(arg);
                ss+=bsplines_data[j]*sin(arg);
            }
            bsplines_moduli[d][i]=sc*sc+ss*ss;
        }
        for (int i = 0; i < ndata; i++) {
            if (bsplines_moduli[d][i]<1.0e-7)
                bsplines_moduli[d][i]=(bsplines_moduli[d][(i-1+ndata)%ndata]+bsplines_moduli[d][(i+1)%ndata])/2;
        }
    }
}

void ReferencePME::update_grid_index_and_fraction(const vector<Vec3>& atomCoordinates, const Vec3 recipBoxVectors[3]) {
    for (int i = 0; i < natoms; i++) {
        /* Index calculation (Look mom, no conditionals!):
         *
         * Both for Cuda and modern CPUs it is nice to avoid conditionals, but we still need to apply periodic boundary conditions.
         * Instead of having loops to add/subtract the box dimension, we do it this way:
         *
         * 1. First add the box size, to make sure this atom coordinate isnt -0.1 or something.
         *    After this we assume all fractional box positions are *positive*.
         *    The reason for this is that we always want to round coordinates _down_ to get
         *    their grid index, and when taking the integer part of -3.4 we would get -3, not -4 as we want.
         *    Since we anyway need the grid indices to fall in the central box, it is more convenient
         *    to first manipulate the coordinates to be positive.
         * 2. Convert to integer grid index
         *    Since we have added a whole box unit in step 1, this index might actually be larger than
         *    the grid dimension. Examples, assuming 10*10*10nm box and grid dimension 100*100*100 (spacing 0.1 nm):
         *
         *    coordinate is { 0.543 , 6.235 , -0.73 }
         *
         *    x[i][d]/box[d]                      becomes   { 0.0543 , 0.6235 , -0.073 }
         *    (x[i][d]/box[d] + 1.0)              becomes   { 1.0543 , 1.6235 , 0.927 }
         *    (x[i][d]/box[d] + 1.0)*ngrid[d]     becomes   { 105.43 , 162.35 , 92.7 }
         *
         *    integer part is now { 105 , 162 , 92 }
         *
         *    The fraction is calculates as t-ti, which becomes { 0.43 , 0.35 , 0.7 }
         *
         * 3. Take the first integer index part (which can be larger than the grid) modulo the grid dimension
         *
         *    Now we get { 5 , 62 , 92 }
         *
         *    Voila, both index and fraction, entirely without conditionals. The one limitation here is that
         *    we only add one box length, so if the particle had a coordinate <=-10.0, we would be screwed.
         *    In principle we can of course add 100.0, but that just moves the problem, it doesnt solve it.
         *    In practice, MD programs will apply PBC to reset particles inside the central box to avoid
         *    numerical problems, so this shouldnt cause any problems.
         *    (And, by adding 100.0 box lengths, we would lose a bit of numerical accuracy here!)
         */
        Vec3 coord = atomCoordinates[i];
        for (int d = 0; d < 3; d++) {
            double t = coord[0]*recipBoxVectors[0][d]+coord[1]*recipBoxVectors[1][d]+coord[2]*recipBoxVectors[2][d];
            t = (t-floor(t))*ngrid[d];
            int ti = (int) t;

            particlefraction[i][d] = t - ti;
            particleindex[i][d]    = ti % ngrid[d];
        }
    }
}


/* Ugly bspline calculation taken from Tom Dardens reference equations.
 * This probably very sub-optimal in Cuda? Separate kernel?
 *
 * In practice, it might help to require order=4 for the cuda port.
 */
void ReferencePME::update_bsplines() {
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j < 3; j++) {
            /* dr is relative offset from lower cell limit */
            double dr = particlefraction[i][j];

            double* data  = &(bsplines_theta[j][i*order]);
            double* ddata = &(bsplines_dtheta[j][i*order]);
            data[order-1] = 0;
            data[1]       = dr;
            data[0]       = 1-dr;

            for (int k = 3; k < order; k++) {
                double div = 1.0/(k-1.0);
                data[k-1] = div*dr*data[k-2];
                for (int l = 1; l < (k-1); l++)
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)*data[k-l-1]);
                data[0] = div*(1-dr)*data[0];
            }

            /* differentiate */
            ddata[0] = -data[0];

            for (int k = 1; k < order; k++)
                ddata[k] = data[k-1]-data[k];

            double div = 1.0/(order-1);
            data[order-1] = div*dr*data[order-2];

            for (int l = 1; l < (order-1); l++)
                data[order-l-1] = div*((dr+l)*data[order-l-2]+(order-l-dr)*data[order-l-1]);
            data[0] = div*(1-dr)*data[0];
        }
    }
}


void ReferencePME::grid_spread_charge(const vector<double>& charges) {
    /* Reset the grid */
    for (int i = 0; i < ngrid[0]*ngrid[1]*ngrid[2]; i++)
        grid[i] = complex<double>(0, 0);

    for (int i = 0; i < natoms; i++) {
        double q = charges[i];

        /* Grid index for the actual atom position */
        int x0index = particleindex[i][0];
        int y0index = particleindex[i][1];
        int z0index = particleindex[i][2];

        /* Bspline factors for this atom in each dimension , calculated from fractional coordinates */
        double* thetax  = &(bsplines_theta[0][i*order]);
        double* thetay  = &(bsplines_theta[1][i*order]);
        double* thetaz  = &(bsplines_theta[2][i*order]);

        /* Loop over norder*norder*norder (typically 5*5*5) neighbor cells.
         *
         * As a neat optimization, we only spread in the forward direction, but apply PBC!
         *
         * Since we are going to do an FFT on the grid, it doesn't matter where the data is,
         * in frequency space the result will be the same.
         *
         * So, the influence function (bsplines) will probably be something like (0.15,0.35,0.35,0.15),
         * with largest weight 2-3 steps forward (you don't need to understand that for the implementation :-)
         * Effectively, you can look at this as translating the entire grid.
         *
         * Why do we do this stupid thing?
         *
         * 1) The loops get much simpler
         * 2) Just looking forward will hopefully get us more cache hits
         * 3) When we parallelize things, we only need to communicate in one direction instead of two!
         */

        for (int ix = 0; ix < order; ix++) {
            /* Calculate index, apply PBC so we spread to index 0/1/2 when a particle is close to the upper limit of the grid */
            int xindex = (x0index + ix) % ngrid[0];

            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ngrid[1];

                for (int iz = 0; iz < order; iz++) {
                    /* Can be optimized, but we keep it simple here */
                    int zindex = (z0index + iz) % ngrid[2];
                    /* Calculate index in the charge grid */
                    int index = xindex*ngrid[1]*ngrid[2] + yindex*ngrid[2] + zindex;
                    /* Add the charge times the bspline spread/interpolation factors to this grid position */
                    grid[index] += q*thetax[ix]*thetay[iy]*thetaz[iz];
                }
            }
        }
    }
}



void ReferencePME::pme_reciprocal_convolution(const Vec3 periodicBoxVectors[3], const Vec3 recipBoxVectors[3], double& energy) {
    int nx = ngrid[0];
    int ny = ngrid[1];
    int nz = ngrid[2];

    double one_4pi_eps = ONE_4PI_EPS0/epsilon_r;
    double factor = M_PI*M_PI/(ewaldcoeff*ewaldcoeff);
    double boxfactor = M_PI*periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2];

    double esum = 0;

    double maxkx = (nx+1)/2;
    double maxky = (ny+1)/2;
    double maxkz = (nz+1)/2;

    for (int kx = 0; kx < nx; kx++) {
        /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
        double mx  = (kx<maxkx) ? kx : (kx-nx);
        double mhx = mx*recipBoxVectors[0][0];
        double bx  = boxfactor*bsplines_moduli[0][kx];

        for (int ky = 0; ky < ny; ky++) {
            /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
            double my  = (ky<maxky) ? ky : (ky-ny);
            double mhy = mx*recipBoxVectors[1][0]+my*recipBoxVectors[1][1];
            double by  = bsplines_moduli[1][ky];

            for (int kz = 0; kz < nz; kz++) {
                /* Pointer to the grid cell in question */
                complex<double>& ptr = grid[kx*ny*nz + ky*nz + kz];

                /* The zero frequency term is undefined due to division by the frequency below.  Set this term to zero;
                 * in the case that the net charge of the system is non-zero, this is equivalent to applying a uniform
                 * neutralizing background charge density.  The contribution to the energy and charge derivatives of
                 * this neutralizing plasma is applied elsewhere.  If this term is not zeroed, however, energies and
                 * forces will be unaffected but charge derivatives for non-neutral systems will be incorrect!
                 */
                if (kx==0 && ky==0 && kz==0) {
                    ptr = 0;
                    continue;
                }

                /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
                double mz        = (kz<maxkz) ? kz : (kz-nz);
                double mhz       = mx*recipBoxVectors[2][0]+my*recipBoxVectors[2][1]+mz*recipBoxVectors[2][2];

                /* Get grid data for this frequency */
                double d1        = ptr.real();
                double d2        = ptr.imag();

                /* Calculate the convolution - see the Essman/Darden paper for the equation! */
                double m2        = mhx*mhx+mhy*mhy+mhz*mhz;
                double bz        = bsplines_moduli[2][kz];
                double denom     = m2*bx*by*bz;

                double eterm     = one_4pi_eps*exp(-factor*m2)/denom;

                /* write back convolution data to grid */
                ptr.real(d1*eterm);
                ptr.imag(d2*eterm);

                double struct2   = (d1*d1+d2*d2);

                /* Long-range PME contribution to the energy for this frequency */
                double ets2      = eterm*struct2;
                esum     += ets2;
            }
        }
    }

    /* The factor 0.5 is nothing special, but it is better to have it here than inside the loop :-) */
    energy = 0.5*esum;
}


void ReferencePME::dpme_reciprocal_convolution(const Vec3 periodicBoxVectors[3], const Vec3 recipBoxVectors[3], double& energy) {
    int nx = ngrid[0];
    int ny = ngrid[1];
    int nz = ngrid[2];

    double boxfactor = -2*M_PI*sqrt(M_PI) / (6.0*periodicBoxVectors[0][0]*periodicBoxVectors[1][1]*periodicBoxVectors[2][2]);

    double esum = 0;

    double maxkx = (nx+1)/2;
    double maxky = (ny+1)/2;
    double maxkz = (nz+1)/2;

    double bfac = M_PI / ewaldcoeff;
    double fac1 = 2.0*M_PI*M_PI*M_PI*sqrt(M_PI);
    double fac2 = ewaldcoeff*ewaldcoeff*ewaldcoeff;
    double fac3 = -2.0*ewaldcoeff*M_PI*M_PI;

    for (int kx = 0; kx < nx; kx++) {
        /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
        double mx  = ((kx<maxkx) ? kx : (kx-nx));
        double mhx = mx*recipBoxVectors[0][0];
        double bx  = bsplines_moduli[0][kx];

        for (int ky = 0; ky < ny; ky++) {
            /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
            double my  = ((ky<maxky) ? ky : (ky-ny));
            double mhy = mx*recipBoxVectors[1][0]+my*recipBoxVectors[1][1];
            double by  = bsplines_moduli[1][ky];

            for (int kz = 0; kz < nz; kz++) {
                /*
                 * Unlike the Coulombic case, there's an m=0 term so all terms are considered here.
                 */

                /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
                double mz  = ((kz<maxkz) ? kz : (kz-nz));
                double mhz = mx*recipBoxVectors[2][0]+my*recipBoxVectors[2][1]+mz*recipBoxVectors[2][2];

                /* Pointer to the grid cell in question */
                complex<double>& ptr = grid[kx*ny*nz + ky*nz + kz];

                /* Get grid data for this frequency */
                double d1 = ptr.real();
                double d2 = ptr.imag();

                /* Calculate the convolution - see the Essman/Darden paper for the equation! */
                double m2    = mhx*mhx+mhy*mhy+mhz*mhz;
                double bz    = bsplines_moduli[2][kz];
                double denom = boxfactor / (bx*by*bz);

                double m = sqrt(m2);
                double m3 = m*m2;
                double b = bfac*m;
                double expfac = -b*b;
                double erfcterm = erfc(b);
                double expterm = exp(expfac);

                double eterm = (fac1*erfcterm*m3 + expterm*(fac2 + fac3*m2)) * denom;

                /* write back convolution data to grid */
                ptr.real(d1*eterm);
                ptr.imag(d2*eterm);

                double struct2 = (d1*d1+d2*d2);

                /* Long-range PME contribution to the energy for this frequency */
                double ets2 = eterm*struct2;
                esum += ets2;
            }
        }
    }
    // Remember the C6 energy is attractive, hence the negative sign.
    energy = 0.5*esum;
}


void ReferencePME::grid_interpolate_force(const Vec3 recipBoxVectors[3], const vector<double>& charges, vector<Vec3>& forces) {
    /* This is almost identical to the charge spreading routine! */

    for (int i = 0; i < natoms; i++) {
        double fx = 0, fy = 0, fz = 0;

        double q = charges[i];

        /* Grid index for the actual atom position */
        int x0index = particleindex[i][0];
        int y0index = particleindex[i][1];
        int z0index = particleindex[i][2];

        /* Bspline factors for this atom in each dimension , calculated from fractional coordinates */
        double* thetax  = &(bsplines_theta[0][i*order]);
        double* thetay  = &(bsplines_theta[1][i*order]);
        double* thetaz  = &(bsplines_theta[2][i*order]);
        double* dthetax = &(bsplines_dtheta[0][i*order]);
        double* dthetay = &(bsplines_dtheta[1][i*order]);
        double* dthetaz = &(bsplines_dtheta[2][i*order]);

        /* See grid_spread_charge() for comments about the order here, and only interpolation in one direction */

        /* Since we will add order^3 (typically 5*5*5=125) terms to the force on each particle, we use temporary fx/fy/fz
         * variables, and only add it to memory forces[] at the end.
         */
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0index + ix) % ngrid[0];
            /* Get both the bspline factor and its derivative with respect to the x coordinate! */
            double tx  = thetax[ix];
            double dtx = dthetax[ix];

            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ngrid[1];
                /* bspline + derivative wrt y */
                double ty  = thetay[iy];
                double dty = dthetay[iy];

                for (int iz = 0; iz < order; iz++) {
                    /* Can be optimized, but we keep it simple here */
                    int zindex = (z0index + iz) % ngrid[2];
                    /* bspline + derivative wrt z */
                    double tz  = thetaz[iz];
                    double dtz = dthetaz[iz];
                    int index  = xindex*ngrid[1]*ngrid[2] + yindex*ngrid[2] + zindex;

                    /* Get the fft+convoluted+ifft:d data from the grid, which must be real by definition */
                    /* Checking that the imaginary part is indeed zero might be a good check :-) */
                    double gridvalue = grid[index].real();

                    /* The d component of the force is calculated by taking the derived bspline in dimension d, normal bsplines in the other two */
                    fx += dtx*ty*tz*gridvalue;
                    fy += tx*dty*tz*gridvalue;
                    fz += tx*ty*dtz*gridvalue;
                }
            }
        }
        /* Update memory force, note that we multiply by charge and some box stuff */
        forces[i][0] -= q*(fx*ngrid[0]*recipBoxVectors[0][0]);
        forces[i][1] -= q*(fx*ngrid[0]*recipBoxVectors[1][0]+fy*ngrid[1]*recipBoxVectors[1][1]);
        forces[i][2] -= q*(fx*ngrid[0]*recipBoxVectors[2][0]+fy*ngrid[1]*recipBoxVectors[2][1]+fz*ngrid[2]*recipBoxVectors[2][2]);
    }
}


void ReferencePME::grid_interpolate_charge_derivatives(const Vec3 recipBoxVectors[3], const vector<double>& charges,
            vector<double>& chargeDerivatives, const vector<int>& chargeIndices) {
    /* This is similar to grid_interpolate_force() */
    
    int nderiv = chargeIndices.size();
    for (int ideriv = 0; ideriv < nderiv; ideriv++) {
        int i = chargeIndices[ideriv];
        double dq = 0;

        /* Grid index for the actual atom position */
        int x0index = particleindex[i][0];
        int y0index = particleindex[i][1];
        int z0index = particleindex[i][2];

        /* Bspline factors for this atom in each dimension , calculated from fractional coordinates */
        double* thetax  = &(bsplines_theta[0][i*order]);
        double* thetay  = &(bsplines_theta[1][i*order]);
        double* thetaz  = &(bsplines_theta[2][i*order]);

        /* See grid_spread_charge() for comments about the order here, and only interpolation in one direction */

        /* Since we will add order^3 (typically 5*5*5=125) terms to the charge
         * derivative on each particle, we use a temporary dq variable, and only
         * add it to memory forces[] at the end.
         */
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0index + ix) % ngrid[0];
            /* Get the bspline factor with respect to the x coordinate */
            double tx  = thetax[ix];

            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ngrid[1];
                /* bspline wrt y */
                double ty  = thetay[iy];

                for (int iz = 0; iz < order; iz++) {
                    /* Can be optimized, but we keep it simple here */
                    int zindex = (z0index + iz) % ngrid[2];
                    /* bspline wrt z */
                    double tz  = thetaz[iz];
                    int index  = xindex*ngrid[1]*ngrid[2] + yindex*ngrid[2] + zindex;

                    /* Get the fft+convoluted+ifft:d data from the grid, which must be real by definition */
                    /* Checking that the imaginary part is indeed zero might be a good check :-) */
                    double gridvalue = grid[index].real();

                    /* The d component of the force is calculated by taking the derived bspline in dimension d, normal bsplines in the other two */
                    dq += tx*ty*tz*gridvalue;
                }
            }
        }

        chargeDerivatives[ideriv] += dq;
    }
}

ReferencePME::ReferencePME(double ewaldcoeff, int natoms, const int ngrid[3], int pme_order, double epsilon_r) :
                    order(pme_order), epsilon_r(epsilon_r), ewaldcoeff(ewaldcoeff), natoms(natoms) {
    for (int d = 0; d < 3; d++) {
        this->ngrid[d] = ngrid[d];
        bsplines_theta[d].resize(pme_order*natoms);
        bsplines_dtheta[d].resize(pme_order*natoms);
    }

    particlefraction.resize(natoms);
    particleindex.resize(natoms);

    /* Allocate charge grid storage */
    grid.resize(ngrid[0]*ngrid[1]*ngrid[2]);

    /* Setup bspline moduli (see Essman paper) */
    calculate_bsplines_moduli();
}

void ReferencePME::exec(const vector<Vec3>& atomCoordinates, vector<Vec3>& forces, const vector<double>& charges,
            const Vec3 periodicBoxVectors[3], double& energy) {
    /* Routine is called with coordinates in x, a box, and charges in q */

    Vec3 recipBoxVectors[3];
    ReferenceForce::invertBoxVectors(periodicBoxVectors, recipBoxVectors);
    
    /* Before we can do the actual interpolation, we need to recalculate and update
     * the indices for each particle in the charge grid (initialized in constructor),
     * and what its fractional offset in this grid cell is.
     */

    /* Update charge grid indices and fractional offsets for each atom.
     * The indices/fractions are stored internally in the pme datatype
     */
    update_grid_index_and_fraction(atomCoordinates, recipBoxVectors);

    /* Calculate bsplines (and their differentials) from current fractional coordinates, store in pme structure */
    update_bsplines();

    /* Spread the charges on grid (using newly calculated bsplines in the pme structure) */
    grid_spread_charge(charges);

    /* do 3d-fft */
    vector<size_t> shape = {(size_t) ngrid[0], (size_t) ngrid[1], (size_t) ngrid[2]};
    vector<size_t> axes = {0, 1, 2};
    vector<ptrdiff_t> stride = {(ptrdiff_t) (ngrid[1]*ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) (ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) sizeof(complex<double>)};
    pocketfft::c2c(shape, stride, stride, axes, true, grid.data(), grid.data(), 1.0, 0);

    /* solve in k-space */
    pme_reciprocal_convolution(periodicBoxVectors, recipBoxVectors, energy);

    /* do 3d-invfft */
    pocketfft::c2c(shape, stride, stride, axes, false, grid.data(), grid.data(), 1.0, 0);

    /* Get the particle forces from the grid and bsplines in the pme structure */
    grid_interpolate_force(recipBoxVectors, charges, forces);
}

void ReferencePME::exec_charge_derivatives(const vector<Vec3>& atomCoordinates, vector<double>& chargeDerivatives,
            const vector<int>& chargeIndices, const vector<double>& charges, const Vec3 periodicBoxVectors[3]) {
    /* Routine is called with coordinates in x, a box, and charges in q */

    Vec3 recipBoxVectors[3];
    ReferenceForce::invertBoxVectors(periodicBoxVectors, recipBoxVectors);
    
    /* Before we can do the actual interpolation, we need to recalculate and update
     * the indices for each particle in the charge grid (initialized in constructor),
     * and what its fractional offset in this grid cell is.
     */

    /* Update charge grid indices and fractional offsets for each atom.
     * The indices/fractions are stored internally in the pme datatype
     */
    update_grid_index_and_fraction(atomCoordinates, recipBoxVectors);

    /* Calculate bsplines (and their differentials) from current fractional coordinates, store in pme structure */
    update_bsplines();

    /* Spread the charges on grid (using newly calculated bsplines in the pme structure) */
    grid_spread_charge(charges);

    /* do 3d-fft */
    vector<size_t> shape = {(size_t) ngrid[0], (size_t) ngrid[1], (size_t) ngrid[2]};
    vector<size_t> axes = {0, 1, 2};
    vector<ptrdiff_t> stride = {(ptrdiff_t) (ngrid[1]*ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) (ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) sizeof(complex<double>)};
    pocketfft::c2c(shape, stride, stride, axes, true, grid.data(), grid.data(), 1.0, 0);

    /* solve in k-space */
    double energy;
    pme_reciprocal_convolution(periodicBoxVectors, recipBoxVectors, energy);

    /* do 3d-invfft */
    pocketfft::c2c(shape, stride, stride, axes, false, grid.data(), grid.data(), 1.0, 0);

    /* Get the charge derivatives from the grid and bsplines in the pme structure */
    grid_interpolate_charge_derivatives(recipBoxVectors, charges, chargeDerivatives, chargeIndices);
}

void ReferencePME::exec_dpme(const vector<Vec3>& atomCoordinates, vector<Vec3>& forces, const vector<double>& c6s,
            const Vec3 periodicBoxVectors[3], double& energy) {
    /* Routine is called with coordinates in x, a box, and charges in q */

    Vec3 recipBoxVectors[3];
    ReferenceForce::invertBoxVectors(periodicBoxVectors, recipBoxVectors);

    /* Before we can do the actual interpolation, we need to recalculate and update
     * the indices for each particle in the charge grid (initialized in constructor),
     * and what its fractional offset in this grid cell is.
     */

    /* Update charge grid indices and fractional offsets for each atom.
     * The indices/fractions are stored internally in the pme datatype
     */
    update_grid_index_and_fraction(atomCoordinates, recipBoxVectors);

    /* Calculate bsplines (and their differentials) from current fractional coordinates, store in pme structure */
    update_bsplines();

    /* Spread the charges on grid (using newly calculated bsplines in the pme structure) */
    grid_spread_charge(c6s);

    /* do 3d-fft */
    vector<size_t> shape = {(size_t) ngrid[0], (size_t) ngrid[1], (size_t) ngrid[2]};
    vector<size_t> axes = {0, 1, 2};
    vector<ptrdiff_t> stride = {(ptrdiff_t) (ngrid[1]*ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) (ngrid[2]*sizeof(complex<double>)),
                                (ptrdiff_t) sizeof(complex<double>)};
    pocketfft::c2c(shape, stride, stride, axes, true, grid.data(), grid.data(), 1.0, 0);

    /* solve in k-space */
    dpme_reciprocal_convolution(periodicBoxVectors, recipBoxVectors, energy);

    /* do 3d-invfft */
    pocketfft::c2c(shape, stride, stride, axes, false, grid.data(), grid.data(), 1.0, 0);

    /* Get the particle forces from the grid and bsplines in the pme structure */
    grid_interpolate_force(recipBoxVectors, c6s, forces);
}

} // namespace OpenMM
