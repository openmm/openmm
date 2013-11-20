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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ReferencePME.h"
#include "fftpack.h"

using std::vector;
using OpenMM::RealVec;

typedef int    ivec[3];


struct pme
{
    int          natoms;
    RealOpenMM       ewaldcoeff;

    t_complex *  grid;                 /* Memory for the grid we spread charges on.
                                        * Element (i,j,k) is accessed as:
                                        * grid[i*ngrid[1]*ngrid[2] + j*ngrid[2] + k]
                                        */
    int          ngrid[3];             /* Total grid dimensions (all data is complex!) */
    fftpack_t    fftplan;              /* Handle to fourier transform setup  */

    int          order;                /* PME interpolation order. Almost always 4 */

    /* Data for bspline interpolation, see the Essman PME paper */
    RealOpenMM *     bsplines_moduli[3];   /* 3 pointers, to x/y/z bspline moduli, each of length ngrid[x/y/z]   */
    RealOpenMM *     bsplines_theta[3];    /* each of x/y/z has length order*natoms */
    RealOpenMM *     bsplines_dtheta[3];   /* each of x/y/z has length order*natoms */

    ivec *       particleindex;        /* Array of length natoms. Each element is
                                        * an ivec (3 ints) that specify the grid
                                        * indices for that particular atom. Updated every step!
                                        */
    rvec *       particlefraction;     /* Array of length natoms. Fractional offset in the grid for
                                        * each atom in all three dimensions.
                                        */

    /* Further explanation of index/fraction:
     *
     * Assume we have a cell of size 10*10*10nm, and a total grid dimension of 100*100*100 cells.
     * In other words, each cell is 0.1*0.1*0.1 nm.
     *
     * If particle i has coordinates { 0.543 , 6.235 , -0.73 }, we will get:
     *
     * particleindex[i]    = { 5 , 62 , 92 }         (-0.73 + 10 = 9.27, we always apply PBC for grid calculations!)
     * particlefraction[i] = { 0.43 , 0.35 , 0.7 }   ( this is the fraction of the cell length where the atom is)
     *
     * (The reason for precaculating / storing these is that it gets a bit more complex for triclinic cells :-)
     *
     * In the current code version we might assume that a particle is not more than a whole box length away from
     * the central cell, i.e., in this case we would assume all coordinates fall in -10 nm < x,y,z < 20 nm.
     */

    RealOpenMM       epsilon_r;             /* Dielectric coefficient to use, typically 1.0 */
};


/* Internal setup routines */



/* Only called once from init_pme(), performance does not matter! */
static void
pme_calculate_bsplines_moduli(pme_t pme)
{
    int       nmax;
    int       i,j,k,l,d;
    int       order;
    int       ndata;
    RealOpenMM *  data;
    RealOpenMM *  ddata;
    RealOpenMM *  bsplines_data;
    RealOpenMM    div;
    RealOpenMM    sc,ss,arg;

    nmax = 0;
    for(d=0;d<3;d++)
    {
        nmax = (pme->ngrid[d] > nmax) ? pme->ngrid[d] : nmax;
        pme->bsplines_moduli[d] = (RealOpenMM *) malloc(sizeof(RealOpenMM)*pme->ngrid[d]);
    }

    order = pme->order;

    /* temp storage in this routine */
    data          = (RealOpenMM *) malloc(sizeof(RealOpenMM)*order);
    ddata         = (RealOpenMM *) malloc(sizeof(RealOpenMM)*order);
    bsplines_data = (RealOpenMM *) malloc(sizeof(RealOpenMM)*nmax);

    data[order-1]=0;
    data[1]=0;
    data[0]=1;

    for(k=3;k<order;k++)
    {
        div=(RealOpenMM) (1.0/(k-1.0));
        data[k-1]=0;
        for(l=1;l<(k-1);l++)
        {
            data[k-l-1]=div*(l*data[k-l-2]+(k-l)*data[k-l-1]);
        }
        data[0]=div*data[0];
    }

    /* differentiate */
    ddata[0]=-data[0];
    for(k=1;k<order;k++)
    {
        ddata[k]=data[k-1]-data[k];
    }

    div=(RealOpenMM) (1.0/(order-1));
    data[order-1]=0;

    for(l=1;l<(order-1);l++)
    {
        data[order-l-1]=div*(l*data[order-l-2]+(order-l)*data[order-l-1]);
    }
    data[0]=div*data[0];

    for(i=0;i<nmax;i++)
    {
        bsplines_data[i]=0;
    }
    for(i=1;i<=order;i++)
    {
        bsplines_data[i]=data[i-1];
    }

    /* Evaluate the actual bspline moduli for X/Y/Z */
    for(d=0;d<3;d++)
    {
        ndata = pme->ngrid[d];
        for(i=0;i<ndata;i++)
        {
            sc=ss=0;
            for(j=0;j<ndata;j++)
            {
                arg=(RealOpenMM) ((2.0*M_PI*i*j)/ndata);
                sc+=bsplines_data[j]*cos(arg);
                ss+=bsplines_data[j]*sin(arg);
            }
            pme->bsplines_moduli[d][i]=sc*sc+ss*ss;
        }
        for(i=0;i<ndata;i++)
        {
            if(pme->bsplines_moduli[d][i]<1.0e-7)
            {
                pme->bsplines_moduli[d][i]=(pme->bsplines_moduli[d][i-1]+pme->bsplines_moduli[d][i+1])/2;
            }
        }
    }

    /* Release temp storage */
    free(data);
    free(ddata);
    free(bsplines_data);
}



static void
pme_update_grid_index_and_fraction(pme_t    pme,
                                   const vector<RealVec>& atomCoordinates,
                                   const RealOpenMM   periodicBoxSize[3])
{
    int    i;
    int    d;
    RealOpenMM t;
    int    ti;

    for(i=0;i<pme->natoms;i++)
    {
        for(d=0;d<3;d++)
        {
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
            RealOpenMM coord = atomCoordinates[i][d];
            coord -= floor(coord/periodicBoxSize[d])*periodicBoxSize[d];
            t  = (coord / periodicBoxSize[d])*pme->ngrid[d];
            ti = (int) t;

            pme->particlefraction[i][d] = t - ti;
            pme->particleindex[i][d]    = ti % pme->ngrid[d];
        }
    }
}


/* Ugly bspline calculation taken from Tom Dardens reference equations.
 * This probably very sub-optimal in Cuda? Separate kernel?
 *
 * In practice, it might help to require order=4 for the cuda port.
 */
static void
pme_update_bsplines(pme_t    pme)
{
    int       i,j,k,l;
    int       order;
    RealOpenMM    dr,div;
    RealOpenMM *  data;
    RealOpenMM *  ddata;

    order = pme->order;

    for(i=0; (i<pme->natoms); i++)
    {
        for(j=0; j<3; j++)
        {
            /* dr is relative offset from lower cell limit */
            dr = pme->particlefraction[i][j];

            data  = &(pme->bsplines_theta[j][i*order]);
            ddata = &(pme->bsplines_dtheta[j][i*order]);
            data[order-1] = 0;
            data[1]       = dr;
            data[0]       = 1-dr;

            for(k=3; k<order; k++)
            {
                div = (RealOpenMM) (1.0/(k-1.0));
                data[k-1] = div*dr*data[k-2];
                for(l=1; l<(k-1); l++)
                {
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)*data[k-l-1]);
                }
                data[0] = div*(1-dr)*data[0];
            }

            /* differentiate */
            ddata[0] = -data[0];

            for(k=1; k<order; k++)
            {
                ddata[k] = data[k-1]-data[k];
            }

            div           = (RealOpenMM) (1.0/(order-1));
            data[order-1] = div*dr*data[order-2];

            for(l=1; l<(order-1); l++)
            {
                data[order-l-1] = div*((dr+l)*data[order-l-2]+(order-l-dr)*data[order-l-1]);
            }
            data[0] = div*(1-dr)*data[0];
        }
    }
}


static void
pme_grid_spread_charge(pme_t pme, const vector<RealOpenMM>& charges)
{
    int       order;
    int       i;
    int       ix,iy,iz;
    int       x0index,y0index,z0index;
    int       xindex,yindex,zindex;
    int       index;
    RealOpenMM    q;
    RealOpenMM *  thetax;
    RealOpenMM *  thetay;
    RealOpenMM *  thetaz;

    order = pme->order;

    /* Reset the grid */
    for(i=0;i<pme->ngrid[0]*pme->ngrid[1]*pme->ngrid[2];i++)
    {
        pme->grid[i].re = pme->grid[i].im = 0;
    }

    for(i=0;i<pme->natoms;i++)
    {
        q = charges[i];

        /* Grid index for the actual atom position */
        x0index = pme->particleindex[i][0];
        y0index = pme->particleindex[i][1];
        z0index = pme->particleindex[i][2];

        /* Bspline factors for this atom in each dimension , calculated from fractional coordinates */
        thetax  = &(pme->bsplines_theta[0][i*order]);
        thetay  = &(pme->bsplines_theta[1][i*order]);
        thetaz  = &(pme->bsplines_theta[2][i*order]);

        /* Loop over norder*norder*norder (typically 4*4*4) neighbor cells.
         *
         * As a neat optimization, we only spread in the forward direction, but apply PBC!
         *
         * Since we are going to do an FFT on the grid, it doesnt matter where the data is,
         * in frequency space the result will be the same.
         *
         * So, the influence function (bsplines) will probably be something like (0.15,0.35,0.35,0.15),
         * with largest weight 2-3 steps forward (you dont need to understand that for the implementation :-)
         * Effectively, you can look at this as translating the entire grid.
         *
         * Why do we do this stupid thing?
         *
         * 1) The loops get much simpler
         * 2) Just looking forward will hopefully get us more cache hits
         * 3) When we parallelize things, we only need to communicate in one direction instead of two!
         */

        for(ix=0;ix<order;ix++)
        {
            /* Calculate index, apply PBC so we spread to index 0/1/2 when a particle is close to the upper limit of the grid */
            xindex = (x0index + ix) % pme->ngrid[0];

            for(iy=0;iy<order;iy++)
            {
                yindex = (y0index + iy) % pme->ngrid[1];

                for(iz=0;iz<order;iz++)
                {
                    /* Can be optimized, but we keep it simple here */
                    zindex               = (z0index + iz) % pme->ngrid[2];
                    /* Calculate index in the charge grid */
                    index                = xindex*pme->ngrid[1]*pme->ngrid[2] + yindex*pme->ngrid[2] + zindex;
                    /* Add the charge times the bspline spread/interpolation factors to this grid position */
                    pme->grid[index].re += q*thetax[ix]*thetay[iy]*thetaz[iz];
                }
            }
        }
    }
}



static void
pme_reciprocal_convolution(pme_t     pme,
                           const RealOpenMM    periodicBoxSize[3],
                           RealOpenMM *  energy,
                           RealOpenMM    pme_virial[3][3])
{
    int kx,ky,kz;
    int nx,ny,nz;
    RealOpenMM mx,my,mz;
    RealOpenMM mhx,mhy,mhz,m2;
    RealOpenMM one_4pi_eps;
    RealOpenMM virxx,virxy,virxz,viryy,viryz,virzz;
    RealOpenMM bx,by,bz;
    RealOpenMM d1,d2;
    RealOpenMM eterm,vfactor,struct2,ets2;
    RealOpenMM esum;
    RealOpenMM factor;
    RealOpenMM denom;
    RealOpenMM boxfactor;
    RealOpenMM maxkx,maxky,maxkz;

    t_complex *ptr;

    nx = pme->ngrid[0];
    ny = pme->ngrid[1];
    nz = pme->ngrid[2];

    one_4pi_eps = (RealOpenMM) (ONE_4PI_EPS0/pme->epsilon_r);
    factor = (RealOpenMM) (M_PI*M_PI/(pme->ewaldcoeff*pme->ewaldcoeff));
    boxfactor = (RealOpenMM) (M_PI*periodicBoxSize[0]*periodicBoxSize[1]*periodicBoxSize[2]);

    esum = 0;
    virxx = 0;
    virxy = 0;
    virxz = 0;
    viryy = 0;
    viryz = 0;
    virzz = 0;

    maxkx = (RealOpenMM) ((nx+1)/2);
    maxky = (RealOpenMM) ((ny+1)/2);
    maxkz = (RealOpenMM) ((nz+1)/2);

    for(kx=0;kx<nx;kx++)
    {
        /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
        mx  = (RealOpenMM) ((kx<maxkx) ? kx : (kx-nx));
        mhx = mx/periodicBoxSize[0];
        bx  = boxfactor*pme->bsplines_moduli[0][kx];

        for(ky=0;ky<ny;ky++)
        {
            /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
            my  = (RealOpenMM) ((ky<maxky) ? ky : (ky-ny));
            mhy = my/periodicBoxSize[1];
            by  = pme->bsplines_moduli[1][ky];

            for(kz=0;kz<nz;kz++)
            {
                /* If the net charge of the system is 0.0, there will not be any DC (direct current, zero frequency) component. However,
                 * we can still handle charged systems through a charge correction, in which case the DC
                 * component should be excluded from recprocal space. We will anyway run into problems below when dividing with the
                 * frequency if it is zero...
                 *
                 * In cuda you could probably work around this by setting something to 0.0 instead, but the short story is that we
                 * should skip the zero frequency case!
                 */
                if (kx==0 && ky==0 && kz==0)
                {
                    continue;
                }

                /* Calculate frequency. Grid indices in the upper half correspond to negative frequencies! */
                mz        = (RealOpenMM) ((kz<maxkz) ? kz : (kz-nz));
                mhz       = mz/periodicBoxSize[2];

                /* Pointer to the grid cell in question */
                ptr       = pme->grid + kx*ny*nz + ky*nz + kz;

                /* Get grid data for this frequency */
                d1        = ptr->re;
                d2        = ptr->im;

                /* Calculate the convolution - see the Essman/Darden paper for the equation! */
                m2        = mhx*mhx+mhy*mhy+mhz*mhz;
                bz        = pme->bsplines_moduli[2][kz];
                denom     = m2*bx*by*bz;

                eterm     = one_4pi_eps*exp(-factor*m2)/denom;

                /* write back convolution data to grid */
                ptr->re   = d1*eterm;
                ptr->im   = d2*eterm;

                struct2   = (d1*d1+d2*d2);

                /* Long-range PME contribution to the energy for this frequency */
                ets2      = eterm*struct2;
                esum     += ets2;

                /* PME long-range contribution to atomic virial. Since it is symmetric, we only calculate half the matrix inside this loop. */
                vfactor   = (factor*m2+1)*2/m2;
                virxx    += ets2*(vfactor*mhx*mhx-1);
                virxy    += ets2*vfactor*mhx*mhy;
                virxz    += ets2*vfactor*mhx*mhz;
                viryy    += ets2*(vfactor*mhy*mhy-1);
                viryz    += ets2*vfactor*mhy*mhz;
                virzz    += ets2*(vfactor*mhz*mhz-1);
            }
        }
    }
    pme_virial[0][0] = (RealOpenMM) (0.25*virxx);
    pme_virial[1][1] = (RealOpenMM) (0.25*viryy);
    pme_virial[2][2] = (RealOpenMM) (0.25*virzz);
    pme_virial[0][1] = pme_virial[1][0] = (RealOpenMM) (0.25*virxy);
    pme_virial[0][2] = pme_virial[2][0] = (RealOpenMM) (0.25*virxz);
    pme_virial[1][2] = pme_virial[2][1] = (RealOpenMM) (0.25*viryz);

    /* The factor 0.5 is nothing special, but it is better to have it here than inside the loop :-) */
    *energy = (RealOpenMM) (0.5*esum);
}


static void
pme_grid_interpolate_force(pme_t pme,
                           const RealOpenMM periodicBoxSize[3],
                           const vector<RealOpenMM>& charges,
                           vector<RealVec>& forces)
{
    int       i;
    int       ix,iy,iz;
    int       x0index,y0index,z0index;
    int       xindex,yindex,zindex;
    int       index;
    int       order;
    RealOpenMM    q;
    RealOpenMM *  thetax;
    RealOpenMM *  thetay;
    RealOpenMM *  thetaz;
    RealOpenMM *  dthetax;
    RealOpenMM *  dthetay;
    RealOpenMM *  dthetaz;
    RealOpenMM    tx,ty,tz;
    RealOpenMM    dtx,dty,dtz;
    RealOpenMM    fx,fy,fz;
    RealOpenMM    gridvalue;
    int       nx,ny,nz;

    nx    = pme->ngrid[0];
    ny    = pme->ngrid[1];
    nz    = pme->ngrid[2];

    order = pme->order;

    /* This is almost identical to the charge spreading routine! */

    for(i=0;i<pme->natoms;i++)
    {
        fx = fy = fz = 0;

        q = charges[i];

        /* Grid index for the actual atom position */
        x0index = pme->particleindex[i][0];
        y0index = pme->particleindex[i][1];
        z0index = pme->particleindex[i][2];

        /* Bspline factors for this atom in each dimension , calculated from fractional coordinates */
        thetax  = &(pme->bsplines_theta[0][i*order]);
        thetay  = &(pme->bsplines_theta[1][i*order]);
        thetaz  = &(pme->bsplines_theta[2][i*order]);
        dthetax = &(pme->bsplines_dtheta[0][i*order]);
        dthetay = &(pme->bsplines_dtheta[1][i*order]);
        dthetaz = &(pme->bsplines_dtheta[2][i*order]);

        /* See pme_grid_spread_charge() for comments about the order here, and only interpolation in one direction */

        /* Since we will add order^3 (typically 4*4*4=64) terms to the force on each particle, we use temporary fx/fy/fz
         * variables, and only add it to memory forces[] at the end.
         */
        for(ix=0;ix<order;ix++)
        {
            xindex = (x0index + ix) % pme->ngrid[0];
            /* Get both the bspline factor and its derivative with respect to the x coordinate! */
            tx     = thetax[ix];
            dtx    = dthetax[ix];

            for(iy=0;iy<order;iy++)
            {
                yindex = (y0index + iy) % pme->ngrid[1];
                /* bspline + derivative wrt y */
                ty     = thetay[iy];
                dty    = dthetay[iy];

                for(iz=0;iz<order;iz++)
                {
                    /* Can be optimized, but we keep it simple here */
                    zindex               = (z0index + iz) % pme->ngrid[2];
                    /* bspline + derivative wrt z */
                    tz                   = thetaz[iz];
                    dtz                  = dthetaz[iz];
                    index                = xindex*pme->ngrid[1]*pme->ngrid[2] + yindex*pme->ngrid[2] + zindex;

                    /* Get the fft+convoluted+ifft:d data from the grid, which must be real by definition */
                    /* Checking that the imaginary part is indeed zero might be a good check :-) */
                    gridvalue            = pme->grid[index].re;

                    /* The d component of the force is calculated by taking the derived bspline in dimension d, normal bsplines in the other two */
                    fx                  += dtx*ty*tz*gridvalue;
                    fy                  += tx*dty*tz*gridvalue;
                    fz                  += tx*ty*dtz*gridvalue;
                }
            }
        }
        /* Update memory force, note that we multiply by charge and some box stuff */
        forces[i][0] -= q*fx*nx/periodicBoxSize[0];
        forces[i][1] -= q*fy*ny/periodicBoxSize[1];
        forces[i][2] -= q*fz*nz/periodicBoxSize[2];
    }
}



/* EXPORTED ROUTINES */

int
pme_init(pme_t *       ppme,
         RealOpenMM        ewaldcoeff,
         int           natoms,
         const int           ngrid[3],
         int           pme_order,
         RealOpenMM        epsilon_r)
{
    pme_t pme;
    int   d;

    pme = (pme_t) malloc(sizeof(struct pme));

    pme->order       = pme_order;
    pme->epsilon_r   = epsilon_r;
    pme->ewaldcoeff  = ewaldcoeff;
    pme->natoms      = natoms;

    for(d=0;d<3;d++)
    {
        pme->ngrid[d]            = ngrid[d];
        pme->bsplines_theta[d]   = (RealOpenMM *)malloc(sizeof(RealOpenMM)*pme_order*natoms);
        pme->bsplines_dtheta[d]  = (RealOpenMM *)malloc(sizeof(RealOpenMM)*pme_order*natoms);
    }

    pme->particlefraction = (rvec *)malloc(sizeof(rvec)*natoms);
    pme->particleindex    = (ivec *)malloc(sizeof(ivec)*natoms);

    /* Allocate charge grid storage */
    pme->grid        = (t_complex *)malloc(sizeof(t_complex)*ngrid[0]*ngrid[1]*ngrid[2]);

    fftpack_init_3d(&pme->fftplan,ngrid[0],ngrid[1],ngrid[2]);

    /* Setup bspline moduli (see Essman paper) */
    pme_calculate_bsplines_moduli(pme);

    *ppme = pme;

    return 0;
}





int pme_exec(pme_t       pme,
             const vector<RealVec>& atomCoordinates,
             vector<RealVec>& forces,
             const vector<RealOpenMM>& charges,
             const RealOpenMM periodicBoxSize[3],
             RealOpenMM* energy,
             RealOpenMM pme_virial[3][3])
{
    /* Routine is called with coordinates in x, a box, and charges in q */

    /* Before we can do the actual interpolation, we need to recalculate and update
     * the indices for each particle in the charge grid (initialized in pme_init()),
     * and what its fractional offset in this grid cell is.
     */

    /* Update charge grid indices and fractional offsets for each atom.
     * The indices/fractions are stored internally in the pme datatype
     */
    pme_update_grid_index_and_fraction(pme,atomCoordinates,periodicBoxSize);

    /* Calculate bsplines (and their differentials) from current fractional coordinates, store in pme structure */
    pme_update_bsplines(pme);

    /* Spread the charges on grid (using newly calculated bsplines in the pme structure) */
    pme_grid_spread_charge(pme, charges);

    /* do 3d-fft */
    fftpack_exec_3d(pme->fftplan,FFTPACK_FORWARD,pme->grid,pme->grid);

    /* solve in k-space */
    pme_reciprocal_convolution(pme,periodicBoxSize,energy,pme_virial);

    /* do 3d-invfft */
    fftpack_exec_3d(pme->fftplan,FFTPACK_BACKWARD,pme->grid,pme->grid);

    /* Get the particle forces from the grid and bsplines in the pme structure */
    pme_grid_interpolate_force(pme,periodicBoxSize,charges,forces);

    return 0;
}



int
pme_destroy(pme_t    pme)
{
    int d;

    free(pme->grid);

    for(d=0;d<3;d++)
    {
        free(pme->bsplines_moduli[d]);
        free(pme->bsplines_theta[d]);
        free(pme->bsplines_dtheta[d]);
    }

    free(pme->particlefraction);
    free(pme->particleindex);

    fftpack_destroy(pme->fftplan);

    /* destroy structure itself */
    free(pme);

    return 0;
}
