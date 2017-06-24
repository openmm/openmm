/*
* This file contains a Fortran to C translation of the 1D transformations
* based on the original FFTPACK, written by paul n swarztrauber
* at the national center for atmospheric research and available
* at www.netlib.org. FFTPACK is in the public domain.
*
* Higher-dimension transforms written by Erik Lindahl, 2008-2009.
* Just as FFTPACK, this file may be redistributed freely, and can be
* considered to be in the public domain. 
*
* Any errors in this (threadsafe, but not threaded) C version
* are probably due to the f2c translator, or hacks by Erik Lindahl,
* rather than FFTPACK. If you find a bug, it would be great if you could
* drop a line to lindahl@cbr.su.se and let me know about it!
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>


#include "fftpack.h"


/** Contents of the FFTPACK fft datatype.
 *
 *  FFTPACK only does 1d transforms, so we use a pointers to another fft for
 *  the transform in the next dimension.
 * Thus, a 3d-structure contains a pointer to a 2d one, which in turns contains
 * a pointer to a 1d. The 1d structure has next==NULL.
 */
struct fftpack
{
    int                    ndim;     /**< Dimensions, including our subdimensions.  */
    int                    n;        /**< Number of points in this dimension.       */
    int                    ifac[15]; /**< 15 bytes needed for cfft and rfft         */
    struct fftpack *       next;     /**< Pointer to next dimension, or NULL.       */
    double *               work;     /**< 1st 4n reserved for cfft, 1st 2n for rfft */
};








static void
fftpack_passf2(int    ido,
               int    l1,
               double cc[],
               double ch[],
               double wa1[],
               int    isign)
{
    int i, k, ah, ac;
    double ti2, tr2;

    if (ido <= 2)
    {
        for (k=0; k<l1; k++)
        {
            ah = k*ido;
            ac = 2*k*ido;
            ch[ah]              = cc[ac]   + cc[ac + ido];
            ch[ah + ido*l1]     = cc[ac]   - cc[ac + ido];
            ch[ah+1]            = cc[ac+1] + cc[ac + ido + 1];
            ch[ah + ido*l1 + 1] = cc[ac+1] - cc[ac + ido + 1];
        }
    }
    else
    {
        for (k=0; k<l1; k++)
        {
            for (i=0; i<ido-1; i+=2)
            {
                ah              = i + k*ido;
                ac              = i + 2*k*ido;
                ch[ah]          = cc[ac] + cc[ac + ido];
                tr2             = cc[ac] - cc[ac + ido];
                ch[ah+1]        = cc[ac+1] + cc[ac + 1 + ido];
                ti2             = cc[ac+1] - cc[ac + 1 + ido];
                ch[ah+l1*ido+1] = wa1[i]*ti2 + isign*wa1[i+1]*tr2;
                ch[ah+l1*ido]   = wa1[i]*tr2 - isign*wa1[i+1]*ti2;
            }
        }
    }
}



static void
fftpack_passf3(int    ido,
               int    l1,
               double cc[],
               double ch[],
               double wa1[],
               double wa2[],
               int    isign)
{
    const double taur = -0.5;
    const double taui = 0.866025403784439;

    int i, k, ac, ah;
    double ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

    if (ido == 2)
    {
        for (k=1; k<=l1; k++)
        {
            ac = (3*k - 2)*ido;
            tr2 = cc[ac] + cc[ac + ido];
            cr2 = cc[ac - ido] + taur*tr2;
            ah = (k - 1)*ido;
            ch[ah] = cc[ac - ido] + tr2;

            ti2 = cc[ac + 1] + cc[ac + ido + 1];
            ci2 = cc[ac - ido + 1] + taur*ti2;
            ch[ah + 1] = cc[ac - ido + 1] + ti2;

            cr3 = isign*taui*(cc[ac] - cc[ac + ido]);
            ci3 = isign*taui*(cc[ac + 1] - cc[ac + ido + 1]);
            ch[ah + l1*ido] = cr2 - ci3;
            ch[ah + 2*l1*ido] = cr2 + ci3;
            ch[ah + l1*ido + 1] = ci2 + cr3;
            ch[ah + 2*l1*ido + 1] = ci2 - cr3;
        }
    }
    else
    {
        for (k=1; k<=l1; k++)
        {
            for (i=0; i<ido-1; i+=2)
            {
                ac = i + (3*k - 2)*ido;
                tr2 = cc[ac] + cc[ac + ido];
                cr2 = cc[ac - ido] + taur*tr2;
                ah = i + (k-1)*ido;
                ch[ah] = cc[ac - ido] + tr2;
                ti2 = cc[ac + 1] + cc[ac + ido + 1];
                ci2 = cc[ac - ido + 1] + taur*ti2;
                ch[ah + 1] = cc[ac - ido + 1] + ti2;
                cr3 = isign*taui*(cc[ac] - cc[ac + ido]);
                ci3 = isign*taui*(cc[ac + 1] - cc[ac + ido + 1]);
                dr2 = cr2 - ci3;
                dr3 = cr2 + ci3;
                di2 = ci2 + cr3;
                di3 = ci2 - cr3;
                ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
                ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
                ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
                ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
            }
        }
    }
}


static void
fftpack_passf4(int    ido,
               int    l1,
               double cc[],
               double ch[],
               double wa1[],
               double wa2[],
               double wa3[],
               int    isign)
{
    int i, k, ac, ah;
    double ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

    if (ido == 2)
    {
        for (k=0; k<l1; k++)
        {
            ac = 4*k*ido + 1;
            ti1 = cc[ac] - cc[ac + 2*ido];
            ti2 = cc[ac] + cc[ac + 2*ido];
            tr4 = cc[ac + 3*ido] - cc[ac + ido];
            ti3 = cc[ac + ido] + cc[ac + 3*ido];
            tr1 = cc[ac - 1] - cc[ac + 2*ido - 1];
            tr2 = cc[ac - 1] + cc[ac + 2*ido - 1];
            ti4 = cc[ac + ido - 1] - cc[ac + 3*ido - 1];
            tr3 = cc[ac + ido - 1] + cc[ac + 3*ido - 1];
            ah = k*ido;
            ch[ah] = tr2 + tr3;
            ch[ah + 2*l1*ido] = tr2 - tr3;
            ch[ah + 1] = ti2 + ti3;
            ch[ah + 2*l1*ido + 1] = ti2 - ti3;
            ch[ah + l1*ido] = tr1 + isign*tr4;
            ch[ah + 3*l1*ido] = tr1 - isign*tr4;
            ch[ah + l1*ido + 1] = ti1 + isign*ti4;
            ch[ah + 3*l1*ido + 1] = ti1 - isign*ti4;
        }
    }
    else
    {
        for (k=0; k<l1; k++)
        {
            for (i=0; i<ido-1; i+=2)
            {
                ac = i + 1 + 4*k*ido;
                ti1 = cc[ac] - cc[ac + 2*ido];
                ti2 = cc[ac] + cc[ac + 2*ido];
                ti3 = cc[ac + ido] + cc[ac + 3*ido];
                tr4 = cc[ac + 3*ido] - cc[ac + ido];
                tr1 = cc[ac - 1] - cc[ac + 2*ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 2*ido - 1];
                ti4 = cc[ac + ido - 1] - cc[ac + 3*ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 3*ido - 1];
                ah = i + k*ido;
                ch[ah] = tr2 + tr3;
                cr3 = tr2 - tr3;
                ch[ah + 1] = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 + isign*tr4;
                cr4 = tr1 - isign*tr4;
                ci2 = ti1 + isign*ti4;
                ci4 = ti1 - isign*ti4;
                ch[ah + l1*ido] = wa1[i]*cr2 - isign*wa1[i + 1]*ci2;
                ch[ah + l1*ido + 1] = wa1[i]*ci2 + isign*wa1[i + 1]*cr2;
                ch[ah + 2*l1*ido] = wa2[i]*cr3 - isign*wa2[i + 1]*ci3;
                ch[ah + 2*l1*ido + 1] = wa2[i]*ci3 + isign*wa2[i + 1]*cr3;
                ch[ah + 3*l1*ido] = wa3[i]*cr4 -isign*wa3[i + 1]*ci4;
                ch[ah + 3*l1*ido + 1] = wa3[i]*ci4 + isign*wa3[i + 1]*cr4;
            }
        }
    }
}


static void
fftpack_passf5(int    ido,
               int    l1,
               double cc[],
               double ch[],
               double wa1[],
               double wa2[],
               double wa3[],
               double wa4[],
               int    isign)
{
    const double tr11 = 0.309016994374947;
    const double ti11 = 0.951056516295154;
    const double tr12 = -0.809016994374947;
    const double ti12 = 0.587785252292473;

    int i, k, ac, ah;
    double ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
        ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

    if (ido == 2)
    {
        for (k = 1; k <= l1; ++k)
        {
            ac = (5*k - 4)*ido + 1;
            ti5 = cc[ac] - cc[ac + 3*ido];
            ti2 = cc[ac] + cc[ac + 3*ido];
            ti4 = cc[ac + ido] - cc[ac + 2*ido];
            ti3 = cc[ac + ido] + cc[ac + 2*ido];
            tr5 = cc[ac - 1] - cc[ac + 3*ido - 1];
            tr2 = cc[ac - 1] + cc[ac + 3*ido - 1];
            tr4 = cc[ac + ido - 1] - cc[ac + 2*ido - 1];
            tr3 = cc[ac + ido - 1] + cc[ac + 2*ido - 1];
            ah = (k - 1)*ido;
            ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
            ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
            cr2 = cc[ac - ido - 1] + tr11*tr2 + tr12*tr3;
            ci2 = cc[ac - ido] + tr11*ti2 + tr12*ti3;
            cr3 = cc[ac - ido - 1] + tr12*tr2 + tr11*tr3;
            ci3 = cc[ac - ido] + tr12*ti2 + tr11*ti3;
            cr5 = isign*(ti11*tr5 + ti12*tr4);
            ci5 = isign*(ti11*ti5 + ti12*ti4);
            cr4 = isign*(ti12*tr5 - ti11*tr4);
            ci4 = isign*(ti12*ti5 - ti11*ti4);
            ch[ah + l1*ido] = cr2 - ci5;
            ch[ah + 4*l1*ido] = cr2 + ci5;
            ch[ah + l1*ido + 1] = ci2 + cr5;
            ch[ah + 2*l1*ido + 1] = ci3 + cr4;
            ch[ah + 2*l1*ido] = cr3 - ci4;
            ch[ah + 3*l1*ido] = cr3 + ci4;
            ch[ah + 3*l1*ido + 1] = ci3 - cr4;
            ch[ah + 4*l1*ido + 1] = ci2 - cr5;
        }
    }
    else
    {
        for (k=1; k<=l1; k++)
        {
            for (i=0; i<ido-1; i+=2)
            {
                ac = i + 1 + (k*5 - 4)*ido;
                ti5 = cc[ac] - cc[ac + 3*ido];
                ti2 = cc[ac] + cc[ac + 3*ido];
                ti4 = cc[ac + ido] - cc[ac + 2*ido];
                ti3 = cc[ac + ido] + cc[ac + 2*ido];
                tr5 = cc[ac - 1] - cc[ac + 3*ido - 1];
                tr2 = cc[ac - 1] + cc[ac + 3*ido - 1];
                tr4 = cc[ac + ido - 1] - cc[ac + 2*ido - 1];
                tr3 = cc[ac + ido - 1] + cc[ac + 2*ido - 1];
                ah = i + (k - 1)*ido;
                ch[ah] = cc[ac - ido - 1] + tr2 + tr3;
                ch[ah + 1] = cc[ac - ido] + ti2 + ti3;
                cr2 = cc[ac - ido - 1] + tr11*tr2 + tr12*tr3;
                ci2 = cc[ac - ido] + tr11*ti2 + tr12*ti3;
                cr3 = cc[ac - ido - 1] + tr12*tr2 + tr11*tr3;
                ci3 = cc[ac - ido] + tr12*ti2 + tr11*ti3;
                cr5 = isign*(ti11*tr5 + ti12*tr4);
                ci5 = isign*(ti11*ti5 + ti12*ti4);
                cr4 = isign*(ti12*tr5 - ti11*tr4);
                ci4 = isign*(ti12*ti5 - ti11*ti4);
                dr3 = cr3 - ci4;
                dr4 = cr3 + ci4;
                di3 = ci3 + cr4;
                di4 = ci3 - cr4;
                dr5 = cr2 + ci5;
                dr2 = cr2 - ci5;
                di5 = ci2 - cr5;
                di2 = ci2 + cr5;
                ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
                ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
                ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
                ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
                ch[ah + 3*l1*ido] = wa3[i]*dr4 - isign*wa3[i+1]*di4;
                ch[ah + 3*l1*ido + 1] = wa3[i]*di4 + isign*wa3[i+1]*dr4;
                ch[ah + 4*l1*ido] = wa4[i]*dr5 - isign*wa4[i+1]*di5;
                ch[ah + 4*l1*ido + 1] = wa4[i]*di5 + isign*wa4[i+1]*dr5;
            }
        }
    }
}


static void
fftpack_passf(int*   nac,
              int    ido,
              int    ip,
              int    l1,
              int    idl1,
              double cc[],
              double ch[],
              double wa[],
              int    isign)
{
    int idij, idlj, idot, ipph, i, j, k, l, jc, lc, ik, nt, idj, idl, inc,idp;
    double wai, war;

    idot = ido / 2;
    nt = ip*idl1;
    ipph = (ip + 1) / 2;
    idp = ip*ido;
    if (ido >= l1)
    {
        for (j=1; j<ipph; j++)
        {
            jc = ip - j;
            for (k=0; k<l1; k++)
            {
                for (i=0; i<ido; i++)
                {
                    ch[i + (k + j*l1)*ido]  = cc[i + (j + k*ip)*ido] + cc[i + (jc + k*ip)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (j + k*ip)*ido] - cc[i + (jc + k*ip)*ido];
                }
            }
        }
        for (k=0; k<l1; k++)
            for (i=0; i<ido; i++)
                ch[i + k*ido] = cc[i + k*ip*ido];
    }
    else
    {
        for (j=1; j<ipph; j++)
        {
            jc = ip - j;
            for (i=0; i<ido; i++)
            {
                for (k=0; k<l1; k++)
                {
                    ch[i + (k + j*l1)*ido] =  cc[i + (j + k*ip)*ido] + cc[i + (jc + k*ip)*ido];
                    ch[i + (k + jc*l1)*ido] = cc[i + (j + k*ip)*ido] - cc[i + (jc + k*ip)*ido];
                }
            }
        }
        for (i=0; i<ido; i++)
            for (k=0; k<l1; k++)
                ch[i + k*ido] = cc[i + k*ip*ido];
    }

    idl = 2 - ido;
    inc = 0;
    for (l=1; l<ipph; l++)
    {
        lc = ip - l;
        idl += ido;
        for (ik=0; ik<idl1; ik++)
        {
            cc[ik + l*idl1] = ch[ik] + wa[idl - 2]*ch[ik + idl1];
            cc[ik + lc*idl1] = isign*wa[idl-1]*ch[ik + (ip-1)*idl1];
        }
        idlj = idl;
        inc += ido;
        for (j=2; j<ipph; j++)
        {
            jc = ip - j;
            idlj += inc;
            if (idlj > idp) idlj -= idp;
            war = wa[idlj - 2];
            wai = wa[idlj-1];
            for (ik=0; ik<idl1; ik++)
            {
                cc[ik + l*idl1] += war*ch[ik + j*idl1];
                cc[ik + lc*idl1] += isign*wai*ch[ik + jc*idl1];
            }
        }
    }
    for (j=1; j<ipph; j++)
        for (ik=0; ik<idl1; ik++)
            ch[ik] += ch[ik + j*idl1];
    for (j=1; j<ipph; j++)
    {
        jc = ip - j;
        for (ik=1; ik<idl1; ik+=2)
        {
            ch[ik - 1 + j*idl1] = cc[ik - 1 + j*idl1] - cc[ik + jc*idl1];
            ch[ik - 1 + jc*idl1] = cc[ik - 1 + j*idl1] + cc[ik + jc*idl1];
            ch[ik + j*idl1] = cc[ik + j*idl1] + cc[ik - 1 + jc*idl1];
            ch[ik + jc*idl1] = cc[ik + j*idl1] - cc[ik - 1 + jc*idl1];
        }
    }
    *nac = 1;
    if (ido == 2)
        return;
    *nac = 0;
    for (ik=0; ik<idl1; ik++)
    {
        cc[ik] = ch[ik];
    }
    for (j=1; j<ip; j++)
    {
        for (k=0; k<l1; k++)
        {
            cc[(k + j*l1)*ido + 0] = ch[(k + j*l1)*ido + 0];
            cc[(k + j*l1)*ido + 1] = ch[(k + j*l1)*ido + 1];
        }
    }
    if (idot <= l1)
    {
        idij = 0;
        for (j=1; j<ip; j++)
        {
            idij += 2;
            for (i=3; i<ido; i+=2)
            {
                idij += 2;
                for (k=0; k<l1; k++)
                {
                    cc[i - 1 + (k + j*l1)*ido] =
                    wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                    isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                        isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
    else
    {
        idj = 2 - ido;
        for (j=1; j<ip; j++)
        {
            idj += ido;
            for (k = 0; k < l1; k++)
            {
                idij = idj;
                for (i=3; i<ido; i+=2)
                {
                    idij += 2;
                    cc[i - 1 + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                        isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
                    cc[i + (k + j*l1)*ido] =
                        wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                        isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
                }
            }
        }
    }
}





static void
fftpack_cfftf1(int    n,
               double c[],
               double ch[],
               double wa[],
               int    ifac[15],
               int    isign)
{
    int idot, i;
    int k1, l1, l2;
    int na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;
    double *cinput, *coutput;
    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 0;

    for (k1=2; k1<=nf+1; k1++)
    {
        ip = ifac[k1];
        l2 = ip*l1;
        ido = n / l2;
        idot = ido + ido;
        idl1 = idot*l1;
        if (na)
        {
            cinput = ch;
            coutput = c;
        }
        else
        {
            cinput = c;
            coutput = ch;
        }
        switch (ip)
        {
            case 4:
                ix2 = iw + idot;
                ix3 = ix2 + idot;
                fftpack_passf4(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], isign);
                na = !na;
                break;
            case 2:
                fftpack_passf2(idot, l1, cinput, coutput, &wa[iw], isign);
                na = !na;
                break;
            case 3:
                ix2 = iw + idot;
                fftpack_passf3(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], isign);
                na = !na;
                break;
            case 5:
                ix2 = iw + idot;
                ix3 = ix2 + idot;
                ix4 = ix3 + idot;
                fftpack_passf5(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4], isign);
                na = !na;
                break;
            default:
                fftpack_passf(&nac, idot, ip, l1, idl1, cinput, coutput, &wa[iw], isign);
                if (nac != 0) na = !na;
        }
        l1 = l2;
        iw += (ip - 1)*idot;
    }
    if (na == 0)
        return;
    for (i=0; i<2*n; i++)
        c[i] = ch[i];
}




static void
fftpack_factorize(int    n,
                  int    ifac[15])
{
    static const int ntryh[4] = { 3,4,2,5 };
    int ntry=3, i, j=0, ib, nf=0, nl=n, nq, nr;

startloop:
    if (j < 4)
        ntry = ntryh[j];
    else
        ntry+= 2;
    j++;
    do
    {
        nq = nl / ntry;
        nr = nl - ntry*nq;
        if (nr != 0) goto startloop;
        nf++;
        ifac[nf + 1] = ntry;
        nl = nq;
        if (ntry == 2 && nf != 1)
        {
            for (i=2; i<=nf; i++)
            {
                ib = nf - i + 2;
                ifac[ib + 1] = ifac[ib];
            }
            ifac[2] = 2;
        }
    }
    while (nl != 1);
    ifac[0] = n;
    ifac[1] = nf;
}


static void
fftpack_cffti1(int    n,
               double wa[],
               int    ifac[15])
{
    const double twopi = 6.28318530717959;
    double arg, argh, argld, fi;
    int idot, i, j;
    int i1, k1, l1, l2;
    int ld, ii, nf, ip;
    int ido, ipm;

    fftpack_factorize(n,ifac);
    nf = ifac[1];
    argh = twopi/n;
    i = 1;
    l1 = 1;
    for (k1=1; k1<=nf; k1++)
    {
        ip = ifac[k1+1];
        ld = 0;
        l2 = l1*ip;
        ido = n / l2;
        idot = ido + ido + 2;
        ipm = ip - 1;
        for (j=1; j<=ipm; j++)
        {
            i1 = i;
            wa[i-1] = 1;
            wa[i] = 0;
            ld += l1;
            fi = 0;
            argld = ld*argh;
            for (ii=4; ii<=idot; ii+=2)
            {
                i+= 2;
                fi+= 1;
                arg = fi*argld;
                wa[i-1] = cos(arg);
                wa[i] = sin(arg);
            }
            if (ip > 5)
            {
                wa[i1-1] = wa[i-1];
                wa[i1] = wa[i];
            }
        }
        l1 = l2;
    }
}








static int
fftpack_transpose_2d(t_complex *          in_data,
                     t_complex *          out_data,
                     int                  nx,
                     int                  ny)
{
    t_complex *  src;
    int          i,j;

    if (nx<2 || ny<2)
    {
        if (in_data != out_data)
        {
            memcpy(out_data,in_data,sizeof(t_complex)*nx*ny);
        }
        return 0;
    }

    if (in_data == out_data)
    {
        src = (t_complex *)malloc(sizeof(t_complex)*nx*ny);
        memcpy(src,in_data,sizeof(t_complex)*nx*ny);
    }
    else
    {
        src = in_data;
    }

    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            out_data[j*nx+i].re = src[i*ny+j].re;
            out_data[j*nx+i].im = src[i*ny+j].im;
        }
    }

    if (src != in_data)
    {
        free(src);
    }
    return 0;
}



static int
fftpack_transpose_2d_nelem(t_complex *          in_data,
                           t_complex *          out_data,
                           int                  nx,
                           int                  ny,
                           int                  nelem)
{
    t_complex *   src;
    int           ncopy;
    int           i,j;

    ncopy = nelem*sizeof(t_complex);

    if (nx<2 || ny<2)
    {
        if (in_data != out_data)
        {
            memcpy(out_data,in_data,nx*ny*ncopy);
        }
        return 0;
    }

    if (in_data == out_data)
    {
        src = (t_complex *)malloc(nx*ny*ncopy);
        memcpy(src,in_data,nx*ny*ncopy);
    }
    else
    {
        src = in_data;
    }

    for (i=0;i<nx;i++)
    {
        for (j=0;j<ny;j++)
        {
            memcpy(out_data + (j*nx+i)*nelem , src + (i*ny+j)*nelem , ncopy);
        }
    }

    if (src != in_data)
    {
        free(src);
    }
    return 0;
}






int
fftpack_init_1d(fftpack_t *        pfft,
                int                nx)
{
    fftpack_t    fft;

    if (pfft==NULL)
    {
        fprintf(stderr,"Fatal error - Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    if ((fft = (struct fftpack *)malloc(sizeof(struct fftpack))) == NULL)
    {
        return ENOMEM;
    }

    fft->next = NULL;
    fft->n    = nx;

    /* Need 4*n storage for 1D complex FFT */
    if ((fft->work = (double *)malloc(sizeof(double)*(4*nx))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }

    if (fft->n>1)
        fftpack_cffti1(nx,fft->work,fft->ifac);

    *pfft = fft;
    return 0;
};




int
fftpack_init_2d(fftpack_t *        pfft,
                int                nx,
                int                ny)
{
    fftpack_t     fft;
    int           rc;

    if (pfft==NULL)
    {
        fprintf(stderr,"Fatal error - Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */
    if ((rc = fftpack_init_1d(&fft,nx)) != 0)
    {
        return rc;
    }

    /* Create Y transform as a link from X */
    if ((rc=fftpack_init_1d(&(fft->next),ny)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;
    return 0;
};



int
fftpack_init_3d(fftpack_t *        pfft,
                int                nx,
                int                ny,
                int                nz)
{
    fftpack_t     fft;
    int           rc;

    if (pfft==NULL)
    {
        fprintf(stderr,"Fatal error - Invalid FFT opaque type pointer.");
        return EINVAL;
    }
    *pfft = NULL;

    /* Create the X transform */

    if ((fft = (struct fftpack *)malloc(sizeof(struct fftpack))) == NULL)
    {
        return ENOMEM;
    }

    fft->n    = nx;

    /* Need 4*nx storage for 1D complex FFT.
     */
    if ((fft->work = (double *)malloc(sizeof(double)*(4*nx))) == NULL)
    {
        free(fft);
        return ENOMEM;
    }

    fftpack_cffti1(nx,fft->work,fft->ifac);

    /* Create 2D Y/Z transforms as a link from X */
    if ((rc=fftpack_init_2d(&(fft->next),ny,nz)) != 0)
    {
        free(fft);
        return rc;
    }

    *pfft = fft;

    return 0;
};


int
fftpack_exec_1d          (fftpack_t                  fft,
                          enum fftpack_direction     dir,
                          t_complex *                in_data,
                          t_complex *                out_data)
{
    int i, n;
    double* p1;
    double* p2;

    n=fft->n;

    if (n==1)
    {
        p1 = (double *)in_data;
        p2 = (double *)out_data;
        p2[0] = p1[0];
        p2[1] = p1[1];
    }

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     */
    if (in_data != out_data)
    {
        p1 = (double *)in_data;
        p2 = (double *)out_data;

        /* n complex = 2*n double elements */
        for (i=0;i<2*n;i++)
        {
            p2[i] = p1[i];
        }
    }

    /* Elements 0   .. 2*n-1 in work are used for ffac values,
     * Elements 2*n .. 4*n-1 are internal FFTPACK work space.
     */

    if (dir == FFTPACK_FORWARD)
    {
        fftpack_cfftf1(n,(double *)out_data,fft->work+2*n,fft->work,fft->ifac, -1);
    }
    else if (dir == FFTPACK_BACKWARD)
    {
        fftpack_cfftf1(n,(double *)out_data,fft->work+2*n,fft->work,fft->ifac, 1);
    }
    else
    {
        fprintf(stderr,"FFT plan mismatch - bad plan or direction.");
        return EINVAL;
    }

    return 0;
}





int
fftpack_exec_2d          (fftpack_t                  fft,
                          enum fftpack_direction     dir,
                          t_complex *                in_data,
                          t_complex *                out_data)
{
    int                i,nx,ny;
    t_complex *    data;

    nx = fft->n;
    ny = fft->next->n;

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 2D there is likely enough data to benefit from memcpy().
     */
    if (in_data != out_data)
    {
        memcpy(out_data,in_data,sizeof(t_complex)*nx*ny);
    }

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* y transforms */
    for (i=0;i<nx;i++)
    {
        fftpack_exec_1d(fft->next,dir,data+i*ny,data+i*ny);
    }

    /* Transpose in-place to get data in place for x transform now */
    fftpack_transpose_2d(data,data,nx,ny);

    /* x transforms */
    for (i=0;i<ny;i++)
    {
        fftpack_exec_1d(fft,dir,data+i*nx,data+i*nx);
    }

    /* Transpose in-place to get data back in original order */
    fftpack_transpose_2d(data,data,ny,nx);

    return 0;
}





int
fftpack_exec_3d     (fftpack_t                  fft,
                     enum fftpack_direction     dir,
                     t_complex *                in_data,
                     t_complex *                out_data)
{
    int              i,nx,ny,nz,rc;
    t_complex *      data;

    nx=fft->n;
    ny=fft->next->n;
    nz=fft->next->next->n;

    /* FFTPACK only does in-place transforms, so emulate out-of-place
     * by copying data to the output array first.
     * For 3D there is likely enough data to benefit from memcpy().
     */
    if (in_data != out_data)
    {
        memcpy(out_data,in_data,sizeof(t_complex)*nx*ny*nz);
    }

    /* Much easier to do pointer arithmetic when base has the correct type */
    data = (t_complex *)out_data;

    /* Perform z transforms */
    for (i=0;i<nx*ny;i++)
        fftpack_exec_1d(fft->next->next,dir,data+i*nz,data+i*nz);

    /* For each X slice, transpose the y & z dimensions inside the slice */
    for (i=0;i<nx;i++)
    {
        fftpack_transpose_2d(data+i*ny*nz,data+i*ny*nz,ny,nz);
    }

    /* Array is now (nx,nz,ny) - perform y transforms */
    for (i=0;i<nx*nz;i++)
    {
        fftpack_exec_1d(fft->next,dir,data+i*ny,data+i*ny);
    }

    /* Transpose back to (nx,ny,nz) */
    for (i=0;i<nx;i++)
    {
        fftpack_transpose_2d(data+i*ny*nz,data+i*ny*nz,nz,ny);
    }

    /* Transpose entire x & y slices to go from
     * (nx,ny,nz) to (ny,nx,nz).
     */
    rc=fftpack_transpose_2d_nelem(data,data,nx,ny,nz);
    if (rc != 0)
    {
        fprintf(stderr,"Fatal error - cannot transpose X & Y/Z in fftpack_exec_3d().");
        return rc;
    }

    /* Then go from (ny,nx,nz) to (ny,nz,nx) */
    for (i=0;i<ny;i++)
    {
        fftpack_transpose_2d(data+i*nx*nz,data+i*nx*nz,nx,nz);
    }

    /* Perform x transforms */
    for (i=0;i<ny*nz;i++)
    {
        fftpack_exec_1d(fft,dir,data+i*nx,data+i*nx);
    }

    /* Transpose back from (ny,nz,nx) to (ny,nx,nz) */
    for (i=0;i<ny;i++)
    {
        fftpack_transpose_2d(data+i*nz*nx,data+i*nz*nx,nz,nx);
    }

    /* Transpose from (ny,nx,nz) to (nx,ny,nz).
     */
    rc = fftpack_transpose_2d_nelem(data,data,ny,nx,nz);
    if (rc != 0)
    {
        fprintf(stderr,"Fatal error - cannot transpose Y/Z & X in fftpack_exec_3d().");
        return rc;
    }

    return 0;
}



void
fftpack_destroy(fftpack_t      fft)
{
    if (fft != NULL)
    {
        free(fft->work);
        if (fft->next != NULL)
            fftpack_destroy(fft->next);
        free(fft);
    }
}


