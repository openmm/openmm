/* See LICENSE below for information on rights to use, modify and distribute
   this code. */

/*
 * hilbert.c - Computes Hilbert space-filling curve coordinates, without
 * recursion, from integer index, and vice versa, and other Hilbert-related
 * calculations.  Also known as Pi-order or Peano scan.
 *
 * Author:      Doug Moore
 *              Dept. of Computational and Applied Math
 *              Rice University
 *              http://www.caam.rice.edu/~dougm
 * Date:        Sun Feb 20 2000
 * Copyright (c) 1998-2000, Rice University
 *
 * Acknowledgement:
 * This implementation is based on the work of A. R. Butz ("Alternative
 * Algorithm for Hilbert's Space-Filling Curve", IEEE Trans. Comp., April,
 * 1971, pp 424-426) and its interpretation by Spencer W. Thomas, University
 * of Michigan (http://www-personal.umich.edu/~spencer/Home.html) in his widely
 * available C software.  While the implementation here differs considerably
 * from his, the first two interfaces and the style of some comments are very
 * much derived from his work. */


#include "hilbert.h"

/* implementation of the hilbert functions */

#define adjust_rotation(rotation,nDims,bits)                            \
do {                                                                    \
      /* rotation = (rotation + 1 + ffs(bits)) % nDims; */              \
      bits &= -bits & nd1Ones;                                          \
      while (bits)                                                      \
        bits >>= 1, ++rotation;                                         \
      if ( ++rotation >= nDims )                                        \
        rotation -= nDims;                                              \
} while (0)

#define ones(T,k) ((((T)2) << (k-1)) - 1)

#define rdbit(w,k) (((w) >> (k)) & 1)

#define rotateRight(arg, nRots, nDims)                                  \
((((arg) >> (nRots)) | ((arg) << ((nDims)-(nRots)))) & ones(bitmask_t,nDims))

#define rotateLeft(arg, nRots, nDims)                                   \
((((arg) << (nRots)) | ((arg) >> ((nDims)-(nRots)))) & ones(bitmask_t,nDims))

#define DLOGB_BIT_TRANSPOSE
static bitmask_t
bitTranspose(unsigned nDims, unsigned nBits, bitmask_t inCoords)
#if defined(DLOGB_BIT_TRANSPOSE)
{
  unsigned const nDims1 = nDims-1;
  unsigned inB = nBits;
  unsigned utB;
  bitmask_t inFieldEnds = 1;
  bitmask_t inMask = ones(bitmask_t,inB);
  bitmask_t coords = 0;

  while ((utB = inB / 2))
    {
      unsigned const shiftAmt = nDims1 * utB;
      bitmask_t const utFieldEnds =
	inFieldEnds | (inFieldEnds << (shiftAmt+utB));
      bitmask_t const utMask =
	(utFieldEnds << utB) - utFieldEnds;
      bitmask_t utCoords = 0;
      unsigned d;
      if (inB & 1)
	{
	  bitmask_t const inFieldStarts = inFieldEnds << (inB-1);
	  unsigned oddShift = 2*shiftAmt;
	  for (d = 0; d < nDims; ++d)
	    {
	      bitmask_t in = inCoords & inMask;
	      inCoords >>= inB;
	      coords |= (in & inFieldStarts) <<	oddShift++;
	      in &= ~inFieldStarts;
	      in = (in | (in << shiftAmt)) & utMask;
	      utCoords |= in << (d*utB);
	    }
	}
      else
	{
	  for (d = 0; d < nDims; ++d)
	    {
	      bitmask_t in = inCoords & inMask;
	      inCoords >>= inB;
	      in = (in | (in << shiftAmt)) & utMask;
	      utCoords |= in << (d*utB);
	    }
	}
      inCoords = utCoords;
      inB = utB;
      inFieldEnds = utFieldEnds;
      inMask = utMask;
    }
  coords |= inCoords;
  return coords;
}
#else
{
  bitmask_t coords = 0;
  unsigned d;
  for (d = 0; d < nDims; ++d)
    {
      unsigned b;
      bitmask_t in = inCoords & ones(bitmask_t,nBits);
      bitmask_t out = 0;
      inCoords >>= nBits;
      for (b = nBits; b--;)
	{
	  out <<= nDims;
	  out |= rdbit(in, b);
	}
      coords |= out << d;
    }
  return coords;
}
#endif

/*****************************************************************
 * hilbert_i2c
 *
 * Convert an index into a Hilbert curve to a set of coordinates.
 * Inputs:
 *  nDims:      Number of coordinate axes.
 *  nBits:      Number of bits per axis.
 *  index:      The index, contains nDims*nBits bits
 *              (so nDims*nBits must be <= 8*sizeof(bitmask_t)).
 * Outputs:
 *  coord:      The list of nDims coordinates, each with nBits bits.
 * Assumptions:
 *      nDims*nBits <= (sizeof index) * (bits_per_byte)
 */
void WINDOWS_EXPORT
hilbert_i2c(unsigned nDims, unsigned nBits, bitmask_t index, bitmask_t coord[])
{
  if (nDims > 1)
    {
      bitmask_t coords;
      halfmask_t const nbOnes = ones(halfmask_t,nBits);
      unsigned d;

      if (nBits > 1)
	{
	  unsigned const nDimsBits = nDims*nBits;
	  halfmask_t const ndOnes = ones(halfmask_t,nDims);
	  halfmask_t const nd1Ones= ndOnes >> 1; /* for adjust_rotation */
	  unsigned b = nDimsBits;
	  unsigned rotation = 0;
	  halfmask_t flipBit = 0;
	  bitmask_t const nthbits = ones(bitmask_t,nDimsBits) / ndOnes;
	  index ^= (index ^ nthbits) >> 1;
	  coords = 0;
	  do
	    {
	      halfmask_t bits = (index >> (b-=nDims)) & ndOnes;
	      coords <<= nDims;
	      coords |= rotateLeft(bits, rotation, nDims) ^ flipBit;
	      flipBit = (halfmask_t)1 << rotation;
	      adjust_rotation(rotation,nDims,bits);
	    } while (b);
	  for (b = nDims; b < nDimsBits; b *= 2)
	    coords ^= coords >> b;
	  coords = bitTranspose(nBits, nDims, coords);
	}
      else
	coords = index ^ (index >> 1);

      for (d = 0; d < nDims; ++d)
	{
	  coord[d] = coords & nbOnes;
	  coords >>= nBits;
	}
    }
  else
    coord[0] = index;
}

/*****************************************************************
 * hilbert_c2i
 *
 * Convert coordinates of a point on a Hilbert curve to its index.
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBits:      Number of bits/coordinate.
 *  coord:      Array of n nBits-bit coordinates.
 * Outputs:
 *  index:      Output index value.  nDims*nBits bits.
 * Assumptions:
 *      nDims*nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */
bitmask_t WINDOWS_EXPORT
hilbert_c2i(unsigned nDims, unsigned nBits, bitmask_t const coord[])
{
  if (nDims > 1)
    {
      unsigned const nDimsBits = nDims*nBits;
      bitmask_t index;
      unsigned d;
      bitmask_t coords = 0;
      for (d = nDims; d--; )
	{
	  coords <<= nBits;
	  coords |= coord[d];
	}

      if (nBits > 1)
	{
	  halfmask_t const ndOnes = ones(halfmask_t,nDims);
	  halfmask_t const nd1Ones= ndOnes >> 1; /* for adjust_rotation */
	  unsigned b = nDimsBits;
	  unsigned rotation = 0;
	  halfmask_t flipBit = 0;
	  bitmask_t const nthbits = ones(bitmask_t,nDimsBits) / ndOnes;
	  coords = bitTranspose(nDims, nBits, coords);
	  coords ^= coords >> nDims;
	  index = 0;
	  do
	    {
	      halfmask_t bits = (coords >> (b-=nDims)) & ndOnes;
	      bits = rotateRight(flipBit ^ bits, rotation, nDims);
	      index <<= nDims;
	      index |= bits;
	      flipBit = (halfmask_t)1 << rotation;
	      adjust_rotation(rotation,nDims,bits);
	    } while (b);
	  index ^= nthbits >> 1;
	}
      else
	index = coords;
      for (d = 1; d < nDimsBits; d *= 2)
	index ^= index >> d;
      return index;
    }
  else
    return coord[0];
}

/*****************************************************************
 * Readers and writers of bits
 */

typedef bitmask_t (*BitReader) (unsigned nDims, unsigned nBytes,
				char const* c, unsigned y);
typedef void (*BitWriter) (unsigned d, unsigned nBytes,
	   char* c, unsigned y, int fold);


#if defined(sparc)
#define __BIG_ENDIAN__
#endif

#if defined(__BIG_ENDIAN__)
#define whichByte(nBytes,y) (nBytes-1-y/8)
#define setBytes(dst,pos,nBytes,val) \
     memset(&dst[pos+1],val,nBytes-pos-1)
#else
#define whichByte(nBytes,y) (y/8)
#define setBytes(dst,pos,nBytes,val) \
     memset(&dst[0],val,pos)
#endif

static bitmask_t
getIntBits(unsigned nDims, unsigned nBytes, char const* c, unsigned y)
{
  unsigned const bit = y%8;
  unsigned const offs = whichByte(nBytes,y);
  unsigned d;
  bitmask_t bits = 0;
  c += offs;
  for (d = 0; d < nDims; ++d)
    {
      bits |= rdbit(*c, bit) << d;
      c += nBytes;
    }
  return bits;
}

#include <string.h>
static void
propogateIntBits(unsigned d, unsigned nBytes,
		 char* c, unsigned y, int fold)
{
  unsigned const byteId = whichByte(nBytes,y);
  unsigned const b = y%8;
  char const bthbit = 1 << b;
  char* const target = &c[d*nBytes];
  target[byteId] ^= bthbit;
  if (!fold)
    {
      char notbit = ((target[byteId] >> b) & 1) - 1;
      if (notbit)
	target[byteId] |= bthbit-1;
      else
	target[byteId] &=  -bthbit;
      setBytes(target,byteId,nBytes,notbit);
    }
}

/* An IEEE double is treated as a 2100 bit number.  In particular, 0 is treated
   as a 1 followed by 2099 zeroes, and negative 0 as a 0 followed by 2099 ones.
   Only 53 bits differ between a number and a zero of the same sign, with the
   position of the 53 determined by the exponent, and the values of the 53 by
   the significand (with implicit leading 1 bit).  Although IEEE 754 uses the
   maximum exponent for NaN's and infinities, this implementation ignores that
   decision, so that infinities and NaN's are treated as very large numbers.
   Note that we do not explicitly construct a 2100 bit bitmask in the IEEE
   routines below. */

enum { IEEEexpBits = 11 };
enum { IEEEsigBits = 52 };
enum { IEEErepBits = (1 << IEEEexpBits) + IEEEsigBits };

typedef union ieee754_double
  {
    double d;

    /* This is the IEEE 754 double-precision format.  */
    struct
      {
#if defined(__BIG_ENDIAN__)
	unsigned int negative:1;
	unsigned int exponent:11;
	/* Together these comprise the mantissa.  */
	unsigned int mantissa0:20;
	unsigned int mantissa1:32;
#else				/* Big endian.  */
	/* Together these comprise the mantissa.  */
	unsigned int mantissa1:32;
	unsigned int mantissa0:20;
	unsigned int exponent:11;
	unsigned int negative:1;
#endif				/* Little endian.  */
      } ieee;
  } ieee754_double;

static bitmask_t
getIEEESignBits(unsigned nDims, double const* c)
{
  unsigned d;
  ieee754_double x;
  bitmask_t bits = 0;
  for (d = 0; d < nDims; ++d)
    {
      x.d = c[d];
      bits |= x.ieee.negative << d;
    }
  return bits;
}

static bitmask_t
getIEEEBits(unsigned nDims,
	    unsigned ignoreMe, /* ignored */
	    char const* cP,
	    unsigned y)
     /* retrieve bits y of elements of double array c, where an expanded IEEE
	double has 2100 bits. */
{
  unsigned d;
  double const* c = (double const*) cP;
  ieee754_double x;
  bitmask_t bits = 0;
  for (x.d = c[d=0]; d < nDims; x.d = c[++d])
    {
      bitmask_t bit = x.ieee.negative;
      unsigned normalized = (x.ieee.exponent != 0);
      unsigned diff = y - (x.ieee.exponent - normalized);
      if (diff <= 52)
	bit ^= 1 & ((diff <  32)? x.ieee.mantissa1 >> diff:
		    (diff <  52)? x.ieee.mantissa0 >> (diff - 32):
		    /* else */    normalized);
      else
	bit ^= (y == IEEErepBits-1);

      bits |= bit << d;
    }
  return bits;
}

static void
propogateIEEEBits(unsigned d, unsigned nBytes,
		  char* cP, unsigned y, int fold)
{
  ieee754_double* x = d + (ieee754_double*) cP;
  unsigned normalized = (x->ieee.exponent != 0);
  unsigned diff = y - (x->ieee.exponent - normalized);
  if (diff < 32)
    {
      unsigned b = 1 << diff;
      unsigned bit = x->ieee.mantissa1 & b;
      x->ieee.mantissa1 &= ~(b-1);
      x->ieee.mantissa1 |= b;
      if (bit)
	--x->ieee.mantissa1;
    }
  else if (diff < 52)
    {
      unsigned b = 1 << (diff - 32);
      unsigned bit = x->ieee.mantissa0 & b;
      x->ieee.mantissa0 &= ~(b-1);
      x->ieee.mantissa0 |= b;
      if (bit)
	--x->ieee.mantissa0;
      x->ieee.mantissa1 = bit?-1: 0;
    }
  else if (diff == 52) /* "flip" the implicit 1 bit */
    {
      if (normalized)
	--x->ieee.exponent;
      else
	x->ieee.exponent = 1;
      x->ieee.mantissa0 = -normalized;
      x->ieee.mantissa1 = -normalized;
    }
  else if (diff < IEEErepBits)
    {
      if (y == IEEErepBits-1)
	{
	  x->ieee.negative ^= 1;
	  x->ieee.exponent = 0;
	}
      else
	x->ieee.exponent = y - 51;
      x->ieee.mantissa0 = 0;
      x->ieee.mantissa1 = 0;
    }
}

static unsigned
getIEEEexptMax(unsigned nDims, double const* c)
{
  unsigned max = 0;
  unsigned d;
  for (d = 0; d < nDims; ++d)
    {
      ieee754_double x;
      x.d = c[d];
      if (max < x.ieee.exponent)
	max = x.ieee.exponent;
    }
  if (max) --max;
  return max;
}

static void
getIEEEinitValues(double const* c1,
		  unsigned y,
		  unsigned nDims,
		  unsigned* rotation,
		  bitmask_t* bits,
		  bitmask_t* index)
{
  bitmask_t const one = 1;
  unsigned d;
  bitmask_t signBits = getIEEESignBits(nDims, c1);
  unsigned signParity, leastZeroBit, strayBit;

  /* compute the odd/evenness of the number of sign bits */
  {
    bitmask_t signPar = signBits;
    for (d = 1; d < nDims; d *= 2)
      signPar ^= signPar >> d;
    signParity = signPar & 1;
  }

  /* find the position of the least-order 0 bit in among signBits and adjust it
     if necessary */
  for (leastZeroBit = 0; leastZeroBit < nDims; ++leastZeroBit)
    if (rdbit(signBits, leastZeroBit) == 0)
      break;
  strayBit = 0;
  if (leastZeroBit == nDims-2)
    strayBit = 1;
  else if (leastZeroBit == nDims)
    leastZeroBit = nDims-1;

  if (y % 2 == 1)
    {
      *rotation = (IEEErepBits - y + 1 + leastZeroBit) % nDims;
      if (y < IEEErepBits-1)
	{
	  *bits = signBits ^ (one << ((*rotation + strayBit) % nDims));
	  *index = signParity;
	}
      else /* y == IEEErepBits-1 */
	{
	  *bits = signBits ^ (ones(bitmask_t,nDims) &~ 1);
	  *index =  signParity ^ (nDims&1);
	}
    }
  else /* y % 2 == 0 */
    if (y < IEEErepBits)
      {
	unsigned shift_amt = (IEEErepBits - y + leastZeroBit) % nDims;
	*rotation = (shift_amt + 2 + strayBit) % nDims;
	*bits = signBits ^ (one << shift_amt);
	*index = signParity ^ 1;
      }
    else /* y == IEEErepBits */
      {
	*rotation = 0;
	*bits = one << (nDims-1);
	*index = 1;
      }
}

/*****************************************************************
 * hilbert_cmp, hilbert_ieee_cmp
 *
 * Determine which of two points lies further along the Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes of storage/coordinate (hilbert_cmp only)
 *  nBits:      Number of bits/coordinate. (hilbert_cmp only)
 *  coord1:     Array of nDims nBytes-byte coordinates (or doubles for ieee_cmp).
 *  coord2:     Array of nDims nBytes-byte coordinates (or doubles for ieee_cmp).
 * Return value:
 *      -1, 0, or 1 according to whether
           coord1<coord2, coord1==coord2, coord1>coord2
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */

static int
hilbert_cmp_work(unsigned nDims, unsigned nBytes, unsigned nBits,
		 unsigned max, unsigned y,
		 char const* c1, char const* c2,
		 unsigned rotation,
		 bitmask_t bits,
		 bitmask_t index,
		 BitReader getBits)
{
  bitmask_t const one = 1;
  bitmask_t const nd1Ones = ones(bitmask_t,nDims) >> 1; /* used in adjust_rotation macro */
  while (y-- > max)
    {
      bitmask_t reflection = getBits(nDims, nBytes, c1, y);
      bitmask_t diff = reflection ^ getBits(nDims, nBytes, c2, y);
      bits ^= reflection;
      bits = rotateRight(bits, rotation, nDims);
      if (diff)
	{
	  unsigned d;
	  diff = rotateRight(diff, rotation, nDims);
	  for (d = 1; d < nDims; d *= 2)
	    {
	      index ^= index >> d;
	      bits  ^= bits  >> d;
	      diff  ^= diff  >> d;
	    }
	  return (((index ^ y ^ nBits) & 1) == (bits < (bits^diff)))? -1: 1;
	}
      index ^= bits;
      reflection ^= one << rotation;
      adjust_rotation(rotation,nDims,bits);
      bits = reflection;
    }
  return 0;
}

int WINDOWS_EXPORT
hilbert_cmp(unsigned nDims, unsigned nBytes, unsigned nBits,
	    void const* c1, void const* c2)
{
  bitmask_t const one = 1;
  bitmask_t bits = one << (nDims-1);
  return hilbert_cmp_work(nDims, nBytes, nBits, 0, nBits,
			  (char const*)c1, (char const*)c2,
			  0, bits, bits, getIntBits);
}

int WINDOWS_EXPORT
hilbert_ieee_cmp(unsigned nDims, double const* c1, double const* c2)
{
  unsigned rotation, max;
  bitmask_t bits, index;
  if (getIEEESignBits(nDims, c1) != getIEEESignBits(nDims, c2))
    max = 2047;
  else
    {
      unsigned max1 = getIEEEexptMax(nDims, c1);
      unsigned max2 = getIEEEexptMax(nDims, c2);
      max = (max1 > max2)? max1: max2;
    }

  getIEEEinitValues(c1, max+53, nDims, &rotation, &bits, &index);
  return hilbert_cmp_work(nDims, 8, 64, max, max+53,
			  (char const*)c1, (char const*)c2,
			  rotation, bits, index, getIEEEBits);
}

/*****************************************************************
 * hilbert_box_vtx
 *
 * Determine the first or last vertex of a box to lie on a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes/coordinate.
 *  nBits:      Number of bits/coordinate.
 *  findMin:    Is it the least vertex sought?
 *  coord1:     Array of nDims nBytes-byte coordinates - one corner of box
 *  coord2:     Array of nDims nBytes-byte coordinates - opposite corner
 * Output:
 *      c1 and c2 modified to refer to selected corner
 *      value returned is log2 of size of largest power-of-two-aligned box that
 *      contains the selected corner and no other corners
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */


static unsigned
hilbert_box_vtx_work(unsigned nDims, unsigned nBytes, unsigned nBits,
		     int findMin,
		     unsigned max, unsigned y,
		     char* c1, char* c2,
		     unsigned rotation,
		     bitmask_t bits,
		     bitmask_t index,
		     BitReader getBits)
{
  bitmask_t const one = 1;
  bitmask_t const ndOnes = ones(bitmask_t,nDims);
  bitmask_t const nd1Ones= ndOnes >> 1;
  bitmask_t bitsFolded = 0;

  while (y--)
    {
      bitmask_t reflection = getBits(nDims, nBytes, c1, y);
      bitmask_t diff = reflection ^ getBits(nDims, nBytes, c2, y);
      if (diff)
	{
	  unsigned d;
	  bitmask_t smear = rotateRight(diff, rotation, nDims) >> 1;
	  bitmask_t digit = rotateRight(bits ^ reflection, rotation, nDims);
	  for (d = 1; d < nDims; d *= 2)
	    {
	      index ^= index >> d;
	      digit ^= (digit  >> d) &~ smear;
	      smear |= smear >> d;
	    }
	  index &= 1;
	  if ((index ^ y ^ findMin) & 1)
	    digit ^= smear+1;
	  digit = rotateLeft(digit, rotation, nDims) & diff;
	  reflection ^= digit;

	  for (d = 0; d < nDims; ++d)
	    if (rdbit(diff, d))
	      {
		int way = rdbit(digit, d);
		char* target = d*nBytes + (way? c1: c2);
		char* const source = 2*d*nBytes + c1 - target + c2;
		memcpy(target, source, nBytes);
	      }

	  bitsFolded |= diff;
	  if (bitsFolded == ndOnes)
	    return y;
	}

      bits ^= reflection;
      bits = rotateRight(bits, rotation, nDims);
      index ^= bits;
      reflection ^= one << rotation;
      adjust_rotation(rotation,nDims,bits);
      bits = reflection;
    }
  return y;
}

unsigned WINDOWS_EXPORT
hilbert_box_vtx(unsigned nDims, unsigned nBytes, unsigned nBits,
		int findMin, void* c1, void* c2)
{
  bitmask_t const one = 1;
  bitmask_t bits = one << (nDims-1);
  return hilbert_box_vtx_work(nDims, nBytes, nBits, findMin,
			      0, nBits, (char*)c1, (char*)c2,
			      0, bits, bits, getIntBits);
}

unsigned WINDOWS_EXPORT
hilbert_ieee_box_vtx(unsigned nDims,
		     int findMin, double* c1, double* c2)
{
  unsigned rotation, max;
  bitmask_t bits, index;
  if (getIEEESignBits(nDims, c1) != getIEEESignBits(nDims, c2))
    max = 2047;
  else
    {
      unsigned max1 = getIEEEexptMax(nDims, c1);
      unsigned max2 = getIEEEexptMax(nDims, c2);
      max = (max1 > max2)? max1: max2;
    }

  getIEEEinitValues(c1, max+53, nDims, &rotation, &bits, &index);

  return hilbert_box_vtx_work(nDims, 8, 64, findMin,
			      max, max+53, (char *)c1, (char *)c2,
			      rotation, bits, index, getIEEEBits);
}

/*****************************************************************
 * hilbert_box_pt
 *
 * Determine the first or last point of a box to lie on a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes/coordinate.
 *  nBits:      Number of bits/coordinate.
 *  findMin:    Is it the least vertex sought?
 *  coord1:     Array of nDims nBytes-byte coordinates - one corner of box
 *  coord2:     Array of nDims nBytes-byte coordinates - opposite corner
 * Output:
 *      c1 and c2 modified to refer to least point
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */
unsigned WINDOWS_EXPORT
hilbert_box_pt_work(unsigned nDims, unsigned nBytes, unsigned nBits,
		    int findMin,
		    unsigned max, unsigned y,
		    char* c1, char* c2,
		    unsigned rotation,
		    bitmask_t bits,
		    bitmask_t index,
		    BitReader getBits,
		    BitWriter propogateBits)
{
  bitmask_t const one = 1;
  bitmask_t const nd1Ones = ones(bitmask_t,nDims) >> 1;
  bitmask_t fold1 = 0, fold2 = 0;
  unsigned smearSum = 0;

  while (y-- > max)
    {
      bitmask_t reflection = getBits(nDims, nBytes, c1, y);
      bitmask_t diff = reflection ^ getBits(nDims, nBytes, c2, y);
      if (diff)
	{
	  bitmask_t smear = rotateRight(diff, rotation, nDims) >> 1;
	  bitmask_t digit = rotateRight(bits ^ reflection, rotation, nDims);
	  unsigned d;
	  for (d = 1; d < nDims; d *= 2)
	    {
	      index ^= index >> d;
	      digit ^= (digit  >> d) &~ smear;
	      smear |= smear >> d;
	    }
	  smearSum += smear;
	  index &= 1;
	  if ((index ^ y ^ findMin) & 1)
	    digit ^= smear+1;
	  digit = rotateLeft(digit, rotation, nDims) & diff;
	  reflection ^= digit;

	  for (d = 0; d < nDims; ++d)
	    if (rdbit(diff, d))
	      {
		int way = rdbit(digit, d);
		char* c = way? c1: c2;
		bitmask_t fold = way? fold1: fold2;
		propogateBits(d, nBytes, c, y, rdbit(fold, d));
	      }
	  diff ^= digit;
	  fold1 |= digit;
	  fold2 |= diff;
	}

      bits ^= reflection;
      bits = rotateRight(bits, rotation, nDims);
      index ^= bits;
      reflection ^= one << rotation;
      adjust_rotation(rotation,nDims,bits);
      bits = reflection;
    }
  return smearSum;
}

unsigned WINDOWS_EXPORT
hilbert_box_pt(unsigned nDims, unsigned nBytes, unsigned nBits,
		int findMin, void* c1, void* c2)
{
  bitmask_t const one = 1;
  bitmask_t bits = one << (nDims-1);
  return hilbert_box_pt_work(nDims, nBytes, nBits, findMin,
			     0, nBits, (char*)c1, (char*)c2,
			     0, bits, bits,
			     getIntBits, propogateIntBits);
}

unsigned WINDOWS_EXPORT
hilbert_ieee_box_pt(unsigned nDims,
		    int findMin, double* c1, double* c2)
{
  unsigned rotation, max;
  bitmask_t bits, index;
  bitmask_t c1Signs = getIEEESignBits(nDims, c1);
  bitmask_t c2Signs = getIEEESignBits(nDims, c2);
  if (c1Signs != c2Signs)
    {
      rotation = 0;
      bits = (bitmask_t)1 << (nDims-1);
      index = 1;
      hilbert_box_pt_work(nDims, 8, 64, findMin,
			  IEEErepBits-1, IEEErepBits, (char *)c1, (char *)c2,
			  rotation, bits, index,
			  getIEEEBits, propogateIEEEBits);
    }

  /* having put everything in the same orthant, start */
  {
    unsigned max1 = getIEEEexptMax(nDims, c1);
    unsigned max2 = getIEEEexptMax(nDims, c2);
    max = (max1 > max2)? max1: max2;
  }

  getIEEEinitValues(c1, max+53, nDims, &rotation, &bits, &index);

  return hilbert_box_pt_work(nDims, 8, 64, findMin,
			     max, max+53, (char *)c1, (char *)c2,
			     rotation, bits, index,
			     getIEEEBits, propogateIEEEBits);
}

/*****************************************************************
 * hilbert_nextinbox
 *
 * Determine the first point of a box after or before a given point to lie on
 * a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes/coordinate.
 *  nBits:      Number of bits/coordinate.
 *  findPrev:   Is it a previous point that you want?
 *  coord1:     Array of nDims nBytes-byte coordinates - one corner of box
 *  coord2:     Array of nDims nBytes-byte coordinates - opposite corner
 *  point:      Array of nDims nBytes-byte coordinates - lower bound on point returned
 *
 * Output:
      if returns 1:
 *      c1 and c2 modified to refer to least point after "point" in box
      else returns 0:
        arguments unchanged; "point" is beyond the last point of the box
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */
int WINDOWS_EXPORT
hilbert_nextinbox(unsigned nDims, unsigned nBytes, unsigned nBits,
		  int findPrev, void* c1V, void* c2V, void const* ptV)
{
  bitmask_t const one = 1;
  unsigned y = nBits;
  bitmask_t const ndOnes = ones(bitmask_t,nDims);
  bitmask_t const nd1Ones = ndOnes >> 1;
  unsigned rotation = 0;
  bitmask_t bits = 0;
  bitmask_t index = 0;
  bitmask_t fold1 = 0, fold2 = 0;
  bitmask_t valu1 = 0, valu2 = 0;
  unsigned p_y;
  bitmask_t p_separator = 0, p_firstSeparator;
  bitmask_t p_cornerdiff, p_reflection;
  bitmask_t p_fold1, p_fold2, p_valu1, p_valu2;

  char* c1 = (char*)c1V;
  char* c2 = (char*)c2V;
  char const* pt = (const char*)ptV;

  while (y-- > 0)
    {
      bitmask_t reflection = getIntBits(nDims, nBytes, pt, y);
      bitmask_t diff = reflection ^ /* planes that separate box and point */
	((getIntBits(nDims, nBytes, c1, y) &~ fold1) | valu1);

      if (diff)
	/* some coordinate planes separate point from box or
	   dividing box or both; smear the bits of diff to reflect that
	   after the first diff dimension, they might as well all be
	   diffing; adjust the diff to reflect the fact that diffed
	   dimensions don't matter. */
	{
	  /* compute (the complement of) a "digit" in the integer index of this
	     point */
	  bitmask_t cornerdiff = (diff ^ reflection) ^ /* separate box crnrs */
	    ((getIntBits(nDims, nBytes, c2, y) &~ fold2) | valu2);
	  bitmask_t separator = diff & ~cornerdiff;
	  /* eventually, the most significant separating cutting plane */
	  bitmask_t firstSeparator;
	  /* bits less significant than the msb of separator are irrelevant;
	     for convenience, call them all separators too */
	  bitmask_t rotSep = rotateRight(separator, rotation, nDims);
	  /* compute the (complement of the) digit of the hilbert code
	     assoc with point */
	  bitmask_t digit = rotateRight(bits ^ reflection, rotation, nDims);
	  unsigned d;
	  for (d = 1; d < nDims; d *= 2)
	    {
	      index ^= index >> d;
	      digit ^= digit >> d;
	      rotSep |= rotSep >> d;
	    }
	  index &= 1;
	  digit &= rotSep;
	  if ((index ^ y ^ findPrev) & 1)
	    digit ^= rotSep;

	  separator = rotateLeft(rotSep, rotation, nDims);
	  rotSep -= rotSep >> 1;
	  firstSeparator = rotateLeft(rotSep, rotation, nDims);
	  /* forget about all the planes that split the box, except those that
	     are more significant than the most significant separator. */
	  cornerdiff &= ~separator;

	  if (cornerdiff && digit)
	    /* some coordinate planes divide the box.  Call the part of the
	       box in the same orthant as the point "here" and the part of
	       the box in the next (or previous) orthant "there".  Remember
	       what the "there" orthant of the box looks like in case it
	       turns out that the curve doesn't reenter the box "here" after
	       (before) passing thru point.  Continue working with the
	       "here" part. If there is no "there" there, skip it */
	    {
	      p_firstSeparator = digit & -digit;
	      p_separator = 2*p_firstSeparator-1;
	      p_separator = rotateLeft(p_separator, rotation, nDims);
	      p_firstSeparator = rotateLeft(p_firstSeparator, rotation, nDims);
	      p_cornerdiff = cornerdiff &~ (p_separator ^ p_firstSeparator);
	      p_y = y;
	      p_reflection = reflection ^ p_firstSeparator;
	      p_fold1 = fold1;
	      p_fold2 = fold2;
	      p_valu1 = valu1;
	      p_valu2 = valu2;
	    }

	  if (digit < rotSep)

	    /* use next box */
	    {
	      if (!p_separator) return 0; /* no next point */
	      separator = p_separator;
	      firstSeparator = p_firstSeparator;
	      y = p_y;
	      cornerdiff = p_cornerdiff;
	      reflection = p_reflection;
	      fold1 = p_fold1;
	      fold2 = p_fold2;
	      valu1 = p_valu1;
	      valu2 = p_valu2;
	    }

	  if (cornerdiff)
	    {
	      /* reduce currbox */
	      bitmask_t corner = diff & cornerdiff;
	      cornerdiff ^= corner;
	      fold1 |= corner;
	      fold2 |= cornerdiff;
	      valu1 |= ~reflection & corner;
	      valu2 |= ~reflection & cornerdiff;
	    }

	  separator ^= firstSeparator;
	  if (firstSeparator)
	    /* we have completely separated the point from a part of the box
	       ahead of it on the curve; almost done */
	    {
	      unsigned byteId = whichByte(nBytes,y);
	      bitmask_t bthbit = one << y%8;
	      for (d = 0; d < nDims; ++d)
		{
		  char lo1, lo2;
		  char* cc1 = &c1[d*nBytes];
		  char* cc2 = &c2[d*nBytes];
		  char const* pnt = &pt[d*nBytes];
		  char hibits = -bthbit;
		  char hipart = pnt[byteId] & hibits;
		  memcpy(cc1, pnt, byteId);
		  memcpy(cc2, pnt, byteId);

		  if (rdbit(separator, d))
		    hibits ^= bthbit;
		  if (rdbit(firstSeparator, d))
		    hipart ^= bthbit;

		  if (rdbit(fold1, d))
		    {
		      lo1 = -rdbit(valu1, d);
		      setBytes(cc1,byteId,nBytes,lo1);
		    }
		  else lo1 = cc1[byteId];
		  cc1[byteId] = hipart | (lo1 &~ hibits);

		  if (rdbit(fold2, d))
		    {
		      lo2 = -rdbit(valu2, d);
		      setBytes(cc2,byteId,nBytes,lo2);
		    }
		  else lo2 = cc2[byteId];
		  cc2[byteId] = hipart | (lo2 &~ hibits);
		}

	      hilbert_box_pt(nDims, nBytes, nBits, !findPrev, c1V, c2V);
	      return 1;
	    }
	}

      bits ^= reflection;
      bits = rotateRight(bits, rotation, nDims);
      index ^= bits;
      reflection ^= one << rotation;
      adjust_rotation(rotation,nDims,bits);
      bits = reflection;
    }

  /* point is in box */
  {
    unsigned d;
    for (d = 0; d < nDims; ++d)
      ((char*)c1)[d] = ((char*)c2)[d] = ((char*)pt)[d];
  }
  return 1;
}



/*****************************************************************
 * hilbert_incr
 *
 * Advance from one point to its successor on a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBits:      Number of bits/coordinate.
 *  coord:      Array of nDims nBits-bit coordinates.
 * Output:
 *  coord:      Next point on Hilbert curve
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */

void WINDOWS_EXPORT
hilbert_incr(unsigned nDims, unsigned nBits, bitmask_t coord[])
{
  bitmask_t const one = 1;
  bitmask_t const ndOnes = ones(bitmask_t,nDims);
  bitmask_t const nd1Ones= ndOnes >> 1;
  unsigned b, d;
  unsigned rotation = 0;
  bitmask_t reflection = 0;
  bitmask_t index = 0;
  unsigned rb = nBits-1;
  bitmask_t rd = ndOnes;

  for (b = nBits; b--;)
    {
      bitmask_t bits = reflection;
      reflection = 0;
      for (d = 0; d < nDims; ++d)
        reflection |= rdbit(coord[d], b) << d;
      bits ^= reflection;
      bits = rotateRight(bits, rotation, nDims);
      index ^= bits;
      for (d = 1; d < nDims; d *= 2)
        index ^= index >> d;
      if (index++ != ndOnes)
        {
          rb = b;
          rd = index & -index;
          rd = rotateLeft(rd, rotation, nDims);

        }
      index &= 1;
      index <<= nDims-1;

      reflection ^= one << rotation;
      adjust_rotation(rotation,nDims,bits);
    }
  for (d = 0; !rdbit(rd, d); ++d) {}
  coord[d] ^= (2 << rb) - 1;
}


/* LICENSE
 *
 * This software is copyrighted by Rice University.  It may be freely copied,
 * modified, and redistributed, provided that the copyright notice is
 * preserved on all copies.
 *
 * There is no warranty or other guarantee of fitness for this software,
 * it is provided solely "as is".  Bug reports or fixes may be sent
 * to the author, who may or may not act on them as he desires.
 *
 * You may include this software in a program or other software product,
 * but must display the notice:
 *
 * Hilbert Curve implementation copyright 1998, Rice University
 *
 * in any place where the end-user would see your own copyright.
 *
 * If you modify this software, you should include a notice giving the
 * name of the person performing the modification, the date of modification,
 * and the reason for such modification.
 */



/* Revision history:

   July 1998: Initial release

   Sept 1998: Second release

   Dec 1998: Fixed bug in hilbert_c2i that allowed a shift by number of bits in
   bitmask to vaporize index, in last bit of the function.  Implemented
   hilbert_incr.

   August 1999: Added argument to hilbert_nextinbox so that you can, optionally,
   find the previous point along the curve to intersect the box, rather than the
   next point.

   Nov 1999: Defined fast bit-transpose function (fast, at least, if the number
   of bits is large), and reimplemented i2c and c2i in terms of it.  Collapsed
   loops in hilbert_cmp, with the intention of reusing the cmp code to compare
   more general bitstreams.

   Feb 2000: Implemented almost all the floating point versions of cmp, etc, so
   that coordinates expressed in terms of double-precision IEEE floating point
   can be ordered.  Still have to do next-in-box, though.

   Oct 2001: Learned that some arbitrary coding choices caused some routines
   to fail in one dimension, and changed those choices.

   version 2001-10-20-05:34

*/

/* What remains is test code that won't be compiled unless you define the
   TEST_HILBERT preprocessor symbol */

#ifdef TEST_HILBERT
#include <stdio.h>
#define abs(x) (((x)>=0)?(x):(-(x)))

int main()
{
#define maxDim (8*sizeof(bitmask_t))
  bitmask_t coord[maxDim], coordPrev[maxDim];
  unsigned nDims, nBits, nPrints, orderCheck, i;
  bitmask_t r, r1;

  for (;;)
    {
      printf( "Enter nDims, nBits, nPrints, orderCheck: " );
      scanf( "%d", &nDims);
      if ( nDims == 0 )
        break;
      scanf( "%d%d%d", &nBits, &nPrints, &orderCheck);
      while ( (i = getchar()) != '\n' && i != EOF )
        ;
      if ( i == EOF )
        break;

      if (nDims*nBits > 8*sizeof(r))
        {
          printf("Product of nDims and nBits not exceed %d.\n", 8*sizeof(r));
          break;
        }

      if (nBits == 0)
        {
          printf("nBits must be positive.\n");
          break;
        }

      if (nPrints > (1ULL << (nDims*nBits)))
        nPrints = 1ULL << (nDims*nBits);

      for (r = 0; r < nPrints; ++r)
        {
	  bitmask_t coord1[maxDim];
	  int miscount = 0;
          hilbert_i2c( nDims, nBits, r, coord );
          printf("%d: ", (unsigned)r);
          for (i = 0; i < nDims; ++i)
	    {
	      int diff = (int)(coord[i] - coordPrev[i]);
	      miscount += abs(diff);
	      coordPrev[i] = coord[i];
	      printf(" %d", (unsigned)coord[i]);
	    }
	  if (r > 0 && miscount != 1)
	    printf(".....error");
	  printf("\n");
          r1 = hilbert_c2i( nDims, nBits, coord );
          if ( r != r1 )
            printf( "r = 0x%x; r1 = 0x%x\n", (unsigned)r, (unsigned)r1);
	  for (i = 0; i < nDims; ++i)
	    coord[i] = coordPrev[i];

	  if (! orderCheck)
	    continue;

	  for (r1 = 0; r1 < r; ++r1 )
	    {
	      unsigned ans;
	      hilbert_i2c( nDims, nBits, r1, coord1 );
	      ans = hilbert_cmp( nDims, sizeof(coord[0]), nBits, coord, coord1);
	      if (ans != 1)
		{
		  int width = (nDims*nBits + 3) / 4;
		  printf( "cmp r = 0x%0*x; r1 = 0x%0*x, ans = %2d\n",
			  width, (unsigned)r,
			  width, (unsigned)r1, ans );
		}
	    }
	  hilbert_i2c( nDims, nBits, r1, coord1 );
	  if (hilbert_cmp( nDims, sizeof(coord[0]), nBits, coord, coord1) != 0)
	    printf( "cmp r = 0x%0*x; r1 = 0x%0*x\n", (nDims*nBits+3)/4, (unsigned)r,
		    (nDims*nBits+3)/4, (unsigned)r1 );

        }
    }
  return 0;
}

#endif

#ifdef TEST_IEEE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int cmp(const void* xv, const void* yv)
{
  double const* x = xv;
  double const* y = yv;
  /* return hilbert_cmp(2, 8, 64, x, y); */
  return hilbert_ieee_cmp(2, x, y);
}

int main()
{
  double *a;
  unsigned i;
  unsigned n;
  printf("How many points? ");
  scanf("%d", &n);
  a = (double*) malloc(2*n*sizeof(double));
  for (i = 0; i < n; ++i)
    a[2*i] = drand48()-0.5, a[2*i+1] = drand48()-0.5;

  qsort(a, n, 2*sizeof(double), cmp);

  for (i = 0; i < n; ++i)
    printf("%8g %8g\n", a[2*i], a[2*i+1]);
  free(a);
  return 0;
}

#endif

#ifdef TEST_CMP
#include <stdio.h>

#define maxDim (8*sizeof(bitmask_t))
int main()
{
  double coord[maxDim];
  unsigned nDims, i, k;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;
  while ( (i = getchar()) != '\n' && i != EOF )
    ;
  if ( i == EOF )
    return 0;

  for (k = 0; k < (1<<nDims); ++k)
    {
      printf("Orth %2d\n", k);
      for (i = 0; i < nDims; ++i)
	coord[i] = ((k>>i)&1)? -1.: 1.;


      hilbert_ieee_cmp( nDims, coord, coord);
    }
  return 0;
}

#endif

#ifdef TEST_VTX
#include <stdio.h>
#include <stdlib.h>

#define maxDim (8*sizeof(bitmask_t))

unsigned g_nDims;

int cmp(void const* c1p, void const* c2p)
{
  return hilbert_cmp(g_nDims, sizeof(unsigned), 8*sizeof(unsigned), c1p, c2p);
}

int main()
{
  unsigned corner0[maxDim], corner1[maxDim];
  unsigned cornerlo[maxDim], cornerhi[maxDim], work[maxDim];
  typedef unsigned array_t[maxDim];
  array_t* array;

  unsigned nDims, i, k;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;
  while ( (i = getchar()) != '\n' && i != EOF )
    ;
  if ( i == EOF )
    return 0;

  printf("Enter one corner (%d coordinates): ", nDims);
  for (k = 0; k < nDims; ++k)
    scanf("%d", &corner0[k]);

  printf("Enter other corner (%d coordinates): ", nDims);
  for (k = 0; k < nDims; ++k)
    scanf("%d", &corner1[k]);


  /* find first corner */
  for (k = 0; k < nDims; ++k)
    {
      cornerlo[k] = corner0[k];
      work[k] = corner1[k];
    }

  hilbert_box_vtx(nDims, sizeof(unsigned), 8*sizeof(unsigned),
		  1, cornerlo, work);
  printf("Predicted lo corner: ");
  for (k = 0; k < nDims; ++k)
    printf("%4u", cornerlo[k]);
  printf("\n");


  /* find last corner */
  for (k = 0; k < nDims; ++k)
    {
      work[k] = corner0[k];
      cornerhi[k] = corner1[k];
    }

  hilbert_box_vtx(nDims, sizeof(unsigned), 8*sizeof(unsigned),
		  0, work, cornerhi);
  printf("Predicted hi corner: ");
  for (k = 0; k < nDims; ++k)
    printf("%4u", cornerhi[k]);
  printf("\n");

  array = (array_t*) malloc(maxDim*sizeof(unsigned) << nDims);
  for (k = 0; k < (1<<nDims); ++k)
    {
      unsigned j;
      unsigned* eltk = &array[k][0];
      for (j = 0; j < nDims; ++j)
	{
	  unsigned* src = ((k>>j)&1)? corner1: corner0;
	  eltk[j] = src[j];
	}
    }

  g_nDims = nDims;
  qsort(array, (1<<nDims), maxDim*sizeof(unsigned), cmp);

  printf("Result of sort\n");
  for (k = 0; k < (1<<nDims); k += (1 << nDims) - 1)
    {
      unsigned j;
      unsigned* eltk = &array[k][0];
      for (j = 0; j < nDims; ++j)
	printf("%4u", eltk[j]);
      printf("\n");
    }
  free((char*)array);
  return 0;
}

#endif

#ifdef TEST_IEEE_VTX
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define maxDim (8*sizeof(bitmask_t))
typedef double key_t;

unsigned g_nDims;

int cmp(void const* c1p, void const* c2p)
{
  return hilbert_ieee_cmp(g_nDims, c1p, c2p);
}

int main()
{
  key_t corner0[maxDim], corner1[maxDim];
  key_t cornerlo[maxDim], cornerhi[maxDim], work[maxDim];
  typedef key_t array_t[maxDim];
  array_t* array;

  unsigned nDims, i, k;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;

  for (i = 0; i < 10000; ++i)
    {
      for (k = 0; k < nDims; ++k)
	{
	  corner0[k] = 2.*drand48() - 1.;
	  corner1[k] = 2.*drand48() - 1.;
	}

      /* find first corner */
      for (k = 0; k < nDims; ++k)
	{
	  cornerlo[k] = corner0[k];
	  work[k] = corner1[k];
	}

      hilbert_ieee_box_vtx(nDims, 1, cornerlo, work);

      /* find last corner */
      for (k = 0; k < nDims; ++k)
	{
	  work[k] = corner0[k];
	  cornerhi[k] = corner1[k];
	}

      hilbert_ieee_box_vtx(nDims, 0, work, cornerhi);

      array = (array_t*) malloc(maxDim*sizeof(key_t) << nDims);
      for (k = 0; k < (1<<nDims); ++k)
	{
	  unsigned j;
	  key_t* eltk = &array[k][0];
	  for (j = 0; j < nDims; ++j)
	    {
	      key_t* src = ((k>>j)&1)? corner1: corner0;
	      eltk[j] = src[j];
	    }
	}

      g_nDims = nDims;
      qsort(array, (1<<nDims), maxDim*sizeof(key_t), cmp);

      for (k = 0; k < (1<<nDims); k += (1 << nDims) - 1)
	{
	  unsigned j;
	  int mismatch = 0;
	  key_t* eltk = &array[k][0];
	  for (j = 0; j < nDims & !mismatch; ++j)
	    {
	      mismatch = (eltk[j] != ((k==0)? cornerlo: cornerhi)[j]);
	    }
	  assert (!mismatch);
	}
      free((char*)array);
    }
  return 0;
}

#endif

#ifdef TEST_PT
#include <stdio.h>
#include <stdlib.h>

#define maxDim (8*sizeof(bitmask_t))

unsigned g_nDims;

int cmp(void const* c1p, void const* c2p)
{
  return hilbert_cmp(g_nDims, sizeof(unsigned), 8*sizeof(unsigned), c1p, c2p);
}

int main()
{
  unsigned point0[maxDim], point1[maxDim];
  unsigned pointlo[maxDim], pointhi[maxDim], work[maxDim];
  typedef unsigned array_t[maxDim];
  array_t* array;

  unsigned nDims, i, k, outvolume = 1, involume = 1;
  unsigned nextItem;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;
  while ( (i = getchar()) != '\n' && i != EOF )
    ;
  if ( i == EOF )
    return 0;

  printf("Enter one point (%d coordinates): ", nDims);
  for (k = 0; k < nDims; ++k)
    scanf("%d", &point0[k]);

  printf("Enter other point (%d coordinates, strictly greater): ", nDims);
  for (k = 0; k < nDims; ++k)
    {
      unsigned diff;
      scanf("%d", &point1[k]);
      diff = point1[k] - point0[k];
      outvolume *= diff + 1;
      involume *= diff - 1;
    }


  /* find first point */
  for (k = 0; k < nDims; ++k)
    {
      pointlo[k] = point0[k];
      work[k] = point1[k];
    }

  hilbert_box_pt(nDims, sizeof(unsigned), 8*sizeof(unsigned),
		 1, pointlo, work);
  printf("Predicted lo point: ");
  for (k = 0; k < nDims; ++k)
    printf("%4u", pointlo[k]);
  printf("\n");


  /* find last point */
  for (k = 0; k < nDims; ++k)
    {
      work[k] = point0[k];
      pointhi[k] = point1[k];
    }

  hilbert_box_pt(nDims, sizeof(unsigned), 8*sizeof(unsigned),
		 0, work, pointhi);
  printf("Predicted hi point: ");
  for (k = 0; k < nDims; ++k)
    printf("%4u", pointhi[k]);
  printf("\n");



  array = (array_t*) malloc(maxDim*sizeof(unsigned) * (outvolume-involume));
  if (array == 0)
    {
      fprintf(stderr, "Out of memory.\n");
      exit(-1);
    }
  nextItem = 0;
  for (k = 0; k < outvolume; ++k)
    {
      unsigned kk = k;
      unsigned j;
      unsigned* eltk = &array[nextItem][0];
      int boundary = 0;

      for (j = 0; j < nDims; ++j)
	{
	  unsigned diff1 = point1[j] - point0[j] + 1;
	  unsigned pos = point0[j] + (kk % diff1);
	  boundary |= (point0[j] == pos || pos == point1[j]);
	  eltk[j] = pos;
	  kk /= diff1;
	}
      if (boundary)
	++nextItem;
    }

  g_nDims = nDims;
  qsort(array, outvolume-involume, maxDim*sizeof(unsigned), cmp);

  printf("Result of sort\n");
  for (k = 0; k < outvolume-involume; k += outvolume-involume-1)
    {
      unsigned j;
      unsigned* eltk = &array[k][0];
      for (j = 0; j < nDims; ++j)
	printf("%4u", eltk[j]);
      printf("\n");
    }
  free((char*)array);
  return 0;
}

#endif

#ifdef TEST_IEEE_PT
#include <stdio.h>
#include <stdlib.h>

#define maxDim (8*sizeof(bitmask_t))

int main()
{
  double point0[maxDim], point1[maxDim];
  double pointlo[maxDim], pointhi[maxDim], work[maxDim];

  unsigned nDims, k, i;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;
  while ( (i = getchar()) != '\n' && i != EOF )
    ;
  if ( i == EOF )
    return 0;

  printf("Enter one point (%d coordinates): ", nDims);
  for (k = 0; k < nDims; ++k)
    scanf("%lf", &point0[k]);

  printf("Enter other point (%d coordinates, strictly greater): ", nDims);
  for (k = 0; k < nDims; ++k)
    scanf("%lf", &point1[k]);

  /* find last point */
  for (k = 0; k < nDims; ++k)
    {
      work[k] = point0[k];
      pointhi[k] = point1[k];
    }

  hilbert_ieee_box_pt(nDims, 0, work, pointhi);
  printf("Predicted hi point: ");
  for (k = 0; k < nDims; ++k)
    printf("%10lg", pointhi[k]);
  printf("\n");

  /* find first point */
  for (k = 0; k < nDims; ++k)
    {
      pointlo[k] = point0[k];
      work[k] = point1[k];
    }

  hilbert_ieee_box_pt(nDims, 1, pointlo, work);
  printf("Predicted lo point: ");
  for (k = 0; k < nDims; ++k)
    printf("%10lg", pointlo[k]);
  printf("\n");

  /* validate by sorting random boundary points */
#define nPts 1000000
  assert(hilbert_ieee_cmp(nDims, pointlo, pointhi) < 0);
  for (i = 0; i < nPts; ++i)
    {
      double pt1[maxDim], pt2[maxDim];
      for (k = 0; k < nDims; ++k)
	{
	  if (i % nDims == k)
	    pt1[k] = point0[k];
	  else
	    pt1[k] = point0[k] + drand48()*(point1[k]-point0[k]);
	}
      for (k = 0; k < nDims; ++k)
	{
	  if (i % nDims == k)
	    pt2[k] = point1[k];
	  else
	    pt2[k] = point0[k] + drand48()*(point1[k]-point0[k]);
	}
      if (hilbert_ieee_cmp(nDims, pt1, pt2) < 0)
	{
	  if (hilbert_ieee_cmp(nDims, pt1, pointlo) < 0)
	    memcpy(pointlo, pt1, maxDim*sizeof(double));
	  if (hilbert_ieee_cmp(nDims, pointhi, pt2) < 0)
	    memcpy(pointhi, pt2, maxDim*sizeof(double));
	}
      else
	{
	  if (hilbert_ieee_cmp(nDims, pt2, pointlo) < 0)
	    memcpy(pointlo, pt2, maxDim*sizeof(double));
	  if (hilbert_ieee_cmp(nDims, pointhi, pt1) < 0)
	    memcpy(pointhi, pt1, maxDim*sizeof(double));
	}
    }

  printf("Sorted hi and lo:\n");
  for (k = 0; k < nDims; ++k)
    printf("%10lg", pointhi[k]);
  printf("\n");
  for (k = 0; k < nDims; ++k)
    printf("%10lg", pointlo[k]);
  printf("\n");

  return 0;
}

#endif

#ifdef TEST_NEXT
#include <stdio.h>

int main()
{
  unsigned i;
  unsigned c1[100], c2[100], pt[100];
  unsigned nDims, nBytes = 4;
  int stat, findPrev;
  printf("Enter nDims: " );
  scanf("%u", &nDims);

  printf("Enter 1st box corner: ");
  for (i = 0; i < nDims; ++i)
    scanf("%u", &c1[i]);
  printf("Enter 2nd box corner: ");
  for (i = 0; i < nDims; ++i)
    scanf("%u", &c2[i]);
  printf("Enter point: ");
  for (i = 0; i < nDims; ++i)
    scanf("%u", &pt[i]);
  printf("Find prev?: ");
  scanf("%d", &findPrev);

  stat =  hilbert_nextinbox(nDims, nBytes, 8*nBytes, findPrev, c1, c2, pt);

  if (stat)
    for (i = 0; i < nDims; ++i)
      printf("%u ", c1[i]);
  else
    printf("No such point");

  printf("\n");
  return 0;
}
#endif

