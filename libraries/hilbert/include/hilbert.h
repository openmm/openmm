/* C header file for Hilbert curve functions */
#if !defined(_hilbert_h_)
#define _hilbert_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _MSC_VER
/* define the bitmask_t type as an integer of sufficient size */
typedef unsigned long long bitmask_t;
/* define the halfmask_t type as an integer of 1/2 the size of bitmask_t */
typedef unsigned int halfmask_t;
#if defined(OPENMM_BUILDING_SHARED_LIBRARY)
    #define WINDOWS_EXPORT __declspec(dllexport)
#else
#define WINDOWS_EXPORT
#endif
#else
#include <stdint.h>
/* define the bitmask_t type as an integer of sufficient size */
typedef uint64_t bitmask_t;
/* define the halfmask_t type as an integer of 1/2 the size of bitmask_t */
typedef uint32_t halfmask_t;
#define WINDOWS_EXPORT
#endif

/*****************************************************************
 * hilbert_i2c
 *
 * Convert an index into a Hilbert curve to a set of coordinates.
 * Inputs:
 *  nDims:      Number of coordinate axes.
 *  nBits:      Number of bits per axis.
 *  index:      The index, contains nDims*nBits bits (so nDims*nBits must be <= 8*sizeof(bitmask_t)).
 * Outputs:
 *  coord:      The list of nDims coordinates, each with nBits bits.
 * Assumptions:
 *      nDims*nBits <= (sizeof index) * (bits_per_byte)
 */

void WINDOWS_EXPORT hilbert_i2c(unsigned nDims, unsigned nBits, bitmask_t index, bitmask_t coord[]);

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

bitmask_t WINDOWS_EXPORT hilbert_c2i(unsigned nDims, unsigned nBits, bitmask_t const coord[]);

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

int WINDOWS_EXPORT hilbert_cmp(unsigned nDims, unsigned nBytes, unsigned nBits, void const* coord1, void const* coord2);
int WINDOWS_EXPORT hilbert_ieee_cmp(unsigned nDims, double const* coord1, double const* coord2);

/*****************************************************************
 * hilbert_box_vtx
 *
 * Determine the first or last vertex of a box to lie on a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes/coordinate.
 *  nBits:      Number of bits/coordinate. (hilbert_cmp only)
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
unsigned WINDOWS_EXPORT
hilbert_box_vtx(unsigned nDims, unsigned nBytes, unsigned nBits,
		int findMin, void* c1, void* c2);
unsigned WINDOWS_EXPORT
hilbert_ieee_box_vtx(unsigned nDims,
		     int findMin, double* c1, double* c2);

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
hilbert_box_pt(unsigned nDims, unsigned nBytes, unsigned nBits,
	       int findMin, void* coord1, void* coord2);
unsigned WINDOWS_EXPORT
hilbert_ieee_box_pt(unsigned nDims,
		    int findMin, double* c1, double* c2);

/*****************************************************************
 * hilbert_nextinbox
 *
 * Determine the first point of a box after a given point to lie on a Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes/coordinate.
 *  nBits:      Number of bits/coordinate.
 *  findPrev:   Is the previous point sought?
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
		  int findPrev, void* coord1, void* coord2,
		  void const* point);

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
hilbert_incr(unsigned nDims, unsigned nBits, bitmask_t coord[]);

#ifdef __cplusplus
}
#endif

#endif /* _hilbert_h_ */
