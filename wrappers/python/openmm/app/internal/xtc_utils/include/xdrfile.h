/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile.h,v 1.6 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */
#ifndef __XXX
#define __XXX 1

//#include "rpc/types.h"
//#include "rpc/xdr.h"


#include <inttypes.h>
enum xdr_op
{
  XDR_ENCODE = 0,
  XDR_DECODE = 1,
  XDR_FREE   = 2
};

typedef struct XDR XDR;

struct XDR
{
  enum xdr_op x_op;
  struct xdr_ops
  {
    int (*x_getlong) (XDR *__xdrs, int32_t *__lp);
    int (*x_putlong) (XDR *__xdrs, int32_t *__lp);
    int (*x_getbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
    int (*x_putbytes) (XDR *__xdrs, char *__addr, unsigned int __len);
    /* two next routines are not 64-bit IO safe - don't use! */
    unsigned int (*x_getpostn) (XDR *__xdrs);
    int (*x_setpostn) (XDR *__xdrs, unsigned int __pos);
    void (*x_destroy) (XDR *__xdrs);
  }
    *x_ops;
  char *x_private;
};

struct XDRFILE
{
    FILE *   fp;       /**< pointer to standard C library file handle */
    XDR *    xdr;      /**< pointer to corresponding XDR handle       */
    char     mode;     /**< r=read, w=write, a=append                 */
    int *    buf1;     /**< Buffer for internal use                   */
    int      buf1size; /**< Current allocated length of buf1          */
    int *    buf2;     /**< Buffer for internal use                   */
    int      buf2size; /**< Current allocated length of buf2          */
};


/*! \file  xdrfile.h
 *  \brief Interface to read/write portabile binary files using XDR.
 *
 * This file provides an interface to read & write portably binary files,
 * using XDR - the external data representation standard defined in RFC 1014. 
 * 
 * There are several advantages to the XDR approach:
 *
 * -# It is portable. And not just portable between big/small integer endian,
 *    but truly portable if you have system XDR routines. For example:
 *       - It doesn't matter if the character representation is ASCII or EBCDIC.
 *       - Some systems are small endian but use big endian order of the two
 *         dword in a double precision floating-point variable. The system XDR
 *         libraries will read/write this correctly.
 *       - Some systems (VAX...) don't use IEEE floating point. Their system
 *         XDR libraries will convert to/from this automatically.
 * -# XDR libraries are required for NFS and lots of other network functions.
 *    This means there isn't a single Unix-like system that doesn't have them.
 * -# There is NO extra metadata whatsoever, and we write plain XDR files.
 *    If you write a float, it will take exactly 4 bytes in the file. 
 *    (All basic datatypes are 4 bytes, double fp 8 bytes). 
 * -# You can read/write the files by calling the system XDR routines directly
 *    too - you don't have to use the routines defined in this file.
 * -# It is no problem if your system doesn't have XDR libraries (MS Windows).
 *    We have written our own versions of the necessary routines that work if
 *    your system uses ASCII for strings and IEEE floating-point. All types
 *    of byte and dword endian for integer and floating-point are supported.
 * -# You can use these routines for any type of data, but since we designed
 *    them for Gromacs we also provide a special routine to write coordinates
 *    with (adjustable) lossy compression. The default precision will give you
 *    three decimals guaranteed accuracy, and reduces the filesize to 1/10th
 *    of normal binary data.
 *
 * We do not support getting or setting positions in XDR files, since it can
 * break in horrible ways for large (64-bit) files, resulting in silent data
 * corruption. Note that it works great to open/read/write 64-bit files if
 * your system supports it; it is just the random access we cannot trust!

#include "rpc/types.h"
#include "rpc/xdr.h"
 *
 * We also provide wrapper routines so this module can be used from FORTRAN -
 * see the file xdrfile_fortran.txt in the Gromacs distribution for 
 * documentation on the FORTRAN interface!
 */


#ifndef _XDRFILE_H_
#define _XDRFILE_H_

#ifdef __cplusplus
extern "C" 
{
#endif

	/*! \brief Abstract datatype for an portable binary file handle 
	 *
	 *  This datatype essentially works just like the standard FILE type in C.
	 *  The actual contents is hidden in the implementation, so you can only
	 *  define pointers to it, for use with the xdrfile routines. 
	 *  
	 *  If you \a really need to see the definition it is in xdrfile.c, but you
	 *  cannot access elements of the structure outside that file.
	 *
	 *  \warning The implementation is completely different from the C standard
	 *  library FILE, so don't even think about using an XDRFILE pointer as an 
	 *  argument to a routine that needs a standard FILE pointer.
	 */
	typedef struct XDRFILE XDRFILE;

	enum { exdrOK, exdrHEADER, exdrSTRING, exdrDOUBLE, 
		   exdrINT, exdrFLOAT, exdrUINT, exdr3DX, exdrCLOSE, exdrMAGIC,
		   exdrNOMEM, exdrENDOFFILE, exdrFILENOTFOUND, exdrDUFF, exdrNR };

	extern char *exdr_message[exdrNR];

#define DIM 3
	typedef float matrix[DIM][DIM];
	typedef float rvec[DIM];
	typedef int   mybool;


	/*! \brief Open a portable binary file, just like fopen()
	 *
	 *  Use this routine much like calls to the standard library function
	 *  fopen(). The only difference is that the returned pointer should only
	 *  be used with routines defined in this header.
	 *
	 *  \param path  Full or relative path (including name) of the file
	 *  \param mode  "r" for reading, "w" for writing, "a" for append.
	 *
	 *  \return Pointer to abstract xdr file datatype, or NULL if an error occurs.
	 *
	 */
	XDRFILE *
	xdrfile_open    (const char *    path, 
					 const char *    mode);


	/*! \brief Close a previously opened portable binary file, just like fclose()
	 *
	 *  Use this routine much like calls to the standard library function
	 *  fopen(). The only difference is that it is used for an XDRFILE handle
	 *  instead of a FILE handle.
	 * 
	 *  \param xfp  Pointer to an abstract XDRFILE datatype
	 *
	 *  \return     0 on success, non-zero on error. 
	 */
	int
	xdrfile_close   (XDRFILE *       xfp);




	/*! \brief Read one or more \a char type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of characters to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of characters read
	 */
	int
	xdrfile_read_char(char *      ptr, 
					  int         ndata, 
					  XDRFILE *   xfp);



	/*! \brief Write one or more \a characters type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of characters to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of characters written
	 */
	int
	xdrfile_write_char(char *      ptr, 
					   int         ndata, 
					   XDRFILE *   xfp);



	/*! \brief Read one or more \a unsigned \a char type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of unsigned characters to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned characters read
	 */
	int
	xdrfile_read_uchar(unsigned char *    ptr, 
					   int		          ndata, 
					   XDRFILE *          xfp);



	/*! \brief Write one or more \a unsigned \a characters type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of unsigned characters to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned characters written
	 */
	int
	xdrfile_write_uchar(unsigned char *   ptr, 
						int               ndata, 
						XDRFILE *         xfp);



	/*! \brief Read one or more \a short type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of shorts to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of shorts read
	 */
	int
	xdrfile_read_short(short *             ptr, 
					   int                 ndata, 
					   XDRFILE *           xfp);



	/*! \brief Write one or more \a short type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of shorts to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of shorts written
	 */
	int
	xdrfile_write_short(short *            ptr, 
						int                ndata, 
						XDRFILE *          xfp);



	/*! \brief Read one or more \a unsigned \a short type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of unsigned shorts to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned shorts read
	 */
	int
	xdrfile_read_ushort(unsigned short *   ptr, 
						int                ndata, 
						XDRFILE *          xfp);



	/*! \brief Write one or more \a unsigned \a short type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of unsigned shorts to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned shorts written
	 */
	int
	xdrfile_write_ushort(unsigned short *     ptr, 
						 int                  ndata, 
						 XDRFILE *            xfp);


	/*! \brief Read one or more \a integer type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of integers to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of integers read
	 *
	 *  The integer data type is assumed to be less than or equal to 32 bits.
	 *
	 *  We do not provide any routines for reading/writing 64-bit integers, since
	 *  - Not all XDR implementations support it
	 *  - Not all machines have 64-bit integers
	 *
	 *  Split your 64-bit data into two 32-bit integers for portability!
	 */
	int
	xdrfile_read_int(int *         ptr, 
					 int           ndata, 
					 XDRFILE *     xfp);



	/*! \brief Write one or more \a integer type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of integers to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of integers written
	 *
	 *  The integer data type is assumed to be less than or equal to 32 bits.
	 *
	 *  We do not provide any routines for reading/writing 64-bit integers, since
	 *  - Not all XDR implementations support it
	 *  - Not all machines have 64-bit integers 
	 *
	 *  Split your 64-bit data into two 32-bit integers for portability!
	 */
	int
	xdrfile_write_int(int *        ptr, 
					  int          ndata, 
					  XDRFILE *    xfp);

	/*! \brief Read one or more \a unsigned \a integers type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of unsigned integers to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned integers read
	 *
	 *  The integer data type is assumed to be less than or equal to 32 bits.
	 *
	 *  We do not provide any routines for reading/writing 64-bit integers, since
	 *  - Not all XDR implementations support it
	 *  - Not all machines have 64-bit integers
	 *
	 *  Split your 64-bit data into two 32-bit integers for portability!
	 */
	int
	xdrfile_read_uint(unsigned int *    ptr, 
					  int               ndata, 
					  XDRFILE *         xfp);



	/*! \brief Write one or more \a unsigned \a integer type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of unsigned integers to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of unsigned integers written
	 *
	 *  The integer data type is assumed to be less than or equal to 32 bits.
	 *
	 *  We do not provide any routines for reading/writing 64-bit integers, since
	 *  - Not all XDR implementations support it
	 *  - Not all machines have 64-bit integers 
	 *
	 *  Split your 64-bit data into two 32-bit integers for portability!
	 */
	int
	xdrfile_write_uint(unsigned int *    ptr, 
					   int               ndata,  
					   XDRFILE *         xfp);



	/*! \brief Read one or more \a float type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of floats to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of floats read
	 */
	int
	xdrfile_read_float(float *           ptr, 
					   int               ndata, 
					   XDRFILE *         xfp);



	/*! \brief Write one or more \a float type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of floats to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of floats written
	 */
	int
	xdrfile_write_float(float *          ptr, 
						int              ndata, 
						XDRFILE *        xfp);



	/*! \brief Read one or more \a double type variable(s) 
	 *
	 *  \param ptr    Pointer to memory where data should be written
	 *  \param ndata  Number of doubles to read
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of doubles read
	 */
	int
	xdrfile_read_double(double *          ptr, 
						int               ndata, 
						XDRFILE *         xfp);



	/*! \brief Write one or more \a double type variable(s)
	 *
	 *  \param ptr    Pointer to memory where data should be read 
	 *  \param ndata  Number of double to write.
	 *  \param xfp    Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return       Number of doubles written
	 */
	int
	xdrfile_write_double(double *        ptr, 
						 int             ndata, 
						 XDRFILE *       xfp);



	/*! \brief Read a string (array of characters)
	 *
	 *  \param ptr     Pointer to memory where data should be written
	 *  \param maxlen  Maximum length of string. If no end-of-string is encountered,
	 *                 one byte less than this is read and end-of-string appended.
	 *  \param xfp     Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return        Number of characters read, including end-of-string
	 */
	int
	xdrfile_read_string(char *          ptr, 
						int             maxlen, 
						XDRFILE *       xfp);



	/*! \brief Write a string (array of characters)
	 *
	 *  \param ptr     Pointer to memory where data should be read
	 *  \param xfp     Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return        Number of characters written, including end-of-string
	 */
	int
	xdrfile_write_string(char *          ptr, 
						 XDRFILE *       xfp);



	/*! \brief Read raw bytes from file (unknown datatype)
	 *
	 *  \param ptr     Pointer to memory where data should be written
	 *  \param nbytes  Number of bytes to read. No conversion whatsoever is done.
	 *  \param xfp     Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return        Number of bytes read from file
	 */
	int
	xdrfile_read_opaque(char *             ptr, 
						int                nbytes, 
						XDRFILE *          xfp);




	/*! \brief Write raw bytes to file (unknown datatype)
	 *
	 *  \param ptr     Pointer to memory where data should be read
	 *  \param nbytes  Number of bytes to write. No conversion whatsoever is done.
	 *  \param xfp     Handle to portable binary file, created with xdrfile_open()
	 *
	 *  \return        Number of bytes written to file
	 */
	int
	xdrfile_write_opaque(char *            ptr, 
						 int               nbytes, 
						 XDRFILE *         xfp);






	/*! \brief Compress coordiates in a float array to XDR file
	 *
	 *  This routine will perform \a lossy compression on the three-dimensional
	 *  coordinate data data specified and store it in the XDR file.
	 *
	 *  The lossy part of the compression consists of multiplying each
	 *  coordinate with the precision argument and then rounding to integers.
	 *  We suggest a default value of 1000.0, which means you are guaranteed
	 *  three decimals of accuracy. The only limitation is that scaled coordinates
	 *  must still fit in an integer variable, so if the precision is 1000.0 the
	 *  coordinate magnitudes must be less than +-2e6.
	 *
	 *  \param ptr        Pointer to coordinates to compress (length 3*ncoord)
	 *  \param ncoord     Number of coordinate triplets in data
	 *  \param precision  Scaling factor for lossy compression. If it is <=0, 
	 *                    the default value of 1000.0 is used.
	 *  \param xfp        Handle to portably binary file
	 *
	 *  \return           Number of coordinate triplets written. 
	 *                    IMPORTANT: Check that this is equal to ncoord - if it is
	 *                    negative, an error occured. This should not happen with 
	 *	   	              normal data, but if your coordinates are NaN or very
	 *                    large (>1e6) it is not possible to use the compression.
	 *
	 *  \warning          The compression algorithm is not part of the XDR standard,
	 *                    and very complicated, so you will need this xdrfile module 
	 *                    to read it later. 
	 */
	int
	xdrfile_compress_coord_float(float *     ptr,
								 int         ncoord,
								 float       precision,
								 XDRFILE *   xfp);




	/*! \brief Decompress coordiates from XDR file to array of floats
	 *
	 *  This routine will decompress three-dimensional coordinate data previously 
	 *  stored in an XDR file and store it in the specified array of floats.
	 *
	 *  The precision used during the earlier compression is read from the file
	 *  and returned - you cannot adjust the accuracy at this stage.
	 *
	 *  \param ptr        Pointer to coordinates to compress (length>= 3*ncoord)
	 *  \param ncoord     Max number of coordinate triplets to read on input, actual
	 *                    number of coordinate triplets read on return. If this
	 *                    is smaller than the number of coordinates in the frame an
	 *                    error will occur.
	 *  \param precision  The precision used in the previous compression will be
	 *                    written to this variable on return.
	 *  \param xfp        Handle to portably binary file
	 *
	 *  \return           Number of coordinate triplets read. If this is negative,
	 *                    an error occured.
	 *
	 *  \warning          Since we cannot count on being able to set/get the 
	 *                    position of large files (>2Gb), it is not possible to
	 *                    recover from errors by re-reading the frame if the 
	 *                    storage area you provided was too small. To avoid this 
	 *                    from happening, we recommend that you store the number of 
	 *                    coordinates triplet as an integer either in a header or
	 *                    just before the compressed coordinate data, so you can 
	 *                    read it first and allocated enough memory.
	 */
	int
	xdrfile_decompress_coord_float(float *     ptr,
								   int *	   ncoord,
								   float *     precision,
								   XDRFILE *   xfp);




	/*! \brief Compress coordiates in a double array to XDR file
	 *
	 *  This routine will perform \a lossy compression on the three-dimensional
	 *  coordinate data data specified and store it in the XDR file. Double will
	 *  NOT give you any extra precision since the coordinates are compressed. This
	 *  routine just avoids allocating a temporary array of floats.
	 *
	 *  The lossy part of the compression consists of multiplying each
	 *  coordinate with the precision argument and then rounding to integers.
	 *  We suggest a default value of 1000.0, which means you are guaranteed
	 *  three decimals of accuracy. The only limitation is that scaled coordinates
	 *  must still fit in an integer variable, so if the precision is 1000.0 the
	 *  coordinate magnitudes must be less than +-2e6.
	 *
	 *  \param ptr        Pointer to coordinates to compress (length 3*ncoord)
	 *  \param ncoord     Number of coordinate triplets in data
	 *  \param precision  Scaling factor for lossy compression. If it is <=0, the
	 *                    default value of 1000.0 is used.
	 *  \param xfp        Handle to portably binary file
	 *
	 *  \return           Number of coordinate triplets written. 
	 *                    IMPORTANT: Check that this is equal to ncoord - if it is
	 *                    negative, an error occured. This should not happen with 
	 *                    normal data, but if your coordinates are NaN or very
	 *                    large (>1e6) it is not possible to use the compression.
	 *
	 *  \warning          The compression algorithm is not part of the XDR standard,
	 *                    and very complicated, so you will need this xdrfile module 
	 *                    to read it later. 
	 */
	int
	xdrfile_compress_coord_double(double *     ptr,
								  int          ncoord,
								  double       precision,
								  XDRFILE *    xfp);




	/*! \brief Decompress coordiates from XDR file to array of doubles
	 *
	 *  This routine will decompress three-dimensional coordinate data previously 
	 *  stored in an XDR file and store it in the specified array of doubles.
	 *  Double will NOT give you any extra precision since the coordinates are 
	 *  compressed. This routine just avoids allocating a temporary array of floats.
	 *
	 *  The precision used during the earlier compression is read from the file
	 *  and returned - you cannot adjust the accuracy at this stage.
	 *
	 *  \param ptr        Pointer to coordinates to compress (length>= 3*ncoord)
	 *  \param ncoord     Max number of coordinate triplets to read on input, actual
	 *                    number of coordinate triplets read on return. If this
	 *                    is smaller than the number of coordinates in the frame an
	 *                    error will occur.
	 *  \param precision  The precision used in the previous compression will be
	 *                    written to this variable on return.
	 *  \param xfp        Handle to portably binary file
	 *
	 *  \return           Number of coordinate triplets read. If this is negative,
	 *                    an error occured.
	 *
	 *  \warning          Since we cannot count on being able to set/get the 
	 *                    position of large files (>2Gb), it is not possible to
	 *                    recover from errors by re-reading the frame if the 
	 *                    storage area you provided was too small. To avoid this 
	 *                    from happening, we recommend that you store the number of 
	 *                    coordinates triplet as an integer either in a header or
	 *                    just before the compressed coordinate data, so you can 
	 *                    read it first and allocated enough memory.
	 */
	int
	xdrfile_decompress_coord_double(double *     ptr,
									int *	     ncoord,
									double *     precision,
									XDRFILE *    xfp);



#ifdef __cplusplus
}
#endif

#endif /* _XDRFILE_H_ */
#endif

