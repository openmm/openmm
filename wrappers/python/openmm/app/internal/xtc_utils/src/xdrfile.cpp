/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: xdrfile.c,v 1.3 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */


/* Get HAVE_RPC_XDR_H, F77_FUNC from config.h if available */
#include <string.h>
#ifdef PLATFORM_Linux
#if defined(__i386__) || defined(__x86_64__)
__asm__(".symver memcpy,memcpy@GLIBC_2.2.5");
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
static void * condfree( void*p) {
        if(p) { free(p); }
	else {
		fprintf( stderr, "Double free detected\n" );
		abort();
	}
        return NULL;
}

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>


/* get fixed-width types if we are using ANSI C99 */
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#elif (defined HAVE_INTTYPES_H)
#  include <inttypes.h>
#endif


#ifdef HAVE_RPC_XDR_H
#  include <rpc/rpc.h>
#  include <rpc/xdr.h>
#endif

#include "xdrfile.h"

/* Default FORTRAN name mangling is: lower case name, append underscore */
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

char *exdr_message[exdrNR] = {
	"OK", 
	"Header",
	"String", 
	"Double",
	"Integer",
	"Float",
	"Unsigned integer",
	"Compressed 3D coordinate",
	"Closing file",
	"Magic number",
	"Not enough memory",
	"End of file",
	"File not found" 
};

/*
 * Declare our own XDR routines statically if no libraries are present.
 * Actual implementation is at the end of this file.
 *
 * We don't want the low-level XDR implementation as part of the Gromacs
 * documentation, so skip it for doxygen too...
 */
#if (!defined HAVE_RPC_XDR_H && !defined DOXYGEN)


/* We need integer types that are guaranteed to be 4 bytes wide.
 * If ANSI C99 headers were included they are already defined
 * as int32_t and uint32_t. Check, and if not define them ourselves.
 * Since it is just our workaround for missing ANSI C99 types, avoid adding
 * it to the doxygen documentation.
 */
#if !(defined INT32_MAX || defined DOXYGEN)
#    if (INT_MAX == 2147483647)
#        define int32_t int
#        define uint32_t unsigned int
#        define INT32_MAX 2147483647
#    elif (LONG_MAX == 2147483647)
#        define int32_t long
#        define uint32_t unsigned long
#        define INT32_MAX 2147483647L
#    else
#        error ERROR: No 32 bit wide integer type found!
#        error Use system XDR libraries instead, or update xdrfile.c
#    endif
#endif

static int  xdr_char        (XDR *xdrs, char *ip);
static int  xdr_u_char      (XDR *xdrs, unsigned char *ip);
static int  xdr_short       (XDR *xdrs, short *ip);
static int  xdr_u_short     (XDR *xdrs, unsigned short *ip);
static int  xdr_int         (XDR *xdrs, int *ip);
static int  xdr_u_int       (XDR *xdrs, unsigned int *ip);
static int  xdr_float       (XDR *xdrs, float *ip);
static int  xdr_double      (XDR *xdrs, double *ip);
static int  xdr_string      (XDR *xdrs, char **ip, unsigned int maxsize);
static int  xdr_opaque      (XDR *xdrs, char *cp, unsigned int cnt);
static void xdrstdio_create (XDR *xdrs, FILE *fp, enum xdr_op xop);

#define xdr_getpos(xdrs)                                \
        (*(xdrs)->x_ops->x_getpostn)(xdrs)
#define xdr_setpos(xdrs, pos)                           \
        (*(xdrs)->x_ops->x_setpostn)(xdrs, pos)
#define xdr_destroy(xdrs)                                       \
        do {                                                    \
                if ((xdrs)->x_ops->x_destroy)                   \
                        (*(xdrs)->x_ops->x_destroy)(xdrs);      \
        } while (0)
#endif /* end of our own XDR declarations */





/** Contents of the abstract XDRFILE data structure.
 *
 *  @internal
 * 
 *  This structure is used to provide an XDR file interface that is
 *  virtual identical to the standard UNIX fopen/fread/fwrite/fclose.
 */



/*************************************************************
 * Implementation of higher-level routines to read/write     * 
 * portable data based on the XDR standard. These should be  *
 * called from C - see further down for Fortran77 wrappers.  *
 *************************************************************/

XDRFILE *
xdrfile_open(const char *path, const char *mode)
{
	char newmode[5];
	enum xdr_op xdrmode;
	XDRFILE *xfp;
  
	/* make sure XDR files are opened in binary mode... */
	if(*mode=='w' || *mode=='W') 
    {
		sprintf(newmode,"wb+");
		xdrmode=XDR_ENCODE;
	} else if(*mode == 'a' || *mode == 'A') 
    {
		sprintf(newmode,"ab+");
		xdrmode = XDR_ENCODE;
	} else if(*mode == 'r' || *mode == 'R')
    {
		sprintf(newmode,"rb");
		xdrmode = XDR_DECODE;
	} else /* cannot determine mode */
		return NULL;
  
	if((xfp=(XDRFILE *)malloc(sizeof(XDRFILE)))==NULL)
		return NULL;
	if((xfp->fp=fopen(path,newmode))==NULL)
    {
		xfp=(XDRFILE *)condfree(xfp);
		return NULL;
	}
	if((xfp->xdr=(XDR *)malloc(sizeof(XDR)))==NULL) 
    {
		fclose(xfp->fp);
		xfp=(XDRFILE *)condfree(xfp);
		return NULL;
	}
	xfp->mode=*mode;
	xdrstdio_create((XDR *)(xfp->xdr),xfp->fp,xdrmode);
	xfp->buf1 = xfp->buf2 = NULL;
	xfp->buf1size = xfp->buf2size = 0;
	return xfp;
}

int 
xdrfile_close(XDRFILE *xfp)
{
	int ret=exdrCLOSE;
	if(xfp) 
    {
		/* flush and destroy XDR stream */
		if(xfp->xdr)
			xdr_destroy((XDR *)(xfp->xdr));
		xfp->xdr=(XDR *)condfree(xfp->xdr);
		/* close the file */
		ret=fclose(xfp->fp);
		if(xfp->buf1size)
			xfp->buf1=(int *)condfree(xfp->buf1);
		if(xfp->buf2size)
			xfp->buf2=(int *)condfree(xfp->buf2);
		xfp=(XDRFILE *)condfree(xfp);
	}
	return ret; /* return 0 if ok */
}



int 
xdrfile_read_int(int *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_int((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_int(int *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_int((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}


int 
xdrfile_read_uint(unsigned int *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_u_int((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_uint(unsigned int *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_u_int((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_read_char(char *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_char((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_char(char *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_char((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}


int 
xdrfile_read_uchar(unsigned char *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_u_char((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_uchar(unsigned char *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_u_char((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_read_short(short *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_short((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_short(short *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_short((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}


int 
xdrfile_read_ushort(unsigned short *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_u_short((XDR *)(xfp->xdr),ptr+i))
		i++;
  
	return i;
}

int 
xdrfile_write_ushort(unsigned short *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
  
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_u_short((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_read_float(float *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_float((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_write_float(float *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_float((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_read_double(double *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_double((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int 
xdrfile_write_double(double *ptr, int ndata, XDRFILE* xfp) 
{
	int i=0;
	/* read write is encoded in the XDR struct */  
	while(i<ndata && xdr_double((XDR *)(xfp->xdr),ptr+i))
		i++;
	return i;
}

int
xdrfile_read_string(char *ptr, int maxlen, XDRFILE* xfp)
{
	int i;
	if(xdr_string((XDR *)(xfp->xdr),&ptr,maxlen)) {
		i=0;
		while(i<maxlen && ptr[i]!=0)
			i++;
		if(i==maxlen)
			return maxlen;
		else
			return i+1;
	} else
		return 0;
}

int
xdrfile_write_string(char *ptr, XDRFILE* xfp)
{
	int len=strlen(ptr)+1;
  
	if(xdr_string((XDR *)(xfp->xdr),&ptr,len)) 
		return len;
	else
		return 0;
}


int
xdrfile_read_opaque(char *ptr, int cnt, XDRFILE* xfp)
{
	if(xdr_opaque((XDR *)(xfp->xdr),ptr,cnt))
		return cnt;
	else 
		return 0;
}


int
xdrfile_write_opaque(char *ptr, int cnt, XDRFILE* xfp)
{
	if(xdr_opaque((XDR *)(xfp->xdr),ptr,cnt))
		return cnt;
	else
		return 0;
}


/* Internal support routines for reading/writing compressed coordinates 
 * sizeofint - calculate smallest number of bits necessary
 * to represent a certain integer.
 */
static int 
sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32) 
    {
		num_of_bits++;
		num <<= 1;
    }
    return num_of_bits;
}


/*
 * sizeofints - calculate 'bitsize' of compressed ints
 *
 * given a number of small unsigned integers and the maximum value
 * return the number of bits needed to read or write them with the
 * routines encodeints/decodeints. You need this parameter when
 * calling those routines. 
 * (However, in some cases we can just use the variable 'smallidx' 
 * which is the exact number of bits, and them we dont need to call
 * this routine).
 */
static int 
sizeofints(int num_of_ints, unsigned int sizes[]) 
{
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++)
    {	
		tmp = 0;
		for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0) 
        {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) 
    {
		num_of_bits++;
		num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}
    

/*
 * encodebits - encode num into buf using the specified number of bits
 *
 * This routines appends the value of num to the bits already present in
 * the array buf. You need to give it the number of bits to use and you had
 * better make sure that this number of bits is enough to hold the value.
 * Num must also be positive.
 */
static void 
encodebits(int buf[], int num_of_bits, int num) 
{
    
    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;
    
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8)
    {
		lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
		cbuf[cnt++] = lastbyte >> lastbits;
		num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
		lastbyte = (lastbyte << num_of_bits) | num;
		lastbits += num_of_bits;
		if (lastbits >= 8) 
        {
			lastbits -= 8;
			cbuf[cnt++] = lastbyte >> lastbits;
		}
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0) 
    {
		cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

/*
 * encodeints - encode a small set of small integers in compressed format
 *
 * this routine is used internally by xdr3dfcoord, to encode a set of
 * small integers to the buffer for writing to a file.
 * Multiplication with fixed (specified maximum) sizes is used to get
 * to one big, multibyte integer. Allthough the routine could be
 * modified to handle sizes bigger than 16777216, or more than just
 * a few integers, this is not done because the gain in compression
 * isn't worth the effort. Note that overflowing the multiplication
 * or the byte buffer (32 bytes) is unchecked and whould cause bad results.
 * THese things are checked in the calling routines, so make sure not
 * to remove those checks...
 */
 
static void 
encodeints(int buf[], int num_of_ints, int num_of_bits,
		   unsigned int sizes[], unsigned int nums[]) 
{

    int i;
    unsigned int bytes[32], num_of_bytes, bytecnt, tmp;

    tmp = nums[0];
    num_of_bytes = 0;
    do 
    {
		bytes[num_of_bytes++] = tmp & 0xff;
		tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++) 
    {
		if (nums[i] >= sizes[i])
        {
			fprintf(stderr,"major breakdown in encodeints - num %u doesn't "
					"match size %u\n", nums[i], sizes[i]);
			abort();
		}
		/* use one step multiply */    
		tmp = nums[i];
		for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) 
        {
			tmp = bytes[bytecnt] * sizes[i] + tmp;
			bytes[bytecnt] = tmp & 0xff;
			tmp >>= 8;
		}
		while (tmp != 0)
        {
			bytes[bytecnt++] = tmp & 0xff;
			tmp >>= 8;
		}
		num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8) 
    {
		for (i = 0; i < num_of_bytes; i++) 
        {
			encodebits(buf, 8, bytes[i]);
		}
		encodebits(buf, num_of_bits - num_of_bytes * 8, 0);
    } 
    else
    {
		for (i = 0; i < num_of_bytes-1; i++)
        {
			encodebits(buf, 8, bytes[i]);
		}
		encodebits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}


/*
 * decodebits - decode number from buf using specified number of bits
 * 
 * extract the number of bits from the array buf and construct an integer
 * from it. Return that value.
 *
 */

static int 
decodebits(int buf[], int num_of_bits) 
{

    int cnt, num; 
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
    
    num = 0;
    while (num_of_bits >= 8)
    {
		lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
		num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
		num_of_bits -=8;
    }
    if (num_of_bits > 0) 
    {
		if (lastbits < num_of_bits) 
        {
			lastbits += 8;
			lastbyte = (lastbyte << 8) | cbuf[cnt++];
		}
		lastbits -= num_of_bits;
		num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num; 
}

/*
 * decodeints - decode 'small' integers from the buf array
 *
 * this routine is the inverse from encodeints() and decodes the small integers
 * written to buf by calculating the remainder and doing divisions with
 * the given sizes[]. You need to specify the total number of bits to be
 * used from buf in num_of_bits.
 *
 */

void 
decodeints(int buf[], int num_of_ints, int num_of_bits,
		   unsigned int sizes[], int nums[])
{

	int bytes[32];
	int i, j, num_of_bytes, p, num;
  
	bytes[1] = bytes[2] = bytes[3] = 0;
	num_of_bytes = 0;
	while (num_of_bits > 8)
    {
		bytes[num_of_bytes++] = decodebits(buf, 8);
		num_of_bits -= 8;
	}
	if (num_of_bits > 0)
    {
		bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
	}
	for (i = num_of_ints-1; i > 0; i--) 
    {
		num = 0;
		for (j = num_of_bytes-1; j >=0; j--) 
        {
			num = (num << 8) | bytes[j];
			switch( sizes[i] ) {
				case 645: p = num / 645; break;
				case 512: p = num / 512; break;
				case 203: p = num / 203; break;
				default: p = num / sizes[i];
			}
//printf("%d\n", sizes[i]);
			bytes[j] = p;

			num = num - p * sizes[i];
		}
		nums[i] = num;
	}
	nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}
    

static const int magicints[] = 
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
    1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003, 
    16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031, 
    131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561, 
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216 
};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

/* Compressed coordinate routines - modified from the original
 * implementation by Frans v. Hoesel to make them threadsafe.
 */
int
xdrfile_decompress_coord_float(float     *ptr,
							   int       *size,
							   float     *precision,
							   XDRFILE*   xfp)
{
	int minint[3], maxint[3], *lip;
	int smallidx,  minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int k, *buf1, *buf2, lsize, flag;
	int smallnum, smaller, larger, i, is_smaller, run;
	float *lfp, inv_precision;
	int tmp, *thiscoord,  prevcoord[3];
	unsigned int bitsize;
	int errval=1;
	int countmjh=0;
 
	sizeint[0] = sizeint[1] = sizeint[2] = 0;
 
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

	if(xfp==NULL || ptr==NULL)
		return -1;
	tmp=xdrfile_read_int(&lsize,1,xfp);
	if(tmp==0)
		return -1; /* return if we could not read size */
	if (*size < lsize) 
    {
		fprintf(stderr, "Requested to decompress %d coords, file contains %d\n",
				*size, lsize);
		return -1;
	}
	*size = lsize;
	size3 = *size * 3;
	if(size3>xfp->buf1size) 
    {
		if((xfp->buf1=(int *)malloc(4096 + sizeof(int)*size3))==NULL) 
        {
			fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
			return -1; 
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(4096 + sizeof(int)*xfp->buf2size))==NULL)
        {
			fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
			return -1;
		}
	}
	/* Dont bother with compression for three atoms or less */
	if(*size<=9) 
    {
		return xdrfile_read_float(ptr,size3,xfp)/3;
		/* return number of coords, not floats */
	}
	/* Compression-time if we got here. Read precision first */
	xdrfile_read_float(precision,1,xfp);
  
	/* avoid repeated pointer dereferencing. */
	buf1=xfp->buf1; 
	buf2=xfp->buf2;
	/* buf2[0-2] are special and do not contain actual data */
	buf2[0] = buf2[1] = buf2[2] = 0;
	xdrfile_read_int(minint,3,xfp);
	xdrfile_read_int(maxint,3,xfp);
  
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;

	if( !sizeint[0] || !sizeint[1] || !sizeint[2] ) { 
		//fprintf( stderr, "XTC file corrupt\n" );
		return exdrDUFF;
	}
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
    {
		bitsizeint[0] = sizeofint(sizeint[0]);
		bitsizeint[1] = sizeofint(sizeint[1]);
		bitsizeint[2] = sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	}
    else 
    {
		bitsize = sizeofints(3, sizeint);
	}
	
	if (xdrfile_read_int(&smallidx,1,xfp) == 0)	
		return 0; /* not sure what has happened here or why we return... */
	tmp=smallidx+8;
	maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
	minidx = maxidx - 8; /* often this equal smallidx */
	tmp = smallidx-1;
	tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
	smaller = magicints[tmp] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	larger = magicints[maxidx];

	if( !sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) {
		//fprintf(stderr, "XTC file corrupt\n");
		return exdrDUFF;
	}
	/* buf2[0] holds the length in bytes */
  
	if (xdrfile_read_int(buf2,1,xfp) == 0)
		return 0;
	if (xdrfile_read_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp) == 0)
		return 0;
	buf2[0] = buf2[1] = buf2[2] = 0;
  
	lfp = ptr;
	inv_precision = 1.0 / * precision;
	run = 0;
	i = 0;
	lip = buf1;
	while ( i < lsize ) 
    {
		thiscoord = (int *)(lip) + i * 3;
    
		if (bitsize == 0) 
        {
			thiscoord[0] = decodebits(buf2, bitsizeint[0]);
			thiscoord[1] = decodebits(buf2, bitsizeint[1]);
			thiscoord[2] = decodebits(buf2, bitsizeint[2]);
		}
        else
        {
			decodeints(buf2, 3, bitsize, sizeint, thiscoord);
		}
    
		i++;
		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];
    
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
    
		flag = decodebits(buf2, 1);
		is_smaller = 0;
		if (flag == 1) 
        {
			run = decodebits(buf2, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
		if (run > 0)
        {
			thiscoord += 3;
			for (k = 0; k < run; k+=3) 
            {
				decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
				i++;
				thiscoord[0] += prevcoord[0] - smallnum;
				thiscoord[1] += prevcoord[1] - smallnum;
				thiscoord[2] += prevcoord[2] - smallnum;
				if (k == 0) {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					*lfp++ = prevcoord[0] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
					*lfp++ = prevcoord[1] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
					*lfp++ = prevcoord[2] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
				} else {
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
				}
				if(countmjh >= size3 ) { return -1; }
				*lfp++ = thiscoord[0] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
				if(countmjh >= size3 ) { return -1; }
				*lfp++ = thiscoord[1] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
				if(countmjh >= size3 ) { return -1; }
				*lfp++ = thiscoord[2] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
			}
		} 
        else
        {
				if(countmjh >= size3 ) { return -1; }
			*lfp++ = thiscoord[0] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
				if(countmjh >= size3 ) { return -1; }
			*lfp++ = thiscoord[1] * inv_precision;
					countmjh++; if( countmjh> 3* *size) { return -1; }
				if(countmjh >= size3 ) { return -1; }
			*lfp++ = thiscoord[2] * inv_precision;		
					countmjh++; if( countmjh> 3* *size) { return -1; }
		}
		smallidx += is_smaller;
		if (is_smaller < 0) 
        {
			smallnum = smaller;
            
			if (smallidx > FIRSTIDX) 
            {
				smaller = magicints[smallidx - 1] /2;
			} 
            else 
            {
				smaller = 0;
			}
		} 
        else if (is_smaller > 0)
        {
			smaller = smallnum;
			smallnum = magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;

		if(!sizesmall[0] || !sizesmall[1] || !sizesmall[2] ) { 
			//fprintf( stderr, "XTC file corrupt\n");
			return exdrDUFF;
		}
	}
	return *size;
}

int
xdrfile_compress_coord_float(float   *ptr,
							 int      size,
							 float    precision,
							 XDRFILE* xfp)
{
	int minint[3], maxint[3], mindiff, *lip, diff;
	int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
	int minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
	int k, *buf1, *buf2;
	int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
	float *lfp, lf;
	int tmp=0, tmpsum=0, *thiscoord=NULL,  prevcoord[3]={0,0,0};
	unsigned int tmpcoord[30];
	unsigned int bitsize;
	int countmjh=0;
	int errval;
  
	if(xfp==NULL)
		return -1;
	size3=3*size;
    
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

	if(size3>xfp->buf1size)
    {
		if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL) 
        {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL)
        {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
	}
	if(xdrfile_write_int(&size,1,xfp)==0)
		return -1; /* return if we could not write size */
	/* Dont bother with compression for three atoms or less */
	if(size<=9) 
    {
		return xdrfile_write_float(ptr,size3,xfp)/3;
		/* return number of coords, not floats */
	}
	/* Compression-time if we got here. Write precision first */
	if (precision <= 0)
		precision = 1000;
	xdrfile_write_float(&precision,1,xfp);
	/* avoid repeated pointer dereferencing. */
	buf1=xfp->buf1; 
	buf2=xfp->buf2;
	/* buf2[0-2] are special and do not contain actual data */
	buf2[0] = buf2[1] = buf2[2] = 0;
	minint[0] = minint[1] = minint[2] = INT_MAX;
	maxint[0] = maxint[1] = maxint[2] = INT_MIN;
	prevrun = -1;
	lfp = ptr;
	lip = buf1;
	mindiff = INT_MAX;
	oldlint1 = oldlint2 = oldlint3 = 0;
	while(lfp < ptr + size3 )
    {
		/* find nearest integer */
		if (*lfp >= 0.0)
			lf = *lfp * precision + 0.5;
		else
			lf = *lfp * precision - 0.5;
		if (fabs(lf) > INT_MAX-2) 
        {
			/* scaling would cause overflow */
			fprintf(stderr,"Internal overflow compressing coordinates.\n");
			errval=0;
		}
		lint1 = lf;
		if (lint1 < minint[0]) minint[0] = lint1;
		if (lint1 > maxint[0]) maxint[0] = lint1;
		*lip++ = lint1;
		lfp++;
	  	countmjh++; if( countmjh> size3) { return -1; }
		if (*lfp >= 0.0)
			lf = *lfp * precision + 0.5;
		else
			lf = *lfp * precision - 0.5;
		if (fabs(lf) > INT_MAX-2)
        {
			/* scaling would cause overflow */
			fprintf(stderr,"Internal overflow compressing coordinates.\n");
			errval=0;
		}
		lint2 = lf;
		if (lint2 < minint[1]) minint[1] = lint2;
		if (lint2 > maxint[1]) maxint[1] = lint2;
		*lip++ = lint2;
		lfp++;
	  countmjh++; if( countmjh> size3) { return -1; }
		if (*lfp >= 0.0)
			lf = *lfp * precision + 0.5;
		else
			lf = *lfp * precision - 0.5;
		if (fabs(lf) > INT_MAX-2) 
        {
			errval=0;      
		}
		lint3 = lf;
		if (lint3 < minint[2]) minint[2] = lint3;
		if (lint3 > maxint[2]) maxint[2] = lint3;
		*lip++ = lint3;
		lfp++;
	  countmjh++; if( countmjh> size3) { return -1; }
		diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
		if (diff < mindiff && lfp > ptr + 3)
			mindiff = diff;
		oldlint1 = lint1;
		oldlint2 = lint2;
		oldlint3 = lint3;
	}  
	xdrfile_write_int(minint,3,xfp);
	xdrfile_write_int(maxint,3,xfp);
  
	if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
		(float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
		(float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
		/* turning value in unsigned by subtracting minint
		 * would cause overflow
		 */
		fprintf(stderr,"Internal overflow compressing coordinates.\n");
		errval=0;
	}
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
  
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
    {
		bitsizeint[0] = sizeofint(sizeint[0]);
		bitsizeint[1] = sizeofint(sizeint[1]);
		bitsizeint[2] = sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	}
    else
    {
		bitsize = sizeofints(3, sizeint);
	}
	lip = buf1;
	luip = (unsigned int *) buf1;
	smallidx = FIRSTIDX;
	while (smallidx < LASTIDX && magicints[smallidx] < mindiff)
    {
		smallidx++;
	}
	xdrfile_write_int(&smallidx,1,xfp);
	tmp=smallidx+8;
	maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
	minidx = maxidx - 8; /* often this equal smallidx */
	tmp=smallidx-1;
	tmp= (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
	smaller = magicints[tmp] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	larger = magicints[maxidx] / 2;
	i = 0;
	while (i < size) 
    {
		is_small = 0;
		thiscoord = (int *)(luip) + i * 3;
		if (smallidx < maxidx && i >= 1 &&
			abs(thiscoord[0] - prevcoord[0]) < larger &&
			abs(thiscoord[1] - prevcoord[1]) < larger &&
			abs(thiscoord[2] - prevcoord[2]) < larger) {
			is_smaller = 1;
		} 
        else if (smallidx > minidx) 
        {
			is_smaller = -1;
		}
        else
        {
			is_smaller = 0;
		}
		if (i + 1 < size) 
        {
			if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
				abs(thiscoord[1] - thiscoord[4]) < smallnum &&
				abs(thiscoord[2] - thiscoord[5]) < smallnum) 
            {
				/* interchange first with second atom for better
				 * compression of water molecules
				 */
				tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
				thiscoord[3] = tmp;
				tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
				thiscoord[4] = tmp;
				tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
				thiscoord[5] = tmp;
				is_small = 1;
			} 
		}
		tmpcoord[0] = thiscoord[0] - minint[0];
		tmpcoord[1] = thiscoord[1] - minint[1];
		tmpcoord[2] = thiscoord[2] - minint[2];
		if (bitsize == 0) 
        {
			encodebits(buf2, bitsizeint[0], tmpcoord[0]);
			encodebits(buf2, bitsizeint[1], tmpcoord[1]);
			encodebits(buf2, bitsizeint[2], tmpcoord[2]);
		} 
        else
        {
			encodeints(buf2, 3, bitsize, sizeint, tmpcoord);
		}
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
		thiscoord = thiscoord + 3;
		i++;

		run = 0;
		if (is_small == 0 && is_smaller == -1)
			is_smaller = 0;
		while (is_small && run < 8*3)
        {
			tmpsum=0;
			for(j=0;j<3;j++) 
            {
				tmp=thiscoord[j] - prevcoord[j];
				tmpsum+=tmp*tmp;
			}
			if (is_smaller == -1 && tmpsum >= smaller * smaller)
            {
				is_smaller = 0;
			}
      
			tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
			tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
			tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;
      
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
      
			i++;
			thiscoord = thiscoord + 3;
			is_small = 0;
			if (i < size &&
				abs(thiscoord[0] - prevcoord[0]) < smallnum &&
				abs(thiscoord[1] - prevcoord[1]) < smallnum &&
				abs(thiscoord[2] - prevcoord[2]) < smallnum)
            {
				is_small = 1;
			}
		}
		if (run != prevrun || is_smaller != 0) 
        {
			prevrun = run;
			encodebits(buf2, 1, 1); /* flag the change in run-length */
			encodebits(buf2, 5, run+is_smaller+1);
		} 
        else 
        {
			encodebits(buf2, 1, 0); /* flag the fact that runlength did not change */
		}
		for (k=0; k < run; k+=3) 
        {
			encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]);	
		}
		if (is_smaller != 0) 
        {
			smallidx += is_smaller;
			if (is_smaller < 0) 
            {
				smallnum = smaller;
				smaller = magicints[smallidx-1] / 2;
			} 
            else 
            {
				smaller = smallnum;
				smallnum = magicints[smallidx] / 2;
			}
			sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
		}   
	}
	if (buf2[1] != 0) buf2[0]++;
	xdrfile_write_int(buf2,1,xfp); /* buf2[0] holds the length in bytes */
	tmp=xdrfile_write_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp);
	if(tmp==(unsigned int)buf2[0])
		return size;
	else
		return -1;
}


int
xdrfile_decompress_coord_double(double     *ptr, 
								int        *size,
								double     *precision,
								XDRFILE*   xfp)
{
	int minint[3], maxint[3], *lip;
	int smallidx,  maxidx, minidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
	int k, *buf1, *buf2, lsize, flag;
	int smallnum, smaller, larger, i, is_smaller, run;
	double *lfp, inv_precision;
	float float_prec, tmpdata[30];
	int tmp, *thiscoord,  prevcoord[3];
	unsigned int bitsize;
	int errval=1;
  
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

	if(xfp==NULL || ptr==NULL)
		return -1;
	tmp=xdrfile_read_int(&lsize,1,xfp);
	if(tmp==0)
		return -1; /* return if we could not read size */
	if (*size < lsize) 
    {
		fprintf(stderr, "Requested to decompress %d coords, file contains %d\n",
				*size, lsize);
		return -1;
	}
	*size = lsize;
	size3 = *size * 3;
	if(size3>xfp->buf1size) 
    {
		if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL) 
        {
			fprintf(stderr,"Cannot allocate memory for decompression coordinates.\n");
			return -1; 
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL)
        {
			fprintf(stderr,"Cannot allocate memory for decompressing coordinates.\n");
			return -1;
		}
	}
	/* Dont bother with compression for three atoms or less */
	if(*size<=9)
    {
		tmp=xdrfile_read_float(tmpdata,size3,xfp);
		for(i=0;i<9*3;i++)
			ptr[i]=tmpdata[i];
		return tmp/3;
		/* return number of coords, not floats */
	}
	/* Compression-time if we got here. Read precision first */
	xdrfile_read_float(&float_prec,1,xfp);
	*precision=float_prec;
	/* avoid repeated pointer dereferencing. */
	buf1=xfp->buf1; 
	buf2=xfp->buf2;
	/* buf2[0-2] are special and do not contain actual data */
	buf2[0] = buf2[1] = buf2[2] = 0;
	xdrfile_read_int(minint,3,xfp);
	xdrfile_read_int(maxint,3,xfp);
  
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;


	if( !sizeint[0] || !sizeint[1] || !sizeint[2] ) { 
		//fprintf( stderr, "XTC file corrupt\n" );
		return exdrDUFF;
	}
	
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
    {
		bitsizeint[0] = sizeofint(sizeint[0]);
		bitsizeint[1] = sizeofint(sizeint[1]);
		bitsizeint[2] = sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	}
    else 
    {
		bitsize = sizeofints(3, sizeint);
	}
	
	if (xdrfile_read_int(&smallidx,1,xfp) == 0)	
		return 0;
	tmp=smallidx+8;
	maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
	minidx = maxidx - 8; /* often this equal smallidx */
	tmp = smallidx-1;
	tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
	smaller = magicints[tmp] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	larger = magicints[maxidx];

	/* buf2[0] holds the length in bytes */
  
	if (xdrfile_read_int(buf2,1,xfp) == 0)
		return 0;
	if (xdrfile_read_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp) == 0)
		return 0;
	buf2[0] = buf2[1] = buf2[2] = 0;
  
	lfp = ptr;
	inv_precision = 1.0 / * precision;
	run = 0;
	i = 0;
	lip = buf1;
	while ( i < lsize ) 
    {
		thiscoord = (int *)(lip) + i * 3;
    
		if (bitsize == 0) 
        {
			thiscoord[0] = decodebits(buf2, bitsizeint[0]);
			thiscoord[1] = decodebits(buf2, bitsizeint[1]);
			thiscoord[2] = decodebits(buf2, bitsizeint[2]);
		} else {
			decodeints(buf2, 3, bitsize, sizeint, thiscoord);
		}
    
		i++;
		thiscoord[0] += minint[0];
		thiscoord[1] += minint[1];
		thiscoord[2] += minint[2];
    
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
    
		flag = decodebits(buf2, 1);
		is_smaller = 0;
		if (flag == 1) 
        {
			run = decodebits(buf2, 5);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}
		if (run > 0) 
        {
			thiscoord += 3;
			for (k = 0; k < run; k+=3) 
            {
				decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
				i++;
				thiscoord[0] += prevcoord[0] - smallnum;
				thiscoord[1] += prevcoord[1] - smallnum;
				thiscoord[2] += prevcoord[2] - smallnum;
				if (k == 0)
                {
					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
					prevcoord[0] = tmp;
					tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
					prevcoord[1] = tmp;
					tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
					prevcoord[2] = tmp;
					*lfp++ = prevcoord[0] * inv_precision;
					*lfp++ = prevcoord[1] * inv_precision;
					*lfp++ = prevcoord[2] * inv_precision;
				}
                else 
                {
					prevcoord[0] = thiscoord[0];
					prevcoord[1] = thiscoord[1];
					prevcoord[2] = thiscoord[2];
				}
				*lfp++ = thiscoord[0] * inv_precision;
				*lfp++ = thiscoord[1] * inv_precision;
				*lfp++ = thiscoord[2] * inv_precision;
			}
		} else {
			*lfp++ = thiscoord[0] * inv_precision;
			*lfp++ = thiscoord[1] * inv_precision;
			*lfp++ = thiscoord[2] * inv_precision;		
		}
		smallidx += is_smaller;
		if (is_smaller < 0) {
			smallnum = smaller;
			if (smallidx > FIRSTIDX) {
				smaller = magicints[smallidx - 1] /2;
			} else {
				smaller = 0;
			}
		} else if (is_smaller > 0) {
			smaller = smallnum;
			smallnum = magicints[smallidx] / 2;
		}
		sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
	}
	return *size;
}

int
xdrfile_compress_coord_double(double   *ptr,
							  int      size,
							  double    precision,
							  XDRFILE* xfp)
{
	int minint[3], maxint[3], mindiff, *lip, diff;
	int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
	int minidx, maxidx;
	unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
	int k, *buf1, *buf2;
	int smallnum, smaller, larger, i, j, is_small, is_smaller, run, prevrun;
	double *lfp;
	float float_prec, lf,tmpdata[30];
	int tmp=0, tmpsum=0, *thiscoord=NULL,  prevcoord[3]={0,0,0};
	unsigned int tmpcoord[30];
	unsigned int bitsize;
	int errval=1;
  
    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

	if(xfp==NULL)
		return -1;
	size3=3*size;
	if(size3>xfp->buf1size) {
		if((xfp->buf1=(int *)malloc(sizeof(int)*size3))==NULL) {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
		xfp->buf1size=size3;
		xfp->buf2size=size3*1.2;
		if((xfp->buf2=(int *)malloc(sizeof(int)*xfp->buf2size))==NULL) {
			fprintf(stderr,"Cannot allocate memory for compressing coordinates.\n");
			return -1;
		}
	}
	if(xdrfile_write_int(&size,1,xfp)==0)
		return -1; /* return if we could not write size */
	/* Dont bother with compression for three atoms or less */
	if(size<=9) {
		for(i=0;i<9*3;i++)
			tmpdata[i]=ptr[i];
		return xdrfile_write_float(tmpdata,size3,xfp)/3;
		/* return number of coords, not floats */
	}
	/* Compression-time if we got here. Write precision first */
	if (precision <= 0)
		precision = 1000;
	float_prec=precision;
	xdrfile_write_float(&float_prec,1,xfp);
	/* avoid repeated pointer dereferencing. */
	buf1=xfp->buf1; 
	buf2=xfp->buf2;
	/* buf2[0-2] are special and do not contain actual data */
	buf2[0] = buf2[1] = buf2[2] = 0;
	minint[0] = minint[1] = minint[2] = INT_MAX;
	maxint[0] = maxint[1] = maxint[2] = INT_MIN;
	prevrun = -1;
	lfp = ptr;
	lip = buf1;
	mindiff = INT_MAX;
	oldlint1 = oldlint2 = oldlint3 = 0;
	while(lfp < ptr + size3 ) {
		/* find nearest integer */
		if (*lfp >= 0.0)
			lf = (float)*lfp * float_prec + 0.5;
		else
			lf = (float)*lfp * float_prec - 0.5;
		if (fabs(lf) > INT_MAX-2) {
			/* scaling would cause overflow */
			fprintf(stderr,"Internal overflow compressing coordinates.\n");
			errval=0;
		}
		lint1 = lf;
		if (lint1 < minint[0]) minint[0] = lint1;
		if (lint1 > maxint[0]) maxint[0] = lint1;
		*lip++ = lint1;
		lfp++;
		if (*lfp >= 0.0)
			lf = (float)*lfp * float_prec + 0.5;
		else
			lf = (float)*lfp * float_prec - 0.5;
		if (fabs(lf) > INT_MAX-2) {
			/* scaling would cause overflow */
			fprintf(stderr,"Internal overflow compressing coordinates.\n");
			errval=0;
		}
		lint2 = lf;
		if (lint2 < minint[1]) minint[1] = lint2;
		if (lint2 > maxint[1]) maxint[1] = lint2;
		*lip++ = lint2;
		lfp++;
		if (*lfp >= 0.0)
			lf = (float)*lfp * float_prec + 0.5;
		else
			lf = (float)*lfp * float_prec - 0.5;
		if (fabs(lf) > INT_MAX-2) {
			errval=0;      
		}
		lint3 = lf;
		if (lint3 < minint[2]) minint[2] = lint3;
		if (lint3 > maxint[2]) maxint[2] = lint3;
		*lip++ = lint3;
		lfp++;
		diff = abs(oldlint1-lint1)+abs(oldlint2-lint2)+abs(oldlint3-lint3);
		if (diff < mindiff && lfp > ptr + 3)
			mindiff = diff;
		oldlint1 = lint1;
		oldlint2 = lint2;
		oldlint3 = lint3;
	}  
	xdrfile_write_int(minint,3,xfp);
	xdrfile_write_int(maxint,3,xfp);
  
	if ((float)maxint[0] - (float)minint[0] >= INT_MAX-2 ||
		(float)maxint[1] - (float)minint[1] >= INT_MAX-2 ||
		(float)maxint[2] - (float)minint[2] >= INT_MAX-2) {
		/* turning value in unsigned by subtracting minint
		 * would cause overflow
		 */
		fprintf(stderr,"Internal overflow compressing coordinates.\n");
		errval=0;
	}
	sizeint[0] = maxint[0] - minint[0]+1;
	sizeint[1] = maxint[1] - minint[1]+1;
	sizeint[2] = maxint[2] - minint[2]+1;
  
	/* check if one of the sizes is to big to be multiplied */
	if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
		bitsizeint[0] = sizeofint(sizeint[0]);
		bitsizeint[1] = sizeofint(sizeint[1]);
		bitsizeint[2] = sizeofint(sizeint[2]);
		bitsize = 0; /* flag the use of large sizes */
	} else {
		bitsize = sizeofints(3, sizeint);
	}
	lip = buf1;
	luip = (unsigned int *) buf1;
	smallidx = FIRSTIDX;
	while (smallidx < LASTIDX && magicints[smallidx] < mindiff) {
		smallidx++;
	}
	xdrfile_write_int(&smallidx,1,xfp);
	tmp=smallidx+8;
	maxidx = (LASTIDX<tmp) ? LASTIDX : tmp;
	minidx = maxidx - 8; /* often this equal smallidx */
	tmp=smallidx-1;
	tmp= (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
	smaller = magicints[tmp] / 2;
	smallnum = magicints[smallidx] / 2;
	sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
	larger = magicints[maxidx] / 2;
	i = 0;
	while (i < size) {
		is_small = 0;
		thiscoord = (int *)(luip) + i * 3;
		if (smallidx < maxidx && i >= 1 &&
			abs(thiscoord[0] - prevcoord[0]) < larger &&
			abs(thiscoord[1] - prevcoord[1]) < larger &&
			abs(thiscoord[2] - prevcoord[2]) < larger) {
			is_smaller = 1;
		} else if (smallidx > minidx) {
			is_smaller = -1;
		} else {
			is_smaller = 0;
		}
		if (i + 1 < size) {
			if (abs(thiscoord[0] - thiscoord[3]) < smallnum &&
				abs(thiscoord[1] - thiscoord[4]) < smallnum &&
				abs(thiscoord[2] - thiscoord[5]) < smallnum) {
				/* interchange first with second atom for better
				 * compression of water molecules
				 */
				tmp = thiscoord[0]; thiscoord[0] = thiscoord[3];
				thiscoord[3] = tmp;
				tmp = thiscoord[1]; thiscoord[1] = thiscoord[4];
				thiscoord[4] = tmp;
				tmp = thiscoord[2]; thiscoord[2] = thiscoord[5];
				thiscoord[5] = tmp;
				is_small = 1;
			} 
		}
		tmpcoord[0] = thiscoord[0] - minint[0];
		tmpcoord[1] = thiscoord[1] - minint[1];
		tmpcoord[2] = thiscoord[2] - minint[2];
		if (bitsize == 0) {
			encodebits(buf2, bitsizeint[0], tmpcoord[0]);
			encodebits(buf2, bitsizeint[1], tmpcoord[1]);
			encodebits(buf2, bitsizeint[2], tmpcoord[2]);
		} else {
			encodeints(buf2, 3, bitsize, sizeint, tmpcoord);
		}
		prevcoord[0] = thiscoord[0];
		prevcoord[1] = thiscoord[1];
		prevcoord[2] = thiscoord[2];
		thiscoord = thiscoord + 3;
		i++;

		run = 0;
		if (is_small == 0 && is_smaller == -1)
			is_smaller = 0;
		while (is_small && run < 8*3) {
			tmpsum=0;
			for(j=0;j<3;j++) {
				tmp=thiscoord[j] - prevcoord[j];
				tmpsum+=tmp*tmp;
			}
			if (is_smaller == -1 && tmpsum >= smaller * smaller) {
				is_smaller = 0;
			}
      
			tmpcoord[run++] = thiscoord[0] - prevcoord[0] + smallnum;
			tmpcoord[run++] = thiscoord[1] - prevcoord[1] + smallnum;
			tmpcoord[run++] = thiscoord[2] - prevcoord[2] + smallnum;
      
			prevcoord[0] = thiscoord[0];
			prevcoord[1] = thiscoord[1];
			prevcoord[2] = thiscoord[2];
      
			i++;
			thiscoord = thiscoord + 3;
			is_small = 0;
			if (i < size &&
				abs(thiscoord[0] - prevcoord[0]) < smallnum &&
				abs(thiscoord[1] - prevcoord[1]) < smallnum &&
				abs(thiscoord[2] - prevcoord[2]) < smallnum) {
				is_small = 1;
			}
		}
		if (run != prevrun || is_smaller != 0) {
			prevrun = run;
			encodebits(buf2, 1, 1); /* flag the change in run-length */
			encodebits(buf2, 5, run+is_smaller+1);
		} else {
			encodebits(buf2, 1, 0); /* flag the fact that runlength did not change */
		}
		for (k=0; k < run; k+=3) {
			encodeints(buf2, 3, smallidx, sizesmall, &tmpcoord[k]);	
		}
		if (is_smaller != 0) {
			smallidx += is_smaller;
			if (is_smaller < 0) {
				smallnum = smaller;
				smaller = magicints[smallidx-1] / 2;
			} else {
				smaller = smallnum;
				smallnum = magicints[smallidx] / 2;
			}
			sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
		}   
	}
	if (buf2[1] != 0) buf2[0]++;
	xdrfile_write_int(buf2,1,xfp); /* buf2[0] holds the length in bytes */
	tmp=xdrfile_write_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp);
	if(tmp==(unsigned int)buf2[0])
		return size;
	else
		return -1; 
}


/* Dont try do document Fortran interface, since
 * Doxygen barfs at the F77_FUNC macro 
 */
#ifndef DOXYGEN 

/*************************************************************
 * Fortran77 interface for reading/writing portable data     *
 * The routine are not threadsafe when called from Fortran   *
 * (as they are when called from C) unless you compile with  *
 * this file with posix thread support.                      *
 * Note that these are not multithread-safe.                 * 
 *************************************************************/
#define MAX_FORTRAN_XDR 1024
static XDRFILE *f77xdr[MAX_FORTRAN_XDR]; /* array of file handles */
static int      f77init = 1;             /* zero array first time */

/* internal to this file: C<-->Fortran string conversion */
static int ftocstr(char *dest, int dest_len, char *src, int src_len);
static int ctofstr(char *dest, int dest_len, char *src);


void
F77_FUNC(xdropen,XDROPEN)(int *fid, char *filename, char *mode,
						  int fn_len, int mode_len)
{
	char cfilename[512];
	char cmode[5];
	int i;
  
	/* zero array at first invocation */
	if(f77init) {
		for(i=0;i<MAX_FORTRAN_XDR;i++)
			f77xdr[i]=NULL;
		f77init=0;
	}
	i=0;
  
	/* nf77xdr is always smaller or equal to MAX_FORTRAN_XDR */
	while(i<MAX_FORTRAN_XDR && f77xdr[i]!=NULL)
		i++;
	if(i==MAX_FORTRAN_XDR) {
		*fid = -1;
	} else if (ftocstr(cfilename, sizeof(cfilename), filename, fn_len)) {
		*fid = -1;
	} else if (ftocstr(cmode, sizeof(cmode), mode,mode_len)) {
		*fid = -1;
	} else {
		f77xdr[i]=xdrfile_open(cfilename,cmode);
		/* return the index in the array as a fortran file handle */
		*fid=i; 
	}
}

void
F77_FUNC(xdrclose,XDRCLOSE)(int *fid)
{
    /* first close it */
    xdrfile_close(f77xdr[*fid]);
    /* the remove it from file handle list */
    f77xdr[*fid]=NULL;
}


void
F77_FUNC(xdrrint,XDRRINT)(int *fid, int *data, int *ndata, int *ret)
{
	*ret = xdrfile_read_int(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwint,XDRWINT)(int *fid, int *data, int *ndata, int *ret)
{
	*ret = xdrfile_write_int(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrruint,XDRRUINT)(int *fid, unsigned int *data, int *ndata, int *ret)
{
	*ret = xdrfile_read_uint(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwuint,XDRWUINT)(int *fid, unsigned int *data, int *ndata, int *ret)
{
	*ret = xdrfile_write_uint(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrrchar,XDRRCHAR)(int *fid, char *ip, int *ndata, int *ret)
{
	*ret = xdrfile_read_char(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwchar,XDRWCHAR)(int *fid, char *ip, int *ndata, int *ret)
{
	*ret = xdrfile_write_char(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrruchar,XDRRUCHAR)(int *fid, unsigned char *ip, int *ndata, int *ret)
{
	*ret = xdrfile_read_uchar(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwuchar,XDRWUCHAR)(int *fid, unsigned char *ip, int *ndata, int *ret)
{
	*ret = xdrfile_write_uchar(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrrshort,XDRRSHORT)(int *fid, short *ip, int *ndata, int *ret)
{
	*ret = xdrfile_read_short(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwshort,XDRWSHORT)(int *fid, short *ip, int *ndata, int *ret)
{
	*ret = xdrfile_write_short(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrrushort,XDRRUSHORT)(int *fid, unsigned short *ip, int *ndata, int *ret)
{
	*ret = xdrfile_read_ushort(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwushort,XDRWUSHORT)(int *fid, unsigned short *ip, int *ndata, int *ret)
{
	*ret = xdrfile_write_ushort(ip,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrrsingle,XDRRSINGLE)(int *fid, float *data, int *ndata, int *ret)
{
	*ret = xdrfile_read_float(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwsingle,XDRWSINGLE)(int *fid, float *data, int *ndata, int *ret)
{
	*ret = xdrfile_write_float(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrrdouble,XDRRDOUBLE)(int *fid, double *data, int *ndata, int *ret)
{
	*ret = xdrfile_read_double(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwdouble,XDRWDOUBLE)(int *fid, double *data, int *ndata, int *ret)
{
	*ret = xdrfile_write_double(data,*ndata,f77xdr[*fid]);
}

static int ftocstr(char *dest, int destlen, char *src, int srclen)
{
    char *p;

    p = src + srclen;
    while ( --p >= src && *p == ' ' );
    srclen = p - src + 1;
    destlen--;
    dest[0] = 0;
    if (srclen > destlen)
		return 1;
    while (srclen--)
		(*dest++ = *src++);
    *dest = '\0';
    return 0;
}


static int ctofstr(char *dest, int destlen, char *src)
{
    while (destlen && *src) {
        *dest++ = *src++;
        destlen--;
    }
    while (destlen--)
        *dest++ = ' ';
    return 0;
}


void
F77_FUNC(xdrrstring,XDRRSTRING)(int *fid, char *str, int *ret, int len)
{
	char *cstr;
  
	if((cstr=(char*)malloc((len+1)*sizeof(char)))==NULL) {
		*ret = 0;
		return;
	}
	if (ftocstr(cstr, len+1, str, len)) {
		*ret = 0;
		cstr=(char *)condfree(cstr);
		return;
	}
  
	*ret = xdrfile_read_string(cstr, len+1,f77xdr[*fid]);
	ctofstr( str, len , cstr);
	cstr=(char *)condfree(cstr);
}

void
F77_FUNC(xdrwstring,XDRWSTRING)(int *fid, char *str, int *ret, int len)
{
	char *cstr;
  
	if((cstr=(char*)malloc((len+1)*sizeof(char)))==NULL) {
		*ret = 0;
		return;
	}
	if (ftocstr(cstr, len+1, str, len)) {
		*ret = 0;
		cstr=(char *)condfree(cstr);
		return;
	}
  
	*ret = xdrfile_write_string(cstr, f77xdr[*fid]);
	ctofstr( str, len , cstr);
	cstr=(char *)condfree(cstr);
}

void
F77_FUNC(xdrropaque,XDRROPAQUE)(int *fid, char *data, int *ndata, int *ret)
{
	*ret = xdrfile_read_opaque(data,*ndata,f77xdr[*fid]);
}

void
F77_FUNC(xdrwopaque,XDRWOPAQUE)(int *fid, char *data, int *ndata, int *ret)
{
	*ret = xdrfile_write_opaque(data,*ndata,f77xdr[*fid]);
}


/* Write single-precision compressed 3d coordinates */
void
F77_FUNC(xdrccs,XDRCCS)(int *fid, float *data, int *ncoord, 
						float *precision, int *ret)
{
    *ret = xdrfile_compress_coord_float(data,*ncoord,*precision,f77xdr[*fid]);
}


/* Read single-precision compressed 3d coordinates */
void
F77_FUNC(xdrdcs,XDRDCS)(int *fid, float *data, int *ncoord, 
						float *precision, int *ret)
{
	*ret = xdrfile_decompress_coord_float(data,ncoord,precision,f77xdr[*fid]);
}


/* Write compressed 3d coordinates from double precision data */
void
F77_FUNC(xdrccd,XDRCCD)(int *fid, double *data, int *ncoord, 
						double *precision, int *ret)
{
	*ret = xdrfile_compress_coord_double(data,*ncoord,*precision,f77xdr[*fid]);
}

/* Read compressed 3d coordinates into double precision data */
void
F77_FUNC(xddcd,XDRDCD)(int *fid, double *data, int *ncoord, 
					   double *precision, int *ret)
{
    *ret = xdrfile_decompress_coord_double(data,ncoord,precision,f77xdr[*fid]);
}







#endif /* DOXYGEN */

/*************************************************************
 * End of higher-level routines - dont change things below!  *
 *************************************************************/




















/*************************************************************
 * The rest of this file contains our own implementation     *
 * of the XDR calls in case you are compiling without them.  *
 * You do NOT want to change things here since it would make *
 * things incompatible with the standard RPC/XDR routines.   *
 *************************************************************/
#ifndef HAVE_RPC_XDR_H

/*
 * What follows is a modified version of the Sun XDR code. For reference
 * we include their copyright and license:
 * 
 * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
 * unrestricted use provided that this legend is included on all tape
 * media and as a part of the software program in whole or part.  Users
 * may copy or modify Sun RPC without charge, but are not authorized
 * to license or distribute it to anyone else except as part of a product or
 * program developed by the user.
 *
 * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
 * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
 *
 * Sun RPC is provided with no support and without any obligation on the
 * part of Sun Microsystems, Inc. to assist in its use, correction,
 * modification or enhancement.
 *
 * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
 * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
 * OR ANY PART THEREOF.
 *
 * In no event will Sun Microsystems, Inc. be liable for any lost revenue
 * or profits or other special, indirect and consequential damages, even if
 * Sun has been advised of the possibility of such damages.
 *
 * Sun Microsystems, Inc.
 * 2550 Garcia Avenue
 * Mountain View, California  94043
 */

/* INT_MAX is defined in limits.h according to ANSI C */
#if (INT_MAX > 2147483647)
#    error Error: Cannot use builtin XDR support when size of int
#    error is larger than 4 bytes. Use your system XDR libraries 
#    error instead, or modify the source code in xdrfile.c
#endif /* Check for 4 byte int type */





typedef int (*xdrproc_t) (XDR *, void *,...);

#define xdr_getlong(xdrs, longp)			\
	(*(xdrs)->x_ops->x_getlong)(xdrs, longp)
#define xdr_putlong(xdrs, longp)			\
	(*(xdrs)->x_ops->x_putlong)(xdrs, longp)
#define xdr_getbytes(xdrs, addr, len)			\
	(*(xdrs)->x_ops->x_getbytes)(xdrs, addr, len)
#define xdr_putbytes(xdrs, addr, len)			\
	(*(xdrs)->x_ops->x_putbytes)(xdrs, addr, len)

#define BYTES_PER_XDR_UNIT 4 

// Apple clang does not support thread_local
#ifdef __APPLE__
static __thread char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};
#else
thread_local char xdr_zero[BYTES_PER_XDR_UNIT] = {0, 0, 0, 0};
#endif


static int32_t
xdr_swapbytes(int32_t x)
{
	int32_t y,i;
	char *px=(char *)&x;
	char *py=(char *)&y;
  
	for(i=0;i<4;i++)
		py[i]=px[3-i];
  
	return y;
}

static int32_t
xdr_htonl(int32_t x)
{
	int s=0x1234;
	if( *((char *)&s)==(char)0x34) {
		/* smallendian,swap bytes */
		return xdr_swapbytes(x);
	} else {
		/* bigendian, do nothing */
		return x;
	}
}

static int32_t
xdr_ntohl(int x)
{
	int s=0x1234;
	if( *((char *)&s)==(char)0x34) {
		/* smallendian, swap bytes */
		return xdr_swapbytes(x);
	} else {
		/* bigendian, do nothing */
		return x;
	}
}

static int
xdr_int (XDR *xdrs, int *ip)
{
	int32_t i32;

	switch (xdrs->x_op)
		{
		case XDR_ENCODE:
			i32 = (int32_t) *ip;
			return xdr_putlong (xdrs, &i32);

		case XDR_DECODE:
			if (!xdr_getlong (xdrs, &i32))
				{
					return 0;
				}
			*ip = (int) i32;
		case XDR_FREE:
			return 1;
		}
	return 0;
}

static int
xdr_u_int (XDR *xdrs, unsigned int *up)
{
	uint32_t ui32;

	switch (xdrs->x_op)
		{
		case XDR_ENCODE:
			ui32 = (uint32_t) * up;
			return xdr_putlong (xdrs, (int32_t *)&ui32);

		case XDR_DECODE:
			if (!xdr_getlong (xdrs, (int32_t *)&ui32))
				{
					return 0;
				}
			*up = (uint32_t) ui32;
		case XDR_FREE:
			return 1;
		}
	return 0;
}

static int
xdr_short (XDR *xdrs, short *sp)
{
	int32_t i32;

	switch (xdrs->x_op)
		{
		case XDR_ENCODE:
			i32 = (int32_t) *sp;
			return xdr_putlong (xdrs, &i32);

		case XDR_DECODE:
			if (!xdr_getlong (xdrs, &i32))
				{
					return 0;
				}
			*sp = (short) i32;
			return 1;

		case XDR_FREE:
			return 1;
		}
	return 0;
}

static int
xdr_u_short (XDR *xdrs, unsigned short *sp)
{
	uint32_t ui32;

	switch (xdrs->x_op)
		{
		case XDR_ENCODE:
			ui32 = (uint32_t) *sp;
			return xdr_putlong (xdrs, (int32_t *)&ui32);

		case XDR_DECODE:
			if (!xdr_getlong (xdrs, (int32_t *)&ui32))
				{
					return 0;
				}
			*sp = (unsigned short) ui32;
			return 1;

		case XDR_FREE:
			return 1;
		}
	return 0;
}

static int
xdr_char (XDR *xdrs, char *cp)
{
	int i;

	i = (*cp);
	if (!xdr_int (xdrs, &i))
		{
			return 0;
		}
	*cp = i;
	return 1;
}

static int
xdr_u_char (XDR *xdrs, unsigned char *cp)
{
	unsigned int u;

	u = (*cp);
	if (!xdr_u_int (xdrs, &u))
		{
			return 0;
		}
	*cp = u;
	return 1;
}

/*
 * XDR opaque data
 * Allows the specification of a fixed size sequence of opaque bytes.
 * cp points to the opaque object and cnt gives the byte length.
 */
static int
xdr_opaque (XDR *xdrs, char *cp, unsigned int cnt)
{
	unsigned int rndup;
	static char crud[BYTES_PER_XDR_UNIT];

	/*
	 * if no data we are done
	 */
	if (cnt == 0)
		return 1;

	/*
	 * round byte count to full xdr units
	 */
	rndup = cnt % BYTES_PER_XDR_UNIT;
	if (rndup > 0)
		rndup = BYTES_PER_XDR_UNIT - rndup;

	switch (xdrs->x_op)
		{
		case XDR_DECODE:
			if (!xdr_getbytes (xdrs, cp, cnt))
				{
					return 0;
				}
			if (rndup == 0)
				return 1;
			return xdr_getbytes (xdrs, (char *)crud, rndup);

		case XDR_ENCODE:
			if (!xdr_putbytes (xdrs, cp, cnt))
				{
					return 0;
				}
			if (rndup == 0)
				return 1;
			return xdr_putbytes (xdrs, xdr_zero, rndup);

		case XDR_FREE:
			return 1;
		}
#undef BYTES_PER_XDR_UNIT
	return 0;
}


/*
 * XDR null terminated ASCII strings
 */
static int
xdr_string (XDR *xdrs, char **cpp, unsigned int maxsize)
{
	char *sp = *cpp;	/* sp is the actual string pointer */
	unsigned int size = 0 ;
	unsigned int nodesize = 0;

	/*
	 * first deal with the length since xdr strings are counted-strings
	 */
	switch (xdrs->x_op)
		{
		case XDR_FREE:
			if (sp == NULL)
				{
					return 1;		/* already free */
				}
			/* fall through... */
		case XDR_ENCODE:
			if (sp == NULL)
				return 0;
			size = strlen (sp);
			break;
		case XDR_DECODE:
			break;
		}
	if (!xdr_u_int (xdrs, &size))
		{
			return 0;
		}
	if (size > maxsize)
		{
			return 0;
		}
	nodesize = size + 1;

	/*
	 * now deal with the actual bytes
	 */
	switch (xdrs->x_op)
		{
		case XDR_DECODE:
			if (nodesize == 0)
				{
					return 1;
				}
			if (sp == NULL)
				*cpp = sp = (char *) malloc (nodesize);
			if (sp == NULL)
				{
					(void) fputs ("xdr_string: out of memory\n", stderr);
					return 0;
				}
			sp[size] = 0;
			/* fall into ... */

		case XDR_ENCODE:
			return xdr_opaque (xdrs, sp, size);

		case XDR_FREE:
			sp=(char *)condfree (sp);
			*cpp = NULL;
			return 1;
		}
	return 0;
}



/* Floating-point stuff */

static int
xdr_float(XDR *xdrs, float *fp)
{
	switch (xdrs->x_op) {

	case XDR_ENCODE:
		if (sizeof(float) == sizeof(int32_t))
			return (xdr_putlong(xdrs, (int32_t *)fp));
		else if (sizeof(float) == sizeof(int)) {
			int32_t tmp = *(int *)fp;
			return (xdr_putlong(xdrs, &tmp));
		}
		break;

	case XDR_DECODE:
		if (sizeof(float) == sizeof(int32_t))
			return (xdr_getlong(xdrs, (int32_t *)fp));
		else if (sizeof(float) == sizeof(int)) {
			int32_t tmp;
			if (xdr_getlong(xdrs, &tmp)) {
				*(int *)fp = tmp;
				return (1);
			}
		}
		break;

	case XDR_FREE:
		return (1);
	}
	return (0);
}


static int
xdr_double(XDR *xdrs, double *dp)
{
    /* Gromacs detects floating-point stuff at compile time, which is faster */
#ifdef GROMACS
#  ifndef FLOAT_FORMAT_IEEE754 
#    error non-IEEE floating point system, or you defined GROMACS yourself...
#  endif
    int LSW;
#  ifdef IEEE754_BIG_ENDIAN_WORD_ORDER 
    int LSW=1;
#  else
    int LSW=0;
#  endif /* Big endian word order */
#else 
    /* Outside Gromacs we rely on dynamic detection of FP order. */
    int LSW; /* Least significant fp word */

    double x=0.987654321; /* Just a number */
    unsigned char ix = *((char *)&x);
    
    /* Possible representations in IEEE double precision: 
     * (S=small endian, B=big endian)
     *  
     * Byte order, Word order, Hex
     *     S           S       b8 56 0e 3c dd 9a ef 3f    
     *     B           S       3c 0e 56 b8 3f ef 9a dd
     *     S           B       dd 9a ef 3f b8 56 0e 3c
     *     B           B       3f ef 9a dd 3c 0e 56 b8
     */ 
    if(ix==0xdd || ix==0x3f)
		LSW=1;  /* Big endian word order */
    else if(ix==0xb8 || ix==0x3c)
		LSW=0;  /* Small endian word order */
    else { /* Catch strange errors */
		fprintf(stderr,"Cannot detect floating-point word order.\n"
				"Do you have a non-IEEE system?\n"
				"Use system XDR libraries or fix xdr_double().\n");
		abort();
    }
#endif /* end of dynamic detection of fp word order */

	switch (xdrs->x_op) {
    
	case XDR_ENCODE:
		if (2*sizeof(int32_t) == sizeof(double)) {
			int32_t *lp = (int32_t *)dp;
			return (xdr_putlong(xdrs, lp+!LSW) &&
					xdr_putlong(xdrs, lp+LSW));
		} else if (2*sizeof(int) == sizeof(double)) {
			int *ip = (int *)dp;
			int32_t tmp[2];
			tmp[0] = ip[!LSW];
			tmp[1] = ip[LSW];
			return (xdr_putlong(xdrs, tmp) &&
					xdr_putlong(xdrs, tmp+1));
		}
		break;
    
	case XDR_DECODE:
		if (2*sizeof(int32_t) == sizeof(double)) {
			int32_t *lp = (int32_t *)dp;
			return (xdr_getlong(xdrs, lp+!LSW) &&
					xdr_getlong(xdrs, lp+LSW));
		} else if (2*sizeof(int) == sizeof(double)) {
			int *ip = (int *)dp;
			int32_t tmp[2];
			if (xdr_getlong(xdrs, tmp+!LSW) &&
				xdr_getlong(xdrs, tmp+LSW)) {
				ip[0] = tmp[0];
				ip[1] = tmp[1];
				return (1);
			}
		}
		break;
    
	case XDR_FREE:
		return (1);
	}
	return (0);
}


static int xdrstdio_getlong (XDR *, int32_t *);
static int xdrstdio_putlong (XDR *, int32_t *);
static int xdrstdio_getbytes (XDR *, char *, unsigned int);
static int xdrstdio_putbytes (XDR *, char *, unsigned int);
static unsigned int xdrstdio_getpos (XDR *);
static int xdrstdio_setpos (XDR *, unsigned int);
static void xdrstdio_destroy (XDR *);

/*
 * Ops vector for stdio type XDR
 */
static const struct XDR::xdr_ops xdrstdio_ops =
	{
		xdrstdio_getlong,		/* deserialize a long int */
		xdrstdio_putlong,		/* serialize a long int */
		xdrstdio_getbytes,       	/* deserialize counted bytes */
		xdrstdio_putbytes,     	/* serialize counted bytes */
		xdrstdio_getpos,		/* get offset in the stream */
		xdrstdio_setpos,		/* set offset in the stream */
		xdrstdio_destroy,		/* destroy stream */
	};

/*
 * Initialize a stdio xdr stream.
 * Sets the xdr stream handle xdrs for use on the stream file.
 * Operation flag is set to op.
 */
static void
xdrstdio_create (XDR *xdrs, FILE *file, enum xdr_op op)
{
	xdrs->x_op = op;

	xdrs->x_ops = (struct XDR::xdr_ops *) &xdrstdio_ops;
	xdrs->x_private = (char *) file;
}

/*
 * Destroy a stdio xdr stream.
 * Cleans up the xdr stream handle xdrs previously set up by xdrstdio_create.
 */
static void
xdrstdio_destroy (XDR *xdrs)
{
	(void) fflush ((FILE *) xdrs->x_private);
	/* xx should we close the file ?? */
}

static int
xdrstdio_getlong (XDR *xdrs, int32_t *lp)
{
	int32_t mycopy;

	if (fread ((char *) & mycopy, 4, 1, (FILE *) xdrs->x_private) != 1)
		return 0;
	*lp = (int32_t) xdr_ntohl (mycopy);
	return 1;
}

static int
xdrstdio_putlong (XDR *xdrs, int32_t *lp)
{
	int32_t mycopy = xdr_htonl (*lp);
	lp = &mycopy;
	if (fwrite ((char *) lp, 4, 1, (FILE *) xdrs->x_private) != 1)
		return 0;
	return 1;
}

static int
xdrstdio_getbytes (XDR *xdrs, char *addr, unsigned int len)
{
	if ((len != 0) && (fread (addr, (int) len, 1,
							  (FILE *) xdrs->x_private) != 1))
		return 0;
	return 1;
}

static int
xdrstdio_putbytes (XDR *xdrs, char *addr, unsigned int len)
{
	if ((len != 0) && (fwrite (addr, (int) len, 1,
							   (FILE *) xdrs->x_private) != 1)) 
		return 0;
	return 1;
}

/* 32 bit fileseek operations */
static unsigned int
xdrstdio_getpos (XDR *xdrs)
{
	return (unsigned int) ftell ((FILE *) xdrs->x_private);
}

static int
xdrstdio_setpos (XDR *xdrs, unsigned int pos)
{
	return fseek ((FILE *) xdrs->x_private, pos, 0) < 0 ? 0 : 1;
}



#endif /* HAVE_RPC_XDR_H not defined */
