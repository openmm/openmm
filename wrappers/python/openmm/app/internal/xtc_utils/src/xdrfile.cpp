/* This file was adapted from GROMACS for OpenMM.  Only parts necessary for
 * reading and writing XTC files are retained.  xdr_int(), xdr_float(), and
 * xdr_opaque() have been re-implemented since the original implementations
 * were adapted from code whose license restricted redistribution.
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

#include "xdrfile.h"

#include <climits>
#include <cmath>
#include <cstdlib>

static_assert(sizeof(int) == 4 && sizeof(float) == 4, "unsupported type sizes for XDR");

// Swaps 4 bytes to/from big-endian order if native order is little-endian.
static void swap_big_endian(char* bytes) {
	int word = 1;
	if (*reinterpret_cast<char*>(&word)) {
		// Native order is little-endian.
		char b0 = bytes[0];
		char b1 = bytes[1];
		char b2 = bytes[2];
		char b3 = bytes[3];
		bytes[0] = b3;
		bytes[1] = b2;
		bytes[2] = b1;
		bytes[3] = b0;
	}
}

static bool xdr_int(FILE* file, xdr_op op, int* ptr) {
	int value;
	switch (op) {
		case XDR_ENCODE:
			value = *ptr;
			swap_big_endian(reinterpret_cast<char*>(&value));
			return fwrite(&value, 4, 1, file) == 1;
		case XDR_DECODE:
			if (fread(ptr, 4, 1, file) != 1) {
				return false;
			}
			swap_big_endian(reinterpret_cast<char*>(ptr));
			return true;
		default:
			return false;
	}
}

static bool xdr_float(FILE* file, xdr_op op, float* ptr) {
	float value;
	switch (op) {
		case XDR_ENCODE:
			value = *ptr;
			swap_big_endian(reinterpret_cast<char*>(&value));
			return fwrite(&value, 4, 1, file) == 1;
		case XDR_DECODE:
			if (fread(ptr, 4, 1, file) != 1) {
				return false;
			}
			swap_big_endian(reinterpret_cast<char*>(ptr));
			return true;
		default:
			return false;
	}
}

static bool xdr_opaque(FILE* file, xdr_op op, char* ptr, int n) {
	char pad[4] = {0, 0, 0, 0};
	int nPad = (n + 3) / 4 * 4 - n;
	switch(op) {
		case XDR_ENCODE:
			return fwrite(ptr, 1, n, file) == n && (!nPad || fwrite(pad, 1, nPad, file) == nPad);
		case XDR_DECODE:
			return fread(ptr, 1, n, file) == n && (!nPad || fread(pad, 1, nPad, file) == nPad);
		default:
			return false;
	}
}

static void * condfree( void*p) {
        if(p) { free(p); }
	else {
		fprintf( stderr, "Double free detected\n" );
		abort();
	}
        return NULL;
}

const char *exdr_message[exdrNR] = {
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
	xfp->op = xdrmode;
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
	while(i<ndata && xdr_int(xfp->fp,xfp->op,ptr+i))
		i++;

	return i;
}

int
xdrfile_write_int(int *ptr, int ndata, XDRFILE* xfp)
{
	int i=0;

	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_int(xfp->fp,xfp->op,ptr+i))
		i++;
	return i;
}

int
xdrfile_read_float(float *ptr, int ndata, XDRFILE* xfp)
{
	int i=0;
	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_float(xfp->fp,xfp->op,ptr+i))
		i++;
	return i;
}

int
xdrfile_write_float(float *ptr, int ndata, XDRFILE* xfp)
{
	int i=0;
	/* read write is encoded in the XDR struct */
	while(i<ndata && xdr_float(xfp->fp,xfp->op,ptr+i))
		i++;
	return i;
}

int
xdrfile_read_opaque(char *ptr, int cnt, XDRFILE* xfp)
{
	if(xdr_opaque(xfp->fp,xfp->op,ptr,cnt))
		return cnt;
	else
		return 0;
}


int
xdrfile_write_opaque(char *ptr, int cnt, XDRFILE* xfp)
{
	if(xdr_opaque(xfp->fp,xfp->op,ptr,cnt))
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
