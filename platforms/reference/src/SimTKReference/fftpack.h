/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#include "../SimTKUtilities/SimTKOpenMMCommon.h"

#ifndef _FFTPACK_H_
#define _FFTPACK_H_


#include <stdio.h>



#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


typedef struct {
	RealOpenMM re;
	RealOpenMM im;
} t_complex;


/*! \brief Datatype for FFT setup 
 *
 *  The fftpack_t type contains all the setup information, e.g. twiddle
 *  factors, necessary to perform an FFT. Internally it is mapped to 
 *  whatever FFT library we are using, or the built-in FFTPACK if no fast
 *  external library is available.
 */
typedef struct fftpack *
fftpack_t;




/*! \brief Specifier for FFT direction. 
 *
 *  The definition of the 1D forward transform from input x[] to output y[] is
 *  \f[
 *  y_{k} = \sum_{j=0}^{N-1} x_{j} \exp{-i 2 \pi j k /N}
 *  \f]
 *
 *  while the corresponding backward transform is
 *
 *  \f[
 *  y_{k} = \sum_{j=0}^{N-1} x_{j} \exp{i 2 \pi j k /N}
 *  \f]
 *
 *  A forward-backward transform pair will this result in data scaled by N.
 *
 */
typedef enum fftpack_direction
{
    FFTPACK_FORWARD,         /*!< Forward complex-to-complex transform  */
    FFTPACK_BACKWARD,        /*!< Backward complex-to-complex transform */
} fftpack_direction;



/*! \brief Setup a 1-dimensional complex-to-complex transform 
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform 
 *
 *  \return status - 0 or a standard error message.
 */
int
fftpack_init_1d        (fftpack_t *       fft,
                        int               nx);



/*! \brief Setup a 2-dimensional complex-to-complex transform 
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform in first dimension
 *  \param ny     Length of transform in second dimension
 *
 *  \return status - 0 or a standard error message.
 *   
 */
int
fftpack_init_2d        (fftpack_t *         fft,
                        int                 nx, 
                        int                 ny);




/*! \brief Setup a 3-dimensional complex-to-complex transform 
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform in first dimension
 *  \param ny     Length of transform in second dimension
 *  \param nz     Length of transform in third dimension
 *
 *  \return status - 0 or a standard error message.
 *   
 */
int
fftpack_init_3d        (fftpack_t *         fft,
                        int                 nx, 
                        int                 ny,
                        int                 nz);




/*! \brief Perform a 1-dimensional complex-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from fftpack_init_1d()
 *  \param dir       Forward or Backward
 *  \param in_data   Input grid data. 
 *  \param out_data  Output grid data. 
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers 
 *       depending on your grid type.
 */
int 
fftpack_exec_1d          (fftpack_t                  setup,
                          enum fftpack_direction     dir,
                          t_complex *                in_data,
                          t_complex *                out_data);

/*! \brief Perform a 2-dimensional complex-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from fftpack_init_1d()
 *  \param dir       Forward or Backward
 *  \param in_data   Input grid data. 
 *  \param out_data  Output grid data. 
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers 
 *       depending on your grid type.
 */
int 
fftpack_exec_2d          (fftpack_t                  setup,
                          enum fftpack_direction     dir,
                          t_complex *                in_data,
                          t_complex *                out_data);


/*! \brief Perform a 3-dimensional complex-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from fftpack_init_1d()
 *  \param dir       Forward or Backward
 *  \param in_data   Input grid data. 
 *  \param out_data  Output grid data. 
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers 
 *       depending on your grid type.
 */
int 
fftpack_exec_3d          (fftpack_t                  setup,
                          enum fftpack_direction     dir,
                          t_complex *                in_data,
                          t_complex *                out_data);


/*! \brief Release an FFT setup structure 
 *
 *  Destroy setup and release all allocated memory.
 *
 *  \param setup Setup returned from fftpack_init_1d(), or one
 *		 of the other initializers.
 *
 */
void
fftpack_destroy          (fftpack_t                 setup);



#ifdef __cplusplus
}
#endif

#endif /* _FFTPACK_H_ */
