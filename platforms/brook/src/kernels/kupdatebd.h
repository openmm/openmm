/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs, Chris Bruns                       *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/*
 * Brownian dynamics integration
 *
 * @param xstrwidth      atom stream width
 * @param gstrwidth      Gaussian stream width
 * @param goffset        Gaussian offset into stream
 * @param forceScale     force scale factor
 * @param noiseAmplitude noise amplitude
 * @param fgauss         random numbers
 * @param pos            atom positions
 * @param force          force
 * @param posp           delta positions
 *
 **/
 
void kintegrate_bd( const float xstrwidth, const float gstrwidth, const float goffset,
                    const float forceScale, const float noiseAmplitude, 
                    ::brook::stream fgauss, ::brook::stream posq, 
                    ::brook::stream force, ::brook::stream posp );

/*
 * Brownian dynamics update
 *
 * @param velocityScale  velocity scale
 * @param posp           atom positions
 * @param posIn          atom positions
 * @param velocity       velocity
 * @param posp           delta positions
 *
 **/

void kupdate_bd( const float velocityScale, ::brook::stream posp, ::brook::stream posIn, 
                 ::brook::stream velocity, ::brook::stream posOut );


/*
 * Brownian dynamics update
 *
 * @param velocityScale  velocity scale
 * @param posp           atom positions
 * @param velocity       velocity
 * @param posp           delta positions
 *
 **/

void kupdate_bd2( const float velocityScale, ::brook::stream posp,
                  ::brook::stream posIn, ::brook::stream velocity, ::brook::stream posOut );


