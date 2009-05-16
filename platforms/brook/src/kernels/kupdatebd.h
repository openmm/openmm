/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Mark Friedrichs, Mike Houston                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
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


