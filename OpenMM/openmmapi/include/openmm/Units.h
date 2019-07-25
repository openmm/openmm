#ifndef OPENMM_UNITS_H_
#define OPENMM_UNITS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include <cmath>

/** \file
 * OpenMM uses the following units everywhere:
 *
 * length: nanometers<br/>
 * time: picoseconds<br/>
 * mass: atomic mass units (daltons)<br/>
 * charge: proton charge<br/>
 * temperature: Kelvin<br/>
 * angle: radians<br/>
 * energy: kJ/mol<br/>
 * force: kJ/mol/nm<br/>
 *
 * Because some programs use other units (e.g. Angstroms for length or kcal/mol for energy), this file defines
 * constants which can be used to convert to and from OpenMM's units.
 */

namespace OpenMM {
    /**
     * The number of nanometers in an Angstrom.
     */
    static const double NmPerAngstrom = 0.1;
    /**
     * The number of Angstroms in a nanometer.
     */
    static const double AngstromsPerNm = 10.0;
    /**
     * The number of picoseconds in a femtosecond.
     */
    static const double PsPerFs = 0.001;
    /**
     * The number of femtoseconds in a picosecond.
     */
    static const double FsPerPs = 1000.0;
    /**
     * The number of kJ in a kcal.
     */
    static const double KJPerKcal = 4.184;
    /**
     * The number of kcal in a kJ.
     */
    static const double KcalPerKJ = 1.0/4.184;
    /**
     * The number of radians in a degree.
     */
    static const double RadiansPerDegree = 3.1415926535897932385/180.0;
    /**
     * The number of degrees in a radian.
     */
    static const double DegreesPerRadian = 180.0/3.1415926535897932385;
    /**
     * L-J sigma per unit van der Waals radius: 2/(2^1/6).
     */
    static const double SigmaPerVdwRadius = 1.7817974362806786095;
    /**
     * van der Waals radius per unit L-J sigma: (2^1/6)/2.
     */
    static const double VdwRadiusPerSigma = .56123102415468649070;

} // namespace OpenMM

#endif /*OPENMM_UNITS_H_*/
