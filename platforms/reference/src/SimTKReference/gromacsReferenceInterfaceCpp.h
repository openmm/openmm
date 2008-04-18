
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef __GromacsReferenceInterfaceCpp_H__
#define __GromacsReferenceInterfaceCpp_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif

/**---------------------------------------------------------------------------------------

   Allocate memory for atom indices participating in bonds and the bond parameters

   @param numberOfBonds            number of bonds
   @param numberOfAtomIndices      (number of atoms)/bond
   @param outputBondIndices        allocated memory for atom indices outputBondIndices[bondIndex][atomIndex]
   @param numberParameters         number of parameters/bond
   @param outputBondParameters     allocated memory for parameters outputBondParameters[bondIndex][parameterIndex]

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsAllocateBondIxnArrays( int numberOfBonds, int numberOfAtomIndices,
                                  int*** outputBondIndices, int numberParameters,
                                   RealOpenMM*** outputBondParameters );

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for harmonic bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

externC
int gromacsGetHarmonicBondParameters( const t_topology* top, int*** atomIndices,
                                      RealOpenMM*** bondParameters );

/**---------------------------------------------------------------------------------------

   Calculate harmonic bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param harmonicBondForces       output array of forces: harmonicBondForces[atomIndex][3]
   @param harmonicBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsHarmonicBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                       RealOpenMM** harmonicBondForces,
                                       RealOpenMM* harmonicBondAtomEnergies,
                                       RealOpenMM* totalEnergy );

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for angle bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

externC
int gromacsGetAngleBondParameters( const t_topology* top, int*** atomIndices,
                                   RealOpenMM*** bondParameters );

/**---------------------------------------------------------------------------------------

   Calculate angle bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param angleBondForces       output array of forces: angleBondForces[atomIndex][3]
   @param angleBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsAngleBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                    RealOpenMM** angleBondForces,
                                    RealOpenMM* angleBondAtomEnergies,
                                    RealOpenMM* totalEnergy ); 

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for properDihedral bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

externC
int gromacsGetProperDihedralBondParameters( const t_topology* top, int*** atomIndices,
                                            RealOpenMM*** bondParameters );

/**---------------------------------------------------------------------------------------

   Calculate properDihedral bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param properDihedralBondForces       output array of forces: angleBondForces[atomIndex][3]
   @param properDihedralBondEnergies     output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsProperDihedralBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                             RealOpenMM** properDihedralBondForces,
                                             RealOpenMM* properDihedralBondAtomEnergies,
                                             RealOpenMM* totalEnergy );

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for RB Dihedral bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
		free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

externC
int gromacsGetRbDihedralBondParameters( const t_topology* top, int*** atomIndices,
                                            RealOpenMM*** bondParameters );

/**---------------------------------------------------------------------------------------

   Calculate RB dihedral bond forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param rbDihedralBondForces     output array of forces: angleBondForces[atomIndex][3]
   @param rbDihedralBondEnergies   output array of bond energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsRbDihedralBondReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                                         RealOpenMM** rbDihedralBondForces,
                                         RealOpenMM* rbDihedralBondAtomEnergies,
                                         RealOpenMM* totalEnergy );

/**---------------------------------------------------------------------------------------

   Get bond atom indices and parameters for LJ 14 bonds

   @param top                      Gromacs t_topology data struct
   @param atomIndices              array of atomIndices[bondIndex][atomIndex]
   @param bondParameters           array of bond parameters[bondIndex][parameterIndex]

   The arrays atomIndices & bondParameters are allocated off the heap;
   The memory for the second index is one block of size (number of bonds)*(number of atoms in bond)

   The memory should be freed by the callee via calls:
      free( atomIndices[0] ); free( atomIndices );

   @return number of bonds

   --------------------------------------------------------------------------------------- */

externC
int gromacsGetLJ14Parameters( const t_topology* top, int*** atomIndices,
                              RealOpenMM*** bondParameters );

/**---------------------------------------------------------------------------------------

   Calculate LJ 14 forces/energies for given configuration and topology

   @param top                      Gromacs t_topology data struct
   @param atomCoordinates          atom configuration
   @param fr                       Gromacs t_forcerec struct
   @param lj14Forces               output array of forces: lj14Forces[atomIndex][3]
   @param lj14AtomEnergies         output array of LJ 14 energies[atomIndex]
   @param totalEnergy              output total energy

   @return 0

   --------------------------------------------------------------------------------------- */

externC
int gromacsLJ14ReferenceForce( const t_topology* top, RealOpenMM** atomCoordinates,
                               t_forcerec *fr, const RealOpenMM* charges, RealOpenMM** lj14Forces,
                               RealOpenMM* lj14AtomEnergies,
                               RealOpenMM* totalEnergy );

#endif
