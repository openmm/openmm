
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

#ifndef __GromacsReferenceParameterPrint_H__
#define __GromacsReferenceParameterPrint_H__

#ifdef __cplusplus
#define externC extern "C"
#else
#define externC extern
#endif


/**---------------------------------------------------------------------------------------
      
   Print bond parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printBondParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* bondFile );
 
/**---------------------------------------------------------------------------------------
      
   Print angle parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printAngleParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* angleFile );

/**---------------------------------------------------------------------------------------
      
   Print rbDihedral parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */
externC int printRbDihedralParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* rbDihedralFile );

/**---------------------------------------------------------------------------------------
      
   Print properDihedral parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printProperDihedralParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* properDihedralFile );

/**---------------------------------------------------------------------------------------
      
   Print LJ 14 parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printLJ14Parameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* lj14File );

/**---------------------------------------------------------------------------------------
      
   Print LJ Coulomb parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printLJCoulombParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* ljCoulombFile );

/**---------------------------------------------------------------------------------------
      
   Print Shake parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC int printShakeParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* shakeFile );

/**---------------------------------------------------------------------------------------
      
   Print Gromacs parameters 
   
   @param baseFileName        base file name
   @param top                 Gromacs parameter data struct
   @param fr                  Gromacs parameter data struct
   @param md                  Gromacs parameter data struct
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC void printReferenceParameters( char* baseFileName, t_topology *top, t_forcerec *fr,
                                           t_mdatoms *md, FILE *log );

/**---------------------------------------------------------------------------------------
      
   Print forces/energies
   
   @param numberOfAtoms       number of atoms
   @param x                   atom coordinates
   @param f                   atom forces
   @param energies            energies by atom
   @param totalEnergy         total energy
   @param fileName            file name

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC void writeGromacsForceFile(  int numberOfAtoms, rvec x[], rvec f[], real* energies,
                                         real totalEnergy, char* fileName );

/**---------------------------------------------------------------------------------------
      
   Print 'state' mass/coordinates/velocities/forces
   
   @param numberOfAtoms       number of atoms
   @param x                   atom coordinates
   @param v                   atom velocities
   @param f                   atom forces
   @param timestep            time step
   @param deltaT              delta t
   @param temperature         temperature
   @param tau                 tau
   @param fileName            file name

   @return 0
      
   --------------------------------------------------------------------------------------- */

externC void writeGromacsStateFile( int numberOfAtoms, real mass[], rvec x[], rvec v[],
                                       rvec f[], int timestep, real deltaT,
                                       real temperature, real tau, char* fileName );


#endif
