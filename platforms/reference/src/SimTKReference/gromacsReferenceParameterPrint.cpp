
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

#include "../SimTKUtilities/SimTKOpenMMGromacsUtilities.h"
#include "gromacsReferenceParameterPrint.h"

#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"

/**---------------------------------------------------------------------------------------
      
   Print bond parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printBondParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* bondFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintBondParameters";

   // ---------------------------------------------------------------------------------------

   t_iatom *atoms;
   int  ii, numberOfBonds, type, i, j, bondIndex;
  
   t_idef *idef   = &top->idef;
   atoms          = idef->il[ F_BONDS ].iatoms;
   numberOfBonds  = idef->il[ F_BONDS ].nr/3; 
   
  (void) fprintf( bondFile, "%d Bond parameters: [atom indices] idealBond k\n", numberOfBonds );

   bondIndex      = 1;
   for( ii = 0; ii < idef->il[ F_BONDS ].nr; ii += 3 ) {

      type = atoms[ ii     ]; 
      i    = atoms[ ii + 1 ]; 
      j    = atoms[ ii + 2 ]; 

      (void) fprintf( bondFile, "%6d %6d %6d %13.6e %13.6e\n",
                                 bondIndex++, i, j, idef->iparams[ type ].harmonic.rA, idef->iparams[ type ].harmonic.krA );

   }
   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print angle parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printAngleParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* angleFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintAngleParameters";

   // ---------------------------------------------------------------------------------------

   t_iatom *atoms;
   int  ii, numberOfAngles, type, i, j, k, angleIndex;
  
   t_idef *idef   = &top->idef;
   atoms          = idef->il[ F_ANGLES ].iatoms;
   numberOfAngles = idef->il[ F_ANGLES ].nr/4; 
   
  (void) fprintf( angleFile, "%d Angle parameters: [atom indices] idealAngle k\n", numberOfAngles );

   angleIndex     = 1;
   for( ii = 0; ii < idef->il[ F_ANGLES ].nr; ii += 4 ) {

      type = atoms[ ii     ]; 
      i    = atoms[ ii + 1 ]; 
      j    = atoms[ ii + 2 ]; 
      k    = atoms[ ii + 3 ]; 

      (void) fprintf( angleFile, "%6d %6d %6d %6d %13.6e %13.6e\n",
                                 angleIndex++, i, j, k, idef->iparams[ type ].harmonic.rA, idef->iparams[ type ].harmonic.krA );

   }
   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print rbDihedral parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printRbDihedralParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* rbDihedralFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintRbDihedralParameters";

   // ---------------------------------------------------------------------------------------

   t_iatom *atoms;
   int  ii, numberOfRbDihedrals, type, i, j, k, l, rbDihedralIndex;
  
   t_idef *idef        = &top->idef;
   atoms               = idef->il[ F_RBDIHS ].iatoms;
   numberOfRbDihedrals = idef->il[ F_RBDIHS ].nr/5; 
   
  (void) fprintf( rbDihedralFile, "%d RbDihedral parameters: [atom indices] RbDihedral parameters\n", numberOfRbDihedrals );

   rbDihedralIndex     = 1;
   for( ii = 0; ii < idef->il[ F_RBDIHS ].nr; ii += 5 ) {

      type = atoms[ ii     ]; 
      i    = atoms[ ii + 1 ]; 
      j    = atoms[ ii + 2 ]; 
      k    = atoms[ ii + 3 ]; 
      l    = atoms[ ii + 4 ]; 

      (void) fprintf( rbDihedralFile, "%6d %6d %6d %6d %6d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                      rbDihedralIndex++, i, j, k, l, 
                      idef->iparams[ type ].rbdihs.rbc[0],
                      idef->iparams[ type ].rbdihs.rbc[1],
                      idef->iparams[ type ].rbdihs.rbc[2],
                      idef->iparams[ type ].rbdihs.rbc[3],
                      idef->iparams[ type ].rbdihs.rbc[4],
                      idef->iparams[ type ].rbdihs.rbc[5] );

   }
   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print properDihedral parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printProperDihedralParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* properDihedralFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintProperDihedralParameters";

   // ---------------------------------------------------------------------------------------

   t_iatom *atoms;
   int  ii, numberOfProperDihedrals, type, i, j, k, l, properDihedralIndex;
  
   t_idef *idef            = &top->idef;
   atoms                   = idef->il[ F_PDIHS ].iatoms;
   numberOfProperDihedrals = idef->il[ F_PDIHS ].nr/5; 
   
  (void) fprintf( properDihedralFile, "%d ProperDihedral parameters: [atom indices] ProperDihedral parameters\n", numberOfProperDihedrals );

   properDihedralIndex     = 1;
   for( ii = 0; ii < idef->il[ F_PDIHS ].nr; ii += 5 ) {

      type = atoms[ ii     ]; 
      i    = atoms[ ii + 1 ]; 
      j    = atoms[ ii + 2 ]; 
      k    = atoms[ ii + 3 ]; 
      l    = atoms[ ii + 4 ]; 

      (void) fprintf( properDihedralFile, "%6d %6d %6d %6d %6d %13.6e %13.6e %3d\n",
                      properDihedralIndex++, i, j, k, l, 
                      idef->iparams[ type ].pdihs.cpA,
                      idef->iparams[ type ].pdihs.phiA,
                      idef->iparams[ type ].pdihs.mult );

   }
   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print LJ 14 parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printLJ14Parameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* lj14File ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintLj14Parameters";

   real c6, c12;
   real oneSixth = ((real) 1.0)/ ( (real) 6.0);

   real zero     = 0.0;

   // ---------------------------------------------------------------------------------------

   t_iatom *atoms;
   int  ii, numberOfLJ14, type, i, j, lj14Index;
  
   t_idef *idef            = &top->idef;
   real* charges           = md->chargeA;

   atoms                   = idef->il[ F_LJ14 ].iatoms;
   numberOfLJ14            = idef->il[ F_LJ14 ].nr/3; 
   
  (void) fprintf( lj14File, "%d LJ 14 parameters: [atom indices] c6 c12 q1 q2 where epsfac=%.5e fudge=%.5e\n", numberOfLJ14,
                  fr->epsfac, fr->fudgeQQ );

   lj14Index     = 1;
   for( ii = 0; ii < idef->il[ F_LJ14 ].nr; ii += 3 ) {

      type = atoms[ ii     ]; 
      i    = atoms[ ii + 1 ]; 
      j    = atoms[ ii + 2 ]; 

      c6   = idef->iparams[ type ].lj14.c6A;
      c12  = idef->iparams[ type ].lj14.c12A;


/*
      if( c12 > zero && c6 > zero ){
         sig  =  c6*c6/c12;
         eps  =  (real) pow( (c12/c6), oneSixth );
      } else {
         sig  = zero;
         eps  = zero;
      }

      coulomb = fr->epsfac*fr->fudgeQQ*charges[i]*charges[j];
*/

      (void) fprintf( lj14File, "%6d %6d %6d %14.6e %14.6e %9.5f %9.5f\n",
                      lj14Index++, i, j, c6, c12, charges[i], charges[j] );

   }

   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print LJ Coulomb parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printLJCoulombParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* ljCoulombFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintLjCoulombParameters";

   real c6, c12;

   real zero      = 0.0;
   real half      = ((real) 1.0)/ ( (real) 2.0);
   real oneSixth  = ((real) 1.0)/ ( (real) 6.0);
   real onetwelth = ((real) 1.0)/ ( (real) 12.0);

   real mass;
   real scaleFactor;

   int errors    = 0;

   // ---------------------------------------------------------------------------------------

   int  ii, jj;
  
   t_idef *idef            = &top->idef;
   real* charges           = md->chargeA;
   t_block *gmx_excl       = &top->atoms.excl;
   int ntypes              = fr->ntype;
   int* types              = md->typeA;
   real* nbfp              = fr->nbfp;
   
  (void) fprintf( ljCoulombFile, "%d epsfac=%.5e fudge=%.5e  LJ Coulomb parameters: { atom c6 c12 q numberOfExclusions Exclsions[] }\n", md->nr,
                  fr->epsfac, fr->fudgeQQ );

   // loop over atoms

   for( ii = 0; ii < md->nr; ii++ ) {

      c6              = nbfp[ types[ii]*2*ntypes + types[ii]*2 ];
      c12             = nbfp[ types[ii]*2*ntypes + types[ii]*2 + 1 ];

      char* atomTypeC = *(top->atoms.atomtype[ii]);

/*
      if( c12 > zero && c6 > zero ){
         sig     =  half*c6*c6/c12;
         eps     =  (real) pow( (c12/c6), onetwelth );
      } else {
         sig     = zero;
         eps     = zero;
      }
*/

      mass = md->massA[ii];
      if ( mass < 1.2f && mass >= 1.0f ){        // hydrogen
         scaleFactor  = 0.85f; 
      } else if( mass > 11.8f && mass < 12.2f ){ // carbon
         scaleFactor  = 0.72f; 
      } else if( mass > 14.0f && mass < 15.0f ){ // nitrogen
         scaleFactor  = 0.79f;
      } else if( mass > 15.5f && mass < 16.5f ){ // oxygen
         scaleFactor  = 0.85f; 
      } else if( mass > 31.5f && mass < 32.5f ){ // sulphur
         scaleFactor  = 0.96f;
      } else if( mass > 29.5f && mass < 30.5f ){ // phosphorus
         scaleFactor  = 0.86f;
      } else {
         errors++;
         (void) fprintf( ljCoulombFile, "Warning: mass for atom %d mass=%.4f not recognized.",
                         ii, mass );
      }   

      (void) fprintf( ljCoulombFile, "%6d %14.6e %14.6e %9.5f %12s %9.4f %d",
                      ii, c6, c12, charges[ii], atomTypeC, scaleFactor, (gmx_excl->index[ii+1] - gmx_excl->index[ii]) );

      for ( jj = gmx_excl->index[ii]; jj < gmx_excl->index[ii+1]; jj++ ) {
         (void) fprintf( ljCoulombFile, " %d", gmx_excl->a[jj] );
      }    
      (void) fprintf( ljCoulombFile, "\n" );

      // symmetry check

      for ( jj = gmx_excl->index[ii]; jj < gmx_excl->index[ii+1]; jj++ ){
         int checkI = gmx_excl->a[jj];
         int hit    = 0;
         for( int kk = gmx_excl->index[checkI]; kk < gmx_excl->index[checkI+1] && !hit; kk++ ){
            if( gmx_excl->a[kk] == ii ){
               hit = 1;
            }
         }
         if( !hit ){
            (void) fprintf( ljCoulombFile, " Missing %d in list for %d\n", ii, checkI );
         }
      }    

   }

   if( errors ){
      (void) fprintf( stdout, "printLJ14Parameters errors -- aborting!\n" );
      (void) fprintf( stderr, "printLJ14Parameters errors -- aborting!\n" );
   }

   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print Shake parameters 
   
   @param idef                Gromacs parameter data struct
   @param fr                  force record?
   @param charges             atom charges
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" int printShakeParameters( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* shakeFile ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nprintShakeParameters";

   // ---------------------------------------------------------------------------------------

   t_idef* idef                = &top->idef;
   t_iatom* atoms              = idef->il[ F_SHAKE ].iatoms;
   int numberOfConstraints     = idef->il[ F_SHAKE ].nr/3;

  (void) fprintf( shakeFile, "# %d atomI atomJ d0 invMassI invMassJ\n", numberOfConstraints );

   // loop over constraints 

   int offset = 0;
   for( int ii = 0; ii < numberOfConstraints; ii++ ) {

      int type       = atoms[ offset     ]; 
      int atomI      = atoms[ offset + 1 ]; 
      int atomJ      = atoms[ offset + 2 ]; 
      offset        += 3;

      // atom indicies, distance constraint, inverse mass

      (void) fprintf( shakeFile, "%6d %6d %6d %14.6e %14.6e %14.6e\n", ii, atomI, atomJ, idef->iparams[type].shake.dA, md->invmass[atomI], md->invmass[atomJ] );
   }

   return 0;
}

/**---------------------------------------------------------------------------------------
      
   Print Gromacs parameters 
   
   @param baseFileName        base file name
   @param top                 Gromacs parameter data struct
   @param fr                  Gromacs parameter data struct
   @param md                  Gromacs parameter data struct
   @param file                FILE* reference

   @return 0
      
   --------------------------------------------------------------------------------------- */

extern "C" void printReferenceParameters( char* baseFileName, t_topology *top, t_forcerec *fr,
                                          t_mdatoms *md, FILE *log ){

   // ---------------------------------------------------------------------------------------

   typedef int printParameterFunction( t_topology *top, t_forcerec *fr, t_mdatoms *md, FILE* file );

   typedef struct {
      char* parameterFileName;
      printParameterFunction* printP;
   } printParameters;

   static const int ParameterCount = 7;
   static printParameters printParameterList[7] = {
      { "HarmonicBondParameter",    printBondParameters             },
      { "AngleBondParameter",       printAngleParameters            },
      { "RbDihedralParameter",      printRbDihedralParameters       },
      { "ProperDihedralParameter",  printProperDihedralParameters   },
      { "LJ14Parameter",            printLJ14Parameters             },
      { "ShakeParameters",          printShakeParameters            },
      { "LJCoulombParameter",       printLJCoulombParameters        }
                                                  };

   // ---------------------------------------------------------------------------------------

   int ii;
   static const int bufferSz = 1024;
   char  parameterFileName[1024];
   FILE* parameterFile;

   // ---------------------------------------------------------------------------------------

   // loop over parameter types 

   for( ii = 0; ii < ParameterCount; ii++ ){

      // open file

#ifdef   WIN32
      (void) sprintf_s( parameterFileName, bufferSz, "%s%s.txt", baseFileName, printParameterList[ii].parameterFileName );
      fopen_s( &parameterFile, parameterFileName, "w" );
#else
      (void) sprintf( parameterFileName, "%s%s.txt", baseFileName, printParameterList[ii].parameterFileName );
      parameterFile = fopen( parameterFileName, "w" );
#endif

      if( parameterFile == NULL ){
         (void) fprintf( log, "Parameter file=<%s> could not be opened.\n", parameterFileName );
      } else {

         // do it

         (*printParameterList[ii].printP)( top, fr, md, parameterFile );
         (void) fclose( parameterFile );
      }
   }

   return;
}

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

extern "C" void writeGromacsForceFile(  int numberOfAtoms, rvec x[], rvec f[], real* energies,
                                         real totalEnergy, char* fileName ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\ngromacsReferenceParameterPrint::writeGromacsForceFile";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* forceFile;
#ifdef   WIN32
   fopen_s( &forceFile, fileName, "w" );
#else
   forceFile = fopen( fileName, "w" );
#endif

   if( forceFile == NULL ){
      std::stringstream message;
      message << methodName << " force file=<" << fileName << " could not be opened.\n";
      SimTKOpenMMLog::printMessage( message );
   } else {

      // do it

      (void) fprintf( forceFile, "# %6d E=%14.6e x[i][3] f[i][3] E[i]\n", numberOfAtoms, totalEnergy );
      for( int ii = 0; ii < numberOfAtoms; ii++ ){
         (void) fprintf( forceFile, "%5d %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ii,
                         x[ii][0], x[ii][1], x[ii][2],
                         f[ii][0], f[ii][1], f[ii][2], (energies ? energies[ii] : 0.0) );
      }
      (void) fclose( forceFile );
   }

   return;
}

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

extern "C" void writeGromacsStateFile( int numberOfAtoms, real mass[], rvec x[], rvec v[],
                                       rvec f[], int timestep, real deltaT,
                                       real temperature, real tau, char* fileName ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\ngromacsReferenceParameterPrint::writeGromacsStateFile";

   // ---------------------------------------------------------------------------------------

   // open file

   FILE* stateFile;
#ifdef   WIN32
   fopen_s( &stateFile, fileName, "w" );
#else
   stateFile = fopen( fileName, "w" );
#endif

   if( stateFile == NULL ){

      std::stringstream message;
      message << methodName << " state file=<" << fileName << " could not be opened.\n";
      SimTKOpenMMLog::printMessage( message );

   } else {

      // do it

      (void) fprintf( stateFile, "# mass coordinates velocities forces\n", numberOfAtoms );
      (void) fprintf( stateFile, "%6d         # Atoms\n", numberOfAtoms );
      (void) fprintf( stateFile, "%6d         # timestep\n", timestep);
      (void) fprintf( stateFile, "%14.6e # delta_t\n", deltaT );
      (void) fprintf( stateFile, "%14.6e # temperature\n", temperature );
      (void) fprintf( stateFile, "%14.6e # tau\n", tau );
      for( int ii = 0; ii < numberOfAtoms; ii++ ){

         (void) fprintf( stateFile, "%6d ", ii );

         if( mass ){
            (void) fprintf( stateFile, "%14.6e   ", mass[ii] );
         }
         if( x ){
            (void) fprintf( stateFile, "%14.6e %14.6e %14.6e    ", x[ii][0], x[ii][1], x[ii][2] );
         }
         if( v ){
            (void) fprintf( stateFile, "%14.6e %14.6e %14.6e    ", v[ii][0], v[ii][1], v[ii][2] );
         }
         if( f ){
            (void) fprintf( stateFile, "%14.6e %14.6e %14.6e    ", f[ii][0], f[ii][1], f[ii][2] );
         }
         (void) fprintf( stateFile, "\n" );
      }
      (void) fclose( stateFile );

   }

   return;
}

