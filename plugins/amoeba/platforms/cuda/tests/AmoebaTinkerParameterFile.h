/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "../../../tests/AssertionUtilities.h"

#include "openmm/CMMotionRemover.h"
#include "openmm/System.h"
#include "openmm/Context.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/VariableLangevinIntegrator.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VariableVerletIntegrator.h"
#include "openmm/BrownianIntegrator.h"

#include "AmoebaHarmonicBondForce.h"
#include "AmoebaHarmonicAngleForce.h"
#include "AmoebaHarmonicInPlaneAngleForce.h"
#include "AmoebaTorsionForce.h"
#include "AmoebaPiTorsionForce.h"
#include "AmoebaStretchBendForce.h"
#include "AmoebaOutOfPlaneBendForce.h"
#include "AmoebaTorsionTorsionForce.h"
#include "AmoebaMultipoleForce.h"
#include "AmoebaGeneralizedKirkwoodForce.h"
#include "AmoebaVdwForce.h"
#include "AmoebaWcaDispersionForce.h"
#include "AmoebaSASAForce.h"
#include "internal/windowsExport.h"

#include <ctime>
#include <vector>
#include <algorithm>
#include <map>
#include <cfloat>
#include <cstring>
#include <cstdlib>
#include <typeinfo>
#include <time.h>

// force enums

#define MAX_PRINT 5

static std::string AMOEBA_HARMONIC_BOND_FORCE                         = "AmoebaHarmonicBond";
static std::string AMOEBA_HARMONIC_ANGLE_FORCE                        = "AmoebaHarmonicAngle"; 
static std::string AMOEBA_HARMONIC_IN_PLANE_ANGLE_FORCE               = "AmoebaHarmonicInPlaneAngle"; 
static std::string AMOEBA_TORSION_FORCE                               = "AmoebaTorsion"; 
static std::string AMOEBA_PI_TORSION_FORCE                            = "AmoebaPiTorsion"; 
static std::string AMOEBA_STRETCH_BEND_FORCE                          = "AmoebaStretchBend";
static std::string AMOEBA_OUT_OF_PLANE_BEND_FORCE                     = "AmoebaOutOfPlaneBend";
static std::string AMOEBA_TORSION_TORSION_FORCE                       = "AmoebaTorsionTorsion";
static std::string AMOEBA_MULTIPOLE_FORCE                             = "AmoebaMultipole";
static std::string AMOEBA_GK_FORCE                                    = "AmoebaGk";
static std::string AMOEBA_VDW_FORCE                                   = "AmoebaVdw";
static std::string AMOEBA_WCA_DISPERSION_FORCE                        = "AmoebaWcaDispersion";
static std::string AMOEBA_SASA_FORCE                                  = "AmoebaSASA";
static std::string ALL_FORCES                                         = "AllForces";

static std::string AMOEBA_MULTIPOLE_ROTATION_MATRICES                 = "AmoebaMultipoleRotationMatrices";
static std::string AMOEBA_MULTIPOLE_ROTATED                           = "AmoebaMultipolesRotated";
static std::string AMOEBA_FIXED_E                                     = "AmoebaFixedE";
static std::string AMOEBA_FIXED_E_GK                                  = "AmoebaFixedE_GK";
static std::string AMOEBA_INDUCDED_DIPOLES                            = "AmoebaInducedDipoles";
static std::string AMOEBA_INDUCDED_DIPOLES_GK                         = "AmoebaInducedDipoles_GK";

static std::string INCLUDE_OBC_CAVITY_TERM                            = "includeObcCavityTerm";
static std::string MUTUAL_INDUCED_MAX_ITERATIONS                      = "mutualInducedMaxIterations";
static std::string MUTUAL_INDUCED_TARGET_EPSILON                      = "mutualInducedTargetEpsilon";

#define AmoebaHarmonicBondIndex                            0
#define AmoebaHarmonicAngleIndex                           1
#define AmoebaHarmonicInPlaneAngleIndex                    2
#define AmoebaTorsionIndex                                 3
#define AmoebaPiTorsionIndex                               4
#define AmoebaStretchBendIndex                             5
#define AmoebaOutOfPlaneBendIndex                          6
#define AmoebaTorsionTorsionIndex                          7
#define AmoebaMultipoleIndex                               8
#define AmoebaVdwIndex                                     9
#define AmoebaWcaDispersionIndex                          10
#define AmoebaObcIndex                                    11
#define SumIndex                                          12
#define AmoebaLastIndex                                   13

#define BOLTZMANN  (1.380658e-23)                         /* (J/K) */
#define AVOGADRO   (6.0221367e23)                         /* ()    */
#define RGAS       (BOLTZMANN*AVOGADRO)                   /* (J/(mol K))  */
#define BOLTZ      (RGAS/1.0e+03)                         /* (kJ/(mol K)) */

#define AngstromToNm 0.1
#define CalToJoule   4.184

const double DegreesToRadians = 3.14159265/180.0;
const double RadiansToDegrees = 180/3.14159265;

using namespace OpenMM;
using namespace std;

// the following are used in parsing parameter file

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

typedef std::vector<StringVector> StringVectorVector;

typedef std::vector<std::vector<double> > VectorOfVectors;
typedef VectorOfVectors::iterator VectorOfVectorsI;
typedef VectorOfVectors::const_iterator VectorOfVectorsCI;

typedef std::map< std::string, VectorOfVectors > MapStringVectorOfVectors;
typedef MapStringVectorOfVectors::iterator MapStringVectorOfVectorsI;
typedef MapStringVectorOfVectors::const_iterator MapStringVectorOfVectorsCI;

typedef std::map< std::string, std::string > MapStringString;
typedef MapStringString::iterator MapStringStringI;
typedef MapStringString::const_iterator MapStringStringCI;

typedef std::map< std::string, int > MapStringInt;
typedef MapStringInt::iterator MapStringIntI;
typedef MapStringInt::const_iterator MapStringIntCI;

typedef std::map< std::string,  std::vector<Vec3> > MapStringVec3;
typedef MapStringVec3::iterator MapStringVec3I;
typedef MapStringVec3::const_iterator MapStringVec3CI;

typedef std::map< std::string, double > MapStringDouble;
typedef MapStringDouble::iterator MapStringDoubleI;
typedef MapStringDouble::const_iterator MapStringDoubleCI;

typedef std::map< std::string, Force*> MapStringForce;
typedef MapStringForce::iterator MapStringForceI;
typedef MapStringForce::const_iterator MapStringForceCI;

// default return value from methods

static const int DefaultReturnValue               = 0;


static const int LengthUnit                       = 0;
static const int EnergyUnit                       = 1;
static const int ForceUnit                        = 2;
static const int LastUnits                        = ForceUnit + 1;

static const int NoUnitsConversion                = 0;
static const int KcalA_To_kJNm                    = 1;

/**---------------------------------------------------------------------------------------
 * Initialize units
 *
 * @param unitType        has w/ force name as key and int as value
 * @param units array
 *
 *
    --------------------------------------------------------------------------------------- */

void setUnits( int unitType, double* units );

/**---------------------------------------------------------------------------------------

    Read parameter file

    @param inputParameterFile   input parameter file name
    @param system               system to which forces based on parameters are to be added
    @param coordinates          Vec3 array containing coordinates on output
    @param velocities           Vec3 array containing velocities on output
    @param inputLog             log file pointer -- may be NULL

    @return number of lines read

    --------------------------------------------------------------------------------------- */

Integrator* readAmoebaParameterFile( const std::string& inputParameterFile, MapStringInt& forceMap, System& system,
                                     std::vector<Vec3>& coordinates, 
                                     std::vector<Vec3>& velocities,
                                     MapStringVec3& forces, MapStringDouble& potentialEnergy,
                                     MapStringVectorOfVectors& supplementary, FILE* inputLog );


/**---------------------------------------------------------------------------------------
 * Get integrator
 * 
 * @param  integratorName       integratorName (VerletIntegrator, BrownianIntegrator, LangevinIntegrator, ...)
 * @param  timeStep             time step
 * @param  friction (ps)        friction
 * @param  temperature          temperature
 * @param  shakeTolerance       Shake tolerance
 * @param  errorTolerance       Error tolerance
 * @param  randomNumberSeed     seed
 *
 * @return DefaultReturnValue or ErrorReturnValue
 *
    --------------------------------------------------------------------------------------- */

Integrator* getIntegrator( std::string& integratorName, double timeStep,
                             double friction, double temperature,
                             double shakeTolerance, double errorTolerance,
                             int randomNumberSeed, FILE* log );

/**---------------------------------------------------------------------------------------
 * Initialize forceMap
 *
 * @param forceMap        has w/ force name as key and int as value
 * @param initialValue    initial value
 *
 *
    --------------------------------------------------------------------------------------- */

void initializeForceMap( MapStringInt& forceMap, int initialValue );

void testUsingAmoebaTinkerParameterFile( const std::string& amoebaTinkerParameterFileName, MapStringInt& forceMap,
                                         double tolerance, FILE* summaryFile, FILE* log );

int OPENMM_EXPORT runTestsUsingAmoebaTinkerParameterFile( MapStringString& argumentMap );
void OPENMM_EXPORT appendInputArgumentsToArgumentMap( int numberOfArguments, char* argv[], MapStringString& argumentMap );
