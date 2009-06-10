// -----------------------------------------------------------------------------
//        OpenMM(tm) example C and Fortran wrapper functions (June 2009)
// -----------------------------------------------------------------------------
// This is the C++ implementation of the C wrappers for the OpenMM workshop.
// The functions here convert between C types and OpenMM's C++ objects
// and then call the appropriate OpenMM methods.
//
// Each C function comes in two forms -- one is intended to be called from C
// main programs, and the other from Fortran main programs. The Fortran one
// typically just translates Fortran naming and argument conventions into C 
// and then calls the C function.
//
// A C main program can use this just by including the OpenMM_CWrapper.h 
// header file that is included here as well. A Fortran 95 program can use
// the "use OpenMM" module which defines an interface to the Fortran-callable
// functions defined here. Fortran 77 programs have to call these directly;
// you can use an integer*8 to hold the pointers.
// -----------------------------------------------------------------------------

#include "OpenMM_CWrapper.h"

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
    #pragma warning(disable:4251)   // no dll interface for some classes
#endif

#include "OpenMM.h"
using namespace OpenMM;


static inline Vec3 toVec3(const OpenMM_Vec3 src) {
    return Vec3(src[0], src[1], src[2]);
}
static inline void fromVec3(const Vec3& src, OpenMM_Vec3 dest) {
    dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2];
}

extern "C" {

    ///////////////////////////////
    // std::vector<OpenMM::Vec3> //
    ///////////////////////////////

// This provides std::vector<Vec3> functionality to the
// C program. It isn't as elegant as in C++ but is still better
// than doing it in C. Also, this allows the C program to
// communicate with OpenMM without having to copy arrays of
// Vec3s to and from std::vectors.

OpenMM_Vec3Array* OpenMM_Vec3Array_create(int n) 
{   return (OpenMM_Vec3Array*)new std::vector<Vec3>(n); }
void openmm_vec3array_create_(OpenMM_Vec3Array*& a, const int& n)
{   a = OpenMM_Vec3Array_create(n); }

int OpenMM_Vec3Array_size(const OpenMM_Vec3Array* a) 
{   return (int)((const std::vector<Vec3>*)a)->size(); }
int openmm_vec3array_size_(const OpenMM_Vec3Array*& a)
{   return OpenMM_Vec3Array_size(a); }

void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* a, int n) 
{   ((std::vector<Vec3>*)a)->resize(n); }
void openmm_vec3array_resize_(OpenMM_Vec3Array* const& a, const int& n)
{   OpenMM_Vec3Array_resize(a, n); }

void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* doomed) 
{   delete ((std::vector<Vec3>*)doomed); }
void openmm_vec3array_destroy_(OpenMM_Vec3Array*& doomed)
{   OpenMM_Vec3Array_destroy(doomed); doomed = 0; }

void OpenMM_Vec3Array_append(OpenMM_Vec3Array* a, const OpenMM_Vec3 v) 
{   ((std::vector<Vec3>*)a)->push_back(*(const Vec3*)v); }
void openmm_vec3array_append_(OpenMM_Vec3Array* const& a, const OpenMM_Vec3 v)
{   OpenMM_Vec3Array_append(a, v); }

// Get a single Vec3 element from the array. Index is 0-relative in C, 1-relative in Fortran.
void OpenMM_Vec3Array_get(const OpenMM_Vec3Array* a, int i0, OpenMM_Vec3 ov3) {
    fromVec3((*(const std::vector<Vec3>*)a)[i0], ov3);
}
void openmm_vec3array_get_(const OpenMM_Vec3Array* const& a, const int& i1, OpenMM_Vec3 ov3)
{   OpenMM_Vec3Array_get(a, i1-1, ov3); }


    /////////////////
    // std::string //
    /////////////////

// This is an interface to std::string primarily for Fortran. You
// can use null-terminated char arrays directly in C.

OpenMM_String* OpenMM_String_create(const char* nullTerminatedInitVal) {
    OpenMM_String* os = (OpenMM_String*)new std::string(nullTerminatedInitVal);
    return os;
}
void openmm_string_create_(OpenMM_String*& os, const char* init, int len) {
    std::string* s = new std::string();
    os = (OpenMM_String*)s;
    if (len > 0) {
        s->resize(len);
        std::strncpy(&(*s)[0], init, len);
    }
}
void OpenMM_String_destroy(OpenMM_String* os) {
    delete ((std::string*)os);
}
void openmm_string_destroy_(OpenMM_String*& os) {OpenMM_String_destroy(os);}
const char* OpenMM_String_getAsC(const OpenMM_String* os) {
    return ((const std::string*)os)->c_str();
}
int OpenMM_String_length(const OpenMM_String* os) {
    return (int)((const std::string*)os)->size();
}
int openmm_string_length_(const OpenMM_String* const& os) {
    return OpenMM_String_length(os);
}
// Copy out as a null-terminated C string.
void OpenMM_String_get(const OpenMM_String* os, char* buf, int buflen) {
    if (buflen <= 0) return;
    const std::string& s = *(const std::string*)os;
    const int minlen = std::min((int)s.size(), buflen);
    for (int i=0; i < minlen; ++i)
        buf[i] = s[i];
    const int nullpos = std::min(minlen, buflen-1);
    buf[nullpos] = '\0';
}
// Copy out as a blank-padded Fortran string.
void openmm_string_get_(const OpenMM_String* const& os, char* buf, int buflen) {
    if (buflen <= 0) return;
    const std::string& s = *(const std::string*)os;
    const int minlen = std::min((int)s.size(), buflen);
    for (int i=0; i < minlen; ++i)
        buf[i] = s[i];
    for (int i=minlen; i < buflen; ++i)
        buf[i] = ' ';
}

// Set string from a null-terminated C string, stripping trailing blanks.
void OpenMM_String_set(OpenMM_String* os, const char* in) {
    std::string& s = *(std::string*)os;
    int len = std::strlen(in);
    s = std::string(in, in+len);
    while (len > 0 && s[len-1]==' ')
        --len;
    s.erase(len);
}
// Set string from a fix-sized Fortran character array, 
// stripping trailing blanks.
void openmm_string_set_(OpenMM_String*& os, const char* in, int len) {
    std::string& s = *(std::string*)os;
    s = std::string(in, in+len);
    while (len > 0 && s[len-1]==' ')
        --len;
    s.erase(len);
}


    //////////////////////
    // OpenMM::Platform //
    //////////////////////
void OpenMM_Platform_loadPluginsFromDirectory(const char* dir) {
    OpenMM::Platform::loadPluginsFromDirectory(std::string(dir));
}
const char* OpenMM_Platform_getDefaultPluginsDirectory() {
    static std::string dir;
    dir = OpenMM::Platform::getDefaultPluginsDirectory();
    const char* out = dir.c_str();
    return dir.c_str();
}

void openmm_platform_loadpluginsfromdirectory_(const OpenMM_String* const& dir)
{    OpenMM_Platform_loadPluginsFromDirectory(OpenMM_String_getAsC(dir)); }

void openmm_platform_getdefaultpluginsdirectory_(OpenMM_String* const& dir)
{    OpenMM_String_set(dir, OpenMM_Platform_getDefaultPluginsDirectory()); }


    ////////////////////
    // OpenMM::System //
    ////////////////////
OpenMM_System* 
OpenMM_System_create() {
    return (OpenMM_System*)new System();
}
void openmm_system_create_(OpenMM_System*& sys) {sys=OpenMM_System_create();}

void OpenMM_System_destroy(OpenMM_System* doomed) {
    delete (System*)doomed;
}
void openmm_system_destroy_(OpenMM_System*& doomed) 
{OpenMM_System_destroy(doomed); doomed=0;}

void OpenMM_System_addForce(OpenMM_System* sys, OpenMM_Force* frc) {
    ((System*)sys)->addForce((NonbondedForce*)frc);
}
void openmm_system_addforce_(OpenMM_System*& sys, OpenMM_Force*& frc) 
{OpenMM_System_addForce(sys,frc);}

void OpenMM_System_addParticle(OpenMM_System* sys, double mass) {
    ((System*)sys)->addParticle(mass);
}
void openmm_system_addparticle_(OpenMM_System*& sys, const double& mass) 
{OpenMM_System_addParticle(sys,mass);}

int OpenMM_System_getNumParticles(const OpenMM_System* sys) {
    return ((const System*)sys)->getNumParticles();
}
int openmm_system_getnumparticles_(const OpenMM_System*& sys)
{return OpenMM_System_getNumParticles(sys);}

    ////////////////////////////
    // OpenMM::NonbondedForce //
    ////////////////////////////
OpenMM_NonbondedForce* OpenMM_NonbondedForce_create() 
{   return (OpenMM_NonbondedForce*)new NonbondedForce(); }
void openmm_nonbondedforce_create_(OpenMM_NonbondedForce*& frc)
{   frc = OpenMM_NonbondedForce_create();}

void OpenMM_NonbondedForce_destroy(OpenMM_NonbondedForce* doomed) 
{   delete (NonbondedForce*)doomed; }
void openmm_nonbondedforce_destroy_(OpenMM_NonbondedForce*& doomed) 
{   OpenMM_NonbondedForce_destroy(doomed); doomed = 0;}

// Fortran only: recast NonbondedForce as a Force.
void openmm_nonbondedforce_asforce_(OpenMM_NonbondedForce* const& nonbond,
									OpenMM_Force*&                force)
{   force = (OpenMM_Force*)nonbond; }

void OpenMM_NonbondedForce_setNonbondedMethod(OpenMM_NonbondedForce* nbf, 
                                              OpenMM_NonbondedForce_NonbondedMethod method) 
{    ((NonbondedForce*)nbf)->setNonbondedMethod(NonbondedForce::NonbondedMethod(method)); }
void openmm_nonbondedforce_setnonbondedmethod_(OpenMM_NonbondedForce*& nbf, const int& method)
{   OpenMM_NonbondedForce_setNonbondedMethod(nbf,OpenMM_NonbondedForce_NonbondedMethod(method));}

void OpenMM_NonbondedForce_setCutoffDistance(OpenMM_NonbondedForce* nbf, double d) 
{   ((NonbondedForce*)nbf)->setCutoffDistance(d); }
void openmm_nonbondedforce_setcutoffdistance_(OpenMM_NonbondedForce*& nbf, const double& d)
{   OpenMM_NonbondedForce_setCutoffDistance(nbf,d);}

void OpenMM_NonbondedForce_setPeriodicBoxVectors(OpenMM_NonbondedForce* nbf, 
                                                 const OpenMM_Vec3 a,const OpenMM_Vec3 b,const OpenMM_Vec3 c) 
{
    ((NonbondedForce*)nbf)->setPeriodicBoxVectors(*(const Vec3*)a, *(const Vec3*)b, *(const Vec3*)c);
}
void openmm_nonbondedforce_setperiodicboxvectors_(OpenMM_NonbondedForce*& nbf, 
                                                  const OpenMM_Vec3 a,const OpenMM_Vec3 b,const OpenMM_Vec3 c) 
{   OpenMM_NonbondedForce_setPeriodicBoxVectors(nbf,a,b,c);}

void OpenMM_NonbondedForce_addParticle(OpenMM_NonbondedForce* nbf, 
                                       double charge, double sigmaInNm, double vdwEnergyInKJ)
{
    ((NonbondedForce*)nbf)->addParticle(charge, sigmaInNm, vdwEnergyInKJ);
}
void openmm_nonbondedforce_addparticle_(OpenMM_NonbondedForce*& nbf, 
                                        const double& charge, const double& sigmaInNm, const double& vdwEnergyInKJ)
{   OpenMM_NonbondedForce_addParticle(nbf,charge, sigmaInNm, vdwEnergyInKJ);}

    //////////////////////////
    // OpenMM::GBSAOBCForce //
    //////////////////////////
OpenMM_GBSAOBCForce* OpenMM_GBSAOBCForce_create() 
{   return (OpenMM_GBSAOBCForce*)new GBSAOBCForce(); }
void openmm_gbsaobcforce_create_(OpenMM_GBSAOBCForce*& frc)
{   frc = OpenMM_GBSAOBCForce_create();}

void OpenMM_GBSAOBCForce_destroy(OpenMM_GBSAOBCForce* doomed) 
{   delete (GBSAOBCForce*)doomed; }
void openmm_gbsaobcforce_destroy_(OpenMM_GBSAOBCForce*& doomed) 
{   OpenMM_GBSAOBCForce_destroy(doomed); doomed = 0;}

// Fortran only: recast NonbondedForce as a Force.
void openmm_gbsaobcforce_asforce_(OpenMM_GBSAOBCForce* const& gbsa,
							      OpenMM_Force*&              force)
{   force = (OpenMM_Force*)gbsa; }

void OpenMM_GBSAOBCForce_setSolventDielectric(OpenMM_GBSAOBCForce* gbsa, double d) 
{   ((GBSAOBCForce*)gbsa)->setSolventDielectric(d); }
void openmm_gbsaobcforce_setsolventdielectric_(OpenMM_GBSAOBCForce*& gbsa, const double& d)
{   OpenMM_GBSAOBCForce_setSolventDielectric(gbsa,d);}

void OpenMM_GBSAOBCForce_setSoluteDielectric(OpenMM_GBSAOBCForce* gbsa, double d) 
{   ((GBSAOBCForce*)gbsa)->setSoluteDielectric(d); }
void openmm_gbsaobcforce_setsolutedielectric_(OpenMM_GBSAOBCForce*& gbsa, const double& d)
{   OpenMM_GBSAOBCForce_setSoluteDielectric(gbsa,d);}

void OpenMM_GBSAOBCForce_addParticle(OpenMM_GBSAOBCForce* gbsa, 
                                     double charge, double radiusInNm, double scalingFactor)
{
    ((GBSAOBCForce*)gbsa)->addParticle(charge, radiusInNm, scalingFactor);
}
void openmm_gbsaobcforce_addparticle_(OpenMM_GBSAOBCForce*& gbsa, 
                                      const double& charge, const double& radiusInNm, const double& scalingFactor)
{   OpenMM_GBSAOBCForce_addParticle(gbsa,charge, radiusInNm, scalingFactor);}



    ////////////////////////
    // OpenMM::Integrator //
    ////////////////////////
void OpenMM_Integrator_step(OpenMM_Integrator* integ, int numSteps) 
{    ((Integrator*)integ)->step(numSteps); }
void openmm_integrator_step_(OpenMM_Integrator* const& integ, int& numSteps) 
{    OpenMM_Integrator_step(integ, numSteps); }

void OpenMM_Integrator_destroy(OpenMM_Integrator* doomed) 
{   delete ((Integrator*)doomed); }
void openmm_integrator_destroy_(OpenMM_Integrator*& doomed)
{   OpenMM_Integrator_destroy(doomed); doomed = 0; }

    // OpenMM::VerletIntegrator
OpenMM_VerletIntegrator* OpenMM_VerletIntegrator_create(double stepSzInPs) 
{   return (OpenMM_VerletIntegrator*)new VerletIntegrator(stepSzInPs); }
void openmm_verletintegrator_create_(OpenMM_VerletIntegrator*& verlet, double& stepSzInPs)
{   verlet = OpenMM_VerletIntegrator_create(stepSzInPs); }

void OpenMM_VerletIntegrator_destroy(OpenMM_VerletIntegrator* doomed) 
{   delete (VerletIntegrator*)doomed; }
void openmm_verletintegrator_destroy_(OpenMM_VerletIntegrator*& doomed)
{   OpenMM_VerletIntegrator_destroy(doomed); doomed = 0; }

// Fortran only: recast VerletIntegrator as an Integrator.
void openmm_verletintegrator_asintegrator_(OpenMM_VerletIntegrator* const& verlet,
										   OpenMM_Integrator*&             integ)
{   integ = (OpenMM_Integrator*)verlet; }

void OpenMM_VerletIntegrator_step(OpenMM_VerletIntegrator* verlet, int numSteps) 
{   ((VerletIntegrator*)verlet)->step(numSteps); }
void openmm_verletintegrator_step_(OpenMM_VerletIntegrator* const& verlet, int& numSteps)
{   OpenMM_VerletIntegrator_step(verlet, numSteps); }


    // OpenMM::LangevinIntegrator
OpenMM_LangevinIntegrator* OpenMM_LangevinIntegrator_create(double temperature, double frictionInPerPs, double stepSzInPs) 
{   return (OpenMM_LangevinIntegrator*)new LangevinIntegrator(temperature, frictionInPerPs, stepSzInPs); }
void openmm_langevinintegrator_create_(OpenMM_LangevinIntegrator*& langevin, double& temperature, double& frictionInPerPs, double& stepSzInPs)
{   langevin = OpenMM_LangevinIntegrator_create(temperature, frictionInPerPs, stepSzInPs); }

void OpenMM_LangevinIntegrator_destroy(OpenMM_LangevinIntegrator* doomed) 
{   delete (LangevinIntegrator*)doomed; }
void openmm_langevinintegrator_destroy_(OpenMM_LangevinIntegrator*& doomed)
{   OpenMM_LangevinIntegrator_destroy(doomed); doomed = 0; }

// Fortran only: recast LangevinIntegrator as an Integrator.
void openmm_langevinintegrator_asintegrator_(OpenMM_LangevinIntegrator* const& langevin,
										     OpenMM_Integrator*&               integ)
{   integ = (OpenMM_Integrator*)langevin; }

void OpenMM_LangevinIntegrator_step(OpenMM_LangevinIntegrator* langevin, int numSteps) 
{   ((LangevinIntegrator*)langevin)->step(numSteps); }
void openmm_langevinintegrator_step_(OpenMM_LangevinIntegrator* const& langevin, int& numSteps)
{   OpenMM_LangevinIntegrator_step(langevin, numSteps); }

    /////////////////////
    // OpenMM::Context //
    /////////////////////
OpenMM_Context* OpenMM_Context_create(OpenMM_System* sys, OpenMM_Integrator* integ) {
    return (OpenMM_Context*)new OpenMM::OpenMMContext(*(System*)sys, *(Integrator*)integ);
}
void openmm_context_create_(OpenMM_Context*& context, OpenMM_System*& sys, OpenMM_Integrator*& integ)
{   context = OpenMM_Context_create(sys, integ); }

void OpenMM_Context_destroy(OpenMM_Context* doomed) {
    delete (OpenMMContext*)doomed;
}
void openmm_context_destroy_(OpenMM_Context*& doomed) 
{    OpenMM_Context_destroy(doomed); }

void OpenMM_Context_setPositions(OpenMM_Context* context, const OpenMM_Vec3Array* positions) {
    ((OpenMMContext*)context)->setPositions(*(const std::vector<Vec3>*)positions);
}
void openmm_context_setpositions_(OpenMM_Context*& context, const OpenMM_Vec3Array* const& positions)
{    OpenMM_Context_setPositions(context, positions); }

void OpenMM_Context_setVelocities(OpenMM_Context* context, const OpenMM_Vec3Array* velocities) {
    ((OpenMMContext*)context)->setVelocities(*(const std::vector<Vec3>*)velocities);
}
void openmm_context_setvelocities_(OpenMM_Context*& context, const OpenMM_Vec3Array* const& velocities)
{    OpenMM_Context_setVelocities(context, velocities); }

// Note that a Context creates the OpenMM::State object, but you have to destroy
// it using OpenMM_State_destroy.
OpenMM_State* OpenMM_Context_createState(const OpenMM_Context* context, int types) {
    return (OpenMM_State*)new State(((OpenMMContext*)context)->getState(types));
}
void openmm_context_createstate_(const OpenMM_Context* const& context, const int& types, OpenMM_State*& state)
{   state=OpenMM_Context_createState(context, types); }

// Return a reference to a static null terminated C string containing the
// Platform name.
const char* OpenMM_Context_getPlatformName(const OpenMM_Context* context) {
    static std::string platform;
    platform = ((const OpenMMContext*)context)->getPlatform().getName();
    return platform.c_str();
}
// Return a blank-padded Fortran string containing the Platform name. There
// is no terminating null.
void openmm_context_getplatformname_(const OpenMM_Context* const& context, char* buf, int len) {
    const std::string name = ((const OpenMMContext*)context)->getPlatform().getName();
    const int minLen = std::min((int)name.size(), len);
    for (int i=0; i<minLen; ++i) buf[i] = name[i];
    for (int i=minLen; i<len; ++i) buf[i] = ' ';
}
	
double OpenMM_Context_getTime(OpenMM_Context* context) {
    return ((OpenMMContext*)context)->getTime();
}
double openmm_context_gettime_(OpenMM_Context* const& context)
{   return OpenMM_Context_getTime(context); }

    ///////////////////
    // OpenMM::State //
    ///////////////////
void OpenMM_State_destroy(OpenMM_State* doomed) 
{   delete (State*)doomed; }
void openmm_state_destroy_(OpenMM_State*& doomed)
{   OpenMM_State_destroy(doomed); doomed=0; }

double OpenMM_State_getTime(const OpenMM_State* state) 
{   return ((const State*)state)->getTime(); }
double openmm_state_gettime_(const OpenMM_State* const& state)
{   return OpenMM_State_getTime(state); }

double OpenMM_State_getPotentialEnergy(const OpenMM_State* state) 
{   return ((const State*)state)->getPotentialEnergy(); }
double openmm_state_getpotentialenergy_(const OpenMM_State* const& state)
{   return OpenMM_State_getPotentialEnergy(state); }

double OpenMM_State_getKineticEnergy(const OpenMM_State* state) 
{   return ((const State*)state)->getKineticEnergy(); }
double openmm_state_getkineticenergy_(const OpenMM_State* const& state)
{   return OpenMM_State_getKineticEnergy(state); }

const OpenMM_Vec3Array* OpenMM_State_getPositions(const OpenMM_State* state) 
{   return (const OpenMM_Vec3Array*)&((const State*)state)->getPositions(); }
void openmm_state_getpositions_(const OpenMM_State* const& state, const OpenMM_Vec3Array*& positions)
{   positions = OpenMM_State_getPositions(state); }

const OpenMM_Vec3Array*  OpenMM_State_getVelocities(const OpenMM_State* state) 
{   return (const OpenMM_Vec3Array*)&((const State*)state)->getVelocities(); }
void openmm_state_getvelocities_(const OpenMM_State* const& state, const OpenMM_Vec3Array*& velocities)
{   velocities = OpenMM_State_getVelocities(state); }


} // extern "C"
