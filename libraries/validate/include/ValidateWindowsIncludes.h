#ifndef OPENMM_VALIDATE_WINDOW_INCLUDE_H_
#define OPENMM_VALIDATE_WINDOW_INCLUDE_H_

/*
 * Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the OpenMMValidate shared library
 *       (dllexport)
 *   (2) this header is being used by a *client* of the OpenMMValidate shared
 *       library (dllimport)
 *   (3) we are building the OpenMMValidate static library, or the client is
 *       being compiled with the expectation of linking with the
 *       OpenMMValidate static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     OpenMMValidate_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol OPENMM_VALIDATE_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. Note that
 * the client symbol is not library dependent, while the library symbols
 * affect only the OpenMMValidate library, meaning that other libraries can
 * be clients of this one. However, we are assuming all-static or all-shared.
 */

#ifdef _MSC_VER
    // We don't want to hear about how sprintf is "unsafe".
    #pragma warning(disable:4996)
    #if defined(OPENMM_VALIDATE_BUILDING_SHARED_LIBRARY)
        #define OPENMM_VALIDATE_EXPORT __declspec(dllexport)
        // Keep MS VC++ quiet about lack of dll export of private members.
        #pragma warning(disable:4251)
    #elif defined(OPENMM_VALIDATE_BUILDING_STATIC_LIBRARY) || defined(OPENMM_VALIDATE_USE_STATIC_LIBRARIES)
		#define OPENMM_VALIDATE_EXPORT
    #else
		#define OPENMM_VALIDATE_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define OPENMM_VALIDATE_EXPORT // Linux, Mac
#endif

#endif // OPENMM_VALIDATE_WINDOW_INCLUDE_H_
