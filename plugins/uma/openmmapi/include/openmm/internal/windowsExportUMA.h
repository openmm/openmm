#ifndef OPENMM_WINDOWSEXPORTUMA_H_
#define OPENMM_WINDOWSEXPORTUMA_H_

/* Global flags for building DLLs on Windows */
#if defined(_MSC_VER)
    #define OPENMM_EXPORT_UMA __declspec(dllimport)
    #define OPENMM_EXPORT_UMA_COMMON __declspec(dllimport)
#else
    #define OPENMM_EXPORT_UMA
    #define OPENMM_EXPORT_UMA_COMMON
#endif

#endif /*OPENMM_WINDOWSEXPORTUMA_H_*/
