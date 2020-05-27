#if defined(__ARM__) || defined(__ARM64__)
    #include "neon_mathfun.h"
#else
    #if !defined(__PNACL__)
        #define USE_SSE2
        #include "sse_mathfun.h"
    #endif
#endif
