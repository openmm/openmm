#if defined(__ANDROID__)
    #include "neon_mathfun.h"
#else
    #if !defined(__PNACL__)
        #include "sse_mathfun.h"
    #endif
#endif
