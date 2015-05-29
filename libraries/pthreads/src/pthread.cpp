#include "pthread.h"

PTHREAD_EXPORT volatile long _pthread_cancelling;

PTHREAD_EXPORT int _pthread_concur;

/* Will default to zero as needed */
PTHREAD_EXPORT pthread_once_t _pthread_tls_once;
PTHREAD_EXPORT DWORD _pthread_tls;

/* Note initializer is zero, so this works */
PTHREAD_EXPORT pthread_rwlock_t _pthread_key_lock;
PTHREAD_EXPORT long _pthread_key_max;
PTHREAD_EXPORT long _pthread_key_sch;
PTHREAD_EXPORT void (**_pthread_key_dest)(void *);
