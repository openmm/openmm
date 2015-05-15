#include "pthread.h"

volatile long _pthread_cancelling;

int _pthread_concur;

/* Will default to zero as needed */
pthread_once_t _pthread_tls_once;
DWORD _pthread_tls;

/* Note initializer is zero, so this works */
pthread_rwlock_t _pthread_key_lock;
long _pthread_key_max;
long _pthread_key_sch;
void (**_pthread_key_dest)(void *);