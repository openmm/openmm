void  kshakeh_fix1 (const float  nit,
		const float  strwidth,
		const float  invmH,
		const float  omega,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream params,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3);


void  kshakeh_fix2 (const float  nit,
		const float  strwidth,
		const float  invmH,
		const float  omega,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream params,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3);



void  kshakeh_update (const float  strwidth,
		::brook::stream invmap,
		::brook::stream posq,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3,
		::brook::stream oposq) ;


void  kshakeh (const float  nit,
		const float  strwidth,
		const float  invmH,
		const float  omega,
		::brook::stream atoms,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream params,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3); 

void  kshakeh_update1_fix1 (
      const float  strwidth,
      const float  sdpc1,
		::brook::stream invmap,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream vPrime,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3,
		::brook::stream oposq); 

void  kshakeh_update1_fix1Old (const float  strwidth,
		::brook::stream invmap,
		::brook::stream posq,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3,
		::brook::stream oposq); 

void  kshakeh_update2_fix1 (const float  strwidth,
		::brook::stream invmap,
		::brook::stream posq,
		::brook::stream posqp,
		::brook::stream cposq0,
		::brook::stream cposq1,
		::brook::stream cposq2,
		::brook::stream cposq3,
		::brook::stream oposq); 
