void  kupdate_md2 (const float  dtinv,
		::brook::stream posqp,
		::brook::stream posq,
		::brook::stream vnew,
		::brook::stream posqnew
	    ); 

void  kupdate_md1_extended (
		const float  dt,
		const float3  lg,
		const float  xi,
		const float3  M0,
		const float3  M1,
		const float3  M2,
		const float3  uold,
		::brook::stream posq,
		::brook::stream v,
		::brook::stream f,
		::brook::stream invmass,
		::brook::stream vnew,
		::brook::stream posqp);

void  kupdate_md_verlet(
		const float  dt,
		::brook::stream posq,
		::brook::stream v,
		::brook::stream f,
		::brook::stream invmass,
		::brook::stream vnew,
		::brook::stream posqp);


void  kupdate_md1_berendsen (
		const float  dt,
		const float3  lg,
		const float3  uold,
		::brook::stream posq,
		::brook::stream v,
		::brook::stream f,
		::brook::stream invmass,
		::brook::stream vnew,
		::brook::stream posqp); 

void  kupdate_md2_fix1 (const float  dtinv,
		::brook::stream posqp,
		::brook::stream posq,
		::brook::stream vnew,
		::brook::stream posqnew); 

void  kupdate_md1_extended_fix1 (const float  dt,
		const float3  lg,
		const float  xi,
		const float3  M0,
		const float3  M1,
		const float3  M2,
		const float3  uold,
		::brook::stream posq,
		::brook::stream v,
		::brook::stream f,
		::brook::stream invmass,
		::brook::stream vnew,
		::brook::stream posqp);

void  kupdate_md1_berendsen_fix1 (const float  dt,
		const float3  lg,
		const float3  uold,
		::brook::stream posq,
		::brook::stream v,
		::brook::stream f,
		::brook::stream invmass,
		::brook::stream vnew,
		::brook::stream posqp); 
