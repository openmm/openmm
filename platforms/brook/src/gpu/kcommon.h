#ifndef __KCOMMON_H__
#define __KCOMMON_H__

void  kgetxyz (::brook::stream instr,
		::brook::stream outstr); 

void  kzerof3 (::brook::stream outstr);
void  kzerof4 (::brook::stream outstr); 
void  kzerof4 (::brook::stream outstr); 

void  ksetf4 (const float4  val, ::brook::stream outstr); 
void kadd3( ::brook::stream instr, ::brook::stream outstr );
void ksetStr3( ::brook::stream instr, ::brook::stream outstr );

#endif // __KCOMMON_H__
