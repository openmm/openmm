
void  kCalculateLinearMomentum( ::brook::stream mass, ::brook::stream velocities, ::brook::stream linearMomentum );
void  kReduceLinearMomentum( ::brook::stream momentum, ::brook::stream linearMomentum );
//void  kSumLinearMomentum( ::brook::stream momentum, float3* linearMomentum );
void  kScale( float scale, ::brook::stream linearMomentumIn, ::brook::stream linearMomentumOut );
void  kRemoveLinearMomentum( ::brook::stream linearMomentum, ::brook::stream velocitiesIn, ::brook::stream velocitiesOut );
void kSum( ::brook::stream array, ::brook::stream sum );

/**---------------------------------------------------------------------------------------

   This kernel calculates the total linear momentum via a reduction

   @param atomStrWidth   atom stream width
   @param numberOfAtoms  number of atoms
   @param momentum       momentum
   @param linearMomentum total momentum

   --------------------------------------------------------------------------------------- */

void kSumLinearMomentum( float atomStrWidth, float numberOfAtoms, ::brook::stream momentum, ::brook::stream linearMomentum );


