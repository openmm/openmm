/**
 * Convert a real4 to a real3 by removing its last element.
 */
inline __device__ real3 ccb_trim(real4 v) {
    return make_real3(v.x, v.y, v.z);
}

/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
inline __device__ real4 ccb_delta(real4 vec1, real4 vec2, bool periodic, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = make_real4(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
    if (periodic)
        APPLY_PERIODIC_TO_DELTA(result);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
__device__ real ccb_computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real3 crossProduct = cross(vec1, vec2);
        real scale = vec1.w*vec2.w;
        angle = ASIN(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0.0f)
            angle = M_PI-angle;
    }
    else
       angle = ACOS(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
inline __device__ real4 ccb_computeCross(real4 vec1, real4 vec2) {
    real3 cp = cross(vec1, vec2);
    return make_real4(cp.x, cp.y, cp.z, cp.x*cp.x+cp.y*cp.y+cp.z*cp.z);
}
