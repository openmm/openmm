/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
real4 ccb_delta(real4 vec1, real4 vec2, bool periodic, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ) {
    real4 result = (real4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0);
    if (periodic)
        APPLY_PERIODIC_TO_DELTA(result);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
real ccb_computeAngle(real4 vec1, real4 vec2) {
    real dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    real cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    real angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        real4 crossProduct = cross(vec1, vec2);
        real scale = vec1.w*vec2.w;
        angle = asin(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0)
            angle = M_PI-angle;
    }
    else
       angle = acos(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
real4 ccb_computeCross(real4 vec1, real4 vec2) {
    real4 result = cross(vec1, vec2);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}
