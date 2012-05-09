/**
 * Compute the difference between two vectors, setting the fourth component to the squared magnitude.
 */
float4 ccb_delta(float4 vec1, float4 vec2) {
    float4 result = (float4) (vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z, 0.0f);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}

/**
 * Compute the angle between two vectors.  The w component of each vector should contain the squared magnitude.
 */
float ccb_computeAngle(float4 vec1, float4 vec2) {
    float dotProduct = vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    float cosine = dotProduct*RSQRT(vec1.w*vec2.w);
    float angle;
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        float4 crossProduct = cross(vec1, vec2);
        float scale = vec1.w*vec2.w;
        angle = asin(SQRT(dot(crossProduct, crossProduct)/scale));
        if (cosine < 0.0f)
            angle = M_PI-angle;
    }
    else
       angle = acos(cosine);
    return angle;
}

/**
 * Compute the cross product of two vectors, setting the fourth component to the squared magnitude.
 */
float4 ccb_computeCross(float4 vec1, float4 vec2) {
    float4 result = cross(vec1, vec2);
    result.w = result.x*result.x + result.y*result.y + result.z*result.z;
    return result;
}
