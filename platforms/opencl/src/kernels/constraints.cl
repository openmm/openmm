__kernel void applyPositionDeltas(__global float4* restrict posq, __global float4* restrict posDelta) {
    for (unsigned int index = get_global_id(0); index < NUM_ATOMS; index += get_global_size(0)) {
        float4 position = posq[index];
        position.xyz += posDelta[index].xyz;
        posq[index] = position;
    }
}
