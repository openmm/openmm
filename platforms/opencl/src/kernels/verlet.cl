/**
 * Perform the first step of verlet integration.
 */

__kernel void integrateVerletPart1(int numAtoms, float dt, __global float4* posq, __global float4* velm, __global float4* force, __global float4* posDelta) {
    int index = get_global_id(0);
    while (index < numAtoms) {
        float4 pos = posq[index];
        float4 velocity = velm[index];
        velocity.xyz += force[index].xyz*dt*velocity.w;
        pos.xyz = velocity.xyz*dt;
        posDelta[index] = pos;
        velm[index] = velocity;
        index += get_global_size(0);
    }
}

/**
 * Perform the second step of verlet integration.
 */

__kernel void integrateVerletPart2(int numAtoms, float dt, __global float4* posq, __global float4* velm, __global float4* posDelta) {
    int index = get_global_id(0);
    float oneOverDt = 1.0f/dt;
    while (index < numAtoms) {
        float4 pos = posq[index];
        float4 delta = posDelta[index];
        float4 velocity = velm[index];
        pos.xyz += delta.xyz;
        velocity.xyz = delta.xyz*oneOverDt;
        posq[index] = pos;
        velm[index] = velocity;
        index += get_global_size(0);
    }
}
