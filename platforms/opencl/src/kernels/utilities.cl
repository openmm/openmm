__kernel void clearBuffer(__global float* buffer, int size) {
    int index = get_global_id(0);
    int step = get_global_size(0);
    __global float4* buffer4 = (__global float4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = (float4) (0.0f);
        index += step;
    }
    if (get_global_id(0) == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0.0f;
}
