// ISPC batch AABB filter
export uniform int bbox_overlap_batch_f32(
    uniform float qminx, uniform float qminy, uniform float qminz,
    uniform float qmaxx, uniform float qmaxy, uniform float qmaxz,
    const uniform float minx[], const uniform float miny[], const uniform float minz[],
    const uniform float maxx[], const uniform float maxy[], const uniform float maxz[],
    uniform int n,
    uniform int out_idx[],
    uniform int out_cap)
{
    uniform int w = 0;
    for (uniform int i = 0; i < n && w < out_cap; ++i) {
        bool hit = (qminx <= maxx[i]) && (minx[i] <= qmaxx) &&
                   (qminy <= maxy[i]) && (miny[i] <= qmaxy) &&
                   (qminz <= maxz[i]) && (minz[i] <= qmaxz);
        if (hit)
            out_idx[w++] = i;
    }
    return w;
}
