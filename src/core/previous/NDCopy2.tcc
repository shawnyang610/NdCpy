#ifndef NDCOPY_TCC
#define NDCOPY_TCC

#include "NDCopy2.h"

template <class T>
using Box = std::pair<T, T>;

bool GetOverlap(const Box<Dims> &b1, const Box<Dims> &b2,
                                     Box<Dims> &overlapBox);
void CopyLocalToGlobal(char *dst, const Box<Dims> &dstBox,
                                            const char *src,
                                            const Box<Dims> &srcBox,
                                            const size_t size,
                                            const Box<Dims> &overlapBox);
template<class T>
int NdCopy2(
        const Buffer &in_buffer,
        const Dims &in_start,
        const Dims &in_count,
        NdCopyFlag in_flag,
        Buffer &out_buffer,
        const Dims &out_start,
        const Dims &out_count,
        NdCopyFlag out_flag
        )
{

    Dims dstSecond(out_start.size());
    for(size_t i=0; i<out_start.size(); ++i){
        dstSecond[i] = out_start[i] + out_count[i];
    }
    Box<Dims> dstBox = {out_start, dstSecond};
    Dims srcSecond(in_start.size());
    for(size_t i=0; i<in_start.size(); ++i){
        srcSecond[i] = in_start[i] + in_count[i];
    }
    Box<Dims> srcBox = {in_start, srcSecond};
    Box<Dims> overlapBox;
    bool t = GetOverlap(dstBox, srcBox, overlapBox);
    CopyLocalToGlobal(out_buffer.data(), dstBox, in_buffer.data(), srcBox, sizeof(T), overlapBox);
    return 0;
}

#endif
