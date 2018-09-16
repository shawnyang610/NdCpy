#ifndef NDCOPY_H
#define NDCOPY_H

#include <vector>
#include <iostream>

using Dims = std::vector<size_t>;
using Buffer = std::vector<char>;

struct NdCopyFlag {
    bool isRowMajor;
    bool isBigEndian;
};

template <typename T>
int NdCopy(const char* in,
           const Dims &inStart,
           const Dims &inCount,
           const bool inIsRowMajor,
           const bool inIsBigEndian,
           char* out,
           const Dims &outStart,
           const Dims &outCount,
           const bool outIsRowMajor,
           const bool outIsBigEndian,
           const bool safeMode=false
           );

template<typename T>
int NdCopy2(
        const Buffer &in_buffer,
        const Dims &in_start,
        const Dims &in_count,
        NdCopyFlag in_flag,
        Buffer &out_buffer,
        const Dims &out_start,
        const Dims &out_count,
        NdCopyFlag out_flag
        );

#endif
