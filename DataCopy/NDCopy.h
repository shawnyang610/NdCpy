#ifndef NDCOPY_H
#define NDCOPY_H

#include <vector>
#include <iostream>

using Dims = std::vector<size_t>;
using Buffer = std::vector<char>;

enum NdCopyFlag{
    RowMajorBigEndian,
    RowMajorSmallEndian,
    ColumnMajorBigEndian,
    ColumnMajorLittleEndian,
};

template<typename T>
int NdCopy(
        const Buffer &in_buffer,
        const Dims &in_start,
        const Dims &in_count,
        NdCopyFlag in_flag,
        Buffer &out_buffer,
        const Dims &out_start,
        const Dims &out_count,
        NdCopyFlag out_flag
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
