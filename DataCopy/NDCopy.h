#include <vector>

using Dims = std::vector<size_t>;
using Buffer = std::vector<char>;

enum NdCopyFlag{
    RowMajorBigEndian,
    RowMajorSmallEndian,
    ColumnMajorBigEndian,
    ColumnMajorSmallEndian,
};

template<class T>
int NdCopy(const Buffer &in_buffer, const Dims &in_start, Dims &in_count, NdCopyFlag in_flag,
           Buffer &out_buffer, const Dims &out_start, Dims &out_count, NdCopyFlag out_flag);

template<class T>
int NdCopy2(const Buffer &in_buffer, const Dims &in_start, Dims &in_count, NdCopyFlag in_flag,
           Buffer &out_buffer, const Dims &out_start, Dims &out_count, NdCopyFlag out_flag);


