#include <vector>

using Dims = std::vector<size_t>;
using Buffer = std::vector<char>;

int NdCopy(const Buffer &input, const Dims &input_start, Dims &input_count,
           Buffer &output, const Dims &output_start, Dims &output_count,
           size_t element_size, bool dimension_reversed=false);

int NdCopy2(const Buffer &input, const Dims &input_start, Dims &input_count,
           Buffer &output, const Dims &output_start, Dims &output_count,
           size_t element_size, bool dimension_reversed=false);


