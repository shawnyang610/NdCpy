# ndcopy -  A highly optimized I/O utility library written in C++ that moves n-dimension, non-contiguous data array between contiguous memories, supports any hardware architectures (Big - Little Endian conversion), and Fortran memory scheme (Row - Column Major conversion)


## Features
 * Highly optimized specificly for high dimensional, contiguous data copying.
 * Copies n-dimensional Data from a source buffer to destination buffer, either
 * can be of any Major and Endianess. Return 1 if no overlap is found.
 * Copying from row-major to row-major of the same endian yields the best speed.
 * Optimizations: first: looks for the largest contiguous data
 * block size and copies the block of data as a whole. Second: by using
 * depth-first traversal for dynamic memory pointer calculations.
 * address calculation for each copied block is reduced to O(1) from O(n).
 * which means the computational cost is drastically reduced for data of higher
 * dimensions.
 * For copying involving column major, or different endianess only the
 * second optimization is applied.
 
## Use case
 * Used as the new "dataman" core function for data copying to replace the old one used
 * in ornladios/ADIOS2, a library runs on super computers for data I/O

## Work in progress
 * All functional features are in production.
 * Convert the library to a professional, maintenable structure. (Work in progress)
 * Adding tests for all functionalities. (Work in progress)
 * Create a website for demo and algorithm analysis for possible further improvements. (Work in progress)
