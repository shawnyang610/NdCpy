# n-dimension-data-copy

## Description
 * Developed as a new optimized data copying core to replace the existing one used
 * in ornladios/ADIOS2 which runs on super computers for data I/O
## Features
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
 
 
