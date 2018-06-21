//
//  NDCopy.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/20/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include "NDCopy.h"
using namespace std;

/*helper functions*/
void getInputEnd(vector<size_t>&input_end, const vector<size_t>&input_start,
                 const vector<size_t>&input_count);
void getOutputEnd(vector<size_t>& output_end, const vector<size_t>&output_start,
                  const vector<size_t>&output_count);
void getOverlapStart(vector<size_t>&overlap_start,const vector<size_t>&input_start,
                     const vector<size_t>&output_start);
void getOverlapEnd(vector<size_t>&overlap_end, vector<size_t>&input_end,
                   vector<size_t>&output_end);
void getOverlapCount(vector<size_t>&overlap_count,vector<size_t>&overlap_start,
                     vector<size_t>&overlap_end);
bool hasOverlap(vector<size_t>&overlap_start, vector<size_t>overlap_end);
//size_t getPosition (vector<size_t>& start_pos, vector<size_t>& overlap_start);
size_t getMinContDimn(const vector<size_t>&input_start,
                      const vector<size_t>&input_count,
                      vector<size_t>&overlap_start, vector<size_t>overlap_count);
size_t getBlockSize(vector<size_t> overlap_count, size_t min_cont_dim,
                    size_t elm_size);
void copy_cat(size_t cur_dim,char*& input_overlap_base,
              char*& output_overlap_base,vector<size_t>in_ovlp_gap_size,
              vector<size_t>out_ovlp_gap_size,
              vector<size_t>&overlap_count,size_t min_cont_dim,size_t block_size);
void getIOStrides(vector<size_t>&io_stride,vector<size_t>&io_count,
                 size_t elm_size);
void getIO_OverlapBase(char*& IO_overlap_base, Buffer io,
                                const vector<size_t>& io_start,
                                vector<size_t>io_elm_counts,
                                vector<size_t>& overlap_start);
void getIO_OvlpGapSize(vector<size_t>&io_ovlp_gap_size,vector<size_t>&in_stride,
                       vector<size_t>&input_count,vector<size_t>&overlap_count);

/*end of helper functions*/


/*calling function*/
template<class T>
int NdCopy(const Buffer &input, const Dims &input_start, Dims &input_count, NdCopyFlag input_flag,
           Buffer &output, const Dims &output_start, Dims &output_count, NdCopyFlag output_flag)
{
    vector<size_t>input_end(input_start.size());
    vector<size_t>output_end(input_start.size());
    vector<size_t>overlap_start(input_start.size());
    vector<size_t>overlap_end(input_start.size());
    vector<size_t>overlap_count(input_start.size());
    vector<size_t>in_stride(input_start.size());
    vector<size_t>out_stride(input_start.size());
    vector<size_t>in_ovlp_gap_size(input_start.size());
    vector<size_t>out_ovlp_gap_size(input_start.size());
    size_t min_cont_dim, block_size;
    char* input_overlap_base=NULL;
    char* output_overlap_base=NULL;
    //main algm starts here
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start,output_count);
    getOverlapStart(overlap_start,input_start, output_start);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getIOStrides(in_stride,input_count,sizeof(T));
    getIOStrides(out_stride, output_count,sizeof(T));
    getIO_OvlpGapSize(in_ovlp_gap_size,in_stride,input_count, overlap_count);
    getIO_OvlpGapSize(out_ovlp_gap_size,out_stride,output_count, overlap_count);
    getIO_OverlapBase(input_overlap_base,input,input_start,in_stride,
                      overlap_start);
    getIO_OverlapBase(output_overlap_base,output,output_start,out_stride,
                      overlap_start);
    min_cont_dim=getMinContDimn(input_start, input_count,overlap_start,
                   overlap_count);
    block_size=getBlockSize(overlap_count, min_cont_dim, sizeof(T));
    copy_cat(0,input_overlap_base,output_overlap_base,in_ovlp_gap_size,
             out_ovlp_gap_size,overlap_count,min_cont_dim,block_size);
    //end of main algm
    return 0;
}
/*end of calling function*/


/* helper function definitions*/
void getInputEnd(vector<size_t>& input_end, const vector<size_t>&input_start,
                 const vector<size_t>&input_count){
    for (size_t i=0;i<input_start.size();i++)
        input_end[i]=input_start[i]+input_count[i]-1;
}
void getOutputEnd(vector<size_t>& output_end, const vector<size_t>&output_start,
                  const vector<size_t>&output_count){
    for (size_t i=0; i<output_start.size();i++)
        output_end[i]=output_start[i]+output_count[i]-1;
}
void getOverlapStart(vector<size_t>&overlap_start,const vector<size_t>&input_start,
                     const vector<size_t>&output_start){
    for (size_t i=0;i<overlap_start.size();i++)
        overlap_start[i]=input_start[i]>output_start[i]?input_start[i]:output_start[i];
}
void getOverlapEnd(vector<size_t>&overlap_end, vector<size_t>&input_end,
                   vector<size_t>&output_end){
    for (size_t i=0; i<overlap_end.size();i++)
        overlap_end[i]=input_end[i]<output_end[i]?input_end[i]:output_end[i];
}
void getOverlapCount(vector<size_t>&overlap_count,vector<size_t>&overlap_start,
                     vector<size_t>&overlap_end){
    for (size_t i=0; i<overlap_count.size(); i++)
        overlap_count[i]=overlap_end[i]-overlap_start[i]+1;
}
bool hasOverlap(vector<size_t>&overlap_start, vector<size_t>overlap_end){
    for (size_t i=0; i<overlap_start.size();i++)
        if (overlap_end[i]<overlap_start[i]) return false;
    return true;
}
void getIOStrides(vector<size_t>&io_stride,vector<size_t>&io_count,
                 size_t elm_size){
    //io_elm_count[i] holds total number of elements under each element
    //of the i'th dimension
    io_stride[io_stride.size()-1]=elm_size;
    if (io_stride.size()>1)
        io_stride[io_stride.size()-2]=io_count[io_stride.size()-1]*elm_size;
    if (io_stride.size()>2) {
        size_t i=io_stride.size()-3;
        while (true){
            io_stride[i]=io_count[i+1]*io_stride[i+1];
            if (i==0) break;
            else i--;
        }
    }
}
void getIO_OverlapBase(char*& IO_overlap_base, Buffer io,
                       const vector<size_t>& io_start,
                       vector<size_t>io_elm_counts,
                       vector<size_t>& overlap_start){
    IO_overlap_base = io.data();
    for (size_t i=0; i<io_start.size();i++)
        IO_overlap_base+=(overlap_start[i]-io_start[i])*io_elm_counts[i];
}
void getIO_OvlpGapSize(vector<size_t>&io_ovlp_gap_size,vector<size_t>&in_stride,
                       vector<size_t>&input_count,vector<size_t>&overlap_count){
    for (size_t i=0;i<io_ovlp_gap_size.size();i++)
        io_ovlp_gap_size[i]=(input_count[i]-overlap_count[i])*in_stride[i];
}
size_t getMinContDimn(const vector<size_t>&input_start, const vector<size_t>&input_count,
                      vector<size_t>&overlap_start,vector<size_t>overlap_count){
    //    note: min_cont_dim is the first index where its input box and overlap box
    //    are not fully match. therefore all data below this branch is continous
    //    and this determins the Biggest continuous block size - Each element of the
    //    current dimension.
    size_t i=input_start.size()-1;
    while (true){
        if (i==0) break;
        if (input_count[i]!=overlap_count[i]) break;
        i--;
    }
    return i;
}
size_t getBlockSize(vector<size_t> overlap_count, size_t min_cont_dim,
                    size_t elm_size){
    size_t res=elm_size;
    for (size_t i=min_cont_dim; i<overlap_count.size();i++)
        res*=overlap_count[i];
    return res;
}

/*
 worst case: this recursive function calls itself for W times.
 W = overlap_area_count[0]*overlap_area_count[1]*...*overlap_area_count[total_dims-1]
 best case: B =1, happens when min_cont_dim is the top dimension
 the computational overhead for copying each block = two additions,
 and another two additions everytime it backtracks up a dimension.
 */
void copy_cat(size_t cur_dim,char*& input_overlap_base,
              char*& output_overlap_base,vector<size_t>in_ovlp_gap_size,
              vector<size_t>out_ovlp_gap_size,
              vector<size_t>&overlap_count,size_t min_cont_dim,size_t block_size)
{
    //note: all elements in and below this node is continuous on input
    //copy the continous data block
    if (cur_dim==min_cont_dim){
        memcpy(output_overlap_base, input_overlap_base,
               block_size);
        input_overlap_base+=block_size;
        output_overlap_base+=block_size;
    }
    //recursively call itself in order, for every element current node has
    //on a deeper level, stops upon reaching min_cont_dim
    if (cur_dim<min_cont_dim)
        for (size_t i=0; i<overlap_count[i];i++)
            copy_cat(cur_dim+1, input_overlap_base, output_overlap_base,
                     in_ovlp_gap_size,out_ovlp_gap_size,overlap_count,
                     min_cont_dim,block_size);
    //the gap between current node and the next needs to be padded so that
    //next continous blocks starts at the correct position for both input and output
    //the size of the gap depends on the depth in dimensions,level backtracked and
    //the difference in element counts between the Input/output and overlap area.
    input_overlap_base+=in_ovlp_gap_size[cur_dim];
    output_overlap_base+=out_ovlp_gap_size[cur_dim];
}
/*end of helper function definitions*/
