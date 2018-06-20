//
//  NDCopy.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/20/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include<cstring>
using namespace std;

/*helper functions*/
void getInputEnd(vector<size_t>&input_end, vector<size_t>&input_start,vector<size_t>&input_count);
void getOutputEnd(vector<size_t>& output_end, vector<size_t>&output_start,vector<size_t>&output_count);
void getOverlapStart(vector<size_t>&overlap_start,vector<size_t>&input_start,vector<size_t>&output_start);
void getOverlapEnd(vector<size_t>&overlap_end, vector<size_t>&input_end,vector<size_t>&output_end);
void getOverlapCount(vector<size_t>&overlap_count,vector<size_t>&overlap_start,vector<size_t>&overlap_end);
bool hasOverlap(vector<size_t>&overlap_start, vector<size_t>overlap_end);
size_t getPosition (vector<size_t>& start_pos, vector<size_t>& overlap_start);
size_t getMinContDimn(vector<size_t>&input_start, vector<size_t>&input_count,
                      vector<size_t>&overlap_start, vector<size_t>overlap_count);
size_t getBlockSize(vector<size_t> overlap_count, size_t min_cont_dim);
void copy_cat(size_t lv, vector<char>&input,vector<char>&output,size_t& in_start_offset,
              size_t& out_start_offset,vector<size_t>in_elem_count,vector<size_t>out_elem_count,
              vector<size_t>&overlap_count,size_t min_cont_dim, size_t block_size, size_t elm_size);
void getElmCount(vector<size_t>&io_elm_count,vector<size_t>&io_count);
size_t getIOStartOffset(vector<size_t>& io_start,vector<size_t>io_elm_counts,
                        vector<size_t>& overlap_start);

/*end of helper functions*/


/*calling function*/
int NdCopy(vector<char> input, vector<size_t> input_start, vector<size_t> input_count,
           vector<char> output, vector<size_t> output_start, vector<size_t> output_count,
           size_t element_size, bool dimension_reversed=false)
{
    vector<size_t>input_end(input_start.size());
    vector<size_t>output_end(input_start.size());
    vector<size_t>overlap_start(input_start.size());
    vector<size_t>overlap_end(input_start.size());
    vector<size_t>overlap_count(input_start.size());
    vector<size_t>input_elm_count(input_start.size());
    vector<size_t>output_elm_count(input_start.size());
    size_t min_cont_dim, block_size, input_overlap_start_offset,
    output_overlap_start_offset;
    //main algm starts here
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start,output_count);
    getOverlapStart(overlap_start,input_start, output_start);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getElmCount(input_elm_count,input_count);
    getElmCount(output_elm_count, output_count);
    input_overlap_start_offset= getIOStartOffset(input_start,input_elm_count,overlap_start);
    output_overlap_start_offset=getIOStartOffset(output_start,output_elm_count,overlap_start);
    min_cont_dim = getMinContDimn(input_start, input_count,overlap_start, overlap_count);
    block_size= getBlockSize(overlap_count, min_cont_dim);
    copy_cat(0,input,output,input_overlap_start_offset, output_overlap_start_offset,
             input_elm_count, output_elm_count,overlap_count,min_cont_dim,block_size, element_size);
    //end of main algm
    return 0;
}
/*end of calling function*/


/* helper function definitions*/
void getInputEnd(vector<size_t>& input_end, vector<size_t>&input_start, vector<size_t>&input_count){
    for (size_t i=0;i<input_start.size();i++)
        input_end[i]=input_start[i]+input_count[i]-1;
};
void getOutputEnd(vector<size_t>& output_end, vector<size_t>&output_start, vector<size_t>&output_count){
    for (size_t i=0; i<output_start.size();i++)
        output_end[i]=output_start[i]+output_count[i]-1;
};
void getOverlapStart(vector<size_t>&overlap_start,vector<size_t>&input_start, vector<size_t>&output_start){
    for (size_t i=0;i<overlap_start.size();i++)
        overlap_start[i]=input_start[i]>output_start[i]?input_start[i]:output_start[i];
}
void getOverlapEnd(vector<size_t>&overlap_end, vector<size_t>&input_end, vector<size_t>&output_end){
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
void getElmCount(vector<size_t>&io_elm_count,vector<size_t>&io_count){
    //io_elm_count[i] holds total number of elements under each element of the i th dimension
    io_elm_count[io_elm_count.size()-1]=1;
    if (io_elm_count.size()>1) io_elm_count[io_elm_count.size()-2]=io_count[io_elm_count.size()-1];
    if (io_elm_count.size()>2) {
        size_t i=io_elm_count.size()-3;
        while (true){
            io_elm_count[i]=io_count[i+1]*io_elm_count[i+1];
            if (i==0) break;
            else i--;
        }
    }
}
size_t getIOStartOffset(vector<size_t>& io_start,vector<size_t>io_elm_counts, vector<size_t>& overlap_start){
    size_t res=1;
    for (size_t i=0; i<io_start.size();i++)
        res=res+(overlap_start[i]-io_start[i])*io_elm_counts[i];
    return res;
}
size_t getMinContDimn(vector<size_t>&input_start, vector<size_t>&input_count,
                      vector<size_t>&overlap_start,vector<size_t>overlap_count){
    //    note: min_cont_dim is the first index where its input box and overlap box
    //    are not fully match. therefore all data below this branch is continous
    //    and this determins the Biggest continuous block size - Each element of the current dimension.
    size_t i=input_start.size()-1;
    while (true){
        if (i==0) break;
        if (input_count[i]!=overlap_count[i]) break;
        i--;
    }
    return i;
}
size_t getBlockSize(vector<size_t> overlap_count, size_t min_cont_dim){
    size_t res=1;
    for (size_t i=min_cont_dim; i<overlap_count.size();i++)
        res*=overlap_count[i];
    return res;
}

/*
 worst case: this recursive function calls itself for W times.
 W = overlap_area_count[0]*overlap_area_count[1]*...*overlap_area_count[total_dims-1]
 best case: B =1, happens when min_cont_dim is the top dimension
 */
void copy_cat(size_t cur_dim, vector<char>&input,vector<char>&output,size_t& in_start_offset,
              size_t& out_start_offset,vector<size_t>in_elem_count,vector<size_t>out_elem_count,
              vector<size_t>&overlap_count,size_t min_cont_dim, size_t block_size, size_t elm_size){
    //note: all elements in and below this node is continuous on input
    //copy the continous data block
    if (cur_dim==min_cont_dim){
        memcpy(&output+out_start_offset*elm_size, &input+in_start_offset*elm_size, block_size*elm_size);
        in_start_offset+=block_size;
        out_start_offset+=block_size;
    }
    //recursively call itself in order, for every element current node has
    //on a deeper level, stops upon reaching min_cont_dim
    if (cur_dim<min_cont_dim)
        for (size_t i=0; i<overlap_count[i];i++)
            copy_cat(cur_dim+1, input, output, in_start_offset, out_start_offset,in_elem_count,
                     out_elem_count,overlap_count, min_cont_dim, block_size, elm_size);
    //the gap between current node and the next needs to be padded so that
    //next continous blocks starts at the correct position for both input and output
    //the size of the gap depends on the depth in dimensions,level backtracked and
    //the difference in element counts between the Input/output and overlap area.
    in_start_offset+=in_elem_count[cur_dim];
    out_start_offset+=out_elem_count[cur_dim];
}
/*end of helper function definitions*/
