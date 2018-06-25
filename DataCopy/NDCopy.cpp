//
//  NDCopy.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/20/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include <cstring>
#include <algorithm>
#include "NDCopy.h"
using namespace std;

/*helper functions*/
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
void flipOutputCount(vector<size_t> &output_count){
  reverse(output_count.begin(), output_count.end());
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
void getIoOverlapBase(char*& IO_overlap_base, Buffer io,
                      const vector<size_t>& io_start,
                      vector<size_t>io_elm_counts,
                      vector<size_t>& overlap_start){
  IO_overlap_base = io.data();
  for (size_t i=0; i<io_start.size();i++)
    IO_overlap_base+=(overlap_start[i]-io_start[i])*io_elm_counts[i];
}
void getIoOvlpGapSize(vector<size_t>&io_ovlp_gap_size,vector<size_t>&in_stride,
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
void copyCat(size_t cur_dim,char*& input_overlap_base,
             char*& output_overlap_base,vector<size_t>in_ovlp_gap_size,
             vector<size_t>out_ovlp_gap_size,
             vector<size_t>&overlap_count,size_t min_cont_dim,size_t block_size)
{
  //note: all elements in and below this node is continuous on input
  //copy the continous data block
  if (cur_dim==min_cont_dim){
      std::memcpy(output_overlap_base, input_overlap_base,
           block_size);
    input_overlap_base+=block_size;
    output_overlap_base+=block_size;
  }
  //recursively call itself in order, for every element current node has
  //on a deeper level, stops upon reaching min_cont_dim
  if (cur_dim<min_cont_dim)
    for (size_t i=0; i<overlap_count[i];i++)
      copyCat(cur_dim+1, input_overlap_base, output_overlap_base,
              in_ovlp_gap_size,out_ovlp_gap_size,overlap_count,
              min_cont_dim,block_size);
  //the gap between current node and the next needs to be padded so that
  //next continous blocks starts at the correct position for both input and output
  //the size of the gap depends on the depth in dimensions,level backtracked and
  //the difference in element counts between the Input/output and overlap area.
  input_overlap_base+=in_ovlp_gap_size[cur_dim];
  output_overlap_base+=out_ovlp_gap_size[cur_dim];
}

string copyMode(NdCopyFlag input_flag, NdCopyFlag output_flag){
  if ((input_flag==RowMajorBigEndian && output_flag==RowMajorBigEndian)||
      (input_flag==ColumnMajorBigEndian && output_flag==ColumnMajorBigEndian)||
      (input_flag==RowMajorSmallEndian && output_flag==RowMajorSmallEndian)||
      (input_flag==ColumnMajorLittleEndian && output_flag==ColumnMajorLittleEndian))
  {
    return "same_maj_same_endian";
  }
  else if ((input_flag==RowMajorBigEndian && output_flag==RowMajorSmallEndian)||
           (input_flag==RowMajorSmallEndian&&output_flag==RowMajorBigEndian)||
           (input_flag==ColumnMajorBigEndian && output_flag==ColumnMajorLittleEndian)
           ||(input_flag==ColumnMajorLittleEndian &&
              output_flag==ColumnMajorBigEndian))
  {
    return "same_maj_dif_endian";
  }
  else if ((input_flag==RowMajorBigEndian && output_flag==ColumnMajorBigEndian)||
           (input_flag==ColumnMajorBigEndian&&output_flag==RowMajorBigEndian)||
           (input_flag==RowMajorSmallEndian && output_flag==ColumnMajorLittleEndian)
           ||(input_flag==ColumnMajorLittleEndian &&
              output_flag==RowMajorSmallEndian))
  {
    return "dif_maj_same_endian";
  }
  else
    return "dif_maj_dif_endian";
}
void getRelativeOvlpHeadPos(vector<size_t>&ioRelOvlpStart,
                            const vector<size_t>&ioStart,
                            vector<size_t>&ovlpStart){
  for(size_t i=0;i<ioStart.size();i++)
    ioRelOvlpStart[i]=ovlpStart[i]-ioStart[i];
}
void getFlippedIOStrides(vector<size_t>&io_stride,
                         vector<size_t>&io_count,size_t elmSize){
  io_count.push_back(elmSize); //add one extra "byte dimension" to count
  reverse(io_count.begin(), io_count.end()); //reverse the order
  //also add one extra "byte dim" to stride
  //stride of the byte dim is always 1 byte
  io_stride.push_back(1);
  if (io_stride.size()>1) {
    size_t i=io_stride.size()-2;
    while (true){
      io_stride[i]=io_count[i+1]*io_stride[i+1];
      if (i==0) break;
      else i--;
    }
  }
  //reverse stride to the same order as input for easy access
  reverse(io_stride.begin(), io_stride.end());
}
void flipCopyByByte(size_t curDim, char* inBase,char* outBase,
                    vector<size_t>&InRelOvlpHeadPos,
                    vector<size_t>&outRelOvlpHeadPos,vector<size_t>&in_stride,
                    vector<size_t>&out_stride,vector<size_t>&overlap_count,
                    size_t elmSize){
  //Note: in_stride and out_stride have 1 extra dimension-the "byte" dim.
  //index in_stride.size()-1 is the last dim which is the byte dim.
  if (curDim==in_stride.size()-1){
    for (size_t i=0; i<elmSize;i++){
      memcpy(outBase+i*out_stride[curDim], inBase+i*in_stride[curDim], 1);
    }
  }
  else {
    for (size_t i=0; i<overlap_count[i];i++){
      flipCopyByByte(curDim+1,
                     inBase+(InRelOvlpHeadPos[curDim]+i)*in_stride[curDim],
                     outBase+(outRelOvlpHeadPos[curDim]+i)*out_stride[curDim],
                     InRelOvlpHeadPos,outRelOvlpHeadPos, in_stride,out_stride,
                     overlap_count, elmSize);
    }
  }
  //todo
}
/*end of helper functions*/


/*calling function*/
template<class T>
int NdCopy(const Buffer &input, const Dims &input_start, Dims &input_count,
           NdCopyFlag input_flag,Buffer &output, const Dims &output_start,
           Dims &output_count, NdCopyFlag output_flag){
  vector<size_t>input_end(input_start.size());
  vector<size_t>output_end(input_start.size());
  vector<size_t>overlap_start(input_start.size());
  vector<size_t>overlap_end(input_start.size());
  vector<size_t>overlap_count(input_start.size());
  vector<size_t>in_stride(input_start.size());
  vector<size_t>out_stride(input_start.size());
  vector<size_t>in_ovlp_gap_size(input_start.size());
  vector<size_t>out_ovlp_gap_size(input_start.size());
  vector<size_t> InRelOvlpHeadPos(input_start.size());
  vector<size_t> outRelOvlpHeadPos(input_start.size());
  size_t min_cont_dim, block_size;
  char* input_overlap_base=nullptr;
  char* output_overlap_base=nullptr;


  if (copyMode(input_flag,output_flag)=="same_maj_same_endian"){
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start,output_count);
    getOverlapStart(overlap_start,input_start, output_start);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getIOStrides(in_stride,input_count,sizeof(T));
    getIOStrides(out_stride, output_count,sizeof(T));
    getIoOvlpGapSize(in_ovlp_gap_size,in_stride,input_count, overlap_count);
    getIoOvlpGapSize(out_ovlp_gap_size,out_stride,output_count, overlap_count);
    getIoOverlapBase(input_overlap_base,input,input_start,in_stride,
                     overlap_start);
    getIoOverlapBase(output_overlap_base,output,output_start,out_stride,
                     overlap_start);
    min_cont_dim=getMinContDimn(input_start, input_count,overlap_start,
                                overlap_count);
    block_size=getBlockSize(overlap_count, min_cont_dim, sizeof(T));
    copyCat(0,input_overlap_base,output_overlap_base,in_ovlp_gap_size,
            out_ovlp_gap_size,overlap_count,min_cont_dim,block_size);
  }
  else if (copyMode(input_flag,output_flag)=="same_maj_dif_endian"){

  }
  else if (copyMode(input_flag,output_flag)=="dif_maj_same_endian"){
    //avg computational overhead is O(1) for each intersecting byte
    //worst case can be O(n) where n is number of dimensions
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start,output_count);
    getOverlapStart(overlap_start,input_start, output_start);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getIOStrides(in_stride,input_count,sizeof(T));
    in_stride.push_back(1); //adds an additional "byte" dimension
    getFlippedIOStrides(out_stride, output_count,sizeof(T));
    getRelativeOvlpHeadPos(InRelOvlpHeadPos,input_start,overlap_start);
    getRelativeOvlpHeadPos(outRelOvlpHeadPos, output_start,overlap_start);
    input_overlap_base=(char*)input.data();
    output_overlap_base=output.data();
    flipCopyByByte(0,input_overlap_base,output_overlap_base,InRelOvlpHeadPos,
                   outRelOvlpHeadPos,in_stride,out_stride,overlap_count,sizeof(T));
  }

  //diff_maj_dif_endian
  else {

  }
  //end of main algm
  return 0;
}
/*end of calling function*/

