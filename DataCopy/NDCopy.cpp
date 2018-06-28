//
//  NDCopy.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/20/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//
#include <algorithm>
#include <cstring>
#include "NDCopy.h"

/*helper functions*/
static void getInputEnd(Dims& input_end, const Dims&input_start,
                 const Dims&input_count){
  for (size_t i=0;i<input_start.size();i++)
    input_end[i]=input_start[i]+input_count[i]-1;
}
static void getOutputEnd(Dims& output_end, const Dims&output_start,
                  const Dims&output_count){
  for (size_t i=0; i<output_start.size();i++)
    output_end[i]=output_start[i]+output_count[i]-1;
}
static void getOverlapStart(Dims&overlap_start,const Dims&input_start,
                     const Dims&output_start){
  for (size_t i=0;i<overlap_start.size();i++)
    overlap_start[i]=input_start[i]>output_start[i]?input_start[i]:output_start[i];
}
static void getOverlapEnd(Dims&overlap_end, Dims&input_end,
                   Dims&output_end){
  for (size_t i=0; i<overlap_end.size();i++)
    overlap_end[i]=input_end[i]<output_end[i]?input_end[i]:output_end[i];
}
static void getOverlapCount(Dims&overlap_count,Dims&overlap_start,
                     Dims&overlap_end){
  for (size_t i=0; i<overlap_count.size(); i++)
    overlap_count[i]=overlap_end[i]-overlap_start[i]+1;
}
static bool hasOverlap(Dims&ovlpHead, Dims&ovlpTail){
  for (size_t i=0; i<ovlpHead.size();i++)
    if (ovlpTail[i]<ovlpHead[i]) return false;
  return true;
}
static void flipOutputCount(Dims&output_count){
  reverse(output_count.begin(), output_count.end());
}
static void getIOStrides(Dims&io_stride,const Dims&io_count,
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

static void getIoOverlapBase(char *& IO_overlap_base, const Buffer& io,
                      const Dims& io_start,
                      Dims &io_stride,
                      Dims & overlap_start){
  IO_overlap_base =(char*)io.data();
  for (size_t i=0; i<io_start.size();i++)
    IO_overlap_base=IO_overlap_base+(overlap_start[i]-io_start[i])*io_stride[i];
}
static void getIoOvlpGapSize(Dims&io_ovlp_gap_size,Dims&in_stride,
                      const Dims&input_count,Dims&overlap_count){
  for (size_t i=0;i<io_ovlp_gap_size.size();i++)
    io_ovlp_gap_size[i]=(input_count[i]-overlap_count[i])*in_stride[i];
}
static size_t getMinContDimn(const Dims& input_count,
                            const Dims output_count,
                             Dims& overlap_count)
{
  //    note: min_cont_dim is the first index where its input box and overlap box
  //    are not fully match. therefore all data below this branch is continous
  //    and this determins the Biggest continuous block size - Each element of the
  //    current dimension.
  size_t i=overlap_count.size()-1;
  while (true){
    if (i==0) break;
    if ((input_count[i]!=overlap_count[i]) ||
        (output_count[i]!=overlap_count[i])) break;
    i--;
  }
  return i;
}
static size_t getBlockSize(Dims& overlap_count, size_t min_cont_dim,
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
static void copyCat(size_t cur_dim,char*& input_overlap_base,
             char*& output_overlap_base,Dims& in_ovlp_gap_size,
             Dims& out_ovlp_gap_size,
             Dims& overlap_count,size_t min_cont_dim,size_t block_size)
{
  //note: all elements in and below this node is continuous on input
  //copy the continous data block
  if (cur_dim==min_cont_dim){
    std::memcpy(output_overlap_base, input_overlap_base, block_size);
    input_overlap_base+=block_size;
    output_overlap_base+=block_size;
  }
  //recursively call itself in order, for every element current node has
  //on a deeper level, stops upon reaching min_cont_dim
  if (cur_dim<min_cont_dim)
    for (size_t i=0; i<overlap_count[cur_dim];i++)
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
static void endianCopyCat(size_t cur_dim,char*& input_overlap_base,
                    char*& output_overlap_base,Dims& in_ovlp_gap_size,
                    Dims& out_ovlp_gap_size,
                    Dims& overlap_count,size_t min_cont_dim,size_t block_size,
                          size_t elmSize, size_t numElmPerBlock)
{
  if (cur_dim==min_cont_dim){
    //only the following block is different from the original copyCat
    //each all bytes of each element in the continuous block needs
    //to be copied in reverse order
    for(size_t i=0;i<numElmPerBlock;i++){
      for(size_t j=0; j<elmSize;j++){
        output_overlap_base[0]=input_overlap_base[elmSize-1-j];
      }
      input_overlap_base+=elmSize;
      output_overlap_base+=elmSize;
    }
  }
  
  if (cur_dim<min_cont_dim)
    for (size_t i=0; i<overlap_count[cur_dim];i++)
      copyCat(cur_dim+1, input_overlap_base, output_overlap_base,
              in_ovlp_gap_size,out_ovlp_gap_size,overlap_count,
              min_cont_dim,block_size);
  input_overlap_base+=in_ovlp_gap_size[cur_dim];
  output_overlap_base+=out_ovlp_gap_size[cur_dim];
}

static std::string copyMode(NdCopyFlag& input_flag, NdCopyFlag& output_flag){
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
static void getRelativeOvlpHeadPos(Dims& ioRelOvlpStart,
                            const Dims& ioStart,
                            Dims& ovlpStart){
  for(size_t i=0;i<ioStart.size();i++)
    ioRelOvlpStart[i]=ovlpStart[i]-ioStart[i];
}
static void getFlippedIOStrides(Dims& io_stride,
                         const Dims& io_count,size_t elmSize){
  Dims flippedCount(io_count);
  flippedCount.push_back(elmSize);   //add one extra "byte dimension" to count
  reverse(flippedCount.begin(), flippedCount.end()); //reverse the order
  //also add one extra "byte dim" to stride
  //stride of the byte dim is always 1 byte
  io_stride.push_back(1);
  if (io_stride.size()>1) {
    size_t i=io_stride.size()-2;
    while (true){
      io_stride[i]=flippedCount[i+1]*io_stride[i+1];
      if (i==0) break;
      else i--;
    }
  }
  //reverse stride to the same order as input for easy access
  reverse(io_stride.begin(), io_stride.end());
}
static void flippedCopyCat(size_t curDim, char* inBase,char* outBase,
                    Dims& InRelOvlpHeadPos,
                    Dims& outRelOvlpHeadPos,Dims& in_stride,
                    Dims& out_stride,Dims& overlap_count,
                    size_t elmSize){
  if (curDim==in_stride.size()){
    memcpy(outBase, inBase, elmSize);
  }
  else {
    for (size_t i=0; i<overlap_count[curDim];i++){
      flippedCopyCat(curDim+1,
                     inBase+(InRelOvlpHeadPos[curDim]+i)*in_stride[curDim],
                     outBase+(outRelOvlpHeadPos[curDim]+i)*out_stride[curDim],
                     InRelOvlpHeadPos,outRelOvlpHeadPos, in_stride,out_stride,
                     overlap_count, elmSize);
    }
  }
}
static void flippedEndianCopyCat(size_t curDim, char* inBase,char* outBase,
                           Dims& InRelOvlpHeadPos,
                           Dims& outRelOvlpHeadPos,Dims& in_stride,
                           Dims& out_stride,Dims& overlap_count,
                           size_t elmSize){
  if (curDim==in_stride.size()){
    //the following for-loop block is the only difference from the original
    //flippedCopyCat
    for (size_t i=0; i<elmSize;i++){
      outBase[i]=inBase[elmSize-1-i];
    }
    
  }
  else {
    for (size_t i=0; i<overlap_count[curDim];i++){
      flippedCopyCat(curDim+1,
                     inBase+(InRelOvlpHeadPos[curDim]+i)*in_stride[curDim],
                     outBase+(outRelOvlpHeadPos[curDim]+i)*out_stride[curDim],
                     InRelOvlpHeadPos,outRelOvlpHeadPos, in_stride,out_stride,
                     overlap_count, elmSize);
    }
  }
}
/*end of helper functions*/


/*calling function*/
template<class T>
int NdCopy(const Buffer &input, const Dims &input_start, const Dims &input_count,
           NdCopyFlag input_flag,Buffer &output, const Dims &output_start,
           const Dims &output_count, NdCopyFlag output_flag){
  Dims input_end(input_start.size());
  Dims output_end(input_start.size());
  Dims overlap_start(input_start.size());
  Dims overlap_end(input_start.size());
  Dims overlap_count(input_start.size());
  Dims in_stride(input_start.size());
  Dims out_stride(input_start.size());
  Dims in_ovlp_gap_size(input_start.size());
  Dims out_ovlp_gap_size(input_start.size());
  Dims InRelOvlpHeadPos(input_start.size());
  Dims outRelOvlpHeadPos(input_start.size());
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

    min_cont_dim=getMinContDimn(input_count,output_count,overlap_count);
    block_size=getBlockSize(overlap_count, min_cont_dim, sizeof(T));
    copyCat(0,input_overlap_base,output_overlap_base,in_ovlp_gap_size,
            out_ovlp_gap_size,overlap_count,min_cont_dim,block_size);
  }
  else if (copyMode(input_flag,output_flag)=="same_maj_dif_endian"){
    std::cout<<"NDCopy copyMode selection:same_maj_dif_endian "<<std::endl;
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
    
    min_cont_dim=getMinContDimn(input_count,output_count,overlap_count);
    block_size=getBlockSize(overlap_count, min_cont_dim, sizeof(T));
    endianCopyCat(0,input_overlap_base,output_overlap_base,in_ovlp_gap_size,
            out_ovlp_gap_size,overlap_count,min_cont_dim,block_size,sizeof(T),
                  block_size/sizeof(T));
  }
  else if (copyMode(input_flag,output_flag)=="dif_maj_same_endian"){
    std::cout<<"NDCopy copyMode selection:dif_maj_same_endian "<<std::endl;

    //avg computational overhead is O(1) for each intersecting byte
    //worst case can be O(n) where n is number of dimensions
    Dims output_start_rev(output_start);
    Dims output_count_rev(output_count);
    std::reverse(output_start_rev.begin(), output_start_rev.end());
    std::reverse(output_count_rev.begin(), output_count_rev.end());
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start_rev,output_count_rev);
    getOverlapStart(overlap_start,input_start, output_start_rev);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getIOStrides(in_stride,input_count,sizeof(T));
    getIOStrides(out_stride, output_count, sizeof(T));
    std::reverse(out_stride.begin(), out_stride.end());
    getRelativeOvlpHeadPos(InRelOvlpHeadPos,input_start,overlap_start);
    getRelativeOvlpHeadPos(outRelOvlpHeadPos, output_start_rev,overlap_start);
    input_overlap_base=(char*)input.data();
    output_overlap_base=output.data();
    flippedCopyCat(0,input_overlap_base,output_overlap_base,InRelOvlpHeadPos,
                   outRelOvlpHeadPos,in_stride,out_stride,overlap_count,sizeof(T));
  }

  //diff_maj_dif_endian
  else {
    std::cout<<"NDCopy copyMode selection:dif_maj_dif_endian "<<std::endl;
    Dims output_start_rev(output_start);
    Dims output_count_rev(output_count);
    std::reverse(output_start_rev.begin(), output_start_rev.end());
    std::reverse(output_count_rev.begin(), output_count_rev.end());
    getInputEnd(input_end,input_start,input_count);
    getOutputEnd(output_end,output_start_rev,output_count_rev);
    getOverlapStart(overlap_start,input_start, output_start_rev);
    getOverlapEnd(overlap_end, input_end, output_end);
    getOverlapCount(overlap_count,overlap_start, overlap_end);
    if (!hasOverlap(overlap_start, overlap_end)) return 1;//no overlap found
    getIOStrides(in_stride,input_count,sizeof(T));
    getIOStrides(out_stride, output_count, sizeof(T));
    std::reverse(out_stride.begin(), out_stride.end());
    getRelativeOvlpHeadPos(InRelOvlpHeadPos,input_start,overlap_start);
    getRelativeOvlpHeadPos(outRelOvlpHeadPos, output_start_rev,overlap_start);
    input_overlap_base=(char*)input.data();
    output_overlap_base=output.data();
    flippedEndianCopyCat(0,input_overlap_base,output_overlap_base,InRelOvlpHeadPos,
                   outRelOvlpHeadPos,in_stride,out_stride,overlap_count,sizeof(T));
  }
  return 0;
}
/*end of calling function*/

