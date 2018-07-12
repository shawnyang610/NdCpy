//
//  NDCopy.cpp
//  DataCopy
//  Created by Shawn Yang on 6/20/18.
//  shawnyang610@gmail.com

#include <algorithm>
#include <cstring>
#include "NDCopy.h"
#include <functional>

//RecurCoreFastMode(): helper function
//Copys n-dimensional Data from input to output in the same major and endianess.
//It looks for the largest contiguous data block size in the overlap (by its helper
//functions) and copies to the output buffer in blocks. the memory address
//calculation complexity for copying each block is minimized to O(1), which is
//independent of the number of dimensions.
static void RecurCoreFastMode(size_t curDim,
                    char*& inOvlpBase,
                    char*& outOvlpBase,
                    Dims& inOvlpGapSize,
                    Dims& outOvlpGapSize,
                    Dims& ovlpCount,
                    size_t minContDim,
                    size_t blockSize)
{
  //note: all elements in and below this node are contiguous on input and output
  //copy the contiguous data block
  if (curDim==minContDim){
    std::memcpy(outOvlpBase, inOvlpBase, blockSize);
    inOvlpBase+=blockSize;
    outOvlpBase+=blockSize;
  }
  //recursively call itself in order, for every element current node has
  //on a deeper level, stops upon reaching minContDim
  //case: curDim<minCountDim
  else{
    for (size_t i=0; i<ovlpCount[curDim];i++)
      RecurCoreFastMode(curDim+1,
              inOvlpBase,
              outOvlpBase,
              inOvlpGapSize,
              outOvlpGapSize,
              ovlpCount,
              minContDim,
              blockSize);
  }
  //the gap between current node and the next needs to be padded so that
  //next contigous block starts at the correct position for both input and output
  //the size of the gap depends on the depth in dimensions,level backtracked and
  //the difference in element counts between the Input/output and overlap area.
  inOvlpBase+=inOvlpGapSize[curDim];
  outOvlpBase+=outOvlpGapSize[curDim];
}

//RecurCoreRevEndianMode(): helper function
//Copys n-dimensional Data from input to output in the same major but in reversed
//endianess. the memory address calculation complexity for copying each element is
//minimized to average O(1), which is independent of the number of dimensions.
static void RecurCoreRevEndianMode(size_t curDim,
                          char*& inOvlpBase,
                          char*& outOvlpBase,
                          Dims& inOvlpGapSize,
                          Dims& outOvlpGapSize,
                          Dims& ovlpCount,
                          size_t minCountDim,
                          size_t blockSize,
                          size_t elmSize,
                          size_t numElmsPerBlock)
{
  if (curDim==minCountDim){
    //only the following block is different from the original copyCat
    //each byte of each element in the continuous block needs
    //to be copied in reverse order
    for(size_t i=0;i<numElmsPerBlock;i++){
      for(size_t j=0; j<elmSize;j++){
        outOvlpBase[j]=inOvlpBase[elmSize-1-j];
      }
      inOvlpBase+=elmSize;
      outOvlpBase+=elmSize;
    }
  }
  //case: curDim<minCountDim
  else{
    for (size_t i=0; i<ovlpCount[curDim];i++)
      RecurCoreRevEndianMode(curDim+1,inOvlpBase,outOvlpBase,inOvlpGapSize,
                    outOvlpGapSize,ovlpCount,minCountDim,blockSize,elmSize,
                    numElmsPerBlock);
  }
  inOvlpBase+=inOvlpGapSize[curDim];
  outOvlpBase+=outOvlpGapSize[curDim];
}

//RecurCoreRevMajorMode(): helper function
//Copys n-dimensional Data from input to output in the same Endianess but in the
//reversed major. the memory address calculation complexity for copying each element is
//minimized to average O(1), which is independent of the number of dimensions.
static void RecurCoreRevMajorMode(size_t curDim,
                           char* inBase,
                           char* outBase,
                           Dims& inRltvOvlpSPos,
                           Dims& outRltvOvlpSPos,
                           Dims& inStride,
                           Dims& outStride,
                           Dims& ovlpCount,
                           size_t elmSize)
{
  if (curDim==inStride.size()){
    std::memcpy(outBase, inBase, elmSize);
  }
  else {
    for (size_t i=0; i<ovlpCount[curDim];i++)
      RecurCoreRevMajorMode(curDim+1,
                     inBase+(inRltvOvlpSPos[curDim]+i)*inStride[curDim],
                     outBase+(outRltvOvlpSPos[curDim]+i)*outStride[curDim],
                     inRltvOvlpSPos,outRltvOvlpSPos,
                     inStride,outStride,
                     ovlpCount,
                     elmSize);
  }
}

//RecurCoreRevMajorRevEndianMode(): helper function
//Copys n-dimensional Data from input to output in the reversed Endianess and Major.
//The memory address calculation complexity for copying each element is
//minimized to average O(1), which is independent of the number of dimensions.
static void RecurCoreRevMajorRevEndianMode(size_t curDim,
                                 char* inBase,
                                 char* outBase,
                                 Dims& inRltvOvlpSPos,
                                 Dims& outRltvOvlpSPos,
                                 Dims& inStride,
                                 Dims& outStride,
                                 Dims& ovlpCount,
                                 size_t elmSize)
{
  if (curDim==inStride.size()){
    //the following for-loop block is the only difference from the original
    //flippedCopyCat
    for (size_t i=0; i<elmSize;i++){
      //memcpy(outBase+i, inBase+elmSize-1-i, 1);
      outBase[i]=inBase[elmSize-1-i];
    }
  }
  else {
    for (size_t i=0; i<ovlpCount[curDim];i++)
      RecurCoreRevMajorRevEndianMode(curDim+1,
                           inBase+(inRltvOvlpSPos[curDim]+i)*inStride[curDim],
                           outBase+(outRltvOvlpSPos[curDim]+i)*outStride[curDim],
                           inRltvOvlpSPos,
                           outRltvOvlpSPos,
                           inStride,
                           outStride,
                           ovlpCount,
                           elmSize);
  }
}


//NdCopy():
//Copys n-dimensional Data from an input buffer of any major and endianess to an
//output buffer of any major and endianess.
//Copying between same major and endianess yields the best speed. The optimization
//is achived by, first: looks for the largest contiguous data block size and copies in
//big chunks. second: By using dynamic, depth-first traversal, The overhead for memory
//address calculation for each block copied is reduced from O(n) to average O(1).
//which means copying speed is drastically improved for data of higher dimensions.
//For copying between buffers of diffenrent majors or endianesses, only optimization
//for memory address calculation is applied.
template<class T>
int NdCopy(const Buffer& in,
           const Dims& inStart,
           const Dims& inCount,
           NdCopyFlag inFlag,
           Buffer& out,
           const Dims& outStart,
           const Dims& outCount,
           NdCopyFlag outFlag)
{
  Dims inEnd(inStart.size());
  Dims outEnd(inStart.size());
  Dims ovlpStart(inStart.size());
  Dims ovlpEnd(inStart.size());
  Dims ovlpCount(inStart.size());
  Dims inStride(inStart.size());
  Dims outStride(inStart.size());
  Dims inOvlpGapSize(inStart.size());
  Dims outOvlpGapSize(inStart.size());
  Dims inRltvOvlpStartPos(inStart.size());
  Dims outRltvOvlpStartPos(inStart.size());
  size_t minContDim, blockSize;
  char* inOvlpBase=nullptr;
  char* outOvlpBase=nullptr;
  bool isSameMajor= inFlag.isRowMajor==outFlag.isRowMajor? true:false;
  bool isSameEndian= inFlag.isBigEndian==outFlag.isBigEndian? true:false;
  auto GetInEnd=[](Dims& inEnd, const Dims&inStart,const Dims&inCount)
  {
    for (size_t i=0;i<inStart.size();i++)
      inEnd[i]=inStart[i]+inCount[i]-1;
  };
  auto GetOutEnd=[](Dims& outEnd,
                    const Dims& outStart,
                    const Dims&output_count)
  {
    for (size_t i=0; i<outStart.size();i++)
      outEnd[i]=outStart[i]+output_count[i]-1;
  };
  auto GetOvlpStart=[](Dims& ovlpStart,
                       const Dims& inStart,
                       const Dims& outStart)
  {
    for (size_t i=0;i<ovlpStart.size();i++)
      ovlpStart[i]=inStart[i]>outStart[i]?inStart[i]:outStart[i];
  };
  auto GetOvlpEnd=[](Dims& ovlpEnd, Dims& inEnd, Dims& outEnd)
  {
    for (size_t i=0; i<ovlpEnd.size();i++)
      ovlpEnd[i]=inEnd[i]<outEnd[i]?inEnd[i]:outEnd[i];
  };
  auto GetOvlpCount=[](Dims& ovlpCount, Dims& ovlpStart, Dims& ovlpEnd)
  {
    for (size_t i=0; i<ovlpCount.size(); i++)
      ovlpCount[i]=ovlpEnd[i]-ovlpStart[i]+1;
  };
  auto HasOvlp=[](Dims& ovlpStart, Dims& ovlpEnd)
  {
    for (size_t i=0; i<ovlpStart.size();i++)
      if (ovlpEnd[i]<ovlpStart[i]) return false;
    return true;
  };

  auto GetIoStrides=[](Dims& ioStride,const Dims& ioCount, size_t elmSize)
  {
    //ioStride[i] holds the total number of elements under each element
    //of the i'th dimension
    ioStride[ioStride.size()-1]=elmSize;
    if (ioStride.size()>1)
      ioStride[ioStride.size()-2]=ioCount[ioStride.size()-1]*elmSize;
    if (ioStride.size()>2) {
      size_t i=ioStride.size()-3;
      while (true){
        ioStride[i]=ioCount[i+1]*ioStride[i+1];
        if (i==0) break;
        else i--;
      }
    }
  };

  auto GetIoOvlpBase=[](char *& ioOvlpBase, const Buffer& io,
                        const Dims& ioStart,
                        Dims & ioStride,
                        Dims & ovlpStart)
  {
    ioOvlpBase =(char*)io.data();
    for (size_t i=0; i<ioStart.size();i++)
      ioOvlpBase=ioOvlpBase+(ovlpStart[i]-ioStart[i])*ioStride[i];
  };
  auto GetIoOvlpGapSize=[](Dims& ioOvlpGapSize,
                           Dims& ioStride,
                           const Dims& ioCount,
                           Dims& ovlpCount)
  {
    for (size_t i=0;i<ioOvlpGapSize.size();i++)
      ioOvlpGapSize[i]=(ioCount[i]-ovlpCount[i])*ioStride[i];
  };
  auto GetMinContDim=[](const Dims& inCount,
                        const Dims outCount,
                        Dims& ovlpCount)
  {
    //    note: minContDim is the first index where its input box and overlap box
    //    are not fully match. therefore all data below this branch is continous
    //    and this determins the Biggest continuous block size - Each element of the
    //    current dimension.
    size_t i=ovlpCount.size()-1;
    while (true){
      if (i==0) break;
      if ((inCount[i]!=ovlpCount[i]) ||
          (outCount[i]!=ovlpCount[i])) break;
      i--;
    }
    return i;
  };
  auto GetBlockSize=[](Dims& ovlpCount, size_t minContDim, size_t elmSize)
  {
    size_t res=elmSize;
    for (size_t i=minContDim; i<ovlpCount.size();i++)
      res*=ovlpCount[i];
    return res;
  };

  auto GetRltvOvlpStartPos=[](Dims& ioRltvOvlpStart,
                              const Dims& ioStart,
                              Dims& ovlpStart)
  {
    for(size_t i=0;i<ioStart.size();i++)
      ioRltvOvlpStart[i]=ovlpStart[i]-ioStart[i];
  };

  //main flow
  if (isSameMajor){
    GetInEnd(inEnd,inStart,inCount);
    GetOutEnd(outEnd,outStart,outCount);
    GetOvlpStart(ovlpStart,inStart, outStart);
    GetOvlpEnd(ovlpEnd, inEnd, outEnd);
    GetOvlpCount(ovlpCount,ovlpStart, ovlpEnd);
    if (!HasOvlp(ovlpStart, ovlpEnd)) return 1;//no overlap found
    GetIoStrides(inStride,inCount,sizeof(T));
    GetIoStrides(outStride, outCount,sizeof(T));
    GetIoOvlpGapSize(inOvlpGapSize,inStride,inCount, ovlpCount);
    GetIoOvlpGapSize(outOvlpGapSize,outStride,outCount, ovlpCount);
    GetIoOvlpBase(inOvlpBase,in,inStart,inStride,
                  ovlpStart);
    GetIoOvlpBase(outOvlpBase,out,outStart,outStride,ovlpStart);
    minContDim=GetMinContDim(inCount,outCount,ovlpCount);
    blockSize=GetBlockSize(ovlpCount, minContDim, sizeof(T));
    if(isSameEndian){
      //Quick Copying Mode:Same Major, Same Endian"
      RecurCoreFastMode(0,
              inOvlpBase,
              outOvlpBase,
              inOvlpGapSize,
              outOvlpGapSize,
              ovlpCount,
              minContDim,
              blockSize);
    }
    else{
      //same major, dif. endian mode
      RecurCoreRevEndianMode(0,
                    inOvlpBase,
                    outOvlpBase,
                    inOvlpGapSize,
                    outOvlpGapSize,
                    ovlpCount,
                    minContDim,
                    blockSize,
                    sizeof(T),
                    blockSize/sizeof(T));
    }
  }
  else {
    //avg computational overhead is O(1) for each intersecting byte
    //worst case can be O(n) where n is number of dimensions
    Dims revOutStart(outStart);
    Dims revOutCount(outCount);
    std::reverse(revOutStart.begin(), revOutStart.end());
    std::reverse(revOutCount.begin(), revOutCount.end());
    GetInEnd(inEnd,inStart,inCount);
    GetOutEnd(outEnd,revOutStart,revOutCount);
    GetOvlpStart(ovlpStart,inStart, revOutStart);
    GetOvlpEnd(ovlpEnd, inEnd, outEnd);
    GetOvlpCount(ovlpCount,ovlpStart, ovlpEnd);
    if (!HasOvlp(ovlpStart, ovlpEnd)) return 1;//no overlap found
    GetIoStrides(inStride,inCount,sizeof(T));
    GetIoStrides(outStride, outCount, sizeof(T));
    std::reverse(outStride.begin(), outStride.end());
    GetRltvOvlpStartPos(inRltvOvlpStartPos,inStart,ovlpStart);
    GetRltvOvlpStartPos(outRltvOvlpStartPos, revOutStart,ovlpStart);
    inOvlpBase=(char*)in.data();
    outOvlpBase=(char*)out.data();
    if (isSameEndian){
      //Copy Mode: Different Major, Same Endian"
      RecurCoreRevMajorMode(0,
                     inOvlpBase,
                     outOvlpBase,
                     inRltvOvlpStartPos,
                     outRltvOvlpStartPos,
                     inStride,
                     outStride,
                     ovlpCount,
                     sizeof(T));
    }
    else{
      //dif major, dif. endian mode
      RecurCoreRevMajorRevEndianMode(0,
                           inOvlpBase,
                           outOvlpBase,
                           inRltvOvlpStartPos,
                           outRltvOvlpStartPos,
                           inStride,
                           outStride,
                           ovlpCount,
                           sizeof(T));
    }
  }
  return 0;
}
