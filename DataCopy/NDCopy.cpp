//
//  NDCopy.cpp
//  DataCopy
//  Created by Shawn Yang on 6/20/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//
#include <algorithm>
#include <cstring>
#include "NDCopy.h"

/*helper functions*/
static void getInEnd(Dims& inEnd, const Dims&inStart,const Dims&inCount)
{
  for (size_t i=0;i<inStart.size();i++)
    inEnd[i]=inStart[i]+inCount[i]-1;
}
static void getOutEnd(Dims& outEnd,
                         const Dims& outStart,
                         const Dims&output_count)
{
  for (size_t i=0; i<outStart.size();i++)
    outEnd[i]=outStart[i]+output_count[i]-1;
}
static void getOvlpStart(Dims& ovlpStart,
                            const Dims& inStart,
                            const Dims& outStart)
{
  for (size_t i=0;i<ovlpStart.size();i++)
    ovlpStart[i]=inStart[i]>outStart[i]?inStart[i]:outStart[i];
}
static void getOvlpEnd(Dims& ovlpEnd, Dims& inEnd, Dims& outEnd)
{
  for (size_t i=0; i<ovlpEnd.size();i++)
    ovlpEnd[i]=inEnd[i]<outEnd[i]?inEnd[i]:outEnd[i];
}
static void getOvlpCount(Dims& ovlpCount, Dims& ovlpStart, Dims& ovlpEnd)
{
  for (size_t i=0; i<ovlpCount.size(); i++)
    ovlpCount[i]=ovlpEnd[i]-ovlpStart[i]+1;
}
static bool hasOvlp(Dims& ovlpStart, Dims& ovlpEnd)
{
  for (size_t i=0; i<ovlpStart.size();i++)
    if (ovlpEnd[i]<ovlpStart[i]) return false;
  return true;
}

static void getIoStrides(Dims& ioStride,const Dims& ioCount, size_t elmSize)
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
}

static void getIoOvlpBase(char *& ioOvlpBase, const Buffer& io,
                      const Dims& ioStart,
                      Dims & ioStride,
                      Dims & ovlpStart)
{
  ioOvlpBase =(char*)io.data();
  for (size_t i=0; i<ioStart.size();i++)
    ioOvlpBase=ioOvlpBase+(ovlpStart[i]-ioStart[i])*ioStride[i];
}
static void getIoOvlpGapSize(Dims& ioOvlpGapSize,
                             Dims& ioStride,
                             const Dims& ioCount,
                             Dims& ovlpCount)
{
  for (size_t i=0;i<ioOvlpGapSize.size();i++)
    ioOvlpGapSize[i]=(ioCount[i]-ovlpCount[i])*ioStride[i];
}
static size_t getMinContDim(const Dims& inCount,
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
}
static size_t getBlockSize(Dims& ovlpCount, size_t minContDim, size_t elmSize)
{
  size_t res=elmSize;
  for (size_t i=minContDim; i<ovlpCount.size();i++)
    res*=ovlpCount[i];
  return res;
}


/*
 worst case: this recursive function calls itself for W times.
 W = overlap_area_count[0]*overlap_area_count[1]*...*overlap_area_count[total_dims-1]
 best case: B =1, happens when min_cont_dim is the top dimension
 the computational overhead for copying each block = two additions,
 and another two additions everytime it backtracks up a dimension.
 */
static void copyCat(size_t curDim,
                    char*& inOvlpBase,
                    char*& outOvlpBase,
                    Dims& inOvlpGapSize,
                    Dims& outOvlpGapSize,
                    Dims& ovlpCount,
                    size_t minContDim,
                    size_t blockSize)
{
  //note: all elements in and below this node are continuous on input and output
  //copy the continous data block
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
      copyCat(curDim+1,
              inOvlpBase,
              outOvlpBase,
              inOvlpGapSize,
              outOvlpGapSize,
              ovlpCount,
              minContDim,
              blockSize);
  }
  //the gap between current node and the next needs to be padded so that
  //next continous blocks starts at the correct position for both input and output
  //the size of the gap depends on the depth in dimensions,level backtracked and
  //the difference in element counts between the Input/output and overlap area.
  inOvlpBase+=inOvlpGapSize[curDim];
  outOvlpBase+=outOvlpGapSize[curDim];
}
static void iterativeCopyCat (
                       char*& inOvlpBase,
                       char*& outOvlpBase,
                       Dims& inOvlpGapSize,
                       Dims& outOvlpGapSize,
                       Dims& ovlpCount,
                       size_t minContDim,
                       size_t blockSize)
{
  Dims pos(ovlpCount.size(),0);
  size_t curDim=0;
  while (true){
    while(curDim!=minContDim){
      pos[curDim]++;
      curDim++;
    }
    std::memcpy(outOvlpBase, inOvlpBase, blockSize);
    inOvlpBase+=blockSize;
    outOvlpBase+=blockSize;
    do {
      if (curDim==0) //more logical but expensive place for the check
        return;
      inOvlpBase+=inOvlpGapSize[curDim];
      outOvlpBase+=outOvlpGapSize[curDim];
      pos[curDim]=0;
      curDim--;
    } while (pos[curDim]==ovlpCount[curDim]);
  }
}
static void endianCopyCat(size_t curDim,
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
    //each all bytes of each element in the continuous block needs
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
      endianCopyCat(curDim+1,inOvlpBase,outOvlpBase,inOvlpGapSize,
                    outOvlpGapSize,ovlpCount,minCountDim,blockSize,elmSize,
                    numElmsPerBlock);
  }
  inOvlpBase+=inOvlpGapSize[curDim];
  outOvlpBase+=outOvlpGapSize[curDim];
}

static void iterativeEndianCopyCat (
                              char*& inOvlpBase,
                              char*& outOvlpBase,
                              Dims& inOvlpGapSize,
                              Dims& outOvlpGapSize,
                              Dims& ovlpCount,
                              size_t minContDim,
                              size_t blockSize,
                              size_t elmSize,
                              size_t numElmsPerBlock)
{
  Dims pos(ovlpCount.size(),0);
  size_t curDim=0;
  while (true){
    while(curDim!=minContDim){
      pos[curDim]++;
      curDim++;
    }
    for(size_t i=0;i<numElmsPerBlock;i++){
      for(size_t j=0; j<elmSize;j++){
        outOvlpBase[j]=inOvlpBase[elmSize-1-j];
      }
      inOvlpBase+=elmSize;
      outOvlpBase+=elmSize;
    }
    do {
      if (curDim==0) //more logical but expensive place for the check
        return;
      inOvlpBase+=inOvlpGapSize[curDim];
      outOvlpBase+=outOvlpGapSize[curDim];
      pos[curDim]=0;
      curDim--;
    } while (pos[curDim]==ovlpCount[curDim]);
  }
}

static void getRltvOvlpStartPos(Dims& ioRltvOvlpStart,
                            const Dims& ioStart,
                            Dims& ovlpStart)
{
  for(size_t i=0;i<ioStart.size();i++)
    ioRltvOvlpStart[i]=ovlpStart[i]-ioStart[i];
}

static void flippedCopyCat(size_t curDim,
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
      flippedCopyCat(curDim+1,
                     inBase+(inRltvOvlpSPos[curDim]+i)*inStride[curDim],
                     outBase+(outRltvOvlpSPos[curDim]+i)*outStride[curDim],
                     inRltvOvlpSPos,outRltvOvlpSPos,
                     inStride,outStride,
                     ovlpCount,
                     elmSize);
  }
}
//performance is 50% slower than the recursive version
static void iterativeFlippedCopyCat(
                           char* inBase,
                           char* outBase,
                           Dims& inRltvOvlpSPos,
                           Dims& outRltvOvlpSPos,
                           Dims& inStride,
                           Dims& outStride,
                           Dims& ovlpCount,
                           size_t elmSize)
{
  size_t curDim=0;
  Dims pos(ovlpCount.size()+1,0);
  std::vector<char*>inAddr(ovlpCount.size()+1);
  inAddr[0]=inBase;
  std::vector<char*>outAddr(ovlpCount.size()+1);
  outAddr[0]=outBase;
  while (true){
    while (curDim!=inStride.size()){
      inAddr[curDim+1]=
      inAddr[curDim]+(inRltvOvlpSPos[curDim]+pos[curDim])*inStride[curDim];
      outAddr[curDim+1]=
      outAddr[curDim]+(outRltvOvlpSPos[curDim]+pos[curDim])*outStride[curDim];
      pos[curDim]++;
      curDim++;
    }
      std::memcpy(outAddr[curDim], inAddr[curDim], elmSize);
    do {
      if (curDim==0)
        return;
      pos[curDim]=0;
      curDim--;
    }while(pos[curDim]==ovlpCount[curDim]);
  }
}

static void flippedEndianCopyCat(size_t curDim,
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
      flippedEndianCopyCat(curDim+1,
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
//performance is 50% slower than the recursive version
static void iterativeFlippedEndianCopyCat(
                                    char* inBase,
                                    char* outBase,
                                    Dims& inRltvOvlpSPos,
                                    Dims& outRltvOvlpSPos,
                                    Dims& inStride,
                                    Dims& outStride,
                                    Dims& ovlpCount,
                                    size_t elmSize)
{
  size_t curDim=0;
  Dims pos(ovlpCount.size()+1,0);
  std::vector<char*>inAddr(ovlpCount.size()+1);
  inAddr[0]=inBase;
  std::vector<char*>outAddr(ovlpCount.size()+1);
  outAddr[0]=outBase;
  while (true){
    while (curDim!=inStride.size()){
      inAddr[curDim+1]=
      inAddr[curDim]+(inRltvOvlpSPos[curDim]+pos[curDim])*inStride[curDim];
      outAddr[curDim+1]=
      outAddr[curDim]+(outRltvOvlpSPos[curDim]+pos[curDim])*outStride[curDim];
      pos[curDim]++;
      curDim++;
    }
    for (size_t i=0; i<elmSize;i++){
      //memcpy(outBase+i, inBase+elmSize-1-i, 1);
      outAddr[curDim][i]=inAddr[curDim][elmSize-1-i];
    }
    do {
      if (curDim==0)
        return;
      pos[curDim]=0;
      curDim--;
    }while(pos[curDim]==ovlpCount[curDim]);
  }
}
class CopyMode{
public:
  bool isSameMajor, isSameEndian;
  CopyMode(NdCopyFlag& in, NdCopyFlag& out){
    isSameMajor= in.isRowMajor==out.isRowMajor? true:false;
    isSameEndian= in.isBigEndian==out.isBigEndian? true:false;
  }
};
/*end of helper functions*/


/*calling function*/
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
  CopyMode copyMode(inFlag, outFlag);
  if (copyMode.isSameMajor){
    getInEnd(inEnd,inStart,inCount);
    getOutEnd(outEnd,outStart,outCount);
    getOvlpStart(ovlpStart,inStart, outStart);
    getOvlpEnd(ovlpEnd, inEnd, outEnd);
    getOvlpCount(ovlpCount,ovlpStart, ovlpEnd);
    if (!hasOvlp(ovlpStart, ovlpEnd)) return 1;//no overlap found
    getIoStrides(inStride,inCount,sizeof(T));
    getIoStrides(outStride, outCount,sizeof(T));
    getIoOvlpGapSize(inOvlpGapSize,inStride,inCount, ovlpCount);
    getIoOvlpGapSize(outOvlpGapSize,outStride,outCount, ovlpCount);
    getIoOvlpBase(inOvlpBase,in,inStart,inStride,
                     ovlpStart);
    getIoOvlpBase(outOvlpBase,out,outStart,outStride,ovlpStart);
    minContDim=getMinContDim(inCount,outCount,ovlpCount);
    blockSize=getBlockSize(ovlpCount, minContDim, sizeof(T));
    if(copyMode.isSameEndian){
//      std::cout<<"Quick Copying Mode:Same Major, Same Endian"<<std::endl;
//      copyCat(0,
//              inOvlpBase,
//              outOvlpBase,
//              inOvlpGapSize,
//              outOvlpGapSize,
//              ovlpCount,
//              minContDim,
//              blockSize);
      
      iterativeCopyCat (inOvlpBase,
                        outOvlpBase,
                        inOvlpGapSize,
                        outOvlpGapSize,
                        ovlpCount,
                        minContDim,
                        blockSize);
    }
    else{//same major, dif. endian mode
//      std::cout<<"Copying Mode:Same Major, Dif. Endian"<<std::endl;
//      endianCopyCat(0,
//                    inOvlpBase,
//                    outOvlpBase,
//                    inOvlpGapSize,
//                    outOvlpGapSize,
//                    ovlpCount,
//                    minContDim,
//                    blockSize,
//                    sizeof(T),
//                    blockSize/sizeof(T));
      iterativeEndianCopyCat (
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
    getInEnd(inEnd,inStart,inCount);
    getOutEnd(outEnd,revOutStart,revOutCount);
    getOvlpStart(ovlpStart,inStart, revOutStart);
    getOvlpEnd(ovlpEnd, inEnd, outEnd);
    getOvlpCount(ovlpCount,ovlpStart, ovlpEnd);
    if (!hasOvlp(ovlpStart, ovlpEnd)) return 1;//no overlap found
    getIoStrides(inStride,inCount,sizeof(T));
    getIoStrides(outStride, outCount, sizeof(T));
    std::reverse(outStride.begin(), outStride.end());
    getRltvOvlpStartPos(inRltvOvlpStartPos,inStart,ovlpStart);
    getRltvOvlpStartPos(outRltvOvlpStartPos, revOutStart,ovlpStart);
    inOvlpBase=(char*)in.data();
    outOvlpBase=(char*)out.data();
    if (copyMode.isSameEndian){
//      std::cout<<"Copy Mode: Different Major, Same Endian"<<std::endl;
//      flippedCopyCat(0,
//                     inOvlpBase,
//                     outOvlpBase,
//                     inRltvOvlpStartPos,
//                     outRltvOvlpStartPos,
//                     inStride,
//                     outStride,
//                     ovlpCount,
//                     sizeof(T));
      
      iterativeFlippedCopyCat(
                     inOvlpBase,
                     outOvlpBase,
                     inRltvOvlpStartPos,
                     outRltvOvlpStartPos,
                     inStride,
                     outStride,
                     ovlpCount,
                     sizeof(T));
    }
    else{ //dif major, dif. endian mode
//      std::cout<<"Copy Mode: Different Major, Different Endian"<<std::endl;
//      flippedEndianCopyCat(0,
//                           inOvlpBase,
//                           outOvlpBase,
//                           inRltvOvlpStartPos,
//                           outRltvOvlpStartPos,
//                           inStride,
//                           outStride,
//                           ovlpCount,
//                           sizeof(T));
      iterativeFlippedEndianCopyCat(
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
/*end of calling function*/




