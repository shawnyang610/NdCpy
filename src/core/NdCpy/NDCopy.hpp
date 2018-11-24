//
//  NDCopy.hpp
//  src
//  Created by Shawn Yang on 6/20/18.
//  shawnyang610@gmail.com

#ifndef NDCOPY_HPP
#define NDCOPY_HPP

#include <algorithm>
#include <cstring>
//#include "NDCopy.h"
#include <vector>
#include <functional>

using Dims = std::vector<size_t>;
using Buffer = std::vector<char>;

template <class T>
int NdCopy(const char *in, const Dims &inStart, const Dims &inCount,
           const bool inIsRowMajor, const bool inIsLittleEndian, char *out,
           const Dims &outStart, const Dims &outCount, const bool outIsRowMajor,
           const bool outIsLittleEndian, const Dims &inMemStart = Dims(),
           const Dims &inMemCount = Dims(), const Dims &outMemStart = Dims(),
           const Dims &outMemCount = Dims(), const bool safeMode = false);

//***************Start of NdCopy() and its 8 helpers ***************
// Author:Shawn Yang, shawnyang610@gmail.com
//
// NdCopyRecurDFSeqPadding(): helper function
// Copys n-dimensional Data from input to output in row major and
// same endianess.
// It looks for the largest contiguous data block size in the overlap (by its
// helper
// functions) and copies to the output buffer in blocks. the memory address
// calculation complexity for copying each block is minimized to O(1), which is
// independent of the number of dimensions.
static void NdCopyRecurDFSeqPadding(size_t curDim, const char *&inOvlpBase,
                                    char *&outOvlpBase, Dims &inOvlpGapSize,
                                    Dims &outOvlpGapSize, Dims &ovlpCount,
                                    size_t &minContDim, size_t &blockSize)
{
  // note: all elements in and below this node are contiguous on input and
  // output
  // copy the contiguous data block
  if (curDim == minContDim)
  {
    std::memcpy(outOvlpBase, inOvlpBase, blockSize);
    inOvlpBase += blockSize + inOvlpGapSize[curDim];
    outOvlpBase += blockSize + outOvlpGapSize[curDim];
  }
  // recursively call itself in order, for every element current node has
  // on a deeper level, stops upon reaching minContDim
  // case: curDim<minCountDim
  else
  {
    for (size_t i = 0; i < ovlpCount[curDim]; i++)
    {
      NdCopyRecurDFSeqPadding(curDim + 1, inOvlpBase, outOvlpBase,
                              inOvlpGapSize, outOvlpGapSize, ovlpCount,
                              minContDim, blockSize);
    }
    // the gap between current node and the next needs to be padded so that
    // next contigous block starts at the correct position for both input
    // and output
    // the size of the gap depends on the depth in dimensions,level
    // backtracked and
    // the difference in element counts between the Input/output and overlap
    // area.
    inOvlpBase += inOvlpGapSize[curDim];
    outOvlpBase += outOvlpGapSize[curDim];
  }
}

// NdCopyRecurDFSeqPaddingRevEndian(): helper function
// Copys n-dimensional Data from input to output in the row major but in
// reversed endianess. the memory address calculation complexity for copying
// each element is minimized to average O(1), which is independent of
// the number of dimensions.

static void
NdCopyRecurDFSeqPaddingRevEndian(size_t curDim, const char *&inOvlpBase,
                                 char *&outOvlpBase, Dims &inOvlpGapSize,
                                 Dims &outOvlpGapSize, Dims &ovlpCount,
                                 size_t minCountDim, size_t blockSize,
                                 size_t elmSize, size_t numElmsPerBlock)
{
  if (curDim == minCountDim)
  {
    // each byte of each element in the continuous block needs
    // to be copied in reverse order
    for (size_t i = 0; i < numElmsPerBlock; i++)
    {
      for (size_t j = 0; j < elmSize; j++)
      {
        outOvlpBase[j] = inOvlpBase[elmSize - 1 - j];
      }
      inOvlpBase += elmSize;
      outOvlpBase += elmSize;
    }
  }
  // case: curDim<minCountDim
  else
  {
    for (size_t i = 0; i < ovlpCount[curDim]; i++)
      NdCopyRecurDFSeqPaddingRevEndian(
                                       curDim + 1, inOvlpBase, outOvlpBase, inOvlpGapSize,
                                       outOvlpGapSize, ovlpCount, minCountDim, blockSize, elmSize,
                                       numElmsPerBlock);
  }
  inOvlpBase += inOvlpGapSize[curDim];
  outOvlpBase += outOvlpGapSize[curDim];
}

// NdCopyRecurDFNonSeqDynamic(): helper function
// Copys n-dimensional Data from input to output in the same Endianess
// used for buffer of Column major
// the memory address calculation complexity for copying each element is
// minimized to average O(1), which is independent of the number of dimensions.
static void NdCopyRecurDFNonSeqDynamic(size_t curDim, const char *inBase,
                                       char *outBase, Dims &inRltvOvlpSPos,
                                       Dims &outRltvOvlpSPos, Dims &inStride,
                                       Dims &outStride, Dims &ovlpCount,
                                       size_t elmSize)
{
  if (curDim == inStride.size())
  {
    std::memcpy(outBase, inBase, elmSize);
  }
  else
  {
    for (size_t i = 0; i < ovlpCount[curDim]; i++)
      NdCopyRecurDFNonSeqDynamic(
                                 curDim + 1,
                                 inBase + (inRltvOvlpSPos[curDim] + i) * inStride[curDim],
                                 outBase + (outRltvOvlpSPos[curDim] + i) * outStride[curDim],
                                 inRltvOvlpSPos, outRltvOvlpSPos, inStride, outStride, ovlpCount,
                                 elmSize);
  }
}

// NdCopyRecurDFNonSeqDynamicRevEndian(): helper function
// Copys n-dimensional Data from input to output in the reversed Endianess and
// Major.
// The memory address calculation complexity for copying each element is
// minimized to average O(1), which is independent of the number of dimensions.

static void NdCopyRecurDFNonSeqDynamicRevEndian(
                                                size_t curDim, const char *inBase, char *outBase, Dims &inRltvOvlpSPos,
                                                Dims &outRltvOvlpSPos, Dims &inStride, Dims &outStride, Dims &ovlpCount,
                                                size_t elmSize)
{
  if (curDim == inStride.size())
  {
    for (size_t i = 0; i < elmSize; i++)
    {
      outBase[i] = inBase[elmSize - 1 - i];
    }
  }
  else
  {
    for (size_t i = 0; i < ovlpCount[curDim]; i++)
      NdCopyRecurDFNonSeqDynamicRevEndian(
                                          curDim + 1,
                                          inBase + (inRltvOvlpSPos[curDim] + i) * inStride[curDim],
                                          outBase + (outRltvOvlpSPos[curDim] + i) * outStride[curDim],
                                          inRltvOvlpSPos, outRltvOvlpSPos, inStride, outStride, ovlpCount,
                                          elmSize);
  }
}

static void NdCopyIterDFSeqPadding(const char *&inOvlpBase, char *&outOvlpBase,
                                   Dims &inOvlpGapSize, Dims &outOvlpGapSize,
                                   Dims &ovlpCount, size_t minContDim,
                                   size_t blockSize)
{
  Dims pos(ovlpCount.size(), 0);
  size_t curDim = 0;
  while (true)
  {
    while (curDim != minContDim)
    {
      pos[curDim]++;
      curDim++;
    }
    std::memcpy(outOvlpBase, inOvlpBase, blockSize);
    inOvlpBase += blockSize;
    outOvlpBase += blockSize;
    do
    {
      if (curDim == 0)
        return;
      inOvlpBase += inOvlpGapSize[curDim];
      outOvlpBase += outOvlpGapSize[curDim];
      pos[curDim] = 0;
      curDim--;
    } while (pos[curDim] == ovlpCount[curDim]);
  }
}

static void NdCopyIterDFSeqPaddingRevEndian(
                                            const char *&inOvlpBase, char *&outOvlpBase, Dims &inOvlpGapSize,
                                            Dims &outOvlpGapSize, Dims &ovlpCount, size_t minContDim, size_t blockSize,
                                            size_t elmSize, size_t numElmsPerBlock)
{
  Dims pos(ovlpCount.size(), 0);
  size_t curDim = 0;
  while (true)
  {
    while (curDim != minContDim)
    {
      pos[curDim]++;
      curDim++;
    }
    for (size_t i = 0; i < numElmsPerBlock; i++)
    {
      for (size_t j = 0; j < elmSize; j++)
      {
        outOvlpBase[j] = inOvlpBase[elmSize - 1 - j];
      }
      inOvlpBase += elmSize;
      outOvlpBase += elmSize;
    }
    do
    {
      if (curDim == 0)
        return;
      inOvlpBase += inOvlpGapSize[curDim];
      outOvlpBase += outOvlpGapSize[curDim];
      pos[curDim] = 0;
      curDim--;
    } while (pos[curDim] == ovlpCount[curDim]);
  }
}
static void NdCopyIterDFDynamic(const char *inBase, char *outBase,
                                Dims &inRltvOvlpSPos, Dims &outRltvOvlpSPos,
                                Dims &inStride, Dims &outStride,
                                Dims &ovlpCount, size_t elmSize)
{
  size_t curDim = 0;
  Dims pos(ovlpCount.size() + 1, 0);
  std::vector<const char *> inAddr(ovlpCount.size() + 1);
  inAddr[0] = inBase;
  std::vector<char *> outAddr(ovlpCount.size() + 1);
  outAddr[0] = outBase;
  while (true)
  {
    while (curDim != inStride.size())
    {
      inAddr[curDim + 1] =
      inAddr[curDim] +
      (inRltvOvlpSPos[curDim] + pos[curDim]) * inStride[curDim];
      outAddr[curDim + 1] =
      outAddr[curDim] +
      (outRltvOvlpSPos[curDim] + pos[curDim]) * outStride[curDim];
      pos[curDim]++;
      curDim++;
    }
    std::memcpy(outAddr[curDim], inAddr[curDim], elmSize);
    do
    {
      if (curDim == 0)
        return;
      pos[curDim] = 0;
      curDim--;
    } while (pos[curDim] == ovlpCount[curDim]);
  }
}

static void NdCopyIterDFDynamicRevEndian(const char *inBase, char *outBase,
                                         Dims &inRltvOvlpSPos,
                                         Dims &outRltvOvlpSPos, Dims &inStride,
                                         Dims &outStride, Dims &ovlpCount,
                                         size_t elmSize)
{
  size_t curDim = 0;
  Dims pos(ovlpCount.size() + 1, 0);
  std::vector<const char *> inAddr(ovlpCount.size() + 1);
  inAddr[0] = inBase;
  std::vector<char *> outAddr(ovlpCount.size() + 1);
  outAddr[0] = outBase;
  while (true)
  {
    while (curDim != inStride.size())
    {
      inAddr[curDim + 1] =
      inAddr[curDim] +
      (inRltvOvlpSPos[curDim] + pos[curDim]) * inStride[curDim];
      outAddr[curDim + 1] =
      outAddr[curDim] +
      (outRltvOvlpSPos[curDim] + pos[curDim]) * outStride[curDim];
      pos[curDim]++;
      curDim++;
    }
    for (size_t i = 0; i < elmSize; i++)
    {
      outAddr[curDim][i] = inAddr[curDim][elmSize - 1 - i];
    }
    do
    {
      if (curDim == 0)
        return;
      pos[curDim] = 0;
      curDim--;
    } while (pos[curDim] == ovlpCount[curDim]);
  }
}


template <class T>
int NdCopy(const char *in, const Dims &inStart, const Dims &inCount,
           const bool inIsRowMajor, const bool inIsLittleEndian, char *out,
           const Dims &outStart, const Dims &outCount, const bool outIsRowMajor,
           const bool outIsLittleEndian, const Dims &inMemStart,
           const Dims &inMemCount, const Dims &outMemStart,
           const Dims &outMemCount, const bool safeMode)

{
    // use values of ioStart and ioCount if ioMemStart and ioMemCount are
    // left as default
    Dims inMemStartNC = inMemStart.empty() ? inStart : inMemStart;
    Dims inMemCountNC = inMemCount.empty() ? inCount : inMemCount;
    Dims outMemStartNC = outMemStart.empty() ? outStart : outMemStart;
    Dims outMemCountNC = outMemCount.empty() ? outCount : outMemCount;

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
    const char *inOvlpBase = nullptr;
    char *outOvlpBase = nullptr;
    auto GetInEnd = [](Dims &inEnd, const Dims &inStart, const Dims &inCount) {
        for (size_t i = 0; i < inStart.size(); i++)
            inEnd[i] = inStart[i] + inCount[i] - 1;
    };
    auto GetOutEnd = [](Dims &outEnd, const Dims &outStart,
                        const Dims &output_count) {
        for (size_t i = 0; i < outStart.size(); i++)
            outEnd[i] = outStart[i] + output_count[i] - 1;
    };
    auto GetOvlpStart = [](Dims &ovlpStart, const Dims &inStart,
                           const Dims &outStart) {
        for (size_t i = 0; i < ovlpStart.size(); i++)
            ovlpStart[i] = inStart[i] > outStart[i] ? inStart[i] : outStart[i];
    };
    auto GetOvlpEnd = [](Dims &ovlpEnd, Dims &inEnd, Dims &outEnd) {
        for (size_t i = 0; i < ovlpEnd.size(); i++)
            ovlpEnd[i] = inEnd[i] < outEnd[i] ? inEnd[i] : outEnd[i];
    };
    auto GetOvlpCount = [](Dims &ovlpCount, Dims &ovlpStart, Dims &ovlpEnd) {
        for (size_t i = 0; i < ovlpCount.size(); i++)
            ovlpCount[i] = ovlpEnd[i] - ovlpStart[i] + 1;
    };
    auto HasOvlp = [](Dims &ovlpStart, Dims &ovlpEnd) {
        for (size_t i = 0; i < ovlpStart.size(); i++)
            if (ovlpEnd[i] < ovlpStart[i])
                return false;
        return true;
    };

    auto GetIoStrides = [](Dims &ioStride, const Dims &ioCount,
                           size_t elmSize) {
        // ioStride[i] holds the total number of elements under each element
        // of the i'th dimension
        ioStride[ioStride.size() - 1] = elmSize;
        if (ioStride.size() > 1)
            ioStride[ioStride.size() - 2] =
                    ioCount[ioStride.size() - 1] * elmSize;
        if (ioStride.size() > 2)
        {
            size_t i = ioStride.size() - 3;
            while (true)
            {
                ioStride[i] = ioCount[i + 1] * ioStride[i + 1];
                if (i == 0)
                    break;
                else
                    i--;
            }
        }
    };

    auto GetInOvlpBase = [](const char *&inOvlpBase, const char *in,
                            const Dims &inStart, Dims &inStride,
                            Dims &ovlpStart) {
        inOvlpBase = in;
        for (size_t i = 0; i < inStart.size(); i++)
            inOvlpBase = inOvlpBase + (ovlpStart[i] - inStart[i]) * inStride[i];
    };
    auto GetOutOvlpBase = [](char *&outOvlpBase, char *out,
                             const Dims &outStart, Dims &outStride,
                             Dims &ovlpStart) {
        outOvlpBase = out;
        for (size_t i = 0; i < outStart.size(); i++)
            outOvlpBase =
                    outOvlpBase + (ovlpStart[i] - outStart[i]) * outStride[i];
    };
    auto GetIoOvlpGapSize = [](Dims &ioOvlpGapSize, Dims &ioStride,
                               const Dims &ioCount, Dims &ovlpCount) {
        for (size_t i = 0; i < ioOvlpGapSize.size(); i++)
            ioOvlpGapSize[i] = (ioCount[i] - ovlpCount[i]) * ioStride[i];
    };
    auto GetMinContDim = [](const Dims &inCount, const Dims outCount,
                            Dims &ovlpCount) {
        //    note: minContDim is the first index where its input box and
        //    overlap box
        //    are not fully match. therefore all data below this branch is
        //    contingous
        //    and this determins the Biggest continuous block size - Each
        //    element of the
        //    current dimension.
        size_t i = ovlpCount.size() - 1;
        while (true)
        {
            if (i == 0)
                break;
            if ((inCount[i] != ovlpCount[i]) || (outCount[i] != ovlpCount[i]))
                break;
            i--;
        }
        return i;
    };
    auto GetBlockSize = [](Dims &ovlpCount, size_t minContDim, size_t elmSize) {
        size_t res = elmSize;
        for (size_t i = minContDim; i < ovlpCount.size(); i++)
            res *= ovlpCount[i];
        return res;
    };

    auto GetRltvOvlpStartPos = [](Dims &ioRltvOvlpStart, const Dims &ioStart,
                                  Dims &ovlpStart) {
        for (size_t i = 0; i < ioStart.size(); i++)
            ioRltvOvlpStart[i] = ovlpStart[i] - ioStart[i];
    };

    // main flow
    // row-major ==> row-major mode
    // algrithm optimizations:
    // 1. contigous data copying
    // 2. mem pointer arithmetics by sequential padding. O(1) overhead/block
    if (inIsRowMajor == true && outIsRowMajor == true)
    {
        GetInEnd(inEnd, inStart, inCount);
        GetOutEnd(outEnd, outStart, outCount);
        GetOvlpStart(ovlpStart, inStart, outStart);
        GetOvlpEnd(ovlpEnd, inEnd, outEnd);
        GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
        if (!HasOvlp(ovlpStart, ovlpEnd))
            return 1; // no overlap found
        GetIoStrides(inStride, inMemCountNC, sizeof(T));
        GetIoStrides(outStride, outMemCountNC, sizeof(T));
        GetIoOvlpGapSize(inOvlpGapSize, inStride, inMemCountNC, ovlpCount);
        GetIoOvlpGapSize(outOvlpGapSize, outStride, outMemCountNC, ovlpCount);
        GetInOvlpBase(inOvlpBase, in, inMemStartNC, inStride, ovlpStart);
        GetOutOvlpBase(outOvlpBase, out, outMemStartNC, outStride, ovlpStart);
        minContDim = GetMinContDim(inMemCountNC, outMemCountNC, ovlpCount);
        blockSize = GetBlockSize(ovlpCount, minContDim, sizeof(T));
        // same endianess mode: most optimized, contiguous data copying
        // algorithm used.
        if (inIsLittleEndian == outIsLittleEndian)
        {
            // most efficient algm
            // warning: number of function stacks used is number of dimensions
            // of data.
            if (!safeMode)
                NdCopyRecurDFSeqPadding(0, inOvlpBase, outOvlpBase,
                                        inOvlpGapSize, outOvlpGapSize,
                                        ovlpCount, minContDim, blockSize);
            else // safeMode
                //      //alternative iterative version, 10% slower then
                //      recursive
                //      //use it when very high demension is used.
                NdCopyIterDFSeqPadding(inOvlpBase, outOvlpBase, inOvlpGapSize,
                                       outOvlpGapSize, ovlpCount, minContDim,
                                       blockSize);
        }
            // different endianess mode
        else
        {
            if (!safeMode)
                NdCopyRecurDFSeqPaddingRevEndian(
                        0, inOvlpBase, outOvlpBase, inOvlpGapSize, outOvlpGapSize,
                        ovlpCount, minContDim, blockSize, sizeof(T),
                        blockSize / sizeof(T));
            else
                NdCopyIterDFSeqPaddingRevEndian(
                        inOvlpBase, outOvlpBase, inOvlpGapSize, outOvlpGapSize,
                        ovlpCount, minContDim, blockSize, sizeof(T),
                        blockSize / sizeof(T));
        }
    }

        // Copying modes involing col-major
        // algorithm optimization:
        // 1. mem ptr arithmetics: O(1) overhead per block, dynamic/non-sequential
        // padding
    else
    {
        Dims revInCount(inCount);
        Dims revOutCount(outCount);
        GetInEnd(inEnd, inStart, inCount);
        GetOutEnd(outEnd, outStart, outCount);
        GetOvlpStart(ovlpStart, inStart, outStart);
        GetOvlpEnd(ovlpEnd, inEnd, outEnd);
        GetOvlpCount(ovlpCount, ovlpStart, ovlpEnd);
        if (!HasOvlp(ovlpStart, ovlpEnd))
            return 1; // no overlap found
        // col-major ==> col-major mode
        if (!inIsRowMajor && !outIsRowMajor)
        {
            std::reverse(revInCount.begin(), revInCount.end());
            GetIoStrides(inStride, revInCount, sizeof(T));
            std::reverse(inStride.begin(), inStride.end());
            std::reverse(revOutCount.begin(), revOutCount.end());
            GetIoStrides(outStride, revOutCount, sizeof(T));
            std::reverse(outStride.begin(), outStride.end());
        }
            // row-major ==> col-major mode
        else if (inIsRowMajor && !outIsRowMajor)
        {
            GetIoStrides(inStride, inCount, sizeof(T));
            std::reverse(revOutCount.begin(), revOutCount.end());
            GetIoStrides(outStride, revOutCount, sizeof(T));
            std::reverse(outStride.begin(), outStride.end());
        }
            // col-major ==> row-major mode
        else if (!inIsRowMajor && outIsRowMajor)
        {
            std::reverse(revInCount.begin(), revInCount.end());
            GetIoStrides(inStride, revInCount, sizeof(T));
            std::reverse(inStride.begin(), inStride.end());
            GetIoStrides(outStride, outCount, sizeof(T));
        }
        GetRltvOvlpStartPos(inRltvOvlpStartPos, inStart, ovlpStart);
        GetRltvOvlpStartPos(outRltvOvlpStartPos, outStart, ovlpStart);
        inOvlpBase = in;
        outOvlpBase = out;
        // Same Endian"
        if (inIsLittleEndian == outIsLittleEndian)
        {
            if (!safeMode)
                NdCopyRecurDFNonSeqDynamic(0, inOvlpBase, outOvlpBase,
                                           inRltvOvlpStartPos,
                                           outRltvOvlpStartPos, inStride,
                                           outStride, ovlpCount, sizeof(T));
            else
                NdCopyIterDFDynamic(inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                                    outRltvOvlpStartPos, inStride, outStride,
                                    ovlpCount, sizeof(T));
        }
            // different Endian"
        else
        {
            if (!safeMode)
                NdCopyRecurDFNonSeqDynamicRevEndian(
                        0, inOvlpBase, outOvlpBase, inRltvOvlpStartPos,
                        outRltvOvlpStartPos, inStride, outStride, ovlpCount,
                        sizeof(T));
            else
                NdCopyIterDFDynamicRevEndian(inOvlpBase, outOvlpBase,
                                             inRltvOvlpStartPos,
                                             outRltvOvlpStartPos, inStride,
                                             outStride, ovlpCount, sizeof(T));
        }
    }
    return 0;
}
//*************** End of NdCopy() and its 8 helpers ***************

#endif