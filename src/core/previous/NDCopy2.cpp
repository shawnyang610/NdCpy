#include <vector>
#include <numeric>
#include <functional>
#include <cstring>

#include "NDCopy2.tcc"

bool IsContinuous(const Box<Dims> &inner,
                                       const Box<Dims> &outer)
{
    for (size_t i = 1; i < inner.first.size(); ++i)
    {
        if (inner.first[i] != outer.first[i])
        {
            return false;
        }
        if (inner.second[i] != outer.second[i])
        {
            return false;
        }
    }
    return true;
}

Dims GetAbsolutePosition(const Dims &inner,
                                              const Dims &outer)
{
    Dims ret;
    size_t size = inner.size();
    ret.resize(size);
    for (int i = 0; i < size; ++i)
    {
        ret[i] = inner[i] + outer[i];
    }
    return ret;
}

Dims GetRelativePosition(const Dims &inner,
                                              const Dims &outer)
{
    Dims ret;
    size_t size = inner.size();
    ret.resize(size);
    for (int i = 0; i < size; ++i)
    {
        ret[i] = inner[i] - outer[i];
    }
    return ret;
}

Dims OneToMulti(const Dims &global, size_t position)
{
    std::vector<size_t> index(global.size());
    for (int i = 1; i < global.size(); ++i)
    {
        size_t s = std::accumulate(global.begin() + i, global.end(), 1,
                                   std::multiplies<size_t>());
        index[i - 1] = position / s;
        position -= index[i - 1] * s;
    }
    index.back() = position;
    return index;
}

size_t MultiToOne(const Dims &global, const Dims &position)
{
    size_t index = 0;
    for (int i = 1; i < global.size(); ++i)
    {
        index += std::accumulate(global.begin() + i, global.end(),
                                 position[i - 1], std::multiplies<size_t>());
    }
    index += position.back();
    return index;
}

bool GetOverlap(const Box<Dims> &b1, const Box<Dims> &b2,
                                     Box<Dims> &overlapBox)
{
    overlapBox.first.resize(b1.first.size());
    overlapBox.second.resize(b1.first.size());

    for (size_t i = 0; i < b1.first.size(); ++i)
    {
        if (b1.first[i] > b2.first[i])
        {
            overlapBox.first[i] = b1.first[i];
        }
        else
        {
            overlapBox.first[i] = b2.first[i];
        }
        if (b1.second[i] < b2.second[i])
        {
            overlapBox.second[i] = b1.second[i];
        }
        else
        {
            overlapBox.second[i] = b2.second[i];
        }
    }

    for (size_t i = 0; i < overlapBox.first.size(); ++i)
    {
        if (overlapBox.first[i] > overlapBox.second[i])
        {
            return false;
        }
    }

    return true;
}

void CopyLocalToGlobal(char *dst, const Box<Dims> &dstBox,
                                            const char *src,
                                            const Box<Dims> &srcBox,
                                            const size_t size,
                                            const Box<Dims> &overlapBox)
{

    size_t dimensions = overlapBox.first.size();
    size_t overlapSize = 1;
    for (int i = 0; i < dimensions; ++i)
    {
        overlapSize =
            overlapSize * (overlapBox.second[i] - overlapBox.first[i]);
    }

    Dims srcCount(dimensions);
    for (int i = 0; i < dimensions; ++i)
    {
        srcCount[i] = srcBox.second[i] - srcBox.first[i];
    }

    Dims dstCount(dimensions);
    for (int i = 0; i < dimensions; ++i)
    {
        dstCount[i] = dstBox.second[i] - dstBox.first[i];
    }

    if (IsContinuous(overlapBox, srcBox) && IsContinuous(overlapBox, dstBox))
    {
        Dims overlapInSrcRelativeLeftBoundary =
            GetRelativePosition(overlapBox.first, srcBox.first);
        Dims overlapInDstRelativeLeftBoundary =
            GetRelativePosition(overlapBox.first, dstBox.first);
        size_t srcStartPtrOffset =
            MultiToOne(srcCount, overlapInSrcRelativeLeftBoundary);
        size_t dstStartPtrOffset =
            MultiToOne(dstCount, overlapInDstRelativeLeftBoundary);
        std::memcpy(dst + dstStartPtrOffset * size,
                    src + srcStartPtrOffset * size, overlapSize * size);
    }
    else
    {
        size_t overlapChunkSize =
            (overlapBox.second.back() - overlapBox.first.back());

        Dims overlapCount(dimensions);
        for (int i = 0; i < dimensions; ++i)
        {
            overlapCount[i] = overlapBox.second[i] - overlapBox.first[i];
        }

        for (size_t i = 0; i < overlapSize; i += overlapChunkSize)
        {
            Dims currentPositionLocal = OneToMulti(overlapCount, i);
            Dims currentPositionGlobal =
                GetAbsolutePosition(currentPositionLocal, overlapBox.first);
            Dims overlapInSrcRelativeCurrentPosition =
                GetRelativePosition(currentPositionGlobal, srcBox.first);
            Dims overlapInDstRelativeCurrentPosition =
                GetRelativePosition(currentPositionGlobal, dstBox.first);
            size_t srcStartPtrOffset =
                MultiToOne(srcCount, overlapInSrcRelativeCurrentPosition);
            size_t dstStartPtrOffset =
                MultiToOne(dstCount, overlapInDstRelativeCurrentPosition);
            std::memcpy(dst + dstStartPtrOffset * size,
                        src + srcStartPtrOffset * size,
                        overlapChunkSize * size);
        }
    }
}

