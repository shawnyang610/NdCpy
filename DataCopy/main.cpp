//
//  main.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/17/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include <numeric>

#include "NDCopy2.tcc"
#include "NDCopy.h"

void PrintDims(const Dims &dims, const std::string &name){
    std::string s;
    s += name + " = [";
    for(auto i : dims){
        s += std::to_string(i) + ", ";
    }
    s.resize(s.size() - 2);
    s += "]";
    std::cout << s << std::endl;
}

template<class T>
void PrintData(const Buffer &buffer, const Dims &count){
    size_t total_elements = accumulate(count.begin(), count.end(), 1,std::multiplies<size_t>());
    for(size_t i = 0; i < total_elements; ++i){
        std::cout << reinterpret_cast<const T *>(buffer.data())[i] << "\t";
        if((i+1)%count.back() == 0){
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

template<class T>
void MakeData(Buffer &buffer, const Dims &count, bool zero){
    size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
    for(size_t i=0; i<size; ++i){
        if(zero){
            reinterpret_cast<T*>(buffer.data())[i] = 0;
        }
        else{
            reinterpret_cast<T*>(buffer.data())[i] = i;
        }
    }
}

template<class T>
void RunTest(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count){

    Buffer input_buffer, output_buffer;

    input_buffer.resize(std::accumulate(input_count.begin(), input_count.end(), sizeof(T), std::multiplies<size_t>()));
    MakeData<int>(input_buffer, input_count, false);

    output_buffer.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
    NdCopy2<int>(
            input_buffer,
            input_start,
            input_count,
            RowMajorBigEndian,
            output_buffer,
            output_start,
            output_count,
            RowMajorBigEndian
            );

    std::cout << "*************** input_buffer ****************" << std::endl;
    PrintData<int>(input_buffer, input_count);
    std::cout << "*************** output_buffer ****************" << std::endl;
    PrintData<int>(output_buffer, output_count);
}


int main(int argc, const char * argv[]) {

    Dims input_start = {10,20};
    Dims input_count = {5,8};

    Dims output_start = {0,0};
    Dims output_count = {30,30};

    RunTest<int>(input_start, input_count, output_start, output_count);

    return 0;
}
