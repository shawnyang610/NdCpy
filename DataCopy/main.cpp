//
//  main.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/17/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include <numeric>
#include <chrono>
#include "NDCopy.cpp"
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
void RunTest(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count, int iters){

    Buffer input_buffer, output_buffer, output_buffer2;

    input_buffer.resize(std::accumulate(input_count.begin(), input_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer2.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));

    MakeData<int>(input_buffer, input_count, false);

    // performance testing begin

    auto start = std::chrono::system_clock::now();
    for(int i=0; i<iters; ++i){
        if(NdCopy<int>(
                    input_buffer,
                    input_start,
                    input_count,
                    RowMajorBigEndian,
                    output_buffer,
                    output_start,
                    output_count,
                    RowMajorBigEndian
                    ))
        {
            std::cout<<"no overlap found"<<std::endl;
        }
    }
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    auto start2 = std::chrono::system_clock::now();
    for(int i=0; i<iters; ++i){
        NdCopy2<int>(
                input_buffer,
                input_start,
                input_count,
                RowMajorBigEndian,
                output_buffer2,
                output_start,
                output_count,
                RowMajorBigEndian
                );
    }
    auto end2 = std::chrono::system_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);

    // performance testing end


    // verify data

    if(output_buffer.size() != output_buffer2.size()){
        std::cout << "Two buffers do not have the same size!" << std::endl;;
    }
    else{
        for(size_t i=0; i<output_buffer.size(); ++i){
            if(output_buffer[i] != output_buffer2[i]){
                std::cout << "Data not correct!" << std::endl;
            }
        }
    }

//    std::cout << "*************** input_buffer ****************" << std::endl;
//    PrintData<int>(input_buffer, input_count);
//    std::cout << "*************** output_buffer ****************" << std::endl;
//    PrintData<int>(output_buffer, output_count);
//    std::cout << "*************** output_buffer2 ****************" << std::endl;
//    PrintData<int>(output_buffer2, output_count);

    std::cout << "Algorithm 1 spent: " << duration.count() << std::endl;
    std::cout << "Algorithm 2 spent: " << duration2.count() << std::endl;
    if(duration2.count() > duration.count()){
        std::cout << "Algorithm 1 is faster than Algorithm 2 by " << ((float)duration2.count() - (float)duration.count()) / (float)duration.count() * 100.0 << "%" << std::endl;
    }
    else{
        std::cout << "Algorithm 2 is faster than Algorithm 1 by " << ((float)duration.count() - (float)duration2.count()) / (float)duration2.count() * 100.0 << "%" << std::endl;
    }
}

int main(int argc, const char * argv[]) {

    int iters = 10;

    if(argc > 1){
        iters = atoi(argv[1]);
    }

    Dims input_start = {0,0,0};
    Dims input_count = {200,400,600};

    Dims output_start = {0,0,0};
    Dims output_count = {800,800,800};

    RunTest<int>(input_start, input_count, output_start, output_count, iters);

//    input_start = {10,20,15,5,2,3,5,5,5};
//    input_count = { 2, 2, 2,4,2,3,5,5,5};
//    output_start = {8,20,15,5,2,3,5,5,5};
//    output_count = {6, 5, 2,4,2,3,5,5,5};
//    RunTest<int>(input_start, input_count, output_start, output_count);
//    RunTest2<int>(input_start, input_count, output_start, output_count);

    return 0;
}
