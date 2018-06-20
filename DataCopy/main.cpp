//
//  main.cpp
//  DataCopy
//
//  Created by Shawn Yang on 6/17/18.
//  Copyright © 2018 Shawn Yang. All rights reserved.
//

#include <iostream>
#include <numeric>

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

Dims GenGlobalShape(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count){
    Dims global_shape(input_count.size());
    for(size_t i=0; i<input_count.size(); ++i){
        global_shape[i] = (input_start[i] + input_count[i] < output_start[i] + output_count[i]) ? output_start[i] + output_count[i] : input_start[i] + input_count[i];
    }
    return global_shape;
}

template<class T>
std::shared_ptr<std::vector<char>> MakeData(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count){
    // TODO: verify data generated
    Dims global_shape = GenGlobalShape(input_start, input_count, output_start, output_count);
    size_t global_size = std::accumulate(global_shape.begin(), global_shape.end(), 1, std::multiplies<size_t>());
    size_t global_bytes = global_size * sizeof(T);
    std::shared_ptr<std::vector<char>> input = std::make_shared<std::vector<char>>(global_size);
    for(size_t i=0; i<global_size; ++i){
        reinterpret_cast<T*>(input->data())[i] = i;
    }
    return input;
}

int main(int argc, const char * argv[]) {

    Dims input_start = {20,40,60};
    Dims input_count = {100,200,300};

    Dims output_start = {100,20,300};
    Dims output_count = {200,100,400};

    // Generate data
    auto input_buffer = MakeData<int>(input_start, input_count, output_start, output_count);

    PrintDims(input_start, "input_start");



    return 0;
}
