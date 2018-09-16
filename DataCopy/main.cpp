
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
            reinterpret_cast<T*>(buffer.data())[i] = i*1.0;
        }
    }
}

template<class T>
void RunTest(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count, int iters){

    Buffer input_buffer, output_buffer, output_buffer2;

    input_buffer.resize(std::accumulate(input_count.begin(), input_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer2.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));

    MakeData<T>(input_buffer, input_count, false);

    // performance testing begin
  NdCopyFlag input_flag, output_flag;
  input_flag.isRowMajor = true;
  input_flag.isBigEndian = true;
  output_flag.isRowMajor = true;
  output_flag.isBigEndian = true;
  bool inIsRowMaj=true;
  bool inIsBigEndian=true;
  bool outIsRowMaj=true;
  bool outIsBigEndian=true;
  bool safeMode=false;

    auto start = std::chrono::system_clock::now();
    for(int i=0; i<iters; ++i){
        if(NdCopy<T>(
                    input_buffer.data(),
                    input_start,
                    input_count,
                       inIsRowMaj,
                       inIsBigEndian,
                    output_buffer.data(),
                    output_start,
                    output_count,
                       outIsRowMaj,
                       outIsBigEndian,
                       safeMode
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
                input_flag,
                output_buffer2,
                output_start,
                output_count,
                output_flag
                );
//      if(NdCopy<T>(
//                     input_buffer.data(),
//                     input_start,
//                     input_count,
//                     inIsRowMaj,
//                     inIsBigEndian,
//                     output_buffer2.data(),
//                     output_start,
//                     output_count,
//                     outIsRowMaj,
//                     outIsBigEndian,
//                     safeMode
//                     ))
//      {
//        std::cout<<"no overlap found"<<std::endl;
//      }
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

  std::cout << "Algorithm 1 spent: Total: " << duration.count()<<";  seconds per iter: "
  <<duration.count()*1e-6/iters <<std::endl;
  std::cout << "Algorithm 2 spent: Total: " << duration2.count()<<";  seconds per iter: "
  <<duration2.count()*1e-6/iters<< std::endl;
    if(duration2.count() > duration.count()){
        std::cout << "Algorithm 1 is faster than Algorithm 2 by " << ((float)duration2.count() - (float)duration.count()) / (float)duration.count() * 100.0 << "%" << std::endl;
    }
    else{
        std::cout << "Algorithm 2 is faster than Algorithm 1 by " << ((float)duration.count() - (float)duration2.count()) / (float)duration2.count() * 100.0 << "%" << std::endl;
    }
}



template<class T>
void RunTestDiffMajorMode( Dims &input_start, Dims &input_count, Dims &output_start,
        Dims &output_count, bool inIsRowMaj, bool outIsRowMaj, bool safeMode)
{
    Buffer input_buffer, output_buffer, output_buffer2;
  if (!inIsRowMaj) std::reverse(input_count.begin(), input_count.end());
  if (!outIsRowMaj) std::reverse(output_count.begin(),output_count.end());

    input_buffer.resize(std::accumulate(input_count.begin(), input_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
    output_buffer2.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));


  bool inIsBigEndian=true;
  bool outIsBigEndian=true;


    MakeData<T>(input_buffer, input_count, false);
  if (!inIsRowMaj) std::reverse(input_count.begin(), input_count.end());
  if (!outIsRowMaj) std::reverse(output_count.begin(),output_count.end());
  auto start = std::chrono::system_clock::now();
  if(NdCopy<T>(
                 input_buffer.data(),
                 input_start,
                 input_count,
                 inIsRowMaj,
                 inIsBigEndian,
                 reinterpret_cast<char*>(output_buffer.data()),
                 output_start,
                 output_count,
                 outIsRowMaj,
                 outIsBigEndian,
                 safeMode
                 ))
    {
        std::cout<<"no overlap found"<<std::endl;
    }
  auto end = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
  
  if (!inIsRowMaj) std::reverse(input_count.begin(), input_count.end());
  if (!outIsRowMaj) std::reverse(output_count.begin(),output_count.end());
    std::cout<<"time spent: "<<duration.count()<<" usec"<<std::endl;
    std::cout << "*************** input_buffer ****************" << std::endl;
    PrintData<T>(input_buffer, input_count);
    std::cout << "*************** output_buffer ****************" << std::endl;
    PrintData<T>(output_buffer, output_count);
}


template<class T>
void MakeDataEndianMode(Buffer &buffer, const Dims &count){
  size_t size = std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>());
  for(size_t i=0; i<size; ++i){
      reinterpret_cast<T*>(buffer.data())[i] = 4278255360;
  }
}
template<class T>
void RunTestEndianMode(const Dims &input_start, const Dims &input_count, const Dims &output_start, const Dims &output_count,
                       bool inIsRowMaj, bool outIsRowMaj,bool inIsBigEndian ,bool outIsBigEndian,bool safeMode){
  Buffer input_buffer, output_buffer, output_buffer2;
  
  input_buffer.resize(std::accumulate(input_count.begin(), input_count.end(), sizeof(T), std::multiplies<size_t>()));
  output_buffer.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
  output_buffer2.resize(std::accumulate(output_count.begin(), output_count.end(), sizeof(T), std::multiplies<size_t>()));
  

  
  MakeDataEndianMode<T>(input_buffer, input_count);
  if(NdCopy<T>(
                input_buffer.data(),
                input_start,
                input_count,
                inIsRowMaj,
                inIsBigEndian,
                reinterpret_cast<char*>(output_buffer.data()),
                output_start,
                output_count,
                outIsRowMaj,
                outIsBigEndian,
                safeMode
                ))
  {
    std::cout<<"no overlap found"<<std::endl;
  }
  std::cout << "*************** input_buffer ****************" << std::endl;
  PrintData<T>(input_buffer, input_count);
  std::cout << "*************** output_buffer ****************" << std::endl;
  PrintData<T>(output_buffer, output_count);
}




void performance_test_pointer_arithmetic_optimization(int iters){
    //memory addr. calc. performance test:
    //both algorithm, only 1 element is copied at a time
    std::cout<<"pointer arithmetic performance test, copies by element."<<std::endl<<"compares ndcopy algorithm against the old."<<std::endl;
  
  std::cout<<"3 dimensional data:"<<std::endl;
  Dims input_start = {1,1,1};
  Dims input_count = {5,1,5};
  Dims output_start = {1,1,1};
  Dims output_count = {5,5,5};
  RunTest<int>(input_start, input_count, output_start, output_count, iters);
  
  std::cout<<"10 dimensional data:"<<std::endl;
    input_start = {1,1,1,1,1,1,1,1,1,1};
    input_count = {5,5,5,5,5,5,5,5,1,5};
    output_start = {1,1,1,1,1,1,1,1,1,1};
    output_count = {5,5,5,5,5,5,5,5,5,5};
    RunTest<int>(input_start, input_count, output_start, output_count, iters);
}

void performance_test_max_cont_block_optimization(int iters){
    std::cout<<"max contiguous block optimized perfomance test."<<std::endl;
    std::cout<<"3 dimensional data:"<<std::endl;
    Dims input_start = {1,1,1};
    Dims input_count = {5,1,5};
    Dims output_start = {1,1,1};
    Dims output_count = {5,5,5};
    RunTest<int>(input_start, input_count, output_start, output_count, iters);
    std::cout<<"10  dimensional data:"<<std::endl;
    input_start = {1,1,1,1,1,1,1,1,1,1};
    input_count = {5,1,5,5,5,5,5,5,5,5};
    output_start = {1,1,1,1,1,1,1,1,1,1};
    output_count = {5,5,5,5,5,5,5,5,5,5};
    RunTest<int>(input_start, input_count, output_start, output_count, iters);
}

void demo_reversed_major_copy(){
    // input:row major, output:col major, same-endian demo
    std::cout<<"copy from row major to col major, 2d data:"<<std::endl;
    Dims input_start = {5,10};
    Dims input_count = {5,10};
    Dims output_start = {5,10};
    Dims output_count = {10,20};
    RunTestDiffMajorMode<int>(input_start, input_count, output_start, output_count, true, false, false);
}

void demo_copy_between_any_majors(){
    std::cout<<"customized copy between any majors, 2d data:"<<std::endl;
    bool inIsRowMaj=true;
    Dims input_start = {3,0}; //in its own format
    Dims input_count = {5,4}; //in its own format
    bool outIsRowMaj=false;
    Dims output_start = {3,0}; //in its own format
    Dims output_count = {10,5}; //in its own format
    bool safeMode=false;
    RunTestDiffMajorMode<size_t>(input_start, input_count,
                                 output_start, output_count,
                                 inIsRowMaj,outIsRowMaj,safeMode);
}

void demo_copy_between_reversed_endians(){
    // endian mode demo
    // unsigned DEC 4278255360 == BIN 11111111 00000000 11111111 00000000
    // unsigned DEC 16711935   == BIN 00000000 11111111 00000000 11111111
    std::cout<<"copy between reversed endians, 2d data:"<<std::endl;
    std::cout<<"unsigned DEC 4278255360 = BIN 11111111 00000000 11111111 00000000"<<std::endl;
    std::cout<<"unsigned DEC 16711935 = BIN 00000000 11111111 00000000 11111111"<<std::endl;

    Dims input_start = {2,4};
    Dims input_count = {3,3};
    Dims output_start = {0,0};
    Dims output_count = {10,10};
    bool inIsRowMaj=true;
    bool inIsBigEnd=true;
    bool outIsRowMaj=true;
    bool outIsBigEnd=false;
    bool safeMode=false;

    RunTestEndianMode<unsigned>(input_start, input_count, output_start,
                                output_count, inIsRowMaj,inIsBigEnd,outIsRowMaj,
                                outIsBigEnd,safeMode);
}

int main(int argc, const char * argv[]) {
    int iters = 1;
    if(argc > 1){
        iters = atoi(argv[1]);
    }
  
  std::cout<<std::endl<<"demo 1:"<<std::endl;
  performance_test_pointer_arithmetic_optimization(iters);

  std::cout<<std::endl<<"demo 2:"<<std::endl;
  performance_test_max_cont_block_optimization(iters);

  std::cout<<std::endl<<"demo 3:"<<std::endl;
  // copy from row-maj to col-maj, same endianess demo
  demo_reversed_major_copy();
  
  std::cout<<std::endl<<"demo 4:"<<std::endl;
  // DEMO input of any major to output of any major,with same endianess
  demo_copy_between_any_majors();
  
  std::cout<<std::endl<<"demo 5:"<<std::endl;
  // DEMO copy between reversed endians
  demo_copy_between_reversed_endians();
  
  
  
  
  //    RunTest<int>(input_start, input_count, output_start, output_count, iters);
    //performance test: diff-maj-same-endian vs same-maj-same-endian:
    //copy 1 element at a time. anticipation: same-maj algm should be slightly faster
    //due to the way mem address is calculated in diff-maj where 2 additional multiplications
    //are used for each element copied, but in both case, avg overhead is still O(1)
  
//  Dims input_start = {1,1,1,1,1,1,1,1,1};
//  Dims input_count = {5,1,5,5,5,5,5,5,5};
//  Dims output_start = {1,1,1,1,1,1,1,1,1};
//  Dims output_count = {5,5,5,5,5,5,5,5,5};
//  //need to comment out data verification part of RunTest before runnning.
//  RunTest<int>(input_start, input_count,output_start, output_count, iters);
//
    return 0;
}
