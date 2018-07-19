
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


int main(int argc, const char * argv[]) {
    int iters = 1;
    if(argc > 1){
        iters = atoi(argv[1]);
    }
  
//    //memory addr. calc. performance test:
//    //both algorithm, only 1 element is copied at a time
//        Dims input_start = {1,1,1,1,1,1,1,1,1,1};
//        Dims input_count = {5,5,5,5,5,5,5,5,1,5};
//        Dims output_start = {1,1,1,1,1,1,1,1,1,1};
//        Dims output_count = {5,5,5,5,5,5,5,5,5,5};
//        RunTest<int>(input_start, input_count, output_start, output_count, iters);

    //largest-continous-block-method performance test
        Dims input_start = {1,1,1,1,1,1,1,1,1,1};
        Dims input_count = {5,1,5,5,5,5,5,5,5,5};
        Dims output_start = {1,1,1,1,1,1,1,1,1,1};
        Dims output_count = {5,5,5,5,5,5,5,5,5,5};
        RunTest<int>(input_start, input_count, output_start, output_count, iters);

  //    input:row major, output:col major, same-endian demo
//    Dims input_start = {5,10};
//    Dims input_count = {5,10};
//    Dims output_start = {10,5};
//    Dims output_count = {20,10};
//    RunTestDiffMajorMode<int>(input_start, input_count, output_start, output_count);
////    RunTest<int>(input_start, input_count, output_start, output_count, iters);

  
//  //DEMO input of any major to output of any major,with same endianess
//  bool inIsRowMaj=true;
//  Dims input_start = {3,0}; //in its own format
//  Dims input_count = {5,4}; //in its own format
//  bool outIsRowMaj=false;
//  Dims output_start = {3,0}; //in its own format
//  Dims output_count = {10,5}; //in its own format
//  bool safeMode=false;
//  RunTestDiffMajorMode<size_t>(input_start, input_count,
//                            output_start, output_count,
//                            inIsRowMaj,outIsRowMaj,safeMode);
  
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
  //endian mode demo
  //unsigned DEC 4278255360 == BIN 11111111 00000000 11111111 00000000
  //unsigned DEC 16711935   == BIN 00000000 11111111 00000000 11111111
//  Dims input_start = {2,4};
//  Dims input_count = {3,3};
//  Dims output_start = {0,0};
//  Dims output_count = {10,10};
//  bool inIsRowMaj=true;
//  bool inIsBigMaj=true;
//  bool outIsRowMaj=true;
//  bool outIsBigMaj=true;
//  bool safeMode=false;
//
//  RunTestEndianMode<unsigned>(input_start, input_count, output_start,
//                            output_count, inIsRowMaj,inIsBigMaj,outIsRowMaj,
//                              outIsBigMaj,safeMode);


  
    return 0;
}
