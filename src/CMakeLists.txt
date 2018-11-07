cmake_minimum_required(VERSION 3.12)
project(DataCopy)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(src
        main.cpp
        core/NdCpy/NDCopy.hpp
        core/previous/NDCopy2.h
        core/previous/NDCopy2.cpp
        core/previous/NDCopy2.tcc
        core/NdCpy/NDCopyCore.cpp tests/test.cpp tests/test.h)
