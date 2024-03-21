cmake_minimum_required(VERSION 3.00)

project(cmake_test CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-Ofast -funroll-loops -fomit-frame-pointer -march=native -mtune=native")

include_directories(include)
include_directories(/home/user9/local_libs/openssl/include/)
include_directories(/usr/local/include)

add_executable(test test.cpp)
target_link_libraries(test pthread crypto)

add_executable(luttest LUT_test.cpp)
target_link_libraries(luttest pthread crypto)

add_executable(face FacialRecognition.cpp)
target_link_libraries(face pthread crypto)

add_executable(mnist MNIST.cpp)
target_link_libraries(mnist pthread crypto)

add_executable(object ObjectClassification.cpp)
target_link_libraries(object pthread crypto)

add_executable(speech SpeechRecognition.cpp)
target_link_libraries(speech pthread crypto)

add_executable(text TextClassification.cpp)
target_link_libraries(text pthread crypto)