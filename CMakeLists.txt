cmake_minimum_required(VERSION 2.8.12)
project(mini-pic)

#Bring the headers, such as Student.h into the project
include_directories(include)

#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.*")

add_executable(mini-pic ${SOURCES})

set_property(TARGET mini-pic PROPERTY CXX_STANDARD 20)

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(omp_libraries OpenMP::OpenMP_CXX)
endif()

target_link_libraries(mini-pic PUBLIC ${lapackblas_libraries} ${omp_libraries})

set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

