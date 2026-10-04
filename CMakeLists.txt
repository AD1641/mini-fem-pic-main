cmake_minimum_required(VERSION 3.20)
project(mini-pic CXX)

# Build optimised unless asked otherwise (e.g. cmake -DCMAKE_BUILD_TYPE=Debug ..)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB SOURCES "src/*.cpp" "src/*.h")

add_executable(mini-pic ${SOURCES})

set_property(TARGET mini-pic PROPERTY CXX_STANDARD 20)

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# Apple's compiler doesn't come with OpenMP, so point it at Homebrew's (brew install libomp)
if(APPLE AND NOT OpenMP_ROOT)
  execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE OpenMP_ROOT
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  set(omp_libraries OpenMP::OpenMP_CXX)
endif()

target_link_libraries(mini-pic PUBLIC ${lapackblas_libraries} ${omp_libraries})

# SYCL: Intel's oneAPI compiler (icpx) supports it directly, anything else uses
# AdaptiveCpp (e.g. on macOS: brew install adaptivecpp)
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_options(mini-pic PRIVATE -fsycl)
  target_link_options(mini-pic PRIVATE -fsycl)
else()
  find_package(AdaptiveCpp CONFIG REQUIRED)
  add_sycl_to_target(TARGET mini-pic SOURCES ${SOURCES})
endif()
