cmake_minimum_required (VERSION 3.16)
project (Harmony NONE)

if (NOT CMAKE_BUILD_TYPE)
 # set(CMAKE_BUILD_TYPE RelWithDebInfo)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT DEFINED BUILD_EXPLICIT_SOLVER)
  set(BUILD_EXPLICIT_SOLVER ON)
endif()

message( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
message( ${CMAKE_MODULE_PATH} )
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

enable_language(C)
enable_language(CXX)

if(BUILD_EXPLICIT_SOLVER)
  find_package(MPI REQUIRED)
  if (MPI_FOUND)
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    #message(STATUS "MPI INCLUDE ${MPI_CXX_INCLUDE_DIRS}")
    #message(STATUS "MPI HOME ${MPI_CXX_LIBRARIES}")
  endif()  
endif()

#add_compile_options(-g -O0)
add_compile_options(-O3 -lm)
set (CMAKE_CXX_STANDARD 17)

if(BUILD_EXPLICIT_SOLVER)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()


if(BUILD_EXPLICIT_SOLVER)
add_subdirectory(Source)
endif()
