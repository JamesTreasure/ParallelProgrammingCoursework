cmake_minimum_required(VERSION 3.6)
project(Parallel)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES distributedMemory.c)
#set(SOURCE_FILES linux.c)
add_executable(Parallel ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(Parallel pthread)
target_link_libraries(Parallel m)