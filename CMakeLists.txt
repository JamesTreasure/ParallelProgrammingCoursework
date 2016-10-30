cmake_minimum_required(VERSION 3.6)
project(Parallel)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES test.c)
#set(SOURCE_FILES old-from-balena.c)
#set(SOURCE_FILES threadedHelloWorld.c)
#set(SOURCE_FILES githubParallel.c)
#set(SOURCE_FILES workSplitTest.c)
add_executable(Parallel ${SOURCE_FILES})