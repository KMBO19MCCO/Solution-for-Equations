#cmake_minimum_required(VERSION 3.23)
#project(Parfenov)
#
#set(CMAKE_CXX_STANDARD 23)
#
#add_executable(Parfenov main.cpp excerptMod.h excerptMod2.h PolynomialRoots.hh)
#
#include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
##target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.a)
#
##find_package(OpenMP REQUIRED)
##include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
##add_executable(nsmfsqe main.cpp Timer.cpp Timer.h)
##link_libraries()
#target_link_libraries(Parfenov PRIVATE ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll)
#file(COPY ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll
#        DESTINATION ${PROJECT_BINARY_DIR})


cmake_minimum_required(VERSION 3.23)
project(Parfenov)

set(CMAKE_CXX_STANDARD 23)

include_directories(${PROJECT_SOURCE_DIR}/excerpt/include)
add_executable(Parfenov main.cpp excerpt/include/excerpt.cpp)
target_link_libraries(Parfenov PRIVATE ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll)
file(COPY ${PROJECT_SOURCE_DIR}/excerpt/libexcerpt.dll
        DESTINATION ${PROJECT_BINARY_DIR})