# SDL2 skeleton for DH2323.
project(DgiLab)

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall")

add_executable(${PROJECT_NAME}
  skeleton.cpp
  SDL2Auxiliary.cpp
  )

find_package (SDL2 REQUIRED)
find_package (GLM REQUIRED)
# find_package (OpenGL REQUIRED)  # You might have to use this instead
                                  # of GLM above.

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SDL_DIR /usr/local/Cellar/sdl2/2.0.14_1/) # 这个SDL开发包的路径，可以通过brew info sdl2查出来
include_directories(${SDL_DIR}/include/) # 添加SDL头文件搜索路径
link_libraries(${SDL_DIR}/lib/) # 增加SDL链接库目录

include_directories(
  ${GLM_INCLUDE_DIRS}
  # ${OpenGL_INCLUDE_DIRS}  # You might have to use this instead of
                            # GLM above.
  ${SDL2_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME}
  ${SDL2_LIBRARIES}
  )
