cmake_minimum_required(VERSION 3.15)
project(sdi-stream)

set(CMAKE_CXX_STANDARD 17)

set(SDIS
        sdis-cache.cpp
        sdis-cache.h
        sdis-skyline.cpp
        sdis-skyline.h
        sdis-stream.h
        timer.cpp
        timer.h
        )

add_executable(rss-count rss-count.cpp rss-count.h ${SDIS})
