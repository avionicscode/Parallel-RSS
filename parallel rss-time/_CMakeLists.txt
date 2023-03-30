cmake_minimum_required(VERSION 3.15)
project(sdi-stream)

set(CMAKE_CXX_STANDARD 11)

set(SDIS
        sdis-cache.cpp
        sdis-cache.h
        sdis-skyline.cpp
        sdis-skyline.h
        sdis-stream.h
        timer.cpp
        timer.h
        types.h
        )

add_executable(rss-count rss-count.cpp rss-count.h ${SDIS})
add_executable(rss-time rss-time.cpp rss-time.h ${SDIS})
set_target_properties(rss-time PROPERTIES COMPILE_DEFINITIONS "WITH_TIME_WINDOW=1")

set(SDISi
        sdis-index.cpp
        sdis-index.h
        sdis-skyline.cpp
        sdis-skyline.h
        sdis-stream.h
        timer.cpp
        timer.h
        types.h
        )

add_executable(rssi-count rssi-count.cpp rssi-count.h ${SDISi})
add_executable(rssi-time rssi-time.cpp rssi-time.h ${SDISi})
set_target_properties(rssi-time PROPERTIES COMPILE_DEFINITIONS "WITH_TIME_WINDOW=1")
