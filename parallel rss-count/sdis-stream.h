#ifndef SDIX_STREAM_H
#define SDIX_STREAM_H

#ifndef BUFFER
// try 8192 bytes 819200 or go to 32kb 32,768
#define BUFFER 409600 // 4kb = 4 * 2 ^ 10   409600
#endif

#include <array>
#include <cstring>
#include <iostream>

/*
I need to generate the stream here and then modify the rss-count to remove the std::ifstream.
run skyline takes only name of function, dimensionality, window.
*/
namespace sdistream
{

  template <class T>
  auto input(T &in, size_t width, double *buffer) -> bool
  {
    auto delim = ", ";
    std::array<char, BUFFER> line{};
    auto data = line.data();
    in.getline(data, BUFFER);
    if (!in.good())
    {
      return false;
    }
    char *ptr = strtok(data, delim);
    size_t n = 0;
    while (ptr != nullptr && n < width)
    {
      buffer[n] = strtod(ptr, nullptr);
      ptr = strtok(nullptr, delim);
      ++n;
    }
    return true;
  }

}

#endif // SDIX_STREAM_H
