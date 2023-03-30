#ifndef SDIX_STREAM_H
#define SDIX_STREAM_H

#ifndef BUFFER
#define BUFFER 409600
#endif

#include <array>
#include <cstring>
#include <iostream>

namespace sdistream
{

  template <class IN>
  auto input(IN &in, size_t width, double *buffer) -> bool
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