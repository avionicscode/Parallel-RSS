#ifndef SDIS_SKYLINE_H
#define SDIS_SKYLINE_H

#ifndef SLICE
#define SLICE 32
#endif

#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <omp.h>

namespace sdistream
{

  class skyline
  {
    friend auto operator<<(std::ostream &, const skyline &) -> std::ostream &;

  public:
    static size_t DT;
    skyline() = default;
    virtual ~skyline() = default;
    // Add a skyline point to d-tree.
    auto add(const size_t &) -> skyline &;
    // The first skyline point dominates the second point.
    auto append(const size_t &, const size_t &) -> skyline &;
    // Check whether a point is skyline point.
    auto contains(const size_t &) -> bool;
    // Get all dominated points of a skyline point.
    auto get(const size_t &) -> std::vector<size_t> &;
    // Move the first skyline point to the second skyline point.
    auto move(const size_t &, const size_t &) -> skyline &;
    // Remove a skyline point.
    auto remove(const size_t &) -> skyline &;
    // The number of skyline points.
    auto size() -> size_t;

  private:
    static auto slice_(size_t) -> size_t;
    size_t count_ = 0;
    std::vector<size_t> empty_;
    std::array<std::unordered_map<size_t, std::vector<size_t>>, SLICE> tree_;
  };

  template <class T>
  auto dominate(T *row1, T *row2, size_t width) -> bool
  {
    //++skyline::DT;
    T *p1 = row1;
    T *p2 = row2;
    bool dominating = false;
    for (size_t i = 0; i < width; ++i, ++p1, ++p2)
    {
      if (*p1 > *p2)
      {
        return false;
      }
      else if (*p1 < *p2 && !dominating)
      {
        dominating = true;
      }
    }
    return dominating;
  }

  template <class T>
  auto par_dominate(T *row1, T *row2, size_t width) -> bool
  {
    //++skyline::DT;
    T *p1 = row1;
    T *p2 = row2;
    bool dominating = false;
    bool should_return = false;
#pragma omp parallel for reduction(||                         \
                                   : dominating) reduction(|| \
                                                           : should_return)
    for (size_t i = 0; i < width; ++i)
    {
      if (should_return)
      {
        continue;
      }
      if (p1[i] > p2[i])
      {
        should_return = true;
      }
      else if (p1[i] < p2[i] && !dominating)
      {
        dominating = true;
      }
    }
    return !should_return && dominating;
  }

}

#endif // SDIS_SKYLINE_H
