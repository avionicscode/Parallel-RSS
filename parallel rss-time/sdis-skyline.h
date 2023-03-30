#ifndef SDIS_SKYLINE_H
#define SDIS_SKYLINE_H

#ifndef SLICE
#define SLICE 32
#endif

#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>

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
    auto add(const double &) -> skyline &;
    // The first skyline point dominates the second point.
    auto append(const double &, const double &) -> skyline &;
    // Check whether a point is skyline point.
    auto contains(const double &) -> bool;
    // Get all dominated points of a skyline point.
    auto get(const double &) -> std::vector<double> &;
    // Move the first skyline point to the second skyline point.
    auto move(const double &, const double &) -> skyline &;
    // Remove a skyline point.
    auto remove(const double &) -> skyline &;
    // The number of skyline points.
    auto size() -> size_t;

  private:
    static auto slice_(double) -> size_t;
    size_t count_ = 0;
    std::vector<double> empty_;
    std::array<std::unordered_map<double, std::vector<double>>, SLICE> tree_;
  };

  template <class V>
  auto dominate(V *row1, V *row2, size_t width) -> bool
  {
    ++skyline::DT;
    V *p1 = row1;
    V *p2 = row2;
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
