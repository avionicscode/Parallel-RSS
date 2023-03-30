#include <cmath>
#include <cstring>
#include "sdis-cache.h"

namespace sdistream
{

  cache_entry::cache_entry(size_t index) : index(index)
  {
  }

  cache_entry::cache_entry(size_t index, double value) : index(index), value(value)
  {
  }

  auto operator==(const cache_entry &row1, const cache_entry &row2) -> bool
  {
    return row1.index == row2.index;
  }

  auto operator<(const cache_entry &row1, const cache_entry &row2) -> bool
  {
    if (row1.value == row2.value)
    {
      return row1.index < row2.index;
    }
    return row1.value < row2.value;
  }

  auto operator<<(std::ostream &out, const cache_entry &row) -> std::ostream &
  {
    out << row.value << ":" << row.index;
    return out;
  }

  auto estimate(const cache_entry &ent, const std::set<cache_entry> &dim) -> double
  {
    if (dim.empty())
    {
      return 0;
    }
    auto &&first = *(dim.begin());
    auto &&last = *(dim.rbegin());
    if (first == last)
    {
      return 1;
    }
    if (ent < first)
    {
      return 0;
    }
    if (last < ent)
    {
      return 1;
    }
    return 1.0 * std::abs(ent.value - first.value) / std::abs(last.value - first.value);
  }

  auto lower_dimension(const cache_entry *entries, const std::set<cache_entry> *indexes, size_t width) -> size_t
  {
    size_t d = 0;
    double lower = 1;
    for (size_t i = 0; i < width; ++i)
    {
      double est = estimate(entries[i], indexes[i]);
      if (est == 0)
      {
        return i;
      }
      else
      {
        if (est < lower)
        {
          lower = est;
          d = i;
        }
      }
    }
    return d;
  }

  auto upper_dimension(const cache_entry *entries, const std::set<cache_entry> *indexes, size_t width) -> size_t
  {
    size_t d = 0;
    double upper = 0;
    for (size_t i = 0; i < width; ++i)
    {
      double est = estimate(entries[i], indexes[i]);
      if (est == 1)
      {
        return i;
      }
      else
      {
        if (est > upper)
        {
          upper = est;
          d = i;
        }
      }
    }
    return d;
  }

  auto par_lower_dimension(const cache_entry *entries, const std::set<cache_entry> *indexes, size_t width) -> size_t
  {
    size_t d = 0;
    double lower = 1;

#pragma omp parallel reduction(min \
                               : d) shared(entries, indexes, width, lower)
    {
      size_t private_d = 0;
      double private_lower = 1;

#pragma omp for
      for (size_t i = 0; i < width; ++i)
      {
        double est = estimate(entries[i], indexes[i]);
        if (est == 0)
        {
          private_d = i;
#pragma omp cancel for
        }
        else
        {
          if (est < private_lower)
          {
            private_lower = est;
            private_d = i;
          }
        }
      }

#pragma omp critical
      {
        if (private_lower < lower)
        {
          lower = private_lower;
          d = private_d;
        }
      }
    }

    return d;
  }

  auto par_upper_dimension(const cache_entry *entries, const std::set<cache_entry> *indexes, size_t width) -> size_t
  {
    size_t d = 0;
    double upper = 0;
    for (size_t i = 0; i < width; ++i)
    {
      double est = estimate(entries[i], indexes[i]);
      if (est == 1)
      {
        return i;
      }
      else
      {
        if (est > upper)
        {
          upper = est;
          d = i;
        }
      }
    }
    return d;
  }
  // -----------------------------------------------------------------------------------
  cache::cache(size_t width, size_t window) : count_(0), width_(width), window_(window)
  {
    cache_ = new double[width_ * window_];
    skyline_ = new bool[window_];
  }

  cache::~cache()
  {
    delete[] cache_;
    delete[] skyline_;
  }

  auto cache::get(size_t index) -> double *
  {
    if (index > count_)
    {
      return nullptr;
    }
    return &cache_[(index % window_) * width_];
  }

  auto cache::put(double *buffer) -> double *
  {
    size_t index = count_ % window_;
    double *base = &cache_[index * width_];
    std::memcpy(base, buffer, sizeof(double) * width_);
    skyline_[index] = false;
    ++count_;
    return base;
  }

  auto cache::put(double *buffer, bool skyline) -> double *
  {
    size_t index = count_ % window_;
    double *base = &cache_[index * width_];
    std::memcpy(base, buffer, sizeof(double) * width_);
    skyline_[index] = skyline;
    ++count_;
    return base;
  }

  auto cache::skyline(size_t index) -> bool &
  {
    return skyline_[index % window_];
  }

}
