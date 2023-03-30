#include "sdis-cache2.h"

namespace sdistream
{

  cache_entry2::cache_entry2(double index) : index(index)
  {
  }

  cache_entry2::cache_entry2(double index, double value) : index(index), value(value)
  {
  }

  auto operator==(const cache_entry2 &row1, const cache_entry2 &row2) -> bool
  {
    return row1.index == row2.index;
  }

  auto operator<(const cache_entry2 &row1, const cache_entry2 &row2) -> bool
  {
    if (row1.value == row2.value)
    {
      return row1.index < row2.index;
    }
    return row1.value < row2.value;
  }

  auto operator<<(std::ostream &out, const cache_entry2 &row) -> std::ostream &
  {
    out << row.value << ":" << row.index;
    return out;
  }

  auto estimate(const cache_entry2 &ent, const std::set<cache_entry2> &dim) -> double
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

  auto lower_dimension(const cache_entry2 *entries, const std::set<cache_entry2> *indexes, size_t width) -> size_t
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

  auto upper_dimension(const cache_entry2 *entries, const std::set<cache_entry2> *indexes, size_t width) -> size_t
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

  //--------------------------------------------------------

  size_t cache2::zero_ = 0;

  auto cache2::timestamp() -> double
  {
    struct timeval t
    {
    };
    gettimeofday(&t, (struct timezone *)nullptr);
    return (double)(t.tv_sec - zero_) + t.tv_usec / 1000000.0;
  }

  cache2::cache2(size_t width, size_t window) : free_(CACHE2), skyline_(width), stamp_(width + 1), width_(width), width2_(width + 2), window_(window)
  {
    cache_ = new double[width2_ * CACHE2]; // Tuple + Skyline flag (-1/1) + Timestamp
    for (size_t i = 0; i < CACHE2; ++i)
    {
      free_[i] = &cache_[width2_ * i];
    }
    struct timeval t
    {
    };
    gettimeofday(&t, (struct timezone *)nullptr);
    zero_ = t.tv_sec;
  }

  cache2::~cache2()
  {
    delete[] cache_;
  }

  void cache2::clean()
  {
    for (auto &&it = list_.begin(); it != list_.end();)
    {
      auto &&early = (*it)[stamp_];
      if (latest_ - early > window_)
      {
        free_.push_back(*it);
        index_[slice_(early)].erase(early);
        it = list_.erase(it);
      }
      else
      {
        break;
      }
    }
  }

  auto cache2::contains(double stamp) -> bool
  {
    return index_[slice_(stamp)].count(stamp);
  }

  auto cache2::expired() -> std::vector<double> &
  {
    expired_.clear();
    for (auto &&it = list_.begin(); it != list_.end(); ++it)
    {
      auto &&early = (*it)[stamp_];
      if (latest_ - early > window_)
      {
        expired_.push_back(early);
      }
      else
      {
        break;
      }
    }
    return expired_;
  }

  auto cache2::get(double stamp) -> double *
  {
    auto &&it = index_[slice_(stamp)].find(stamp);
    if (it == index_[slice_(stamp)].end())
    {
      return nullptr;
    }
    return it->second;
  }

  auto cache2::put(double *buffer) -> double
  {
    latest_ = timestamp();
    auto &&row = free_.back();
    std::memcpy(row, buffer, sizeof(double) * width_);
    free_.pop_back();
    index_[slice_(latest_)].insert(std::make_pair(latest_, row));
    list_.push_back(row);
    row[stamp_] = latest_;
    ++count_;
    return latest_;
  }

  auto cache2::put(double *buffer, bool skyline) -> double
  {
    latest_ = timestamp();
    auto &&row = free_.back();
    std::memcpy(row, buffer, sizeof(double) * width_);
    free_.pop_back();
    index_[slice_(latest_)].insert(std::make_pair(latest_, row));
    list_.push_back(row);
    row[stamp_] = latest_;
    ++count_;
    return latest_;
  }

  auto cache2::size() -> size_t
  {
    return CACHE2 - free_.size();
  }

  auto cache2::slice_(double index) -> size_t
  {
    return (int)index % BLOCK;
  }

}
