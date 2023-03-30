
#include "sdis-cache.h"

namespace sdistream
{

  cache_entry::cache_entry(double index) : index(index)
  {
  }

  cache_entry::cache_entry(double index, double value) : index(index), value(value)
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

  //--------------------------------------------------------

  size_t cache::zero_ = 0;

  auto cache::timestamp() -> double
  {
    struct timeval t
    {
    };
    gettimeofday(&t, (struct timezone *)nullptr);
    return (double)(t.tv_sec - zero_) + t.tv_usec / 1000000.0;
  }

  cache::cache(size_t width, size_t window) : free_(CACHE), skyline_(width), stamp_(width + 1), width_(width), width2_(width + 2), window_(window)
  {
    cache_ = new double[width2_ * CACHE]; // Tuple + Skyline flag (-1/1) + Timestamp
    for (size_t i = 0; i < CACHE; ++i)
    {
      free_[i] = &cache_[width2_ * i];
    }
    struct timeval t
    {
    };
    gettimeofday(&t, (struct timezone *)nullptr);
    zero_ = t.tv_sec;
  }

  cache::~cache()
  {
    delete[] cache_;
  }

  void cache::clean()
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

  auto cache::contains(double stamp) -> bool
  {
    return index_[slice_(stamp)].count(stamp);
  }

  auto cache::expired() -> std::vector<double> &
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

  auto cache::get(double stamp) -> double *
  {
    auto &&it = index_[slice_(stamp)].find(stamp);
    if (it == index_[slice_(stamp)].end())
    {
      return nullptr;
    }
    return it->second;
  }

  auto cache::put(double *buffer) -> double
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

  auto cache::put(double *buffer, bool skyline) -> double
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

  auto cache::size() -> size_t
  {
    return CACHE - free_.size();
  }

  auto cache::slice_(double index) -> size_t
  {
    return (int)index % BLOCK;
  }

}
