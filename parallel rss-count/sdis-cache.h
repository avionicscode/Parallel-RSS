#ifndef SDIS_CACHE_H
#define SDIS_CACHE_H

#include <iostream>
#include <set>
#include <omp.h>

namespace sdistream
{

  struct cache_entry
  {
    size_t index = 0;
    double value = 0;
    explicit cache_entry(size_t);
    cache_entry(size_t, double);

    inline auto set(size_t n) -> cache_entry &
    {
      index = n;
      return *this;
    }
    cache_entry() = default;
  };

  auto operator==(const cache_entry &, const cache_entry &) -> bool;
  auto operator<(const cache_entry &, const cache_entry &) -> bool;
  auto operator<<(std::ostream &, const cache_entry &) -> std::ostream &;

  auto estimate(const cache_entry &, const std::set<cache_entry> &) -> double;
  auto lower_dimension(const cache_entry *, const std::set<cache_entry> *, size_t) -> size_t;
  auto upper_dimension(const cache_entry *, const std::set<cache_entry> *, size_t) -> size_t;
  auto par_lower_dimension(const cache_entry *, const std::set<cache_entry> *, size_t) -> size_t;
  auto par_upper_dimension(const cache_entry *, const std::set<cache_entry> *, size_t) -> size_t;

  class cache
  {
  public:
    cache() = default;
    cache(size_t, size_t);
    virtual ~cache();
    auto get(size_t) -> double *;
    auto put(double *) -> double *;
    auto put(double *, bool) -> double *;
    auto skyline(size_t) -> bool &;

  private:
    double *cache_ = nullptr;
    size_t count_ = 0;
    bool *skyline_ = nullptr;
    size_t width_ = 0;
    size_t window_ = 0;
  };

}

#endif // SDIS_CACHE_H
