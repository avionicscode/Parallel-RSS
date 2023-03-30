#ifndef SDIS_CACHE_H
#define SDIS_CACHE_H

#define CACHE 10000

#define BLOCK 16

#include <iostream>
#include <set>
#include <array>
#include <list>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <ctime>
#include <cstring>
#include <omp.h>

namespace sdistream
{
  struct cache_entry
  {
    double index = 0;
    double value = 0;
    explicit cache_entry(double);
    cache_entry(double, double);
    inline auto set(double n) -> cache_entry &
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

  class cache
  {
  public:
    static auto timestamp() -> double;
    cache() = default;
    cache(size_t, size_t);
    virtual ~cache();
    void clean();
    auto contains(double) -> bool;
    auto expired() -> std::vector<double> &;
    auto get(double) -> double *;
    auto put(double *) -> double;
    auto put(double *, bool) -> double;
    auto size() -> size_t;

  private:
    static auto slice_(double) -> size_t;
    static size_t zero_;
    double *cache_ = nullptr;
    size_t count_ = 0;
    std::vector<double> expired_;
    std::vector<double *> free_;
    std::array<std::unordered_map<double, double *>, BLOCK> index_;
    double latest_ = 0;
    std::list<double *> list_;
    size_t skyline_ = 0;
    size_t stamp_ = 0;
    size_t width_ = 0;
    size_t width2_ = 0;
    double window_ = 0;
  };

}

#endif // SDIS_CACHE_H
