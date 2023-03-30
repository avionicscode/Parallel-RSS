#ifndef SDIS_CACHE2_H
#define SDIS_CACHE2_H

#define CACHE2 10000

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
  struct cache_entry2
  {
    double index = 0;
    double value = 0;
    explicit cache_entry2(double);
    cache_entry2(double, double);
    inline auto set(double n) -> cache_entry2 &
    {
      index = n;
      return *this;
    }
    cache_entry2() = default;
  };

  auto operator==(const cache_entry2 &, const cache_entry2 &) -> bool;
  auto operator<(const cache_entry2 &, const cache_entry2 &) -> bool;
  auto operator<<(std::ostream &, const cache_entry2 &) -> std::ostream &;

  auto estimate(const cache_entry2 &, const std::set<cache_entry2> &) -> double;
  auto lower_dimension(const cache_entry2 *, const std::set<cache_entry2> *, size_t) -> size_t;
  auto upper_dimension(const cache_entry2 *, const std::set<cache_entry2> *, size_t) -> size_t;

  class cache2
  {
  public:
    static auto timestamp() -> double;
    cache2() = default;
    cache2(size_t, size_t);
    virtual ~cache2();
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
