#ifndef SDIS_PRSS_COUNT_H
#define SDIS_PRSS_COUNT_H

#ifndef POST_WINDOW_COUNT
#define POST_WINDOW_COUNT 2000
#endif

#include <unordered_set>
#include "sdis-cache.h"
#include "sdis-skyline.h"
#include "sdis-stream.h"

namespace sdistream
{
  static skyline Skyline2;

  template <class T>
  void par_skyline_update(T &in, size_t width, size_t window)
  {
    cache cache2(width, window); // Tuple cache.
    size_t count2 = 0;
    std::unordered_set<size_t> deal2;
    auto entries2 = new cache_entry[width];           // Index entry buffer.
    auto entries_remove2 = new cache_entry[width];    // Index entry of the tuple to remove.
    auto entries_update2 = new cache_entry[width];    // Index entry of the non-skyline tuple to update while removing a tuple.
    size_t index2 = 0;                                // Index ID of the incoming tuple.
    auto indexes2 = new std::set<cache_entry>[width]; // Dimensional indexes.
    auto tuple2 = new double[width];                  // Tuple input buffer.
    // Add the first tuple.
    if (!input(in, width, tuple2))
    {
      delete[] tuple2;
      return;
    }

#pragma omp parallel for
    for (size_t i = 0; i < width; ++i)
    {
      entries2[i].index = index2;
      entries2[i].value = tuple2[i];
      indexes2[i].insert(entries2[i]);
    }

    cache2.put(tuple2, true);
    Skyline2.add(index2);
    ++index2;

    // Process incoming tuples.
    while (input(in, width, tuple2))
    {

      if (count2 >= POST_WINDOW_COUNT)
      {
        break;
      }

#pragma omp parallel for
      for (size_t i = 0; i < width; ++i)
      {
        entries2[i].index = index2;
        entries2[i].value = tuple2[i];
      }
      // Remove the expired tuple.
      if (index2 >= window)
      {
        ++count2;
        // Build index entry of the tuple to remove.
        auto &&index_remove2 = index2 - window;
        auto &&tuple_remove2 = cache2.get(index_remove2);
#pragma omp parallel for
        for (size_t i = 0; i < width; ++i)
        {
          entries_remove2[i].index = index_remove2;
          entries_remove2[i].value = tuple_remove2[i];
        }
        // The expired tuple is a skyline tuple.
        if (cache2.skyline(index_remove2))
        {
          deal2.clear();
          for (auto &&index_update2 : Skyline2.get(index_remove2))
          {
            // Ignore tuples that have already been removed.
            if (index_update2 < index_remove2)
            {
              continue;
            }
            deal2.insert(index_update2);
            auto &&tuple_update2 = cache2.get(index_update2); // Green warm.
#pragma omp parallel for
            for (size_t i = 0; i < width; ++i)
            {
              entries_update2[i].index = index_update2;
              entries_update2[i].value = tuple_update2[i];
            }
            auto &&lower_bound_dimension2 = par_lower_dimension(entries_update2, indexes2, width);
            auto &&lower_bound_entry2 = entries_update2[lower_bound_dimension2];
            auto &&lower_bound_index2 = indexes2[lower_bound_dimension2];
            auto &&lower2 = lower_bound_index2.begin();
            bool dominated2 = false;

            while (lower2->value <= lower_bound_entry2.value && lower2 != lower_bound_index2.end())
            {
              // If the lower tuple is not in skyline set or is the expired tuple,
              // ignore it.

              if (!cache2.skyline(lower2->index) || lower2->index == index_remove2)
              {
                ++lower2;
                continue;
              }
              // If current tuple is dominated ALSO by the lower tuple, do break.
              if (par_dominate<double>(cache2.get(lower2->index), tuple_update2, width))
              {
                Skyline2.append(lower2->index, index_update2);
                dominated2 = true;
                break;
              }
              // If current tuple is not dominated by the lower tuple, even the
              // lower tuple has the same value as the current tuple, the reverse
              // dominance checking is not necessary. Just comment the following
              // useless code.
              // if (lower->value == lower_bound_entry.value && dominate<double>(tuple, cache.get(lower->index), width)) {
              //  Skyline.move(lower->index, index);
              //}
              ++lower2;
            }
            if (!dominated2)
            {
              cache2.skyline(index_update2) = true;
              Skyline2.add(index_update2);
            }
            // Dominance tree entries do not respect dimensional indexing order,
            // a local BNL must be applied to fix this problem.
            for (auto &&x : deal2)
            {
              if (x != index_update2 && cache2.skyline(x))
              {
                if (par_dominate<double>(tuple_update2, cache2.get(x), width))
                {
                  cache2.skyline(x) = false;
                  Skyline2.move(x, index_update2);
                }
              }
            }
          }
          cache2.skyline(index_remove2) = false; // Not really necessary.
          Skyline2.remove(index_remove2);
        }
// Remove expired tuple from all dimensional indexes.
#pragma omp parallel for
        for (size_t i = 0; i < width; ++i)
        {
          indexes2[i].erase(entries_remove2[i]);
        }
      }
      // Do lower-bound dominance checking.
      bool dominated2 = false;
      auto &&lower_bound_dimension2 = par_lower_dimension(entries2, indexes2, width);
      auto &&lower_bound_entry2 = entries2[lower_bound_dimension2];
      auto &&lower_bound_index2 = indexes2[lower_bound_dimension2];
      auto &&lower2 = lower_bound_index2.begin();
      while (lower2->value <= lower_bound_entry2.value && lower2 != lower_bound_index2.end())
      {
        // Only compare the incoming tuple with skyline tuples.
        if (!cache2.skyline(lower2->index))
        {
          ++lower2;
          continue;
        }
        // If the incoming tuple is dominated by a lower skyline tuple, do break.
        // The skyline flag of the incoming tuple will be set while adding it
        // to the cache.
        if (par_dominate<double>(cache2.get(lower2->index), tuple2, width))
        {
          dominated2 = true;
          Skyline2.append(lower2->index, index2);
          break;
        }
        // If the incoming tuple is not dominated by the lower tuple, however the
        // lower tuple has the same value as the current tuple, then do reverse
        // dominance checking.
        if (lower2->value == lower_bound_entry2.value && par_dominate<double>(tuple2, cache2.get(lower2->index), width))
        {
          cache2.skyline(lower2->index) = false;
          Skyline2.move(lower2->index, index2);
        }
        ++lower2;
      }
      // Do upper-bound dominance checking.
      if (!dominated2)
      {
        Skyline2.add(index2);
        auto &&upper_bound_dimension2 = par_upper_dimension(entries2, indexes2, width);
        auto &&upper_bound_entry2 = entries2[upper_bound_dimension2];
        auto &&upper_bound_index2 = indexes2[upper_bound_dimension2];
        auto &&upper_repeat2 = std::set<cache_entry>::reverse_iterator(upper_bound_index2.lower_bound(upper_bound_entry2));
        // For repeating dimensional values.
        while (upper_repeat2 != upper_bound_index2.rend())
        {
          if (!cache2.skyline(upper_repeat2->index))
          {
            ++upper_repeat2;
            continue;
          }
          if (upper_repeat2->value < upper_bound_entry2.value)
          {
            break;
          }
          // A tuple with repeat dimensional value is dominated by the incoming
          // tuple.
          if (par_dominate<double>(tuple2, cache2.get(upper_repeat2->index), width))
          {
            cache2.skyline(upper_repeat2->index) = false;
            Skyline2.move(upper_repeat2->index, index2);
          }
          ++upper_repeat2;
        }
        // Find all upper skyline tuples that are dominated by the
        // incoming tuple.
        auto &&upper2 = upper_bound_index2.upper_bound(upper_bound_entry2); // this set can be parallel.

        while (upper2 != upper_bound_index2.end())
        {
          if (!Skyline2.contains(upper2->index))
          {
            // if (!cache.skyline(upper->index)) {
            ++upper2;
            continue;
          }
          if (par_dominate<double>(tuple2, cache2.get(upper2->index), width))
          {
            cache2.skyline(upper2->index) = false;
            Skyline2.move(upper2->index, index2);
          }
          ++upper2;
        }
      }
// Add the incoming tuple to all dimensional indexes.
#pragma omp parallel for
      for (size_t i = 0; i < width; ++i)
      {
        indexes2[i].insert(entries2[i]);
      }
      // Finally, replace the expired tuple by the incoming tuple.
      cache2.put(tuple2, !dominated2);
      ++index2;
    }

    delete[] entries2;
    delete[] entries_remove2;
    delete[] entries_update2;
    delete[] indexes2;
    delete[] tuple2;
  }
}

#endif // SDIS_PRSS_COUNT_H