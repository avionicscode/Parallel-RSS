#ifndef SDIS_PRSS_TIME_H
#define SDIS_PRSS_TIME_H

#define POST_WINDOW_COUNT2 2000

#include <unordered_set>
#include "sdis-cache2.h"
#include "sdis-skyline.h"
#include "sdis-stream.h"

namespace sdistream
{

  static skyline Skyline2;

  template <class IN>
  void par_skyline_update(IN &in, size_t width, size_t window)
  {
    cache2 cache2(width, window); // Tuple cache2.
    size_t count = 0;
    std::unordered_set<double> deal;
    bool display = false;
    auto entries = new cache_entry2[width];           // Index entry buffer.
    auto entries_remove = new cache_entry2[width];    // Index entry buffer of the tuple to remove.
    auto entries_update = new cache_entry2[width];    // Index entry buffer of the non-skyline tuple to update while removing a tuple.
    double index = 0;                                 // Index ID of the incoming tuple.
    auto indexes = new std::set<cache_entry2>[width]; // Dimensional indexes.
    std::set<double> remove;
    double start = 0;
    auto tuple = new double[width]; // Tuple input buffer.

    // Add the first tuple.
    if (!input(in, width, tuple))
    {
      delete[] tuple;
      return;
    }

    index = cache2.put(tuple, true);
#pragma omp parallel for
    for (size_t i = 0; i < width; ++i)
    {
      entries[i].index = index;
      entries[i].value = tuple[i];
      indexes[i].insert(entries[i]);
    }

    Skyline2.add(index);
    start = index;

    // Process incoming tuples.
    while (input(in, width, tuple))
    {
      if (count >= POST_WINDOW_COUNT2)
      {
        break;
      }

      index = cache2.put(tuple);
      if (!display && index - start > window)
      {
        count = 0;
        display = true;
      }
#pragma omp parallel for
      for (size_t i = 0; i < width; ++i)
      {
        entries[i].index = index;
        entries[i].value = tuple[i];
      }
      // Remove expired tuples.
      auto &&expired = cache2.expired();
      if (!expired.empty())
      {
        remove.clear();
        // #pragma omp parallel for
        for (auto &&x : expired)
        {
          remove.insert(x);
        }
        for (auto &&index_remove : remove)
        {
          // Build index entry of the tuple to remove.
          auto &&tuple_remove = cache2.get(index_remove);
          // If tuple does not exist (should not happen), continue with the next one.
          if (!tuple_remove)
          {
            continue;
          }
// Remove expired tuple from all dimensional indexes.
#pragma omp parallel for
          for (size_t i = 0; i < width; ++i)
          {
            entries_remove[i].index = index_remove;
            entries_remove[i].value = tuple_remove[i];
            indexes[i].erase(entries_remove[i]);
          }
          // The expired tuple is a skyline tuple.
          // Can it be optimized ???????
          if (Skyline2.contains(index_remove))
          {
            deal.clear();
            for (auto &&index_update : Skyline2.get(index_remove))
            {
              // Ignore tuples that have already been removed.
              if (!cache2.contains(index_update) || index_update < index_remove)
              {
                continue;
              }
              deal.insert(index_update);
              auto &&tuple_update = cache2.get(index_update); // Green warm.
#pragma omp parallel for
              for (size_t i = 0; i < width; ++i)
              {
                entries_update[i].index = index_update;
                entries_update[i].value = tuple_update[i];
              }
              auto &&lower_bound_dimension = lower_dimension(entries_update, indexes, width);
              auto &&lower_bound_entry = entries_update[lower_bound_dimension];
              auto &&lower_bound_index = indexes[lower_bound_dimension];
              auto &&lower = lower_bound_index.begin();
              bool dominated = false;
              while (lower->value <= lower_bound_entry.value && lower != lower_bound_index.end())
              {
                // If the lower tuple is not in skyline set or is the expired tuple,
                // ignore it.
                if (!Skyline2.contains(lower->index) || lower->index == index_remove)
                {
                  ++lower;
                  continue;
                }
                // If current tuple is dominated ALSO by the lower tuple, do break.
                if (par_dominate<double>(cache2.get(lower->index), tuple_update, width))
                {
                  Skyline2.append(lower->index, index_update);
                  dominated = true;
                  break;
                }
                // If current tuple is not dominated by the lower tuple, even the
                // lower tuple has the same value as the current tuple, the reverse
                // dominance checking is not necessary. Just comment the following
                // code.
                // if (lower->value == lower_bound_entry.value && par_dominate<double>(tuple, cache2.get(lower->index), width)) {
                //  Skyline2.move(lower->index, index);
                //}
                ++lower;
              }
              if (!dominated)
              {
                Skyline2.add(index_update);
              }
              // Dominance tree entries do not respect dimensional indexing order,
              // a local BNL must be applied to fix this problem.
              for (auto &&x : deal)
              {
                if (x != index_update && Skyline2.contains(x) > 0)
                {
                  if (par_dominate<double>(tuple_update, cache2.get(x), width))
                  {
                    Skyline2.move(x, index_update);
                  }
                }
              }
            }
            Skyline2.remove(index_remove);
          }
        }
        cache2.clean();
      }
      // Do lower-bound dominance checking.
      bool dominated = false;
      auto &&lower_bound_dimension = lower_dimension(entries, indexes, width);
      auto &&lower_bound_entry = entries[lower_bound_dimension];
      auto &&lower_bound_index = indexes[lower_bound_dimension];
      auto &&lower = lower_bound_index.begin();
      while (lower->index < index && lower->value <= lower_bound_entry.value && lower != lower_bound_index.end())
      {
        // Only compare the incoming tuple with skyline tuples.
        if (!Skyline2.contains(lower->index))
        {
          ++lower;
          continue;
        }
        // If the incoming tuple is dominated by a lower skyline tuple, do break.
        // The skyline flag of the incoming tuple will be set while adding it
        // to the cache2.
        if (par_dominate<double>(cache2.get(lower->index), tuple, width))
        {
          dominated = true;
          Skyline2.append(lower->index, index);
          break;
        }
        // If the incoming tuple is not dominated by the lower tuple, however the
        // lower tuple has the same value as the current tuple, then do reverse
        // dominance checking.
        if (lower->value == lower_bound_entry.value && par_dominate<double>(tuple, cache2.get(lower->index), width))
        {
          Skyline2.move(lower->index, index);
        }
        ++lower;
      }
      // Do upper-bound dominance checking.
      if (!dominated)
      {
        Skyline2.add(index);
        auto &&upper_bound_dimension = upper_dimension(entries, indexes, width);
        auto &&upper_bound_entry = entries[upper_bound_dimension];
        auto &&upper_bound_index = indexes[upper_bound_dimension];
        auto &&upper_repeat = std::set<cache_entry2>::reverse_iterator(upper_bound_index.lower_bound(upper_bound_entry));
        // For repeating dimensional values.
        while (upper_repeat->index < index && upper_repeat != upper_bound_index.rend())
        {
          if (!Skyline2.contains(upper_repeat->index))
          {
            ++upper_repeat;
            continue;
          }
          if (upper_repeat->value < upper_bound_entry.value)
          {
            break;
          }
          // A tuple with repeat dimensional value is dominated by the incoming
          // tuple.
          if (par_dominate<double>(tuple, cache2.get(upper_repeat->index), width))
          {
            Skyline2.move(upper_repeat->index, index);
          }
          ++upper_repeat;
        }
        // Find all upper skyline tuples that are dominated by the
        // incoming tuple.
        auto &&upper = upper_bound_index.upper_bound(upper_bound_entry);
        while (upper != upper_bound_index.end())
        {
          if (!Skyline2.contains(upper->index))
          {
            ++upper;
            continue;
          }
          if (par_dominate<double>(tuple, cache2.get(upper->index), width))
          {
            Skyline2.move(upper->index, index);
          }
          ++upper;
        }
      }
// Add the incoming tuple to all dimensional indexes.
#pragma omp parallel for
      for (size_t i = 0; i < width; ++i)
      {
        indexes[i].insert(entries[i]);
      }

      if (display)
      {
        ++count;
      }
    }

    delete[] entries;
    delete[] entries_remove;
    delete[] entries_update;
    delete[] indexes;
    delete[] tuple;
  }

}

#endif // SDIS_RSS_TIME_H
