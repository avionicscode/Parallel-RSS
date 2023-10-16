/*
Test parallel RSS with big window size.
*/
#include <iostream>
#include <unordered_set>
#include <array>
#include <unordered_map>
#include <vector>
#include <cstring>
#include <set>
#include <cmath>
#include <omp.h>

#define SLICE 32

const size_t width = 8;
const size_t window = 100000;
const size_t datasize = 140000;

double data[datasize][width];

// just for demonstration purpose, PRSS uses data from csv.
void InitData()
{
    for (size_t i = 0; i < datasize; i++)
    {
        for (size_t j = 0; j < width; j++)
        {

            data[i][j] = 2 * i + 1;
        }
    }
}

size_t skyline_count_ = 0;
std::vector<size_t> empty_;
std::array<std::unordered_map<size_t, std::vector<size_t>>, SLICE> tree_;

size_t cache_count_ = 0;
double *cache_ = new double[width * window];
bool *skyline_ = new bool[window];

static size_t skyline_slice_(size_t index)
{
    return (int)index % SLICE;
}

static void skyline_add(const size_t &s)
{
    auto &&it = tree_[skyline_slice_(s)].find(s);
    if (it == tree_[skyline_slice_(s)].end())
    {
        std::vector<size_t> v;
        tree_[skyline_slice_(s)].insert(std::make_pair(s, v));
        ++skyline_count_;
    }
}

static void skyline_append(const size_t &s, const size_t &p)
{
    auto &&it = tree_[skyline_slice_(s)].find(s);
    if (it == tree_[skyline_slice_(s)].end())
    {
        std::vector<size_t> v;
        v.push_back(p);
        tree_[skyline_slice_(s)].insert(std::make_pair(s, v));
    }
    else
    {
        it->second.push_back(p);
    }
}

static bool skyline_contains(const size_t &p)
{
    return tree_[skyline_slice_(p)].count(p);
}

static std::vector<size_t> &skyline_get(const size_t &s)
{
    auto &&it = tree_[skyline_slice_(s)].find(s);
    if (it == tree_[skyline_slice_(s)].end())
    {
        return empty_;
    }
    return it->second;
}

static void skyline_move(const size_t &s1, const size_t &s2)
{
    auto &&it1 = tree_[skyline_slice_(s1)].find(s1);
    if (it1 != tree_[skyline_slice_(s1)].end())
    {
        auto &&it2 = tree_[skyline_slice_(s2)].find(s2);
        if (it2 != tree_[skyline_slice_(s2)].end())
        {
            it2->second.push_back(s1);
            it2->second.insert(it2->second.end(), it1->second.begin(), it1->second.end());
            tree_[skyline_slice_(s1)].erase(it1);
            --skyline_count_;
        }
    }
}

static void skyline_remove(const size_t &s)
{
    tree_[skyline_slice_(s)].erase(s);
    --skyline_count_;
}

static size_t skyline_size()
{
    return skyline_count_;
}

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

void destroyCache()
{
    delete[] cache_;
    delete[] skyline_;
}

double *cache_get(size_t index)
{
    if (index > cache_count_)
    {
        return nullptr;
    }
    return &cache_[(index % window) * width];
}

double *cache_put(double *buffer, bool skyline)
{
    size_t index = cache_count_ % window;
    double *base = &cache_[index * width];
    std::memcpy(base, buffer, sizeof(double) * width);
    skyline_[index] = skyline;
    ++cache_count_;
    return base;
}

bool &cache_skyline(size_t index)
{
    return skyline_[index % window];
}

template <class T>
auto dominate(T *row1, T *row2, size_t width) -> bool
{
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

void run_skyline_for_window(int start, int end)
{
    std::unordered_set<size_t> deal;
    auto entries = new cache_entry[width];
    auto entries_remove = new cache_entry[width];
    auto entries_update = new cache_entry[width];
    size_t index = start;
    auto indexes = new std::set<cache_entry>[width];
    auto tuple = new double[width]; // tuple = data[j][i]

    for (size_t i = 0; i < width; ++i)
    {
        tuple[i] = data[index][i];
        entries[i].index = index;
        entries[i].value = data[index][i];
        indexes[i].insert(entries[i]);
    }

    cache_put(tuple, true);
    skyline_add(index);
    ++index;
    for (; index < end; ++index)
    {

        for (size_t i = 0; i < width; ++i)
        {
            tuple[i] = data[index][i];
            entries[i].index = index;
            entries[i].value = data[index][i];
        }

        if (index >= window)
        {
            auto &&index_remove = index - window;
            auto &&tuple_remove = cache_get(index_remove);
            for (size_t i = 0; i < width; ++i)
            {
                entries_remove[i].index = index_remove;
                entries_remove[i].value = tuple_remove[i];
            }

            if (cache_skyline(index_remove))
            {
                deal.clear();
                for (auto &&index_update : skyline_get(index_remove))
                {
                    if (index_update < index_remove)
                    {
                        continue;
                    }
                    deal.insert(index_update);
                    auto &&tuple_update = cache_get(index_update);
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
                        if (!cache_skyline(lower->index) || lower->index == index_remove)
                        {
                            ++lower;
                            continue;
                        }

                        if (dominate<double>(cache_get(lower->index), tuple_update, width))
                        {
                            skyline_append(lower->index, index_update);
                            dominated = true;
                            break;
                        }

                        ++lower;
                    }
                    if (!dominated)
                    {
                        cache_skyline(index_update) = true;
                        skyline_add(index_update);
                    }

                    for (auto &&x : deal)
                    {
                        if (x != index_update && cache_skyline(x))
                        {
                            if (dominate<double>(tuple_update, cache_get(x), width))
                            {
                                cache_skyline(x) = false;
                                skyline_move(x, index_update);
                            }
                        }
                    }
                }
                cache_skyline(index_remove) = false;
                skyline_remove(index_remove);
            }

            for (size_t i = 0; i < width; ++i)
            {
                indexes[i].erase(entries_remove[i]);
            }
        }

        bool dominated = false;

        auto &&lower_bound_dimension = lower_dimension(entries, indexes, width);
        auto &&lower_bound_entry = entries[lower_bound_dimension];
        auto &&lower_bound_index = indexes[lower_bound_dimension];
        auto &&lower = lower_bound_index.begin();
        while (lower->value <= lower_bound_entry.value && lower != lower_bound_index.end())
        {

            if (!cache_skyline(lower->index))
            {
                ++lower;
                continue;
            }

            if (dominate<double>(cache_get(lower->index), tuple, width))
            {
                dominated = true;
                skyline_append(lower->index, index);
                break;
            }

            if (lower->value == lower_bound_entry.value && dominate<double>(tuple, cache_get(lower->index), width))
            {
                cache_skyline(lower->index) = false;
                skyline_move(lower->index, index);
            }
            ++lower;
        }

        if (!dominated)
        {
            skyline_add(index);
            auto &&upper_bound_dimension = upper_dimension(entries, indexes, width);
            auto &&upper_bound_entry = entries[upper_bound_dimension];
            auto &&upper_bound_index = indexes[upper_bound_dimension];
            auto &&upper_repeat = std::set<cache_entry>::reverse_iterator(upper_bound_index.lower_bound(upper_bound_entry));

            while (upper_repeat != upper_bound_index.rend())
            {
                if (!cache_skyline(upper_repeat->index))
                {
                    ++upper_repeat;
                    continue;
                }
                if (upper_repeat->value < upper_bound_entry.value)
                {
                    break;
                }

                if (dominate<double>(tuple, cache_get(upper_repeat->index), width))
                {
                    cache_skyline(upper_repeat->index) = false;
                    skyline_move(upper_repeat->index, index);
                }
                ++upper_repeat;
            }

            auto &&upper = upper_bound_index.upper_bound(upper_bound_entry);

            while (upper != upper_bound_index.end())
            {
                if (!skyline_contains(upper->index))
                {
                    ++upper;
                    continue;
                }
                if (dominate<double>(tuple, cache_get(upper->index), width))
                {
                    cache_skyline(upper->index) = false;
                    skyline_move(upper->index, index);
                }
                ++upper;
            }
        }

        for (size_t i = 0; i < width; ++i)
        {
            indexes[i].insert(entries[i]);
        }

        cache_put(tuple, !dominated);
        ++index;
    }

    delete[] entries;
    delete[] entries_remove;
    delete[] entries_update;
    delete[] indexes;
    delete[] tuple;
}

int main()
{
    InitData();
    int num_cores = 4; // Change this to the number of cores you want to use.

    double startTime = omp_get_wtime();

#pragma omp parallel num_threads(num_cores)
    {
        int tid = omp_get_thread_num();
        int window_size = window / num_cores;
        int start = tid * window_size;
        int end = start + window_size;
        if (tid == num_cores - 1)
        {
            // Last thread should handle the remaining part of the window.
            end = window;
        }
        run_skyline_for_window(start, end);
    }

    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime;
    std::cout << "Time Skyline: " << secsElapsed << std::endl;

    destroyCache();

    return 0;
}