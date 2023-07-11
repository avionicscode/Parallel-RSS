#include "sdis-skyline.h"

namespace sdistream
{

  size_t skyline::DT = 0;

  auto operator<<(std::ostream &out, const skyline &skyline) -> std::ostream &
  {
    for (auto &&tree : skyline.tree_)
    {
      for (auto &&s : tree)
      {
        out << s.first << std::endl;
      }
    }
    return out;
  }

  auto skyline::add(const double &s) -> skyline &
  {
    auto &&it = tree_[slice_(s)].find(s);
    if (it != tree_[slice_(s)].end())
    {
      return *this;
    }
    std::vector<double> v;
    tree_[slice_(s)].insert(std::make_pair(s, v));
    ++count_;
    return *this;
  }

  auto skyline::append(const double &s, const double &p) -> skyline &
  {
    auto &&it = tree_[slice_(s)].find(s);
    if (it == tree_[slice_(s)].end())
    {
      std::vector<double> v;
      v.push_back(p);
      tree_[slice_(s)].insert(std::make_pair(s, v));
      return *this;
    }
    it->second.push_back(p);
    return *this;
  }

  auto skyline::contains(const double &p) -> bool
  {
    return tree_[slice_(p)].count(p);
  }

  auto skyline::get(const double &s) -> std::vector<double> &
  {
    auto &&it = tree_[slice_(s)].find(s);
    if (it == tree_[slice_(s)].end())
    {
      return empty_;
    }
    return it->second;
  }

  auto skyline::move(const double &s1, const double &s2) -> skyline &
  {
    auto &&it1 = tree_[slice_(s1)].find(s1);
    if (it1 == tree_[slice_(s1)].end())
    {
      return *this;
    }
    auto &&it2 = tree_[slice_(s2)].find(s2);
    if (it2 == tree_[slice_(s2)].end())
    {
      return *this;
    }
    it2->second.push_back(s1);
    it2->second.insert(it2->second.end(), it1->second.begin(), it1->second.end());
    tree_[slice_(s1)].erase(it1);
    --count_;
    return *this;
  }

  auto skyline::remove(const double &s) -> skyline &
  {
    tree_[slice_(s)].erase(s);
    --count_;
    return *this;
  }

  auto skyline::size() -> size_t
  {
    return count_;
  }

  auto skyline::slice_(double index) -> size_t
  {
    return (int)index % SLICE; // produces the remainder of index / SLICE.
  }

}
