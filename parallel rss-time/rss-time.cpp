#include <fstream>
#include <iostream>
#include "rss-time.h"
#include "prss-time.h"
#include <omp.h>

using namespace sdistream;

auto run_skyline(const char *name, size_t dimensionality, size_t window, const char *stream) -> bool
{
  std::cerr << " Sequential Running..." << std::endl;
  if (stream)
  {
    std::ifstream in(stream);
    if (!in.good())
    {
      std::cerr << "Cannot open stream " << stream << std::endl;
      return false;
    }

    double startTime = omp_get_wtime();
    skyline_update<std::ifstream>(in, dimensionality, window);
    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime;
    std::cout << "Sequential Skyline: " << secsElapsed << std::endl;

    in.close();
  }
  else
  {
    skyline_update<std::istream>(std::cin, dimensionality, window);
  }
  return true;
}

auto par_run_skyline(const char *name, size_t dimensionality, size_t window, const char *stream) -> bool
{
  std::cerr << " Parallel  Running..." << std::endl;
  if (stream)
  {
    std::ifstream in(stream);
    if (!in.good())
    {
      std::cerr << "Cannot open stream " << stream << std::endl;
      return false;
    }

    double startTime = omp_get_wtime();
    par_skyline_update<std::ifstream>(in, dimensionality, window);
    double stopTime = omp_get_wtime();
    double secsElapsed = stopTime - startTime;
    std::cout << "Parallel  Skyline: " << secsElapsed << std::endl;

    in.close();
  }
  else
  {
    par_skyline_update<std::istream>(std::cin, dimensionality, window);
  }
  return true;
}

auto main(int argc, char **argv) -> int
{
  if (argc < 3)
  {
    std::cout << "Usage: rss-time DIMENSIONALITY WINDOW [STREAM]" << std::endl;
    return 0;
  }
  size_t dimensionality = strtoul(argv[1], nullptr, 10);
  size_t window = strtoul(argv[2], nullptr, 10);
  const char *stream = argc > 3 ? argv[3] : nullptr;
  run_skyline("RSS-TIME", dimensionality, window, stream);
  par_run_skyline("RSS-TIME", dimensionality, window, stream);
  return 0;
}