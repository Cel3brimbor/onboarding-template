#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstddef>
#include <functional>
#include <thread>
#include <vector>

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid {

//attributes
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<double> data_;

public:
  //constructor
  Grid(std::size_t rows, std::size_t cols)
  {
    rows_ = rows;
    cols_ = cols;
    data_.assign(rows * cols, 0.0);
  }

  //getters
  double& operator()(std::size_t i, std::size_t j)
  {
    return data_[i*cols_ + j];
  }
  double  operator()(std::size_t i, std::size_t j) const
  {
    return data_[i*cols_ + j];
  }

  std::size_t get_row_size() const
  {
    return rows_;
  }

  std::size_t get_col_size() const
  {
    return cols_;
  }
};  

// Helper function to process a row slice of the grid
void process_rows(std::size_t row_begin, std::size_t row_end, const Grid& old_grid, Grid& new_grid)
{
  const std::size_t rows = old_grid.get_row_size();
  const std::size_t cols = old_grid.get_col_size();

  for(std::size_t i = row_begin; i < row_end; i++)
  {
    if(i == 0 || i == rows - 1)
    {
      for(std::size_t j = 0; j < cols; j++)
      {
        new_grid(i, j) = old_grid(i, j);
      }
    }
    else 
    {
      new_grid(i, 0) = old_grid(i, 0);
      new_grid(i, cols - 1) = old_grid(i, cols - 1);

      for(std::size_t j = 1; j < cols - 1; j++)
      {
        new_grid(i, j) = 0.5 * old_grid(i, j) +
                        0.125 * (old_grid(i - 1, j) + old_grid(i + 1, j) +
                                 old_grid(i, j - 1) + old_grid(i, j + 1));
      }
    }
  }
}

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid)
{
  const std::size_t rows = old_grid.get_row_size();

  unsigned hardware_threads = std::thread::hardware_concurrency();
  
  if (hardware_threads == 0) //safety
  {
    std::cout<<"No threads detected. Falling back to 1 thread."<<std::endl;
    hardware_threads = 1;
  }

  const unsigned nthreads = static_cast<unsigned>(std::min<std::size_t>(hardware_threads, rows)); //cap threads

  if (rows == 0 || nthreads <= 1) //safety
  {
    process_rows(0, rows, old_grid, new_grid);
    return;
  }

  std::vector<std::thread> thread_pool;
  thread_pool.reserve(nthreads - 1);

  const std::size_t chunk = (rows + nthreads - 1) / nthreads;
  std::size_t start = 0;

  for (unsigned t = 0; t < nthreads && start < rows; t++)
  {
    const std::size_t end = std::min(start + chunk, rows);
    if (t == nthreads - 1)
    {
      process_rows(start, end, old_grid, new_grid);
    }
    else
    {
      thread_pool.emplace_back(process_rows, start, end, std::cref(old_grid), std::ref(new_grid));
    }
    start = end;
  }

  for (auto& worker_thread : thread_pool)
  {
    worker_thread.join();
  }
}