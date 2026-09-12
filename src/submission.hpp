#pragma once

#include <cmath>
#include <cstddef>
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
    //data_ = std::vector<double> (rows*cols, 0.0);
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

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid& old_grid, Grid& new_grid)
{
  const std::size_t rows = old_grid.get_row_size();
  const std::size_t cols = old_grid.get_col_size();

  //copy boundary rows
  // for(std::size_t j = 0; j < cols; j++)
  // {
  //   new_grid(0,j) = old_grid(0,j);
  //   new_grid(rows - 1, j) = old_grid(rows - 1, j);
  // }

  // //copy boundary columns
  // for(std::size_t i = 0; i < rows; i++)
  // {
  //   new_grid(i,0) = old_grid(i,0);
  //   new_grid(i,cols-1) = old_grid(i,cols-1);
  // }

  //apply formula to inner cells
  for(std::size_t i = 0; i < rows; i++)
  {
    if(i == 0 || i == rows-1)
    {
      for(std::size_t j = 0; j < cols; j++)
      {
        new_grid(i,j) = old_grid(i,j);
      }
    }
    else 
    {
      for(std::size_t j = 0; j < cols; j++)
      {
        new_grid(i,0) = old_grid(i,0);

        new_grid(i,j) = 0.5 * old_grid(i,j) +
                        0.125 * (old_grid(i-1,j) + old_grid(i+1,j) +
                                 old_grid(i,j-1) + old_grid(i,j+1));
        
        new_grid(i,cols-1) = old_grid(i,cols-1);
      }
    }
  }
}
