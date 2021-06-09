//include/deal.II-translator/base/cuda_size_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_cuda_size_h
#define dealii_cuda_size_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * 定义启动CUDA内核时的块的大小。这个数字可以根据代码所运行的架构来改变。
   *
   */
  constexpr int block_size = 512;

  /**
   * 定义一个线程所处理的数据块的大小。这个数字可以根据代码运行的架构而改变。
   *
   */
  constexpr int chunk_size = 1;

  /**
   * 定义一个经线中的线程数量。
   *
   */
  constexpr int warp_size = 32;

  /**
   * 定义可以使用 CUDAWrappers::MatrixFree.
   * 求解的最大有限元度，改变这个数字将影响正在使用的恒定内存量。
   *
   */
  constexpr unsigned int mf_max_elem_degree = 10;

  /**
   * 定义有效 CUDAWrappers::MatrixFree 对象的最大数量。
   * 改变这个数字将影响正在使用的恒定内存的数量。
   *
   */
  constexpr unsigned int mf_n_concurrent_objects = 5;
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif


