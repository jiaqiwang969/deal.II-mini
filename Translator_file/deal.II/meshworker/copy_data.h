//include/deal.II-translator/meshworker/copy_data_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_meshworker_copy_data_h
#define dealii_meshworker_copy_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * 帮助复制数据结构。    这个类对于 WorkStream::run() 和
   * MeshWorker::mesh_loop()
   * 函数来说是一个很好的默认投放的CopyData对象。
   * 它是（局部）全矩阵、向量和局部自由度索引向量的数组，其大小由相应的模板参数决定。
   * 特别是，你可以指定以下模板参数
   *
   *
   *
   *
   *
   * -  @tparam  n_matrices: 矩阵阵列的大小
   *
   *
   *
   *
   *
   *
   * -  @tparam  n_vectors: 向量阵列的大小
   *
   *
   *
   *
   *
   *
   * -  @tparam  n_dof_indices: 本地dof指数阵列的大小
   *
   */
  template <int n_matrices    = 1,
            int n_vectors     = n_matrices,
            int n_dof_indices = n_matrices>
  struct CopyData
  {
    /**
     * 用相同的 @p size. 初始化所有东西
     * 这通常是局部自由度的数量。
     *
     */
    explicit CopyData(const unsigned int size);

    /**
     * 对于每个对象，指定它们应该有的尺寸。
     *
     */
    explicit CopyData(
      const ndarray<unsigned int, n_matrices, 2> &   matrix_sizes,
      const std::array<unsigned int, n_vectors> &    vector_sizes,
      const std::array<unsigned int, n_dof_indices> &dof_indices_sizes);

    /**
     * 深度复制构造函数。
     *
     */
    CopyData(const CopyData<n_matrices, n_vectors, n_dof_indices> &other) =
      default;

    /**
     * 一个本地矩阵的数组。
     *
     */
    std::array<FullMatrix<double>, n_matrices> matrices;

    /**
     * 一个本地向量的数组。
     *
     */
    std::array<Vector<double>, n_vectors> vectors;

    /**
     * 一个本地自由度指数的数组。
     *
     */
    std::array<std::vector<types::global_dof_index>, n_dof_indices>
      local_dof_indices;
  };


#ifndef DOXYGEN
  //
  // Template definitions
  //
  template <int n_matrices, int n_vectors, int n_dof_indices>
  CopyData<n_matrices, n_vectors, n_dof_indices>::CopyData(
    const unsigned int size)
  {
    for (auto &m : matrices)
      m.reinit({size, size});
    for (auto &v : vectors)
      v.reinit(size);
    for (auto &d : local_dof_indices)
      d.resize(size);
  }



  template <int n_matrices, int n_vectors, int n_dof_indices>
  CopyData<n_matrices, n_vectors, n_dof_indices>::CopyData(
    const ndarray<unsigned int, n_matrices, 2> &   matrix_sizes,
    const std::array<unsigned int, n_vectors> &    vector_sizes,
    const std::array<unsigned int, n_dof_indices> &dof_indices_sizes)
  {
    for (unsigned int i = 0; i < n_matrices; ++i)
      matrices[i].reinit(matrix_sizes[i++]);

    for (unsigned int i = 0; i < n_vectors; ++i)
      vectors[i].reinit(vector_sizes[i++]);

    for (unsigned int i = 0; i < n_dof_indices; ++i)
      local_dof_indices[i].resize(dof_indices_sizes[i++]);
  }

#endif // DOXYGEN
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


