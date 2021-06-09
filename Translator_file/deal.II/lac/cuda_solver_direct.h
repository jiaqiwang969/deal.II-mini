//include/deal.II-translator/lac/cuda_solver_direct_0.txt
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

#ifndef dealii_cuda_solver_direct_h
#define dealii_cuda_solver_direct_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE
#  include <deal.II/base/cuda.h>

#  include <deal.II/lac/cuda_sparse_matrix.h>
#  include <deal.II/lac/cuda_vector.h>
#  include <deal.II/lac/solver_control.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * 直接求解器。这些求解器在下面调用cuSOLVER。
   * @note  此模板的实例化提供给<tt>  @<float@></tt>  和 <tt>
   * @<double@></tt>.  。
   * @ingroup CUDAWrappers
   *
   */
  template <typename Number>
  class SolverDirect
  {
  public:
    /**
     * 用于SolverDirect额外设置的结构。
     *
     */
    struct AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的解算器。
       *
       */
      explicit AdditionalData(const std::string &solver_type = "LU_dense");

      /**
       * 设置解算器类型。可能的情况是。        <ul>   <li>  "Cholesky"，在设备上执行Cholesky分解  </li>   <li>  "LU_dense"，将稀疏矩阵转换为密集矩阵并使用LU分解  </li>   <li>  "LU_host"，在主机上使用LU分解  </li>   </ul>  。
       *
       */
      std::string solver_type;
    };

    /**
     * 构造函数。接受求解器控制对象并创建求解器。
     *
     */
    SolverDirect(const Utilities::CUDA::Handle &handle,
                 SolverControl &                cn,
                 const AdditionalData &         data = AdditionalData());

    /**
     * 解构器。
     *
     */
    virtual ~SolverDirect() = default;

    /**
     * 解决线性系统<tt>Ax=b</tt>。
     *
     */
    void
    solve(const SparseMatrix<Number> &                       A,
          LinearAlgebra::CUDAWrappers::Vector<Number> &      x,
          const LinearAlgebra::CUDAWrappers::Vector<Number> &b);

    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;

  private:
    /**
     * 处理
     *
     */
    const Utilities::CUDA::Handle &cuda_handle;

    /**
     * 对控制迭代求解器收敛的对象的引用。事实上，对于这些CUDA封装器来说，cuSOLVER和cuSPARSE本身就是这样做的，但是我们在开始求解过程之前从这个对象中复制数据，之后再将数据复制回这个对象中。
     *
     */
    SolverControl &solver_control;

    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };
} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif

#endif


