//include/deal.II-translator/lac/generic_linear_algebra_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_generic_linear_algebra_h
#define dealii_generic_linear_algebra_h

#include <deal.II/base/config.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN


/**
 * 一个命名空间，其中deal.II线性代数类被别名为通用名称。还有类似的命名空间LinearAlgebraPETSc和LinearAlgebraTrilinos，用于别名与PETSc和Trilinos库接口的类。
 *
 *
 */
namespace LinearAlgebraDealII
{
  /**
   * 用于向量类型的类型化定义
   *
   */
  using Vector = Vector<double>;

  /**
   * 用于块向量类型的类型定义
   *
   */
  using BlockVector = BlockVector<double>;

  /**
   * 用于稀疏矩阵类型的类型定义
   *
   */
  using SparseMatrix = SparseMatrix<double>;

  /**
   * 描述由多个块组成的稀疏矩阵的类型定义。
   *
   */
  using BlockSparseMatrix = BlockSparseMatrix<double>;

  /**
   * 用于SSOR预处理程序的类型定义
   *
   */
  using PreconditionSSOR = PreconditionSSOR<SparseMatrix>;
} // namespace LinearAlgebraDealII


DEAL_II_NAMESPACE_CLOSE


#ifdef DEAL_II_WITH_PETSC

#  include <deal.II/lac/petsc_block_sparse_matrix.h>
#  include <deal.II/lac/petsc_precondition.h>
#  include <deal.II/lac/petsc_solver.h>
#  include <deal.II/lac/petsc_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，其中PETSc线性代数类的封装器被别名为通用名称。还有类似的命名空间
 * LinearAlgebraDealII 和 LinearAlgebraTrilinos，用于别名 deal.II
 * 自己的类和与 Trilinos 接口的类。
 *
 *
 */
namespace LinearAlgebraPETSc
{
  /**
   * 用于CG求解器类型的类型化定义。
   *
   */
  using SolverCG = PETScWrappers::SolverCG;

  /**
   * 用于GMRES求解器类型的类型定义。
   *
   */
  using SolverGMRES = PETScWrappers::SolverGMRES;

  /**
   * 一个带有平行PETSc线性代数对象通用名称别名的命名空间。
   *
   */
  namespace MPI
  {
    /**
     * 用于向量类型的类型定义。
     *
     */
    using Vector = PETScWrappers::MPI::Vector;

    /**
     * 用于描述由多个块组成的向量的类型的类型定义。
     *
     */
    using BlockVector = PETScWrappers::MPI::BlockVector;

    /**
     * 用于描述稀疏矩阵类型的类型定义。
     *
     */
    using SparseMatrix = PETScWrappers::MPI::SparseMatrix;

    /**
     * 用于描述由多个块组成的稀疏矩阵的类型的类型定义。
     *
     */
    using BlockSparseMatrix = PETScWrappers::MPI::BlockSparseMatrix;

    /**
     * 用于压缩块疏散模式的类型定义。
     *
     */
    using BlockCompressedSparsityPattern = dealii::BlockDynamicSparsityPattern;

    /**
     * AMG预处理器类型的类型定义。
     *
     */
    using PreconditionAMG = PETScWrappers::PreconditionBoomerAMG;

    /**
     * 用于不完全Cholesky预处理的类型定义。
     *
     */
    using PreconditionIC = PETScWrappers::PreconditionICC;

    /**
     * 不完全LU分解预处理器的类型定义。
     *
     */
    using PreconditionILU = PETScWrappers::PreconditionILU;

    /**
     * 不完全雅各比分解预处理程序的类型定义。
     *
     */
    using PreconditionJacobi = PETScWrappers::PreconditionJacobi;

    /**
     * SSOR预处理程序的类型定义。
     *
     */
    using PreconditionSSOR = PETScWrappers::PreconditionSSOR;

  } // namespace MPI

} // namespace LinearAlgebraPETSc
DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_PETSC

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/lac/trilinos_block_sparse_matrix.h>
#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_solver.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，其中Trilinos线性代数类的包装器被别名为通用名称。还有类似的命名空间
 * LinearAlgebraDealII 和 LinearAlgebraPETSc，用于别名 deal.II
 * 自己的类和与 PETSc 接口的类。
 *
 *
 */
namespace LinearAlgebraTrilinos
{
  /**
   * 用于CG求解器类型的类型化定义。
   *
   */
  using SolverCG = TrilinosWrappers::SolverCG;

  /**
   * 用于GMRES求解器类型的类型定义。
   *
   */
  using SolverGMRES = TrilinosWrappers::SolverGMRES;

  /**
   * 一个带有别名的命名空间，用于并行Trilinos线性代数对象的通用名称。
   *
   */
  namespace MPI
  {
    /**
     * 用于向量类型的类型定义。
     *
     */
    using Vector = TrilinosWrappers::MPI::Vector;

    /**
     * 用于描述由多个块组成的向量的类型的类型定义。
     *
     */
    using BlockVector = TrilinosWrappers::MPI::BlockVector;

    /**
     * 用于描述稀疏矩阵类型的类型定义。
     *
     */
    using SparseMatrix = TrilinosWrappers::SparseMatrix;

    /**
     * 用于描述由多个块组成的稀疏矩阵的类型的类型定义。
     *
     */
    using BlockSparseMatrix = TrilinosWrappers::BlockSparseMatrix;

    /**
     * 用于压缩块稀疏模式的类型的类型定义。
     *
     */
    using BlockCompressedSparsityPattern =
      TrilinosWrappers::BlockSparsityPattern;

    /**
     * AMG预处理器类型的类型定义。
     *
     */
    using PreconditionAMG = TrilinosWrappers::PreconditionAMG;

    /**
     * 用于不完全Cholesky预处理程序的类型定义。
     *
     */
    using PreconditionIC = TrilinosWrappers::PreconditionIC;

    /**
     * 不完全LU分解预处理程序的类型定义。
     *
     */
    using PreconditionILU = TrilinosWrappers::PreconditionILU;

    /**
     * 不完全雅各比分解预处理程序的类型定义。
     *
     */
    using PreconditionJacobi = TrilinosWrappers::PreconditionJacobi;

    /**
     * SSOR预处理程序的类型定义
     *
     */
    using PreconditionSSOR = TrilinosWrappers::PreconditionSSOR;


  } // namespace MPI

} // namespace LinearAlgebraTrilinos

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_TRILINOS



#endif


