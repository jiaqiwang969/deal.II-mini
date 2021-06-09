//include/deal.II-translator/multigrid/mg_transfer_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#ifndef dealii_mg_transfer_h
#define dealii_mg_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  template <typename VectorType>
  struct MatrixSelector
  {
    using Sparsity = ::dealii::SparsityPattern;
    using Matrix   = ::dealii::SparseMatrix<typename VectorType::value_type>;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &                   matrix,
           Sparsity &                 sparsity,
           int                        level,
           const SparsityPatternType &sp,
           const DoFHandler<dim, spacedim> &)
    {
      sparsity.copy_from(sp);
      (void)level;
      matrix.reinit(sparsity);
    }
  };

#ifdef DEAL_II_WITH_TRILINOS
  template <typename Number>
  struct MatrixSelector<LinearAlgebra::distributed::Vector<Number>>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                              level,
           const SparsityPatternType &      sp,
           const DoFHandler<dim, spacedim> &dh)
    {
      const MPI_Comm communicator = dh.get_communicator();

      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };

  template <>
  struct MatrixSelector<dealii::TrilinosWrappers::MPI::Vector>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                              level,
           const SparsityPatternType &      sp,
           const DoFHandler<dim, spacedim> &dh)
    {
      const MPI_Comm communicator = dh.get_communicator();

      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };

#  ifdef DEAL_II_WITH_MPI
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number>
  struct MatrixSelector<dealii::LinearAlgebra::TpetraWrappers::Vector<Number>>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                              level,
           const SparsityPatternType &      sp,
           const DoFHandler<dim, spacedim> &dh)
    {
      const MPI_Comm communicator = dh.get_communicator();

      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };
#    endif

  template <>
  struct MatrixSelector<dealii::LinearAlgebra::EpetraWrappers::Vector>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                              level,
           const SparsityPatternType &      sp,
           const DoFHandler<dim, spacedim> &dh)
    {
      const MPI_Comm communicator = dh.get_communicator();

      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };
#  endif

#else
  // ! DEAL_II_WITH_TRILINOS
  template <typename Number>
  struct MatrixSelector<LinearAlgebra::distributed::Vector<Number>>
  {
    using Sparsity = ::dealii::SparsityPattern;
    using Matrix   = ::dealii::SparseMatrix<Number>;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &,
           Sparsity &,
           int,
           const SparsityPatternType &,
           const DoFHandler<dim, spacedim> &)
    {
      AssertThrow(
        false,
        ExcNotImplemented(
          "ERROR: MGTransferPrebuilt with LinearAlgebra::distributed::Vector currently "
          "needs deal.II to be configured with Trilinos."));
    }
  };

#endif

#ifdef DEAL_II_WITH_PETSC
  template <>
  struct MatrixSelector<dealii::PETScWrappers::MPI::Vector>
  {
    using Sparsity = ::dealii::DynamicSparsityPattern;
    using Matrix   = ::dealii::PETScWrappers::MPI::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = true;

    template <typename SparsityPatternType, int dim, int spacedim>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                              level,
           const SparsityPatternType &      sp,
           const DoFHandler<dim, spacedim> &dh)
    {
      const MPI_Comm communicator = dh.get_communicator();

      // Reinit PETSc matrix
      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator);
    }
  };
#endif
} // namespace internal

/* MGTransferBase在mg_base.h中定义。

* 
*
*/

 /*!@addtogroup mg */ 
 /*@{*/ 



/**
 * 实现全局向量和多网格层次之间的转移，用于派生类MGTransferPrebuilt和其他类。
 *
 *
 */
template <typename VectorType>
class MGLevelGlobalTransfer : public MGTransferBase<VectorType>
{
public:
  /**
   * 将对象重置为默认构造函数后的状态。
   *
   */
  void
  clear();

  /**
   * 从全局网格上的矢量转移到为活动自由度分别定义在各层的矢量。特别是，对于一个全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，所有其他的层次对象都没有被触动。
   *
   */
  template <int dim, class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<VectorType> &      dst,
             const InVector &                 src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将MGVector的活动部分的数据复制到<tt>Vector<number></tt>的相应位置。为了保持结果的一致性，受限自由度被设置为零。
   *
   */
  template <int dim, class OutVector, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &dof_handler,
               OutVector &                      dst,
               const MGLevelObject<VectorType> &src) const;

  /**
   * 将一个多级向量添加到一个法向量。
   * 与前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, class OutVector, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &dof_handler,
                   OutVector &                      dst,
                   const MGLevelObject<VectorType> &src) const;

  /**
   * 如果这个对象对BlockVector对象进行操作，我们需要描述各个矢量分量如何被映射到矢量的块上。
   * 例如，对于斯托克斯系统，我们有dim+1的速度和压力的矢量分量，但我们可能想使用只有两个块的块状矢量，将所有的速度放在一个块中，压力变量放在另一个块中。
   * 默认情况下，如果不调用这个函数，块向量的块数与有限元的向量分量一样多。然而，这可以通过调用这个函数来改变，该函数用一个数组来描述如何将矢量分量分组为块。该参数的含义与给
   * DoFTools::count_dofs_per_component 函数的参数相同。
   *
   */
  void
  set_component_to_block_map(const std::vector<unsigned int> &map);

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 为调试目的打印拷贝索引字段。
   *
   */
  void
  print_indices(std::ostream &os) const;

protected:
  /**
   * @p fill  copy_indices*的内部函数。由派生类调用。
   *
   */
  template <int dim, int spacedim>
  void
  fill_and_communicate_copy_indices(
    const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 多级向量的大小。
   *
   */
  std::vector<types::global_dof_index> sizes;

  /**
   * copy_to_mg()和copy_from_mg()函数的映射。这里只有本地拥有的索引对
   * 数据的组织方式如下：每层一个向量。这些向量的每个元素首先包含全局索引，然后是级别索引。
   *
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices;

  /**
   * 用于copy_to_mg()函数的额外自由度。这些是全局自由度为本地所有，而层次自由度则不是。
   * 数据的组织与 @p copy_indices_mine. 相同。
   *
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices_global_mine;

  /**
   * 用于copy_from_mg()函数的额外自由度。这些是水平自由度是本地拥有的，而全局自由度不是。
   * 数据的组织与 @p copy_indices_mine. 的情况一样。
   *
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices_level_mine;

  /**
   * 这个变量存储了从全局到层次向量的复制操作是否实际上是对最细层次的普通复制。这意味着网格没有自适应细化，最细的多网格层次上的编号与全局的情况相同。
   *
   */
  bool perform_plain_copy;

  /**
   * 存储给set_component_to_block_map()函数的向量。
   *
   */
  std::vector<unsigned int> component_to_block_map;

  /**
   * 水平系统的mg_constrained_dofs。
   *
   */
  SmartPointer<const MGConstrainedDoFs, MGLevelGlobalTransfer<VectorType>>
    mg_constrained_dofs;

private:
  /**
   * 调用这个函数是为了确保build()已经被调用。
   *
   */
  template <int dim, int spacedim>
  void
  assert_built(const DoFHandler<dim, spacedim> &dof_handler) const;
};



/**
 * 实现全局向量和多网格层次之间的转移，用于派生类MGTransferPrebuilt和其他类。这个类是针对
 * LinearAlgebra::distributed::Vector
 * 的情况的特殊化，与PETScWrappers和TrilinosWrappers命名空间中的%并行向量相比，需要一些不同的调用例程。
 *
 *
 */
template <typename Number>
class MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>
  : public MGTransferBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * 将对象重置为默认构造函数后的状态。
   *
   */
  void
  clear();

  /**
   * 从全局网格上的矢量转移到为活动自由度分别定义在每一层的矢量。特别是，对于一个全局细化的网格，只有
   * @p dst 中最细的层次被填充为 @p src.
   * 的普通拷贝，所有其他的层次对象都没有被触动。
   *
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
             const LinearAlgebra::distributed::Vector<Number2> &src) const;

  /**
   * 从多级向量转移到普通向量。
   * 将MGVector的活动部分的数据复制到<tt>Vector<number></tt>的相应位置。为了保持结果的一致性，受限自由度被设置为零。
   *
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_from_mg(
    const DoFHandler<dim, spacedim> &            dof_handler,
    LinearAlgebra::distributed::Vector<Number2> &dst,
    const MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &src) const;

  /**
   * 将一个多级向量添加到一个法向量。
   * 与前面的函数一样工作，但可能不适合连续元素。
   *
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_from_mg_add(
    const DoFHandler<dim, spacedim> &            dof_handler,
    LinearAlgebra::distributed::Vector<Number2> &dst,
    const MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &src) const;

  /**
   * 如果这个对象对BlockVector对象进行操作，我们需要描述各个矢量分量如何被映射到矢量的块上。
   * 例如，对于斯托克斯系统，我们有dim+1的速度和压力的矢量分量，但我们可能想使用只有两个块的块状矢量，将所有速度放在一个块中，压力变量放在另一个块中。
   * 默认情况下，如果不调用这个函数，块向量的块数与有限元的向量分量一样多。然而，这可以通过调用这个函数来改变，该函数用一个数组来描述如何将矢量分量分组为块。该参数的含义与给
   * DoFTools::count_dofs_per_component 函数的参数相同。
   *
   */
  void
  set_component_to_block_map(const std::vector<unsigned int> &map);

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 为调试目的打印拷贝索引字段。
   *
   */
  void
  print_indices(std::ostream &os) const;

protected:
  /**
   * 内部函数，根据标志 @p solution_transfer.
   * 执行残差或解决方案的转移。
   *
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
             const LinearAlgebra::distributed::Vector<Number2> &        src,
             const bool solution_transfer) const;

  /**
   * @p fill  copy_indices*的内部函数。由派生类调用。
   *
   */
  template <int dim, int spacedim>
  void
  fill_and_communicate_copy_indices(
    const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 多级向量的大小。
   *
   */
  std::vector<types::global_dof_index> sizes;

  /**
   * copy_to_mg()和copy_from_mg()函数的映射。这里只存储本地拥有的索引对。
   * 数据的组织方式如下：每层有一个表。这个表有两行。第一行包含全局索引，第二行是级别索引。
   *
   */
  std::vector<Table<2, unsigned int>> copy_indices;

  /**
   * 与上述相同，但用于传输解决方案向量。
   *
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices;

  /**
   * 用于copy_to_mg()函数的额外自由度。这些是全局自由度是本地拥有的，而水平自由度则不是。
   * 数据的组织与 @p copy_indices. 相同。
   *
   */
  std::vector<Table<2, unsigned int>> copy_indices_global_mine;

  /**
   * 与上述相同，但用于转移解向量。
   *
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices_global_mine;

  /**
   * 用于copy_from_mg()函数的额外自由度。这些是水平自由度是本地拥有的，而全局自由度不是。
   * 数据的组织与 @p copy_indices. 一样。
   *
   */
  std::vector<Table<2, unsigned int>> copy_indices_level_mine;

  /**
   * 与上述相同，但用于转移解向量。
   *
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices_level_mine;

  /**
   * 这个变量存储了从全局到水平向量的复制操作是否实际上是对最细水平的普通复制。这意味着网格没有自适应细化，最细的多网格层次上的编号与全局情况下的编号相同。
   *
   */
  bool perform_plain_copy;

  /**
   * 这个变量存储了从全局到层次向量的复制操作，除了在最细层次内对自由度进行重新编号外，实际上是对最细层次的普通复制。这意味着该网格没有自适应细化。
   *
   */
  bool perform_renumbered_plain_copy;

  /**
   * 存储给set_component_to_block_map()函数的向量。
   *
   */
  std::vector<unsigned int> component_to_block_map;

  /**
   * 水平系统的mg_constrained_dofs。
   *
   */
  SmartPointer<
    const MGConstrainedDoFs,
    MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>>
    mg_constrained_dofs;

  /**
   * 在函数copy_to_mg中，我们需要访问全局向量的重影项，以便插入到关卡向量中。这个向量是由这些条目填充的。
   *
   */
  mutable LinearAlgebra::distributed::Vector<Number> ghosted_global_vector;

  /**
   * 和上面一样，但在处理解向量时使用。
   *
   */
  mutable LinearAlgebra::distributed::Vector<Number>
    solution_ghosted_global_vector;

  /**
   * 在函数copy_from_mg中，我们访问所有带有某些鬼魂条目的水平向量，以便将结果插入全局向量中。
   *
   */
  mutable MGLevelObject<LinearAlgebra::distributed::Vector<Number>>
    ghosted_level_vector;

  /**
   * 与上述相同，但在处理解向量时使用。
   *
   */
  mutable MGLevelObject<LinearAlgebra::distributed::Vector<Number>>
    solution_ghosted_level_vector;

private:
  /**
   * 调用这个函数是为了确保build()已经被调用。
   *
   */
  template <int dim, int spacedim>
  void
  assert_built(const DoFHandler<dim, spacedim> &dof_handler) const;
};



/**
 * MGTransferBase接口的实现，该接口的转移操作是在构建该类的对象为矩阵时预先构建的。这是快速的方法，因为它只需要通过循环所有单元格并将结果存储在每一层的矩阵中来构建一次操作，但需要额外的内存。
 * 参见MGTransferBase，了解哪一个转移类最适合你的需要。
 *
 *
 */
template <typename VectorType>
class MGTransferPrebuilt : public MGLevelGlobalTransfer<VectorType>
{
public:
  /**
   * 没有约束矩阵的构造函数。只在不连续的有限元或没有局部细化的情况下使用这个构造函数。
   *
   */
  MGTransferPrebuilt() = default;

  /**
   * 带约束的构造器。相当于默认的构造函数，后面加上initialize_constraints()。
   *
   */
  MGTransferPrebuilt(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 解构器。
   *
   */
  virtual ~MGTransferPrebuilt() override = default;

  /**
   * 初始化将在build()中使用的约束。
   *
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * 将对象重置为默认构造函数之后的状态。
   *
   */
  void
  clear();

  /**
   * 实际构建转移操作所需的信息。在使用prolongate()或restrict_and_add()之前需要调用。
   *
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 使用底层有限元的嵌入矩阵，将一个向量从<tt>to_level-1</tt>延长到<tt>to_level</tt>。之前<tt>dst</tt>的内容被覆盖。
   * @arg
   * src是一个向量，其元素数与所涉及的较粗层次上的自由度相同。
   * @arg  dst有多少个元素，就有多少个更细层次的自由度。
   *
   */
  virtual void
  prolongate(const unsigned int to_level,
             VectorType &       dst,
             const VectorType & src) const override;

  /**
   * 使用 @p prolongate
   * 方法的转置操作，将一个矢量从<tt>from_level</tt>级限制到<tt>from_level-1</tt>级。如果<tt>from_level</tt>层的单元所覆盖的区域小于<tt>from_level-1</tt>层的区域（局部细化），那么<tt>dst</tt>中的一些自由度是有效的，将不会被改变。对于其他自由度，限制的结果将被添加。
   * @arg
   * src是一个向量，其元素数量与所涉及的更细层次上的自由度相同。
   * @arg  dst的元素数与较粗层次上的自由度相同。
   *
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType &       dst,
                   const VectorType & src) const override;

  /**
   * 有限元不提供延长矩阵。
   *
   */
  DeclException0(ExcNoProlongation);

  /**
   * 在使用这个对象之前，你必须调用build()。
   *
   */
  DeclException0(ExcMatricesNotBuilt);

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 为调试目的打印所有的矩阵。
   *
   */
  void
  print_matrices(std::ostream &os) const;

private:
  /**
   * 转移矩阵的稀疏性模式。
   *
   */
  std::vector<
    std::shared_ptr<typename internal::MatrixSelector<VectorType>::Sparsity>>
    prolongation_sparsities;

  /**
   * 实际的延长矩阵。列指数属于母单元的道夫指数，即粗级别。而行指数属于子单元，即细级别。
   *
   */
  std::vector<
    std::shared_ptr<typename internal::MatrixSelector<VectorType>::Matrix>>
    prolongation_matrices;

  /**
   * 细化边上的自由度，不包括边界上的自由度。
   *
   */
  std::vector<std::vector<bool>> interface_dofs;
};


 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


