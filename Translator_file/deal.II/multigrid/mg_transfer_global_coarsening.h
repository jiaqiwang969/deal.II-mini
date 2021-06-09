//include/deal.II-translator/multigrid/mg_transfer_global_coarsening_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_mg_transfer_global_coarsening_h
#define dealii_mg_transfer_global_coarsening_h

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace internal
{
  class MGTwoLevelTransferImplementation;
}

namespace RepartitioningPolicyTools
{
  template <int dim, int spacedim>
  class Base;
}
#endif



/**
 * 全局粗化的实用函数。
 *
 *
 */
namespace MGTransferGlobalCoarseningTools
{
  /**
   * 常见的多项式粗化序列。
   * @note
   * 这些高达9度的多项式粗化序列在MGTwoLevelTransfer中被预先编译了。也请参见。
   * MGTwoLevelTransfer::fast_polynomial_transfer_supported()
   *
   */
  enum class PolynomialCoarseningSequenceType
  {
    /**
     * 通过整数除法得到半数多项式的度数。例如，对于度数=7，将得到以下序列：。7
     *
     * -> 3
     *
     * -> 1
     *
     */
    bisect,
    /**
     * 将多项式的度数减少一个。例如，对于度数=7的情况，将产生以下序列：7
     *
     * -> 6
     *
     * -> 5
     *
     * -> 4
     *
     * -> 3
     *
     * -> 2
     *
     * -> 1
     *
     */
    decrease_by_one,
    /**
     * 将多项式的度数减少到1。例如，对于度数=7的情况，将产生以下序列：7
     *
     * -> 1
     *
     */
    go_to_one
  };

  /**
   * 对于给定的 @p degree 和多项式粗化序列 @p p_sequence,
   * ，确定下一个粗化程度。
   *
   */
  unsigned int
  create_next_polynomial_coarsening_degree(
    const unsigned int                      degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * 对于给定的 @p max_degree 和多项式粗化序列 @p p_sequence,
   * ，确定多项式度数的完整序列，按升序排序。
   *
   */
  std::vector<unsigned int>
  create_polynomial_coarsening_sequence(
    const unsigned int                      max_degree,
    const PolynomialCoarseningSequenceType &p_sequence);

  /**
   * 对于一个给定的三角形 @p tria,
   * 通过对所提供的三角形重复进行全局粗化来确定几何粗化序列。
   * @note
   * 为方便起见，在返回向量的最后一项中存储了对输入三角形的引用。
   * @note  目前，没有为 parallel::fullydistributed::Triangulation.
   * 实现。
   * @note  返回三角形的类型与输入三角形的类型相同。
   *
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &tria);

  /**
   * 类似于上面的函数，但也取一个 @p policy
   * ，用于在更粗的层次上重新划分三角形。如果 @p
   * preserve_fine_triangulation
   * 被设置，输入的三角形不被改变，否则三角形被粗化。如果
   * @p repartition_fine_triangulation
   * 被设置，最细层的三角形也会被重新划分。如果标志被设置为true/false，则输入的三角图被简单地用作最细的三角图。
   * @note
   * 为方便起见，在返回向量的最后一项中存储了对输入三角形的引用。
   * @note  返回三角形的类型是
   * parallel::fullydistributed::Triangulation.  。
   * @note  目前，只对 parallel::distributed::Triangulation. 实现。
   *
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    Triangulation<dim, spacedim> &                        tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool preserve_fine_triangulation,
    const bool repartition_fine_triangulation);

  /**
   * 类似于上面的函数，但吸收了 @p tria
   * 的恒定版本，因此不允许直接用于粗化，需要在内部创建一个时间性的副本。
   *
   */
  template <int dim, int spacedim>
  std::vector<std::shared_ptr<const Triangulation<dim, spacedim>>>
  create_geometric_coarsening_sequence(
    const Triangulation<dim, spacedim> &                  tria,
    const RepartitioningPolicyTools::Base<dim, spacedim> &policy,
    const bool repartition_fine_triangulation = false);

} // namespace MGTransferGlobalCoarseningTools



/**
 * 用于在两个多网格层次之间进行p-或全局粗化的转移的类。
 *
 *
 */
template <int dim, typename VectorType>
class MGTwoLevelTransfer
{
public:
  /**
   * 进行延长。
   *
   */
  void
  prolongate(VectorType &dst, const VectorType &src) const;

  /**
   * 执行限制。
   *
   */
  void
  restrict_and_add(VectorType &dst, const VectorType &src) const;

  /**
   * 对解向量进行从细级到粗级的插值。这个功能与限制不同，在限制中，加权残差被转移到一个更粗的层次（延长矩阵的转置）。
   *
   */
  void
  interpolate(VectorType &dst, const VectorType &src) const;
};



/**
 * 用于在两个多网格层次之间转移p-或全局粗化的类。对
 * LinearAlgebra::distributed::Vector. 的专门化。
 *
 */
template <int dim, typename Number>
class MGTwoLevelTransfer<dim, LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * 在给定的DoFHandler对象之间设置全局粗化（  @p
   * dof_handler_fine  和  @p dof_handler_coarse).
   * 转移只能在活动层上进行。
   *
   */
  void
  reinit_geometric_transfer(const DoFHandler<dim> &          dof_handler_fine,
                            const DoFHandler<dim> &          dof_handler_coarse,
                            const AffineConstraints<Number> &constraint_fine,
                            const AffineConstraints<Number> &constraint_coarse);

  /**
   * 在给定的DoFHandler对象之间设置多项式粗化 (  @p
   * dof_handler_fine  和  @p dof_handler_coarse).
   * 多项式转移只能在活动层  (`numbers::invalid_unsigned_int`)
   * 或粗格层上执行。
   * @note
   * 函数polynomial_transfer_supported()可以用来检查是否支持给定的多项式粗化策略。
   *
   */
  void
  reinit_polynomial_transfer(
    const DoFHandler<dim> &          dof_handler_fine,
    const DoFHandler<dim> &          dof_handler_coarse,
    const AffineConstraints<Number> &constraint_fine,
    const AffineConstraints<Number> &constraint_coarse,
    const unsigned int mg_level_fine   = numbers::invalid_unsigned_int,
    const unsigned int mg_level_coarse = numbers::invalid_unsigned_int);

  /**
   * 检查是否有 @p fe_degree_fine 和 @p fe_degree_coarse
   * 之间的快速模板化版本的多项式转移。
   * @note 目前，多项式的粗化策略。1) go-to-one, 2) bisect, and 3)
   * decrease-by-one 是预先编译好的9度以下的模板。
   *
   */
  static bool
  fast_polynomial_transfer_supported(const unsigned int fe_degree_fine,
                                     const unsigned int fe_degree_coarse);

  /**
   * 进行延长。
   *
   */
  void
  prolongate(LinearAlgebra::distributed::Vector<Number> &      dst,
             const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * 执行限制。
   *
   */
  void
  restrict_and_add(LinearAlgebra::distributed::Vector<Number> &      dst,
                   const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * 执行从细级到粗级的解向量的插值。这个功能与限制不同，在限制中，加权残差被转移到一个更粗的层次（延长矩阵的转置）。
   *
   */
  void
  interpolate(LinearAlgebra::distributed::Vector<Number> &      dst,
              const LinearAlgebra::distributed::Vector<Number> &src) const;

private:
  /**
   * 一个多栅格转移方案。一个multrigrid转移类可以有不同的转移方案，以实现p-adaptivity（每个多项式度数对有一个转移方案）和实现全局粗化（子单元和父单元之间的转移有一个转移方案，以及，未被细化的单元有一个转移方案）。
   *
   */
  struct MGTransferScheme
  {
    /**
     * 粗化单元的数量。
     *
     */
    unsigned int n_coarse_cells;

    /**
     * 一个粗单元的自由度数。
     *
     */
    unsigned int dofs_per_cell_coarse;

    /**
     * 精细单元的自由度数。
     *
     */
    unsigned int dofs_per_cell_fine;

    /**
     * 粗单元的有限元的多项式程度。
     *
     */
    unsigned int degree_coarse;

    /**
     * 精细单元的有限元的多项式程度。
     *
     */
    unsigned int degree_fine;

    /**
     * 连续元素的权重。
     *
     */
    std::vector<Number> weights;

    /**
     * 非张量积元素的延长矩阵。
     *
     */
    AlignedVector<VectorizedArray<Number>> prolongation_matrix;

    /**
     * 张量-乘积元素的一维延长矩阵。
     *
     */
    AlignedVector<VectorizedArray<Number>> prolongation_matrix_1d;

    /**
     * 非张量产品元素的限制矩阵。
     *
     */
    AlignedVector<VectorizedArray<Number>> restriction_matrix;

    /**
     * 张量-产品元素的一维限制矩阵。
     *
     */
    AlignedVector<VectorizedArray<Number>> restriction_matrix_1d;

    /**
     * 粗大单元的DoF指数，以MPI等级的本地指数表示。
     *
     */
    std::vector<unsigned int> level_dof_indices_coarse;

    /**
     * 精细单元的DoF指数，以MPI等级的本地指数表示。
     *
     */
    std::vector<unsigned int> level_dof_indices_fine;
  };

  /**
   * 传输方案。
   *
   */
  std::vector<MGTransferScheme> schemes;

  /**
   * 标志着精细单元上的有限元是否是连续的。如果是，共享一个顶点/线的DoF的倍数以及约束必须通过权重来考虑。
   *
   */
  bool fine_element_is_continuous;

  /**
   * 中间向量所需的分区器。
   *
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_fine;

  /**
   * 中间向量需要的分区器。
   *
   */
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_coarse;

  /**
   * 收集精细单元的所有自由度所需的内部向量。它只有在精细级自由度指数触及局部活动自由度以外的自由度时才被初始化（我们总是假设可以通过延长/限制函数中的给定向量来访问这些自由度），否则它的大小将被保留为零。
   *
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_fine;

  /**
   * 内部向量，实际延长/限制是在其上进行的。
   *
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_coarse;

  /**
   * 用于执行手动约束_coarse.distribution()的内部向量，这对于可接受的性能是必需的。
   *
   */
  mutable LinearAlgebra::distributed::Vector<Number> vec_coarse_constraints;

  /**
   * 用于在MPI本地索引中手动执行constraint_coarse.distribution()的约束条目索引（出于性能原因）。
   *
   */
  std::vector<unsigned int> constraint_coarse_distribute_indices;

  /**
   * 用于在MPI-本地索引中手动执行constraint_coarse.distribution()的约束条目值（出于性能原因）。
   *
   */
  std::vector<Number> constraint_coarse_distribute_values;

  /**
   * 指向用于手动执行constraint_coarse.distribution()的约束条目的指针。
   *
   */
  std::vector<unsigned int> constraint_coarse_distribute_ptr;

  /**
   * 用于执行手动 constraint_coarse.distribut_local_to_global()
   * 的约束条目索引。
   *
   */
  std::vector<unsigned int> distribute_local_to_global_indices;

  /**
   * 用于执行手动 constraint_coarse.distribut_local_to_global()
   * 的约束条目值。
   *
   */
  std::vector<Number> distribute_local_to_global_values;

  /**
   * 用于执行手动
   * constraint_coarse.distribut_local_to_global()的约束条目的指针。
   *
   */
  std::vector<unsigned int> distribute_local_to_global_ptr;

  /**
   * 组件的数量。
   *
   */
  unsigned int n_components;

  friend class internal::MGTwoLevelTransferImplementation;
};



/**
 * MGTransferBase的实现。与其他多网格传输运算符相比，用户可以在每一层之间提供单独的MGTwoLevelTransfer类型的传输运算符。
 * 该类目前只适用于基于FE_Q和FE_DGQ元素的张量-产品有限元。涉及这些元素之一的多个组件的系统，以及具有不同元素或其他元素的系统，目前没有实现。
 *
 *
 */
template <int dim, typename VectorType>
class MGTransferGlobalCoarsening : public dealii::MGTransferBase<VectorType>
{
public:
  /**
   * 值类型。
   *
   */
  using Number = typename VectorType::value_type;

  /**
   * 构造函数取一个转移运算符的集合（在 @p transfer)
   * 中，最粗大的层次保持为空，如果在PreconditionMG的背景下使用，则有一个可选的函数在函数调用copy_to_mg()内初始化内部层次向量。
   *
   */
  MGTransferGlobalCoarsening(
    const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer,
    const std::function<void(const unsigned int, VectorType &)>
      &initialize_dof_vector = {});

  /**
   * 进行延长。
   *
   */
  void
  prolongate(const unsigned int to_level,
             VectorType &       dst,
             const VectorType & src) const override;

  /**
   * 执行限制。
   *
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType &       dst,
                   const VectorType & src) const override;

  /**
   * 初始化内部向量，并将 @p src
   * 向量复制到最细的多网格层。
   * @note DoFHandler在这里不需要，但被接口要求。
   *
   */
  template <class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<VectorType> &      dst,
             const InVector &                 src) const;

  /**
   * 初始化内部向量，并将最细多网格层上的值复制到 @p dst
   * 向量。
   * @note DoFHandler在这里不需要，但被接口所需要。
   *
   */
  template <class OutVector, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &dof_handler,
               OutVector &                      dst,
               const MGLevelObject<VectorType> &src) const;

  /**
   * 将精细网格场 @p src 插值到 @p dof_handler
   * 中的每个多网格层次，并将结果存储在 @p dst.
   * 中。这个函数与限制不同，在限制中，加权残差被转移到更粗糙的层次（延长矩阵的转置）。
   * 参数 @p dst
   * 必须根据三角化的层数以正确的大小进行初始化。
   * 如果 @p dst
   * 的内向量是空的或者有不正确的局部拥有的大小，它将被调整为每一层的局部相关自由度。
   * @note DoFHandler在这里不需要，但被接口所需要。
   *
   */
  template <class InVector, int spacedim>
  void
  interpolate_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
                    MGLevelObject<VectorType> &      dst,
                    const InVector &                 src) const;

private:
  /**
   * 两级转移运算符的集合。
   *
   */
  const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer;

  /**
   * 用于初始化内部级别向量的%函数。
   *
   */
  const std::function<void(const unsigned int, VectorType &)>
    initialize_dof_vector;
};



#ifndef DOXYGEN

 /* ----------------------- Inline functions --------------------------------- */ 



template <int dim, typename VectorType>
MGTransferGlobalCoarsening<dim, VectorType>::MGTransferGlobalCoarsening(
  const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfer,
  const std::function<void(const unsigned int, VectorType &)>
    &initialize_dof_vector)
  : transfer(transfer)
  , initialize_dof_vector(initialize_dof_vector)
{}



template <int dim, typename VectorType>
void
MGTransferGlobalCoarsening<dim, VectorType>::prolongate(
  const unsigned int to_level,
  VectorType &       dst,
  const VectorType & src) const
{
  this->transfer[to_level].prolongate(dst, src);
}



template <int dim, typename VectorType>
void
MGTransferGlobalCoarsening<dim, VectorType>::restrict_and_add(
  const unsigned int from_level,
  VectorType &       dst,
  const VectorType & src) const
{
  this->transfer[from_level].restrict_and_add(dst, src);
}



template <int dim, typename VectorType>
template <class InVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::copy_to_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  MGLevelObject<VectorType> &      dst,
  const InVector &                 src) const
{
  (void)dof_handler;

  Assert(
    initialize_dof_vector,
    ExcMessage(
      "To be able to use this function, a function to initialize an internal "
      "DoF vector has to be provided in the constructor of "
      "MGTransferGlobalCoarsening."));

  for (unsigned int level = dst.min_level(); level <= dst.max_level(); ++level)
    initialize_dof_vector(level, dst[level]);

  dst[dst.max_level()].copy_locally_owned_data_from(src);
}



template <int dim, typename VectorType>
template <class OutVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::copy_from_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  OutVector &                      dst,
  const MGLevelObject<VectorType> &src) const
{
  (void)dof_handler;

  dst.copy_locally_owned_data_from(src[src.max_level()]);
}



template <int dim, typename VectorType>
template <class InVector, int spacedim>
void
MGTransferGlobalCoarsening<dim, VectorType>::interpolate_to_mg(
  const DoFHandler<dim, spacedim> &dof_handler,
  MGLevelObject<VectorType> &      dst,
  const InVector &                 src) const
{
  (void)dof_handler;

  Assert(
    initialize_dof_vector,
    ExcMessage(
      "To be able to use this function, a function to initialize an internal "
      "DoF vector has to be provided in the constructor of "
      "MGTransferGlobalCoarsening."));

  const unsigned int min_level = transfer.min_level();
  const unsigned int max_level = transfer.max_level();

  AssertDimension(min_level, dst.min_level());
  AssertDimension(max_level, dst.max_level());

  for (unsigned int level = min_level; level <= max_level; ++level)
    initialize_dof_vector(level, dst[level]);

  dst[transfer.max_level()].copy_locally_owned_data_from(src);

  for (unsigned int l = max_level; l > min_level; --l)
    this->transfer[l].interpolate(dst[l - 1], dst[l]);
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif


