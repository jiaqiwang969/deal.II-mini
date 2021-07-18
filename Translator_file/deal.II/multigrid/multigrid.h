//include/deal.II-translator/multigrid/multigrid_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_multigrid_h
#define dealii_multigrid_h


#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

#ifdef signals
#  error \
    "The name 'signals' is already defined. You are most likely using the QT library \
and using the 'signals' keyword. You can either #include the Qt headers (or any conflicting headers) \
*after* the deal.II headers or you can define the 'QT_NO_KEYWORDS' macro and use the 'Q_SIGNALS' macro."
#endif

 /*!@addtogroup mg */ 
 /*@{*/ 

namespace mg
{
  /**
   * 一个包含 boost::signal
   * 对象的结构，用于多网格求解器的可选处理。
   * 每个信号都被调用两次，一次是在执行动作之前，一次是在执行动作之后。这两个函数调用的不同之处在于布尔参数
   * @p before, ，第一次为真，第二次为假。
   *
   */
  struct Signals
  {
    /**
     * 这个信号在调用 @p before 之前（ @p before
     * 为真）和之后（ @p before 为假）被触发，
     * MGTransfer::copy_to_mg 将给定的向量转移到一个多级向量。
     *
     */
    boost::signals2::signal<void(const bool before)> transfer_to_mg;

    /**
     * 这个信号在调用 @p before 为真之前和 @p before
     * 为假之后被触发， MGTransfer::copy_from_mg
     * 将给定的多级向量转为正常向量。
     *
     */
    boost::signals2::signal<void(const bool before)> transfer_to_global;

    /**
     * 这个信号在调用 @p before 之前（ @p before
     * 为真）和之后（ @p level.
     * 为假）被触发，粗解将用`defect[leve]`完成，并以`solution[level]`返回，用户可以使用这个信号检查。
     *
     */
    boost::signals2::signal<void(const bool before, const unsigned int level)>
      coarse_solve;

    /**
     * 这个信号在调用 @p before 之前（ @p before
     * 为真）和之后（ @p before
     * 为假）被触发，该信号将一个矢量从 @p level
     * 限制到下一个更粗的矢量（ @p level  ）。
     *
     * - 1).         向量``defect[level-1]`将在这两个触发器之间被更新，用户可以使用这个信号来检查。
     *
     */
    boost::signals2::signal<void(const bool before, const unsigned int level)>
      restriction;

    /**
     * 这个信号在调用 @p before 之前（ @p before
     * 为真）和之后（ @p before 为假）被触发，
     * MGTransfer::prolongate()
     * 将一个向量从下一个较粗的向量延长到 @p level （ @p
     * level  ）。
     *
     * - 1).
     *
     */
    boost::signals2::signal<void(const bool before, const unsigned int level)>
      prolongation;

    /**
     * 这个信号在（ @p before 为真）和（ @p before 为假）通过
     * MGPreSmoother::apply() 对 @p level.
     * 的预平滑步骤的调用之前和之后被触发。平滑的结果将被存储在`solution[level]`中，用户可以通过这个信号来检查。
     *
     */
    boost::signals2::signal<void(const bool before, const unsigned int level)>
      pre_smoother_step;

    /**
     * 这个信号在 @p before 为真之前和 @p before
     * 为假之后被触发，通过 MGPostSmoother::apply() 对 @p level.
     * 的平滑后步骤进行调用。
     *
     */
    boost::signals2::signal<void(const bool before, const unsigned int level)>
      post_smoother_step;
  };
} // namespace mg

/**
 * 多网格方法的实现。该实现支持连续和非连续元素，并遵循 @ref mg_paper "Janssen和Kanschat的多栅格论文 "
 * 中描述的程序。
 * cycle()是在最细级别上启动多网格循环的函数。根据构造函数选择的循环类型（见枚举Cycle），该函数会触发level_v_step()或level_step()中的一个循环，其中后者可以做不同类型的循环。
 * 使用这个类，预计右手边已经从一个生活在本地最好水平上的向量转换为一个多级向量。这是一个不简单的操作，通常由PreconditionMG类自动启动，由MGTransferBase派生的类执行。
 *
 *
 * @note
 * 这个类的接口仍然非常笨拙。特别是，在使用它之前，你必须要设置相当多的辅助性对象。不幸的是，似乎只有以一种不可接受的方式来限制这个类的灵活性，才能避免这种情况。
 *
 *
 */
template <typename VectorType>
class Multigrid : public Subscriptor
{
public:
  /**
   * 实现的循环类型列表。
   *
   */
  enum Cycle
  {
    /// The V-cycle
    v_cycle,
    /// The W-cycle
    w_cycle,
    /// The F-cycle
    f_cycle
  };

  using vector_type       = VectorType;
  using const_vector_type = const VectorType;

  /**
   * 构造函数。<tt>transfer</tt>是一个执行延长和限制的对象。对于[minlevel,
   * maxlevel]中的级别，矩阵必须包含有效的矩阵。默认情况下，maxlevel被设置为最大的有效级别。
   * 这个函数已经初始化了向量，这些向量将在以后的计算过程中使用。因此，你应该尽可能晚地创建这种类型的对象。
   *
   */
  Multigrid(const MGMatrixBase<VectorType> &    matrix,
            const MGCoarseGridBase<VectorType> &coarse,
            const MGTransferBase<VectorType> &  transfer,
            const MGSmootherBase<VectorType> &  pre_smooth,
            const MGSmootherBase<VectorType> &  post_smooth,
            const unsigned int                  minlevel = 0,
            const unsigned int maxlevel = numbers::invalid_unsigned_int,
            Cycle              cycle    = v_cycle);

  /**
   * 根据#minlevel和#maxlevel来重新启动这个类。
   *
   */
  void
  reinit(const unsigned int minlevel, const unsigned int maxlevel);

  /**
   * 执行一个多栅格循环。循环的类型由构造函数参数cycle来选择。可用的类型见枚举循环。
   *
   */
  void
  cycle();

  /**
   * 执行V-cycle算法的一个步骤。
   * 这个函数假定，多级向量#defect是由外部缺陷修正方案的残余物填充的。这通常是由PreconditionMG)处理的。在vcycle()之后，结果是在多级向量#solution中。如果你想自己使用这些向量，请参见MGTools命名空间的<tt>copy_*_mg</tt>。
   * 这个函数的实际工作是在level_v_step()中完成的。
   *
   */
  void
  vcycle();

  /**
   * 设置额外的矩阵来修正细化边缘的残差计算。由于我们只对网格细化部分的内部进行平滑处理，所以缺少对细化边缘的耦合。这种耦合是由这两个矩阵提供的。
   * @note
   * 虽然使用了<tt>edge_out.vmult</tt>，对于第二个参数，我们使用<tt>edge_in.Tvmult</tt>。因此，<tt>edge_in</tt>应该以转置的形式进行组合。这为<tt>edge_in</tt>节省了第二个稀疏模式。特别是，对于对称运算符来说，两个参数都可以指向同一个矩阵，这样可以节省其中一个的装配。
   *
   */
  void
  set_edge_matrices(const MGMatrixBase<VectorType> &edge_out,
                    const MGMatrixBase<VectorType> &edge_in);

  /**
   * 设置额外的矩阵来纠正细化边缘的残差计算。这些矩阵源于不连续的Galerkin方法（见FE_DGQ等），它们对应于两层之间细化边缘的边缘通量。
   * @note
   * 虽然使用了<tt>edge_down.vmult</tt>，对于第二个参数，我们使用<tt>edge_up.Tvmult</tt>。因此，<tt>edge_up</tt>应该以转置的形式进行组合。这为<tt>edge_up</tt>节省了第二个稀疏模式。特别是，对于对称运算符来说，两个参数都可以指向同一个矩阵，这样可以节省其中一个的装配。
   *
   */
  void
  set_edge_flux_matrices(const MGMatrixBase<VectorType> &edge_down,
                         const MGMatrixBase<VectorType> &edge_up);

  /**
   * 返回多网格的最佳水平。
   *
   */
  unsigned int
  get_maxlevel() const;

  /**
   * 返回多网格的最粗层次。
   *
   */
  unsigned int
  get_minlevel() const;

  /**
   * 设置执行多层次方法的最高层次。默认情况下，这是三角法的最细层次。接受的值不小于当前的#minlevel。
   *
   */
  void
  set_maxlevel(const unsigned int);

  /**
   * 设置执行多层次方法的最粗层次。默认情况下，该值为零。接受的是不大于当前#maxlevel的非负值。
   * 如果<tt>relative</tt>是<tt>true</tt>，那么这个函数决定了使用的级别数，也就是说，它将#minlevel设置为#maxlevel-<tt>level</tt>。
   * @note
   * 最粗层次上的网格必须覆盖整个领域。在#minlevel上不能有悬空的节点。
   * @note
   * 如果#minlevel被设置为非零值，不要忘记调整你的粗略网格解算器
   *
   */
  void
  set_minlevel(const unsigned int level, bool relative = false);

  /**
   * 机会#cycle_type在cycle()中使用。
   *
   */
  void set_cycle(Cycle);

  /**
   * 连接一个函数到 mg::Signals::coarse_solve. 。
   *
   */
  boost::signals2::connection
  connect_coarse_solve(
    const std::function<void(const bool, const unsigned int)> &slot);

  /**
   * 将一个函数连接到 mg::Signals::restriction. 。
   *
   */
  boost::signals2::connection
  connect_restriction(
    const std::function<void(const bool, const unsigned int)> &slot);

  /**
   * 将一个函数连接到 mg::Signals::prolongation. 。
   *
   */
  boost::signals2::connection
  connect_prolongation(
    const std::function<void(const bool, const unsigned int)> &slot);

  /**
   * 将一个函数连接到 mg::Signals::pre_smoother_step. 。
   *
   */
  boost::signals2::connection
  connect_pre_smoother_step(
    const std::function<void(const bool, const unsigned int)> &slot);

  /**
   * 将一个函数连接到 mg::Signals::post_smoother_step. 。
   *
   */
  boost::signals2::connection
  connect_post_smoother_step(
    const std::function<void(const bool, const unsigned int)> &slot);

private:
  /**
   * 多重网格算法使用的各种动作的信号。
   *
   */
  mg::Signals signals;

  /**
   * V型循环多网格方法。<tt>level</tt>是函数开始的层次。它通常会从外部调用最高层，但随后会递归调用自己的<tt>level-1</tt>，除非我们在#minlevel上，粗网格求解器完全解决了这个问题。
   *
   */
  void
  level_v_step(const unsigned int level);

  /**
   * 实际的W-循环或F-循环多网格方法。<tt>level</tt>是函数开始的级别。它通常会从外部调用最高层，但随后会递归调用自己的<tt>level-1</tt>，除非我们在#minlevel上，粗网格求解器完全解决了问题。
   *
   */
  void
  level_step(const unsigned int level, Cycle cycle);

  /**
   * 由cycle()方法执行的循环类型。
   *
   */
  Cycle cycle_type;

  /**
   * 用于粗略网格求解的水平。
   *
   */
  unsigned int minlevel;

  /**
   * 最高级别的单元格。
   *
   */
  unsigned int maxlevel;

public:
  /**
   * 循环的输入向量。包含外方法投射到多层次向量的缺陷。
   *
   */
  MGLevelObject<VectorType> defect;

  /**
   * 多网格步骤后的解决方案更新。
   *
   */
  MGLevelObject<VectorType> solution;

private:
  /**
   * 辅助向量。
   *
   */
  MGLevelObject<VectorType> t;

  /**
   * 用于W-和F-循环的辅助向量。在V-循环中未被初始化。
   *
   */
  MGLevelObject<VectorType> defect2;


  /**
   * 每个层次的矩阵。
   *
   */
  SmartPointer<const MGMatrixBase<VectorType>, Multigrid<VectorType>> matrix;

  /**
   * 每一级的矩阵。
   *
   */
  SmartPointer<const MGCoarseGridBase<VectorType>, Multigrid<VectorType>>
    coarse;

  /**
   * 用于网格转移的对象。
   *
   */
  SmartPointer<const MGTransferBase<VectorType>, Multigrid<VectorType>>
    transfer;

  /**
   * 预平滑对象。
   *
   */
  SmartPointer<const MGSmootherBase<VectorType>, Multigrid<VectorType>>
    pre_smooth;

  /**
   * 后平滑对象。
   *
   */
  SmartPointer<const MGSmootherBase<VectorType>, Multigrid<VectorType>>
    post_smooth;

  /**
   * 从细化部分的内部到细化边缘的边缘矩阵。
   * @note 只有<tt>vmult</tt>被用于这些矩阵。
   *
   */
  SmartPointer<const MGMatrixBase<VectorType>> edge_out;

  /**
   * 从细化边缘到细化部分内部的转置边缘矩阵。
   * @note 只有<tt>Tvmult</tt>被用于这些矩阵。
   *
   */
  SmartPointer<const MGMatrixBase<VectorType>> edge_in;

  /**
   * 从细到粗的边缘矩阵。
   * @note  只有<tt>vmult</tt>用于这些矩阵。
   *
   */
  SmartPointer<const MGMatrixBase<VectorType>, Multigrid<VectorType>> edge_down;

  /**
   * 从粗到细的转置边缘矩阵。
   * @note  只有<tt>Tvmult</tt>用于这些矩阵。
   *
   */
  SmartPointer<const MGMatrixBase<VectorType>, Multigrid<VectorType>> edge_up;

  template <int dim, class OtherVectorType, class TRANSFER>
  friend class PreconditionMG;
};


/**
 * 多级预处理器。在这里，我们收集了多级预处理所需的所有信息，并为LAC迭代方法提供标准接口。
 * 此外，它需要函数<tt>void copy_to_mg(const
 * VectorType&)</tt>来存储多级方法右侧的 @p src ，<tt>void
 * copy_from_mg(VectorType&)</tt>来存储 @p dst. 中的v循环结果。
 * 如果VectorType实际上是一个块向量，并且TRANSFER对象支持为每个块使用单独的DoFHandler，这个类也允许被初始化为每个块的单独DoFHandler。
 *
 *
 */
template <int dim, typename VectorType, class TRANSFER>
class PreconditionMG : public Subscriptor
{
public:
  /**
   * 构造函数。参数是多网格对象、前平滑器、后平滑器和粗网格求解器。
   *
   */
  PreconditionMG(const DoFHandler<dim> &dof_handler,
                 Multigrid<VectorType> &mg,
                 const TRANSFER &       transfer);

  /**
   * 在块向量的每个分量都使用其自己的DoFHandler的情况下，与上述相同。
   *
   */
  PreconditionMG(const std::vector<const DoFHandler<dim> *> &dof_handler,
                 Multigrid<VectorType> &                     mg,
                 const TRANSFER &                            transfer);

  /**
   * 其他类所需的假函数。
   *
   */
  bool
  empty() const;

  /**
   * 预设条件运算符。调用传递给构造器的 @p MG 对象的 @p
   * vcycle 函数。    这是LAC迭代求解器使用的运算器。
   *
   */
  template <class OtherVectorType>
  void
  vmult(OtherVectorType &dst, const OtherVectorType &src) const;

  /**
   * 预处理运算符。调用传递给构造函数的 @p MG 对象的 @p
   * vcycle 函数。
   *
   */
  template <class OtherVectorType>
  void
  vmult_add(OtherVectorType &dst, const OtherVectorType &src) const;

  /**
   * 转置的预设条件运算符。    未实现，但可能需要定义。
   *
   */
  template <class OtherVectorType>
  void
  Tvmult(OtherVectorType &dst, const OtherVectorType &src) const;

  /**
   * 转置的预设条件运算符。    未实现，但可能需要定义。
   *
   */
  template <class OtherVectorType>
  void
  Tvmult_add(OtherVectorType &dst, const OtherVectorType &src) const;

  /**
   * 返回这个预处理程序的范围空间的划分，即从矩阵-向量乘积得到的向量的划分。默认情况下，会返回第一个DoFHandler对象的相应信息。
   *
   */
  IndexSet
  locally_owned_range_indices(const unsigned int block = 0) const;

  /**
   * 返回该预处理程序的域空间的划分，即该矩阵必须与之相乘的向量的划分。
   * 默认情况下，会返回第一个DoFHandler对象的相应信息。
   *
   */
  IndexSet
  locally_owned_domain_indices(const unsigned int block = 0) const;

  /**
   * 返回与该预处理程序一起使用的MPI通信器对象。
   *
   */
  MPI_Comm
  get_mpi_communicator() const;

  /**
   * 连接一个函数到 mg::Signals::transfer_to_mg. 。
   *
   */
  boost::signals2::connection
  connect_transfer_to_mg(const std::function<void(bool)> &slot);

  /**
   * 将一个函数连接到 mg::Signals::transfer_to_global. 。
   *
   */
  boost::signals2::connection
  connect_transfer_to_global(const std::function<void(bool)> &slot);

private:
  /**
   * 关联 @p DoFHandler. 。
   *
   */
  std::vector<SmartPointer<const DoFHandler<dim>,
                           PreconditionMG<dim, VectorType, TRANSFER>>>
    dof_handler_vector;

  /**
   * 为没有SmartPointer包装的DoFHandler对象的指针进行存储。
   *
   */
  std::vector<const DoFHandler<dim> *> dof_handler_vector_raw;

  /**
   * 多网格对象。
   *
   */
  SmartPointer<Multigrid<VectorType>, PreconditionMG<dim, VectorType, TRANSFER>>
    multigrid;

  /**
   * 用于网格传输的对象。
   *
   */
  SmartPointer<const TRANSFER, PreconditionMG<dim, VectorType, TRANSFER>>
    transfer;

  /**
   * 指示该对象是用一个DoFHandler初始化还是每个区块用一个的标志。
   *
   */
  const bool uses_dof_handler_vector;

  /**
   * 该对象使用的信号
   *
   */
  mg::Signals signals;
};

 /*@}*/ 

#ifndef DOXYGEN
 /* --------------------------- inline functions --------------------- */ 


template <typename VectorType>
Multigrid<VectorType>::Multigrid(const MGMatrixBase<VectorType> &    matrix,
                                 const MGCoarseGridBase<VectorType> &coarse,
                                 const MGTransferBase<VectorType> &  transfer,
                                 const MGSmootherBase<VectorType> &  pre_smooth,
                                 const MGSmootherBase<VectorType> &post_smooth,
                                 const unsigned int                min_level,
                                 const unsigned int                max_level,
                                 Cycle                             cycle)
  : cycle_type(cycle)
  , matrix(&matrix, typeid(*this).name())
  , coarse(&coarse, typeid(*this).name())
  , transfer(&transfer, typeid(*this).name())
  , pre_smooth(&pre_smooth, typeid(*this).name())
  , post_smooth(&post_smooth, typeid(*this).name())
  , edge_out(nullptr, typeid(*this).name())
  , edge_in(nullptr, typeid(*this).name())
  , edge_down(nullptr, typeid(*this).name())
  , edge_up(nullptr, typeid(*this).name())
{
  if (max_level == numbers::invalid_unsigned_int)
    maxlevel = matrix.get_maxlevel();
  else
    maxlevel = max_level;
  reinit(min_level, maxlevel);
}



template <typename VectorType>
inline unsigned int
Multigrid<VectorType>::get_maxlevel() const
{
  return maxlevel;
}



template <typename VectorType>
inline unsigned int
Multigrid<VectorType>::get_minlevel() const
{
  return minlevel;
}


 /* --------------------------- inline functions --------------------- */ 


namespace internal
{
  namespace PreconditionMGImplementation
  {
    template <int dim,
              typename VectorType,
              class TRANSFER,
              typename OtherVectorType>
    typename std::enable_if<TRANSFER::supports_dof_handler_vector>::type
    vmult(
      const std::vector<const dealii::DoFHandler<dim> *> &dof_handler_vector,
      dealii::Multigrid<VectorType> &                     multigrid,
      const TRANSFER &                                    transfer,
      OtherVectorType &                                   dst,
      const OtherVectorType &                             src,
      const bool                          uses_dof_handler_vector,
      const typename dealii::mg::Signals &signals,
      int)
    {
      signals.transfer_to_mg(true);
      if (uses_dof_handler_vector)
        transfer.copy_to_mg(dof_handler_vector, multigrid.defect, src);
      else
        transfer.copy_to_mg(*dof_handler_vector[0], multigrid.defect, src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      if (uses_dof_handler_vector)
        transfer.copy_from_mg(dof_handler_vector, dst, multigrid.solution);
      else
        transfer.copy_from_mg(*dof_handler_vector[0], dst, multigrid.solution);
      signals.transfer_to_global(false);
    }

    template <int dim,
              typename VectorType,
              class TRANSFER,
              typename OtherVectorType>
    void
    vmult(
      const std::vector<const dealii::DoFHandler<dim> *> &dof_handler_vector,
      dealii::Multigrid<VectorType> &                     multigrid,
      const TRANSFER &                                    transfer,
      OtherVectorType &                                   dst,
      const OtherVectorType &                             src,
      const bool                          uses_dof_handler_vector,
      const typename dealii::mg::Signals &signals,
      ...)
    {
      (void)uses_dof_handler_vector;
      Assert(!uses_dof_handler_vector, ExcInternalError());

      signals.transfer_to_mg(true);
      transfer.copy_to_mg(*dof_handler_vector[0], multigrid.defect, src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      transfer.copy_from_mg(*dof_handler_vector[0], dst, multigrid.solution);
      signals.transfer_to_global(false);
    }

    template <int dim,
              typename VectorType,
              class TRANSFER,
              typename OtherVectorType>
    typename std::enable_if<TRANSFER::supports_dof_handler_vector>::type
    vmult_add(
      const std::vector<const dealii::DoFHandler<dim> *> &dof_handler_vector,
      dealii::Multigrid<VectorType> &                     multigrid,
      const TRANSFER &                                    transfer,
      OtherVectorType &                                   dst,
      const OtherVectorType &                             src,
      const bool                          uses_dof_handler_vector,
      const typename dealii::mg::Signals &signals,
      int)
    {
      signals.transfer_to_mg(true);
      if (uses_dof_handler_vector)
        transfer.copy_to_mg(dof_handler_vector, multigrid.defect, src);
      else
        transfer.copy_to_mg(*dof_handler_vector[0], multigrid.defect, src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      if (uses_dof_handler_vector)
        transfer.copy_from_mg_add(dof_handler_vector, dst, multigrid.solution);
      else
        transfer.copy_from_mg_add(*dof_handler_vector[0],
                                  dst,
                                  multigrid.solution);
      signals.transfer_to_global(false);
    }

    template <int dim,
              typename VectorType,
              class TRANSFER,
              typename OtherVectorType>
    void
    vmult_add(
      const std::vector<const dealii::DoFHandler<dim> *> &dof_handler_vector,
      dealii::Multigrid<VectorType> &                     multigrid,
      const TRANSFER &                                    transfer,
      OtherVectorType &                                   dst,
      const OtherVectorType &                             src,
      const bool                          uses_dof_handler_vector,
      const typename dealii::mg::Signals &signals,
      ...)
    {
      (void)uses_dof_handler_vector;
      Assert(!uses_dof_handler_vector, ExcInternalError());

      signals.transfer_to_mg(true);
      transfer.copy_to_mg(*dof_handler_vector[0], multigrid.defect, src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      transfer.copy_from_mg_add(*dof_handler_vector[0],
                                dst,
                                multigrid.solution);
      signals.transfer_to_global(false);
    }
  } // namespace PreconditionMGImplementation
} // namespace internal

template <int dim, typename VectorType, class TRANSFER>
PreconditionMG<dim, VectorType, TRANSFER>::PreconditionMG(
  const DoFHandler<dim> &dof_handler,
  Multigrid<VectorType> &mg,
  const TRANSFER &       transfer)
  : dof_handler_vector(1, &dof_handler)
  , dof_handler_vector_raw(1, &dof_handler)
  , multigrid(&mg)
  , transfer(&transfer)
  , uses_dof_handler_vector(false)
{}

template <int dim, typename VectorType, class TRANSFER>
PreconditionMG<dim, VectorType, TRANSFER>::PreconditionMG(
  const std::vector<const DoFHandler<dim> *> &dof_handler,
  Multigrid<VectorType> &                     mg,
  const TRANSFER &                            transfer)
  : dof_handler_vector(dof_handler.size())
  , dof_handler_vector_raw(dof_handler.size())
  , multigrid(&mg)
  , transfer(&transfer)
  , uses_dof_handler_vector(true)
{
  for (unsigned int i = 0; i < dof_handler.size(); ++i)
    {
      dof_handler_vector[i]     = dof_handler[i];
      dof_handler_vector_raw[i] = dof_handler[i];
    }
}

template <int dim, typename VectorType, class TRANSFER>
inline bool
PreconditionMG<dim, VectorType, TRANSFER>::empty() const
{
  return false;
}

template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult(
  OtherVectorType &      dst,
  const OtherVectorType &src) const
{
  internal::PreconditionMGImplementation::vmult(dof_handler_vector_raw,
                                                *multigrid,
                                                *transfer,
                                                dst,
                                                src,
                                                uses_dof_handler_vector,
                                                this->signals,
                                                0);
}


template <int dim, typename VectorType, class TRANSFER>
IndexSet
PreconditionMG<dim, VectorType, TRANSFER>::locally_owned_range_indices(
  const unsigned int block) const
{
  AssertIndexRange(block, dof_handler_vector.size());
  return dof_handler_vector[block]->locally_owned_dofs();
}


template <int dim, typename VectorType, class TRANSFER>
IndexSet
PreconditionMG<dim, VectorType, TRANSFER>::locally_owned_domain_indices(
  const unsigned int block) const
{
  AssertIndexRange(block, dof_handler_vector.size());
  return dof_handler_vector[block]->locally_owned_dofs();
}



template <int dim, typename VectorType, class TRANSFER>
MPI_Comm
PreconditionMG<dim, VectorType, TRANSFER>::get_mpi_communicator() const
{
  // currently parallel GMG works with parallel triangulations only,
  // so it should be a safe bet to use it to query MPI communicator:
  const Triangulation<dim> &tria = dof_handler_vector[0]->get_triangulation();
  const parallel::TriangulationBase<dim> *ptria =
    dynamic_cast<const parallel::TriangulationBase<dim> *>(&tria);
  Assert(ptria != nullptr, ExcInternalError());
  return ptria->get_communicator();
}



template <int dim, typename VectorType, class TRANSFER>
boost::signals2::connection
PreconditionMG<dim, VectorType, TRANSFER>::connect_transfer_to_mg(
  const std::function<void(bool)> &slot)
{
  return this->signals.transfer_to_mg.connect(slot);
}



template <int dim, typename VectorType, class TRANSFER>
boost::signals2::connection
PreconditionMG<dim, VectorType, TRANSFER>::connect_transfer_to_global(
  const std::function<void(bool)> &slot)
{
  return this->signals.transfer_to_global.connect(slot);
}



template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult_add(
  OtherVectorType &      dst,
  const OtherVectorType &src) const
{
  internal::PreconditionMGImplementation::vmult_add(dof_handler_vector_raw,
                                                    *multigrid,
                                                    *transfer,
                                                    dst,
                                                    src,
                                                    uses_dof_handler_vector,
                                                    this->signals,
                                                    0);
}


template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::Tvmult(OtherVectorType &,
                                                  const OtherVectorType &) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::Tvmult_add(
  OtherVectorType &,
  const OtherVectorType &) const
{
  Assert(false, ExcNotImplemented());
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


