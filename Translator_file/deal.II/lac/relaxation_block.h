//include/deal.II-translator/lac/relaxation_block_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_relaxation_block_h
#define dealii_relaxation_block_h

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/precondition_block_base.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 用于实现重叠、乘法施瓦兹松弛方法和平滑器的基类。
 * 该类使用PreconditionBlockBase提供的基础设施。它增加了用块列表初始化和做松弛步骤的函数。实际的松弛方法与SolverRelaxation和MGSmootherRelaxation所期望的接口在派生类中。
 * 这个类允许比PreconditionBlock更通用的松弛方法，因为索引集可以是任意的和重叠的，而那里只允许大小相等的连续的、不相交的集合。作为一个缺点，这个类不能作为预处理程序使用，因为它的实现依赖于高斯-塞德尔过程的直接实现。
 * 并行计算要求你在 AdditionalData::temp_ghost_vector.
 * 中指定一个初始化的幽灵向量。
 *
 *
 * @ingroup Preconditioners
 *
 *
 */
template <typename MatrixType,
          typename InverseNumberType = typename MatrixType::value_type,
          typename VectorType        = Vector<double>>
class RelaxationBlock : protected PreconditionBlockBase<InverseNumberType>
{
private:
  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 逆矩阵的数值类型。
   *
   */
  using value_type = InverseNumberType;

public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 块状松弛方法的参数。除了像#relaxation这样的典型控制参数外，这个对象还包含#block_list中的块结构和#order中的块的可选排序。
   *
   */
  class AdditionalData : public Subscriptor
  {
  public:
    /**
     * 构造函数。
     *
     */
    AdditionalData(
      const double relaxation      = 1.,
      const bool   invert_diagonal = true,
      const bool   same_diagonal   = false,
      const typename PreconditionBlockBase<InverseNumberType>::Inversion
                   inversion = PreconditionBlockBase<InverseNumberType>::gauss_jordan,
      const double threshold         = 0.,
      VectorType * temp_ghost_vector = nullptr);

    /**
     * 从索引到块的映射。这个模式的每一行都列举了构成要反转的对角线块的索引。
     *
     */
    SparsityPattern block_list;

    /**
     * 放松参数。
     *
     */
    double relaxation;

    /**
     * 在初始化过程中反转对角线。另外，对角线块在使用时也会被反转。虽然提前反转块需要更多的内存，但通常可以节省大量的计算。
     * 参见#same_diagonal关于如何避免内存开销的问题。
     *
     */
    bool invert_diagonal;

    /**
     * 假设所有对角线块都相等，以节省内存。如果这个标志为真，那么只有矩阵的第一个对角线块被反转和存储。然后再用于所有其他块。
     * \注
     * 如果你的块不相等，特别是它们的大小不同，则避免设置为真。
     *
     */
    bool same_diagonal;

    /**
     * 选择块的反转方法。
     *
     */
    typename PreconditionBlockBase<InverseNumberType>::Inversion inversion;

    /**
     * 如果#反转是SVD，我们可以计算块的Penrose-Moore反转。为了做到这一点，我们可以在这里指定一个阈值，低于这个阈值的奇异值将被视为零，从而不被反转。
     * 将此参数设置为大于0的值比阈值优先，即如果你想使用阈值，kernel_size必须为0。
     * 这个参数在调用 LAPACKFullMatrix::compute_inverse_svd().
     * 时使用。
     *
     */
    double threshold = 0.;

    /**
     * 如果#inversion是SVD，我们可以计算块的Penrose-Moore逆。为了做到这一点，我们可以在这里指定不被反转而被视为零的内核的大小。将这个参数设置为大于零的值比阈值优先，也就是说，如果你想使用阈值，kernel_size必须是零。
     * 这个参数在调用 LAPACKFullMatrix::compute_inverse_svd().
     * 时使用。
     *
     */
    unsigned int kernel_size = 0;

    /**
     * 块应该被遍历的顺序。这个向量可以启动几种执行模式。          <ol>   <li>  如果向量的长度为零，那么放松方法将从第一个到最后一个块执行。 </li>   <li>  如果长度为1，那么内向量的大小必须与块的数量相同。放宽方法是按照这个向量中给出的顺序来应用的。 </li>   <li>  如果外向量的长度大于1，那么松弛方法将被多次应用，每次都按照相应索引的内向量给出的顺序。例如，这种模式可用于ADI方法和类似的方向扫频。 </li>   </ol> 。
     *
     */
    std::vector<std::vector<unsigned int>> order;

    /**
     * 临时的鬼魂向量，在执行并行MPI计算时用于放松方法。用户需要让它指向一个初始化的向量，该向量包含所有出现在
     * @p block_list sa
     * ghost值中的指数。通常情况下，这就是本地活动层DoF的集合。当VectorType是一个像Vector<double>这样的串行向量类型时未被使用。
     *
     */
    mutable VectorType *temp_ghost_vector;

    /**
     * 返回这个对象中分配的内存。
     *
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * 初始化矩阵和附加信息。在第二步，可以计算对角线块的倒数。
   * 请注意，与其他预处理程序不同，AdditionalData定义了相当大的对象，因此该对象不是复制的，而是存储一个指针。因此，
   * <code>additional_data</code>
   * 的寿命急于超过这个对象的寿命。
   *
   */
  void
  initialize(const MatrixType &A, const AdditionalData &parameters);

  /**
   * 删除逆对角线块矩阵（如果存在的话），将块大小设置为0，从而使该类处于调用构造函数后的直接状态。
   *
   */
  void
  clear();

  /**
   * 在 @p inverse. 中存储对角线块的逆值
   * 这需要花费一些额外的内存
   *
   * - 对于DG方法来说，大约是用于矩阵的1/3（对于双倍反转）或1/6（对于浮动反转）。
   *
   * 但它使预处理的速度大大加快。
   * 在调用<tt>clear(...)</tt>之前，不允许两次调用这个函数（会产生一个错误），因为在第二次调用时，已经存在逆矩阵。
   * 在这个函数被调用后，通过 @p use_matrix
   * 函数给出的矩阵的锁被释放，也就是说，你可以覆盖或删除它。
   * 你可能想这样做，以防你用这个矩阵作为另一个矩阵的前提条件。
   *
   */
  void
  invert_diagblocks();

protected:
  /**
   * 执行一个块状松弛步骤。    根据参数 @p dst 和 @p pref,
   * ，这将执行一个SOR步骤（两者都引用同一个向量）或一个Jacobi步骤（两者都是不同的向量）。对于雅可比步骤，调用函数必须在此后将
   * @p dst 复制到 @p prev 。
   *
   */
  void
  do_step(VectorType &      dst,
          const VectorType &prev,
          const VectorType &src,
          const bool        backward) const;

  /**
   * 指向矩阵的指针。确保只要这个类需要，矩阵就存在，即直到调用
   * @p invert_diagblocks,
   * 或（如果不应该存储逆矩阵）直到派生类的预处理 @p
   * vmult 函数的最后一次调用。
   *
   */
  SmartPointer<const MatrixType,
               RelaxationBlock<MatrixType, InverseNumberType, VectorType>>
    A;

  /**
   * 控制信息。
   *
   */
  SmartPointer<const AdditionalData,
               RelaxationBlock<MatrixType, InverseNumberType, VectorType>>
    additional_data;

private:
  /**
   * 计算（反）一个区块的范围。
   *
   */
  void
  block_kernel(const size_type block_begin, const size_type block_end);
};


/**
 * 块状雅可比（加法施瓦茨）方法，可能有重叠的块。
 * 该类实现了 @ref ConceptRelaxationType "放松概念 "
 * 所期望的step()和Tstep()函数。它们对AdditionalData的块列表中提供的块执行加法施瓦茨方法。与PreconditionBlockJacobi不同的是，这些块可以是不同大小的、非连续的和重叠的。另一方面，这个类并没有实现Solver对象所期望的预处理程序接口。
 *
 *
 * @ingroup Preconditioners
 *
 *
 */
template <typename MatrixType,
          typename InverseNumberType = typename MatrixType::value_type,
          typename VectorType        = Vector<double>>
class RelaxationBlockJacobi
  : public virtual Subscriptor,
    protected RelaxationBlock<MatrixType, InverseNumberType, VectorType>
{
public:
  /**
   * 默认构造函数。
   *
   */
  //    RelaxationBlockJacobi();

  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 使类型公开化。
   *
   */
  using typename RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    AdditionalData;

  /**
   * 公开初始化函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::initialize;

  /**
   * 让基类的函数再次公开。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::clear;

  /**
   * 再次公开基类的函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::size;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    inverse_householder;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse_svd;
  /**
   * 再次公开基类的功能。
   *
   */
  using PreconditionBlockBase<InverseNumberType>::log_statistics;
  /**
   * 执行雅可比迭代的一个步骤。
   *
   */
  void
  step(VectorType &dst, const VectorType &rhs) const;

  /**
   * 执行雅各比迭代的一个步骤。
   *
   */
  void
  Tstep(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现vmult()操作，对于这个类来说，在调用step()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现一个转置vmult操作，对于这个类，在调用Tstep()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  Tvmult(VectorType &dst, const VectorType &rhs) const;
};


/**
 * 块状Gauss-Seidel方法，可能有重叠的块。
 * 该类实现了 @ref ConceptRelaxationType "放松概念 "
 * 所期望的step()和Tstep()函数。它们对AdditionalData的块列表中提供的块执行乘法施瓦茨方法。
 * 与PreconditionBlockSOR不同的是，这些块可以是不同大小的、不连续的和重叠的。另一方面，该类没有实现Solver对象所期望的预处理程序接口。
 *
 *
 * @ingroup Preconditioners
 *
 *
 */
template <typename MatrixType,
          typename InverseNumberType = typename MatrixType::value_type,
          typename VectorType        = Vector<double>>
class RelaxationBlockSOR
  : public virtual Subscriptor,
    protected RelaxationBlock<MatrixType, InverseNumberType, VectorType>
{
public:
  /**
   * 默认构造函数。
   *
   */
  //    RelaxationBlockSOR();

  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 使类型公开化。
   *
   */
  using typename RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    AdditionalData;

  /**
   * 公开初始化函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::initialize;

  /**
   * 让基类的函数再次公开。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::clear;

  /**
   * 再次公开基类的函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::size;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    inverse_householder;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse_svd;
  /**
   * 再次公开基类的功能。
   *
   */
  using PreconditionBlockBase<InverseNumberType>::log_statistics;
  /**
   * 执行SOR迭代的一个步骤。
   *
   */
  void
  step(VectorType &dst, const VectorType &rhs) const;

  /**
   * 执行一步转置的SOR迭代。
   *
   */
  void
  Tstep(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现vmult()操作，对于这个类，在调用step()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现一个转置vmult操作，对于这个类，在调用Tstep()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  Tvmult(VectorType &dst, const VectorType &rhs) const;
};


/**
 * 对称块高斯-赛德尔方法，可能有重叠的块。
 * 该类实现了 @ref ConceptRelaxationType "放松概念 "
 * 所期望的step()和Tstep()函数。它们以对称的方式对AdditionalData的块列表中提供的块执行乘法施瓦茨方法。与PreconditionBlockSSOR不同的是，这些块可以是不同大小的、不连续的和重叠的。另一方面，这个类并没有实现Solver对象所期望的预处理接口。
 *
 *
 * @ingroup Preconditioners
 *
 *
 */
template <typename MatrixType,
          typename InverseNumberType = typename MatrixType::value_type,
          typename VectorType        = Vector<double>>
class RelaxationBlockSSOR
  : public virtual Subscriptor,
    protected RelaxationBlock<MatrixType, InverseNumberType, VectorType>
{
public:
  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 使类型公开化。
   *
   */
  using typename RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    AdditionalData;

  /**
   * 公开初始化函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::initialize;

  /**
   * 让基类的函数再次公开。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::clear;

  /**
   * 再次公开基类的函数。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::size;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::
    inverse_householder;
  /**
   * 再次公开基类的功能。
   *
   */
  using RelaxationBlock<MatrixType, InverseNumberType, VectorType>::inverse_svd;
  /**
   * 再次公开基类的功能。
   *
   */
  using PreconditionBlockBase<InverseNumberType>::log_statistics;
  /**
   * 执行SSOR迭代的一个步骤。
   *
   */
  void
  step(VectorType &dst, const VectorType &rhs) const;

  /**
   * 执行一步转置的SSOR迭代。
   *
   */
  void
  Tstep(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现vmult()操作，对于这个类来说，在调用step()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  vmult(VectorType &dst, const VectorType &rhs) const;

  /**
   * 实现一个转置vmult操作，对于这个类，在调用Tstep()方法之前，首先将dst()向量设置为零。
   *
   */
  void
  Tvmult(VectorType &dst, const VectorType &rhs) const;
};


DEAL_II_NAMESPACE_CLOSE

#endif


