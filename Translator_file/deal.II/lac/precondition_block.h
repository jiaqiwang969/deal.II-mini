//include/deal.II-translator/lac/precondition_block_0.txt
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

#ifndef dealii_precondition_block_h
#define dealii_precondition_block_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/precondition_block_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/


/**
 * 实际块状预处理程序的基类。该类假设<tt>MatrixType</tt>由对角线上的
 * @p blocksize
 * 的可逆块组成，并提供矩阵对角线块的反转。对于这个类来说，矩阵不一定是块状对角线；相反，它适用于任意结构的矩阵，其最小属性是在对角线上有可逆块。但是，该矩阵必须能够访问单个矩阵条目。因此，BlockMatrixArray和类似的类不是一个可能的矩阵类模板参数。
 * 这个类所使用的块状矩阵结构是给定的，例如，对于传输方程的DG方法。对于一个下游的编号，矩阵甚至已经得到了一个块状的左下角矩阵结构，也就是说，矩阵在对角线块的上方是空的。
 *
 *
 * @note
 * 该类旨在用于矩阵，其结构是由来自不相交的单元的局部贡献给出的，例如用于DG方法。它不适用于块结构由不同的物理变量产生的问题，如
 * step-22 中考虑的斯托克斯方程。
 * 对于所有对角线块上方和下方为空的矩阵（即所有块状对角线矩阵），
 * @p BlockJacobi
 * 预调节器是一个直接求解器。对于所有只在对角线块上方为空的矩阵（例如通过下游编号的DG方法得到的矩阵），
 * @p BlockSOR 是一个直接求解器。 @p PreconditionBlock
 * 的第一个实现是假设矩阵有相同块大小的块。如果需要的话，矩阵内不同的块大小仍然必须被实现。
 * 第一个模板参数表示稀疏矩阵中的数字表示类型，第二个表示数字表示类型，在这个类中通过<tt>invert_diagblocks()</tt>存储倒置的对角块矩阵。如果你不想把块反转作为一个精确的求解器，而是作为一个预处理器，你可能想用比原始矩阵更低的精度来存储反转的块；例如，<tt>number==double,
 * inverse_type=float</tt>可能是一个可行的选择。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename MatrixType,
          typename inverse_type = typename MatrixType::value_type>
class PreconditionBlock : public virtual Subscriptor,
                          protected PreconditionBlockBase<inverse_type>
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
  using value_type = inverse_type;

public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 块状预处理程序的参数。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造器。由于没有合理的默认参数，必须给出块的大小。
     *
     */
    AdditionalData(const size_type block_size,
                   const double    relaxation      = 1.,
                   const bool      invert_diagonal = true,
                   const bool      same_diagonal   = false);

    /**
     * 松弛参数。
     *
     */
    double relaxation;

    /**
     * 区块大小。
     *
     */
    size_type block_size;

    /**
     * 在初始化过程中反转对角线。
     *
     */
    bool invert_diagonal;

    /**
     * 假设所有的对角线块都相等以节省内存。
     *
     */
    bool same_diagonal;
    /**
     * 选择块的反转方法。
     *
     */
    typename PreconditionBlockBase<inverse_type>::Inversion inversion;

    /**
     * 如果#反转是SVD，那么低于这个阈值的奇异值将被视为零，从而不被反转。这个参数在调用
     * LAPACKFullMatrix::compute_inverse_svd(). 时使用。
     *
     */
    double threshold;
  };


  /**
   * 构造函数。
   *
   */
  PreconditionBlock(bool store_diagonals = false);

  /**
   * 解构器。
   *
   */
  ~PreconditionBlock() override = default;

  /**
   * 初始化矩阵和块大小。
   * 我们将矩阵和块的大小存储在预处理对象中。在第二步中，可以计算对角线块的倒数。
   * 此外，可以提供派生类的松弛参数。
   *
   */
  void
  initialize(const MatrixType &A, const AdditionalData parameters);

protected:
  /**
   * 初始化矩阵和块的大小，以便进行渗透式预处理。
   * 除了其他initialize()函数的参数外，我们还交出两个索引向量，其中包括包络和它的逆向。关于这些向量的含义，见PreconditionBlockSOR。
   * 在第二步中，可以计算对角线块的倒数。
   * 请确保你使用invert_permuted_diagblocks()来产生一致的数据。
   * 此外，还可以为派生类提供一个放松参数。
   *
   */
  void
  initialize(const MatrixType &            A,
             const std::vector<size_type> &permutation,
             const std::vector<size_type> &inverse_permutation,
             const AdditionalData          parameters);

  /**
   * 根据向量的大小，设置行的permutation或者块的permutation。
   * 如果置换向量的大小等于线性系统的维度，则假定行是单独置换的。在这种情况下，set_permutation()必须在initialize()之前调用，因为对角线块是由矩阵的permuted条目建立的。
   * 如果置换向量的大小不等于系统的维度，对角线块将从未置换的条目中计算出来。
   * 相反，松弛方法step()和Tstep()会按照排列向量给出的顺序应用这些块。如果这个向量的长度不等于块的数量，它们将抛出一个异常。
   * @note
   * 块的排列只能应用于松弛运算符step()和Tstep()，不能应用于预处理运算符vmult()和Tvmult()。
   * @note
   * 在initialize()之前调用set_permutation()是安全的，而其他顺序只允许用于块的permutation。
   *
   */
  void
  set_permutation(const std::vector<size_type> &permutation,
                  const std::vector<size_type> &inverse_permutation);

  /**
   * 替换invert_diagblocks()用于包络预处理。
   *
   */
  void
  invert_permuted_diagblocks();

public:
  /**
   * 如果存在的话，删除反对角线块矩阵，将块大小设置为0，从而使该类处于调用构造函数后的直接状态。
   *
   */
  void
  clear();

  /**
   * 检查对象是否为空。
   *
   */
  bool
  empty() const;

  /**
   * 对条目进行只读访问。这个功能只有在反对角线块被存储的情况下才能实现。
   *
   */
  value_type
  el(size_type i, size_type j) const;

  /**
   * 在 @p inverse. 中存储对角线块的倒数
   * 这需要花费一些额外的内存
   *
   * - 对于DG方法来说，大约是用于矩阵的1/3（用于双倍反演）或1/6（用于浮点反演）。
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

  /**
   * 在正向编号中执行一个块状松弛步骤。    根据参数 @p
   * dst 和 @p pref,
   * ，这将执行一个SOR步骤（都引用同一个向量）或一个Jacobi步骤（都是不同的向量）。对于雅可比步骤，调用函数必须将
   * @p dst 复制到 @p pref 之后。
   * @note
   * 如果设置了一个包络，它将自动被这个函数所尊重。
   *
   */
  template <typename number2>
  void
  forward_step(Vector<number2> &      dst,
               const Vector<number2> &prev,
               const Vector<number2> &src,
               const bool             transpose_diagonal) const;

  /**
   * 在后退编号中执行一个块状放松步骤。    根据参数 @p
   * dst 和 @p pref,
   * ，这将执行一个SOR步骤（都引用同一个向量）的雅可比步骤（都是不同的向量）。对于雅可比步骤，调用函数必须将
   * @p dst 复制到 @p pref 之后。
   * @note
   * 如果设置了一个包络，它将自动被这个函数所尊重。
   *
   */
  template <typename number2>
  void
  backward_step(Vector<number2> &      dst,
                const Vector<number2> &prev,
                const Vector<number2> &src,
                const bool             transpose_diagonal) const;


  /**
   * 返回块的大小。
   *
   */
  size_type
  block_size() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */

  /**
   * 对于非重叠块预处理，块大小必须除以矩阵大小。如果不是，就会抛出这个异常。
   *
   */
  DeclException2(ExcWrongBlockSize,
                 int,
                 int,
                 << "The blocksize " << arg1 << " and the size of the matrix "
                 << arg2 << " do not match.");

  /**
   * 异常
   *
   */
  DeclException0(ExcInverseMatricesAlreadyExist);

  //@}

protected:
  /**
   * 块的大小。每个对角线块都被假定为相同的大小。
   *
   */
  size_type blocksize;

  /**
   * 指向矩阵的指针。确保只要这个类需要，矩阵就一直存在，即直到调用
   * @p invert_diagblocks,
   * 或（如果不应该存储逆矩阵）直到最后一次调用派生类的pre-onditoining
   * @p vmult 函数。
   *
   */
  SmartPointer<const MatrixType, PreconditionBlock<MatrixType, inverse_type>> A;
  /**
   * 派生类要使用的松弛参数。
   *
   */
  double relaxation;

  /**
   * 包容向量
   *
   */
  std::vector<size_type> permutation;

  /**
   * 反向的包络向量
   *
   */
  std::vector<size_type> inverse_permutation;
};



/**
 * 块状雅可比预处理。对矩阵的要求见PreconditionBlock。该类满足 @ref ConceptRelaxationType 的 "放松概念"
 * 。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt> @<float@> 和 @<double@></tt>;
 * 其他人可以在应用程序中生成（见手册中的 @ref
 * Instantiations 部分）。
 *
 *
 */
template <typename MatrixType,
          typename inverse_type = typename MatrixType::value_type>
class PreconditionBlockJacobi
  : public virtual Subscriptor,
    private PreconditionBlock<MatrixType, inverse_type>
{
private:
  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 符合标准的迭代器。
   *
   */
  class const_iterator
  {
  private:
    /**
     * 迭代器的访问器类
     *
     */
    class Accessor
    {
    public:
      /**
       * 构造器。因为我们只使用访问器进行读取访问，一个常数矩阵指针就足够了。
       *
       */
      Accessor(const PreconditionBlockJacobi<MatrixType, inverse_type> *matrix,
               const size_type                                          row);

      /**
       * 这个对象所代表的元素的行号。
       *
       */
      size_type
      row() const;

      /**
       * 这个对象所代表的元素的列号。
       *
       */
      size_type
      column() const;

      /**
       * 这个矩阵条目的值。
       *
       */
      inverse_type
      value() const;

    protected:
      /**
       * 访问的矩阵。
       *
       */
      const PreconditionBlockJacobi<MatrixType, inverse_type> *matrix;

      /**
       * 在这里保存块的大小，以便进一步参考。
       *
       */
      size_type bs;

      /**
       * 当前的区块编号。
       *
       */
      size_type a_block;

      /**
       * 块内的迭代器。
       *
       */
      typename FullMatrix<inverse_type>::const_iterator b_iterator;

      /**
       * 当前块的结束。
       *
       */
      typename FullMatrix<inverse_type>::const_iterator b_end;

      // Make enclosing class a friend.
      friend class const_iterator;
    };

  public:
    /**
     * 构造函数。
     *
     */
    const_iterator(
      const PreconditionBlockJacobi<MatrixType, inverse_type> *matrix,
      const size_type                                          row);

    /**
     * 前缀增量。
     *
     */
    const_iterator &
    operator++();

    /**
     * 去引用操作符。
     *
     */
    const Accessor &operator*() const;

    /**
     * 解除引用操作符。
     *
     */
    const Accessor *operator->() const;

    /**
     * 比较。真，如果两个迭代器都指向同一个矩阵位置。
     *
     */
    bool
    operator==(const const_iterator &) const;
    /**
     * <tt>==</tt>的倒数。
     *
     */
    bool
    operator!=(const const_iterator &) const;

    /**
     * 比较运算符。如果第一行数字较小，或者行数字相等且第一个索引较小，则结果为真。
     *
     */
    bool
    operator<(const const_iterator &) const;

  private:
    /**
     * 存储一个访问器类的对象。
     *
     */
    Accessor accessor;
  };

  /**
   * 从私有基类中导入函数
   *
   */
  using typename PreconditionBlock<MatrixType, inverse_type>::AdditionalData;
  using PreconditionBlock<MatrixType, inverse_type>::initialize;
  using PreconditionBlock<MatrixType, inverse_type>::clear;
  using PreconditionBlock<MatrixType, inverse_type>::empty;
  using PreconditionBlock<MatrixType, inverse_type>::el;
  using PreconditionBlock<MatrixType, inverse_type>::invert_diagblocks;
  using PreconditionBlock<MatrixType, inverse_type>::block_size;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  using PreconditionBlock<MatrixType, inverse_type>::set_permutation;

  /**
   * 执行块状雅可比预处理。
   * 如果存在逆矩阵，该函数将自动使用逆矩阵，如果不存在，那么BlockJacobi将需要很多时间在每个预处理步骤中反转对角线块矩阵。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 与 @p vmult, 相同，因为雅可比是对称的。
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &, const Vector<number2> &) const;
  /**
   * 执行块状雅可比预处理，添加到 @p dst.
   * 这个函数将自动使用反矩阵（如果存在），如果没有，那么BlockJacobi将需要很多时间在每个预处理步骤中反转对角块状矩阵。
   *
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 与 @p vmult_add, 相同，因为雅可比是对称的。
   *
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 执行雅可比迭代的一个步骤。
   *
   */
  template <typename number2>
  void
  step(Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * 执行雅可比迭代的一个步骤。
   *
   */
  template <typename number2>
  void
  Tstep(Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * 迭代器从第一条开始。
   *
   */
  const_iterator
  begin() const;

  /**
   * 最后的迭代器。
   *
   */
  const_iterator
  end() const;

  /**
   * 从第 @p r. 行的第一个条目开始的迭代器。
   *
   */
  const_iterator
  begin(const size_type r) const;

  /**
   * 行 @p r. 的最终迭代器
   *
   */
  const_iterator
  end(const size_type r) const;


private:
  /**
   * 预处理程序的实际实现。    根据 @p adding,
   * ，预处理的结果被添加到目标向量中。
   *
   */
  template <typename number2>
  void
  do_vmult(Vector<number2> &, const Vector<number2> &, bool adding) const;

  friend class Accessor;
  friend class const_iterator;
};



/**
 * 块状SOR预处理。该类满足了 @ref ConceptRelaxationType 的 "放松概念"
 * 。 函数 @p vmult 和 @p Tvmult
 * 根据PreconditionBlock中的块，执行一个（转置的）块-SOR步骤。对角线块之外的元素可以任意分布。
 * 对矩阵的要求见PreconditionBlock。该类中使用的块必须是连续的和不重叠的。一个重叠的施瓦兹松弛方法可以在RelaxationBlockSOR中找到；不过该类不提供预处理。
 * <h3>Permutations</h3>
 * 可选地，源向量的条目可以按照#set_permutation设置的包络向量中的索引顺序来处理（或者对于Tvmult()来说是相反顺序）。反向排列用于将元素存储回这个向量中。在调用set_permutation()的非零大小的向量后，这个功能会自动启用。
 *
 *
 * @note  对角线块，就像矩阵一样，是不被置换的!
 * 因此，互换向量只能交换整个块。它不能改变块内部的顺序或交换块之间的单个索引。
 * <h3>Instantiations</h3>
 *
 * @note
 * 该模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他的可以在应用程序中生成（见手册中的 @ref
 * Instantiations 部分）。
 *
 *
 */
template <typename MatrixType,
          typename inverse_type = typename MatrixType::value_type>
class PreconditionBlockSOR
  : public virtual Subscriptor,
    protected PreconditionBlock<MatrixType, inverse_type>
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 默认构造函数。
   *
   */
  PreconditionBlockSOR();

  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 从受保护的基类中导入类型和函数。
   *
   */
  using typename PreconditionBlock<MatrixType, inverse_type>::AdditionalData;
  using PreconditionBlock<MatrixType, inverse_type>::initialize;
  using PreconditionBlock<MatrixType, inverse_type>::clear;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlock<MatrixType, inverse_type>::invert_diagblocks;
  using PreconditionBlock<MatrixType, inverse_type>::set_permutation;
  using PreconditionBlockBase<inverse_type>::log_statistics;

  /**
   * 执行块SOR预处理。
   * 如果存在逆矩阵，该函数将自动使用逆矩阵，如果不存在，那么BlockSOR将在每个预处理步骤中浪费很多时间来反转对角线块矩阵。
   * 对于对角线块以上为空的矩阵，BlockSOR是一个直接求解器。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 执行块SOR预处理。    警告：这个函数执行正常的 @p vmult
   * ，不加。它存在的原因是BlockMatrixArray默认需要添加版本。另一方面，添加需要一个额外的辅助向量，这并不可取。
   * @see  vmult
   *
   */
  template <typename number2>
  void
  vmult_add(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 向后应用vmult()。
   * 在目前的实现中，这不是vmult()的转置。它是一个应用于整个矩阵的转置的高斯-赛德尔算法，但是被反转的对角线块没有转置。因此，如果对角线块是对称的，它就是转置的。
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 执行后向块SOR预处理。    警告：这个函数执行正常的 @p
   * vmult
   * ，不加。它存在的原因是BlockMatrixArray默认需要添加版本。另一方面，添加需要一个额外的辅助向量，这并不可取。
   * @see  vmult
   *
   */
  template <typename number2>
  void
  Tvmult_add(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 执行SOR迭代的一个步骤。
   *
   */
  template <typename number2>
  void
  step(Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * 执行一步转置的SOR迭代。
   *
   */
  template <typename number2>
  void
  Tstep(Vector<number2> &dst, const Vector<number2> &rhs) const;

protected:
  /**
   * 由PreconditionBlockSSOR使用的构造函数。
   *
   */
  PreconditionBlockSOR(bool store);

  /**
   * 实现由vmult()和vmult_add()调用的正向置换循环。
   * 如果#permutation是由set_permutation()设置的，它将自动被这个函数遵守。
   * 参数 @p adding 还没有任何功能。
   *
   */
  template <typename number2>
  void
  forward(Vector<number2> &,
          const Vector<number2> &,
          const bool transpose_diagonal,
          const bool adding) const;

  /**
   * 实现由Tvmult()和Tvmult_add()调用的后向替换循环。
   * 如果一个#permutation是由set_permutation()设置的，它将自动被这个函数兑现。
   * 参数 @p adding 还没有任何功能。
   *
   */
  template <typename number2>
  void
  backward(Vector<number2> &,
           const Vector<number2> &,
           const bool transpose_diagonal,
           const bool adding) const;
};


/**
 * 块状SSOR预处理。该类满足了 @ref ConceptRelaxationType 的 "放松概念"
 * 。 函数 @p vmult 和 @p Tvmult
 * 根据PreconditionBlockSOR中的实现，执行一个块-SSOR步骤。
 * 这个类需要存储对角线块和它们的逆向值。
 * 关于矩阵的要求见PreconditionBlock。该类中使用的块必须是连续的和不重叠的。一个重叠的施瓦兹松弛方法可以在RelaxationBlockSSOR中找到；不过该类不提供预处理。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他的可以在应用程序中生成（见手册中 @ref Instantiations
 * 一节）。
 *
 *
 */
template <typename MatrixType,
          typename inverse_type = typename MatrixType::value_type>
class PreconditionBlockSSOR
  : public virtual Subscriptor,
    private PreconditionBlockSOR<MatrixType, inverse_type>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 定义矩阵的数字类型。
   *
   */
  using number = typename MatrixType::value_type;

  /**
   * 构造函数。
   *
   */
  PreconditionBlockSSOR();

  // Keep AdditionalData accessible
  using typename PreconditionBlockSOR<MatrixType, inverse_type>::AdditionalData;

  // The following are the
  // functions of the base classes
  // which we want to keep
  // accessible.
  /**
   * 使初始化函数公开可用。
   *
   */
  using PreconditionBlockSOR<MatrixType, inverse_type>::initialize;
  using PreconditionBlockSOR<MatrixType, inverse_type>::clear;
  using PreconditionBlockBase<inverse_type>::size;
  using PreconditionBlockBase<inverse_type>::inverse;
  using PreconditionBlockBase<inverse_type>::inverse_householder;
  using PreconditionBlockBase<inverse_type>::inverse_svd;
  using PreconditionBlockBase<inverse_type>::log_statistics;
  using PreconditionBlockSOR<MatrixType, inverse_type>::set_permutation;
  using PreconditionBlockSOR<MatrixType, inverse_type>::empty;
  using PreconditionBlockSOR<MatrixType, inverse_type>::el;
  using PreconditionBlockSOR<MatrixType, inverse_type>::invert_diagblocks;

  /**
   * 执行块SSOR预处理。
   * 如果存在逆矩阵，该函数将自动使用逆矩阵，如果不存在，那么BlockSOR将在每个预处理步骤中浪费很多时间来反转对角线块矩阵。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 与vmult()相同
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &, const Vector<number2> &) const;

  /**
   * 执行SOR迭代的一个步骤。
   *
   */
  template <typename number2>
  void
  step(Vector<number2> &dst, const Vector<number2> &rhs) const;

  /**
   * 执行一步转置的SOR迭代。
   *
   */
  template <typename number2>
  void
  Tstep(Vector<number2> &dst, const Vector<number2> &rhs) const;
};

 /*@}*/ 
//---------------------------------------------------------------------------

#ifndef DOXYGEN

template <typename MatrixType, typename inverse_type>
inline bool
PreconditionBlock<MatrixType, inverse_type>::empty() const
{
  if (A == nullptr)
    return true;
  return A->empty();
}


template <typename MatrixType, typename inverse_type>
inline inverse_type
PreconditionBlock<MatrixType, inverse_type>::el(size_type i, size_type j) const
{
  const size_type    bs = blocksize;
  const unsigned int nb = i / bs;

  const FullMatrix<inverse_type> &B = this->inverse(nb);

  const size_type ib = i % bs;
  const size_type jb = j % bs;

  if (jb + nb * bs != j)
    {
      return 0.;
    }

  return B(ib, jb);
}

//---------------------------------------------------------------------------

template <typename MatrixType, typename inverse_type>
inline PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
  Accessor::Accessor(
    const PreconditionBlockJacobi<MatrixType, inverse_type> *matrix,
    const size_type                                          row)
  : matrix(matrix)
  , bs(matrix->block_size())
  , a_block(row / bs)
  , b_iterator(&matrix->inverse(0), 0, 0)
  , b_end(&matrix->inverse(0), 0, 0)
{
  // This is the end accessor, which
  // does not have a valid block.
  if (a_block == matrix->size())
    return;

  const size_type r = row % bs;

  b_iterator = matrix->inverse(a_block).begin(r);
  b_end      = matrix->inverse(a_block).end();

  AssertIndexRange(a_block, matrix->size());
}


template <typename MatrixType, typename inverse_type>
inline typename PreconditionBlockJacobi<MatrixType, inverse_type>::size_type
PreconditionBlockJacobi<MatrixType,
                        inverse_type>::const_iterator::Accessor::row() const
{
  Assert(a_block < matrix->size(), ExcIteratorPastEnd());

  return bs * a_block + b_iterator->row();
}


template <typename MatrixType, typename inverse_type>
inline typename PreconditionBlockJacobi<MatrixType, inverse_type>::size_type
PreconditionBlockJacobi<MatrixType,
                        inverse_type>::const_iterator::Accessor::column() const
{
  Assert(a_block < matrix->size(), ExcIteratorPastEnd());

  return bs * a_block + b_iterator->column();
}


template <typename MatrixType, typename inverse_type>
inline inverse_type
PreconditionBlockJacobi<MatrixType,
                        inverse_type>::const_iterator::Accessor::value() const
{
  Assert(a_block < matrix->size(), ExcIteratorPastEnd());

  return b_iterator->value();
}


template <typename MatrixType, typename inverse_type>
inline PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
  const_iterator(
    const PreconditionBlockJacobi<MatrixType, inverse_type> *matrix,
    const size_type                                          row)
  : accessor(matrix, row)
{}


template <typename MatrixType, typename inverse_type>
inline
  typename PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator &
  PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
  operator++()
{
  Assert(*this != accessor.matrix->end(), ExcIteratorPastEnd());

  ++accessor.b_iterator;
  if (accessor.b_iterator == accessor.b_end)
    {
      ++accessor.a_block;

      if (accessor.a_block < accessor.matrix->size())
        {
          accessor.b_iterator =
            accessor.matrix->inverse(accessor.a_block).begin();
          accessor.b_end = accessor.matrix->inverse(accessor.a_block).end();
        }
    }
  return *this;
}


template <typename MatrixType, typename inverse_type>
inline const typename PreconditionBlockJacobi<MatrixType, inverse_type>::
  const_iterator::Accessor &
    PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
    operator*() const
{
  return accessor;
}


template <typename MatrixType, typename inverse_type>
inline const typename PreconditionBlockJacobi<MatrixType, inverse_type>::
  const_iterator::Accessor *
    PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
    operator->() const
{
  return &accessor;
}


template <typename MatrixType, typename inverse_type>
inline bool
PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
operator==(const const_iterator &other) const
{
  if (accessor.a_block == accessor.matrix->size() &&
      accessor.a_block == other.accessor.a_block)
    return true;

  if (accessor.a_block != other.accessor.a_block)
    return false;

  return (accessor.row() == other.accessor.row() &&
          accessor.column() == other.accessor.column());
}


template <typename MatrixType, typename inverse_type>
inline bool
PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
operator!=(const const_iterator &other) const
{
  return !(*this == other);
}


template <typename MatrixType, typename inverse_type>
inline bool
PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator::
operator<(const const_iterator &other) const
{
  return (accessor.row() < other.accessor.row() ||
          (accessor.row() == other.accessor.row() &&
           accessor.column() < other.accessor.column()));
}


template <typename MatrixType, typename inverse_type>
inline
  typename PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator
  PreconditionBlockJacobi<MatrixType, inverse_type>::begin() const
{
  return const_iterator(this, 0);
}


template <typename MatrixType, typename inverse_type>
inline
  typename PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator
  PreconditionBlockJacobi<MatrixType, inverse_type>::end() const
{
  return const_iterator(this, this->size() * this->block_size());
}


template <typename MatrixType, typename inverse_type>
inline
  typename PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator
  PreconditionBlockJacobi<MatrixType, inverse_type>::begin(
    const size_type r) const
{
  AssertIndexRange(r, this->A->m());
  return const_iterator(this, r);
}



template <typename MatrixType, typename inverse_type>
inline
  typename PreconditionBlockJacobi<MatrixType, inverse_type>::const_iterator
  PreconditionBlockJacobi<MatrixType, inverse_type>::end(
    const size_type r) const
{
  AssertIndexRange(r, this->A->m());
  return const_iterator(this, r + 1);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


