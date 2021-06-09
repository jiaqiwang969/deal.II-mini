//include/deal.II-translator/lac/sparse_vanka_0.txt
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

#ifndef dealii_sparse_vanka_h
#define dealii_sparse_vanka_h



#include <deal.II/base/config.h>

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/smartpointer.h>

#include <map>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;
template <typename number>
class Vector;

template <typename number>
class SparseVanka;
template <typename number>
class SparseBlockVanka;
#endif

/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/

/**
 * 点对点的Vanka预处理。该类在一个点-明智的基础上进行Vanka预处理。Vanka预处理用于鞍点问题，如Stokes问题或优化中出现的拉格朗日乘数和牛顿方法矩阵有一个零块的问题。对于这些矩阵，应用Jacobi或Gauss-Seidel方法是不可能的，因为在拉格朗日乘子的行中，一些对角线元素是零。Vanka的方法是为每个朗格朗日乘数变量解决一个小的（通常是不确定的）方程组（我们也将把斯托克斯方程中的压力称为朗格朗日乘数，因为它可以被解释为朗格朗日乘数）。
 * 这个类的对象是通过传递拉格朗日乘子的自由度指数的向量来构造的。在实际的预处理方法中，这些行是按照在矩阵中出现的顺序来遍历的。由于这是一个类似于Gauß-Seidel的程序，记得事先要有一个好的排序（对于传输主导的问题，如果选择流入边界上的点作为重新编号的起点，Cuthill-McKee算法是一个很好的手段）。
 * 对于每个选定的自由度，一个局部方程组由自由度本身和所有其他数值立即耦合建立，即自由度
 * @p i 的局部方程组考虑的自由度集合是 @p i 本身和所有 @p
 * j
 * ，使元素<tt>（i,j）</tt>是考虑中的稀疏矩阵的非零条目。元素<tt>(j,i)</tt>不被考虑。我们现在从刚才描述的自由度集合中的行和列中挑选出所有的矩阵条目，并将其放入一个局部矩阵中，随后将其反转。这个系统对于每个自由度可能有不同的大小，例如取决于计算网格上各个节点的局部邻域。
 * 右手边是以同样的方式建立的，即通过复制所有与目前考虑的那个耦合的条目，但它被所有与上述集合的自由度耦合的自由度所增强（即与目前那个耦合的二阶DoFs）。其原因是，要解决的局部问题在二阶耦合的自由度上有迪里希特边界条件，所以我们必须考虑到它们，但在实际求解前要消除它们；这种消除是通过修改右手边来完成的，最终这些自由度不再出现在矩阵和解矢量中。
 * 这个局部系统被解决了，其值被更新到目标向量中。
 * 备注：凡卡方法是一种非对称性的预处理方法。
 *
 *  <h3>Example of Use</h3>
 * 这个小例子取自一个做参数优化的程序。拉格朗日乘数是所使用的有限元的第三部分。该系统用GMRES方法求解。
 *
 * @code
 *  // tag the Lagrange multiplier variable
 *  vector<bool> signature(3);
 *  signature[0] = signature[1] = false;
 *  signature[2] = true;
 *
 *  // tag all dofs belonging to the Lagrange multiplier
 *  vector<bool> selected_dofs (dof.n_dofs(), false);
 *  DoFTools::extract_dofs(dof, signature, p_select);
 *  // create the Vanka object
 *  SparseVanka<double> vanka (global_matrix, selected_dofs);
 *
 *  // create the solver
 *  SolverGMRES<> gmres(control,memory,504);
 *
 *  // solve
 *  gmres.solve (global_matrix, solution, right_hand_side,
 *               vanka);
 * @endcode
 *
 *  <h4>Implementor's remark</h4>
 * 目前，局部矩阵是这样建立的：与局部拉格朗日乘数相关的自由度是第一个。因此，通常局部矩阵的左上角条目是零。我（W.B.）不清楚这是否会给局部矩阵的反演带来一些问题。也许有人会想检查一下这个问题。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他可以在应用程序中生成（见手册中的 @ref Instantiations
 * 部分）。
 *
 *
 */
template <typename number>
class SparseVanka
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 构造函数。什么都不做。
   * 在使用此对象作为预处理（vmult()）之前，调用initialize()函数。
   *
   */
  SparseVanka();

  /**
   * 构造函数也需要两个被废弃的输入。      @deprecated
   * 最后两个参数的使用已被废弃。它们目前被忽略。
   *
   */
  DEAL_II_DEPRECATED
  SparseVanka(const SparseMatrix<number> &M,
              const std::vector<bool> &   selected,
              const bool                  conserve_memory,
              const unsigned int n_threads = MultithreadInfo::n_threads());

  /**
   * 构造函数。获取用于预处理的矩阵和一个比特向量，其条目
   * @p true
   * 为所有要更新的行。对这个向量的引用将被存储，所以它必须比Vanka对象的持续时间长。矩阵的情况也是如此。
   * 这里传递的矩阵 @p M
   * 可能是也可能不是这个对象应作为预处理程序的同一矩阵。特别是，可以想象预处理器只为一个矩阵建立一次，但也用于非线性过程的后续步骤，其中矩阵在每一步中都会有轻微变化。
   *
   */
  SparseVanka(const SparseMatrix<number> &M, const std::vector<bool> &selected);

  /**
   * 解构器。删除所有分配的矩阵。
   *
   */
  ~SparseVanka();

  /**
   * SparseVanka的参数。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造函数。关于参数的描述，见下文。
     *
     */
    explicit AdditionalData(const std::vector<bool> &selected);

    /**
     * 构造函数。关于参数的描述，请见下文。
     * @deprecated  该构造函数的使用已被废弃。
     *
     * - 第二个和第三个参数被忽略。
     *
     */
    DEAL_II_DEPRECATED
    AdditionalData(const std::vector<bool> &selected,
                   const bool               conserve_memory,
                   const unsigned int n_threads = MultithreadInfo::n_threads());

    /**
     * 我们要处理的那些自由度的索引。
     *
     */
    const std::vector<bool> &selected;
  };


  /**
   * 如果使用默认的构造函数，那么这个函数需要在这个类的对象被用作预处理之前被调用。
   * 关于可能的参数的更多细节，请参阅类的文档和
   * SparseVanka::AdditionalData 类的文档。
   * 这个函数被调用后，预处理程序就可以被使用了（使用派生类的
   * <code>vmult</code> 函数）。
   *
   */
  void
  initialize(const SparseMatrix<number> &M,
             const AdditionalData &      additional_data);

  /**
   * 做预处理。这个函数接收 @p src 中的残差并返回 @p dst.
   * 中的结果更新向量。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * 应用转置预处理。这个函数接收 @p 中的残差，并返回 @p
   * dst. 中的结果更新向量。
   *
   */
  template <typename number2>
  void
  Tvmult(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * 返回共域（或范围）空间的维度。注意，矩阵的维度是
   * $m \times n$  .
   * @note
   * 只有在预处理程序已经被初始化的情况下才可以调用这个函数。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维数。请注意，矩阵的维度是 $m \times n$
   * 。
   * @note
   * 只有在预处理程序被初始化的情况下才可以调用这个函数。
   *
   */
  size_type
  n() const;

protected:
  /**
   * 将那些在 @p dof_mask, 中具有 @p true
   * 值的自由度所对应的倒数应用于 @p src
   * 向量，并将结果移入 @p dst.
   * 实际上，只有允许指数的值被写入 @p
   * ]dst，所以如果给定的掩码将所有的值设置为0，那么这个函数的应用只做一般文档中所宣布的事情。提供掩码的原因是，在派生类中，我们可能只想将预处理程序应用于矩阵的一部分，以便将应用并行化。那么，为了消除线程之间的相互依赖，只写到
   * @p dst, 的一些片断是很重要的。
   * 如果传递一个空指针而不是指向 @p dof_mask
   * 的指针（如默认值），那么就假定我们将在所有自由度上工作。这相当于用<tt>向量<bool>(n_dofs,true)</tt>来调用该函数。
   * 这个类的 @p vmult 当然是用一个空指针来调用这个函数。
   *
   */
  template <typename number2>
  void
  apply_preconditioner(Vector<number2> &              dst,
                       const Vector<number2> &        src,
                       const std::vector<bool> *const dof_mask = nullptr) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 指向矩阵的指针。
   *
   */
  SmartPointer<const SparseMatrix<number>, SparseVanka<number>> matrix;

  /**
   * 我们要处理的那些自由度的索引。
   *
   */
  const std::vector<bool> *selected;

  /**
   * 反矩阵的数组，每个自由度有一个。只有那些在 @p
   * selected. 中被标记的元素才会被使用。
   *
   */
  mutable std::vector<SmartPointer<FullMatrix<float>, SparseVanka<number>>>
    inverses;

  /**
   * 范围空间的维度。
   *
   */
  size_type _m;

  /**
   * 域空间的维度。
   *
   */
  size_type _n;

  /**
   * 计算所有选定的对角线元素的反值。
   *
   */
  void
  compute_inverses();

  /**
   * 计算在<tt>[begin,end)</tt>范围内的位置的反演。在非多线程模式下，<tt>compute_inverses()</tt>用整个范围来调用这个函数，但在多线程模式下，这个函数的几个副本被生成。
   *
   */
  void
  compute_inverses(const size_type begin, const size_type end);

  /**
   * 计算位于 @p row. 位置的块的逆值
   * 由于向量被经常使用，它只在这个函数的调用者中生成一次，并传递给这个首先清除它的函数。与本函数每次重新创建向量的情况相比，重复使用向量使整个过程明显加快。
   *
   */
  void
  compute_inverse(const size_type row, std::vector<size_type> &local_indices);

  // Make the derived class a friend. This seems silly, but is actually
  // necessary, since derived classes can only access non-public members
  // through their @p this pointer, but not access these members as member
  // functions of other objects of the type of this base class (i.e. like
  // <tt>x.f()</tt>, where @p x is an object of the base class, and @p f one
  // of it's non-public member functions).
  //
  // Now, in the case of the @p SparseBlockVanka class, we would like to take
  // the address of a function of the base class in order to call it through
  // the multithreading framework, so the derived class has to be a friend.
  template <typename T>
  friend class SparseBlockVanka;
};



/**
 * 稀疏Vanka预处理器的块状版本。这个类将矩阵分成块，只在对角线块上工作，这当然会降低作为预处理程序的效率，但完全可以并行。构造函数需要一个参数，将矩阵分成多少个块，然后让底层的类来做这个工作。矩阵的划分有几种方式，下面将详细介绍。
 * 如果你没有一个多处理器系统，这个类可能是无用的，因为那时每个预处理步骤的工作量与
 * @p SparseVanka
 * 类相同，但预处理特性更差。另一方面，如果你有一个多处理器系统，较差的预处理质量（导致线性求解器的更多迭代）通常会被并行化带来的应用速度的提高所平衡，从而导致解决线性系统的耗时减少。应该注意的是，作为预处理程序的质量会随着块数的增加而降低，所以块数可能有一个最佳值（就每次线性求解的壁挂时间而言）。
 * 为了方便编写可移植的代码，如果将矩阵细分的块数设置为1，那么这个类的作用就像
 * @p SparseVanka
 * 类。因此，你可能想把块的数量设置为等于你拥有的处理器的数量。
 * 请注意，如果<tt>deal.II</tt>被配置为多线程使用，并且被生成的线程数等于块的数量，那么并行化就会完成。这是合理的，因为你不会想把块的数量设置得过大，因为如前所述，这降低了预处理的特性。
 *
 *  <h3>Splitting the matrix into blocks</h3>
 * 将矩阵分割成块的方式总是使块的大小不一定相等，但块与块之间要解决的局部系统的选定自由度的数量是相等的。这种细分策略的原因是多线程的负载平衡。有几种可能性可以将矩阵实际分割成块，由传递给构造函数的标志
 * @p
 * blocking_strategy来选择。在下文中，我们将用块来表示自由度的索引列表；算法将在每个块上单独工作，也就是说，一个块的自由度所对应的局部系统的解将只用于更新属于同一块的自由度，而不会用于更新其他块的自由度。一个区块可以是一个连续的指数列表，如下面的第一个方案，也可以是一个不连续的指数列表。当然，我们假设每两个块的交点是空的，所有块的联合等于区间<tt>[0,N)</tt>，其中
 * @p N 是方程组的自由度数。
 * <ul>   <li>   @p index_intervals:  这里，我们选择区块为区间<tt>[a_i,a_{i+1</tt>)}，即连续的自由度通常也在同一区块内。这是一个合理的策略，如果自由度使用Cuthill-McKee算法进行重新编号，其中空间上相邻的自由度有相邻的指数。在这种情况下，矩阵中的耦合通常也被限制在对角线附近，我们可以简单地将矩阵切割成块。
 * 区间的边界，即上面的 @p a_i
 * ，是这样选择的，即我们将在其中工作的自由度数量（即通常与拉格朗日乘数对应的自由度）在每个区块中大致相同；然而，这并不意味着区块的大小相等，因为区块还包括没有解决局部系统的其他自由度。在极端情况下，考虑到所有拉格朗日乘数都被排序到自由度指数范围的末端，那么第一个区块将非常大，因为它包括所有其他自由度和一些拉格朗日乘数，而所有其他区块则相当小，只包括朗格朗日乘数。因此，这一策略不仅取决于拉格朗日做功因子的排序，也取决于其他做功因子的排序。因此有必要指出，如果自由度是按分量编号的，即所有拉格朗日乘法器都是一整块的，这几乎使作为预处理的能力失去作用。
 * <li>   @p adaptive:
 * 在拉格朗日自由度被分组的情况下，这个策略就比较聪明了，就像上面的例子。它的工作原理如下：它首先将Lagrange
 * DoFs分组为块，使用与上述相同的策略。然而，它不是将其他DoF分组到具有最接近DoF索引的Lagrange
 * DoF块中，而是决定将每个非Lagrange
 * DoF放入最经常写入该非Lagrange DoF的Lagrange
 * DoF块中。这使得它甚至可以将拉格朗日
 * DoFs排序到最后，并且仍然将空间上相邻的非拉格朗日
 * DoFs关联到各自的拉格朗日
 * DoFs所在的相同块中，因为它们相互耦合，而空间上相距遥远的
 * DoFs却不耦合。
 * 与局部系统的反演和应用预处理程序相比，对非拉格朗日DoF进行排序的额外计算量并不大，所以如果你想按分量对自由度进行排序，这种策略可能是合理的。如果自由度不按分量排序，上述两种策略的结果没有太大差别。然而，与第一种策略不同的是，如果自由度按分量重新编号，第二种策略的性能不会变差。  </ul>
 *
 *  <h3>Typical results</h3>
 * 作为一个典型的测试案例，我们使用了一个来自优化的非线性问题，这导致了一系列的鞍点问题，每个问题都使用GMRES和Vanka作为预处理程序来解决。该方程有大约850个自由度。使用非阻塞版本
 * @p SparseVanka （或 @p
 * SparseBlockVanka与<tt>n_blocks==1</tt>），在每个非线性步骤中，需要以下迭代次数来解决线性系统。
 *
 * @verbatim
 * 101 68 64 53 35 21
 * @endverbatim
 *
 * 如果有四个块，我们需要以下的迭代次数
 *
 * @verbatim
 * 124 88 83 66 44 28
 * @endverbatim
 * 可以看出，需要更多的迭代次数。然而，就计算时间而言，第一个版本需要72秒的上墙时间（和79秒的CPU时间，这比上墙时间多，因为程序的其他一些部分也被并行化了），而第二个版本在四处理器的机器上需要53秒的上墙时间（和110秒的CPU时间）。在这两种情况下，总时间都是由线性求解器主导的。在这种情况下，如果壁挂时间比CPU时间更重要，那么使用阻塞版本的预处理程序是值得的。
 * 上述阻断版的结果是用第一种阻断策略得到的，自由度没有按组件编号。使用第二种策略对迭代次数没有太大变化（每一步最多增加一次），而且当自由度按分量排序时也没有变化，而第一种策略则明显恶化了。
 *
 *
 */
template <typename number>
class SparseBlockVanka : public SparseVanka<number>
{
public:
  /**
   * 申报容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 枚举不同的方法，通过这些方法将自由度分配到我们要工作的区块上。
   *
   */
  enum BlockingStrategy
  {
    /**
     * 按索引间隔的区块。
     *
     */
    index_intervals,
    /**
     * 用自适应策略的区块。
     *
     */
    adaptive
  };

  /**
   * 构造函数。将除 @p n_blocks 以外的所有参数传递给基类。
   * @deprecated
   * 这个构造函数已被废弃。传递给最后两个参数的值被忽略了。
   *
   */
  DEAL_II_DEPRECATED
  SparseBlockVanka(const SparseMatrix<number> &M,
                   const std::vector<bool> &   selected,
                   const unsigned int          n_blocks,
                   const BlockingStrategy      blocking_strategy,
                   const bool                  conserve_memory,
                   const unsigned int n_threads = MultithreadInfo::n_threads());

  /**
   * 构造函数。将除 @p n_blocks 以外的所有参数传递给基类。
   *
   */
  SparseBlockVanka(const SparseMatrix<number> &M,
                   const std::vector<bool> &   selected,
                   const unsigned int          n_blocks,
                   const BlockingStrategy      blocking_strategy);

  /**
   * 应用预处理程序。
   *
   */
  template <typename number2>
  void
  vmult(Vector<number2> &dst, const Vector<number2> &src) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 存储块的数量。
   *
   */
  const unsigned int n_blocks;

  /**
   * 在这个字段中，我们为每个块预先计算哪些自由度属于它。因此，如果<tt>dof_masks[i][j]==true</tt>，那么DoF
   * @p j 就属于块 @p i.
   * 当然，没有其他<tt>dof_masks[l][j]</tt>可以为<tt>l！=i</tt>的
   * @p true
   * 。这个计算是在构造函数中完成的，以避免每次调用预处理程序时重新计算。
   *
   */
  std::vector<std::vector<bool>> dof_masks;

  /**
   * 计算字段的内容  @p dof_masks.
   * 这个函数是在构造函数中调用的。
   *
   */
  void
  compute_dof_masks(const SparseMatrix<number> &M,
                    const std::vector<bool> &   selected,
                    const BlockingStrategy      blocking_strategy);
};

 /*@}*/ 
 /* ---------------------------------- Inline functions ------------------- */ 

#ifndef DOXYGEN

template <typename number>
inline typename SparseVanka<number>::size_type
SparseVanka<number>::m() const
{
  Assert(_m != 0, ExcNotInitialized());
  return _m;
}

template <typename number>
inline typename SparseVanka<number>::size_type
SparseVanka<number>::n() const
{
  Assert(_n != 0, ExcNotInitialized());
  return _n;
}

template <typename number>
template <typename number2>
inline void
SparseVanka<number>::Tvmult(Vector<number2> &  /*dst*/ ,
                            const Vector<number2> &  /*src*/ ) const
{
  AssertThrow(false, ExcNotImplemented());
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


