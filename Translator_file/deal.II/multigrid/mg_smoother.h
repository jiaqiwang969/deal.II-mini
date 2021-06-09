//include/deal.II-translator/multigrid/mg_smoother_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_mg_smoother_h
#define dealii_mg_smoother_h


#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/* MGSmootherBase在mg_base.h中定义。

* 
*
*/

 /*!@addtogroup mg */ 
 /*@{*/ 

/**
 * 一个用于平滑器处理平滑信息的基类。虽然没有添加到MGSmootherBase中的抽象接口，但这个类存储了平滑步骤的数量和类型的信息，反过来可以被派生类使用。
 *
 *
 */
template <typename VectorType>
class MGSmoother : public MGSmootherBase<VectorType>
{
public:
  /**
   * 构造函数。
   *
   */
  MGSmoother(const unsigned int steps     = 1,
             const bool         variable  = false,
             const bool         symmetric = false,
             const bool         transpose = false);

  /**
   * 在最好的级别上修改平滑步骤的数量。
   *
   */
  void
  set_steps(const unsigned int);

  /**
   * 开启/关闭变量平滑。
   *
   */
  void
  set_variable(const bool);

  /**
   * 开启/关闭对称平滑。
   *
   */
  void
  set_symmetric(const bool);

  /**
   * 开启/关闭换位平滑。该效果被set_symmetric()所覆盖。
   *
   */
  void
  set_transpose(const bool);

  /**
   * 将 @p debug 设置为非零值以获得记录在 @p
   * deallog的调试信息。增加以获得更多信息
   *
   */
  void
  set_debug(const unsigned int level);

protected:
  /**
   * 一个用于临时向量的内存对象。
   * 该对象被标记为可变的，因为我们将需要用它来分配临时向量，也是在常量的函数中。
   *
   */
  mutable GrowingVectorMemory<VectorType> vector_memory;

  /**
   * 最细级别上的平滑步骤数。如果没有选择#variable
   * smoothing，这就是所有级别的步数。
   *
   */
  unsigned int steps;

  /**
   * 可变平滑：每当进入下一个更粗的级别时，平滑步骤的数量将增加一倍。
   *
   */
  bool variable;

  /**
   * 对称平滑：在平滑迭代中，在松弛方法和其转置之间交替进行。
   *
   */
  bool symmetric;

  /**
   * 使用松弛方法的转置而不是该方法本身。
   * 如果选择了#symmetric smoothing，这就没有影响。
   *
   */
  bool transpose;

  /**
   * 如果这不是零，则将调试信息输出到 @p deallog 。
   *
   */
  unsigned int debug;
};


/**
 * 平滑化不做任何事情。这个类对于很多应用来说是没有用的，除了测试一些多网格程序。另外，有些应用可能在没有平滑的情况下得到收敛，然后这个类会给你带来最便宜的多网格。
 *
 *
 */
template <typename VectorType>
class MGSmootherIdentity : public MGSmootherBase<VectorType>
{
public:
  /**
   * @p Multigrid. 接口的实现
   * 这个函数什么都不做，通过与这个函数的定义比较，这意味着平滑算子等于空算子。
   *
   */
  virtual void
  smooth(const unsigned int level, VectorType &u, const VectorType &rhs) const;

  virtual void
  clear();
};


namespace mg
{
  /**
   * 使用松弛类进行平滑处理。    一个松弛类是一个满足 @ref ConceptRelaxationType "松弛概念 "
   * 的对象。
   * 这个类在每个层次上执行平滑操作。该操作可以由几个参数控制。首先，松弛参数
   * @p  omega用于底层松弛方法。  @p steps
   * 是在最细级别上的松弛步骤的数量（如果 @p variable
   * 关闭，则在所有级别上）。如果 @p variable 是 @p true,
   * ，则每一个较粗的层次上的平滑步数是双倍的。这导致了一种具有W-循环复杂性的方法，但节省了网格转移。这是Bramble等人提出的方法。选项
   * @p symmetric
   * 按照Bramble的建议，在每个步骤中交替使用平滑器和其转置。
   * @p transpose
   * 使用<tt>Tstep</tt>的转置平滑操作，而不是放松方案的常规<tt>step</tt>。
   * 如果你使用块矩阵，第二个 @p initialize
   * 函数提供了提取单个块进行平滑的可能性。在这种情况下，多棱镜方法必须只用于与该单一块相关的矢量。
   *
   */
  template <class RelaxationType, typename VectorType>
  class SmootherRelaxation : public MGLevelObject<RelaxationType>,
                             public MGSmoother<VectorType>
  {
  public:
    /**
     * 构造函数。设置平滑参数。
     *
     */
    SmootherRelaxation(const unsigned int steps     = 1,
                       const bool         variable  = false,
                       const bool         symmetric = false,
                       const bool         transpose = false);

    /**
     * 对矩阵进行初始化。这个函数用每个级别的相同的平滑器来初始化平滑运算器。
     * @p additional_data 是一个 @p  RelaxationType::AdditionalData
     * 类型的对象，被交给松弛方法的初始化函数。
     *
     */
    template <typename MatrixType2>
    void
    initialize(const MGLevelObject<MatrixType2> &             matrices,
               const typename RelaxationType::AdditionalData &additional_data =
                 typename RelaxationType::AdditionalData());

    /**
     * 初始化每个层次的矩阵和附加数据。
     * 如果两个对象的最小或最大级别不同，则利用最大的共同范围。这样，即使矩阵是为所有级别生成的，平滑也可以限制在某些级别上。
     *
     */
    template <typename MatrixType2, class DATA>
    void
    initialize(const MGLevelObject<MatrixType2> &matrices,
               const MGLevelObject<DATA> &       additional_data);

    /**
     * 清空所有的向量。
     *
     */
    void
    clear() override;

    /**
     * 实际的平滑方法。
     *
     */
    virtual void
    smooth(const unsigned int level,
           VectorType &       u,
           const VectorType & rhs) const override;

    /**
     * 平滑化的应用变体，在调用平滑函数之前将向量u设置为零。这个函数等同于以下代码
     * @code
     * u = 0;
     * smooth(level, u, rhs);
     * @endcode
     * 在多网格预处理接口中，apply()方法用于预平滑操作，因为解向量中以前的内容需要为新进来的残差所覆盖。另一方面，所有后续的操作都需要平滑已经存在于向量
     * @p u 中的内容，给定的右手边，这是由smooth()完成的。
     *
     */
    virtual void
    apply(const unsigned int level,
          VectorType &       u,
          const VectorType & rhs) const override;

    /**
     * 这个对象使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;
  };
} // namespace mg

/**
 * 使用满足 @ref ConceptRelaxationType "放松概念 "
 * 的求解器进行平滑处理。
 * 该类在每个层次上执行平滑操作。该操作可以由几个参数控制。首先，松弛参数
 * @p omega 被用于底层松弛方法中。  @p steps
 * 是最细层次上的松弛步数（如果 @p variable
 * 关闭，则在所有层次上）。如果 @p variable 是 @p true,
 * ，则每一个较粗的层次上的平滑步数是双倍的。这导致了一种具有W-循环复杂性的方法，但节省了网格转移。这是Bramble
 * at al.提出的方法。 选项 @p symmetric
 * 按照Bramble的建议，在每个步骤中交替使用平滑器和其转置。
 * @p transpose
 * 使用<tt>Tstep</tt>的转置平滑操作，而不是放松方案的常规<tt>step</tt>。
 * 如果你使用块矩阵，第二个 @p initialize
 * 函数提供了提取单个块进行平滑的可能性。在这种情况下，多棱镜方法必须只用于与该单块相关的向量。
 * 该库包含<tt>SparseMatrix<.></tt>和<tt>Vector<.></tt>的实例化，其中模板参数是
 * @p  float和 @p double.
 * 的所有组合，其他实例化可以通过包括文件mg_smoother.templates.h来创建。
 *
 *
 */
template <typename MatrixType, class RelaxationType, typename VectorType>
class MGSmootherRelaxation : public MGSmoother<VectorType>
{
public:
  /**
   * 构造函数。设置平滑参数。
   *
   */
  MGSmootherRelaxation(const unsigned int steps     = 1,
                       const bool         variable  = false,
                       const bool         symmetric = false,
                       const bool         transpose = false);

  /**
   * 对矩阵进行初始化。这个函数存储指向水平矩阵的指针，并且用每个水平的相同的平滑器来初始化平滑运算器。
   * @p additional_data 是一个 @p RelaxationType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2>
  void
  initialize(const MGLevelObject<MatrixType2> &             matrices,
             const typename RelaxationType::AdditionalData &additional_data =
               typename RelaxationType::AdditionalData());

  /**
   * 为矩阵进行初始化。这个函数存储指向水平矩阵的指针，并且用每个水平的相应平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p RelaxationType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const MGLevelObject<DATA> &       additional_data);

  /**
   * 对单块矩阵进行初始化。在这个块状矩阵中，每一级都会选择
   * @p block_row 和 @p block_col 所指示的块。
   * 这个函数存储了指向水平矩阵的指针，并以每个水平的相同平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p RelaxationType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const DATA &                      additional_data,
             const unsigned int                block_row,
             const unsigned int                block_col);

  /**
   * 对单块矩阵进行初始化。在这个块状矩阵中，每一级都会选择
   * @p block_row 和 @p block_col 所指示的块。
   * 这个函数存储指向水平矩阵的指针，并为每个水平用相应的平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p RelaxationType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const MGLevelObject<DATA> &       additional_data,
             const unsigned int                block_row,
             const unsigned int                block_col);

  /**
   * 清空所有的向量。
   *
   */
  void
  clear();

  /**
   * 实际的平滑方法。
   *
   */
  virtual void
  smooth(const unsigned int level, VectorType &u, const VectorType &rhs) const;

  /**
   * 平滑化的应用变体，在调用平滑函数之前将向量u设置为零。这个函数等同于以下代码
   * @code
   * u = 0;
   * smooth(level, u, rhs);
   * @endcode
   * 在多网格预处理接口中，apply()方法用于预平滑操作，因为解向量中以前的内容需要为新进来的残差所覆盖。另一方面，所有后续的操作都需要平滑已经存在于向量
   * @p u 中的内容，给定的右手边，这是由smooth()完成。
   *
   */
  virtual void
  apply(const unsigned int level, VectorType &u, const VectorType &rhs) const;

  /**
   * 包含放松方法的对象。
   *
   */
  MGLevelObject<RelaxationType> smoothers;

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;


private:
  /**
   * 指向矩阵的指针。
   *
   */
  MGLevelObject<LinearOperator<VectorType>> matrices;
};



/**
 * 使用预处理器类的平滑器。
 * 该类在每一级上执行平滑操作。该操作可以由几个参数控制。首先，放松参数
 * @p omega 用于底层放松方法。  @p steps
 * 是在最细级别上的松弛步骤的数量（如果 @p variable
 * 关闭，则在所有级别上）。如果 @p variable 是 @p true,
 * ，则每一个较粗的层次上的平滑步数是双倍的。这导致了一种具有W-循环复杂性的方法，但节省了网格转移。这是Bramble
 * at al.提出的方法。 选项 @p symmetric
 * 按照Bramble的建议，在每个步骤中交替使用平滑器和其转置。
 * @p transpose
 * 使用<tt>Tvmult</tt>的转置平滑操作，而不是放松方案的常规<tt>vmult</tt>。
 * 如果你使用块矩阵，第二个 @p initialize
 * 函数提供了提取单个块进行平滑的可能性。在这种情况下，多棱镜方法必须只用于与该单块相关的向量。
 * 该库包含<tt>SparseMatrix<.></tt>和<tt>Vector<.></tt>的实例化，其中模板参数是
 * @p  float和 @p double.
 * 的所有组合，其他实例化可以通过包括文件mg_smoother.templates.h来创建。
 *
 *
 */
template <typename MatrixType, typename PreconditionerType, typename VectorType>
class MGSmootherPrecondition : public MGSmoother<VectorType>
{
public:
  /**
   * 构造函数。设置平滑参数。
   *
   */
  MGSmootherPrecondition(const unsigned int steps     = 1,
                         const bool         variable  = false,
                         const bool         symmetric = false,
                         const bool         transpose = false);

  /**
   * 对矩阵进行初始化。这个函数存储指向水平矩阵的指针，并且用每个水平的相同的平滑器来初始化平滑运算器。
   * @p additional_data 是一个 @p  PreconditionerType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const typename PreconditionerType::AdditionalData &
               additional_data = typename PreconditionerType::AdditionalData());

  /**
   * 为矩阵进行初始化。这个函数存储指向水平矩阵的指针，并且用每个水平的相应平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p  PreconditionerType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const MGLevelObject<DATA> &       additional_data);

  /**
   * 对单块矩阵进行初始化。在这个块状矩阵中， @p block_row
   * 和 @p block_col 所指示的块在每一级都被选择。
   * 这个函数存储了指向水平矩阵的指针，并以每个水平的相同平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p  PreconditionerType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const DATA &                      additional_data,
             const unsigned int                block_row,
             const unsigned int                block_col);

  /**
   * 对单块矩阵进行初始化。在这个块状矩阵中， @p block_row
   * 和 @p block_col 所指示的块在每一级都被选择。
   * 这个函数存储指向水平矩阵的指针，并为每个水平用相应的平滑器初始化平滑运算器。
   * @p additional_data 是一个 @p  PreconditionerType::AdditionalData
   * 类型的对象，被交给放松方法的初始化函数。
   *
   */
  template <typename MatrixType2, class DATA>
  void
  initialize(const MGLevelObject<MatrixType2> &matrices,
             const MGLevelObject<DATA> &       additional_data,
             const unsigned int                block_row,
             const unsigned int                block_col);

  /**
   * 清空所有的向量。
   *
   */
  void
  clear() override;

  /**
   * 实际的平滑方法。
   *
   */
  virtual void
  smooth(const unsigned int level,
         VectorType &       u,
         const VectorType & rhs) const override;

  /**
   * 平滑化的应用变体，在调用平滑函数之前将向量u设置为零。这个函数等同于以下代码
   * @code
   * u = 0;
   * smooth(level, u, rhs);
   * @endcode
   * 在多网格预处理接口中，apply()方法用于预平滑操作，因为解向量中以前的内容需要为新进来的残差所覆盖。另一方面，所有后续的操作都需要平滑已经存在于向量
   * @p u 中的内容，给定的右手边，这是由smooth()完成。
   *
   */
  virtual void
  apply(const unsigned int level,
        VectorType &       u,
        const VectorType & rhs) const override;

  /**
   * 包含放松方法的对象。
   *
   */
  MGLevelObject<PreconditionerType> smoothers;

  /**
   * 这个对象所使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;


private:
  /**
   * 指向矩阵的指针。
   *
   */
  MGLevelObject<LinearOperator<VectorType>> matrices;
};

 /*@}*/ 

 /* ------------------------------- Inline functions --------------------------
 */ 

#ifndef DOXYGEN

template <typename VectorType>
inline void
MGSmootherIdentity<VectorType>::smooth(const unsigned int,
                                       VectorType &,
                                       const VectorType &) const
{}

template <typename VectorType>
inline void
MGSmootherIdentity<VectorType>::clear()
{}

//---------------------------------------------------------------------------

template <typename VectorType>
inline MGSmoother<VectorType>::MGSmoother(const unsigned int steps,
                                          const bool         variable,
                                          const bool         symmetric,
                                          const bool         transpose)
  : steps(steps)
  , variable(variable)
  , symmetric(symmetric)
  , transpose(transpose)
  , debug(0)
{}


template <typename VectorType>
inline void
MGSmoother<VectorType>::set_steps(const unsigned int s)
{
  steps = s;
}


template <typename VectorType>
inline void
MGSmoother<VectorType>::set_debug(const unsigned int s)
{
  debug = s;
}


template <typename VectorType>
inline void
MGSmoother<VectorType>::set_variable(const bool flag)
{
  variable = flag;
}


template <typename VectorType>
inline void
MGSmoother<VectorType>::set_symmetric(const bool flag)
{
  symmetric = flag;
}


template <typename VectorType>
inline void
MGSmoother<VectorType>::set_transpose(const bool flag)
{
  transpose = flag;
}

//----------------------------------------------------------------------//

namespace mg
{
  template <class RelaxationType, typename VectorType>
  inline SmootherRelaxation<RelaxationType, VectorType>::SmootherRelaxation(
    const unsigned int steps,
    const bool         variable,
    const bool         symmetric,
    const bool         transpose)
    : MGSmoother<VectorType>(steps, variable, symmetric, transpose)
  {}


  template <class RelaxationType, typename VectorType>
  inline void
  SmootherRelaxation<RelaxationType, VectorType>::clear()
  {
    MGLevelObject<RelaxationType>::clear_elements();
  }


  template <class RelaxationType, typename VectorType>
  template <typename MatrixType2>
  inline void
  SmootherRelaxation<RelaxationType, VectorType>::initialize(
    const MGLevelObject<MatrixType2> &             m,
    const typename RelaxationType::AdditionalData &data)
  {
    const unsigned int min = m.min_level();
    const unsigned int max = m.max_level();

    this->resize(min, max);

    for (unsigned int i = min; i <= max; ++i)
      (*this)[i].initialize(m[i], data);
  }


  template <class RelaxationType, typename VectorType>
  template <typename MatrixType2, class DATA>
  inline void
  SmootherRelaxation<RelaxationType, VectorType>::initialize(
    const MGLevelObject<MatrixType2> &m,
    const MGLevelObject<DATA> &       data)
  {
    const unsigned int min = std::max(m.min_level(), data.min_level());
    const unsigned int max = std::min(m.max_level(), data.max_level());

    this->resize(min, max);

    for (unsigned int i = min; i <= max; ++i)
      (*this)[i].initialize(m[i], data[i]);
  }


  template <class RelaxationType, typename VectorType>
  inline void
  SmootherRelaxation<RelaxationType, VectorType>::smooth(
    const unsigned int level,
    VectorType &       u,
    const VectorType & rhs) const
  {
    unsigned int maxlevel = this->max_level();
    unsigned int steps2   = this->steps;

    if (this->variable)
      steps2 *= (1 << (maxlevel - level));

    bool T = this->transpose;
    if (this->symmetric && (steps2 % 2 == 0))
      T = false;
    if (this->debug > 0)
      deallog << 'S' << level << ' ';

    for (unsigned int i = 0; i < steps2; ++i)
      {
        if (T)
          (*this)[level].Tstep(u, rhs);
        else
          (*this)[level].step(u, rhs);
        if (this->symmetric)
          T = !T;
      }
  }


  template <class RelaxationType, typename VectorType>
  inline void
  SmootherRelaxation<RelaxationType, VectorType>::apply(
    const unsigned int level,
    VectorType &       u,
    const VectorType & rhs) const
  {
    unsigned int maxlevel = this->max_level();
    unsigned int steps2   = this->steps;

    if (this->variable)
      steps2 *= (1 << (maxlevel - level));

    bool T = this->transpose;
    if (this->symmetric && (steps2 % 2 == 0))
      T = false;
    if (this->debug > 0)
      deallog << 'S' << level << ' ';

    if (T)
      (*this)[level].Tvmult(u, rhs);
    else
      (*this)[level].vmult(u, rhs);
    if (this->symmetric)
      T = !T;
    for (unsigned int i = 1; i < steps2; ++i)
      {
        if (T)
          (*this)[level].Tstep(u, rhs);
        else
          (*this)[level].step(u, rhs);
        if (this->symmetric)
          T = !T;
      }
  }


  template <class RelaxationType, typename VectorType>
  inline std::size_t
  SmootherRelaxation<RelaxationType, VectorType>::memory_consumption() const
  {
    return sizeof(*this) - sizeof(MGLevelObject<RelaxationType>) +
           MGLevelObject<RelaxationType>::memory_consumption() +
           this->vector_memory.memory_consumption();
  }
} // namespace mg


//----------------------------------------------------------------------//

template <typename MatrixType, class RelaxationType, typename VectorType>
inline MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::
  MGSmootherRelaxation(const unsigned int steps,
                       const bool         variable,
                       const bool         symmetric,
                       const bool         transpose)
  : MGSmoother<VectorType>(steps, variable, symmetric, transpose)
{}



template <typename MatrixType, class RelaxationType, typename VectorType>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::clear()
{
  smoothers.clear_elements();

  unsigned int i = matrices.min_level(), max_level = matrices.max_level();
  for (; i <= max_level; ++i)
    matrices[i] = LinearOperator<VectorType>();
}


template <typename MatrixType, class RelaxationType, typename VectorType>
template <typename MatrixType2>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &             m,
  const typename RelaxationType::AdditionalData &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] =
        linear_operator<VectorType>(LinearOperator<VectorType>(), m[i]);
      smoothers[i].initialize(m[i], data);
    }
}

template <typename MatrixType, class RelaxationType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const MGLevelObject<DATA> &       data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert(data.min_level() == min, ExcDimensionMismatch(data.min_level(), min));
  Assert(data.max_level() == max, ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] =
        linear_operator<VectorType>(LinearOperator<VectorType>(), m[i]);
      smoothers[i].initialize(m[i], data[i]);
    }
}

template <typename MatrixType, class RelaxationType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const DATA &                      data,
  const unsigned int                row,
  const unsigned int                col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] = linear_operator<VectorType>(LinearOperator<VectorType>(),
                                                m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}

template <typename MatrixType, class RelaxationType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const MGLevelObject<DATA> &       data,
  const unsigned int                row,
  const unsigned int                col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert(data.min_level() == min, ExcDimensionMismatch(data.min_level(), min));
  Assert(data.max_level() == max, ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] = linear_operator<VectorType>(LinearOperator<VectorType>(),
                                                m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}


template <typename MatrixType, class RelaxationType, typename VectorType>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::smooth(
  const unsigned int level,
  VectorType &       u,
  const VectorType & rhs) const
{
  unsigned int maxlevel = smoothers.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  for (unsigned int i = 0; i < steps2; ++i)
    {
      if (T)
        smoothers[level].Tstep(u, rhs);
      else
        smoothers[level].step(u, rhs);
      if (this->symmetric)
        T = !T;
    }
}


template <typename MatrixType, class RelaxationType, typename VectorType>
inline void
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::apply(
  const unsigned int level,
  VectorType &       u,
  const VectorType & rhs) const
{
  unsigned int maxlevel = smoothers.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  if (T)
    smoothers[level].Tvmult(u, rhs);
  else
    smoothers[level].vmult(u, rhs);
  if (this->symmetric)
    T = !T;
  for (unsigned int i = 1; i < steps2; ++i)
    {
      if (T)
        smoothers[level].Tstep(u, rhs);
      else
        smoothers[level].step(u, rhs);
      if (this->symmetric)
        T = !T;
    }
}



template <typename MatrixType, class RelaxationType, typename VectorType>
inline std::size_t
MGSmootherRelaxation<MatrixType, RelaxationType, VectorType>::
  memory_consumption() const
{
  return sizeof(*this) + matrices.memory_consumption() +
         smoothers.memory_consumption() +
         this->vector_memory.memory_consumption();
}


//----------------------------------------------------------------------//

template <typename MatrixType, typename PreconditionerType, typename VectorType>
inline MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::
  MGSmootherPrecondition(const unsigned int steps,
                         const bool         variable,
                         const bool         symmetric,
                         const bool         transpose)
  : MGSmoother<VectorType>(steps, variable, symmetric, transpose)
{}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::clear()
{
  smoothers.clear_elements();

  unsigned int i = matrices.min_level(), max_level = matrices.max_level();
  for (; i <= max_level; ++i)
    matrices[i] = LinearOperator<VectorType>();
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
template <typename MatrixType2>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &                 m,
  const typename PreconditionerType::AdditionalData &data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] =
        linear_operator<VectorType>(LinearOperator<VectorType>(), m[i]);
      smoothers[i].initialize(m[i], data);
    }
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const MGLevelObject<DATA> &       data)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert(data.min_level() == min, ExcDimensionMismatch(data.min_level(), min));
  Assert(data.max_level() == max, ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a rich
      // enough interface to populate reinit_(domain|range)_vector. Thus,
      // apply an empty LinearOperator exemplar.
      matrices[i] =
        linear_operator<VectorType>(LinearOperator<VectorType>(),
                                    Utilities::get_underlying_value(m[i]));
      smoothers[i].initialize(Utilities::get_underlying_value(m[i]), data[i]);
    }
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const DATA &                      data,
  const unsigned int                row,
  const unsigned int                col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data);
    }
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
template <typename MatrixType2, class DATA>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::initialize(
  const MGLevelObject<MatrixType2> &m,
  const MGLevelObject<DATA> &       data,
  const unsigned int                row,
  const unsigned int                col)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  Assert(data.min_level() == min, ExcDimensionMismatch(data.min_level(), min));
  Assert(data.max_level() == max, ExcDimensionMismatch(data.max_level(), max));

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      matrices[i] = &(m[i].block(row, col));
      smoothers[i].initialize(m[i].block(row, col), data[i]);
    }
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::smooth(
  const unsigned int level,
  VectorType &       u,
  const VectorType & rhs) const
{
  unsigned int maxlevel = matrices.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  typename VectorMemory<VectorType>::Pointer r(this->vector_memory);
  typename VectorMemory<VectorType>::Pointer d(this->vector_memory);

  r->reinit(u, true);
  d->reinit(u, true);

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  for (unsigned int i = 0; i < steps2; ++i)
    {
      if (T)
        {
          if (this->debug > 0)
            deallog << 'T';
          matrices[level].Tvmult(*r, u);
          r->sadd(-1., 1., rhs);
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].Tvmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      else
        {
          if (this->debug > 0)
            deallog << 'N';
          matrices[level].vmult(*r, u);
          r->sadd(-1., rhs);
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].vmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      u += *d;
      if (this->symmetric)
        T = !T;
    }
  if (this->debug > 0)
    deallog << std::endl;
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
inline void
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::apply(
  const unsigned int level,
  VectorType &       u,
  const VectorType & rhs) const
{
  unsigned int maxlevel = matrices.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;
  if (this->debug > 0)
    deallog << 'S' << level << ' ';

  // first step where we overwrite the result
  if (this->debug > 2)
    deallog << ' ' << rhs.l2_norm() << ' ';
  if (this->debug > 0)
    deallog << (T ? 'T' : 'N');
  if (T)
    smoothers[level].Tvmult(u, rhs);
  else
    smoothers[level].vmult(u, rhs);
  if (this->debug > 1)
    deallog << ' ' << u.l2_norm() << ' ';
  if (this->symmetric)
    T = !T;

  typename VectorMemory<VectorType>::Pointer r(this->vector_memory);
  typename VectorMemory<VectorType>::Pointer d(this->vector_memory);

  if (steps2 > 1)
    {
      r->reinit(u, true);
      d->reinit(u, true);
    }

  for (unsigned int i = 1; i < steps2; ++i)
    {
      if (T)
        {
          if (this->debug > 0)
            deallog << 'T';
          matrices[level].Tvmult(*r, u);
          r->sadd(-1., 1., rhs);
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].Tvmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      else
        {
          if (this->debug > 0)
            deallog << 'N';
          matrices[level].vmult(*r, u);
          r->sadd(-1., rhs);
          if (this->debug > 2)
            deallog << ' ' << r->l2_norm() << ' ';
          smoothers[level].vmult(*d, *r);
          if (this->debug > 1)
            deallog << ' ' << d->l2_norm() << ' ';
        }
      u += *d;
      if (this->symmetric)
        T = !T;
    }
  if (this->debug > 0)
    deallog << std::endl;
}



template <typename MatrixType, typename PreconditionerType, typename VectorType>
inline std::size_t
MGSmootherPrecondition<MatrixType, PreconditionerType, VectorType>::
  memory_consumption() const
{
  return sizeof(*this) + matrices.memory_consumption() +
         smoothers.memory_consumption() +
         this->vector_memory.memory_consumption();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


