//include/deal.II-translator/multigrid/mg_block_smoother_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_mg_block_smoother_h
#define dealii_mg_block_smoother_h


#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_smoother.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/* MGSmootherBase在mg_base.h中定义。

* 
*
*/

 /*!@addtogroup mg */ 
 /*@{*/ 

/**
 * 用于块状向量的通用平滑器类。该类通过初始化一个矩阵和一个平滑器对象，给了选择块状平滑器的完全自由。因此，每一级的平滑器对象必须由手工构建。
 *
 *
 */
template <typename MatrixType, class RelaxationType, typename number>
class MGSmootherBlock : public MGSmoother<BlockVector<number>>
{
public:
  /**
   * 构造器。
   *
   */
  MGSmootherBlock(const unsigned int steps     = 1,
                  const bool         variable  = false,
                  const bool         symmetric = false,
                  const bool         transpose = false,
                  const bool         reverse   = false);

  /**
   * 对矩阵进行初始化。参数<tt>matrices</tt>可以是任何具有<tt>get_minlevel()</tt>和<tt>get_maxlevel()</tt>函数以及<tt>operator[]</tt>的对象，返回对
   * @p MatrixType.
   * 的引用。参数<tt>smoothers</tt>使用同样的约定，这样<tt>operator[]</tt>返回在单一层次上做块平滑的对象。
   * 这个函数存储了指向每个级别的级别矩阵和平滑运算器的指针。
   *
   */
  template <class MGMatrixType, class MGRelaxationType>
  void
  initialize(const MGMatrixType &matrices, const MGRelaxationType &smoothers);

  /**
   * 清空所有的向量。
   *
   */
  void
  clear();

  /**
   * 开启/关闭反转。这与转置()是相互排斥的。
   *
   */
  void
  set_reverse(const bool);

  /**
   * 对 @p Multigrid. 接口的实现
   * 这个函数什么都不做，通过与这个函数的定义比较，这意味着平滑运算符等于空运算符。
   *
   */
  virtual void
  smooth(const unsigned int         level,
         BlockVector<number> &      u,
         const BlockVector<number> &rhs) const;

  /**
   * 这个对象使用的内存。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 指向矩阵的指针。
   *
   */
  MGLevelObject<LinearOperator<BlockVector<number>>> matrices;

  /**
   * 指向矩阵的指针。
   *
   */
  MGLevelObject<LinearOperator<BlockVector<number>>> smoothers;

  /**
   * 反转？
   *
   */
  bool reverse;

  /**
   * 用于辅助向量的内存。
   *
   */
  SmartPointer<VectorMemory<BlockVector<number>>,
               MGSmootherBlock<MatrixType, RelaxationType, number>>
    mem;
};

 /**@}*/ 

//---------------------------------------------------------------------------

#ifndef DOXYGEN

template <typename MatrixType, class RelaxationType, typename number>
inline MGSmootherBlock<MatrixType, RelaxationType, number>::MGSmootherBlock(
  const unsigned int steps,
  const bool         variable,
  const bool         symmetric,
  const bool         transpose,
  const bool         reverse)
  : MGSmoother<BlockVector<number>>(steps, variable, symmetric, transpose)
  , reverse(reverse)
  , mem(&this->vector_memory)
{}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::clear()
{
  unsigned int i = matrices.min_level(), max_level = matrices.max_level();
  for (; i <= max_level; ++i)
    {
      smoothers[i] = LinearOperator<BlockVector<number>>();
      matrices[i]  = LinearOperator<BlockVector<number>>();
    }
}


template <typename MatrixType, class RelaxationType, typename number>
template <class MGMatrixType, class MGRelaxationType>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::initialize(
  const MGMatrixType &    m,
  const MGRelaxationType &s)
{
  const unsigned int min = m.min_level();
  const unsigned int max = m.max_level();

  matrices.resize(min, max);
  smoothers.resize(min, max);

  for (unsigned int i = min; i <= max; ++i)
    {
      // Workaround: Unfortunately, not every "m[i]" object has a
      // rich enough interface to populate reinit_(domain|range)_vector.
      // Thus, apply an empty LinearOperator exemplar.
      matrices[i] = linear_operator<BlockVector<number>>(
        LinearOperator<BlockVector<number>>(), m[i]);
      smoothers[i] = linear_operator<BlockVector<number>>(matrices[i], s[i]);
    }
}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::set_reverse(
  const bool flag)
{
  reverse = flag;
}


template <typename MatrixType, class RelaxationType, typename number>
inline std::size_t
MGSmootherBlock<MatrixType, RelaxationType, number>::memory_consumption() const
{
  return sizeof(*this) + matrices.memory_consumption() +
         smoothers.memory_consumption() +
         this->vector_memory.memory_consumption();
}


template <typename MatrixType, class RelaxationType, typename number>
inline void
MGSmootherBlock<MatrixType, RelaxationType, number>::smooth(
  const unsigned int         level,
  BlockVector<number> &      u,
  const BlockVector<number> &rhs) const
{
  LogStream::Prefix prefix("Smooth");

  unsigned int maxlevel = matrices.max_level();
  unsigned int steps2   = this->steps;

  if (this->variable)
    steps2 *= (1 << (maxlevel - level));

  typename VectorMemory<BlockVector<number>>::Pointer r(*this->mem);
  typename VectorMemory<BlockVector<number>>::Pointer d(*this->mem);
  r->reinit(u);
  d->reinit(u);

  bool T = this->transpose;
  if (this->symmetric && (steps2 % 2 == 0))
    T = false;

  for (unsigned int i = 0; i < steps2; ++i)
    {
      if (T)
        {
          matrices[level].vmult(*r, u);
          r->sadd(-1., 1., rhs);
          smoothers[level].Tvmult(*d, *r);
        }
      else
        {
          matrices[level].vmult(*r, u);
          r->sadd(-1., 1., rhs);
          smoothers[level].vmult(*d, *r);
        }
      u += *d;
      if (this->symmetric)
        T = !T;
    }
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


