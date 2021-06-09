//include/deal.II-translator/lac/matrix_iterator_0.txt
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

#ifndef dealii_matrix_iterator_h
#define dealii_matrix_iterator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 用于常数和非常数矩阵的迭代器。
 * 这个迭代器是从实际的矩阵类型中抽象出来的，可以用于任何具有所需ACCESSOR类型的矩阵。
 *
 *
 */
template <class ACCESSOR>
class MatrixIterator
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 为我们要操作的矩阵类型（包括常数）提供类型定义。
   *
   */
  using MatrixType = typename ACCESSOR::MatrixType;

  /**
   * 构造函数。为给定的<tt>行</tt>和其中的<tt>索引</tt>创建一个进入矩阵<tt>matrix</tt>的迭代器。
   *
   */
  MatrixIterator(MatrixType *    matrix,
                 const size_type row   = 0,
                 const size_type index = 0);

  /**
   * 从另一个矩阵迭代器复制。主要是为了允许从一个非常量的迭代器初始化，这个函数只要求从另一个迭代器的访问器转换到这个访问器对象就可以了。
   *
   */
  template <class OtherAccessor>
  MatrixIterator(const MatrixIterator<OtherAccessor> &other);

  /**
   * 前缀增量。
   *
   */
  MatrixIterator &
  operator++();

  /**
   * 后缀增量。
   *
   */
  MatrixIterator
  operator++(int);

  /**
   * 撤消运算符。
   *
   */
  const ACCESSOR &operator*() const;

  /**
   * 解除引用操作符。
   *
   */
  const ACCESSOR *operator->() const;

  /**
   * 比较。真，如果两个访问器都相等。
   *
   */
  bool
  operator==(const MatrixIterator &) const;

  /**
   * <tt>==</tt>的倒数。
   *
   */
  bool
  operator!=(const MatrixIterator &) const;

  /**
   * 比较运算符。如果第一个行号较小，或者行号相等且第一个索引较小，则结果为真。
   * 这个函数只有在两个迭代器都指向同一个矩阵时才有效。
   *
   */
  bool
  operator<(const MatrixIterator &) const;

  /**
   * 比较运算符。与上述运算符的工作方式相同，只是反过来了。
   *
   */
  bool
  operator>(const MatrixIterator &) const;

private:
  /**
   * 存储一个访问器类的对象。
   *
   */
  ACCESSOR accessor;

  // Allow other iterators access to private data.
  template <class OtherAccessor>
  friend class MatrixIterator;
};


//----------------------------------------------------------------------//

template <class ACCESSOR>
inline MatrixIterator<ACCESSOR>::MatrixIterator(MatrixType *    matrix,
                                                const size_type r,
                                                const size_type i)
  : accessor(matrix, r, i)
{}


template <class ACCESSOR>
template <class OtherAccessor>
inline MatrixIterator<ACCESSOR>::MatrixIterator(
  const MatrixIterator<OtherAccessor> &other)
  : accessor(other.accessor)
{}


template <class ACCESSOR>
inline MatrixIterator<ACCESSOR> &
MatrixIterator<ACCESSOR>::operator++()
{
  accessor.advance();
  return *this;
}


template <class ACCESSOR>
inline MatrixIterator<ACCESSOR>
MatrixIterator<ACCESSOR>::operator++(int)
{
  const MatrixIterator iter = *this;
  accessor.advance();
  return iter;
}


template <class ACCESSOR>
inline const ACCESSOR &MatrixIterator<ACCESSOR>::operator*() const
{
  return accessor;
}


template <class ACCESSOR>
inline const ACCESSOR *MatrixIterator<ACCESSOR>::operator->() const
{
  return &accessor;
}


template <class ACCESSOR>
inline bool
MatrixIterator<ACCESSOR>::operator==(const MatrixIterator &other) const
{
  return (accessor == other.accessor);
}


template <class ACCESSOR>
inline bool
MatrixIterator<ACCESSOR>::operator!=(const MatrixIterator &other) const
{
  return !(*this == other);
}


template <class ACCESSOR>
inline bool
MatrixIterator<ACCESSOR>::operator<(const MatrixIterator &other) const
{
  Assert(&accessor.get_matrix() == &other.accessor.get_matrix(),
         ExcInternalError());

  return (accessor < other.accessor);
}


template <class ACCESSOR>
inline bool
MatrixIterator<ACCESSOR>::operator>(const MatrixIterator &other) const
{
  return (other < *this);
}

DEAL_II_NAMESPACE_CLOSE

#endif


