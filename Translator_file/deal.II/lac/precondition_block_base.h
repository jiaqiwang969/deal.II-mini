//include/deal.II-translator/lac/precondition_block_base_0.txt
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

#ifndef dealii_precondition_block_base_h
#define dealii_precondition_block_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/householder.h>
#include <deal.II/lac/lapack_full_matrix.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class FullMatrix;
template <typename number>
class Vector;
#endif

/**
 * 一个存储块状预处理程序和块状松弛方法的逆对角线块的类。
 * 这个类为基于反转矩阵对角线上的块的预处理和松弛方法做记账工作。它允许我们存储所有的对角线块和它们的倒数，或者为每个条目存储相同的块，并且它可以跟踪选择。因此，在初始化它并正确填充对角线上的逆块后，派生类可以使用inverse()，其整数参数指的是块数。
 * 此外，它还允许存储原始对角线块，而不仅仅是逆向的。例如，这些都是用于SSOR预处理程序的中间步骤。
 *
 *
 */
template <typename number>
class PreconditionBlockBase
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 选择一个反转块的方法，从而选择反转的数据类型。
   *
   */
  enum Inversion
  {
    /**
     * 使用 FullMatrix::inverse().
     * 中实现的标准高斯-雅各比方法。
     *
     */
    gauss_jordan,
    /**
     * 使用Householder类的QR分解。
     *
     */
    householder,
    /**
     * 使用LAPACKFullMatrix的奇异值分解。
     *
     */
    svd
  };

  /**
   * 构造函数初始化默认值。
   *
   */
  PreconditionBlockBase(bool      store_diagonals = false,
                        Inversion method          = gauss_jordan);

  /**
   * 虚拟解构器
   *
   */
  ~PreconditionBlockBase() = default;

  /**
   * 删除对角线块矩阵（如果存在的话），从而使该类处于调用构造函数后的直接状态。
   *
   */
  void
  clear();

  /**
   * 用给定的块大小调整到这个对角线块的数量。如果<tt>compress</tt>为真，那么将只存储一个块。
   *
   */
  void
  reinit(unsigned int nblocks,
         size_type    blocksize,
         bool         compress,
         Inversion    method = gauss_jordan);

  /**
   * 告诉类，计算反相。
   *
   */
  void
  inverses_computed(bool are_they);

  /**
   * 矩阵是否只使用一个对角线块？
   *
   */
  bool
  same_diagonal() const;

  /**
   * 检查，是否应该存储对角线块（而不是它们的逆向）。
   *
   */
  bool
  store_diagonals() const;

  /**
   * 如果反向块可以使用，则返回true。
   *
   */
  bool
  inverses_ready() const;

  /**
   * 块的数量。
   *
   */
  unsigned int
  size() const;

  /**
   * 与位于<tt>i</tt>位置的逆块相乘。
   *
   */
  template <typename number2>
  void
  inverse_vmult(size_type              i,
                Vector<number2> &      dst,
                const Vector<number2> &src) const;

  /**
   * 与位置<tt>i</tt>的换位逆块相乘。
   *
   */
  template <typename number2>
  void
  inverse_Tvmult(size_type              i,
                 Vector<number2> &      dst,
                 const Vector<number2> &src) const;

  /**
   * 如果反转是#gauss_jordan，则访问反转对角线块。
   *
   */
  FullMatrix<number> &
  inverse(size_type i);

  /**
   * 如果反演是#householder，则访问反对角线块。
   *
   */
  Householder<number> &
  inverse_householder(size_type i);

  /**
   * 如果反转是#householder，则访问反对角线块。
   *
   */
  LAPACKFullMatrix<number> &
  inverse_svd(size_type i);

  /**
   * 访问逆对角线区块。
   *
   */
  const FullMatrix<number> &
  inverse(size_type i) const;

  /**
   * 如果反转是#householder，则访问逆对角线区块。
   *
   */
  const Householder<number> &
  inverse_householder(size_type i) const;

  /**
   * 如果反转是#householder，则访问逆对角线区块。
   *
   */
  const LAPACKFullMatrix<number> &
  inverse_svd(size_type i) const;

  /**
   * 访问对角线块。
   *
   */
  FullMatrix<number> &
  diagonal(size_type i);

  /**
   * 访问对角线块。
   *
   */
  const FullMatrix<number> &
  diagonal(size_type i) const;

  /**
   * 向 @p deallog. 打印一些关于求逆者的统计数据
   * 输出取决于#求逆者。对于svd来说，它是最丰富的，在那里我们可以获得关于极值奇异值和条件数的统计数据。
   *
   */
  void
  log_statistics() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 你试图访问一个对角线块（而不是它的逆向），但你决定不存储对角线块。
   *
   */
  DeclException0(ExcDiagonalsNotStored);

  /**
   * 你正在访问一个对角线块，假设它具有某种类型。
   * 但是，用于反转对角线块的方法并没有使用这种类型
   *
   */
  DeclException0(ExcInverseNotAvailable);

protected:
  /**
   * 用于倒置块的方法。
   *
   */
  Inversion inversion;

private:
  /**
   * 对角线块的（反转）数量，如果只存储一个。
   *
   */
  unsigned int n_diagonal_blocks;

  /**
   * 如果使用反转#gauss_jordan，则以<tt>FullMatrix<number></tt>矩阵的形式存储对角线块矩阵的反转矩阵。
   * 与<tt>number=float</tt>相比，使用<tt>number=double</tt>可以节省内存，但可能会引入数值不稳定。
   *
   */
  std::vector<FullMatrix<number>> var_inverse_full;

  /**
   * 如果使用反转#householder，将对角块矩阵的反矩阵存储为<tt>Householder</tt>矩阵。与<tt>number=float</tt>相比，使用<tt>number=double</tt>可以节省内存，但可能会引入数值不稳定。
   *
   */
  std::vector<Householder<number>> var_inverse_householder;

  /**
   * 如果使用反转#svd，将对角块矩阵的反转矩阵存储为<tt>LAPACKFullMatrix</tt>矩阵。与<tt>number=float</tt>相比，使用<tt>number=double</tt>可以节省内存，但可能会引入数值不稳定。
   *
   */
  std::vector<LAPACKFullMatrix<number>> var_inverse_svd;

  /**
   * 存储原始对角线块。    被封锁的SSOR方法所使用。
   *
   */
  std::vector<FullMatrix<number>> var_diagonal;


  /**
   * 如果要使用#var_diagonal这个字段，则为真。
   *
   */
  bool var_store_diagonals;

  /**
   * 如果只存储一个逆数，则为真。
   *
   */
  bool var_same_diagonal;

  /**
   * 逆矩阵是可以使用的。由父类通过inverses_computed()设置。
   *
   */
  bool var_inverses_ready;
};

//----------------------------------------------------------------------//

template <typename number>
inline PreconditionBlockBase<number>::PreconditionBlockBase(bool      store,
                                                            Inversion method)
  : inversion(method)
  , n_diagonal_blocks(0)
  , var_store_diagonals(store)
  , var_same_diagonal(false)
  , var_inverses_ready(false)
{}


template <typename number>
inline void
PreconditionBlockBase<number>::clear()
{
  if (var_inverse_full.size() != 0)
    var_inverse_full.erase(var_inverse_full.begin(), var_inverse_full.end());
  if (var_inverse_householder.size() != 0)
    var_inverse_householder.erase(var_inverse_householder.begin(),
                                  var_inverse_householder.end());
  if (var_inverse_svd.size() != 0)
    var_inverse_svd.erase(var_inverse_svd.begin(), var_inverse_svd.end());
  if (var_diagonal.size() != 0)
    var_diagonal.erase(var_diagonal.begin(), var_diagonal.end());
  var_same_diagonal  = false;
  var_inverses_ready = false;
  n_diagonal_blocks  = 0;
}

template <typename number>
inline void
PreconditionBlockBase<number>::reinit(unsigned int n,
                                      size_type    b,
                                      bool         compress,
                                      Inversion    method)
{
  inversion          = method;
  var_same_diagonal  = compress;
  var_inverses_ready = false;
  n_diagonal_blocks  = n;

  if (compress)
    {
      switch (inversion)
        {
          case gauss_jordan:
            var_inverse_full.resize(1);
            var_inverse_full[0].reinit(b, b);
            break;
          case householder:
            var_inverse_householder.resize(1);
            break;
          case svd:
            var_inverse_svd.resize(1);
            var_inverse_svd[0].reinit(b, b);
            break;
          default:
            Assert(false, ExcNotImplemented());
        }

      if (store_diagonals())
        {
          var_diagonal.resize(1);
          var_diagonal[0].reinit(b, b);
        }
    }
  else
    {
      // set the arrays to the right
      // size. we could do it like this:
      // var_inverse = vector<>(nblocks,FullMatrix<>())
      // but this would involve copying many
      // FullMatrix objects.
      //
      // the following is a neat trick which
      // avoids copying
      if (store_diagonals())
        {
          std::vector<FullMatrix<number>> tmp(n, FullMatrix<number>(b));
          var_diagonal.swap(tmp);
        }

      switch (inversion)
        {
          case gauss_jordan:
            {
              std::vector<FullMatrix<number>> tmp(n, FullMatrix<number>(b));
              var_inverse_full.swap(tmp);
              break;
            }
          case householder:
            var_inverse_householder.resize(n);
            break;
          case svd:
            {
              std::vector<LAPACKFullMatrix<number>> tmp(
                n, LAPACKFullMatrix<number>(b));
              var_inverse_svd.swap(tmp);
              break;
            }
          default:
            Assert(false, ExcNotImplemented());
        }
    }
}


template <typename number>
inline unsigned int
PreconditionBlockBase<number>::size() const
{
  return n_diagonal_blocks;
}

template <typename number>
template <typename number2>
inline void
PreconditionBlockBase<number>::inverse_vmult(size_type              i,
                                             Vector<number2> &      dst,
                                             const Vector<number2> &src) const
{
  const size_type ii = same_diagonal() ? 0U : i;

  switch (inversion)
    {
      case gauss_jordan:
        AssertIndexRange(ii, var_inverse_full.size());
        var_inverse_full[ii].vmult(dst, src);
        break;
      case householder:
        AssertIndexRange(ii, var_inverse_householder.size());
        var_inverse_householder[ii].vmult(dst, src);
        break;
      case svd:
        AssertIndexRange(ii, var_inverse_svd.size());
        var_inverse_svd[ii].vmult(dst, src);
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
}


template <typename number>
template <typename number2>
inline void
PreconditionBlockBase<number>::inverse_Tvmult(size_type              i,
                                              Vector<number2> &      dst,
                                              const Vector<number2> &src) const
{
  const size_type ii = same_diagonal() ? 0U : i;

  switch (inversion)
    {
      case gauss_jordan:
        AssertIndexRange(ii, var_inverse_full.size());
        var_inverse_full[ii].Tvmult(dst, src);
        break;
      case householder:
        AssertIndexRange(ii, var_inverse_householder.size());
        var_inverse_householder[ii].Tvmult(dst, src);
        break;
      case svd:
        AssertIndexRange(ii, var_inverse_svd.size());
        var_inverse_svd[ii].Tvmult(dst, src);
        break;
      default:
        Assert(false, ExcNotImplemented());
    }
}


template <typename number>
inline const FullMatrix<number> &
PreconditionBlockBase<number>::inverse(size_type i) const
{
  if (same_diagonal())
    return var_inverse_full[0];

  AssertIndexRange(i, var_inverse_full.size());
  return var_inverse_full[i];
}


template <typename number>
inline const Householder<number> &
PreconditionBlockBase<number>::inverse_householder(size_type i) const
{
  if (same_diagonal())
    return var_inverse_householder[0];

  AssertIndexRange(i, var_inverse_householder.size());
  return var_inverse_householder[i];
}


template <typename number>
inline const LAPACKFullMatrix<number> &
PreconditionBlockBase<number>::inverse_svd(size_type i) const
{
  if (same_diagonal())
    return var_inverse_svd[0];

  AssertIndexRange(i, var_inverse_svd.size());
  return var_inverse_svd[i];
}


template <typename number>
inline const FullMatrix<number> &
PreconditionBlockBase<number>::diagonal(size_type i) const
{
  Assert(store_diagonals(), ExcDiagonalsNotStored());

  if (same_diagonal())
    return var_diagonal[0];

  AssertIndexRange(i, var_diagonal.size());
  return var_diagonal[i];
}


template <typename number>
inline FullMatrix<number> &
PreconditionBlockBase<number>::inverse(size_type i)
{
  Assert(var_inverse_full.size() != 0, ExcInverseNotAvailable());

  if (same_diagonal())
    return var_inverse_full[0];

  AssertIndexRange(i, var_inverse_full.size());
  return var_inverse_full[i];
}


template <typename number>
inline Householder<number> &
PreconditionBlockBase<number>::inverse_householder(size_type i)
{
  Assert(var_inverse_householder.size() != 0, ExcInverseNotAvailable());

  if (same_diagonal())
    return var_inverse_householder[0];

  AssertIndexRange(i, var_inverse_householder.size());
  return var_inverse_householder[i];
}


template <typename number>
inline LAPACKFullMatrix<number> &
PreconditionBlockBase<number>::inverse_svd(size_type i)
{
  Assert(var_inverse_svd.size() != 0, ExcInverseNotAvailable());

  if (same_diagonal())
    return var_inverse_svd[0];

  AssertIndexRange(i, var_inverse_svd.size());
  return var_inverse_svd[i];
}


template <typename number>
inline FullMatrix<number> &
PreconditionBlockBase<number>::diagonal(size_type i)
{
  Assert(store_diagonals(), ExcDiagonalsNotStored());

  if (same_diagonal())
    return var_diagonal[0];

  AssertIndexRange(i, var_diagonal.size());
  return var_diagonal[i];
}


template <typename number>
inline bool
PreconditionBlockBase<number>::same_diagonal() const
{
  return var_same_diagonal;
}


template <typename number>
inline bool
PreconditionBlockBase<number>::store_diagonals() const
{
  return var_store_diagonals;
}


template <typename number>
inline void
PreconditionBlockBase<number>::inverses_computed(bool x)
{
  var_inverses_ready = x;
}


template <typename number>
inline bool
PreconditionBlockBase<number>::inverses_ready() const
{
  return var_inverses_ready;
}


template <typename number>
inline void
PreconditionBlockBase<number>::log_statistics() const
{
  deallog << "PreconditionBlockBase: " << size() << " blocks; ";

  if (inversion == svd)
    {
      unsigned int kermin = 100000000, kermax = 0;
      double       sigmin = 1.e300, sigmax = -1.e300;
      double       kappamin = 1.e300, kappamax = -1.e300;

      for (size_type b = 0; b < size(); ++b)
        {
          const LAPACKFullMatrix<number> &matrix = inverse_svd(b);
          size_type                       k      = 1;
          while (k <= matrix.n_cols() &&
                 matrix.singular_value(matrix.n_cols() - k) == 0)
            ++k;
          const double s0 = matrix.singular_value(0);
          const double sm = matrix.singular_value(matrix.n_cols() - k);
          const double co = sm / s0;

          if (kermin > k)
            kermin = k - 1;
          if (kermax < k)
            kermax = k - 1;
          if (s0 < sigmin)
            sigmin = s0;
          if (sm > sigmax)
            sigmax = sm;
          if (co < kappamin)
            kappamin = co;
          if (co > kappamax)
            kappamax = co;
        }
      deallog << "dim ker [" << kermin << ':' << kermax << "] sigma [" << sigmin
              << ':' << sigmax << "] kappa [" << kappamin << ':' << kappamax
              << ']' << std::endl;
    }
  else if (inversion == householder)
    {}
  else if (inversion == gauss_jordan)
    {}
  else
    {
      Assert(false, ExcNotImplemented());
    }
}


template <typename number>
inline std::size_t
PreconditionBlockBase<number>::memory_consumption() const
{
  std::size_t mem = sizeof(*this);
  for (size_type i = 0; i < var_inverse_full.size(); ++i)
    mem += MemoryConsumption::memory_consumption(var_inverse_full[i]);
  for (size_type i = 0; i < var_diagonal.size(); ++i)
    mem += MemoryConsumption::memory_consumption(var_diagonal[i]);
  return mem;
}


DEAL_II_NAMESPACE_CLOSE

#endif


