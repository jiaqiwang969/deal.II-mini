//include/deal.II-translator/multigrid/mg_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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

#ifndef dealii_mg_matrix_h
#define dealii_mg_matrix_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_base.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mg */ 
 /*@{*/ 

namespace mg
{
  /**
   * 多层次矩阵。该矩阵存储了LinearOperator对象的MGLevelObject。它实现了MGMatrixBase中定义的接口，因此它可以作为Multigrid中的一个矩阵使用。
   *
   */
  template <typename VectorType = Vector<double>>
  class Matrix : public MGMatrixBase<VectorType>
  {
  public:
    /**
     * 空对象的默认构造函数。
     *
     */
    Matrix() = default;

    /**
     * 构造函数通过调用initialize()设置指向<tt>M</tt>中的矩阵的指针。
     *
     */
    template <typename MatrixType>
    Matrix(const MGLevelObject<MatrixType> &M);

    /**
     * 初始化该对象，使水平乘法使用<tt>M</tt>中的矩阵。
     *
     */
    template <typename MatrixType>
    void
    initialize(const MGLevelObject<MatrixType> &M);

    /**
     * 重置该对象。
     *
     */
    void
    reset();

    /**
     * 在一个层次上访问矩阵。
     *
     */
    const LinearOperator<VectorType> &operator[](unsigned int level) const;

    virtual void
    vmult(const unsigned int level,
          VectorType &       dst,
          const VectorType & src) const override;
    virtual void
    vmult_add(const unsigned int level,
              VectorType &       dst,
              const VectorType & src) const override;
    virtual void
    Tvmult(const unsigned int level,
           VectorType &       dst,
           const VectorType & src) const override;
    virtual void
    Tvmult_add(const unsigned int level,
               VectorType &       dst,
               const VectorType & src) const override;
    virtual unsigned int
    get_minlevel() const override;
    virtual unsigned int
    get_maxlevel() const override;

    /**
     * 此对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;

  private:
    MGLevelObject<LinearOperator<VectorType>> matrices;
  };

} // namespace mg


/**
 * 从块状矩阵中选择多级矩阵。该类实现了MGMatrixBase定义的接口。
 * 模板参数 @p MatrixType
 * 应该是一个块状矩阵类，如BlockSparseMatrix或 @p
 * BlockSparseMatrixEZ。然后，该类存储一个指向该矩阵类的MGLevelObject的指针。在每个
 * @p vmult, 中，初始化时选择的块将与提供的向量相乘。
 *
 *
 */
template <typename MatrixType, typename number>
class MGMatrixSelect : public MGMatrixBase<Vector<number>>
{
public:
  /**
   * 构造函数。  @p row 和 @p col 是选择块的坐标。
   * 另一个参数被移交给 @p SmartPointer 构造函数。
   *
   */
  MGMatrixSelect(const unsigned int         row    = 0,
                 const unsigned int         col    = 0,
                 MGLevelObject<MatrixType> *matrix = 0);

  /**
   * 设置要使用的矩阵对象。矩阵对象必须作为 @p
   * MGMatrixSelect
   * 对象存在更长时间，因为只存储了一个指针。
   *
   */
  void
  set_matrix(MGLevelObject<MatrixType> *M);

  /**
   * 选择用于乘法的块。
   *
   */
  void
  select_block(const unsigned int row, const unsigned int col);

  /**
   * 矩阵-向量-乘法的某一层次。
   *
   */
  virtual void
  vmult(const unsigned int    level,
        Vector<number> &      dst,
        const Vector<number> &src) const;

  /**
   * 在某一层次上添加矩阵-向量-乘法。
   *
   */
  virtual void
  vmult_add(const unsigned int    level,
            Vector<number> &      dst,
            const Vector<number> &src) const;

  /**
   * 在某一层次上进行矩阵-向量-乘法的转置。
   *
   */
  virtual void
  Tvmult(const unsigned int    level,
         Vector<number> &      dst,
         const Vector<number> &src) const;

  /**
   * 在某一层次上增加转置矩阵-向量-乘法。
   *
   */
  virtual void
  Tvmult_add(const unsigned int    level,
             Vector<number> &      dst,
             const Vector<number> &src) const;

private:
  /**
   * 指向每一层的矩阵对象的指针。
   *
   */
  SmartPointer<MGLevelObject<MatrixType>, MGMatrixSelect<MatrixType, number>>
    matrix;
  /**
   * 所选块的行坐标。
   *
   */
  unsigned int row;
  /**
   * 所选区块的列坐标。
   *
   */
  unsigned int col;
};

 /*@}*/ 

 /*----------------------------------------------------------------------*/ 

namespace mg
{
  template <typename VectorType>
  template <typename MatrixType>
  inline void
  Matrix<VectorType>::initialize(const MGLevelObject<MatrixType> &p)
  {
    matrices.resize(p.min_level(), p.max_level());
    for (unsigned int level = p.min_level(); level <= p.max_level(); ++level)
      {
        // Workaround: Unfortunately, not every "p[level]" object has a
        // rich enough interface to populate reinit_(domain|range)_vector.
        // Thus, apply an empty LinearOperator exemplar.
        matrices[level] =
          linear_operator<VectorType>(LinearOperator<VectorType>(),
                                      Utilities::get_underlying_value(
                                        p[level]));
      }
  }



  template <typename VectorType>
  inline void
  Matrix<VectorType>::reset()
  {
    matrices.resize(0, 0);
  }



  template <typename VectorType>
  template <typename MatrixType>
  inline Matrix<VectorType>::Matrix(const MGLevelObject<MatrixType> &p)
  {
    initialize(p);
  }



  template <typename VectorType>
  inline const LinearOperator<VectorType> &Matrix<VectorType>::
                                           operator[](unsigned int level) const
  {
    return matrices[level];
  }



  template <typename VectorType>
  void
  Matrix<VectorType>::vmult(const unsigned int level,
                            VectorType &       dst,
                            const VectorType & src) const
  {
    matrices[level].vmult(dst, src);
  }



  template <typename VectorType>
  void
  Matrix<VectorType>::vmult_add(const unsigned int level,
                                VectorType &       dst,
                                const VectorType & src) const
  {
    matrices[level].vmult_add(dst, src);
  }



  template <typename VectorType>
  void
  Matrix<VectorType>::Tvmult(const unsigned int level,
                             VectorType &       dst,
                             const VectorType & src) const
  {
    matrices[level].Tvmult(dst, src);
  }



  template <typename VectorType>
  void
  Matrix<VectorType>::Tvmult_add(const unsigned int level,
                                 VectorType &       dst,
                                 const VectorType & src) const
  {
    matrices[level].Tvmult_add(dst, src);
  }



  template <typename VectorType>
  unsigned int
  Matrix<VectorType>::get_minlevel() const
  {
    return matrices.min_level();
  }



  template <typename VectorType>
  unsigned int
  Matrix<VectorType>::get_maxlevel() const
  {
    return matrices.max_level();
  }



  template <typename VectorType>
  inline std::size_t
  Matrix<VectorType>::memory_consumption() const
  {
    return sizeof(*this) + matrices->memory_consumption();
  }
} // namespace mg


 /*----------------------------------------------------------------------*/ 

template <typename MatrixType, typename number>
MGMatrixSelect<MatrixType, number>::MGMatrixSelect(const unsigned int row,
                                                   const unsigned int col,
                                                   MGLevelObject<MatrixType> *p)
  : matrix(p, typeid(*this).name())
  , row(row)
  , col(col)
{}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::set_matrix(MGLevelObject<MatrixType> *p)
{
  matrix = p;
}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::select_block(const unsigned int brow,
                                                 const unsigned int bcol)
{
  row = brow;
  col = bcol;
}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::vmult(const unsigned int    level,
                                          Vector<number> &      dst,
                                          const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).vmult(dst, src);
}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::vmult_add(const unsigned int    level,
                                              Vector<number> &      dst,
                                              const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).vmult_add(dst, src);
}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::Tvmult(const unsigned int    level,
                                           Vector<number> &      dst,
                                           const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).Tvmult(dst, src);
}



template <typename MatrixType, typename number>
void
MGMatrixSelect<MatrixType, number>::Tvmult_add(const unsigned int    level,
                                               Vector<number> &      dst,
                                               const Vector<number> &src) const
{
  Assert(matrix != 0, ExcNotInitialized());

  const MGLevelObject<MatrixType> &m = *matrix;
  m[level].block(row, col).Tvmult_add(dst, src);
}

DEAL_II_NAMESPACE_CLOSE

#endif


