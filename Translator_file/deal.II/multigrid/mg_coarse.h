//include/deal.II-translator/multigrid/mg_coarse_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_mg_coarse_h
#define dealii_mg_coarse_h


#include <deal.II/base/config.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/linear_operator.h>

#include <deal.II/multigrid/mg_base.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mg */ 
 /*@{*/ 

/**
 * 仅使用平滑器的粗略网格求解器。这是一个小的封装器，将平滑器转化为粗网格求解器。
 *
 *
 */
template <class VectorType = Vector<double>>
class MGCoarseGridApplySmoother : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * 默认构造函数。
   *
   */
  MGCoarseGridApplySmoother();

  /**
   * 构造函数。存储一个指向平滑器的指针供以后使用。
   *
   */
  MGCoarseGridApplySmoother(const MGSmootherBase<VectorType> &coarse_smooth);

  /**
   * 清除该指针。
   *
   */
  void
  clear();

  /**
   * 初始化新的数据。
   *
   */
  void
  initialize(const MGSmootherBase<VectorType> &coarse_smooth);

  /**
   * 实现抽象函数。
   *
   */
  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

private:
  /**
   * 对平滑器的引用。
   *
   */
  SmartPointer<const MGSmootherBase<VectorType>,
               MGCoarseGridApplySmoother<VectorType>>
    coarse_smooth;
};



/**
 * 迭代求解器的粗网格多格运算器。
 * 该类为具有给定矩阵和预调节器的deal.II迭代求解器提供了一个包装，作为粗网格算子。
 *
 *
 */
template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
class MGCoarseGridIterativeSolver : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * 默认构造函数。
   *
   */
  MGCoarseGridIterativeSolver();

  /**
   * 构造器。只存储这些对象的引用，所以它们的寿命需要超过这个类中的用量。
   *
   */
  MGCoarseGridIterativeSolver(SolverType &              solver,
                              const MatrixType &        matrix,
                              const PreconditionerType &precondition);

  /**
   * 用新的数据进行初始化，更多细节见相应的构造函数。
   *
   */
  void
  initialize(SolverType &              solver,
             const MatrixType &        matrix,
             const PreconditionerType &precondition);

  /**
   * 清除所有的指针。
   *
   */
  void
  clear();

  /**
   * 抽象函数的实现。用矩阵、向量和预调节器调用给定求解器的求解方法。
   *
   */
  virtual void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

private:
  /**
   * 对解算器的引用。
   *
   */
  SmartPointer<SolverType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    solver;

  /**
   * 对矩阵的参考。
   *
   */
  SmartPointer<const MatrixType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    matrix;

  /**
   * 对预处理程序的引用。
   *
   */
  SmartPointer<const PreconditionerType,
               MGCoarseGridIterativeSolver<VectorType,
                                           SolverType,
                                           MatrixType,
                                           PreconditionerType>>
    preconditioner;
};



/**
 * 在Householder类中通过QR分解实现的粗网格求解器。
 * 初始化时，计算矩阵的QR分解。然后，运算器()使用
 * Householder::least_squares() 来计算逆的作用。
 *
 *
 */
template <typename number = double, class VectorType = Vector<number>>
class MGCoarseGridHouseholder : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * 构造函数，采取粗略的网格矩阵。
   *
   */
  MGCoarseGridHouseholder(const FullMatrix<number> *A = nullptr);

  /**
   * 为一个新的矩阵进行初始化。
   *
   */
  void
  initialize(const FullMatrix<number> &A);

  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const override;

private:
  /**
   * 用于QR因子化的矩阵。
   *
   */
  Householder<number> householder;
};

/**
 * 使用LAPACK矩阵的奇异值分解的粗略网格求解器。
 * 初始化时，计算矩阵的奇异值分解。然后，运算器（）使用
 *
 *
 */
template <typename number = double, class VectorType = Vector<number>>
class MGCoarseGridSVD : public MGCoarseGridBase<VectorType>
{
public:
  /**
   * 构造函数留下一个未初始化的对象。
   *
   */
  MGCoarseGridSVD() = default;

  /**
   * 对一个新的矩阵进行初始化。这将重设尺寸为
   *
   */
  void
  initialize(const FullMatrix<number> &A, const double threshold = 0);

  void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const;

  /**
   * 将奇异值写到 @p deallog. 。
   *
   */
  void
  log() const;

private:
  /**
   * 用于奇异值分解的矩阵。
   *
   */
  LAPACKFullMatrix<number> matrix;
};

 /*@}*/ 

#ifndef DOXYGEN
 /* ------------------ Functions for MGCoarseGridApplySmoother -----------*/ 
template <class VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother()
  : coarse_smooth(nullptr)
{}

template <class VectorType>
MGCoarseGridApplySmoother<VectorType>::MGCoarseGridApplySmoother(
  const MGSmootherBase<VectorType> &coarse_smooth)
  : coarse_smooth(nullptr)
{
  initialize(coarse_smooth);
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::initialize(
  const MGSmootherBase<VectorType> &coarse_smooth_)
{
  coarse_smooth =
    SmartPointer<const MGSmootherBase<VectorType>,
                 MGCoarseGridApplySmoother<VectorType>>(&coarse_smooth_,
                                                        typeid(*this).name());
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::clear()
{
  coarse_smooth = nullptr;
}


template <class VectorType>
void
MGCoarseGridApplySmoother<VectorType>::operator()(const unsigned int level,
                                                  VectorType &       dst,
                                                  const VectorType & src) const
{
  coarse_smooth->smooth(level, dst, src);
}

 /* ------------------ Functions for MGCoarseGridIterativeSolver ------------ */ 

template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::MGCoarseGridIterativeSolver()
  : solver(0, typeid(*this).name())
  , matrix(0, typeid(*this).name())
  , preconditioner(0, typeid(*this).name())
{}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::
  MGCoarseGridIterativeSolver(SolverType &              solver,
                              const MatrixType &        matrix,
                              const PreconditionerType &preconditioner)
  : solver(&solver, typeid(*this).name())
  , matrix(&matrix, typeid(*this).name())
  , preconditioner(&preconditioner, typeid(*this).name())
{}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<
  VectorType,
  SolverType,
  MatrixType,
  PreconditionerType>::initialize(SolverType &              solver_,
                                  const MatrixType &        matrix_,
                                  const PreconditionerType &preconditioner_)
{
  solver         = &solver_;
  matrix         = &matrix_;
  preconditioner = &preconditioner_;
}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::clear()
{
  solver         = 0;
  matrix         = 0;
  preconditioner = 0;
}



template <class VectorType,
          class SolverType,
          class MatrixType,
          class PreconditionerType>
void
MGCoarseGridIterativeSolver<VectorType,
                            SolverType,
                            MatrixType,
                            PreconditionerType>::
operator()(const unsigned int  /*level*/ ,
           VectorType &      dst,
           const VectorType &src) const
{
  Assert(solver != nullptr, ExcNotInitialized());
  Assert(matrix != nullptr, ExcNotInitialized());
  Assert(preconditioner != nullptr, ExcNotInitialized());
  solver->solve(*matrix, dst, src, *preconditioner);
}



 /* ------------------ Functions for MGCoarseGridHouseholder ------------ */ 

template <typename number, class VectorType>
MGCoarseGridHouseholder<number, VectorType>::MGCoarseGridHouseholder(
  const FullMatrix<number> *A)
{
  if (A != nullptr)
    householder.initialize(*A);
}



template <typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::initialize(
  const FullMatrix<number> &A)
{
  householder.initialize(A);
}



template <typename number, class VectorType>
void
MGCoarseGridHouseholder<number, VectorType>::
operator()(const unsigned int  /*level*/ ,
           VectorType &      dst,
           const VectorType &src) const
{
  householder.least_squares(dst, src);
}

//---------------------------------------------------------------------------



template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::initialize(const FullMatrix<number> &A,
                                                double threshold)
{
  matrix.reinit(A.n_rows(), A.n_cols());
  matrix = A;
  matrix.compute_inverse_svd(threshold);
}


template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::operator()(const unsigned int  /*level*/ ,
                                                VectorType &      dst,
                                                const VectorType &src) const
{
  matrix.vmult(dst, src);
}


template <typename number, class VectorType>
void
MGCoarseGridSVD<number, VectorType>::log() const
{
  const unsigned int n = std::min(matrix.n_rows(), matrix.n_cols());

  for (unsigned int i = 0; i < n; ++i)
    deallog << ' ' << matrix.singular_value(i);
  deallog << std::endl;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


