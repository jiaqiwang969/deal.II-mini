//include/deal.II-translator/lac/petsc_matrix_free_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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

#ifndef dealii_petsc_matrix_free_h
#  define dealii_petsc_matrix_free_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC
#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_matrix_base.h>
#    include <deal.II/lac/petsc_vector.h>
DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  /**
   * 实现一个基于PETSc <tt>MatShell</tt>
   * 矩阵类型的并行矩阵类。这个基类只实现了PETSc矩阵对象的接口，而所有的功能都包含在矩阵-向量乘法中，必须在派生类中重新实现。
   * 这个接口是对 dealii::MatrixFree
   * 类的补充，以实现用户定义的矩阵类与PETSc解算器和功能。参见
   * dealii::MatrixFree  类和  step-37  以及  step-48  的文档。
   * 与命名空间PETScWrappers和 PETScWrappers::MPI,
   * 中的其他矩阵类类似，MatrixFree类提供了通常的矩阵-向量乘法<tt>vmult(VectorBase
   * &dst, const VectorBase
   * &src)</tt>，它是纯虚拟的，必须在派生类中重新实现。
   * 除了通常的接口，这个类还有一个矩阵向量乘法<tt>vmult(Vec
   * &dst, const Vec &src)</tt>，采取PETSc
   * Vec对象，它将被<tt>matrix_free_mult(Mat A, Vec src, Vec
   * dst)</tt>注册为这个PETSc矩阵对象的矩阵-向量乘法。基类中vmult函数的默认实现是用
   * PETScWrappers::VectorBase
   * 类包装给定的PETSc向量，然后用常规接口调用常规vmult函数。
   * @ingroup PETScWrappers
   * @ingroup Matrix1
   *
   */
  class MatrixFree : public MatrixBase
  {
  public:
    /**
     * 默认构造函数。创建一个空的矩阵对象。
     *
     */
    MatrixFree();

    /**
     * 创建一个尺寸为 @p m 乘以 @p n
     * 的矩阵对象，通过提供的 @p communicator. 进行通信。
     * 关于 @p local_rows 和 @p local_columns 参数的含义，请参阅
     * PETScWrappers::MPI::SparseMatrix 类文档。
     * 与其他PETSc矩阵一样，无矩阵对象也需要有一个大小，并能有效地并行执行矩阵向量乘法
     * @p local_rows 和 @p local_columns. ，但与 PETSc::SparseMatrix
     * 类相反，PETSc无矩阵对象不需要对非零项进行任何估计，也没有<tt>is_symmetric</tt>选项。
     *
     */
    MatrixFree(const MPI_Comm &   communicator,
               const unsigned int m,
               const unsigned int n,
               const unsigned int local_rows,
               const unsigned int local_columns);

    /**
     * 创建一个尺寸为 @p m 乘以 @p n
     * 的矩阵对象，并通过提供的 @p communicator. 进行通信。
     * 与其他PETSc矩阵一样，无矩阵对象也需要有一个尺寸，并能有效地并行执行矩阵向量乘法，但与
     * PETSc::SparseMatrix
     * 类相反，PETSc无矩阵对象不需要对非零项进行任何估计，也没有<tt>is_symmetric</tt>选项。
     *
     */
    MatrixFree(const MPI_Comm &                 communicator,
               const unsigned int               m,
               const unsigned int               n,
               const std::vector<unsigned int> &local_rows_per_process,
               const std::vector<unsigned int> &local_columns_per_process,
               const unsigned int               this_process);

    /**
     * 串行情况下的构造函数。与<tt>MatrixFree()</tt>函数相同，见上文，<tt>communicator
     * = MPI_COMM_WORLD</tt>。
     *
     */
    MatrixFree(const unsigned int m,
               const unsigned int n,
               const unsigned int local_rows,
               const unsigned int local_columns);

    /**
     * 串行情况下的构造函数。与<tt>MatrixFree()</tt>函数相同，见上文，<tt>communicator
     * = MPI_COMM_WORLD</tt>。
     *
     */
    MatrixFree(const unsigned int               m,
               const unsigned int               n,
               const std::vector<unsigned int> &local_rows_per_process,
               const std::vector<unsigned int> &local_columns_per_process,
               const unsigned int               this_process);

    /**
     * 扔掉当前的矩阵，并生成一个具有相同属性的矩阵，就像它是由这个类的构造函数创建的，参数列表与当前函数相同。
     *
     */
    void
    reinit(const MPI_Comm &   communicator,
           const unsigned int m,
           const unsigned int n,
           const unsigned int local_rows,
           const unsigned int local_columns);

    /**
     * 扔掉当前的矩阵，并生成一个具有相同属性的矩阵，就像它是由这个类的构造函数创建的一样，参数列表与当前函数相同。
     *
     */
    void
    reinit(const MPI_Comm &                 communicator,
           const unsigned int               m,
           const unsigned int               n,
           const std::vector<unsigned int> &local_rows_per_process,
           const std::vector<unsigned int> &local_columns_per_process,
           const unsigned int               this_process);

    /**
     * 用<tt>communicator = MPI_COMM_WORLD</tt>调用上面的 @p reinit()
     * 函数。
     *
     */
    void
    reinit(const unsigned int m,
           const unsigned int n,
           const unsigned int local_rows,
           const unsigned int local_columns);

    /**
     * 用<tt>communicator = MPI_COMM_WORLD</tt>调用上面的 @p reinit()
     * 函数。
     *
     */
    void
    reinit(const unsigned int               m,
           const unsigned int               n,
           const std::vector<unsigned int> &local_rows_per_process,
           const std::vector<unsigned int> &local_columns_per_process,
           const unsigned int               this_process);

    /**
     * 释放所有内存并返回到与调用默认构造函数后一样的状态。
     *
     */
    void
    clear();

    /**
     * 返回一个对该矩阵使用的MPI通信器对象的引用。
     *
     */
    const MPI_Comm &
    get_mpi_communicator() const override;

    /**
     * 矩阵-向量乘法：让<i>dst =
     * M*src</i>与<i>M</i>是这个矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布式矩阵（类型为
     * PETScWrappers::MPI::SparseMatrix),
     * ，那么两个向量也必须是分布式向量。反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    virtual void
    vmult(VectorBase &dst, const VectorBase &src) const = 0;

    /**
     * 矩阵-向量乘法：让<i>dst =
     * M<sup>T</sup>*src</i>与<i>M</i>为这个矩阵。这个函数与 @p
     * vmult() 的作用相同，但需要转置的矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布的矩阵，那么两个向量也必须是分布式向量。
     * 反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    virtual void
    Tvmult(VectorBase &dst, const VectorBase &src) const = 0;

    /**
     * 加法
     * 矩阵-向量乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>为该矩阵。
     * 源和目的不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布的矩阵，那么两个向量也必须是分布式向量。
     * 反之，如果矩阵不是分布式的，那么两个向量都不能是。
     *
     */
    virtual void
    vmult_add(VectorBase &dst, const VectorBase &src) const = 0;

    /**
     * 加法
     * 矩阵-向量乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数与
     * @p vmult_add() 的作用相同，但需要转置的矩阵。
     * 来源和目的地不能是同一个向量。
     * 注意，如果当前对象代表一个平行分布的矩阵，那么两个向量也必须是分布式向量。
     * 反之，如果矩阵不是分布式的，那么两个向量都不可以是。
     *
     */
    virtual void
    Tvmult_add(VectorBase &dst, const VectorBase &src) const = 0;

    /**
     * @p matrix_free_mult(). 调用的矩阵-向量乘法
     * 这个函数可以在派生类中重新实现以提高效率。默认实现将给定的向量复制到
     * PETScWrappers::*::Vector 中，并调用<tt>vmult(VectorBase &dst,
     * const VectorBase
     * &src)</tt>，它是纯虚的，必须在派生类中重新实现。
     *
     */
    virtual void
    vmult(Vec &dst, const Vec &src) const;

  private:
    /**
     * 将用于该并行矩阵-自由对象的通信器对象的副本。
     *
     */
    MPI_Comm communicator;

    /**
     * 注册为PETSc例程所调用的该无矩阵对象的矩阵-向量乘法的回调函数。这个函数必须是静态的，并且接受一个PETSc矩阵
     * @p A, 和向量 @p src 和 @p dst, ，其中<i>dst =
     * A*src</i>源和目的地不能是同一个向量。
     * 这个函数调用<tt>vmult(Vec &dst, const Vec
     * &src)</tt>，应该在派生类中重新实现。
     *
     */
    static int
    matrix_free_mult(Mat A, Vec src, Vec dst);

    /**
     * 为各自的 @p reinit()
     * 函数和匹配的构造函数做实际工作，即创建一个矩阵对象。摆脱之前的矩阵是留给调用者的。
     *
     */
    void
    do_reinit(const unsigned int m,
              const unsigned int n,
              const unsigned int local_rows,
              const unsigned int local_columns);
  };



  // -------- template and inline functions ----------

  inline const MPI_Comm &
  MatrixFree::get_mpi_communicator() const
  {
    return communicator;
  }
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
 /*---------------------------- petsc_matrix_free.h --------------------------*/ 


