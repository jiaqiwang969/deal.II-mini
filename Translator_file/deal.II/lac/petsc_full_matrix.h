//include/deal.II-translator/lac/petsc_full_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_petsc_full_matrix_h
#  define dealii_petsc_full_matrix_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_matrix_base.h>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  /*!   @addtogroup  PETScWrappers  
     * @{   
*
*/

  /**
   * 实现一个基于PETSc的连续密集矩阵类。所有的功能实际上都在基类中，除了生成连续密集矩阵的调用。这是可能的，因为PETSc只在一个抽象的矩阵类型上工作，并在内部根据实际的矩阵类型分配给做实际工作的函数（很像使用虚拟函数）。只有创建特定类型矩阵的函数不同，并在这个特定的类中实现。
   * @ingroup Matrix1
   *
   */
  class FullMatrix : public MatrixBase
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = types::global_dof_index;


    /**
     * 默认构造函数。创建一个空矩阵。
     *
     */
    FullMatrix();


    /**
     * 创建一个尺寸为 @p m 乘以 @p n. 的完整矩阵。
     *
     */
    FullMatrix(const size_type m, const size_type n);


    /**
     * 扔掉现在的矩阵，并生成一个具有相同属性的矩阵，就像它是由这个类的构造函数创建的，参数列表与现在的函数相同。
     *
     */
    void
    reinit(const size_type m, const size_type n);


    /**
     * 返回一个对与本矩阵一起使用的MPI通信器对象的引用。由于这是一个连续的矩阵，它返回MPI_COMM_SELF通信器。
     *
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const override;

  private:
    /**
     * 为各自的reinit()函数和匹配的构造函数做实际工作，即创建一个矩阵。摆脱之前的矩阵是留给调用者的。
     *
     */
    void
    do_reinit(const size_type m, const size_type n);
  };

   /*@}*/ 
} // namespace PETScWrappers


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
 /*---------------------------- petsc_full_matrix.h --------------------------*/ 


