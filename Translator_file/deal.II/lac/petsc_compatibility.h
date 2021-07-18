//include/deal.II-translator/lac/petsc_compatibility_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

/* 与其到处使用ifdefs，不如尝试将PETSc的旧版本函数包裹在一个地方。

* 
*
*/
#ifndef dealii_petsc_compatibility_h
#define dealii_petsc_compatibility_h

#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>

#ifdef DEAL_II_WITH_PETSC

#  include <petscconf.h>
#  include <petscksp.h>
#  include <petscmat.h>
#  include <petscpc.h>

#  include <string>

DEAL_II_NAMESPACE_OPEN

namespace PETScWrappers
{
  /**
   * 在全局PETSc数据库中设置一个选项。这个函数只是用版本检查来包装
   * PetscOptionsSetValue（这个函数的签名在 PETSc 3.7.0
   * 中改变了）。
   *
   */
  inline void
  set_option_value(const std::string &name, const std::string &value)
  {
#  if DEAL_II_PETSC_VERSION_LT(3, 7, 0)
    const PetscErrorCode ierr =
      PetscOptionsSetValue(name.c_str(), value.c_str());
#  else
    const PetscErrorCode ierr =
      PetscOptionsSetValue(nullptr, name.c_str(), value.c_str());
#  endif
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  /**
   * 销毁一个PETSc矩阵。这个函数封装了MatDestroy，并进行了版本检查（这个函数的签名在PETSc
   * 3.2.0中改变）。      @warning
   * 由于这个函数的主要意图是在PETSc封装器中启用RAII语义，如果发生错误，这个函数将不会抛出一个异常，而只是返回MatDestroy给出的错误代码。
   *
   */
  inline PetscErrorCode
  destroy_matrix(Mat &matrix)
  {
    // PETSc will check whether or not matrix is nullptr.
    return MatDestroy(&matrix);
  }



  /**
   * 销毁一个Krylov Subspace (KSP)
   * PETSc求解器。这个函数用一个版本检查来包装KSPDestroy（这个函数的签名在PETSc
   * 3.2.0中改变了）。      @warning
   * 由于该函数的主要意图是在PETSc封装器中启用RAII语义，因此如果发生错误，该函数将不会抛出一个异常，而只是返回MatDestroy给出的错误代码。
   *
   */
  inline PetscErrorCode
  destroy_krylov_solver(KSP &krylov_solver)
  {
    // PETSc will check whether or not matrix is nullptr.
    return KSPDestroy(&krylov_solver);
  }



  /**
   * 设置一个PETSc矩阵选项。这个函数用一个版本检查来包装MatSetOption。
   * @warning
   * 在3.0.0之前的PETSc版本中，参数option_value被忽略，因为相应的函数不接受这个参数。
   *
   */
  inline void
  set_matrix_option(Mat &           matrix,
                    const MatOption option_name,
                    const PetscBool option_value = PETSC_FALSE)
  {
    const PetscErrorCode ierr = MatSetOption(matrix, option_name, option_value);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
  }



  /**
   * 告诉PETSc，我们不打算向矩阵添加新的条目。在调试模式下产生错误。
   *
   */
  inline void
  close_matrix(Mat &matrix)
  {
#  ifdef DEBUG
    set_matrix_option(matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
#  else
    set_matrix_option(matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
#  endif
  }



  /**
   * 告诉PETSc保留SparsityPattern条目，即使我们用clear_rows()删除了某行，并调用MatZeroRows()。否则，以后就不能再向该行写入数据了。
   *
   */
  inline void
  set_keep_zero_rows(Mat &matrix)
  {
    set_matrix_option(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
  }
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_PETSC
#endif // dealii_petsc_compatibility_h


