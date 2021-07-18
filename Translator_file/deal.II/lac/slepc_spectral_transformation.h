//include/deal.II-translator/lac/slepc_spectral_transformation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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


#ifndef dealii_slepc_spectral_transformation_h
#  define dealii_slepc_spectral_transformation_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_SLEPC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/petsc_solver.h>

#    include <petscksp.h>

#    include <slepceps.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#    ifndef DOXYGEN
namespace PETScWrappers
{
  // forward declarations
  class SolverBase;
} // namespace PETScWrappers
#    endif

namespace SLEPcWrappers
{
  // forward declaration
  class SolverBase;

  /**
   * 使用SLEPc求解器的谱系转换类的基类，这些求解器是根据传递给谱系转换的标志来选择的。
   * <code>SLEPcWrappers::TransformationXXX</code>, where <code>XXX</code>
   * 是你最喜欢的变换类型，然后可以在应用程序代码中以下列方式实现
   * <code>XXX=INVERT</code>  与求解器对象  <code>eigensolver</code>  :
   * @code
   * // Set a transformation, this one shifts the eigenspectrum by 3.142..
   * SLEPcWrappers::TransformationShift::AdditionalData
   * additional_data(3.142);
   * SLEPcWrappers::TransformationShift shift(mpi_communicator,additional_data);
   * eigensolver.set_transformation(shift);
   * @endcode
   * 之后像往常一样调用 <code>solve()</code> 函数。
   * @code
   * SolverControl solver_control (1000, 1e-9);
   * SolverArnoldi system (solver_control, mpi_communicator);
   * eigensolver.solve (A, B, lambda, x, size_of_spectrum);
   * @endcode
   *
   * @note  这些选项也可以在命令行中设置。
   * @ingroup SLEPcWrappers
   *
   */
  class TransformationBase
  {
  protected:
    /**
     * 构造器。
     *
     */
    TransformationBase(const MPI_Comm &mpi_communicator);

  public:
    /**
     * 解构器。
     *
     */
    virtual ~TransformationBase();

    /**
     * 设置一个标志，以表明变换后的矩阵是如何被存储在光谱变换中的。
     * 可能的值由SLEPc库中的枚举器STMatMode给出
     * http://www.grycap.upv.es/slepc/documentation/current/docs/manualpages/ST/STMatMode.html
     *
     */
    void
    set_matrix_mode(const STMatMode mode);

    /**
     * 设置在eigensolver内求解线性代数方程组时要使用的解算器。
     *
     */
    void
    set_solver(const PETScWrappers::SolverBase &solver);

  protected:
    /**
     * SLEPc谱系转换对象。
     *
     */
    ST st;

    // Make the solver class a friend, since it needs to set spectral
    // transformation object.
    friend class SolverBase;
  };

  /**
   * 使用SLEPc转换接口的实现。
   * @ingroup SLEPcWrappers
   *
   */
  class TransformationShift : public TransformationBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将移位参数设置为零。
       *
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * 移位参数。
       *
       */
      const double shift_parameter;
    };


    /**
     * 构造函数。
     *
     */
    TransformationShift(const MPI_Comm &      mpi_communicator,
                        const AdditionalData &data = AdditionalData());


  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 一个使用SLEPc移位和反转的转换接口的实现。
   * @ingroup SLEPcWrappers
   *
   */
  class TransformationShiftInvert : public TransformationBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将移位参数设置为零。
       *
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * 移位参数。
       *
       */
      const double shift_parameter;
    };


    /**
     * 构造函数。
     *
     */
    TransformationShiftInvert(const MPI_Comm &      mpi_communicator,
                              const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;

    // Make the solver class a friend, since it may need to set target
    // equal the provided shift value.
    friend class SolverBase;
  };

  /**
   * 使用SLEPc频谱折叠的转换接口的实现。这种变换类型在SLEPc
   * 3.5.0中已经被删除，因此不能在更新的版本中使用。
   * @ingroup SLEPcWrappers
   *
   */
  class TransformationSpectrumFolding : public TransformationBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将移位参数设置为零。
       *
       */
      AdditionalData(const double shift_parameter = 0);

      /**
       * 移位参数。
       *
       */
      const double shift_parameter;
    };


    /**
     * 构造函数。
     *
     */
    TransformationSpectrumFolding(
      const MPI_Comm &      mpi_communicator,
      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };

  /**
   * 使用SLEPc Cayley的转换接口的实现。
   * @ingroup SLEPcWrappers
   *
   */
  class TransformationCayley : public TransformationBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。需要两个移位参数
       *
       */
      AdditionalData(const double shift_parameter     = 0,
                     const double antishift_parameter = 0);

      /**
       * 移位参数。
       *
       */
      const double shift_parameter;

      /**
       * 反移位参数。
       *
       */
      const double antishift_parameter;
    };


    /**
     * 构造函数。
     *
     */
    TransformationCayley(const MPI_Comm &      mpi_communicator,
                         const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_SLEPC

 /*--------------------   slepc_spectral_transformation.h   ------------------*/ 

#endif

 /*--------------------   slepc_spectral_transformation.h   ------------------*/ 


