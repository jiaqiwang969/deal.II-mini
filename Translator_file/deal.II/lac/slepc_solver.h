//include/deal.II-translator/lac/slepc_solver_0.txt
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


#ifndef dealii_slepc_solver_h
#  define dealii_slepc_solver_h

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_SLEPC

#    include <deal.II/lac/petsc_matrix_base.h>
#    include <deal.II/lac/slepc_spectral_transformation.h>
#    include <deal.II/lac/solver_control.h>

#    include <petscconf.h>
#    include <petscksp.h>

#    include <slepceps.h>

#    include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 * 使用SLEPc求解器的求解器类的基础命名空间，这些求解器是根据传递给特征值问题求解器上下文的标志来选择的。派生类设置正确的标志来设置正确的求解器。
 * SLEPc求解器旨在用于解决广义的特征谱问题 $(A-\lambda B)x=0$
 * ，用于 $x\neq0$ ；其中 $A$ 是一个系统矩阵， $B$
 * 是一个质量矩阵， $\lambda, x$
 * 分别为一组特征值和特征向量。重点是适合相关矩阵是稀疏的问题的方法和技术。SLEPc库提供的大多数方法都是投影方法或其他具有类似属性的方法；并提供了包装器来连接处理这两个问题集的SLEPc求解器。
 * SLEPcWrappers可以通过以下方式在应用程序代码中实现。
 *
 * @code
 * SolverControl solver_control (1000, 1e-9);
 * SolverArnoldi system (solver_control, mpi_communicator);
 * system.solve (A, B, lambda, x, size_of_spectrum);
 * @endcode
 * 用于广义特征值问题  $Ax=B\lambda x$  ，其中变量  <code>const
 * unsigned int size_of_spectrum</code>
 * 告诉SLEPc要解决的特征向量/特征值对的数量。额外的选项和求解器参数可以在调用
 * <code>solve()</code>
 * 之前传递给SLEPc求解器。例如，如果一般特征谱问题的矩阵不是隐性的，而且只想要低级特征值，那么在调用
 * <code>solve()</code>  之前可以实现以下代码。
 *
 * @code
 * system.set_problem_type (EPS_NHEP);
 * system.set_which_eigenpairs (EPS_SMALLEST_REAL);
 * @endcode
 * 这些选项也可以在命令行中设置。 也可参见
 * <code>step-36</code> 中的实战案例。
 * 对于谱系转换与Krylov类型求解器或Davidson类型求解器结合使用的情况，可以额外指定使用哪种线性求解器和预处理器。这可以通过以下方式实现
 *
 * @code
 * PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
 * data.symmetric_operator = true;
 * PETScWrappers::PreconditionBoomerAMG preconditioner(mpi_communicator, data);
 * SolverControl linear_solver_control (dof_handler.n_dofs(),
 *                                    1e-12, false, false);
 * PETScWrappers::SolverCG linear_solver(linear_solver_control,
 *                                     mpi_communicator);
 * linear_solver.initialize(preconditioner);
 * SolverControl solver_control (100, 1e-9,false,false);
 * SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control,
 *                                            mpi_communicator);
 * SLEPcWrappers::TransformationShift spectral_transformation(mpi_communicator);
 * spectral_transformation.set_solver(linear_solver);
 * eigensolver.set_transformation(spectral_transformation);
 * eigensolver.solve(stiffness_matrix, mass_matrix,
 *                 eigenvalues, eigenfunctions, eigenfunctions.size());
 * @endcode
 *
 * 为了支持这种使用情况，与PETSc包装器不同，该命名空间中的类的编写方式是在构造函数中初始化基础SLEPc对象。这样做也避免了对不同设置（如目标特征值或问题类型）的缓存；相反，这些设置在调用包装器类的相应函数时直接应用。
 * 上述方法的另一种实现方式是在应用程序代码中直接使用API的内部功能。这样一来，调用序列需要调用SolverBase的几个函数，而不是只有一个。这种自由是为了使用SLEPcWrappers，它需要对特征值问题求解器上下文有更多的处理。也请看API，例如。
 *
 * @code
 * template <typename OutputVector>
 * void
 * SolverBase::solve (const PETScWrappers::MatrixBase &A,
 *                  const PETScWrappers::MatrixBase &B,
 *                  std::vector<PetscScalar>        &eigenvalues,
 *                  std::vector<OutputVector>       &eigenvectors,
 *                  const unsigned int               n_eigenpairs)
 * {
 * ...
 * }
 * @endcode
 * 以此为例说明如何做到这一点。
 *
 *
 * @ingroup SLEPcWrappers
 *
 */
namespace SLEPcWrappers
{
  /**
   * 使用SLEPc求解器的求解器类的基类。由于SLEPc中的求解器是根据传递给通用求解器对象的标志来选择的，所以基本上所有实际的求解器调用都发生在这个类中，而派生类只是设置正确的标志来选择一个或另一个求解器，或者为单个求解器设置某些参数。
   * 关于如何使用这个类及其派生类的例子，包括如何为要计算特征值的矩阵提供预处理程序，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   *
   */
  class SolverBase
  {
  public:
    /**
     * 构造函数。接受MPI通信器，在该通信器上发生并行计算。
     *
     */
    SolverBase(SolverControl &cn, const MPI_Comm &mpi_communicator);

    /**
     * 解构器。
     *
     */
    virtual ~SolverBase();

    /**
     * 解决特征系统的复合方法  $Ax=\lambda x$
     * 。送入的特征向量必须至少有一个元素，我们可以在调整大小时作为模板使用，因为我们不知道使用的具体向量类的参数（即MPI向量的local_dofs）。然而，在复制特征向量时，至少要使用两倍于<tt>特征向量</tt>的内存大小（而且可能更多）。为了避免这样做，这里执行的相当标准的调用序列被使用。设置用于求解的矩阵；实际求解系统；收集求解结果。
     * @note
     * 注意收敛的特征向量的数量可能大于请求的特征向量的数量；这是由于特征问题求解器上下文的舍弃错误（成功）。如果发现这种情况，我们就不去管比请求的更多的特征对，而是通过忽略任何额外的东西来处理可能超过指定的数量。
     * 默认情况下，将计算一个特征向量/特征值对。
     *
     */
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          std::vector<PetscScalar> &       eigenvalues,
          std::vector<OutputVector> &      eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * 同上，但这里是解决系统的复合方法  $A x=\lambda B x$
     * ，用于实数矩阵、向量和值  $A, B, x, \lambda$  。
     *
     */
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          const PETScWrappers::MatrixBase &B,
          std::vector<PetscScalar> &       eigenvalues,
          std::vector<OutputVector> &      eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * 同上，但这里是解决系统 $A x=\lambda B x$ 与实数矩阵 $A,
     * B$ 和虚数特征对 $x, \lambda$ 的复合方法。
     *
     */
    template <typename OutputVector>
    void
    solve(const PETScWrappers::MatrixBase &A,
          const PETScWrappers::MatrixBase &B,
          std::vector<double> &            real_eigenvalues,
          std::vector<double> &            imag_eigenvalues,
          std::vector<OutputVector> &      real_eigenvectors,
          std::vector<OutputVector> &      imag_eigenvectors,
          const unsigned int               n_eigenpairs = 1);

    /**
     * 设置求解器的初始向量空间。
     * 默认情况下，SLEPc随机地初始化起始向量或初始子空间。
     *
     */
    template <typename Vector>
    void
    set_initial_space(const std::vector<Vector> &initial_space);

    /**
     * 设置要使用的谱系变换。
     *
     */
    void
    set_transformation(SLEPcWrappers::TransformationBase &this_transformation);

    /**
     * 设置要计算的频谱中的目标特征值。默认情况下，不设置目标。
     *
     */
    void
    set_target_eigenvalue(const PetscScalar &this_target);

    /**
     * 指明要计算的频谱的哪一部分。默认情况下，计算的是最大量级的特征值。
     * @note 关于其他允许的值，见SLEPc文档。
     *
     */
    void
    set_which_eigenpairs(EPSWhich set_which);

    /**
     * 指定特征谱问题的类型。这可以用来利用构成标准/广义同位素问题的矩阵的已知对称性。
     * 默认情况下，假定是一个非赫米特问题。
     * @note 关于其他允许的值，见SLEPc文档。
     *
     */
    void
    set_problem_type(EPSProblemType set_problem);

    /**
     * 利用SLEPc提供的信息，与deal.II自己的SolverControl对象进行核对，看是否已经达到收敛。
     *
     */
    void
    get_solver_state(const SolverControl::State state);

    /**
     * 异常。标准异常。
     *
     */
    DeclException0(ExcSLEPcWrappersUsageError);

    /**
     * 异常。带有错误号码的SLEPc错误。
     *
     */
    DeclException1(ExcSLEPcError,
                   int,
                   << "    An error with error number " << arg1
                   << " occurred while calling a SLEPc function");

    /**
     * 异常。在特征向量的数量上收敛失败。
     *
     */
    DeclException2(ExcSLEPcEigenvectorConvergenceMismatchError,
                   int,
                   int,
                   << "    The number of converged eigenvectors is " << arg1
                   << " but " << arg2 << " were requested. ");

    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;

  protected:
    /**
     * 对控制迭代求解器收敛的对象的引用。
     *
     */
    SolverControl &solver_control;

    /**
     * 将用于求解器的MPI通信器对象的副本。
     *
     */
    const MPI_Comm mpi_communicator;

    /**
     * 解决 <code>n_eigenpairs</code> 特征状态的线性系统。
     * 参数 <code>n_converged</code>
     * 包含已经收敛的特征态的实际数量；这可以少于或多于n_eigenpairs，取决于使用的SLEPc
     * eigensolver。
     *
     */
    void
    solve(const unsigned int n_eigenpairs, unsigned int *n_converged);

    /**
     * 访问已解决的特征向量问题的解决方案的实部，对索引解决方案，
     * $\text{index}\,\in\,0\dots \mathrm{n\_converged}-1$  。
     *
     */
    void
    get_eigenpair(const unsigned int         index,
                  PetscScalar &              eigenvalues,
                  PETScWrappers::VectorBase &eigenvectors);

    /**
     * 访问已解决的特征向量问题的解的实部和虚部，对索引解，
     * $\text{index}\,\in\,0\dots \mathrm{n\_converged}-1$  。
     *
     */
    void
    get_eigenpair(const unsigned int         index,
                  double &                   real_eigenvalues,
                  double &                   imag_eigenvalues,
                  PETScWrappers::VectorBase &real_eigenvectors,
                  PETScWrappers::VectorBase &imag_eigenvectors);

    /**
     * 初始化线性系统的求解器  $Ax=\lambda x$  .
     * (注意：在调用solve ()之前需要这样做)
     *
     */
    void
    set_matrices(const PETScWrappers::MatrixBase &A);

    /**
     * 同上，但在这里初始化线性系统的求解器  $A x=\lambda B
     * x$  。
     *
     */
    void
    set_matrices(const PETScWrappers::MatrixBase &A,
                 const PETScWrappers::MatrixBase &B);

  protected:
    /**
     * 特征值问题求解器的对象。
     *
     */
    EPS eps;

  private:
    /**
     * 收敛的原因。
     *
     */
    EPSConvergedReason reason;


    /**
     * 一个可以在SLEPc中作为回调使用的函数，以检查收敛情况。
     * @note 这个函数目前没有使用。
     *
     */
    static int
    convergence_test(EPS         eps,
                     PetscScalar real_eigenvalue,
                     PetscScalar imag_eigenvalue,
                     PetscReal   residual_norm,
                     PetscReal * estimated_error,
                     void *      solver_control);
  };



  /**
   * 一个使用SLEPc
   * Krylov-Schur求解器的求解器接口的实现。使用方法。所有光谱，所有问题类型，复杂。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverKrylovSchur : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {};

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在这个上下文中计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverKrylovSchur(SolverControl &       cn,
                      const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                      const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 一个使用SLEPc Arnoldi求解器的求解器接口的实现。
   * 使用方法。所有频谱，所有问题类型，复杂。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverArnoldi : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将延迟正交化的选项设置为false，即不做。
       *
       */
      AdditionalData(const bool delayed_reorthogonalization = false);

      /**
       * 延迟正交化的标志。
       *
       */
      bool delayed_reorthogonalization;
    };

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在这个上下文中计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverArnoldi(SolverControl &       cn,
                  const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                  const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 一个使用SLEPc Lanczos求解器的求解器接口的实现。
   * 使用方法。所有光谱，所有问题类型，复杂。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverLanczos : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {
      /**
       * 在Lanczos迭代过程中使用的正交化类型。
       *
       */
      EPSLanczosReorthogType reorthog;

      /**
       * 构造函数。默认情况下，在Lanczos迭代过程中使用的正交化类型为完全。
       *
       */
      AdditionalData(
        const EPSLanczosReorthogType r = EPS_LANCZOS_REORTHOG_FULL);
    };

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在该上下文上计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverLanczos(SolverControl &       cn,
                  const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                  const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 一个使用SLEPc Power求解器的求解器接口的实现。
   * 使用方法。仅限频谱的最大值，所有问题类型，复杂。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverPower : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，以便在需要时将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {};

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在这个上下文中计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverPower(SolverControl &       cn,
                const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用SLEPc
   * Davidson求解器的求解器接口的一个实现。用法。所有问题类型。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverGeneralizedDavidson : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {
      /**
       * 在搜索子空间中使用双重扩展。
       *
       */
      bool double_expansion;

      /**
       * 构造函数。默认将double_expansion设置为false。
       *
       */
      AdditionalData(bool double_expansion = false);
    };

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在该上下文上计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverGeneralizedDavidson(
      SolverControl &       cn,
      const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
      const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用SLEPc
   * Jacobi-Davidson求解器的求解器接口的一个实现。用法。所有问题类型。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档，以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverJacobiDavidson : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {};

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在这个上下文中计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverJacobiDavidson(SolverControl & cn,
                         const MPI_Comm &mpi_communicator = PETSC_COMM_SELF,
                         const AdditionalData &data       = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 一个使用SLEPc LAPACK直接求解器的求解器接口的实现。
   * 关于如何使用这个类和它的兄弟类的例子，包括如何为要计算特征值的矩阵提供前置条件，请参见SolverBase类的文档以及SLEPcWrappers命名空间的文档中的广泛讨论。
   * @ingroup SLEPcWrappers
   *
   */
  class SolverLAPACK : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，如果需要的话，可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {};

    /**
     * SLEPc求解器将希望有一个MPI通信器上下文，在这个上下文中计算被并行化。默认情况下，这与PETScWrappers的行为相同，但你可以改变它。
     *
     */
    SolverLAPACK(SolverControl &       cn,
                 const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                 const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  // --------------------------- inline and template functions -----------
  /**
   * 在此声明的目的是为了使其有可能采取 std::vector
   * 不同的PETScWrappers向量类型
   *
   */
  // todo: The logic of these functions can be simplified without breaking
  // backward compatibility...

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    std::vector<PetscScalar> &       eigenvalues,
                    std::vector<OutputVector> &      eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged > n_eigenpairs)
      n_converged = n_eigenpairs;
    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));

    AssertThrow(eigenvectors.size() != 0, ExcSLEPcWrappersUsageError());
    eigenvectors.resize(n_converged, eigenvectors.front());
    eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
  }

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    const PETScWrappers::MatrixBase &B,
                    std::vector<PetscScalar> &       eigenvalues,
                    std::vector<OutputVector> &      eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Guard against incompatible matrix sizes:
    AssertThrow(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
    AssertThrow(A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));

    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A, B);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged >= n_eigenpairs)
      n_converged = n_eigenpairs;

    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));
    AssertThrow(eigenvectors.size() != 0, ExcSLEPcWrappersUsageError());

    eigenvectors.resize(n_converged, eigenvectors.front());
    eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index, eigenvalues[index], eigenvectors[index]);
  }

  template <typename OutputVector>
  void
  SolverBase::solve(const PETScWrappers::MatrixBase &A,
                    const PETScWrappers::MatrixBase &B,
                    std::vector<double> &            real_eigenvalues,
                    std::vector<double> &            imag_eigenvalues,
                    std::vector<OutputVector> &      real_eigenvectors,
                    std::vector<OutputVector> &      imag_eigenvectors,
                    const unsigned int               n_eigenpairs)
  {
    // Guard against incompatible matrix sizes:
    AssertThrow(A.m() == B.m(), ExcDimensionMismatch(A.m(), B.m()));
    AssertThrow(A.n() == B.n(), ExcDimensionMismatch(A.n(), B.n()));

    // and incompatible eigenvalue/eigenvector sizes
    AssertThrow(real_eigenvalues.size() == imag_eigenvalues.size(),
                ExcDimensionMismatch(real_eigenvalues.size(),
                                     imag_eigenvalues.size()));
    AssertThrow(real_eigenvectors.size() == imag_eigenvectors.size(),
                ExcDimensionMismatch(real_eigenvectors.size(),
                                     imag_eigenvectors.size()));

    // Panic if the number of eigenpairs wanted is out of bounds.
    AssertThrow((n_eigenpairs > 0) && (n_eigenpairs <= A.m()),
                ExcSLEPcWrappersUsageError());

    // Set the matrices of the problem
    set_matrices(A, B);

    // and solve
    unsigned int n_converged = 0;
    solve(n_eigenpairs, &n_converged);

    if (n_converged >= n_eigenpairs)
      n_converged = n_eigenpairs;

    AssertThrow(n_converged == n_eigenpairs,
                ExcSLEPcEigenvectorConvergenceMismatchError(n_converged,
                                                            n_eigenpairs));
    AssertThrow((real_eigenvectors.size() != 0) &&
                  (imag_eigenvectors.size() != 0),
                ExcSLEPcWrappersUsageError());

    real_eigenvectors.resize(n_converged, real_eigenvectors.front());
    imag_eigenvectors.resize(n_converged, imag_eigenvectors.front());
    real_eigenvalues.resize(n_converged);
    imag_eigenvalues.resize(n_converged);

    for (unsigned int index = 0; index < n_converged; ++index)
      get_eigenpair(index,
                    real_eigenvalues[index],
                    imag_eigenvalues[index],
                    real_eigenvectors[index],
                    imag_eigenvectors[index]);
  }

  template <typename Vector>
  void
  SolverBase::set_initial_space(const std::vector<Vector> &this_initial_space)
  {
    std::vector<Vec> vecs(this_initial_space.size());

    for (unsigned int i = 0; i < this_initial_space.size(); i++)
      {
        Assert(this_initial_space[i].l2_norm() > 0.0,
               ExcMessage("Initial vectors should be nonzero."));
        vecs[i] = this_initial_space[i];
      }

    // if the eigensolver supports only a single initial vector, but several
    // guesses are provided, then all except the first one will be discarded.
    // One could still build a vector that is rich in the directions of all
    // guesses, by taking a linear combination of them. (TODO: make function
    // virtual?)

    const PetscErrorCode ierr =
      EPSSetInitialSpace(eps, vecs.size(), vecs.data());
    AssertThrow(ierr == 0, ExcSLEPcError(ierr));
  }

} // namespace SLEPcWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_SLEPC

 /*----------------------------   slepc_solver.h  ---------------------------*/ 

#endif

 /*----------------------------   slepc_solver.h  ---------------------------*/ 


