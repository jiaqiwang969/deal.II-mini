//include/deal.II-translator/lac/petsc_solver_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2021 by the deal.II authors
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

#ifndef dealii_petsc_solver_h
#  define dealii_petsc_solver_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/solver_control.h>

#    include <petscksp.h>

#    include <memory>

#    ifdef DEAL_II_WITH_SLEPC
#      include <deal.II/lac/slepc_spectral_transformation.h>
#    endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#    ifndef DOXYGEN
#      ifdef DEAL_II_WITH_SLEPC
namespace SLEPcWrappers
{
  // forward declarations
  class TransformationBase;
} // namespace SLEPcWrappers
#      endif
#    endif

namespace PETScWrappers
{
  // forward declarations
#    ifndef DOXYGEN
  class MatrixBase;
  class VectorBase;
  class PreconditionBase;
#    endif


  /**
   * 使用PETSc求解器的求解器类的基类。由于PETSc中的求解器是根据传递给通用求解器对象的标志来选择的，所以基本上所有实际的求解器调用都发生在这个类中，派生类只是设置正确的标志来选择一个或另一个求解器，或者为单个求解器设置某些参数。
   * 作为选择，用户可以创建一个从SolverBase类派生出来的求解器，并且可以通过SolverControl来设置解决线性方程组所需的默认参数。这些默认选项可以通过指定形式为
   * @p 的命令行参数而被覆盖。
   *
   * - sp_*。  例如， @p  * * - sp_*。
   *
   * - sp_monitor_true_residual在每次迭代时打印出真正的残差准则（无条件的）， @p  * - sp_monitor_true_residual在每次迭代时打印出真正的残差准则。
   *
   * - sp_view提供了关于当前上下文中使用的线性求解器和预处理器的信息。解算器的类型也可以在运行时通过指定 @p 来改变。
   *
   * - sp_type{richardson, cg, gmres, fgmres, ...}来动态测试最佳求解器以及使用 @p 的合适的前置条件集。
   *
   * - c_type{jacobi, bjacobi, ilu, lu, ..}。还有其他几个命令行选项可用于修改PETSc线性求解器的行为，可从<a
   * href="http://www.mcs.anl.gov/petsc">documentation and manual
   * pages</a>中获得。
   * @note
   * 在带有预设条件器的求解器对象上重复调用solve()必须谨慎使用。先决条件在第一次调用solve()时被初始化，随后的调用会重复使用求解器和先决条件对象。这样做是出于性能考虑。可以通过调用reset()来重置求解器和前置条件器。
   * PETSc的一个问题是
   *
   * - 特别是在MPI模式下
   *
   * 它经常不产生非常有用的错误信息。为了节省其他用户寻找一个难以追踪的错误的时间，这里有一种情况和得到的错误信息：当你没有给你的求解器的构造器指定一个MPI通信器时。在这种情况下，你会从你的每个并行进程中得到以下形式的错误。
   * @verbatim
   * [1]PETSC ERROR: PCSetVector() line 1173 in src/ksp/pc/interface/precon.c
   * [1]PETSC ERROR:   Arguments must have same communicators!
   * [1]PETSC ERROR:   Different communicators in the two objects: Argument #
   * 1 and 2! [1]PETSC ERROR: KSPSetUp() line 195 in
   * src/ksp/ksp/interface/itfunc.c
   * @endverbatim
   * 这个错误是由于没有指定MPI通信器而导致的，人们可以花很长的时间来弄清楚到底出了什么问题。请注意，通信器
   * @em
   * 必须与我们要解决的线性系统中的矩阵和所有向量相匹配。使情况恶化的是，解算器类的默认参数，
   * @p
   * PETSC_COMM_SELF，是顺序情况下的适当参数（这就是为什么它是默认参数），所以这个错误只在并行模式下显示出来。
   * @ingroup PETScWrappers
   *
   */
  class SolverBase
  {
  public:
    /**
     * 构造函数。接受求解器控制对象和MPI通信器，并行计算将在该通信器上进行。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，以便用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverBase(SolverControl &cn, const MPI_Comm &mpi_communicator);

    /**
     * 解构器。
     *
     */
    virtual ~SolverBase() = default;

    /**
     * 解决线性系统<tt>Ax=b</tt>。根据派生类提供的信息和作为预处理程序传递的对象，会选择PETSc的一个线性求解器和预处理程序。
     * 由于性能的原因，重复调用solve()不会重构预处理程序。参见类的文档。
     *
     */
    void
    solve(const MatrixBase &      A,
          VectorBase &            x,
          const VectorBase &      b,
          const PreconditionBase &preconditioner);


    /**
     * 重置所包含的前置条件器和求解器对象。更多细节见类描述。
     *
     */
    virtual void
    reset();


    /**
     * 为求解器对象设置一个前缀名称。在用命令行选项定制PETSc
     * KSP对象时很有用。
     *
     */
    void
    set_prefix(const std::string &prefix);


    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;

    /**
     * 用前置条件器初始化求解器。该函数旨在与SLEPc谱系变换类一起使用。
     *
     */
    void
    initialize(const PreconditionBase &preconditioner);

  protected:
    /**
     * 对控制迭代求解器收敛的对象的引用。事实上，对于这些PETSc封装器来说，PETSc本身就这样做了，但是我们在开始求解过程之前从这个对象中复制数据，之后再把数据复制回它。
     *
     */
    SolverControl &solver_control;

    /**
     * 将用于求解器的MPI通信器对象的拷贝。
     *
     */
    const MPI_Comm mpi_communicator;

    /**
     * %函数，它接收一个Krylov Subspace
     * Solver上下文对象，并设置派生类所要求的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const = 0;

    /**
     * 解算器前缀名称，用于限定当前上下文中PETSc
     * KSP对象的特定选项。注意：前缀名称的开头不能有连字符（-）。所有运行时选项的第一个字符是自动的连字符。
     *
     */
    std::string prefix_name;

  private:
    /**
     * 一个在PETSc中用作回调的函数，用于检查收敛情况。
     * 它采用PETSc提供的信息，并与deal.II自己的SolverControl对象进行核对，以确定是否已经达到收敛。
     *
     */
    static PetscErrorCode
    convergence_test(KSP                 ksp,
                     const PetscInt      iteration,
                     const PetscReal     residual_norm,
                     KSPConvergedReason *reason,
                     void *              solver_control);

    /**
     * 一个包含PETSc求解器和预处理器对象的结构。
     * 如果在前一个求解步骤中使用了相同的预处理程序，那么这个对象在后续的求解器调用中会被保留下来。如果设置先决条件的成本很高，例如在ILU的情况下，这可以节省一些计算时间。
     * 这个类的实际声明很复杂，因为PETSc在2.1.6和2.2.0版本之间完全不兼容地改变了它的求解器接口
     * :-(
     * 这种类型的对象是明确创建的，但是当周围的求解器对象超出范围时，或者当我们给这个对象的指针分配一个新值时，就会被销毁。因此，各自的Destroy函数被写入该对象的析构器中，即使该对象没有构造函数。
     *
     */
    struct SolverData
    {
      /**
       * 毁灭函数
       *
       */
      ~SolverData();

      /**
       * Krylov子空间求解器的对象。
       *
       */
      KSP ksp;
    };

    /**
     * 指向存储求解器上下文的对象的指针。如果有必要，这将在主求解器例程中重新创建。
     *
     */
    std::unique_ptr<SolverData> solver_data;

#    ifdef DEAL_II_WITH_SLEPC
    // Make the transformation class a friend, since it needs to set the KSP
    // solver.
    friend class SLEPcWrappers::TransformationBase;
#    endif
  };



  /**
   * 使用PETSc Richardson解算器的解算器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverRichardson : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将阻尼参数设置为1。
       *
       */
      explicit AdditionalData(const double omega = 1);

      /**
       * 松弛参数。
       *
       */
      double omega;
    };

    /**
     * 构造函数。与deal.II自己的求解器相比，不需要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverRichardson(SolverControl &       cn,
                     const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                     const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 使用PETSc
   * Chebyshev（或3.3版之前的Chebychev）求解器接口的一个实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverChebychev : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverChebychev(SolverControl &       cn,
                    const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                    const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 一个使用PETSc CG求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverCG(SolverControl &       cn,
             const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
             const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 一个使用PETSc BiCG求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverBiCG : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverBiCG(SolverControl &       cn,
               const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
               const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 使用PETSc GMRES求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，设置临时向量的数量为30，即每30次迭代做一次重启。
       *
       */
      AdditionalData(const unsigned int restart_parameter     = 30,
                     const bool         right_preconditioning = false);

      /**
       * 临时向量的最大数量。
       *
       */
      unsigned int restart_parameter;

      /**
       * 右侧预处理的标志。
       *
       */
      bool right_preconditioning;
    };

    /**
     * 构造函数。与deal.II自己的求解器相比，不需要给出一个向量存储对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverGMRES(SolverControl &       cn,
                const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 使用PETSc BiCGStab求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverBicgstab(SolverControl &       cn,
                   const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                   const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };

  /**
   * 一个使用PETSc CG Squared求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverCGS(SolverControl &       cn,
              const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
              const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 使用PETSc TFQMR求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverTFQMR(SolverControl &       cn,
                const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 一个使用PETSc
   * TFQMR-2求解器（在PETSc中称为TCQMR）的求解器接口的实现。请注意，这个求解器在PETSc
   * 2.1.6以下的版本中有一个严重的错误，即它不检查收敛性，总是返回一个错误代码。因此，在PETSc
   * 2.1.6及以前的版本中，这个类会以一个错误的方式终止，表明未能收敛。不过这个问题应该在以后的PETSc版本中得到修正。
   * @ingroup PETScWrappers
   *
   */
  class SolverTCQMR : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，可以向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给出一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverTCQMR(SolverControl &       cn,
                const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 一个使用PETSc CR求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverCR : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverCR(SolverControl &       cn,
             const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
             const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 为这个特定的求解器存储一份标志。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };



  /**
   * 一个使用PETSc最小二乘法求解器的求解器接口的实现。
   * @ingroup PETScWrappers
   *
   */
  class SolverLSQR : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给出一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverLSQR(SolverControl &       cn,
               const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
               const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };


  /**
   * 一个使用PETSc PREONLY求解器的求解器接口的实现。
   * 实际上这不是一个真正的求解算法。solve()只应用一次预处理程序，并立即返回。它的唯一目的是提供一个求解器对象，当预设条件器应该作为一个真正的求解器使用时。它与完整的LU分解预处理<tt>
   * PreconditionLU
   * </tt>结合起来非常有用，它与这个求解器类一起成为一个直接求解器。
   * @ingroup PETScWrappers
   *
   */
  class SolverPreOnly : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造器。与deal.II自己的求解器相比，没有必要给一个矢量内存对象。然而，PETSc求解器希望有一个MPI通信器的上下文，在这个上下文上计算被并行化。默认情况下，这里使用
     * @p PETSC_COMM_SELF
     * ，但你可以改变这个。注意，对于单处理器（非MPI）版本，这个参数没有任何影响。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的调谐标志。
     * 请注意，这里使用的通信器必须与系统矩阵、解决方案和右侧对象中使用的通信器相匹配，才能用这个求解器来完成。否则，PETSc将产生难以追踪的错误，见SolverBase类的文档。
     *
     */
    SolverPreOnly(SolverControl &       cn,
                  const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                  const AdditionalData &data             = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;

    /**
     * %函数，接收一个Krylov Subspace
     * Solver上下文对象，并设置适合该类的求解器类型。
     *
     */
    virtual void
    set_solver_type(KSP &ksp) const override;
  };

  /**
   * 一个通过PETSc使用稀疏直接MUMPS求解器接口的实现。这个类具有所有其他求解器类的常规接口，但它当然不同，因为它没有实现迭代求解器。因此，像SolverControl对象这样的东西在这里没有特别的意义。
   * MUMPS允许在这个矩阵中使用对称性。在这个类中，这是由set_symmetric_mode()函数实现的。如果你的矩阵是对称的，你可以按以下方式使用这个类。
   * @code
   * SolverControl cn;
   * PETScWrappers::SparseDirectMUMPS solver(cn, mpi_communicator);
   * solver.set_symmetric_mode(true);
   * solver.solve(system_matrix, solution, system_rhs);
   * @endcode
   *
   * @note
   * 该类内部调用KSPSetFromOptions，因此你能够使用MATSOLVERMUMPS包的所有PETSc参数。见http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MATSOLVERMUMPS.html
   * @ingroup PETScWrappers
   *
   */
  class SparseDirectMUMPS : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构可以将额外的数据输送到求解器中。
     *
     */
    struct AdditionalData
    {};
    /**
     * 构造函数
     *
     */
    SparseDirectMUMPS(SolverControl &       cn,
                      const MPI_Comm &      mpi_communicator = PETSC_COMM_SELF,
                      const AdditionalData &data = AdditionalData());

    /**
     * 解决线性系统的方法。
     *
     */
    void
    solve(const MatrixBase &A, VectorBase &x, const VectorBase &b);

    /**
     * 如果系统矩阵是对称的，该方法允许使用LDL^T分解，而不是更昂贵的LU分解。该参数表示矩阵是否对称。
     *
     */
    void
    set_symmetric_mode(const bool flag);

  protected:
    /**
     * 为这个特定的求解器存储一份标志的副本。
     *
     */
    const AdditionalData additional_data;

    virtual void
    set_solver_type(KSP &ksp) const override;

  private:
    /**
     * 一个在PETSc中用作回调的函数，用于检查收敛性。它接收来自PETSc的信息，并与deal.II自己的SolverControl对象进行核对，看是否已经达到收敛。
     *
     */
    static PetscErrorCode
    convergence_test(KSP                 ksp,
                     const PetscInt      iteration,
                     const PetscReal     residual_norm,
                     KSPConvergedReason *reason,
                     void *              solver_control);

    /**
     * 一个包含PETSc求解器和预处理器对象的结构。
     * 由于这里没有使用基础中的solve成员函数，所以位于基础中的私有SolverData结构也不能使用。
     *
     */
    struct SolverDataMUMPS
    {
      /**
       * 破坏器
       *
       */
      ~SolverDataMUMPS();

      KSP ksp;
      PC  pc;
    };

    std::unique_ptr<SolverDataMUMPS> solver_data;

    /**
     * 标志指定被分解的矩阵是否是对称的。它影响到使用的预处理程序的类型（PCLU或PCCHOLESKY）。
     *
     */
    bool symmetric_mode;
  };
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

 /*----------------------------   petsc_solver.h ---------------------------*/ 

#endif
 /*----------------------------   petsc_solver.h ---------------------------*/ 


