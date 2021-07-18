//include/deal.II-translator/lac/trilinos_solver_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_trilinos_solver_h
#  define dealii_trilinos_solver_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_TRILINOS

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/la_parallel_vector.h>
#    include <deal.II/lac/solver_control.h>
#    include <deal.II/lac/vector.h>

#    include <Amesos.h>
#    include <AztecOO.h>
#    include <Epetra_LinearProblem.h>
#    include <Epetra_Operator.h>

#    include <memory>


DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // forward declarations
#    ifndef DOXYGEN
  class SparseMatrix;
  class PreconditionBase;
#    endif


  /**
   * 使用Trilinos求解器的求解器类的基类。由于Trilinos中的求解器是根据传递给通用求解器对象的标志来选择的，基本上所有实际的求解器调用都发生在这个类中，派生类只是设置正确的标志来选择一个或另一个求解器，或者为单个求解器设置某些参数。关于Trilinos求解器包AztecOO的一般讨论，我们参考了<a
   * href="https://trilinos.org/docs/dev/packages/aztecoo/doc/html/index.html">AztecOO
   * user guide</a>。
   * 这个求解器类也可以作为一个独立的类来使用，通过标志<tt>solver_name</tt>来设置各自的Krylov方法。这可以在运行时进行（例如，在从ParameterList解析求解器时），与deal.II类SolverSelector类似。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverBase
  {
  public:
    /**
     * 枚举对象，在派生类的构造函数中设置，并告诉Trilinos要使用哪个求解器。这个选项也可以在用户程序中设置，所以当解算器应该在运行时设置时，人们可能会使用这个基类而不是某个专门的派生类。目前启用的选项有。
     *
     */
    enum SolverName
    {
      /**
       * 使用共轭梯度（CG）算法。
       *
       */
      cg,
      /**
       * 使用共轭梯度平方（CGS）算法。
       *
       */
      cgs,
      /**
       * 使用广义最小残差（GMRES）算法。
       *
       */
      gmres,
      /**
       * 使用双共轭梯度稳定化（BICGStab）算法。
       *
       */
      bicgstab,
      /**
       * 使用无转置的准最小残差（TFQMR）方法。
       *
       */
      tfqmr
    } solver_name;

    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */

    struct AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式，并在派生类为GMRES的情况下放入重启参数。
       * TODO:
       * 找到一个更好的方法来设置GMRES重启参数，因为在基类中为所有求解器设置一个求解器的特定选项是相当不优雅的。
       *
       */
      explicit AdditionalData(const bool         output_solver_details = false,
                              const unsigned int gmres_restart_parameter = 30);

      /**
       * 启用/禁用求解器细节的输出（每次迭代的残差等）。
       *
       */
      const bool output_solver_details;

      /**
       * GMRES求解器的重新启动参数。
       *
       */
      const unsigned int gmres_restart_parameter;
    };

    /**
     * 构造函数。接受求解器控制对象并创建求解器。
     *
     */
    SolverBase(SolverControl &       cn,
               const AdditionalData &data = AdditionalData());

    /**
     * 第二个构造函数。这个构造函数接收一个指定求解器名称的枚举对象，并设置适当的Krylov方法。
     *
     */
    SolverBase(const enum SolverName solver_name,
               SolverControl &       cn,
               const AdditionalData &data = AdditionalData());

    /**
     * 解构器。
     *
     */
    virtual ~SolverBase() = default;

    /**
     * 解决线性系统<tt>Ax=b</tt>。根据派生类提供的信息和作为预处理程序传递的对象，选择Trilinos的线性求解器和预处理程序中的一个。
     *
     */
    void
    solve(const SparseMatrix &    A,
          MPI::Vector &           x,
          const MPI::Vector &     b,
          const PreconditionBase &preconditioner);

    /**
     * 解决线性系统<tt>Ax=b</tt>，其中<tt>A</tt>是一个算子。
     * 这个函数可以用来进行无矩阵计算。根据派生类提供的信息和作为预处理程序传递的对象，选择Trilinos的线性求解器和预处理程序中的一个。
     *
     */
    void
    solve(const Epetra_Operator & A,
          MPI::Vector &           x,
          const MPI::Vector &     b,
          const PreconditionBase &preconditioner);

    /**
     * 求解线性系统<tt>Ax=b</tt>，其中<tt>A</tt>和其 @p
     * preconditioner 都是一个算子。    当<tt>A</tt>和 @p
     * preconditioner
     * 都是由TrilinosPayload派生的LinearOperators时，就可以使用这个函数。
     * 根据派生类提供的信息和作为预处理程序传递的对象，会选择Trilinos的线性求解器和预处理程序之一。
     *
     */
    void
    solve(const Epetra_Operator &A,
          MPI::Vector &          x,
          const MPI::Vector &    b,
          const Epetra_Operator &preconditioner);

    /**
     * 求解线性系统<tt>Ax=b</tt>，其中<tt>A</tt>是一个算子，向量
     * @p x 和 @p b 是Trilinos本地的向量类型。
     * 当<tt>A</tt>是一个从TrilinosPayload派生的LinearOperators时，可以使用这个函数。
     * 根据派生类提供的信息和作为预处理程序传递的对象，会选择Trilinos的线性求解器和预处理程序之一。
     *
     */
    void
    solve(const Epetra_Operator &   A,
          Epetra_MultiVector &      x,
          const Epetra_MultiVector &b,
          const PreconditionBase &  preconditioner);

    /**
     * 求解线性系统<tt>Ax=b</tt>，其中<tt>A</tt>及其 @p
     * preconditioner 都是一个算子，而向量 @p x 和 @p b
     * 是Trilinos的本地向量类型。    当<tt>A</tt>和 @p
     * preconditioner
     * 都是源自TrilinosPayload的LinearOperators时，可以使用这个函数。
     * 根据派生类提供的信息和作为预处理程序传递的对象，会选择Trilinos的线性求解器和预处理程序之一。
     *
     */
    void
    solve(const Epetra_Operator &   A,
          Epetra_MultiVector &      x,
          const Epetra_MultiVector &b,
          const Epetra_Operator &   preconditioner);



    /**
     * 解决线性系统<tt>Ax=b</tt>。根据派生类提供的信息和作为预处理程序传递的对象，选择Trilinos的线性求解器和预处理程序中的一个。
     * 该类按照TrilinosWrappers的格式处理矩阵，但可以将deal.II向量作为参数。由于deal.II是串行向量（非分布式），这个函数只在矩阵为本地所有的情况下做你期望的事情。否则，会产生一个异常。
     *
     */
    void
    solve(const SparseMatrix &          A,
          dealii::Vector<double> &      x,
          const dealii::Vector<double> &b,
          const PreconditionBase &      preconditioner);

    /**
     * 解决线性系统<tt>Ax=b</tt>，其中<tt>A</tt>是一个运算符。
     * 这个函数可用于无矩阵计算。根据派生类提供的信息和作为预处理程序传递的对象，会选择Trilinos的线性求解器和预处理程序之一。该类根据TrilinosWrappers的格式处理矩阵，但可以将deal.II向量作为参数。
     * 由于deal.II是串行向量（非分布式），这个函数只在矩阵为本地所有的情况下做你期望的事情。否则，会产生一个异常。
     *
     */
    void
    solve(Epetra_Operator &             A,
          dealii::Vector<double> &      x,
          const dealii::Vector<double> &b,
          const PreconditionBase &      preconditioner);

    /**
     * 为deal.II的平行分布向量求解线性系统<tt>Ax=b</tt>。根据派生类提供的信息和作为预处理器传递的对象，选择Trilinos的线性求解器和预处理器中的一个。
     *
     */
    void
    solve(const SparseMatrix &                                      A,
          dealii::LinearAlgebra::distributed::Vector<double> &      x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b,
          const PreconditionBase &preconditioner);

    /**
     * 解决线性系统<tt>Ax=b</tt>，其中<tt>A</tt>是一个算子。
     * 这个函数可以用来进行无矩阵计算。根据派生类提供的信息和作为预处理程序传递的对象，选择Trilinos的线性求解器和预处理程序中的一个。
     *
     */
    void
    solve(Epetra_Operator &                                         A,
          dealii::LinearAlgebra::distributed::Vector<double> &      x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b,
          const PreconditionBase &preconditioner);


    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;

    /**
     * 异常情况
     *
     */
    DeclException1(ExcTrilinosError,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while calling a Trilinos function");

  protected:
    /**
     * 对控制迭代求解器收敛性的对象的引用。事实上，对于这些Trilinos包装器来说，Trilinos本身就这样做了，但是我们在开始求解过程之前从这个对象中复制数据，之后再把数据复制回这个对象中。
     *
     */
    SolverControl &solver_control;

  private:
    /**
     * 解决函数用于正确设置Epetra_LinearProblem，一旦完成，这个函数就能解决线性问题。
     *
     */
    template <typename Preconditioner>
    void
    do_solve(const Preconditioner &preconditioner);

    /**
     * 一个设置解算器将应用的预处理程序的函数。
     *
     */
    template <typename Preconditioner>
    void
    set_preconditioner(AztecOO &solver, const Preconditioner &preconditioner);

    /**
     * 一个收集Trilinos稀疏矩阵、右手向量和解向量的结构，它被传递给Trilinos求解器。
     *
     */
    std::unique_ptr<Epetra_LinearProblem> linear_problem;

    /**
     * 一个包含Trilinos对象的结构，可以查询线性求解器并确定是否达到收敛标准。
     *
     */
    std::unique_ptr<AztecOO_StatusTest> status_test;

    /**
     * 一个包含Trilinos求解器和预处理器对象的结构。
     *
     */
    AztecOO solver;

    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };


  // provide a declaration for two explicit specializations
  template <>
  void
  SolverBase::set_preconditioner(AztecOO &               solver,
                                 const PreconditionBase &preconditioner);

  template <>
  void
  SolverBase::set_preconditioner(AztecOO &              solver,
                                 const Epetra_Operator &preconditioner);



  /**
   * 使用Trilinos CG求解器的求解器接口的实现。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverCG : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */

    struct AdditionalData : public SolverBase::AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式。
       *
       */
      explicit AdditionalData(const bool output_solver_details = false);
    };

    /**
     * 构造器。与deal.II自己的求解器相比，不需要给出一个矢量内存对象。
     * 最后一个参数是一个结构，带有额外的、与求解器相关的标志，用于调整。
     *
     */
    SolverCG(SolverControl &cn, const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用Trilinos CGS求解器的求解器接口的实现。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverCGS : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData : public SolverBase::AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式。
       *
       */
      explicit AdditionalData(const bool output_solver_details = false);
    };

    /**
     * 构造器。与deal.II自己的求解器相比，不需要给出一个矢量内存对象。
     * 最后一个参数是一个结构，带有额外的、与求解器相关的标志，用于调整。
     *
     */
    SolverCGS(SolverControl &cn, const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用Trilinos GMRES求解器的求解器接口的实现。
   *
   */
  class SolverGMRES : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData : public SolverBase::AdditionalData
    {
      /**
       * 构造函数。默认情况下，设置临时向量的数量为30，即每30次迭代做一次重启。
       *
       */
      explicit AdditionalData(const bool         output_solver_details = false,
                              const unsigned int restart_parameter     = 30);
    };

    /**
     * 构造函数。与deal.II自己的求解器相比，不需要给出一个向量存储对象。
     * 最后一个参数是一个结构，包含额外的、与求解器相关的标志，用于调谐。
     *
     */
    SolverGMRES(SolverControl &       cn,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用Trilinos BiCGStab求解器的求解器接口的实现。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverBicgstab : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData : public SolverBase::AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式。
       *
       */
      explicit AdditionalData(const bool output_solver_details = false);
    };

    /**
     * 构造器。与deal.II自己的求解器相比，不需要给出一个矢量内存对象。
     * 最后一个参数是一个结构，带有额外的、与求解器相关的标志，用于调整。
     *
     */
    SolverBicgstab(SolverControl &       cn,
                   const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * 使用Trilinos TFQMR求解器的求解器接口的实现。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverTFQMR : public SolverBase
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */
    struct AdditionalData : public SolverBase::AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式。
       *
       */
      explicit AdditionalData(const bool output_solver_details = false);
    };

    /**
     * 构造器。与deal.II自己的求解器相比，不需要给出一个矢量内存对象。
     * 最后一个参数是一个结构，带有额外的、与求解器相关的标志，用于调整。
     *
     */
    SolverTFQMR(SolverControl &       cn,
                const AdditionalData &data = AdditionalData());

  protected:
    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };



  /**
   * Trilinos直接求解器的一个实现（使用Amesos包）。  数据域
   * AdditionalData::solver_type
   * 可以用来指定求解器的类型。它允许使用内置解算器Amesos_Klu以及第三方解算器Amesos_Superludist或Amesos_Mumps。
   * 关于如何安装Trilinos以使用KLU以外的直接求解器的说明，请参见deal.II
   * ReadMe文件中链接的Trilinos安装说明。
   * @ingroup TrilinosWrappers
   *
   */
  class SolverDirect
  {
  public:
    /**
     * 标准化的数据结构，用于向求解器输送额外的数据。
     *
     */

    struct AdditionalData
    {
      /**
       * 将附加数据字段设置为所需的输出格式。
       *
       */
      explicit AdditionalData(const bool         output_solver_details = false,
                              const std::string &solver_type = "Amesos_Klu");

      /**
       * 启用/禁用求解器细节的输出（每次迭代的残差等）。
       *
       */
      bool output_solver_details;

      /**
       * 设置求解器类型（用于支持Triminos Amesos软件包的第三方求解器）。可能的情况是。        <ul>   <li>  "Amesos_Lapack"  </li>   <li>  "Amesos_Scalapack"  </li>   <li>  "Amesos_Klu"  </li>   <li>  ] "Amesos_Umfpack"  </li>   <li>  "Amesos_Pardiso"  </li>   <li>  "Amesos_Taucs"  </li>   <li>  "Amesos_Superlu"  </li>  ]  <li>  "Amesos_Superludist"  </li>   <li>  "Amesos_Dscpack"  </li>   <li>  "Amesos_Mumps"  </li>   </ul>  注意，这些求解器在 deal.II 的可用性取决于配置 Trilinos 时设置了哪些求解器。
       *
       */
      std::string solver_type;
    };

    /**
     * 构造函数。接受求解器控制对象并创建求解器。
     *
     */
    SolverDirect(SolverControl &       cn,
                 const AdditionalData &data = AdditionalData());

    /**
     * 解构器。
     *
     */
    virtual ~SolverDirect() = default;

    /**
     * 初始化矩阵<tt>A</tt>的直接求解器，并用从附加数据结构中选择的包为它创建一个因式分解。请注意，这里不需要预处理程序，也不调用solve()。
     *
     */
    void
    initialize(const SparseMatrix &A);

    /**
     * 根据initialize()中设置的包，解决线性系统<tt>Ax=b</tt>。注意在这个调用过程中，矩阵没有被重构。
     *
     */
    void
    solve(MPI::Vector &x, const MPI::Vector &b);

    /**
     * 根据initialize()中设置的软件包，为deal.II自己的平行向量求解线性系统<tt>Ax=b</tt>。注意在这个调用过程中，矩阵没有被重构。
     *
     */
    void
    solve(dealii::LinearAlgebra::distributed::Vector<double> &      x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b);

    /**
     * 解决线性系统<tt>Ax=b</tt>。用从附加数据结构中选择的包创建一个矩阵的因式分解，并执行求解。注意，这里不需要预处理程序。
     *
     */
    void
    solve(const SparseMatrix &A, MPI::Vector &x, const MPI::Vector &b);

    /**
     * 解决线性系统<tt>Ax=b</tt>。这个类对Trilinos矩阵起作用，但需要deal.II串行向量作为参数。由于这些向量不是分布式的，这个函数只在矩阵是串行的情况下（即本地拥有）做你所期望的事情。否则，将抛出一个异常。
     *
     */
    void
    solve(const SparseMatrix &          A,
          dealii::Vector<double> &      x,
          const dealii::Vector<double> &b);

    /**
     * 为deal.II自己的并行向量求解线性系统<tt>Ax=b</tt>。用从附加数据结构中选择的包创建一个矩阵的因式分解，并执行求解。注意，这里不需要预处理程序。
     *
     */
    void
    solve(const SparseMatrix &                                      A,
          dealii::LinearAlgebra::distributed::Vector<double> &      x,
          const dealii::LinearAlgebra::distributed::Vector<double> &b);

    /**
     * 访问控制收敛的对象。
     *
     */
    SolverControl &
    control() const;

    /**
     * 异常情况
     *
     */
    DeclException1(ExcTrilinosError,
                   int,
                   << "An error with error number " << arg1
                   << " occurred while calling a Trilinos function");

  private:
    /**
     * 实际执行解决线性系统的操作，包括因式分解和前向、后向替换。
     *
     */
    void
    do_solve();

    /**
     * 对控制迭代求解器收敛的对象的引用。事实上，对于这些Trilinos包装器来说，Trilinos本身就是这样做的，但是我们在开始求解过程之前从这个对象中复制数据，之后再将数据复制回对象中。
     *
     */
    SolverControl &solver_control;

    /**
     * 一个收集Trilinos稀疏矩阵、右手边向量和求解向量的结构，它被传递给Trilinos求解器。
     *
     */
    std::unique_ptr<Epetra_LinearProblem> linear_problem;

    /**
     * 一个包含Trilinos求解器和预处理器对象的结构。
     *
     */
    std::unique_ptr<Amesos_BaseSolver> solver;

    /**
     * 存储这个特定求解器的标志的副本。
     *
     */
    const AdditionalData additional_data;
  };

} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_TRILINOS

 /*----------------------------   trilinos_solver.h ---------------------------*/ 

#endif
 /*----------------------------   trilinos_solver.h ---------------------------*/ 


