//include/deal.II-translator/lac/petsc_precondition_0.txt
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

#ifndef dealii_petsc_precondition_h
#  define dealii_petsc_precondition_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/subscriptor.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/lac/exceptions.h>

#    include <petscpc.h>

DEAL_II_NAMESPACE_OPEN



namespace PETScWrappers
{
  // forward declarations
#    ifndef DOXYGEN
  class MatrixBase;
  class VectorBase;
  class SolverBase;
#    endif

  /**
   * 使用PETSc功能的预处理器类的基类。这个层次中的类并不做很多事情，只是提供了一个函数，在求解器的预处理上下文上设置预处理和某些参数。这些类在这里基本上只是为了实现与已经用于deal.II求解器和预处理器类类似的接口。
   * 请注意，派生类只提供与PETSc相关功能的接口。PETSc并没有为所有的矩阵类型实现所有的预处理程序。特别是，有些预处理程序不能用于并行作业，例如ILU预处理程序。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionBase : public Subscriptor
  {
  public:
    /**
     * 构造器。
     *
     */
    PreconditionBase();

    /**
     * 解构器。
     *
     */
    virtual ~PreconditionBase();

    /**
     * 销毁预处理程序，留下一个像刚刚调用构造函数后的对象。
     *
     */
    void
    clear();

    /**
     * 对给定的src向量应用一次预处理程序。
     *
     */
    void
    vmult(VectorBase &dst, const VectorBase &src) const;

    /**
     * 对给定的src向量应用一次转置的预处理程序。
     *
     */
    void
    Tvmult(VectorBase &dst, const VectorBase &src) const;


    /**
     * 给予对底层PETSc对象的访问权。
     *
     */
    const PC &
    get_pc() const;

  protected:
    /**
     * PETSc预处理器对象
     *
     */
    PC pc;

    /**
     * 一个指向作为预处理器的矩阵的指针。
     *
     */
    Mat matrix;

    /**
     * 用于创建PETSc预处理对象的内部函数。如果调用两次就会失败。
     *
     */
    void
    create_pc();

    /**
     * 转换运算符，以获得代表该预处理程序的矩阵的表示。我们在实际求解器中使用它，我们需要将这个矩阵传递给PETSc求解器。
     *
     */
    operator Mat() const;

    // Make the solver class a friend, since it needs to call the conversion
    // operator.
    friend class SolverBase;
  };



  /**
   * 一个实现了使用PETSc雅可比预处理程序接口的类。
   * 参见基类 @ref PreconditionBase
   * 中的注释，了解该预处理器在什么情况下可以或不可以工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionJacobi : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {};

    /**
     * 空的构造器。在使用这个对象之前，你需要调用初始化（）。
     *
     */
    PreconditionJacobi() = default;


    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionJacobi(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 与上述相同，但不设置矩阵以形成预处理程序。
     * 旨在与SLEPc对象一起使用。
     *
     */
    PreconditionJacobi(
      const MPI_Comm &      communicator,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在创建没有参数的预处理程序时才会使用。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;

    /**
     * 在不知道特定矩阵的情况下初始化预处理器对象。这个函数在底层的PETSc对象被创建后为其设置适当的参数。
     *
     */
    void
    initialize();
  };



  /**
   * 一个实现接口的类，用于使用PETSc块状雅可比预处理程序。PETSc将 "块状雅可比 "
   * 定义为一个预处理程序，它查看矩阵的一些对角线块，然后定义一个预处理程序，其中预处理程序矩阵的块状结构只与这些对角线块相同，预处理程序的每个对角线块都是原矩阵相应块的逆的近似值。
   * 矩阵的块状结构是由MPI并行作业中各个处理器的自由度关联决定的。如果你在一个顺序作业(或只有一个进程的MPI作业)中使用这个预处理程序，那么整个矩阵就是唯一的块。
   * 默认情况下，PETSc使用矩阵的每个对角线块的ILU(0)分解来进行预处理。这一点可以改变，在PETSc手册的相关章节中有解释，但在这里没有实现。
   * 关于这个预处理程序在什么时候可以工作，什么时候不可以工作，请看基类
   * @ref PreconditionBase 中的注释。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionBlockJacobi : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {};

    /**
     * 空的构造器。在使用这个对象之前，你需要调用初始化（）。
     *
     */
    PreconditionBlockJacobi() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionBlockJacobi(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 与上述相同，但不设置矩阵以形成预处理程序。
     * 旨在与SLEPc对象一起使用。
     *
     */
    PreconditionBlockJacobi(
      const MPI_Comm &      communicator,
      const AdditionalData &additional_data = AdditionalData());


    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在创建没有参数的预处理程序时才会使用。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;

    /**
     * 在不知道特定矩阵的情况下初始化预处理程序对象。这个函数在底层的PETSc对象被创建后为其设置了适当的参数。
     *
     */
    void
    initialize();
  };



  /**
   * 一个实现接口的类，以使用PETSc SOR预处理器。
   * @note 只在与 PETScWrappers::SparseMatrix. 串联时工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将阻尼参数设置为1。
       *
       */
      AdditionalData(const double omega = 1);

      /**
       * 松弛参数。
       *
       */
      double omega;
    };

    /**
     * 空的构造函数。在使用这个对象之前，你需要调用初始化()。
     *
     */
    PreconditionSOR() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，则取额外的标志。
     *
     */
    PreconditionSOR(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现了使用PETSc SSOR预处理程序的接口的类。
   * @note 只在与 PETScWrappers::SparseMatrix. 串行时工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionSSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将阻尼参数设置为1。
       *
       */
      AdditionalData(const double omega = 1);

      /**
       * 松弛参数。
       *
       */
      double omega;
    };

    /**
     * 空的构造函数。在使用这个对象之前，你需要调用初始化（）。
     *
     */
    PreconditionSSOR() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，则取额外的标志。
     *
     */
    PreconditionSSOR(const MatrixBase &    matrix,
                     const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现使用PETSc Eisenstat预处理的接口的类，它在每个处理器所拥有的对角线块上实现SSOR。
   * 参见基类 @ref PreconditionBase
   * 中的注释，了解该预处理器在什么情况下可能工作或不工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionEisenstat : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将阻尼参数设置为1。
       *
       */
      AdditionalData(const double omega = 1);

      /**
       * 松弛参数。
       *
       */
      double omega;
    };

    /**
     * 空的构造函数。在使用这个对象之前，你需要调用初始化()。
     *
     */
    PreconditionEisenstat() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionEisenstat(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现接口的类，用于使用PETSc不完全Cholesky预处理器。
   * @note  只与 PETScWrappers::SparseMatrix. 串行工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionICC : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将填充参数设置为零。
       *
       */
      AdditionalData(const unsigned int levels = 0);

      /**
       * 填充参数。
       *
       */
      unsigned int levels;
    };

    /**
     * 空的构造函数。在使用这个对象之前，你需要调用initialize()。
     *
     */
    PreconditionICC() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，则取额外的标志。
     *
     */
    PreconditionICC(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现了使用PETSc ILU预处理程序的接口的类。
   * @note  仅在与 PETScWrappers::SparseMatrix. 的串行中工作。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionILU : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将填充参数设置为零。
       *
       */
      AdditionalData(const unsigned int levels = 0);

      /**
       * 填充参数。
       *
       */
      unsigned int levels;
    };

    /**
     * 空的构造函数。在使用这个对象之前，你需要调用initialize()。
     *
     */
    PreconditionILU() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionILU(const MatrixBase &    matrix,
                    const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现使用PETSc LU预处理的接口的类 (  @p PCLU).
   * 与PreconditionILU等类不同，这个类通常（取决于设置）执行矩阵的精确因式分解，所以没有必要用迭代求解器来包装它。该类通常与SolverPreOnly一起使用，以获得一个直接求解器。或者，你也可以直接使用
   * PreconditionBase::vmult() 。
   * @note
   * 这不是一个并行的预处理程序，所以它只适用于单处理器的串行计算。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionLU : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。(默认值取自PETSc库的函数PCCreate_LU)。
       *
       */
      AdditionalData(const double pivoting   = 1.e-6,
                     const double zero_pivot = 1.e-12,
                     const double damping    = 0.0);

      /**
       * 决定在LU分解过程中何时进行枢轴化。0.0表示没有透视，1.0表示完全透视。更多细节请参考PETSc手册。
       *
       */
      double pivoting;

      /**
       * 小枢轴被宣布为零的大小。更多细节请参考PETSc手册。
       *
       */
      double zero_pivot;

      /**
       * 在因式分解过程中，这个量被添加到矩阵的对角线上。
       *
       */
      double damping;
    };

    /**
     * 空构造函数。在使用这个对象之前，你需要调用初始化()。
     *
     */
    PreconditionLU() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，则取额外的标志。
     *
     */
    PreconditionLU(const MatrixBase &    matrix,
                   const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现接口的类，用于使用HYPRE套件中的BoomerAMG代数多网格预处理器。请注意，PETSc必须与HYPRE一起配置（例如，使用\--download-hypre=1）。
   * 该预处理程序确实支持并行分布式计算。见 step-40
   * 中的一个例子。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionBoomerAMG : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 定义了BoomerAMG的可用松弛类型。
       *
       */
      enum class RelaxationType
      {
        Jacobi,
        sequentialGaussSeidel,
        seqboundaryGaussSeidel,
        SORJacobi,
        backwardSORJacobi,
        symmetricSORJacobi,
        l1scaledSORJacobi,
        GaussianElimination,
        l1GaussSeidel,
        backwardl1GaussSeidel,
        CG,
        Chebyshev,
        FCFJacobi,
        l1scaledJacobi,
        None
      };

      /**
       * 构造函数。请注意，BoomerAMG提供了比这里所揭示的更多的选项来设置。
       *
       */
      AdditionalData(
        const bool           symmetric_operator               = false,
        const double         strong_threshold                 = 0.25,
        const double         max_row_sum                      = 0.9,
        const unsigned int   aggressive_coarsening_num_levels = 0,
        const bool           output_details                   = false,
        const RelaxationType relaxation_type_up   = RelaxationType::SORJacobi,
        const RelaxationType relaxation_type_down = RelaxationType::SORJacobi,
        const RelaxationType relaxation_type_coarse =
          RelaxationType::GaussianElimination,
        const unsigned int n_sweeps_coarse = 1,
        const double       tol             = 0.0,
        const unsigned int max_iter        = 1,
        const bool         w_cycle         = false);

      /**
       * 如果你有一个对称的系统矩阵，并且你想使用像CG这样假设有对称预处理的求解器，请将这个标志设置为
       * "true"。
       *
       */
      bool symmetric_operator;

      /**
       * 节点被认为是强连接的阈值。参见HYPRE_BoomerAMGSetStrongThreshold（）。推荐值是2D问题的0.25和3D问题的0.5，但这取决于问题。
       *
       */
      double strong_threshold;

      /**
       * 如果设置为小于1.0的值，那么矩阵的对角线主导部分将被视为没有强连接节点。如果由对角线条目加权的行和大于给定值，则被认为是对角线主导的。这个功能可以通过设置数值为1.0来关闭。这是默认的，因为有些矩阵可能导致只有对角线主导的条目，因此没有构建多网格层次。在BoomerAMG中，这个默认值是0.9。当你尝试这样做的时候，请检查所创建的层次的合理数量。
       *
       */
      double max_row_sum;

      /**
       * 积极粗化的层数。增加这个值可以减少构建时间和内存需求，但是可能会降低效果。
       *
       */
      unsigned int aggressive_coarsening_num_levels;

      /**
       * 当构建预处理程序时，将此标志设置为 "true
       * "会产生来自HYPRE的调试输出。
       *
       */
      bool output_details;

      /**
       * 选择松弛类型了。
       *
       */
      RelaxationType relaxation_type_up;

      /**
       * 选择松弛类型向下。
       *
       */
      RelaxationType relaxation_type_down;

      /**
       * 选择粗放型的放松。
       *
       */
      RelaxationType relaxation_type_coarse;

      /**
       * 选择粗放网格上的扫描次数。
       *
       */
      unsigned int n_sweeps_coarse;

      /**
       * 选择BommerAMG公差。
       *
       */
      double tol;

      /**
       * 选择BommerAMG的最大周期数。
       *
       */
      unsigned int max_iter;

      /**
       * 定义是否应该使用w型循环而不是标准设置的v型循环。
       *
       */
      bool w_cycle;
    };

    /**
     * 空构造函数。在使用这个对象之前，你需要调用initialize()。
     *
     */
    PreconditionBoomerAMG() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionBoomerAMG(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 与上述相同，但不设置矩阵以形成预处理程序。
     * 目的是与SLEPc对象一起使用。
     *
     */
    PreconditionBoomerAMG(
      const MPI_Comm &      communicator,
      const AdditionalData &additional_data = AdditionalData());


    /**
     * 初始化预处理对象并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时自动调用，只有在创建没有参数的预处理程序时才会使用。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  protected:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;

    /**
     * 在不知道特定矩阵的情况下初始化预处理程序对象。这个函数在底层PETSc对象被创建后为其设置了适当的参数。
     *
     */
    void
    initialize();
  };



  /**
   * 一个实现接口的类，用于使用HYPRE套件中的ParaSails稀疏近似反预处理器。请注意，PETSc必须与HYPRE一起配置（例如，使用\--download-hypre=1）。
   * ParaSails使用最小二乘法来计算稀疏的近似逆。所使用的稀疏模式是稀疏矩阵的幂的模式。ParaSails还使用了一种后过滤技术，以减少应用预处理程序的成本。
   * ParaSails使用因子化SPD预处理器解决对称正定（SPD）问题，也可以使用非因子化预处理器解决一般（非对称和/或不确定）问题。问题类型必须在
   * @p AdditionalData.
   * 中设置，预处理程序确实支持并行分布式计算。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionParaSails : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。
       *
       */
      AdditionalData(const unsigned int symmetric      = 1,
                     const unsigned int n_levels       = 1,
                     const double       threshold      = 0.1,
                     const double       filter         = 0.05,
                     const bool         output_details = false);

      /**
       * 这个参数指定了要解决的问题的类型。        <ul>   <li>   @p 0:  非对称和/或不确定问题，以及非对称预处理  <li>   @p 1:  SPD问题，以及SPD（因子）预处理  <li>   @p 2:  非对称、确定问题，以及SPD（因子）预处理  </ul>  默认为<tt>symmetric = 1</tt>。
       *
       */
      unsigned int symmetric;

      /**
       * 用于近似求逆的稀疏模式是幂级<tt>B^m</tt>的模式，其中<tt>B</tt>已经从给定的矩阵<tt>A</tt>中稀疏化，<tt>n_level</tt>等于<tt>m+1</tt>。
       * 默认值是<tt>n_levels = 1</tt>。
       *
       */
      unsigned int n_levels;

      /**
       * 稀疏化是通过丢弃幅度小于<tt>thresh</tt>的非零点来实现的。较低的<tt>thresh</tt>值会导致更精确，但也更昂贵的预处理程序。
       * 默认值是<tt>thresh = 0.1</tt>。设置<tt>thresh <
       * 0</tt>会自动选择一个阈值，这样<tt>thresh</tt>就代表了被放弃的非零元素的比例。例如，如果<tt>thresh
       * =
       *
       * - .9</tt>，那么<tt>B</tt>将包含给定矩阵<tt>A</tt>中约百分之十的非零元素。
       *
       */
      double threshold;

      /**
       * 过滤是一个后处理过程，<tt>filter</tt>代表在创建近似的反稀疏模式后被丢弃的非零元素的一部分。默认值是<tt>filter
       * = 0.05</tt>。设置<tt>filter <
       * 0</tt>会自动选择一个值，这样<tt>-filter</tt>代表被丢弃的非零元素的比例。例如，如果<tt>thresh
       * =
       *
       *
       *
       *
       *
       *
       * - .9</tt>，那么在计算的近似逆中约有90%的条目被放弃。
       *
       */
      double filter;

      /**
       * 当构建预处理程序时，将此标志设置为 "true
       * "会产生来自HYPRE的输出。
       *
       */
      bool output_details;
    };



    /**
     * 空的构造器。在使用这个对象之前，你需要调用初始化()。
     *
     */
    PreconditionParaSails() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志，也要取。
     *
     */
    PreconditionParaSails(
      const MatrixBase &    matrix,
      const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时被自动调用，只有在你创建没有参数的预处理程序时才会用到。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  private:
    /**
     * 为这个特定的先决条件器存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };



  /**
   * 一个实现非预处理方法的类。
   * @ingroup PETScWrappers
   *
   */
  class PreconditionNone : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。
     *
     */
    struct AdditionalData
    {};

    /**
     * 空的构造器。在使用这个对象之前，你需要调用初始化（）。
     *
     */
    PreconditionNone() = default;

    /**
     * 构造函数。取用于形成预处理程序的矩阵，如果有额外的标志的话。该矩阵在计算中被完全忽略。
     *
     */
    PreconditionNone(const MatrixBase &    matrix,
                     const AdditionalData &additional_data = AdditionalData());

    /**
     * 初始化预处理对象，并计算在求解器中应用它所需的所有数据。这个函数在调用具有相同参数的构造函数时自动调用，只有在创建没有参数的预处理程序时才会用到。矩阵在计算中被完全忽略。
     *
     */
    void
    initialize(const MatrixBase &    matrix,
               const AdditionalData &additional_data = AdditionalData());

  private:
    /**
     * 为这个特定的预处理程序存储一份标志的副本。
     *
     */
    AdditionalData additional_data;
  };

  /**
   * 向后兼容的别名。    @deprecated  使用
   * PETScWrappers::PreconditionBase 代替。
   *
   */
  using PreconditionerBase DEAL_II_DEPRECATED = PreconditionBase;
} // namespace PETScWrappers



DEAL_II_NAMESPACE_CLOSE


#  endif // DEAL_II_WITH_PETSC

#endif
 /*--------------------------- petsc_precondition.h --------------------------*/ 


