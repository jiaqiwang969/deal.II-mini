//include/deal.II-translator/lac/trilinos_precondition_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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

#ifndef dealii_trilinos_precondition_h
#  define dealii_trilinos_precondition_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_TRILINOS

#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/la_parallel_vector.h>
#    include <deal.II/lac/trilinos_vector.h>

#    include <memory>

#    ifdef DEAL_II_WITH_MPI
#      include <Epetra_MpiComm.h>
#    else
#      include <Epetra_SerialComm.h>
#    endif
#    include <Epetra_Map.h>
#    include <Epetra_MultiVector.h>
#    include <Epetra_RowMatrix.h>
#    include <Epetra_Vector.h>
#    include <Teuchos_ParameterList.hpp>

// forward declarations
#    ifndef DOXYGEN
class Ifpack_Preconditioner;
class Ifpack_Chebyshev;
namespace ML_Epetra
{
  class MultiLevelPreconditioner;
}
#    endif

DEAL_II_NAMESPACE_OPEN

// forward declarations
#    ifndef DOXYGEN
template <typename number>
class SparseMatrix;
template <typename number>
class Vector;
class SparsityPattern;
#    endif

/*!   @addtogroup  TrilinosWrappers  
     * @{ 

 
*
*/

namespace TrilinosWrappers
{
  // forward declarations
  class SparseMatrix;
  class BlockSparseMatrix;
  class SolverBase;

  /**
   * 所有基于Trilinos稀疏矩阵的预处理器的基类。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionBase : public Subscriptor
  {
  public:
    /**
     * 声明容器大小的类型。
     *
     */
    using size_type = dealii::types::global_dof_index;

    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。
     *
     */
    struct AdditionalData
    {};

    /**
     * 构造函数。不做任何事情。派生类的<tt>initialize</tt>函数将不得不从给定的稀疏矩阵中创建预处理程序。
     *
     */
    PreconditionBase();

    /**
     * 复制构造函数。
     *
     */
    PreconditionBase(const PreconditionBase &);

    /**
     * 解构器。
     *
     */
    ~PreconditionBase() override = default;

    /**
     * 销毁预处理程序，留下一个像刚刚调用构造函数后的对象。
     *
     */
    void
    clear();

    /**
     * 返回与该矩阵一起使用的MPI通信器对象。
     *
     */
    MPI_Comm
    get_mpi_communicator() const;

    /**
     * 设置一个内部标志，使矩阵进行的所有操作，即乘法，都以转置的顺序进行。然而，这并不能直接将矩阵重塑为转置的形式，所以在使用这个标志时要注意。
     * @note
     * 连续调用此函数的任何偶数次，都将使对象返回到其原始状态。
     *
     */
    void
    transpose();

    /**
     * 应用预处理程序。
     *
     */
    virtual void
    vmult(MPI::Vector &dst, const MPI::Vector &src) const;

    /**
     * 应用转置预处理程序。
     *
     */
    virtual void
    Tvmult(MPI::Vector &dst, const MPI::Vector &src) const;

    /**
     * 在deal.II数据结构上应用预处理程序，而不是在Trilinos包装类中提供的数据结构。
     *
     */
    virtual void
    vmult(dealii::Vector<double> &dst, const dealii::Vector<double> &src) const;

    /**
     * 在deal.II数据结构上应用转置预处理程序，而不是在Trilinos包装类中提供的那些。
     *
     */
    virtual void
    Tvmult(dealii::Vector<double> &      dst,
           const dealii::Vector<double> &src) const;

    /**
     * 在deal.II并行数据结构上应用预处理程序，而不是在Trilinos包装类中提供的结构。
     *
     */
    virtual void
    vmult(dealii::LinearAlgebra::distributed::Vector<double> &      dst,
          const dealii::LinearAlgebra::distributed::Vector<double> &src) const;

    /**
     * 在deal.II并行数据结构上应用转置预处理程序，而不是在Trilinos包装类中提供的数据结构。
     *
     */
    virtual void
    Tvmult(dealii::LinearAlgebra::distributed::Vector<double> &      dst,
           const dealii::LinearAlgebra::distributed::Vector<double> &src) const;

    /**
     * @name  访问底层Trilinos数据
     *
     */
    //@{
    /**
     * 从一个未初始化的对象中调用这个函数将导致一个异常。
     *
     */
    Epetra_Operator &
    trilinos_operator() const;
    //@}

    /**
     * @name  分区器
     *
     */
    //@{

    /**
     * 返回这个矩阵的域空间的分区，即这个矩阵要与之相乘的向量的分区。
     *
     */
    IndexSet
    locally_owned_domain_indices() const;

    /**
     * 返回该矩阵的范围空间的划分，即由矩阵-向量乘积产生的向量的划分。
     *
     */
    IndexSet
    locally_owned_range_indices() const;

    //@}

    /**
     * @addtogroup  Exceptions
     *
     */
    //@{
    /**
     * 异常情况。
     *
     */
    DeclException1(ExcNonMatchingMaps,
                   std::string,
                   << "The sparse matrix the preconditioner is based on "
                   << "uses a map that is not compatible to the one in vector "
                   << arg1 << ". Check preconditioner and matrix setup.");
    //@}

    friend class SolverBase;

  protected:
    /**
     * 这是一个指向预处理器对象的指针，在应用预处理器时使用。
     *
     */
    Teuchos::RCP<Epetra_Operator> preconditioner;

    /**
     * 内部通信模式，以防矩阵需要从deal.II格式中复制。
     *
     */
#    ifdef DEAL_II_WITH_MPI
    Epetra_MpiComm communicator;
#    else
    Epetra_SerialComm communicator;
#    endif

    /**
     * 内部Trilinos图，以防矩阵需要从deal.II格式中复制。
     *
     */
    std::shared_ptr<Epetra_Map> vector_distributor;
  };


  /**
   * 一个用于Trilinos矩阵的（点）雅可比预处理的封装类。这个预处理器既可以在串行中工作，也可以在并行中工作，这取决于它所基于的矩阵。
   * AdditionalData数据结构允许设置预处理器选项。
   * 对于雅可比预处理程序，这些选项是阻尼参数<tt>omega</tt>和<tt>min_diagonal</tt>参数，可以用来使预处理程序工作，即使矩阵的对角线上包含一些零元素。默认设置是阻尼参数为1，对角线增量为0。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionJacobi : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。参数<tt>omega</tt>指定雅可比预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在一些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵不太远的地方工作。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将阻尼参数设置为1，并且不修改对角线。
       *
       */
      AdditionalData(const double       omega        = 1,
                     const double       min_diagonal = 0,
                     const unsigned int n_sweeps     = 1);

      /**
       * 这指定了雅可比预处理程序中的松弛参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当雅可比预处理器用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以零。
       *
       */
      double min_diagonal;

      /**
       * 设置在vmult()操作中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取预处理对象应建立的稀疏矩阵，如果有额外的标志（阻尼参数等）。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * Trilinos矩阵的（点）SSOR预处理器的封装类。这个预处理程序既可以在串行中工作，也可以在并行中工作，这取决于它所基于的矩阵。
   * AdditionalData数据结构允许设置预处理器选项。
   * 对于SSOR预处理程序，这些选项是阻尼/松弛参数<tt>omega</tt>，<tt>min_diagonal</tt>参数，可以用来使预处理程序工作，即使矩阵在对角线上包含一些零元素，以及<tt>overlap</tt>参数，决定在不同MPI进程的矩阵分区之间是否应该有重叠以及有多少。默认设置是松弛参数为1，对角线增强为0，重叠为0。
   * 请注意，SSOR预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>与SSOR <em> 近似求解 </em>
   * 作为内部求解器，基于外部并行分区。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionSSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。参数<tt>omega</tt>指定了SSOR预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在一些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵相差不大的矩阵上工作。最后，<tt>overlap</tt>管理预处理程序并行运行时分区的重叠，形成所谓的加法施瓦茨预处理程序。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将阻尼参数设置为1，我们不修改对角线，并且没有重叠（即在并行情况下，我们运行BlockJacobi预处理程序，其中每个块被SSOR近似倒置）。
       *
       */
      AdditionalData(const double       omega        = 1,
                     const double       min_diagonal = 0,
                     const unsigned int overlap      = 0,
                     const unsigned int n_sweeps     = 1);

      /**
       * 这指定了SSOR预处理程序中的（过度）放松参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当SSOR预处理程序用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以对角线元素。
       *
       */
      double min_diagonal;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;

      /**
       * 设置在vmult()操作过程中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取出预处理对象应建立的稀疏矩阵，以及额外的标志（阻尼参数、并行计算中的重叠等），如果有的话。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * Trilinos矩阵的（点）SOR预处理的封装类。这个预处理程序既可以在串行中工作，也可以在并行中工作，这取决于它所基于的矩阵。
   * AdditionalData数据结构允许设置预处理器选项。
   * 对于SOR预处理程序，这些选项是阻尼/松弛参数<tt>omega</tt>，<tt>min_diagonal</tt>参数，可以用来使预处理程序工作，即使矩阵在对角线上包含一些零元素，以及<tt>overlap</tt>参数，决定在不同MPI进程的矩阵分区之间是否应该有重叠以及有多少。默认设置是松弛参数为1，对角线增强为0，重叠为0。
   * 请注意，SOR预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，有一个SOR <em> 近似求解 </em>
   * 作为内部求解器，基于外部并行划分。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。参数<tt>omega</tt>指定了SOR预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在一些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵相差不大的矩阵上工作。最后，<tt>overlap</tt>管理预处理程序并行运行时分区的重叠，形成所谓的加法施瓦茨预处理程序。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，将阻尼参数设置为1，我们不修改对角线，并且没有重叠（即在并行情况下，我们运行BlockJacobi预处理程序，其中每个块都被一个SOR近似倒置。
       *
       */
      AdditionalData(const double       omega        = 1,
                     const double       min_diagonal = 0,
                     const unsigned int overlap      = 0,
                     const unsigned int n_sweeps     = 1);

      /**
       * 这指定了SOR预处理程序中的（过度）放松参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当SOR预处理程序用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以对角线元素。
       *
       */
      double min_diagonal;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;

      /**
       * 设置在vmult()操作过程中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取出预处理对象应建立的稀疏矩阵，以及额外的标志（阻尼参数、并行计算中的重叠等），如果有的话。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * Trilinos矩阵的块状雅可比预处理器的封装类。
   * 与PreconditionSOR不同的是，在PreconditionSOR中，每一行都是单独处理的，这个方案收集特定大小的块，并同时对所有这些行进行全矩阵反演。Trilinos允许选择几种策略来选择哪些行构成一个块，包括
   * "线性"（即把矩阵的局部范围划分为块大小的片断）、"贪婪
   * "或 "metis"。请注意，术语 <em> 块Jacobi </em>
   * 并不涉及MPI设置中的可能块，而是从每个处理器本地的稀疏矩阵中提取的密集矩阵的小块。
   * AdditionalData数据结构允许设置预处理程序选项。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionBlockJacobi : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。参数<tt>block_size</tt>设置小块的大小。建议选择这个参数不要太大（最多几百个），因为这个实现对块使用密集的矩阵。参数<tt>block_creation_type</tt>允许传递给Ifpack寻找块的策略。参数<tt>omega</tt>指定了SOR预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在某些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵不太远的地方工作。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，使用块大小为1，使用行的线性细分，将阻尼参数设置为1，并且不修改对角线。
       *
       */
      AdditionalData(const unsigned int block_size          = 1,
                     const std::string &block_creation_type = "linear",
                     const double       omega               = 1,
                     const double       min_diagonal        = 0,
                     const unsigned int n_sweeps            = 1);

      /**
       * 这指定了块的大小。
       *
       */
      unsigned int block_size;

      /**
       * 传递给Ifpack块放松（变量'partitioner:
       * type'）的块的创建策略，该字符串为给定值。
       * 在Ifpack中可用的类型包括
       * "线性"（即把矩阵的局部范围划分为块大小的片状），"贪婪""metis"。
       * 完整的列表请参见Ifpack的文档。
       *
       */
      std::string block_creation_type;

      /**
       * 这指定了雅可比预处理程序中的（过度）放松参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当块状雅可比预处理器用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以对角线元素。
       *
       */
      double min_diagonal;

      /**
       * 设置在vmult()操作中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取预处理对象应建立的稀疏矩阵，如果有额外的标志（阻尼参数等）。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 一个用于特里诺斯矩阵的块状SSOR预处理的封装类。与PreconditionSSOR不同的是，在PreconditionSSOR中，每一行都是单独处理的（从点到面），这个方案收集一定大小的块，并同时对所有这些行进行全矩阵逆变。Trilinos允许选择几种策略来选择哪些行构成一个块，包括
   * "线性"（即把矩阵的局部范围划分为块大小的片断）、"贪婪
   * "或 "metis"。
   * AdditionalData数据结构允许设置预处理选项。
   * 请注意，这个预处理程序的并行应用实际上是一个块-雅各比预处理程序，其（外部）块大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，以块SSOR <em> 近似求解 </em>
   * 为内部求解器，基于外部并行分区。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionBlockSSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加标志。参数<tt>block_size</tt>设置小块的大小。建议选择这个参数不要太大（最多几百个），因为这个实现对块使用密集的矩阵。参数<tt>block_creation_type</tt>允许传递给Ifpack寻找块的策略。参数<tt>omega</tt>指定了SSOR预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在某些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵相差不大的矩阵上工作。最后，<tt>overlap</tt>管理预处理程序并行运行时分区的重叠，形成所谓的加法施瓦茨预处理程序。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，使用块大小为1，使用行的线性细分，将阻尼参数设置为1，我们不修改对角线，并且没有重叠（即在并行情况下，我们运行BlockJacobi预处理程序，每个块大约由一个块SOR来反转）。
       *
       */
      AdditionalData(const unsigned int block_size          = 1,
                     const std::string &block_creation_type = "linear",
                     const double       omega               = 1,
                     const double       min_diagonal        = 0,
                     const unsigned int overlap             = 0,
                     const unsigned int n_sweeps            = 1);

      /**
       * 这指定了块的大小。
       *
       */
      unsigned int block_size;

      /**
       * 创建块的策略传递给Ifpack块放松（变量'partitioner:
       * type'），以这个字符串为给定值。
       * 在Ifpack中可用的类型包括
       * "线性"（即把矩阵的局部范围划分为块大小的片状），"贪婪""metis"。
       * 完整的列表请参见Ifpack的文档。
       *
       */
      std::string block_creation_type;

      /**
       * 这指定了SOR预处理程序中的（过度）放松参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当SSOR预处理器用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以对角线元素。
       *
       */
      double min_diagonal;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;

      /**
       * 设置在vmult()操作过程中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取出预处理对象应建立的稀疏矩阵，以及额外的标志（阻尼参数、并行计算中的重叠等），如果有的话。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * Trilinos矩阵的块状SOR预处理器的封装类。与PreconditionSOR不同的是，在PreconditionSOR中，每一行都是单独处理的，这个方案收集一定大小的块，并同时对所有这些行进行全矩阵逆变。Trilinos允许选择几种策略来选择哪些行构成一个块，包括
   * "线性"（即把矩阵的局部范围划分为块大小的片断）、"贪婪
   * "或 "metis"。
   * AdditionalData数据结构允许设置预处理选项。
   * 请注意，这个预处理程序的并行应用实际上是一个块-雅各比预处理程序，其（外部）块大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，有一个块SOR <em> 近似求解 </em>
   * 作为内部求解器，基于外部并行分区。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionBlockSOR : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的标志。参数<tt>block_size</tt>设置小块的大小。建议选择这个参数不要太大（最多几百个），因为这个实现对块使用密集的矩阵。参数<tt>block_creation_type</tt>允许传递给Ifpack寻找块的策略。参数<tt>omega</tt>指定了SOR预处理程序中的松弛参数。参数<tt>min_diagonal</tt>可用于在某些对角线元素为零时也可以应用预处理程序。在默认的应用中，这意味着我们要除以0，所以通过设置参数<tt>min_diagonal</tt>到一个小的非零值，SOR将在一个与我们要处理的矩阵相差不大的矩阵上工作。最后，<tt>overlap</tt>管理预处理程序并行运行时分区的重叠，形成所谓的加法施瓦茨预处理程序。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。默认情况下，使用块大小为1，使用行的线性细分，将阻尼参数设置为1，我们不修改对角线，并且没有重叠（即在并行情况下，我们运行BlockJacobi预处理程序，每个块大约由一个块SOR来反转）。
       *
       */
      AdditionalData(const unsigned int block_size          = 1,
                     const std::string &block_creation_type = "linear",
                     const double       omega               = 1,
                     const double       min_diagonal        = 0,
                     const unsigned int overlap             = 0,
                     const unsigned int n_sweeps            = 1);

      /**
       * 这指定了块的大小。
       *
       */
      unsigned int block_size;

      /**
       * 创建块的策略传递给Ifpack块放松（变量'partitioner:
       * type'），以这个字符串为给定值。
       * 在Ifpack中可用的类型包括
       * "线性"（即把矩阵的局部范围划分为块大小的片状），"贪婪""metis"。
       * 完整的列表请参见Ifpack的文档。
       *
       */
      std::string block_creation_type;

      /**
       * 这指定了SOR预处理程序中的（过度）放松参数。
       *
       */
      double omega;

      /**
       * 这指定了对角线元素应该具有的最小值。
       * 当SOR预处理器用于对角线元素为零的矩阵时，这可能是必要的。在这种情况下，直接应用是不可能的，因为我们要除以对角线元素。
       *
       */
      double min_diagonal;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;

      /**
       * 设置在vmult()操作过程中应该应用多少次给定的操作。
       *
       */
      unsigned int n_sweeps;
    };

    /**
     * 取出预处理对象应建立的稀疏矩阵，以及额外的标志（阻尼参数、并行计算中的重叠等），如果有的话。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 用于 @em
   * 对称Trilinos矩阵的不完全Cholesky因子化（IC）预处理的封装类。这个预处理程序既可以串行工作也可以并行工作，这取决于它所基于的矩阵。一般来说，不完全因式分解并不包含所有会出现在完全因式分解中的填充元素（这是直接求解的基础）。Trilinos允许设置填充元素的数量，由附加数据参数<tt>ic_fill</tt>控制，因此可以逐步选择只对稀疏矩阵结构进行因式分解（<tt>ic_fill=0</tt>）到完全因式分解（<tt>ic_fill</tt>在10到50之间，取决于PDE问题的空间维度和有限元基函数的程度；通常，更多需要的填充元素需要将该参数设置为更高的整数值）。
   * AdditionalData数据结构允许设置预处理器选项。
   * 除了填充参数外，这些选项还有一些用于扰动的选项（详见AdditionalData结构的文档），以及一个参数<tt>overlap</tt>，该参数决定了不同MPI进程上的矩阵分区之间是否应该有重叠以及重叠的程度。
   * 默认设置为：额外填充为0，绝对增容为0，相对增容为1，重叠为0。
   * 请注意，IC预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，有一个IC <em> 近似求解 </em>
   * 作为内部求解器，基于（外部）并行分区的。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionIC : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加参数。Trilinos
     * IC分解允许一些填充，所以它实际上是一个阈值不完整的Cholesky分解。填充量，以及该预处理程序使用的内存量，由参数<tt>ic_fill</tt>控制，该参数将其指定为双数。在形成预处理程序时，对于某些问题，不好的条件（或者只是运气不好）会导致预处理程序的条件很差。因此，将对角线扰动添加到原始矩阵中，并为这个稍好的矩阵形成预处理程序，会有所帮助。<tt>ic_atol</tt>是一个绝对扰动，在形成prec之前加到对角线上，<tt>ic_rtol</tt>是一个缩放系数
     * $rtol \geq 1$
     * 。最后一个参数指定了预处理程序并行运行时分区的重叠情况。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，将下降容忍度设置为0，额外填充物的级别设置为0（只使用矩阵结构，不产生任何额外的填充物），容忍度级别分别为0和1，在并行执行的情况下，重叠度为0。在并行情况下，IC的块状应用中的这种重叠使预处理器成为所谓的加性施瓦兹预处理器。
       *
       */
      AdditionalData(const unsigned int ic_fill = 0,
                     const double       ic_atol = 0.,
                     const double       ic_rtol = 1.,
                     const unsigned int overlap = 0);

      /**
       * 这指定了除了稀疏矩阵结构之外的额外填充元素的数量。当<tt>ic_fill</tt>较大时，这意味着将添加许多填充元素，从而使IC预处理程序更接近于直接稀疏Cholesky分解。但是请注意，这将极大地增加内存需求，特别是当预处理程序在三维中使用时。
       *
       */
      unsigned int ic_fill;

      /**
       * 这指定了将被添加到矩阵对角线上的绝对扰动量，有时这有助于获得更好的预处理。
       *
       */
      double ic_atol;

      /**
       * 这指定了矩阵的对角线将被缩放的系数，这有时有助于获得更好的预处理程序。
       *
       */
      double ic_rtol;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;
    };

    /**
     * 初始化函数。接收预处理程序应该计算的矩阵，如果有额外的标志，则接收额外的标志。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 一个用于特里诺斯矩阵的不完全LU分解（ILU(k)）预处理的封装类。这个预处理程序既可以串行工作也可以并行工作，这取决于它所基于的矩阵。一般来说，不完全因式分解并不包含所有会出现在完全因式分解中的填充元素（这是直接求解的基础）。Trilinos允许设置填充元素的数量，由附加数据参数<tt>ilu_fill</tt>控制，因此可以逐步选择只对稀疏矩阵结构进行因式分解（<tt>ilu_fill=0</tt>）到完全因式分解（<tt>ilu_fill</tt>在10到50之间，取决于PDE问题的空间维度和有限元基函数的程度；通常，更多需要填充的元素需要将该参数设置为更高的整数值）。
   * AdditionalData数据结构允许设置预处理器选项。
   * 详情见AdditionalData结构的文档。
   * 请注意，ILU预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，以ILU <em> 近似求解 </em>
   * 为内部求解器，基于（外部）并行划分。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionILU : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加参数。      <ul>   <li>   @p ilu_fill:  这指定了除原始稀疏矩阵结构外的额外填充元素的数量。如果 $k$ 是 @p 填充， $A^{k+1}$ 的稀疏模式被用于存储高斯消除的结果。这在文献中被称为ILU(  $k$  ) 。 当 @p fill 较大时，预调节器更接近于（直接）稀疏LU分解。但是请注意，这将极大地增加内存需求，特别是当预处理程序在三维中使用时。          <li>   @p ilu_atol  和  @p ilu_rtol:  这两个参数允许对矩阵的对角线进行扰动，这有时可以帮助获得更好的预处理程序，特别是在条件不好的情况下。    在因式分解之前，对角线条目 $a_{ii}$ 被 $\alpha sign(a_{ii}) + \beta a_{ii}$ 取代，其中 $\alpha\geq 0$ 是绝对阈值 @p ilu_atol ， $\beta\geq 1$ 是相对阈值 @p ilu_rtol.  默认值（ $\alpha = 0$  ， $\beta = 1$  ）因此使用原始，未修改的对角线条目。建议值是 $10^{-5}$ 到 $10^{-2}$ 的顺序，对于 @p ilu_rtol.   <li>   @p overlap:  这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该是多大。    0的重叠对应于每个处理器上的块状对角线分解，1的重叠将额外包括一个行j，如果在自己的行中的第j列有一个非零条目。更高的重叠数以递归的方式相应地工作。增加 @p overlap 将增加通信和存储成本。根据IFPACK文档，1的重叠通常是有效的，超过3的值很少需要。          </ul>
     *
     */
    struct AdditionalData
    {
      /**
       * 具有所有参数默认值的构造函数。
       *
       */
      AdditionalData(const unsigned int ilu_fill = 0,
                     const double       ilu_atol = 0.,
                     const double       ilu_rtol = 1.,
                     const unsigned int overlap  = 0);

      /**
       * 额外的填充，见上面的类文件。
       *
       */
      unsigned int ilu_fill;

      /**
       * 添加到对角线条目的扰动量。详见上面的类文件。
       *
       */
      double ilu_atol;

      /**
       * 对角线条目的缩放演员。详见上面的类文件。
       *
       */
      double ilu_rtol;

      /**
       * 处理器之间的重叠。详见上面的类文件。
       *
       */
      unsigned int overlap;
    };

    /**
     * 初始化函数。接受用于形成预处理程序的矩阵，如果有额外的标志，则接受额外的标志。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 用于Trilinos矩阵的阈值不完全LU分解（ILU-T）预处理的封装类。这个预处理程序既可以串行工作也可以并行工作，这取决于它所基于的矩阵。一般来说，不完全因式分解并不包含所有会出现在完全因式分解中的填充元素（那是直接求解的基础）。对于ILU-T预处理程序，参数<tt>ilut_drop</tt>让用户指定哪些元素应该被放弃（即不应该成为不完全分解的一部分）。Trilinos首先计算一行的完整分解，然后跳过那些低于阈值的元素。这是与非阈值ILU预处理程序的主要区别，其中参数<tt>ilut_fill</tt>控制着不完全因式分解结构。这个参数在这里也是可用的，但在这里只提供一些额外的信息。
   * AdditionalData数据结构允许设置预处理程序选项。
   * 除了填充参数外，这些选项还有一些用于扰动的选项（详见AdditionalData结构的文档），以及一个参数<tt>overlap</tt>，该参数决定了不同MPI进程上的矩阵分区之间是否应该有重叠以及有多少重叠。默认设置是额外填充为0，绝对增容为0，相对增容为1，重叠为0。
   * 请注意，ILU-T预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，以ILU <em> 近似求解 </em>
   * 为内部求解器，基于（外部）并行分区。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionILUT : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加参数。Trilinos
     * ILU-T分解允许一些填充，所以它实际上是一个不完整的LU因子化门槛。填充量以及该预处理程序使用的内存量由参数<tt>ilut_drop</tt>和<tt>ilut_fill</tt>控制，这些参数指定了关于哪些值应该形成不完全因式分解的阈值以及额外填充的程度。在形成预处理程序时，对于某些问题，不好的条件（或者只是运气不好）会导致预处理程序的条件很差。因此，将对角线扰动添加到原始矩阵中，并为这个稍好的矩阵形成预处理器，会有帮助。<tt>ilut_atol</tt>是一个绝对扰动，在形成prec之前加到对角线上，<tt>ilu_rtol</tt>是一个缩放系数
     * $rtol \geq 1$
     * 。最后一个参数指定了预处理程序并行运行时分区的重叠情况。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，没有元素会被丢弃，额外的填充水平被设置为零（只使用矩阵结构，不产生任何额外的填充，除了不丢弃大元素而产生的填充），容忍度分别为0和1，在并行执行的情况下，重叠度为零。在并行情况下，ILU的块应用中的这种重叠使预处理程序成为所谓的加性施瓦茨预处理程序。
       *
       */
      AdditionalData(const double       ilut_drop = 0.,
                     const unsigned int ilut_fill = 0,
                     const double       ilut_atol = 0.,
                     const double       ilut_rtol = 1.,
                     const unsigned int overlap   = 0);

      /**
       * 这指定了在形成不完整的LU分解时应该放弃的元素的相对大小，并带有阈值。
       *
       */
      double ilut_drop;

      /**
       * 这指定了除了稀疏矩阵结构之外的额外填充元素的数量。当<tt>ilu_fill</tt>较大时，这意味着许多填充元素将被添加，从而使ILU预处理程序更接近于（直接）稀疏LU分解。但是请注意，这将极大地增加内存需求，特别是当预处理程序用于三维时。
       *
       */
      unsigned int ilut_fill;

      /**
       * 这指定了将被添加到矩阵对角线上的绝对扰动量，这有时可以帮助获得更好的预处理程序。
       *
       */
      double ilut_atol;

      /**
       * 这指定了矩阵的对角线将被缩放的系数，这有时有助于获得更好的预处理程序。
       *
       */
      double ilut_rtol;

      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;
    };

    /**
     * 初始化函数。接受用于形成预处理程序的矩阵，如果有额外的标志，则接受额外的标志。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 稀疏直接LU分解的包装类，适用于Trilinos矩阵的并行块。当以串行方式运行时，这相当于对矩阵进行直接求解。
   * AdditionalData数据结构允许设置预处理选项。
   * 请注意，块状直接求解预处理的并行应用实际上是一个块状Jacobi预处理，块状大小等于本地矩阵大小。说得更专业一点，这个并行操作是一个<a
   * href="http://en.wikipedia.org/wiki/Additive_Schwarz_method">additive
   * Schwarz method</a>，有一个 <em> 精确求解 </em>
   * 作为内部求解器，基于（外部）并行分区的。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionBlockwiseDirect : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送附加参数。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。
       *
       */
      AdditionalData(const unsigned int overlap = 0);


      /**
       * 这决定了并行应用中每个处理器上的局部矩阵部分的重叠应该有多大。
       *
       */
      unsigned int overlap;
    };

    /**
     * 初始化函数。接受用于形成预处理程序的矩阵，如果有额外的标志，则接受额外的标志。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 一个用于特里诺斯矩阵的切比雪夫预处理的封装类。
   * AdditionalData数据结构允许设置预处理器选项。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionChebyshev : public PreconditionBase
  {
  public:
    /**
     * 标准化的数据结构，用于向预处理程序输送额外的参数。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造器。
       *
       */
      AdditionalData(const unsigned int degree           = 1,
                     const double       max_eigenvalue   = 10.,
                     const double       eigenvalue_ratio = 30.,
                     const double       min_eigenvalue   = 1.,
                     const double       min_diagonal     = 1e-12,
                     const bool         nonzero_starting = false);

      /**
       * 这决定了切比雪夫多项式的程度。多项式的度数给出了应用vmult()操作所要进行的矩阵-向量乘积的数量。
       *
       */
      unsigned int degree;

      /**
       * 这设置了矩阵的最大特征值，为了使切比雪夫预处理程序有适当的性能，需要适当设置。
       *
       */
      double max_eigenvalue;

      /**
       * 这设置了最大特征值和最小特征值之间的比率。
       *
       */
      double eigenvalue_ratio;

      /**
       * 这设置了最小特征值，这是一个可选的参数，只在内部用于检查我们是否使用了一个身份矩阵。
       *
       */
      double min_eigenvalue;

      /**
       * 这设置了一个阈值，低于这个阈值的对角线元素将不会在切比雪夫算法中被倒置。
       *
       */
      double min_diagonal;

      /**
       * 当这个标志被设置为<tt>true</tt>时，它使方法<tt>vmult(dst,
       * src)</tt>使用向量<tt>dst</tt>中的非零数据，将切比雪夫修正附加到它。这在某些情况下是有用的（例如用于高频误差平滑时），但不是求解器类期望的预处理工作方式（在预处理应用中忽略<tt>dst</tt>中的内容）。用户在触摸这个标志时，应该真正知道他们在做什么。
       *
       */
      bool nonzero_starting;
    };

    /**
     * 初始化函数。接受用于形成预处理程序的矩阵，如果有额外的标志，则接受额外的标志。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());
  };



  /**
   * 这个类实现了一个基于Trilinos
   * ML实现的代数多网格（AMG）预处理，它是一个黑盒预处理，对许多基于PDE的线性问题很有效。
   * 这个类的作用有两个方面。
   * 当调用initialize()函数时，基于我们希望预处理器所基于的矩阵，创建一个ML预处理器对象。调用相应的
   * <code>vmult</code>
   * 函数确实会调用Trilinos包中的相应操作，在那里它被称为
   * <code>ApplyInverse</code>  。这个类的使用在 step-31
   * 的教程程序中解释。
   * 由于我们要使用的Trilinos对象在很大程度上依赖于Epetra对象，我们建议将这个类与Trilinos（Epetra）稀疏矩阵和向量结合起来使用。也支持与
   * dealii::SparseMatrix
   * 类的矩阵和相应的向量一起使用，但这需要生成矩阵的副本，这样会比较慢，而且占用（更多）内存。在进行这样的复制操作时，我们仍然可以从预处理矩阵中的一些条目为零的事实中获益，因此可以忽略不计。
   * 该实现能够区分椭圆问题的矩阵和对流主导的问题。我们使用Trilinos
   * ML为椭圆问题提供的标准选项，但我们使用Chebyshev平滑器而不是对称的Gauss-Seidel平滑器。
   * 对于大多数椭圆问题，Chebyshev比Gauss-Seidel（SSOR）提供了更好的高频阻尼（在代数意义上），而且速度更快（Chebyshev只需要一些矩阵-向量乘积，而SSOR需要进行替换，费用更高）。此外，Chebyshev是完全并行的，因为它在许多处理器上使用时不会退化。另一方面，SSOR在许多处理器上变得更像雅可比。
   * 为了使这个类的功能正常，我们建议使用Trilinos
   * v9.0及以上版本。旧版本在使用每行有许多非零项的矩阵（即源于高阶有限元离散化的矩阵）时，可能在生成粗略的矩阵结构方面有问题。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionAMG : public PreconditionBase
  {
  public:
    /**
     * 一个数据结构，用于控制代数多重网格的设置细节。这里详述的标志会被传递给Trilinos
     * ML实现。一个当前类型的结构被传递给PreconditionAMG的构造函数。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，我们假装在标量方程上用线性有限元处理椭圆问题。
       * 利用 DoFTools::extract_constant_modes() 函数， @p constant_modes
       * 向量可以按以下方式对给定的场进行初始化。
       * @code
       * #include <deal.II/dofs/dof_tools.h>
       * ...
       *
       * DoFHandler<...> dof_handler;
       * FEValuesExtractors::Type... field_extractor;
       * ...
       *
       * TrilinosWrappers::PreconditionAMG::AdditionalData data;
       * DoFTools::extract_constant_modes(
       *   dof_handler,
       *   dof_handler.get_fe_collection().component_mask(field_extractor),
       *   data.constant_modes );
       * @endcode
       *
       *
       */
      AdditionalData(const bool         elliptic              = true,
                     const bool         higher_order_elements = false,
                     const unsigned int n_cycles              = 1,
                     const bool         w_cyle                = false,
                     const double       aggregation_threshold = 1e-4,
                     const std::vector<std::vector<bool>> &constant_modes =
                       std::vector<std::vector<bool>>(0),
                     const unsigned int smoother_sweeps  = 2,
                     const unsigned int smoother_overlap = 0,
                     const bool         output_details   = false,
                     const char *       smoother_type    = "Chebyshev",
                     const char *       coarse_type      = "Amesos-KLU");

      /**
       * 填入一个 @p parameter_list
       * ，可用于初始化AMG预处理器。             @p matrix 与 @p
       * constant_modes
       * 结合使用，用于配置预处理程序的无效空间设置。
       * @p distributed_constant_modes 由该函数初始化，在
       * PreconditionAMG::initialize() 被调用之前必须保持其范围。
       * @note
       * 设置参数反映了这个对象中的当前设置，各种选项既可以通过成员变量的状态直接设置（例如
       * "平滑器：类型"），也可以间接设置（例如
       * "聚合：类型"）。如果你希望对AMG预处理程序的配置进行精细控制，那么你可以使用这个函数创建参数列表（方便地设置运算符的空位），改变相关设置，并使用修正后的参数列表作为
       * PreconditionAMG::initialize(),
       * 的参数，而不是AdditionalData对象本身。
       * 关于有哪些选项可供修改，请参见<a
       * href="https://trilinos.org/docs/dev/packages/ml/doc/html/index.html">
       * Trilinos ML package</a>的文档。
       * @note
       * 任何与此数据结构设置的参数不冲突的用户定义的参数将被保留。
       *
       */
      void
      set_parameters(
        Teuchos::ParameterList &             parameter_list,
        std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
        const Epetra_RowMatrix &             matrix) const;

      /**
       * 填入一个参数列表，可用于初始化AMG预处理程序。
       * @note
       * 任何与此数据结构设置的参数不冲突的用户定义的参数将被保留。
       *
       */
      void
      set_parameters(
        Teuchos::ParameterList &             parameter_list,
        std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
        const SparseMatrix &                 matrix) const;

      /**
       * 根据 @p constant_modes 变量的状态，为输入 @p matrix 配置
       * @p parameter_list 中的无效空间设置。
       *
       */
      void
      set_operator_null_space(
        Teuchos::ParameterList &             parameter_list,
        std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
        const Epetra_RowMatrix &             matrix) const;

      /**
       * 根据 @p constant_modes 变量的状态，为输入 @p matrix 配置
       * @p parameter_list 中的无效空间设置。
       *
       */
      void
      set_operator_null_space(
        Teuchos::ParameterList &             parameter_list,
        std::unique_ptr<Epetra_MultiVector> &distributed_constant_modes,
        const SparseMatrix &                 matrix) const;

      /**
       * 决定AMG预处理是否应该针对椭圆问题（ML选项平滑聚合SA，使用切比雪夫平滑器）或非椭圆问题（ML选项非对称平滑聚合NSSA，平滑器为SSOR与欠放松）进行优化。
       *
       */
      bool elliptic;

      /**
       * 决定预处理程序所依据的矩阵是由线性元素还是高阶元素生成。
       *
       */
      bool higher_order_elements;

      /**
       * 定义预处理程序应执行多少个多重网格循环。
       *
       */
      unsigned int n_cycles;

      /**
       * 定义是否应使用w循环而不是标准设置的v循环。
       *
       */
      bool w_cycle;

      /**
       * 这个阈值告诉AMG设置应该如何进行粗化。在ML使用的AMG中，所有与暂定的粗放级点强烈耦合的点形成一个集合体。术语
       * <em> 强耦合 </em>
       * 由变量<tt>aggregation_threshold</tt>控制，意味着所有不小于<tt>aggregation_threshold</tt>的对角线元素都做强耦合。
       *
       */
      double aggregation_threshold;

      /**
       * 指定矩阵的恒定模式（近空空间）。这个参数告诉AMG我们是在标量方程（近空空间只由1组成，默认值为OK）还是在矢量值方程上工作。对于有<tt>n_component</tt>的向量值方程问题，提供的 @p constant_modes 应满足以下要求。        <ul>   <li>  n_component.size() == <tt>n_component</tt>  </li>   <li>  n_component[*].size() == n_dof_local 或者 n_component[*].size() == n_dof_global  </li>  ]  <li>  n_component[<tt>ic</tt>][<tt>id</tt>] =="<tt>id</tt>  <em>  th  </em>  DoF是对应于组件<tt>ic</tt>  </li>   </ul>  的。
       *
       */
      std::vector<std::vector<bool>> constant_modes;

      /**
       * 决定应该进行多少次平滑器的扫频。当标志<tt>elliptic</tt>被设置为<tt>true</tt>，即对于椭圆或几乎椭圆的问题，切比雪夫平滑器的多项式程度被设置为<tt>smoother_sweeps</tt>。扫频指的是在切比雪夫情况下进行的矩阵-向量乘积的数量。在非椭圆情况下，<tt>smoother_sweeps</tt>设置后平滑的SSOR松弛扫频次数。
       *
       */
      unsigned int smoother_sweeps;

      /**
       * 确定并行运行时SSOR/Chebyshev误差平滑器的重叠度。
       *
       */
      unsigned int smoother_overlap;

      /**
       * 如果这个标志被设置为<tt>true</tt>，那么来自ML预处理程序的内部信息将被打印到屏幕上。这在调试预处理程序时很有用。
       *
       */
      bool output_details;

      /**
       * 决定在AMG循环中使用哪个平滑器。smoother_type的可能性有以下几种。        <ul>   <li>  "Aztec"  </li>   <li>  "IFPACK"  </li>   <li>  "Jacobi"  </li>   <li>  ] "ML对称高斯-塞德尔"  </li>   <li>  "对称高斯-塞德尔"  </li>   <li>  "ML高斯-塞德尔"  </li>   <li>  ] "高斯-赛德尔"  </li>   <li>  "高斯-赛德尔块"  </li>   <li>  "对称块高斯-赛德尔"  </li>   <li>  ] "切比雪夫"  </li>   <li>  "MLS"  </li>   <li>  "Hiptmair"  </li>   <li>  "Amesos-KLU"  </li>   <li>  ] "Amesos-Superlu"  </li>   <li>  "Amesos-UMFPACK"  </li>   <li>  "Amesos-Superludist"  </li>   <li>  ] "Amesos-MUMPS"  </li>   <li>  "用户定义"  </li>   <li>  "SuperLU"  </li>   <li>  "IFPACK-Chebyshev"  </li>   <li>  ] "自己"  </li>   <li>  "无为"  </li>   <li>  "IC"  </li>   <li>  "ICT"  </li>   <li>  "ILU"  </li>  ]  <li>  "ILUT"  </li>   <li>  "Block Chebyshev"  </li>   <li>  "IFPACK-Block Chebyshev"  </li>   </ul>
       *
       */
      const char *smoother_type;

      /**
       * 决定在最粗糙的层次上使用哪种求解器。可以使用与平滑类型相同的设置。
       *
       */
      const char *coarse_type;
    };

    /**
     * 解构器。
     *
     */
    ~PreconditionAMG() override;


    /**
     * 让Trilinos为给定矩阵的线性系统的解计算一个多层次的层次。该函数使用
     * TrilinosWrappers::SparseMatrix. 中指定的矩阵格式。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());

    /**
     * 让Trilinos为给定矩阵的线性系统的解计算一个多级层次结构。与上面的其他初始化函数不同，这个函数使用了一个抽象的接口，指向Epetra_RowMatrix类型的对象，允许用户向ML预处理程序传递相当一般的对象。
     * 这个初始化例程在需要预处理的运算符不是
     * TrilinosWrappers::SparseMatrix
     * 对象的情况下非常有用，但仍然可以获得本地拥有的矩阵行的每个条目的副本（方法ExtractMyRowCopy），并实现矩阵-向量乘积（方法Multiply或Apply）。一个例子是提供比使用矩阵条目更快的矩阵-向量乘法的运算符（无矩阵方法）。这些实现可以与只执行矩阵-向量乘积的切比雪夫平滑器有益地结合起来。接口类Epetra_RowMatrix非常灵活，可以实现这种实现方式。
     *
     */
    void
    initialize(const Epetra_RowMatrix &matrix,
               const AdditionalData &  additional_data = AdditionalData());

    /**
     * 让Trilinos为给定矩阵的线性系统的解计算一个多层次的层次。该函数使用
     * TrilinosWrappers::SparseMatrix.
     * 中指定的矩阵格式，该函数与上述函数类似，但允许用户设置Trilinos
     * ML预处理程序的所有选项。为了了解ML的所有选项，我们参考了<a
     * href="https://trilinos.org/docs/dev/packages/ml/doc/html/index.html">ML
     * documentation
     * page</a>。特别是在需要解决矢量值问题的情况下，用户需要遵循ML指令。
     *
     */
    void
    initialize(const SparseMatrix &          matrix,
               const Teuchos::ParameterList &ml_parameters);

    /**
     * 让Trilinos计算一个多级层次的线性系统的解决方案，并给出矩阵。与上面的其他初始化函数不同，这个函数使用了一个抽象的接口，指向Epetra_RowMatrix类型的对象，允许用户向ML预处理程序传递相当一般的对象。
     *
     */
    void
    initialize(const Epetra_RowMatrix &      matrix,
               const Teuchos::ParameterList &ml_parameters);

    /**
     * 让Trilinos为给定矩阵的线性系统的解计算一个多级层次结构。这个函数需要一个deal.II矩阵，并将其内容复制到Trilinos矩阵中，所以这个函数可以说是相当低效的。
     *
     */
    template <typename number>
    void
    initialize(const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
               const AdditionalData &additional_data = AdditionalData(),
               const double          drop_tolerance  = 1e-13,
               const ::dealii::SparsityPattern *use_this_sparsity = nullptr);

    /**
     * 当预处理器的基础矩阵条目发生变化，但矩阵的稀疏模式保持不变时，该函数可用于更快地重新计算预处理器的构造。这个函数的作用是利用已经生成的粗化结构，根据平滑聚合策略计算AMG延长和限制，然后建立整个多级层次结构。这个函数可以比初始化函数快得多，因为粗化模式通常是设置AMG
     * ML预处理程序时最困难的事情。
     *
     */
    void
    reinit();

    /**
     * 销毁预处理程序，留下一个像刚刚调用构造函数后的对象。
     *
     */
    void
    clear();

    /**
     * 打印该类的内存消耗估计值。
     *
     */
    size_type
    memory_consumption() const;

  private:
    /**
     * 将deal.II矩阵复制成Trilinos格式。
     *
     */
    std::shared_ptr<SparseMatrix> trilinos_matrix;
  };



#    if defined(DOXYGEN) || defined(DEAL_II_TRILINOS_WITH_MUELU)
  /**
   * 该类实现了基于Trilinos
   * MueLu实现的代数多网格（AMG）预处理，这是一个黑盒预处理，对许多基于PDE的线性问题都很有效。PreconditionerAMGMueLu的接口与PreconditionerAMG的接口相同，只是Higher_order_elements参数在PreconditionerAMGMueLu中不存在。
   * @note
   * 你需要配置支持MueLU的Trilinos才能使这个预处理程序工作。
   * @note  目前不支持64位指数。      @warning
   * 这个接口不应该被认为是稳定的。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionAMGMueLu : public PreconditionBase
  {
  public:
    /**
     * 一个数据结构，用于控制代数多重网格的设置细节。这里详述的标志被传递给Trilinos
     * MueLu的实现。当前类型的结构被传递给PreconditionAMGMueLu的构造函数。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。默认情况下，我们假装在标量方程上用线性有限元处理椭圆问题。
       *
       */
      AdditionalData(const bool         elliptic              = true,
                     const unsigned int n_cycles              = 1,
                     const bool         w_cyle                = false,
                     const double       aggregation_threshold = 1e-4,
                     const std::vector<std::vector<bool>> &constant_modes =
                       std::vector<std::vector<bool>>(0),
                     const unsigned int smoother_sweeps  = 2,
                     const unsigned int smoother_overlap = 0,
                     const bool         output_details   = false,
                     const char *       smoother_type    = "Chebyshev",
                     const char *       coarse_type      = "Amesos-KLU");

      /**
       * 决定AMG预处理是否应该针对椭圆问题（MueLu选项平滑聚合SA，使用Chebyshev平滑器）或非椭圆问题（MueLu选项非对称平滑聚合NSSA，平滑器为SSOR与欠放松）进行优化。
       *
       */
      bool elliptic;

      /**
       * 定义预处理程序应执行多少个多重网格循环。
       *
       */
      unsigned int n_cycles;

      /**
       * 定义是否应使用w循环而不是标准设置的v循环。
       *
       */
      bool w_cycle;

      /**
       * 这个阈值告诉AMG设置应该如何进行粗化。在MueLu使用的AMG中，所有与暂定粗略级点强烈耦合的点形成一个集合体。术语
       * <em> 强耦合 </em>
       * 由变量<tt>aggregation_threshold</tt>控制，意味着所有不小于<tt>aggregation_threshold</tt>的对角线元素都做强耦合。
       *
       */
      double aggregation_threshold;

      /**
       * 指定矩阵的恒定模式（近空空间）。这个参数告诉AMG我们是在标量方程（近空空间只由1组成）还是在矢量值方程上工作。
       *
       */
      std::vector<std::vector<bool>> constant_modes;

      /**
       * 决定应该执行多少次平滑器扫频。当标志<tt>elliptic</tt>被设置为<tt>true</tt>，即对于椭圆或几乎椭圆的问题，切比雪夫平滑器的多项式程度被设置为<tt>smoother_sweeps</tt>。扫频指的是在切比雪夫情况下进行的矩阵-向量乘积的数量。在非椭圆情况下，<tt>smoother_sweeps</tt>设置后平滑的SSOR松弛扫频次数。
       *
       */
      unsigned int smoother_sweeps;

      /**
       * 确定并行运行时SSOR/Chebyshev误差平滑器的重叠度。
       *
       */
      unsigned int smoother_overlap;

      /**
       * 如果这个标志被设置为<tt>true</tt>，那么来自ML预处理程序的内部信息将被打印到屏幕上。这在调试预处理程序时很有用。
       *
       */
      bool output_details;

      /**
       * 决定在AMG循环中使用哪个平滑器。smoother_type的可能性有以下几种。        <ul>   <li>  "Aztec"  </li>   <li>  "IFPACK"  </li>   <li>  "Jacobi"  </li>   <li>  ] "ML对称高斯-塞德尔"  </li>   <li>  "对称高斯-塞德尔"  </li>   <li>  "ML高斯-塞德尔"  </li>   <li>  ] "高斯-赛德尔"  </li>   <li>  "高斯-赛德尔块"  </li>   <li>  "对称块高斯-赛德尔"  </li>   <li>  ] "切比雪夫"  </li>   <li>  "MLS"  </li>   <li>  "Hiptmair"  </li>   <li>  "Amesos-KLU"  </li>   <li>  ] "Amesos-Superlu"  </li>   <li>  "Amesos-UMFPACK"  </li>   <li>  "Amesos-Superludist"  </li>   <li>  ] "Amesos-MUMPS"  </li>   <li>  "用户定义"  </li>   <li>  "SuperLU"  </li>   <li>  "IFPACK-Chebyshev"  </li>   <li>  ] "自我"  </li>   <li>  "无为"  </li>   <li>  "IC"  </li>   <li>  "ICT"  </li>   <li>  "ILU"  </li>  ]  <li>  "ILUT"  </li>   <li>  "Block Chebyshev"  </li>   <li>  "IFPACK-Block Chebyshev"  </li>   </ul>
       *
       */
      const char *smoother_type;

      /**
       * 决定在最粗糙的层次上使用哪种求解器。可以使用与平滑类型相同的设置。
       *
       */
      const char *coarse_type;
    };

    /**
     * 构造函数。
     *
     */
    PreconditionAMGMueLu();

    /**
     * 销毁器。
     *
     */
    virtual ~PreconditionAMGMueLu() override = default;

    /**
     * 让Trilinos计算一个多级层次的线性系统的解决方案，并给出矩阵。该函数使用
     * TrilinosWrappers::SparseMatrix. 中指定的矩阵格式。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());

    /**
     * 让Trilinos为给定矩阵的线性系统的解计算一个多级层次结构。与上面的其他初始化函数不同，这个函数使用一个Epetra_CrsMatrixCrs类型的对象。
     *
     */
    void
    initialize(const Epetra_CrsMatrix &matrix,
               const AdditionalData &  additional_data = AdditionalData());

    /**
     * 让Trilinos计算一个多级层次的线性系统的解决方案，给定的矩阵。该函数使用
     * TrilinosWrappers::SparseMatrix. 中指定的矩阵格式
     * 该函数与上述函数类似，但允许用户设置Trilinos
     * ML预处理程序的大部分选项。为了了解ML的所有选项，我们参考了<a
     * href="https://trilinos.org/docs/dev/packages/ml/doc/html/index.html">ML
     * documentation
     * page</a>。不是所有的ML选项都有相应的MueLu选项。
     *
     */
    void
    initialize(const SparseMatrix &    matrix,
               Teuchos::ParameterList &muelu_parameters);

    /**
     * 让Trilinos为给定矩阵的线性系统的解计算出一个多级层次结构。与上面的其他初始化函数不同，这个函数使用一个Epetra_CrsMatrix类型的对象。
     *
     */
    void
    initialize(const Epetra_CrsMatrix &matrix,
               Teuchos::ParameterList &muelu_parameters);

    /**
     * 让Trilinos计算一个多级层次的线性系统的解决方案，给定的矩阵。这个函数接收一个deal.ii矩阵，并将其内容复制到Trilinos矩阵中，所以这个函数可以说是相当低效的。
     *
     */
    template <typename number>
    void
    initialize(const ::dealii::SparseMatrix<number> &deal_ii_sparse_matrix,
               const AdditionalData &additional_data = AdditionalData(),
               const double          drop_tolerance  = 1e-13,
               const ::dealii::SparsityPattern *use_this_sparsity = nullptr);

    /**
     * 销毁预处理程序，留下一个像刚刚调用构造函数后的对象。
     *
     */
    void
    clear();

    /**
     * 打印这个类的内存消耗估计值。
     *
     */
    size_type
    memory_consumption() const;

  private:
    /**
     * 将deal.II矩阵复制成Trilinos格式。
     *
     */
    std::shared_ptr<SparseMatrix> trilinos_matrix;
  };
#    endif



  /**
   * 一个用于Trilinos矩阵的身份预处理的封装类。
   * @ingroup TrilinosWrappers
   * @ingroup Preconditioners
   *
   */
  class PreconditionIdentity : public PreconditionBase
  {
  public:
    /**
     * 这个函数的出现只是为了提供一个预处理程序的接口，以交给平滑器。
     * 这个功能什么都不做。
     *
     */
    struct AdditionalData
    {};

    /**
     * 矩阵参数被忽略，这里只是为了与更复杂的预处理程序兼容。
     * @note
     * 当这个预处理程序要被包裹在一个没有典范母体的LinearOperator中时，必须调用这个函数。
     *
     */
    void
    initialize(const SparseMatrix &  matrix,
               const AdditionalData &additional_data = AdditionalData());

    /**
     * 应用预处理程序，即dst = src。
     *
     */
    void
    vmult(MPI::Vector &dst, const MPI::Vector &src) const override;

    /**
     * 应用传输调节器，即dst = src。
     *
     */
    void
    Tvmult(MPI::Vector &dst, const MPI::Vector &src) const override;

    /**
     * 在deal.II数据结构上应用前置条件器，而不是Trilinos包装类中提供的数据结构，即dst
     * = src。
     *
     */
    void
    vmult(dealii::Vector<double> &      dst,
          const dealii::Vector<double> &src) const override;

    /**
     * 对deal.II数据结构应用转置预处理，而不是Trilinos包装类中提供的数据结构，即dst
     * = src。
     *
     */
    void
    Tvmult(dealii::Vector<double> &      dst,
           const dealii::Vector<double> &src) const override;

    /**
     * 在deal.II并行数据结构上应用预处理程序，而不是Trilinos包装类中提供的数据结构，即dst
     * = src。
     *
     */
    void
    vmult(LinearAlgebra::distributed::Vector<double> &              dst,
          const dealii::LinearAlgebra::distributed::Vector<double> &src)
      const override;

    /**
     * 在deal.II并行数据结构上应用转置预处理程序，而不是Trilinos包装类中提供的数据结构，即dst
     * = src。
     *
     */
    void
    Tvmult(LinearAlgebra::distributed::Vector<double> &              dst,
           const dealii::LinearAlgebra::distributed::Vector<double> &src)
      const override;
  };



  // ----------------------- inline and template functions --------------------


#    ifndef DOXYGEN


  inline void
  PreconditionBase::transpose()
  {
    // This only flips a flag that tells
    // Trilinos that any vmult operation
    // should be done with the
    // transpose. However, the matrix
    // structure is not reset.
    int ierr;

    if (!preconditioner->UseTranspose())
      {
        ierr = preconditioner->SetUseTranspose(true);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
    else
      {
        ierr = preconditioner->SetUseTranspose(false);
        AssertThrow(ierr == 0, ExcTrilinosError(ierr));
      }
  }


  inline void
  PreconditionBase::vmult(MPI::Vector &dst, const MPI::Vector &src) const
  {
    Assert(dst.trilinos_partitioner().SameAs(
             preconditioner->OperatorRangeMap()),
           ExcNonMatchingMaps("dst"));
    Assert(src.trilinos_partitioner().SameAs(
             preconditioner->OperatorDomainMap()),
           ExcNonMatchingMaps("src"));

    const int ierr = preconditioner->ApplyInverse(src.trilinos_vector(),
                                                  dst.trilinos_vector());
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }

  inline void
  PreconditionBase::Tvmult(MPI::Vector &dst, const MPI::Vector &src) const
  {
    Assert(dst.trilinos_partitioner().SameAs(
             preconditioner->OperatorRangeMap()),
           ExcNonMatchingMaps("dst"));
    Assert(src.trilinos_partitioner().SameAs(
             preconditioner->OperatorDomainMap()),
           ExcNonMatchingMaps("src"));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse(src.trilinos_vector(),
                                                  dst.trilinos_vector());
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }


  // For the implementation of the <code>vmult</code> function with deal.II
  // data structures we note that invoking a call of the Trilinos
  // preconditioner requires us to use Epetra vectors as well. We do this by
  // providing a view, i.e., feed Trilinos with a pointer to the data, so we
  // avoid copying the content of the vectors during the iteration (this
  // function is only useful when used in serial anyway). In the declaration
  // of the right hand side, we need to cast the source vector (that is
  // <code>const</code> in all deal.II calls) to non-constant value, as this
  // is the way Trilinos wants to have them.
  inline void
  PreconditionBase::vmult(dealii::Vector<double> &      dst,
                          const dealii::Vector<double> &src) const
  {
    AssertDimension(dst.size(),
                    preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension(src.size(),
                    preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst(View,
                           preconditioner->OperatorDomainMap(),
                           dst.begin());
    Epetra_Vector tril_src(View,
                           preconditioner->OperatorRangeMap(),
                           const_cast<double *>(src.begin()));

    const int ierr = preconditioner->ApplyInverse(tril_src, tril_dst);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }


  inline void
  PreconditionBase::Tvmult(dealii::Vector<double> &      dst,
                           const dealii::Vector<double> &src) const
  {
    AssertDimension(dst.size(),
                    preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension(src.size(),
                    preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst(View,
                           preconditioner->OperatorDomainMap(),
                           dst.begin());
    Epetra_Vector tril_src(View,
                           preconditioner->OperatorRangeMap(),
                           const_cast<double *>(src.begin()));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse(tril_src, tril_dst);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }



  inline void
  PreconditionBase::vmult(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    AssertDimension(dst.locally_owned_size(),
                    preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension(src.locally_owned_size(),
                    preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst(View,
                           preconditioner->OperatorDomainMap(),
                           dst.begin());
    Epetra_Vector tril_src(View,
                           preconditioner->OperatorRangeMap(),
                           const_cast<double *>(src.begin()));

    const int ierr = preconditioner->ApplyInverse(tril_src, tril_dst);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
  }

  inline void
  PreconditionBase::Tvmult(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const
  {
    AssertDimension(dst.locally_owned_size(),
                    preconditioner->OperatorDomainMap().NumMyElements());
    AssertDimension(src.locally_owned_size(),
                    preconditioner->OperatorRangeMap().NumMyElements());
    Epetra_Vector tril_dst(View,
                           preconditioner->OperatorDomainMap(),
                           dst.begin());
    Epetra_Vector tril_src(View,
                           preconditioner->OperatorRangeMap(),
                           const_cast<double *>(src.begin()));

    preconditioner->SetUseTranspose(true);
    const int ierr = preconditioner->ApplyInverse(tril_src, tril_dst);
    AssertThrow(ierr == 0, ExcTrilinosError(ierr));
    preconditioner->SetUseTranspose(false);
  }

#    endif

} // namespace TrilinosWrappers


 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_TRILINOS

#endif // trilinos_precondition_h
 /*------------------------- trilinos_precondition.h -------------------------*/ 


