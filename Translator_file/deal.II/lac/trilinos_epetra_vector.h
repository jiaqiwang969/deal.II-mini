//include/deal.II-translator/lac/trilinos_epetra_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#ifndef dealii_trilinos_epetra_vector_h
#define dealii_trilinos_epetra_vector_h


#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)

#  include <deal.II/base/index_set.h>
#  include <deal.II/base/subscriptor.h>

#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/vector_operation.h>
#  include <deal.II/lac/vector_space_vector.h>
#  include <deal.II/lac/vector_type_traits.h>

#  include <Epetra_FEVector.h>
#  include <mpi.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  // Forward declaration
  template <typename Number>
  class ReadWriteVector;

  /**
   * 一个命名空间，用于为Tridinos的Epetra向量提供包装器的类。
   * 该命名空间为Tridinos的Epetra包（https://trilinos.github.io/epetra.html）中的Epetra_FEVector类提供包装器。
   *
   */
  namespace EpetraWrappers
  {
    /**
     * 该类实现了对Trilinos分布式向量类Epetra_FEVector的包装器。这个类是由
     * LinearAlgebra::VectorSpaceVector
     * 类派生的。但是请注意，Epetra只适用于Number =
     * double。这个类要求Trilinos在编译时支持MPI。
     * @ingroup TrilinosWrappers
     * @ingroup Vectors
     *
     */
    class Vector : public VectorSpaceVector<double>, public Subscriptor
    {
    public:
      /**
       * 构造函数。创建一个维数为0的向量。
       *
       */
      Vector();

      /**
       * 复制构造函数。将维数和分区设置为给定的向量，并复制所有的元素。
       *
       */
      Vector(const Vector &V);

      /**
       * 这个构造函数接收一个IndexSet，它定义了如何在MPI处理器之间分配各个组件。由于它还包括关于向量大小的信息，这就是我们生成一个%并行向量所需要的全部内容。
       *
       */
      explicit Vector(const IndexSet &parallel_partitioner,
                      const MPI_Comm &communicator);

      /**
       * Reinit功能。这个功能销毁了旧的向量内容，并根据输入的分区生成一个新的向量。标志<tt>omit_zeroing_entries</tt>决定了向量是否应该被填入零（false）或不被触动（true）。
       *
       */
      void
      reinit(const IndexSet &parallel_partitioner,
             const MPI_Comm &communicator,
             const bool      omit_zeroing_entries = false);

      /**
       * 将尺寸改为向量 @p V. 的尺寸， @p V
       * 的元素不会被复制。
       *
       */
      virtual void
      reinit(const VectorSpaceVector<double> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * 复制函数。这个函数接收一个Vector并复制所有的元素。该向量将具有与
       * @p  V相同的平行分布。
       *
       */
      Vector &
      operator=(const Vector &V);

      /**
       * 将向量的所有元素设置为标量 @p s. ，只有在 @p s
       * 等于零的情况下才允许这个操作。
       *
       */
      virtual Vector &
      operator=(const double s) override;

      /**
       * 从输入向量 @p V.  VectorOperation::values  @p operation
       * 中导入向量的IndexSet中存在的所有元素，用于决定 @p
       * V
       * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
       *
       */
      virtual void
      import(const ReadWriteVector<double> &V,
             VectorOperation::values        operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {}) override;

      /**
       * 将整个向量乘以一个固定系数。
       *
       */
      virtual Vector &
      operator*=(const double factor) override;

      /**
       * 将整个向量除以一个固定的因子。
       *
       */
      virtual Vector &
      operator/=(const double factor) override;

      /**
       * 将向量 @p V 添加到现在的向量中。
       *
       */
      virtual Vector &
      operator+=(const VectorSpaceVector<double> &V) override;

      /**
       * 从现在的向量中减去向量 @p V 。
       *
       */
      virtual Vector &
      operator-=(const VectorSpaceVector<double> &V) override;

      /**
       * 返回两个向量的标量乘积。这些向量需要有相同的布局。
       *
       */
      virtual double
      operator*(const VectorSpaceVector<double> &V) const override;

      /**
       * 将 @p a 添加到所有分量中。注意， @p is
       * 是一个标量，不是一个向量。
       *
       */
      virtual void
      add(const double a) override;

      /**
       * 一个向量的倍数的简单加法，即<tt>*this +=
       * a*V</tt>。向量需要有相同的布局。
       *
       */
      virtual void
      add(const double a, const VectorSpaceVector<double> &V) override;

      /**
       * 一个向量的多个加法，即：<tt>*this> +=
       * a*V+b*W</tt>。向量需要有相同的布局。
       *
       */
      virtual void
      add(const double                     a,
          const VectorSpaceVector<double> &V,
          const double                     b,
          const VectorSpaceVector<double> &W) override;

      /**
       * 缩放和一个向量的倍数的简单相加，即<tt>*this=s*(*this)+a*V</tt>。
       *
       */
      virtual void
      sadd(const double                     s,
           const double                     a,
           const VectorSpaceVector<double> &V) override;

      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。这些向量需要有相同的布局。
       *
       */
      virtual void
      scale(const VectorSpaceVector<double> &scaling_factors) override;

      /**
       * 赋值 <tt>*this = a*V</tt>.
       *
       */
      virtual void
      equ(const double a, const VectorSpaceVector<double> &V) override;

      /**
       * 返回向量是否只包含值为0的元素。
       *
       */
      virtual bool
      all_zero() const override;

      /**
       * 返回这个向量的元素的平均值。
       *
       */
      virtual double
      mean_value() const override;

      /**
       * 返回该向量的l<sub>1</sub>准则（即所有处理器中所有条目的绝对值之和）。
       *
       */
      virtual double
      l1_norm() const override;

      /**
       * 返回向量的l<sub>2</sub>准则（即所有处理器中所有条目的平方之和的平方根）。
       *
       */
      virtual double
      l2_norm() const override;

      /**
       * 返回向量的最大规范（即所有条目和所有处理器之间的最大绝对值）。
       *
       */
      virtual double
      linfty_norm() const override;

      /**
       * 执行一个向量加法和随后的内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
       * @code
       * this->add(a, V);
       * return_value =this W;
       * @endcode
       * 这个函数存在的原因是这个操作比单独调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
       * @p this,   @p V,   @p W,
       * ，而调用单独的方法意味着要加载调用向量 @p this
       * 两次。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
       * @p W 等于 @p this). ，则为50\%）。
       * 向量需要有相同的布局。
       * 对于复值向量，第二步中的标量乘积被实现为
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
       *
       */
      virtual double
      add_and_dot(const double                     a,
                  const VectorSpaceVector<double> &V,
                  const VectorSpaceVector<double> &W) override;
      /**
       * 这个函数总是返回false，它的存在只是为了向后兼容。
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 返回向量的全局大小，等于所有处理器中本地拥有的索引数之和。
       *
       */
      virtual size_type
      size() const override;

      /**
       * 返回向量的本地大小，即本地拥有的索引数。
       *
       */
      size_type
      locally_owned_size() const;

      /**
       * 返回与此对象一起使用的MPI通信器对象。
       *
       */
      MPI_Comm
      get_mpi_communicator() const;

      /**
       * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相干的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
       * @code
       * vec.locally_owned_elements() == complete_index_set(vec.size())
       * @endcode
       *
       *
       */
      virtual ::dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * 返回一个对底层Trilinos Epetra_FEVector类的常量引用。
       *
       */
      const Epetra_FEVector &
      trilinos_vector() const;

      /**
       * 返回一个（可修改的）对底层Tridinos
       * Epetra_FEVector类的引用。
       *
       */
      Epetra_FEVector &
      trilinos_vector();

      /**
       * 将向量打印到输出流中  @p out.  。
       *
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * 以字节为单位返回这个类的内存消耗。
       *
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * 向量有不同的分区，即它们的IndexSet对象不代表相同的指数。
       *
       */
      DeclException0(ExcDifferentParallelPartitioning);

      /**
       * 试图在两个不兼容的向量类型之间执行操作。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * 在Trilinos中被错误抛出的异常。
       * @ingroup Exceptions
       *
       */
      DeclException1(ExcTrilinosError,
                     int,
                     << "An error with error number " << arg1
                     << " occurred while calling a Trilinos function");

    private:
      /**
       * 为IndexSet  @p source_index_set  和基于通信器 @p mpi_comm.
       * 的当前向量之间的通信创建通信模式。
       *
       */
      void
      create_epetra_comm_pattern(const IndexSet &source_index_set,
                                 const MPI_Comm &mpi_comm);

      /**
       * 指向实际Epetra向量对象的指针。
       *
       */
      std::unique_ptr<Epetra_FEVector> vector;

      /**
       * 最后导入的向量元素的索引集。
       *
       */
      ::dealii::IndexSet source_stored_elements;

      /**
       * Source_stored_elements
       * IndexSet和当前向量之间的通信模式。
       *
       */
      std::shared_ptr<const CommunicationPattern> epetra_comm_pattern;
    };


    inline bool
    Vector::has_ghost_elements() const
    {
      return false;
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra


/**
 * 将 dealii::LinearAlgebra::EpetraWrappers::Vector
 * 声明为分布式向量。
 *
 *
 */
template <>
struct is_serial_vector<LinearAlgebra::EpetraWrappers::Vector> : std::false_type
{};

DEAL_II_NAMESPACE_CLOSE

#endif

#endif


