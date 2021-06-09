//include/deal.II-translator/lac/la_parallel_block_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_la_parallel_block_vector_h
#define dealii_la_parallel_block_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declarations
#ifndef DOXYGEN
#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class BlockVector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class BlockVector;
  }
} // namespace TrilinosWrappers
#  endif
#endif

namespace LinearAlgebra
{
  namespace distributed
  {
    /*!   @addtogroup Vectors  
     * @{     
*
*/


    /**
     * 一个基于分布式deal.II向量的块向量的实现。虽然基类提供了大部分的接口，但这个类处理了向量的实际分配，并提供了底层向量类型的特定函数。
     * @note  这个模板的实例提供给<tt>  @<float@>  和  @<double@></tt>;  其他可以在应用程序中生成（见手册中的 @ref Instantiations 部分）。          @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    template <typename Number>
    class BlockVector : public BlockVectorBase<Vector<Number>>,
                        public VectorSpaceVector<Number>
    {
    public:
      /**
       * 在update_ghost_values()和compress()调用中分割通信的块大小。
       * 当有太多的消息/请求未完成时，大多数常见的MPI实现会变得很慢。即使消息很小，比如只有1
       * kB，我们也应该用 @p communication_block_size
       * 收集足够的数据，以覆盖典型的infiniband延迟，这些延迟大约是几微秒。在5GB/s的吞吐量下发送20
       * kB需要4微秒，所以我们应该达到带宽主导的制度，这已经很好了。
       *
       */
      static constexpr unsigned int communication_block_size = 20;

      /**
       * 对基类进行类型化定义，以便更简单地访问它自己的别名。
       *
       */
      using BaseClass = BlockVectorBase<Vector<Number>>;

      /**
       * 类型化底层向量的类型。
       *
       */
      using BlockType = typename BaseClass::BlockType;

      /**
       * 从基类中导入别名。
       *
       */
      using value_type      = typename BaseClass::value_type;
      using real_type       = typename BaseClass::real_type;
      using pointer         = typename BaseClass::pointer;
      using const_pointer   = typename BaseClass::const_pointer;
      using reference       = typename BaseClass::reference;
      using const_reference = typename BaseClass::const_reference;
      using size_type       = typename BaseClass::size_type;
      using iterator        = typename BaseClass::iterator;
      using const_iterator  = typename BaseClass::const_iterator;

      /**
       * @name  1: 基本操作
       *
       */
      //@{

      /**
       * 构造函数。有三种方法来使用这个构造函数。首先，没有任何参数，它生成一个没有块的对象。给定一个参数，它初始化<tt>num_blocks</tt>块，但这些块的大小为零。第三种变体最后将所有块初始化为相同的大小<tt>block_size</tt>。
       * 如果你打算使用不同大小的块，请参考下面的其他构造函数。
       *
       */
      explicit BlockVector(const size_type num_blocks = 0,
                           const size_type block_size = 0);

      /**
       * 复制-构造器。尺寸设置为V的尺寸，所有的组件都是从V中复制的。
       *
       */
      BlockVector(const BlockVector<Number> &V);

      /**
       * 拷贝构造函数获取另一种数据类型的BlockVector。如果没有从<tt>OtherNumber</tt>到<tt>Number</tt>的转换路径，这将失败。注意，当你复制到一个数据元素精度较低的BlockVector时，可能会失去精度。
       * 旧版本的gcc不尊重模板构造函数上的 @p explicit
       * 关键字。在这种情况下，很容易不小心写出效率很低的代码，因为编译器开始执行隐藏的转换。为了避免这种情况，如果我们在配置过程中检测到有损坏的编译器，这个功能就会被禁用。
       *
       */
      template <typename OtherNumber>
      explicit BlockVector(const BlockVector<OtherNumber> &v);

      /**
       * 构造函数。设置块的数量为<tt>block_sizes.size()</tt>，并以<tt>block_sizes[i]</tt>零元素初始化每个块。
       *
       */
      BlockVector(const std::vector<size_type> &block_sizes);

      /**
       * 构建一个块向量，用IndexSet来表示本地范围，并为每个块建立ghost条目。
       *
       */
      BlockVector(const std::vector<IndexSet> &local_ranges,
                  const std::vector<IndexSet> &ghost_indices,
                  const MPI_Comm &             communicator);

      /**
       * 与上述相同，但假设鬼魂索引为空。
       *
       */
      BlockVector(const std::vector<IndexSet> &local_ranges,
                  const MPI_Comm &             communicator);

      /**
       * 解构器。
       * @note
       * 我们需要明确地提供一个析构器，否则链接器可能会认为它是未使用的而丢弃它，尽管在不同的部分需要。英特尔的编译器很容易出现这种行为。
       *
       */
      virtual ~BlockVector() override = default;

      /**
       * 复制操作符：用给定的标量值填充向量的所有分量。
       *
       */
      virtual BlockVector &
      operator=(const value_type s) override;

      /**
       * 对于相同类型的参数的复制操作符。如果有必要，可以调整现在的向量的大小。
       *
       */
      BlockVector &
      operator=(const BlockVector &V);

      /**
       * 不同类型的模板参数的复制操作。如果有必要的话，调整当前向量的大小。
       *
       */
      template <class Number2>
      BlockVector &
      operator=(const BlockVector<Number2> &V);

      /**
       * 将一个普通向量复制到一个块向量中。
       *
       */
      BlockVector &
      operator=(const Vector<Number> &V);

#ifdef DEAL_II_WITH_PETSC
      /**
       * 将一个PETSc向量的内容复制到调用向量中。这个函数假定向量的布局已经被初始化为匹配。
       * 这个运算符只有在deal.II被配置为PETSc时才可用。
       *
       */
      BlockVector<Number> &
      operator=(const PETScWrappers::MPI::BlockVector &petsc_vec);
#endif

#ifdef DEAL_II_WITH_TRILINOS
      /**
       * 将一个Trilinos向量的内容复制到调用向量中。这个函数假定向量的布局已经被初始化为匹配。
       * 这个运算符只有在deal.II被配置为Trilinos时才可用。
       *
       */
      BlockVector<Number> &
      operator=(const TrilinosWrappers::MPI::BlockVector &trilinos_vec);
#endif

      /**
       * 重新初始化BlockVector，使其包含<tt>num_blocks</tt>每个大小为<tt>block_size</tt>的块。
       * 如果第二个参数是默认值，那么块向量就会分配指定数量的块，但让它们的大小为零。然后你需要重新初始化各个区块，并调用collect_sizes()来更新区块系统对其各个区块大小的认识。
       * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
       *
       */
      void
      reinit(const size_type num_blocks,
             const size_type block_size           = 0,
             const bool      omit_zeroing_entries = false);

      /**
       * 重新初始化BlockVector，使其包含<tt>block_sizes.size()</tt>块。每个块都被重新初始化为<tt>block_sizes[i]</tt>的尺寸。
       * 如果块的数量与调用此函数前相同，所有的向量都保持不变，并且为每个向量调用
       * reinit()。
       * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
       * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新其缓存的向量大小。如果你在其中一个块上调用reinit()，那么对这个对象的后续操作可能会产生不可预测的结果，因为它们可能会被路由到错误的块上。
       *
       */
      void
      reinit(const std::vector<size_type> &N,
             const bool                    omit_zeroing_entries = false);

      /**
       * 将尺寸改为向量<tt>V</tt>的尺寸。这与其他reinit()函数的情况相同。
       * <tt>V</tt>的元素不会被复制，也就是说，这个函数与调用<tt>reinit
       * (V.size(), omit_zeroing_entries)</tt>相同。
       * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新它的向量大小缓存。如果你调用其中一个块的reinit()，那么这个对象的后续操作可能会产生不可预测的结果，因为它们可能被路由到错误的块。
       *
       */
      template <typename Number2>
      void
      reinit(const BlockVector<Number2> &V,
             const bool                  omit_zeroing_entries = false);

      /**
       * 这个函数将积累在数据缓冲区中的鬼魂索引的数据复制到拥有的处理器中。关于参数 @p operation, 的含义，见术语表中的 @ref GlossCompress "压缩分布式向量和矩阵 "
       * 条目。            这个函数有两个变体。如果用参数 @p
       * 调用 VectorOperation::add
       * ，则将所有积累在幽灵元素中的数据添加到拥有处理器的相应元素中，并在之后清除幽灵阵列。如果用参数
       * @p   VectorOperation::insert,
       * 调用，则执行设置操作。因为在一个有鬼魂元素的向量中设置元素是不明确的（因为可以同时设置鬼魂位置和拥有位置上的元素），这个操作假设所有的数据都在拥有处理器上正确设置。因此在调用
       * compress(VectorOperation::insert),
       * 时，所有的ghost条目都被简单地清零（使用zero_ghost_values()）。在调试模式下，会进行一个检查，确保数据集在处理器之间实际上是一致的，也就是说，每当发现一个非零的鬼魂元素，就会与拥有处理器上的值进行比较，如果这些元素不一致就会抛出一个异常。
       *
       */
      virtual void
      compress(::dealii::VectorOperation::values operation) override;

      /**
       * 用存储在拥有处理器各自位置上的值来填充ghost索引的数据域。在从ghost中读取数据之前需要这个函数。该函数是
       * @p const
       * ，即使ghost数据被改变。需要这样做是为了让具有 @p
       * 常量向量的函数在不创建暂存器的情况下进行数据交换。
       *
       */
      void
      update_ghost_values() const;

      /**
       * 这个方法将ghost
       * dofs上的条目清零，但不触及本地拥有的DoF。
       * 调用这个方法后，禁止对向量中的鬼魂元素进行读取访问，并且会抛出一个异常。在这种状态下，只允许对鬼魂元素进行写入访问。
       * @deprecated  使用zero_out_ghost_values()代替。
       *
       */
      DEAL_II_DEPRECATED void
      zero_out_ghosts() const;


      /**
       * 这个方法将ghost
       * dofs上的条目清零，但不涉及本地拥有的DoF。
       * 调用该方法后，禁止对向量中的ghost元素进行读取访问，并且会抛出一个异常。在这种状态下，只允许对鬼魂元素进行写入访问。
       *
       */
      void
      zero_out_ghost_values() const;

      /**
       * 返回这个向量是否包含鬼魂元素。
       *
       */
      bool
      has_ghost_elements() const;

      /**
       * 这是一个集体添加操作，将存储在 @p values
       * 中的一整组值添加到 @p indices. 指定的向量成分中。
       *
       */
      template <typename OtherNumber>
      void
      add(const std::vector<size_type> &       indices,
          const ::dealii::Vector<OtherNumber> &values);

      /**
       * 缩放和简单的向量加法，即<tt>*this = s*(*this)+V</tt>。
       *
       */
      void
      sadd(const Number s, const BlockVector<Number> &V);

      /**
       * 返回向量是否只包含数值为0的元素。
       * 这个函数主要是用于内部一致性检查，在非调试模式下应该很少使用，因为它耗费了不少时间。
       *
       */
      virtual bool
      all_zero() const override;

      /**
       * 计算向量中所有条目的平均值。
       *
       */
      virtual Number
      mean_value() const override;

      /**
       * $l_p$
       * -向量的正负值。元素绝对值的p次方之和的p次根。
       *
       */
      real_type
      lp_norm(const real_type p) const;

      /**
       * 把这个向量的内容和另一个向量<tt>v</tt>交换。我们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换两个向量的数据指针，因此不需要分配临时存储和移动数据。
       * 限制：现在这个函数只在两个向量有相同数量的块的情况下工作。如果需要，也应该交换块的数量。
       * 这个函数类似于所有C++标准容器的swap()函数。此外，还有一个全局函数swap(u,v)，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数类比。
       *
       */
      void
      swap(BlockVector<Number> &v);
      //@}

      /**
       * @name  2：VectorSpaceVector的实现
       *
       */
      //@{

      /**
       * 将维度改为向量V的维度，V的元素不会被复制。
       *
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * 将整个向量乘以一个固定的因子。
       *
       */
      virtual BlockVector<Number> &
      operator*=(const Number factor) override;

      /**
       * 用整个向量除以一个固定的因子。
       *
       */
      virtual BlockVector<Number> &
      operator/=(const Number factor) override;

      /**
       * 将向量 @p V 添加到现在的向量中。
       *
       */
      virtual BlockVector<Number> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * 从现在的向量中减去向量 @p V 。
       *
       */
      virtual BlockVector<Number> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * 从输入向量 @p V.  VectorOperation::values  @p operation
       * 中导入向量的IndexSet中存在的所有元素，用于决定 @p
       * V
       * 中的元素是否应该添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
       *
       */
      virtual void
      import(const LinearAlgebra::ReadWriteVector<Number> &V,
             VectorOperation::values                       operation,
             std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
               communication_pattern = {}) override;

      /**
       * 返回两个向量的标量乘积。
       *
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * 计算这个向量的每个块和 @p V
       * 之间的标量积，并将结果存储在一个完整的矩阵 @p
       * matrix. 中。这个函数通过形成 $A_{ij}=U_i \cdot V_j$
       * 来计算结果，其中 $U_i$ 和 $V_j$ 分别表示 $i$
       * 的第1个块（不是元素！）以及 $j$ 的第1个块 $V$
       * ，。如果 @p symmetric 是 <code>true</code>
       * ，假定内积的结果是一个方形对称矩阵，几乎一半的标量积可以避免。
       * 显然，这个函数只有在两个向量的所有块都是相同大小的情况下才能使用。
       * @note
       * 在内部，将调用一个全局还原来累积本地拥有的自由度之间的标量积。
       *
       */
      template <typename FullMatrixType>
      void
      multivector_inner_product(FullMatrixType &           matrix,
                                const BlockVector<Number> &V,
                                const bool symmetric = false) const;

      /**
       * 使用度量张量 @p matrix. 计算这个向量的每个块和 @p V
       * 之间的标量积，这个函数计算 $ \sum_{ij} A^{ij} U_i \cdot
       * V_j$ 的结果，其中 $U_i$ 和 $V_j$ 分别表示 $U$
       * 的第1个块（不是元素）和 $j$ 的第2个块。如果 @p
       * symmetric 是 <code>true</code> ，则假定 $U_i \cdot V_j$ 和
       * $A^{ij}$ 是对称矩阵，几乎一半的标量积可以避免。
       * 显然，只有当两个向量的所有块都是相同大小时才能使用这个函数。
       * @note
       * 在内部，将调用一个全局还原来累积本地拥有的自由度之间的标量积。
       *
       */
      template <typename FullMatrixType>
      Number
      multivector_inner_product_with_metric(const FullMatrixType &     matrix,
                                            const BlockVector<Number> &V,
                                            const bool symmetric = false) const;

      /**
       * 设置该向量的每个块，如下所示。        $V^i = s V^i + b
       * \sum_{j} U_j A^{ji}$  其中 $V^i$ 和 $U_j$ 分别表示 $V$ 的第
       * $i$ 块（非元素）和 $j$ 的第 $U$ 块，。
       * 显然，这个函数只有在两个向量的所有块都是相同大小时才能使用。
       *
       */
      template <typename FullMatrixType>
      void
      mmult(BlockVector<Number> & V,
            const FullMatrixType &matrix,
            const Number          s = Number(0.),
            const Number          b = Number(1.)) const;

      /**
       * 将 @p a 添加到所有组件中。注意， @p a
       * 是一个标量，而不是一个向量。
       *
       */
      virtual void
      add(const Number a) override;

      /**
       * 向量的倍数的简单加法，即<tt>*this += a*V</tt>。
       *
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 缩放向量的多重加法，即：<tt>*this += a*V+b*W</tt>。
       *
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * 一个集体加法操作。这个函数将存储在 @p values
       * 中的一整套数值添加到 @p indices.
       * 指定的向量成分中。
       *
       */
      virtual void
      add(const std::vector<size_type> &indices,
          const std::vector<Number> &   values);

      /**
       * 缩放和简单的向量倍数相加，即<tt>*this =
       * s*(*this)+a*V</tt>。
       *
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
       *
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * 赋值 <tt>*this = a*V</tt>.
       *
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * 返回向量的l<sub>1</sub>规范（即所有处理器中所有条目的绝对值之和）。
       *
       */
      virtual real_type
      l1_norm() const override;

      /**
       * 返回向量的 $l_2$
       * 准则（即所有处理器中所有条目的平方之和的平方根）。
       *
       */
      virtual real_type
      l2_norm() const override;

      /**
       * 返回向量的 $l_2$ 规范的平方。
       *
       */
      real_type
      norm_sqr() const;

      /**
       * 返回向量的最大规范（即所有条目和所有处理器中的最大绝对值）。
       *
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * 执行一个向量加法和后续内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
       * @code
       * this->add(a, V);
       * return_value =this W;
       * @endcode
       * 这个函数存在的原因是这个操作比单独调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
       * @p this,   @p V,   @p W,
       * ，而调用单独的方法意味着要加载调用向量 @p this
       * 两次。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
       * @p W 等于 @p this).
       * ，则减少50\%）对于复值向量，第二步中的标量乘法被实现为
       * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
       *
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * 返回向量的全局大小，等于所有处理器中本地拥有的索引数之和。
       *
       */
      virtual size_type
      size() const override;

      /**
       * 返回一个索引集，描述这个向量的哪些元素为当前处理器所拥有。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
       * @code
       * vec.locally_owned_elements() == complete_index_set(vec.size())
       * @endcode
       *
       *
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * 将向量打印到输出流  @p out.  。
       *
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 3,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * 返回该类的内存消耗，单位是字节。
       *
       */
      virtual std::size_t
      memory_consumption() const override;
      //@}

      /**
       * @addtogroup  Exceptions  @{
       *
       */

      /**
       * 试图在两个不兼容的向量类型之间执行操作。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcVectorTypeNotCompatible);

      /**
       * 异常情况
       *
       */
      DeclException0(ExcIteratorRangeDoesNotMatchVectorSize);
      //@}
    };

     /*@}*/ 

  } // end of namespace distributed

} // end of namespace LinearAlgebra


/**
 * 全局函数，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  BlockVector
 *
 *
 */
template <typename Number>
inline void
swap(LinearAlgebra::distributed::BlockVector<Number> &u,
     LinearAlgebra::distributed::BlockVector<Number> &v)
{
  u.swap(v);
}


/**
 * 将 dealii::LinearAlgebra::distributed::BlockVector
 * 声明为分布式向量。
 *
 *
 */
template <typename Number>
struct is_serial_vector<LinearAlgebra::distributed::BlockVector<Number>>
  : std::false_type
{};


DEAL_II_NAMESPACE_CLOSE

#ifdef DEAL_II_MSVC
#  include <deal.II/lac/la_parallel_block_vector.templates.h>
#endif

#endif


