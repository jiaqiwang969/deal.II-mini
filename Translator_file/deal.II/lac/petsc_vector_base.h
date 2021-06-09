//include/deal.II-translator/lac/petsc_vector_base_0.txt
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

#ifndef dealii_petsc_vector_base_h
#  define dealii_petsc_vector_base_h


#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_PETSC

#    include <deal.II/base/index_set.h>
#    include <deal.II/base/subscriptor.h>

#    include <deal.II/lac/exceptions.h>
#    include <deal.II/lac/vector.h>
#    include <deal.II/lac/vector_operation.h>

#    include <petscvec.h>

#    include <utility>
#    include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#    ifndef DOXYGEN
template <typename number>
class Vector;

namespace PETScWrappers
{
  class VectorBase;
}
#    endif

/**
 * 一个命名空间，PETSc对象的封装类位于其中。
 *
 *
 * @ingroup PETScWrappers
 *
 * @ingroup Vectors
 *
 */
namespace PETScWrappers
{
  /**
   * @cond 内部
   *
   */

  /**
   * 一个命名空间，用于PETScWrapper成员的内部实现细节。
   * @ingroup PETScWrappers
   *
   */
  namespace internal
  {
    /**
     * 由于对PETSc向量的访问只通过函数进行，而不是通过获得一个向量元素的引用，所以我们需要一个包装器类，它的作用就像一个引用一样，基本上把所有的访问（读和写）重定向到这个类的成员函数。
     * 这个类实现了这样一个封装器：它用一个向量和其中的一个元素进行初始化，并有一个转换操作符来提取这个元素的标量值。它也有各种赋值运算符，用于向这一个元素写入。
     * @ingroup PETScWrappers
     *
     */
    class VectorReference
    {
    public:
      /**
       * 声明容器大小的类型。
       *
       */
      using size_type = types::global_dof_index;

    private:
      /**
       * 构造函数。它是私有的，以便只允许实际的向量类来创建它。
       *
       */
      VectorReference(const VectorBase &vector, const size_type index);

    public:
      /*复制构造函数。     
*
*/
      VectorReference(const VectorReference &vector) = default;

      /**
       * 这看起来像一个复制操作符，但做的事情与平常不同。特别是，它并不复制这个引用的成员变量。相反，它处理的情况是，我们有两个向量
       * @p v 和 @p w,
       * ，像<tt>v(i)=w(i)</tt>中那样分配元素。这里，赋值的左手和右手都有数据类型VectorReference，但我们真正的意思是赋值这两个引用所代表的向量元素。这个操作符实现了这个操作。还要注意的是，这使得我们可以使赋值运算符成为常数。
       *
       */
      const VectorReference &
      operator=(const VectorReference &r) const;

      /**
       * 与上面的函数相同，但用于非const引用对象。这个函数是需要的，因为编译器可能会自动为非常量对象生成一个拷贝操作符。
       *
       */
      VectorReference &
      operator=(const VectorReference &r);

      /**
       * 将向量的引用元素设置为<tt>s</tt>。
       *
       */
      const VectorReference &
      operator=(const PetscScalar &s) const;

      /**
       * 将<tt>s</tt>添加到矢量的引用元素中。
       *
       */
      const VectorReference &
      operator+=(const PetscScalar &s) const;

      /**
       * 从向量的参考元素中减去<tt>s</tt>。
       *
       */
      const VectorReference &
      operator-=(const PetscScalar &s) const;

      /**
       * 用<tt>s</tt>乘以矢量中的参考元素。
       *
       */
      const VectorReference &
      operator*=(const PetscScalar &s) const;

      /**
       * 将向量中的参考元素除以<tt>s</tt>。
       *
       */
      const VectorReference &
      operator/=(const PetscScalar &s) const;

      /**
       * 返回参考元素数值的实数部分。
       *
       */
      PetscReal
      real() const;

      /**
       * 返回参考元素值的虚数部分。
       * @note
       * 这个操作对实数来说没有定义，会产生一个异常。
       *
       */
      PetscReal
      imag() const;

      /**
       * 将引用转换为实际值，即返回向量中被引用元素的值。
       *
       */
      operator PetscScalar() const;
      /**
       * 异常情况
       *
       */
      DeclException3(
        ExcAccessToNonlocalElement,
        int,
        int,
        int,
        << "You tried to access element " << arg1
        << " of a distributed vector, but only elements in range [" << arg2
        << "," << arg3 << "] are stored locally and can be accessed."
        << "\n\n"
        << "A common source for this kind of problem is that you "
        << "are passing a 'fully distributed' vector into a function "
        << "that needs read access to vector elements that correspond "
        << "to degrees of freedom on ghost cells (or at least to "
        << "'locally active' degrees of freedom that are not also "
        << "'locally owned'). You need to pass a vector that has these "
        << "elements as ghost entries.");
      /**
       * 异常情况。
       *
       */
      DeclException2(ExcWrongMode,
                     int,
                     int,
                     << "You tried to do a "
                     << (arg1 == 1 ? "'set'" : (arg1 == 2 ? "'add'" : "???"))
                     << " operation but the vector is currently in "
                     << (arg2 == 1 ? "'set'" : (arg2 == 2 ? "'add'" : "???"))
                     << " mode. You first have to call 'compress()'.");

    private:
      /**
       * 指向我们所参考的向量。
       *
       */
      const VectorBase &vector;

      /**
       * 向量中被引用的元素的索引。
       *
       */
      const size_type index;

      // Make the vector class a friend, so that it can create objects of the
      // present type.
      friend class ::dealii::PETScWrappers::VectorBase;
    };
  } // namespace internal
  /**
   * @endcond
   *
   */


  /**
   * 所有在PETSc向量类型之上实现的向量类的基类。由于在PETSc中，所有的向量类型（即顺序和平行的）都是通过填充一个抽象对象的内容来建立的，而这个抽象对象只能通过一个独立于实际向量类型的指针来引用，所以我们可以在这个基类中实现几乎所有的向量功能。因此，这个类也可以作为任何类型的PETSc
   * <code>Vec</code>
   * 对象的deal.II兼容包装器使用。派生类将只需要提供创建一种或另一种矢量的功能。
   * 这个类的接口是以deal.II中现有的Vector类为模型的。它有几乎相同的成员函数，并且通常是可交换的。然而，由于PETSc只支持单一的标量类型（要么是double，要么是float，要么是一个复杂的数据类型），所以它没有模板化，只能与你的PETSc安装时定义的数据类型
   * @p PetscScalar 一起工作。
   * 请注意，只有在向量装配后调用了 @p VecAssemblyBegin 和 @p
   * VecAssemblyEnd
   * 这两个函数，PETSc才能保证操作符合你的期望。因此，你需要在实际使用矢量之前调用
   * Vector::compress() 。
   * @ingroup PETScWrappers
   *
   */
  class VectorBase : public Subscriptor
  {
  public:
    /**
     * 声明一些在所有容器中使用的标准类型。这些类型与<tt>C++</tt>标准库<tt>vector<...></tt>类中的类型平行。
     *
     */
    using value_type      = PetscScalar;
    using real_type       = PetscReal;
    using size_type       = types::global_dof_index;
    using reference       = internal::VectorReference;
    using const_reference = const internal::VectorReference;

    /**
     * 默认构造函数。它不做任何事情，派生类将不得不初始化数据。
     *
     */
    VectorBase();

    /**
     * 复制构造函数。将维度设置为给定的向量，并复制所有元素。
     *
     */
    VectorBase(const VectorBase &v);

    /**
     * 从一个PETSc
     * Vec对象初始化一个向量。注意，我们没有复制向量，也没有获得所有权，所以我们没有在析构函数中销毁PETSc对象。
     *
     */
    explicit VectorBase(const Vec &v);

    /**
     * 删除了复制赋值运算符，以避免意外的使用带来意外的行为。
     *
     */
    VectorBase &
    operator=(const VectorBase &) = delete;

    /**
     * 解构器。
     *
     */
    virtual ~VectorBase() override;

    /**
     * 释放所有内存并返回到与调用默认构造函数后相同的状态。
     *
     */
    virtual void
    clear();

    /**
     * 压缩PETSc对象的底层表示，即刷新矢量对象的缓冲区（如果它有的话）。这个函数在逐一写入矢量元素后，在对其进行任何其他操作之前是必要的。        更多信息请参见  @ref GlossCompress  "压缩分布式对象"
     * 。
     *
     */
    void
    compress(const VectorOperation::values operation);

    /**
     * 将向量的所有分量设置为给定的数字  @p s.
     * 只需将其传递给各个块对象，但我们仍然需要声明这个函数，以使讨论中给出的关于使构造函数显式的例子发挥作用。
     * 由于将标量分配给向量的语义并不立即明确，这个操作符实际上应该只在你想将整个向量设置为零时才使用。这样就可以使用直观的符号<tt>v=0</tt>。赋予其他的值是被弃用的，将来可能会被禁止使用。
     *
     */
    VectorBase &
    operator=(const PetscScalar s);

    /**
     * 检验是否相等。这个函数假定现在的向量和要比较的向量已经有相同的大小，因为比较不同大小的向量反正没有什么意义。
     *
     */
    bool
    operator==(const VectorBase &v) const;

    /**
     * 测试不平等。这个函数假定现在的向量和要比较的向量已经有相同的大小，因为比较不同大小的向量反正没有什么意义。
     *
     */
    bool
    operator!=(const VectorBase &v) const;

    /**
     * 返回向量的全局尺寸。
     *
     */
    size_type
    size() const;

    /**
     * 返回向量的局部尺寸，即存储在当前MPI进程中的元素数量。对于顺序向量，这个数字与size()相同，但对于并行向量，它可能更小。
     * 要想知道到底哪些元素是存储在本地的，可以使用local_range()或local_owned_elements()。
     * @deprecated 用local_owned_size()代替。
     *
     */
    DEAL_II_DEPRECATED
    size_type
    local_size() const;

    /**
     * 返回向量的本地维度，即存储在当前MPI进程中的元素数量。对于顺序向量，这个数字与size()相同，但对于并行向量，它可能更小。
     * 要想知道哪些元素确切地存储在本地，可以使用local_range()或local_owned_elements()。
     *
     */
    size_type
    locally_owned_size() const;

    /**
     * 返回一对指数，表明该向量的哪些元素被存储在本地。第一个数字是存储的第一个元素的索引，第二个数字是本地存储的最后一个元素之后的索引。如果这是一个顺序向量，那么结果将是一对（0,N），否则将是一对（i,i+n），其中<tt>n=locally_owned_size()</tt>。
     *
     */
    std::pair<size_type, size_type>
    local_range() const;

    /**
     * 返回 @p index 是否在本地范围内，另见local_range()。
     *
     */
    bool
    in_local_range(const size_type index) const;

    /**
     * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。请注意，这个索引集不包括这个向量可能在本地存储为幽灵元素，但实际上是由另一个处理器拥有的元素。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。
     * 很明显，如果一个向量只在一个处理器上创建，那么结果将满足
     * @code
     * vec.locally_owned_elements() == complete_index_set (vec.size())
     * @endcode
     *
     *
     */
    IndexSet
    locally_owned_elements() const;

    /**
     * 如果向量包含鬼魂元素，则返回。          @see   @ref GlossGhostedVector  "含有鬼魂元素的向量"
     *
     */
    bool
    has_ghost_elements() const;

    /**
     * 这个函数的存在只是为了与 @p
     * LinearAlgebra::distributed::Vector
     * 类兼容，并不做任何事情：这个类以不同的方式实现鬼魂值的更新，与底层的PETSc向量对象更加匹配。
     *
     */
    void
    update_ghost_values() const;

    /**
     * 提供对一个给定元素的访问，包括读和写。
     *
     */
    reference
    operator()(const size_type index);

    /**
     * 提供对一个元素的只读访问。
     *
     */
    PetscScalar
    operator()(const size_type index) const;

    /**
     * 提供对一个给定元素的访问，包括读和写。
     * 与operator()完全相同。
     *
     */
    reference operator[](const size_type index);

    /**
     * 提供对一个元素的只读访问。
     * 与operator()完全相同。
     *
     */
    PetscScalar operator[](const size_type index) const;

    /**
     * 一个集体的设置操作：这个函数允许一次性设置整个元素集，而不是设置一个向量的单个元素。
     * 要设置的元素的索引在第一个参数中说明，相应的值在第二个参数中说明。
     *
     */
    void
    set(const std::vector<size_type> &  indices,
        const std::vector<PetscScalar> &values);

    /**
     * 与通过operator()获取向量中的单个元素不同，这个函数允许一次性获取一整组元素。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
     * 如果当前的向量被称为 @p v,
     * ，那么这个函数就等同于代码
     * @code
     * for (unsigned int i=0; i<indices.size(); ++i)
     *   values[i] = v[indices[i]];
     * @endcode
     * @pre  @p indices 和 @p values 数组的大小必须是一致的。
     *
     */
    void
    extract_subvector_to(const std::vector<size_type> &indices,
                         std::vector<PetscScalar> &    values) const;

    /**
     * 这个函数不是通过operator()获得向量的单个元素，而是允许一次获得整个元素集。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
     * 如果当前的向量被称为 @p v,
     * ，那么这个函数就等同于代码
     * @code
     * ForwardIterator indices_p = indices_begin;
     * OutputIterator  values_p  = values_begin;
     * while (indices_p != indices_end)
     * {
     *  values_p = v[*indices_p];
     *   ++indices_p;
     *   ++values_p;
     * }
     * @endcode
     * @pre  必须能够写进从 @p values_begin
     * 开始的尽可能多的内存位置，因为有 @p indices_begin 和
     * @p indices_end. 之间的迭代器。
     *
     */
    template <typename ForwardIterator, typename OutputIterator>
    void
    extract_subvector_to(const ForwardIterator indices_begin,
                         const ForwardIterator indices_end,
                         OutputIterator        values_begin) const;

    /**
     * 一个集体的添加操作。这个函数将存储在 @p values
     * 中的一整组值添加到 @p indices. 指定的向量成分中。
     *
     */
    void
    add(const std::vector<size_type> &  indices,
        const std::vector<PetscScalar> &values);

    /**
     * 这是第二次集体添加操作。作为区别，这个函数需要一个deal.II的数值向量。
     *
     */
    void
    add(const std::vector<size_type> &       indices,
        const ::dealii::Vector<PetscScalar> &values);

    /**
     * 取一个<tt>n_elements</tt>连续存储的地址，并将其添加到向量中。处理上述其他两个<tt>add()</tt>函数未涵盖的所有情况。
     *
     */
    void
    add(const size_type    n_elements,
        const size_type *  indices,
        const PetscScalar *values);

    /**
     * 返回两个向量的标量乘积。这两个向量必须有相同的大小。
     * 对于复值向量，这将得到  $\left(v^\ast,vec\right)$  。
     *
     */
    PetscScalar operator*(const VectorBase &vec) const;

    /**
     * 返回 $l_2$ -norm的平方。
     *
     */
    real_type
    norm_sqr() const;

    /**
     * 返回这个向量的元素的平均值。
     *
     */
    PetscScalar
    mean_value() const;

    /**
     * 该向量的 $l_1$ -norm。绝对值的总和。
     * @note
     * 在3.7.0以前的复值PETSc中，这个规范被实现为复数向量元素的实部和虚部的绝对值之和。
     *
     */
    real_type
    l1_norm() const;

    /**
     * $l_2$  - 矢量的规范。 各元素的平方根之和。
     *
     */
    real_type
    l2_norm() const;

    /**
     * 元素绝对值的p次方之和的p次根。
     *
     */
    real_type
    lp_norm(const real_type p) const;

    /**
     * $l_\infty$
     * -向量的规范。返回具有最大绝对值的向量元素的值。
     *
     */
    real_type
    linfty_norm() const;

    /**
     * 执行一个矢量加法和后续内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用
     * @code
     * this->add(a, V);
     * return_value =this W;
     * @endcode
     * 这个函数存在的原因是为了与deal.II自己的向量类兼容，后者可以用较少的内存传输实现这个功能。然而，对于PETSc向量来说，这样的组合操作是不被原生支持的，因此其代价完全等同于单独调用这两个方法。
     * 对于复值向量，第二步中的标量乘积被实现为
     * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  .
     *
     */
    PetscScalar
    add_and_dot(const PetscScalar a, const VectorBase &V, const VectorBase &W);

    /**
     * 返回具有最大负值的向量元素的值。          @deprecated
     * 为了提高与其他继承自VectorSpaceVector的类的兼容性，这个函数已经被废弃。如果你需要使用这个功能，那么请使用PETSc函数VecMin代替。
     *
     */
    DEAL_II_DEPRECATED
    real_type
    min() const;

    /**
     * 返回具有最大正值的向量元素的值。          @deprecated
     * 这个函数已经被废弃，以提高与其他继承自VectorSpaceVector的类的兼容性。如果你需要使用这个功能，那么请使用PETSc函数VecMax代替。
     *
     */
    DEAL_II_DEPRECATED
    real_type
    max() const;

    /**
     * 返回向量是否只包含值为0的元素。这是一个集体操作。这个函数很昂贵，因为可能所有的元素都要被检查。
     *
     */
    bool
    all_zero() const;

    /**
     * 如果向量没有负的条目，即所有条目都是零或正，则返回
     * @p true
     * 。例如，这个函数用于检查细化指标是否真的都是正的（或零）。
     * @deprecated
     * 这个函数已经被废弃，以改善与其他继承自VectorSpaceVector的类的兼容性。
     *
     */
    DEAL_II_DEPRECATED
    bool
    is_non_negative() const;

    /**
     * 将整个向量乘以一个固定的因子。
     *
     */
    VectorBase &
    operator*=(const PetscScalar factor);

    /**
     * 将整个向量除以一个固定的因子。
     *
     */
    VectorBase &
    operator/=(const PetscScalar factor);

    /**
     * 将给定的向量添加到当前的向量中。
     *
     */
    VectorBase &
    operator+=(const VectorBase &V);

    /**
     * 从现在的向量中减去给定的向量。
     *
     */
    VectorBase &
    operator-=(const VectorBase &V);

    /**
     * 将 @p s 加到所有组件上。注意 @p s
     * 是一个标量而不是一个向量。
     *
     */
    void
    add(const PetscScalar s);

    /**
     * 一个向量的倍数的简单加法，即<tt>*this += a*V</tt>。
     *
     */
    void
    add(const PetscScalar a, const VectorBase &V);

    /**
     * 缩放向量的多重加法，即：<tt>*this += a*V+b*W</tt>。
     *
     */
    void
    add(const PetscScalar a,
        const VectorBase &V,
        const PetscScalar b,
        const VectorBase &W);

    /**
     * 缩放和简单的向量相加，即<tt>*this = s*(*this)+V</tt>。
     *
     */
    void
    sadd(const PetscScalar s, const VectorBase &V);

    /**
     * 缩放和简单加法，即：<tt>*this = s*(*this)+a*V</tt>。
     *
     */
    void
    sadd(const PetscScalar s, const PetscScalar a, const VectorBase &V);

    /**
     * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新分配）。
     *
     */
    void
    scale(const VectorBase &scaling_factors);

    /**
     * 赋值 <tt>*this = a*V</tt>.
     *
     */
    void
    equ(const PetscScalar a, const VectorBase &V);

    /**
     * 使用PETSc内部矢量查看器函数<tt>VecView</tt>打印PETSc矢量对象的值。默认格式是打印矢量的内容，包括矢量元素的索引。对于其他有效的视图格式，请参考http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecView.html
     *
     */
    void
    write_ascii(const PetscViewerFormat format = PETSC_VIEWER_DEFAULT);

    /**
     * 打印到一个流。  @p precision
     * 表示打印数值所需的精度， @p scientific
     * 是否应使用科学符号。如果 @p across 是 @p true
     * ，那么向量将被打印在一行中，而如果 @p false
     * 则元素被打印在单独的一行中。
     *
     */
    void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true,
          const bool         across     = true) const;

    /**
     * 交换这个向量和另一个向量的内容  @p v.
     * 人们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。
     * 这个函数类似于所有C++标准容器的 @p swap
     * 函数。此外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
     *
     */
    void
    swap(VectorBase &v);

    /**
     * 转换操作符，以获得对底层PETSc类型的访问。如果你这样做，你就切断了这个类可能需要的一些信息，所以这个转换操作符应该只在你知道你要做什么的情况下使用。特别是，它应该只用于对向量的只读操作。
     *
     */
    operator const Vec &() const;

    /**
     * 对内存消耗的估计（这个类没有实现）。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 返回一个对与此对象一起使用的MPI通信器对象的引用。
     *
     */
    virtual const MPI_Comm &
    get_mpi_communicator() const;

  protected:
    /**
     * 一个PETSc中的通用向量对象。实际的类型，一个连续的向量，在构造函数中被设置。
     *
     */
    Vec vector;

    /**
     * 表示这个向量是否有与之相关的鬼魂索引。这意味着并行程序中至少有一个进程有至少一个幽灵索引。
     *
     */
    bool ghosted;

    /**
     * 这个向量包含鬼魂值的全局索引。这个向量中的位置表示本地编号，在PETSc中使用。
     *
     */
    IndexSet ghost_indices;

    /**
     * 存储最后一个动作是写操作还是加操作。这个变量是
     * @p mutable
     * ，这样访问器类就可以写到它，即使它们引用的向量对象是常量。
     *
     */
    mutable VectorOperation::values last_action;

    // Make the reference class a friend.
    friend class internal::VectorReference;

    /**
     * 指定该向量是否是PETSc
     * Vec的所有者。如果它是由这个类创建的，这就是真的，并决定它是否在析构器中被销毁。
     *
     */
    bool obtained_ownership;

    /**
     * 集合设置或添加操作。这个函数由集体 @p set 和 @p add
     * 调用， @p add_values 标志设置为相应的值。
     *
     */
    void
    do_set_add_operation(const size_type    n_elements,
                         const size_type *  indices,
                         const PetscScalar *values,
                         const bool         add_values);
  };



  // ------------------- inline and template functions --------------

  /**
   * 全局函数 @p swap
   * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
   * @relatesalso   PETScWrappers::VectorBase .
   *
   */
  inline void
  swap(VectorBase &u, VectorBase &v)
  {
    u.swap(v);
  }

#    ifndef DOXYGEN
  namespace internal
  {
    inline VectorReference::VectorReference(const VectorBase &vector,
                                            const size_type   index)
      : vector(vector)
      , index(index)
    {}


    inline const VectorReference &
    VectorReference::operator=(const VectorReference &r) const
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<PetscScalar>(r);

      return *this;
    }



    inline VectorReference &
    VectorReference::operator=(const VectorReference &r)
    {
      // as explained in the class
      // documentation, this is not the copy
      // operator. so simply pass on to the
      // "correct" assignment operator
      *this = static_cast<PetscScalar>(r);

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      const PetscInt petsc_i = index;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      vector.last_action = VectorOperation::insert;

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator+=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::add) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::add, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::add;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the addend is
      // zero
      if (value == PetscScalar())
        return *this;

      // use the PETSc function to add something
      const PetscInt       petsc_i = index;
      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &value, ADD_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));


      return *this;
    }



    inline const VectorReference &
    VectorReference::operator-=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::add) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::add, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::add;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the addend is
      // zero
      if (value == PetscScalar())
        return *this;

      // use the PETSc function to
      // add something
      const PetscInt       petsc_i     = index;
      const PetscScalar    subtractand = -value;
      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &subtractand, ADD_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator*=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::insert;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the factor is
      // one
      if (value == 1.)
        return *this;

      const PetscInt    petsc_i   = index;
      const PetscScalar new_value = static_cast<PetscScalar>(*this) * value;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline const VectorReference &
    VectorReference::operator/=(const PetscScalar &value) const
    {
      Assert((vector.last_action == VectorOperation::insert) ||
               (vector.last_action == VectorOperation::unknown),
             ExcWrongMode(VectorOperation::insert, vector.last_action));

      Assert(!vector.has_ghost_elements(), ExcGhostsPresent());

      vector.last_action = VectorOperation::insert;

      // we have to do above actions in any
      // case to be consistent with the MPI
      // communication model (see the
      // comments in the documentation of
      // PETScWrappers::MPI::Vector), but we
      // can save some work if the factor is
      // one
      if (value == 1.)
        return *this;

      const PetscInt    petsc_i   = index;
      const PetscScalar new_value = static_cast<PetscScalar>(*this) / value;

      const PetscErrorCode ierr =
        VecSetValues(vector, 1, &petsc_i, &new_value, INSERT_VALUES);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      return *this;
    }



    inline PetscReal
    VectorReference::real() const
    {
#      ifndef PETSC_USE_COMPLEX
      return static_cast<PetscScalar>(*this);
#      else
      return PetscRealPart(static_cast<PetscScalar>(*this));
#      endif
    }



    inline PetscReal
    VectorReference::imag() const
    {
#      ifndef PETSC_USE_COMPLEX
      return PetscReal(0);
#      else
      return PetscImaginaryPart(static_cast<PetscScalar>(*this));
#      endif
    }

  } // namespace internal

  inline bool
  VectorBase::in_local_range(const size_type index) const
  {
    PetscInt             begin, end;
    const PetscErrorCode ierr =
      VecGetOwnershipRange(static_cast<const Vec &>(vector), &begin, &end);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    return ((index >= static_cast<size_type>(begin)) &&
            (index < static_cast<size_type>(end)));
  }


  inline IndexSet
  VectorBase::locally_owned_elements() const
  {
    IndexSet is(size());

    // PETSc only allows for contiguous local ranges, so this is simple
    const std::pair<size_type, size_type> x = local_range();
    is.add_range(x.first, x.second);
    return is;
  }



  inline bool
  VectorBase::has_ghost_elements() const
  {
    return ghosted;
  }



  inline void
  VectorBase::update_ghost_values() const
  {}



  inline internal::VectorReference
  VectorBase::operator()(const size_type index)
  {
    return internal::VectorReference(*this, index);
  }



  inline PetscScalar
  VectorBase::operator()(const size_type index) const
  {
    return static_cast<PetscScalar>(internal::VectorReference(*this, index));
  }



  inline internal::VectorReference VectorBase::operator[](const size_type index)
  {
    return operator()(index);
  }



  inline PetscScalar VectorBase::operator[](const size_type index) const
  {
    return operator()(index);
  }

  inline const MPI_Comm &
  VectorBase::get_mpi_communicator() const
  {
    static MPI_Comm comm;
    PetscObjectGetComm(reinterpret_cast<PetscObject>(vector), &comm);
    return comm;
  }

  inline void
  VectorBase::extract_subvector_to(const std::vector<size_type> &indices,
                                   std::vector<PetscScalar> &    values) const
  {
    Assert(indices.size() <= values.size(),
           ExcDimensionMismatch(indices.size(), values.size()));
    extract_subvector_to(indices.begin(), indices.end(), values.begin());
  }

  template <typename ForwardIterator, typename OutputIterator>
  inline void
  VectorBase::extract_subvector_to(const ForwardIterator indices_begin,
                                   const ForwardIterator indices_end,
                                   OutputIterator        values_begin) const
  {
    const PetscInt n_idx = static_cast<PetscInt>(indices_end - indices_begin);
    if (n_idx == 0)
      return;

    // if we are dealing
    // with a parallel vector
    if (ghosted)
      {
        // there is the possibility
        // that the vector has
        // ghost elements. in that
        // case, we first need to
        // figure out which
        // elements we own locally,
        // then get a pointer to
        // the elements that are
        // stored here (both the
        // ones we own as well as
        // the ghost elements). in
        // this array, the locally
        // owned elements come
        // first followed by the
        // ghost elements whose
        // position we can get from
        // an index set
        PetscInt       begin, end;
        PetscErrorCode ierr = VecGetOwnershipRange(vector, &begin, &end);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        Vec locally_stored_elements = nullptr;
        ierr = VecGhostGetLocalForm(vector, &locally_stored_elements);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscInt lsize;
        ierr = VecGetSize(locally_stored_elements, &lsize);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscScalar *ptr;
        ierr = VecGetArray(locally_stored_elements, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        for (PetscInt i = 0; i < n_idx; ++i)
          {
            const unsigned int index = *(indices_begin + i);
            if (index >= static_cast<unsigned int>(begin) &&
                index < static_cast<unsigned int>(end))
              {
                // local entry
                *(values_begin + i) = *(ptr + index - begin);
              }
            else
              {
                // ghost entry
                const unsigned int ghostidx =
                  ghost_indices.index_within_set(index);

                AssertIndexRange(ghostidx + end - begin, lsize);
                *(values_begin + i) = *(ptr + ghostidx + end - begin);
              }
          }

        ierr = VecRestoreArray(locally_stored_elements, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        ierr = VecGhostRestoreLocalForm(vector, &locally_stored_elements);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
    // if the vector is local or the
    // caller, then simply access the
    // element we are interested in
    else
      {
        PetscInt       begin, end;
        PetscErrorCode ierr = VecGetOwnershipRange(vector, &begin, &end);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        PetscScalar *ptr;
        ierr = VecGetArray(vector, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));

        for (PetscInt i = 0; i < n_idx; ++i)
          {
            const unsigned int index = *(indices_begin + i);

            Assert(index >= static_cast<unsigned int>(begin) &&
                     index < static_cast<unsigned int>(end),
                   ExcInternalError());

            *(values_begin + i) = *(ptr + index - begin);
          }

        ierr = VecRestoreArray(vector, &ptr);
        AssertThrow(ierr == 0, ExcPETScError(ierr));
      }
  }

#    endif // DOXYGEN
} // namespace PETScWrappers

DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_PETSC

#endif
 /*---------------------------- petsc_vector_base.h --------------------------*/ 


