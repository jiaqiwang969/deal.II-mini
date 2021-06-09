//include/deal.II-translator/lac/vector_0.txt
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

#ifndef dealii_vector_h
#define dealii_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/differentiation/ad/ad_number_traits.h>

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <boost/serialization/split_member.hpp>

#include <algorithm>
#include <initializer_list>
#include <iosfwd>
#include <iterator>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declarations
#ifndef DOXYGEN
#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class VectorBase;
}
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace TrilinosWrappers
#  endif

template <typename number>
class LAPACKFullMatrix;

template <typename>
class BlockVector;

namespace parallel
{
  namespace internal
  {
    class TBBPartitioner;
  }
} // namespace parallel
#endif


/*!   @addtogroup Vectors  
     * @{ 

 
*
*/

/**
 * 一个表示数字元素向量的类。至于其他类，在 @ref
 * Vectors 组中，这个类有大量的成员函数。这些包括
 *
 *
 *
 * - 用于初始化向量或改变其大小的函数。
 *
 *
 *
 * - 计算向量属性的函数，如各种规范。
 *
 *
 * - 允许从向量的单个元素中读取或写入的函数。
 *
 *
 * - 实现向量代数运算的函数，例如向量的加法；以及
 *
 *
 *
 * - 允许输入和输出向量存储的数据的函数。
 * 与C++标准库类 `std::vector`,
 * 相比，该类打算实现的不是简单的允许访问其元素的数组，而实际上是一个矢量，是适合数值计算的
 * "矢量空间 "这一数学概念的成员。
 *
 *
 * @note
 * 这个模板的实例提供给<tt>  @<float@>,   @<double@>,
 * @<std::complex@<float@>@>,   @<std::complex@<double@>@></tt>;
 * 其他的可以在应用程序中生成（见手册中 @ref Instantiations
 * 部分）。
 *
 *
 */
template <typename Number>
class Vector : public Subscriptor
{
public:
  // The assertion in vector.templates.h for whether or not a number is
  // finite is not compatible for AD number types.
  static_assert(
    !Differentiation::AD::is_ad_number<Number>::value,
    "The Vector class does not support auto-differentiable numbers.");

  /**
   * 声明所有容器中使用的标准类型。这些类型与<tt>C++</tt>标准库中的<tt>vector<...></tt>类相似。
   *
   */
  using value_type      = Number;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;
  using size_type       = types::global_dof_index;

  /**
   * 声明一个类型，该类型持有与该类的模板参数相同精度的实值数。如果这个类的模板参数是一个实数数据类型，那么real_type就等于模板参数。
   * 如果模板参数是一个 std::complex
   * 类型，那么real_type等于复数的基础类型。
   * 这个别名被用来表示规范的返回类型。
   *
   */
  using real_type = typename numbers::NumberTraits<Number>::real_type;

  /**
   * @name  基本对象处理
   *
   */
  //@{
  /**
   * 构造函数。创建一个维数为0的向量。
   *
   */
  Vector();

  /**
   * 复制构造函数。将维数设置为给定的向量的维数，并复制所有元素。
   * 我们希望这个构造函数是显式的，但是标准的容器坚持隐式使用它。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector(const Vector<Number> &v);

  /**
   * 移动构造函数。通过窃取向量的内部数据创建一个新的向量
   * @p v.  。
   *
   */
  Vector(Vector<Number> &&v) noexcept = default;

  /**
   * 复制构造函数获取另一种数据类型的向量。
   * 如果没有从 @p OtherNumber 到 @p Number.
   * 的转换路径，这个构造函数将无法编译。
   * 当复制到一个数据元素精度较低的向量时，你可能会失去精度。
   *
   */
  template <typename OtherNumber>
  explicit Vector(const Vector<OtherNumber> &v);

  /**
   * 拷贝构造函数取一个 `std::initializer_list`.
   * 类型的对象，该构造函数可用于使用大括号封闭的数字列表来初始化向量，如下面的例子。
   * @code
   * Vector<double> v({1,2,3});
   * @endcode
   * 这将创建一个大小为3的向量，其（双精度）元素的值为1.0、2.0和3.0。
   * 如果没有从 @p OtherNumber 到 @p Number.
   * 的转换路径，这个构造函数将无法编译。
   * 当复制到一个数据元素精度较低的向量时，可能会失去精度。
   *
   */
  template <typename OtherNumber>
  explicit Vector(const std::initializer_list<OtherNumber> &v);

#ifdef DEAL_II_WITH_PETSC
  /**
   * 另一个复制构造函数：从一个PETSc向量类中复制数值。这个复制构造函数只有在配置时检测到PETSc时才可用。
   * 请注意，由于MPI中使用的通信模型，当 <code>v</code>
   * 是一个分布式向量时，只有当所有进程同时进行这一操作时才能成功。不可能只有一个进程获得一个并行向量的副本，而其他作业做其他事情。
   *
   */
  explicit Vector(const PETScWrappers::VectorBase &v);
#endif

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 另一个拷贝构造函数：从一个Trilinos包装向量中拷贝数值。
   * 这个拷贝构造函数只有在配置时检测到Trilinos时才可用。
   * @note
   * 由于MPI中使用的通信模型，这个操作只有在所有知道 @p
   * v 的进程（即 <code>v.get_mpi_communicator()</code>
   * 给出的进程）同时进行时才能成功。这意味着，除非你使用一个分裂的MPI通信器，那么通常不可能只有一个或一个子集的进程获得一个并行向量的副本，而其他工作则做其他事情。换句话说，调用这个函数是一个
   * "集体操作"，需要由所有共同分享 @p v.
   * 的MPI进程来执行。
   *
   */
  explicit Vector(const TrilinosWrappers::MPI::Vector &v);
#endif

  /**
   * 构造函数。设置维度为 @p n ，并将所有元素初始化为0。
   * 构造函数是明确的，以避免像这样的意外。
   * <tt>v=0;</tt>。据推测，用户希望将向量的每一个元素都设置为零，但相反，发生的是这样的调用。
   * <tt>v=向量 @<number@>(0);</tt>,
   * 即向量被一个长度为零的向量所取代。
   *
   */
  explicit Vector(const size_type n);

  /**
   * 用迭代器指向的给定范围的值初始化向量。这个函数是为了类似于
   * @p std::vector 类而存在的。
   *
   */
  template <typename InputIterator>
  Vector(const InputIterator first, const InputIterator last);

  /**
   * 销毁器，去分配内存。虚化，以使派生类的行为正确。
   *
   */
  virtual ~Vector() override = default;

  /**
   * 这个函数什么也不做，只是为了与并行向量类兼容而存在。
   * 对于并行向量封装类来说，这个函数压缩了向量的底层表示，即刷新了向量对象的缓冲区（如果它有的话）。这个函数在逐一写入向量元素之后，在对其进行其他操作之前是必要的。
   * 然而，对于这个类的实现，它是不重要的，因此是一个空函数。
   *
   */
  void
  compress(::dealii::VectorOperation::values operation =
             ::dealii::VectorOperation::unknown) const;

  /**
   * 将向量的维度改为 @p N.
   * 如果可能的话，这个向量的保留内存保持不变，以使事情更快；这可能会浪费一些内存，所以要记住这一点。
   * 然而，如果<tt>N==0</tt>所有的内存都被释放了，也就是说，如果你想调整向量的大小并释放不需要的内存，你必须先调用<tt>reinit(0)</tt>，然后再调用<tt>reinit(N)</tt>。这种引用的行为与标准库容器的行为类似。
   * 如果 @p omit_zeroing_entries
   * 是假的，那么向量就会被零所填充。
   * 否则，元素就会留下一个未指定的状态。
   * 这个函数是虚拟的，以便允许派生类单独处理内存。
   *
   */
  virtual void
  reinit(const size_type N, const bool omit_zeroing_entries = false);

  /**
   * 和上面一样，但在调整大小时将保留向量的值。
   * 如果我们的新尺寸更大，我们将有\f[ \mathbf V \rightarrow
   * \left( \begin{array}{c} \mathbf V   \\ \mathbf 0 \end{array} \right)
   * \f]，而如果想要的尺寸更小，则有\f[ \left( \begin{array}{c}
   * \mathbf V_1   \\ \mathbf V_2 \end{array} \right) \rightarrow \mathbf V_1
   * \f]。
   *
   */
  void
  grow_or_shrink(const size_type N);

  /**
   * 将<a href="https://en.wikipedia.org/wiki/Givens_rotation">Givens
   * rotation</a> @p csr （余弦、正弦和半径的三要素，见
   * Utilities::LinearAlgebra::givens_rotation()) ）应用于 @p i'th 和 @p
   * k'th 单位向量跨越的平面内的向量。
   *
   */
  void
  apply_givens_rotation(const std::array<Number, 3> &csr,
                        const size_type              i,
                        const size_type              k);

  /**
   * 将维度改为向量 @p V. 的维度，与其他 @p reinit
   * 函数的情况相同。     @p V
   * 的元素不会被复制，也就是说，这个函数与调用<tt>reinit
   * (V.size(), omit_zeroing_entries)</tt>相同。
   *
   */
  template <typename Number2>
  void
  reinit(const Vector<Number2> &V, const bool omit_zeroing_entries = false);

  /**
   * 交换这个向量和另一个向量的内容  @p v.
   * 人们可以用一个临时变量和复制数据元素来做这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。
   * 这个函数类似于所有C++标准容器的 @p swap
   * 函数。另外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
   * 这个函数是虚拟的，以便允许派生类单独处理内存。
   *
   */
  virtual void
  swap(Vector<Number> &v);

  /**
   * 将向量的所有分量设置为给定的数字  @p s.
   * 由于将标量分配给向量的语义不是很明确，这个操作符实际上只应该在你想将整个向量设置为0时使用。这允许使用直观的符号<tt>v=0</tt>。
   * 赋值其他值是不允许的，将来可能会被禁止。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator=(const Number s);

  /**
   * 复制给定的向量。如果有必要，可以调整当前向量的大小。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator=(const Vector<Number> &v);

  /**
   * 移动给定的向量。该操作符用向量 @p v
   * 的内部数据替换当前向量，并将 @p v
   * 重置为新的默认构建后的状态。
   *
   */
  Vector<Number> &
  operator=(Vector<Number> &&v) noexcept = default;

  /**
   * 复制给定的向量。如果有必要的话，调整当前向量的大小。
   * @dealiiOperationIsMultithreaded
   *
   */
  template <typename Number2>
  Vector<Number> &
  operator=(const Vector<Number2> &v);

  /**
   * 用于将块状向量分配给普通向量的复制操作。
   *
   */
  Vector<Number> &
  operator=(const BlockVector<Number> &v);

#ifdef DEAL_II_WITH_PETSC
  /**
   * 另一个复制操作符：从PETSc包装的向量类中复制数值。这个操作符只有在配置时检测到PETSc时才可用。
   * 请注意，由于MPI中使用的通信模型，当 <code>v</code>
   * 是一个分布式向量时，只有当所有进程同时进行该操作时，该操作才能成功。不可能只有一个进程获得一个并行向量的副本，而其他作业做其他事情。
   *
   */
  Vector<Number> &
  operator=(const PETScWrappers::VectorBase &v);
#endif


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 另一个拷贝操作符：从一个（顺序的或平行的，取决于底层编译器）Trilinos包装向量类中拷贝数值。这个操作符只有在配置时检测到Trilinos时才可用。
   * @note  由于MPI中使用的通信模型，该操作只有在所有了解
   * @p v 的进程（即 <code>v.get_mpi_communicator()</code>
   * 给出的进程）同时进行时才能成功。这意味着，除非你使用一个分裂的MPI通信器，那么通常不可能只有一个或一个子集的进程获得一个并行向量的副本，而其他工作做其他事情。换句话说，调用这个函数是一个
   * "集体操作"，需要由所有共同分享 @p v.
   * 的MPI进程来执行。
   *
   */
  Vector<Number> &
  operator=(const TrilinosWrappers::MPI::Vector &v);
#endif

  /**
   * 检验是否相等。这个函数假设现在的向量和要比较的向量已经有相同的大小，因为比较不同大小的向量反正没有什么意义。
   *
   */
  template <typename Number2>
  bool
  operator==(const Vector<Number2> &v) const;

  /**
   * 测试不平等。这个函数假定现在的向量和要比较的向量已经有相同的大小，因为比较不同大小的向量反正没有什么意义。
   *
   */
  template <typename Number2>
  bool
  operator!=(const Vector<Number2> &v) const;

  //@}


  /**
   * @name  标量积、规范和相关操作
   *
   */
  //@{

  /**
   * 返回两个向量的标量乘积。 返回类型是 @p this
   * 向量的基本类型，所以返回类型和计算结果的准确性取决于这个向量的参数顺序。
   * 对于复数向量，标量积的实现方式为
   * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
   * @dealiiOperationIsMultithreaded
   * 该算法使用成对求和，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  template <typename Number2>
  Number operator*(const Vector<Number2> &V) const;

  /**
   * 返回 $l_2$ -norm的平方。      @dealiiOperationIsMultithreaded
   * 该算法使用成对求和法，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  real_type
  norm_sqr() const;

  /**
   * 这个向量的元素的平均值。      @dealiiOperationIsMultithreaded
   * 该算法使用成对求和法，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  Number
  mean_value() const;

  /**
   * $l_1$  - 矢量的规范。绝对值的总和。
   * @dealiiOperationIsMultithreaded
   * 该算法使用成对求和，每次运行中的求和顺序相同，这使得每次运行的结果完全可以重复。
   *
   */
  real_type
  l1_norm() const;

  /**
   * $l_2$  - 矢量的规范。元素平方之和的平方根。
   * @dealiiOperationIsMultithreaded
   * 该算法使用成对求和，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  real_type
  l2_norm() const;

  /**
   * $l_p$  - 矢量的规范。元素绝对值的p次方之和的p次根。
   * @dealiiOperationIsMultithreaded
   * 该算法使用成对求和，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  real_type
  lp_norm(const real_type p) const;

  /**
   * 元素的最大绝对值。
   *
   */
  real_type
  linfty_norm() const;

  /**
   * 执行一个矢量加法和随后的内积的组合操作，返回内积的值。换句话说，这个函数的结果与用户调用的
   * @code
   * this->add(a, V);
   * return_value =this W;
   * @endcode
   * 这个函数存在的原因是这个操作比单独调用这两个函数涉及的内存转移要少。这个方法只需要加载三个向量，
   * @p this,   @p V,   @p W,
   * ，而调用单独的方法意味着要加载调用向量 @p this
   * 两次。由于大多数向量操作都有内存传输限制，这就使时间减少了25\%（如果
   * @p W 等于 @p this).
   * ，则减少50\%）对于复值向量，第二步的标量乘法被实现为
   * $\left<v,w\right>=\sum_i v_i \bar{w_i}$  。
   * @dealiiOperationIsMultithreaded
   * 该算法使用成对求和，每次运行的求和顺序相同，这使得每次运行的结果完全可重复。
   *
   */
  Number
  add_and_dot(const Number a, const Vector<Number> &V, const Vector<Number> &W);

  //@}


  /**
   * @name  数据访问
   *
   */
  //@{

  /**
   * 返回一个指向底层数据缓冲区的指针。
   *
   */
  pointer
  data();

  /**
   * 返回一个指向底层数据缓冲区的常量指针。
   *
   */
  const_pointer
  data() const;

  /**
   * 使 @p Vector
   * 类有点像C++标准库中的<tt>vector<>/tt>类，返回该向量元素的起点和终点的迭代器。
   *
   */
  iterator
  begin();

  /**
   * 返回到向量开始的常数迭代器。
   *
   */
  const_iterator
  begin() const;

  /**
   * 返回一个迭代器，指向超过数组末端的元素。
   *
   */
  iterator
  end();

  /**
   * 返回一个恒定的迭代器，指向超过数组末端的元素。
   *
   */
  const_iterator
  end() const;

  /**
   * 访问 @p ith 组件的值。
   *
   */
  Number
  operator()(const size_type i) const;

  /**
   * 访问 @p ith 组件作为一个可写的引用。
   *
   */
  Number &
  operator()(const size_type i);

  /**
   * 访问 @p ith 组件的值。    与operator()完全相同。
   *
   */
  Number operator[](const size_type i) const;

  /**
   * 访问 @p ith 组件作为可写引用。    与operator()完全相同。
   *
   */
  Number &operator[](const size_type i);

  /**
   * 与通过operator()获取向量的单个元素不同，这个函数允许一次性获取整个元素集。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
   * 如果当前的向量被称为 @p v,
   * ，那么这个函数就等同于代码
   * @code
   * for (unsigned int i = 0; i < indices.size(); ++i)
   *   values[i] = v[indices[i]];
   * @endcode
   * @pre  @p indices 和 @p values 数组的大小必须是一致的。
   *
   */
  template <typename OtherNumber>
  void
  extract_subvector_to(const std::vector<size_type> &indices,
                       std::vector<OtherNumber> &    values) const;

  /**
   * 这个函数不是通过operator()获得向量的单个元素，而是允许一次获得整个元素集。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
   * 如果当前的向量被称为 @p v,
   * ，那么这个函数就等同于代码
   * @code
   * ForwardIterator indices_p = indices_begin;
   * OutputIterator  values_p  = values_begin;
   * while (indices_p != indices_end)
   *   {
   *    values_p = v[*indices_p];
   *     ++indices_p;
   *     ++values_p;
   *   }
   * @endcode
   * @pre  必须能够写进从 @p values_begin
   * 开始的尽可能多的内存位置，因为有 @p indices_begin 和 @p
   * indices_end. 之间的迭代器。
   *
   */
  template <typename ForwardIterator, typename OutputIterator>
  void
  extract_subvector_to(ForwardIterator       indices_begin,
                       const ForwardIterator indices_end,
                       OutputIterator        values_begin) const;
  //@}


  /**
   * @name  修改向量
   *
   */
  //@{

  /**
   * 将给定的向量添加到当前的向量中。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator+=(const Vector<Number> &V);

  /**
   * 从现在的向量中减去给定的向量。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator-=(const Vector<Number> &V);

  /**
   * 一个集体加法操作。这个函数将存储在 @p values
   * 中的一整套数值添加到 @p indices. 指定的向量成分中。
   *
   */
  template <typename OtherNumber>
  void
  add(const std::vector<size_type> &  indices,
      const std::vector<OtherNumber> &values);

  /**
   * 这是第二次集体添加操作。作为区别，这个函数需要一个deal.II的数值向量。
   *
   */
  template <typename OtherNumber>
  void
  add(const std::vector<size_type> &indices, const Vector<OtherNumber> &values);

  /**
   * 取一个<tt>n_elements</tt>连续存储的地址，并将其添加到向量中。处理上述其他两个<tt>add()</tt>函数未涵盖的所有情况。
   *
   */
  template <typename OtherNumber>
  void
  add(const size_type    n_elements,
      const size_type *  indices,
      const OtherNumber *values);

  /**
   * 将 @p s 添加到所有组件中。注意， @p s
   * 是一个标量，而不是一个矢量。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  add(const Number s);

  /**
   * 缩放向量的多重加法，即<tt>*this += a*V+b*W</tt>。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  add(const Number          a,
      const Vector<Number> &V,
      const Number          b,
      const Vector<Number> &W);

  /**
   * 矢量的倍数的简单加法，即<tt>*this += a*V</tt>。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  add(const Number a, const Vector<Number> &V);

  /**
   * 缩放和简单的向量相加，即<tt>*this = s*(*this)+V</tt>。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  sadd(const Number s, const Vector<Number> &V);

  /**
   * 缩放和简单加法，即<tt>*this = s*(*this)+a*V</tt>。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  sadd(const Number s, const Number a, const Vector<Number> &V);

  /**
   * 将向量的每个元素按一个常数缩放。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator*=(const Number factor);

  /**
   * 用给定值的倒数来缩放向量的每个元素。
   * @dealiiOperationIsMultithreaded
   *
   */
  Vector<Number> &
  operator/=(const Number factor);

  /**
   * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟对角线缩放矩阵的乘法（和立即重新赋值）。
   * @dealiiOperationIsMultithreaded
   *
   */
  void
  scale(const Vector<Number> &scaling_factors);

  /**
   * 用参数中的相应元素来缩放这个向量的每个元素。这个函数主要是为了模拟用对角线缩放矩阵进行乘法（和立即重新分配）。
   *
   */
  template <typename Number2>
  void
  scale(const Vector<Number2> &scaling_factors);

  /**
   * 赋值 <tt>*this = a*u</tt>.       @dealiiOperationIsMultithreaded
   *
   */
  void
  equ(const Number a, const Vector<Number> &u);

  /**
   * 赋值 <tt>*this = a*u</tt>.
   *
   */
  template <typename Number2>
  void
  equ(const Number a, const Vector<Number2> &u);

  /**
   * 这个函数什么也不做，只是为了与 @p 并行向量类（如
   * LinearAlgebra::distributed::Vector 类）兼容而存在。
   *
   */
  void
  update_ghost_values() const;
  //@}


  /**
   * @name  输入和输出
   *
   */
  //@{
  /**
   * 打印到一个流。  @p precision 表示打印数值所需的精度，
   * @p scientific 是否应使用科学符号。如果 @p across 是 @p true
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
   * 将整个向量写到文件中。这是以二进制模式进行的，所以输出结果既不能被人类阅读，也不能（可能）被其他使用不同操作系统或数字格式的计算机阅读。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从一个文件中读取一个矢量en块。这是用上述函数的逆运算来完成的，所以它的速度相当快，因为位流没有被解释。
   * 如果有必要，矢量会被调整大小。
   * 一个原始形式的错误检查被执行，它将识别最直截了当的尝试，将一些数据解释为存储在文件中的位流向量，但不会超过。
   *
   */
  void
  block_read(std::istream &in);

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * @}
   *
   */

  /**
   * @name  对象的信息
   *
   */
  //@{

  /**
   * 如果给定的全局索引在这个处理器的本地范围内，返回真。
   * 由于这不是一个分布式矢量，该方法总是返回真。
   *
   */
  bool
  in_local_range(const size_type global_index) const;

  /**
   * 返回一个索引集，描述这个向量的哪些元素是由当前处理器拥有的。请注意，这个索引集不包括这个向量可能作为鬼魂元素存储在本地，但实际上是由另一个处理器拥有的元素。因此，如果这是一个分布式向量，在不同处理器上返回的索引集将形成不相交的集合，加起来就是完整的索引集。很明显，如果一个向量只在一个处理器上创建，那么结果将满足
   * @code
   * vec.locally_owned_elements() == complete_index_set (vec.size())
   * @endcode
   * 由于当前的数据类型不支持跨不同处理器的并行数据存储，所以返回的索引集是完整的索引集。
   *
   */
  IndexSet
  locally_owned_elements() const;

  /**
   * 返回向量的尺寸。
   *
   */
  size_type
  size() const;

  /**
   * 返回向量的局部维度。因为这个向量不支持分布式数据，所以这个值总是与size()相同。
   * @note  这个函数的存在是为了与 LinearAlgebra::ReadWriteVector.
   * 兼容。
   *
   */
  size_type
  locally_owned_size() const;

  /**
   * 返回向量是否只包含值为0的元素。这个函数主要用于内部一致性检查，在非调试模式下应该很少使用，因为它使用了相当多的时间。
   *
   */
  bool
  all_zero() const;

  /**
   * 如果向量没有负条目，即所有条目都是零或正，则返回
   * @p true
   * 。例如，这个函数被用来检查细化指标是否真的都是正的（或零）。
   * 这个函数显然只有在这个类的模板参数是一个实数类型时才有意义。如果它是一个复杂的类型，那么就会抛出一个异常。
   *
   */
  bool
  is_non_negative() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 这个函数的存在是为了与 @p 并行向量类（例如，
   * LinearAlgebra::distributed::Vector 类）兼容。
   * 总是返回false，因为这个实现是串行的。
   *
   */
  bool
  has_ghost_elements() const;
  //@}

private:
  /**
   * 这个向量所拥有的元素的数组。
   *
   */
  AlignedVector<Number> values;

  /**
   * 在初始化或重新初始化结束时使用的便利函数。根据循环分区器的当前状态和向量的长度，将其重置（如果需要）到正确的状态。
   *
   */
  void
  maybe_reset_thread_partitioner();

  /**
   * reinit函数的实际实现。
   *
   */
  void
  do_reinit(const size_type new_size,
            const bool      omit_zeroing_entries,
            const bool      reset_partitioner);

  /**
   * 对于带有TBB的并行循环，该成员变量存储循环的亲和力信息。
   *
   */
  mutable std::shared_ptr<parallel::internal::TBBPartitioner>
    thread_loop_partitioner;

  // Make all other vector types friends.
  template <typename Number2>
  friend class Vector;
};

 /*@}*/ 
 /*----------------------- Inline functions ----------------------------------*/ 


#ifndef DOXYGEN


//------------------------ declarations for explicit specializations
template <>
Vector<int>::real_type
Vector<int>::lp_norm(const real_type) const;


//------------------------ inline functions

template <typename Number>
inline Vector<Number>::Vector()
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  Vector<Number>::reinit(0);
}



template <typename Number>
template <typename OtherNumber>
Vector<Number>::Vector(const std::initializer_list<OtherNumber> &v)
  : Vector(v.begin(), v.end())
{}



template <typename Number>
template <typename InputIterator>
Vector<Number>::Vector(const InputIterator first, const InputIterator last)
{
  // allocate memory. do not initialize it, as we will copy over to it in a
  // second
  reinit(std::distance(first, last), true);
  std::copy(first, last, begin());
}



template <typename Number>
inline Vector<Number>::Vector(const size_type n)
{
  // virtual functions called in constructors and destructors never use the
  // override in a derived class
  // for clarity be explicit on which function is called
  Vector<Number>::reinit(n, false);
}



template <typename Number>
inline typename Vector<Number>::size_type
Vector<Number>::size() const
{
  return values.size();
}



template <typename Number>
inline typename Vector<Number>::size_type
Vector<Number>::locally_owned_size() const
{
  return values.size();
}



template <typename Number>
inline bool
Vector<Number>::in_local_range(const size_type) const
{
  return true;
}



template <typename Number>
inline typename Vector<Number>::pointer
Vector<Number>::data()
{
  return values.data();
}



template <typename Number>
inline typename Vector<Number>::const_pointer
Vector<Number>::data() const
{
  return values.data();
}



template <typename Number>
inline typename Vector<Number>::iterator
Vector<Number>::begin()
{
  return values.begin();
}



template <typename Number>
inline typename Vector<Number>::const_iterator
Vector<Number>::begin() const
{
  return values.begin();
}



template <typename Number>
inline typename Vector<Number>::iterator
Vector<Number>::end()
{
  return values.end();
}



template <typename Number>
inline typename Vector<Number>::const_iterator
Vector<Number>::end() const
{
  return values.end();
}



template <typename Number>
inline Number
Vector<Number>::operator()(const size_type i) const
{
  AssertIndexRange(i, size());
  return values[i];
}



template <typename Number>
inline Number &
Vector<Number>::operator()(const size_type i)
{
  AssertIndexRange(i, size());
  return values[i];
}



template <typename Number>
inline Number Vector<Number>::operator[](const size_type i) const
{
  return operator()(i);
}



template <typename Number>
inline Number &Vector<Number>::operator[](const size_type i)
{
  return operator()(i);
}



template <typename Number>
template <typename OtherNumber>
inline void
Vector<Number>::extract_subvector_to(const std::vector<size_type> &indices,
                                     std::vector<OtherNumber> &    values) const
{
  for (size_type i = 0; i < indices.size(); ++i)
    values[i] = operator()(indices[i]);
}



template <typename Number>
template <typename ForwardIterator, typename OutputIterator>
inline void
Vector<Number>::extract_subvector_to(ForwardIterator       indices_begin,
                                     const ForwardIterator indices_end,
                                     OutputIterator        values_begin) const
{
  while (indices_begin != indices_end)
    {
      *values_begin = operator()(*indices_begin);
      indices_begin++;
      values_begin++;
    }
}



template <typename Number>
inline Vector<Number> &
Vector<Number>::operator/=(const Number factor)
{
  AssertIsFinite(factor);
  Assert(factor != Number(0.), ExcZero());

  this->operator*=(Number(1.) / factor);
  return *this;
}



template <typename Number>
template <typename OtherNumber>
inline void
Vector<Number>::add(const std::vector<size_type> &  indices,
                    const std::vector<OtherNumber> &values)
{
  Assert(indices.size() == values.size(),
         ExcDimensionMismatch(indices.size(), values.size()));
  add(indices.size(), indices.data(), values.data());
}



template <typename Number>
template <typename OtherNumber>
inline void
Vector<Number>::add(const std::vector<size_type> &indices,
                    const Vector<OtherNumber> &   values)
{
  Assert(indices.size() == values.size(),
         ExcDimensionMismatch(indices.size(), values.size()));
  add(indices.size(), indices.data(), values.values.begin());
}



template <typename Number>
template <typename OtherNumber>
inline void
Vector<Number>::add(const size_type    n_indices,
                    const size_type *  indices,
                    const OtherNumber *values)
{
  for (size_type i = 0; i < n_indices; ++i)
    {
      AssertIndexRange(indices[i], size());
      Assert(
        numbers::is_finite(values[i]),
        ExcMessage(
          "The given value is not finite but either infinite or Not A Number (NaN)"));

      this->values[indices[i]] += values[i];
    }
}



template <typename Number>
template <typename Number2>
inline bool
Vector<Number>::operator!=(const Vector<Number2> &v) const
{
  return !(*this == v);
}



template <typename Number>
inline void Vector<Number>::compress(::dealii::VectorOperation::values) const
{}


template <typename Number>
inline bool
Vector<Number>::has_ghost_elements() const
{
  return false;
}

template <typename Number>
inline void
Vector<Number>::update_ghost_values() const
{}



// Moved from vector.templates.h as an inline function by Luca Heltai
// on 2009/04/12 to prevent strange compiling errors, after making
// swap virtual.
template <typename Number>
inline void
Vector<Number>::swap(Vector<Number> &v)
{
  values.swap(v.values);
  std::swap(thread_loop_partitioner, v.thread_loop_partitioner);
}



template <typename Number>
template <class Archive>
inline void
Vector<Number>::save(Archive &ar, const unsigned int) const
{
  // forward to serialization function in the base class.
  ar &static_cast<const Subscriptor &>(*this);
  ar &values;
}



template <typename Number>
template <class Archive>
inline void
Vector<Number>::load(Archive &ar, const unsigned int)
{
  // the load stuff again from the archive
  ar &static_cast<Subscriptor &>(*this);
  ar &values;
  maybe_reset_thread_partitioner();
}

#endif


/*!   @addtogroup Vectors  
     * @{ 

* 
*
*/


/**
 * 全局函数 @p swap
 * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  Vectors
 *
 *
 */
template <typename Number>
inline void
swap(Vector<Number> &u, Vector<Number> &v)
{
  u.swap(v);
}


/**
 * 输出操作符将一个向量写入流中。该操作符逐一输出向量的元素，条目之间有一个空格。每个条目都会根据输出流上设置的标志进行格式化。
 * @relatesalso  Vectors
 *
 *
 */
template <typename number>
inline std::ostream &
operator<<(std::ostream &out, const Vector<number> &v)
{
  Assert(v.size() != 0, ExcEmptyObject());
  AssertThrow(out, ExcIO());

  for (typename Vector<number>::size_type i = 0; i < v.size() - 1; ++i)
    out << v(i) << ' ';
  out << v(v.size() - 1);

  AssertThrow(out, ExcIO());

  return out;
}

 /*@}*/ 


/**
 * 声明  dealii::Vector<  数字 > 作为串行矢量。
 * @relatesalso  Vectors
 *
 *
 */
template <typename Number>
struct is_serial_vector<Vector<Number>> : std::true_type
{};


DEAL_II_NAMESPACE_CLOSE

#endif


