//include/deal.II-translator/base/function_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_function_h
#define dealii_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_time.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <functional>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
template <int rank, int dim, typename Number>
class TensorFunction;
#endif

/**
 * 该类是一个一般函数的模型，给定一个评估函数的点，返回一个有一个或多个分量的值的向量。
 * 该类的目的是代表标量和向量值的函数。为此，我们把标量函数看作是向量值函数的一个特例，在前者的情况下，只有一个单分量的返回向量。由于处理向量的成本相对较高，这个类的接口有只要求向量值结果的单个分量的函数（这是在你知道你的函数是标量值的情况下通常需要的），也有可以要求整个结果向量的函数，其分量与函数对象代表的分量一样多。因此，对函数对象的访问是通过以下方法进行的。
 *
 * @code
 * // access to one component at one point
 * double
 * value(const Point<dim> & p,
 *     const unsigned int component = 0) const;
 *
 * // return all components at one point
 * void
 * vector_value(const Point<dim> &p,
 *            Vector<double>   &value) const;
 * @endcode
 *
 * 为了提高效率，还有一些函数一次返回一个或所有的分量在一个列表中的点。
 *
 * @code
 * // access to one component at several points
 * void
 * value_list(const std::vector<Point<dim>> &point_list,
 *          std::vector<double>           &value_list,
 *          const unsigned int             component = 0) const;
 *
 * // return all components at several points
 * void
 * vector_value_list(const std::vector<Point<dim>> &point_list,
 *                 std::vector<Vector<double>>   &value_list) const;
 * @endcode
 *
 * 此外，还有一些函数返回函数的梯度，甚至是在一个或几个点上的高阶导数。
 * 你通常只会重载那些你需要的函数；一次返回几个值的函数（value_list(),
 * vector_value_list(),
 * 和梯度类似物）会调用那些只返回一个值的函数（value(),
 * vector_value(),
 * 和梯度类似物），而那些被调用但没有重载的函数会产生一个异常。
 * 相反，在一个或几个点上返回所有分量的函数（即vector_value(),
 * vector_value_list()），将 <em> 而不是 </em>
 * 重复调用在一个点上返回一个分量的函数，对每个点和分量调用一次。原因是效率问题：这相当于调用了太多的虚拟函数。如果你有矢量值的函数，你应该同时为所有分量提供虚拟函数的重载。
 * 还要注意的是，除非只有非常少的调用次数，否则你应该重载所有的函数集（只返回一个值，以及那些返回整个数组的函数），因为对一个点值的评估成本往往比虚拟函数调用本身要低。
 * 对时间相关函数的支持可以在基类FunctionTime中找到。
 *
 *  <h3>Functions that return tensors</h3>
 * 如果你要处理的函数有一些先验已知的成分（例如，<tt>dim</tt>元素），你可以考虑使用TensorFunction类来代替。特别是，如果你返回的对象具有张量的属性，也就是说，它们是dim-dimensional向量或dim-by-dim矩阵，那么这一点是正确的。另一方面，像
 * VectorTools::interpolate 或 VectorTools::interpolate_boundary_values
 * 这样的函数肯定只想要当前类型的对象。你可以使用
 * VectorFunctionFromTensorFunction 类将前者转换为后者。
 *
 *  <h3>Functions that return vectors of other data types</h3>
 * 大多数时候，你的函数会有  $f : \Omega \rightarrow {\mathbb
 * R}^{n_\text{components}}$
 * 的形式。然而，在有些情况下，你希望函数在不同的数域上返回向量（或标量），例如，返回复数或复数向量的函数。
 * $f : \Omega \rightarrow {\mathbb C}^{n_\text{components}}$  .
 * 在这种情况下，你可以为这个类的第二个模板参数选择一个不同于默认的
 * @p double
 * 的值：它描述了用于你的返回值的每个组成部分的标量类型。它的默认值是
 * @p double, ，但在上面的例子中，它可以被设置为
 * <code>std::complex@<double@></code>  。  step-58 是一个例子。
 * @tparam  dim 函数的域 $\Omega$
 * 所在的范围空间的空间维度。因此，该函数将在类型为 @p
 * Point<dim>. 的对象上被评估。  @tparam  RangeNumberType
 * 作为该函数的范围（或图像）的向量空间的标量类型。如上所述，当前类型的对象代表从
 * ${\mathbb R}^\text{dim}$ 到 $S^{n_\text{components}}$ 的函数，其中
 * $S$ 是向量空间的基本标量类型。 $S$ 的类型是由 @p
 * RangeNumberType 的模板参数给出的。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class Function : public FunctionTime<
                   typename numbers::NumberTraits<RangeNumberType>::real_type>,
                 public Subscriptor
{
public:
  /**
   * 将模板参数的值作为一个静态成员常量导出。
   * 有时对一些表达式模板编程很有用。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 矢量组件的数量。
   *
   */
  const unsigned int n_components;

  /**
   * 用于表示时间的标量值实数类型。
   *
   */
  using time_type = typename FunctionTime<
    typename numbers::NumberTraits<RangeNumberType>::real_type>::time_type;

  /**
   * 构造函数。可以取一个成分数的初始值（默认为1，即一个标量函数），以及时间变量，默认为0。
   *
   */
  explicit Function(const unsigned int n_components = 1,
                    const time_type    initial_time = 0.0);

  /**
   * 复制构造函数。
   *
   */
  Function(const Function &f) = default;

  /**
   * 虚拟析构器；在这种情况下绝对必要。
   * 这个析构器被声明为纯虚拟的，这样这个类的对象就不能被创建。由于所有其他的虚拟函数都有一个伪实现，以避免派生类的开销，它们不能是抽象的。因此，我们可以生成这个类的对象，因为这个类的函数没有一个是抽象的。
   * 我们通过使这个类的析构器成为抽象虚函数来规避这个问题。这就保证了至少有一个成员函数是抽象的，因此也就不能创建Function类型的对象。
   * 然而，派生类没有必要明确地实现一个析构器：每个类都有一个析构器，要么明确地实现，要么由编译器隐含地生成，这就解决了任何派生类的抽象性问题，即使它们没有明确声明的析构器。
   * 尽管如此，由于派生类想要调用基类的析构器，这个析构器被实现了（尽管它是纯虚的）。
   *
   */
  virtual ~Function() override = 0;

  /**
   * 赋值运算符。这在这里只是为了让你可以在容器中拥有派生类的对象，或者以其他方式分配它们。如果你所赋值的对象与被赋值的对象有不同的组件数量，它将引发一个异常。
   *
   */
  Function &
  operator=(const Function &f);

  /**
   * 返回函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
   *
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * 返回一个向量值函数在某一点的所有分量。
   * <tt>values</tt>应事先有正确的大小，即#n_components。
   * 默认实现将为每个分量调用value()。
   *
   */
  virtual void
  vector_value(const Point<dim> &p, Vector<RangeNumberType> &values) const;

  /**
   * 将<tt>values</tt>设置为函数中指定分量在<tt>points</tt>的点值。
   * 假设<tt>values</tt>已经有合适的大小，即与<tt>points</tt>数组的大小相同。
   * 默认情况下，这个函数为每个点分别重复调用value()，以填充输出数组。
   *
   */
  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<RangeNumberType> & values,
             const unsigned int             component = 0) const;

  /**
   * 将<tt>values</tt>设为函数在<tt>points</tt>的点值。
   * 假设<tt>values</tt>已经有了合适的大小，即与<tt>points</tt>数组的大小相同，并且所有元素都是与本函数的分量相同的向量。
   * 默认情况下，这个函数对每个点单独重复调用vector_value()，以填充输出数组。
   *
   */
  virtual void
  vector_value_list(const std::vector<Point<dim>> &       points,
                    std::vector<Vector<RangeNumberType>> &values) const;

  /**
   * 对于函数的每个分量，填充一个值的向量，每个点一个。
   * Function中这个函数的默认实现是为每个分量调用value_list()。为了提高性能，可以在派生类中重新实现以加快性能。
   *
   */
  virtual void
  vector_values(const std::vector<Point<dim>> &            points,
                std::vector<std::vector<RangeNumberType>> &values) const;

  /**
   * 返回函数的指定分量在给定点的梯度。
   *
   */
  virtual Tensor<1, dim, RangeNumberType>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * 返回函数的所有分量在给定点的梯度。
   *
   */
  virtual void
  vector_gradient(
    const Point<dim> &                            p,
    std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const;

  /**
   * 将<tt>gradients</tt>设为函数的指定分量在<tt>点的梯度</tt>。
   * 假设<tt>gradients</tt>已经有正确的大小，即与<tt>points</tt>数组的大小相同。
   *
   */
  virtual void
  gradient_list(const std::vector<Point<dim>> &               points,
                std::vector<Tensor<1, dim, RangeNumberType>> &gradients,
                const unsigned int component = 0) const;

  /**
   * 对于函数的每个分量，填充一个梯度值的向量，每个点一个。
   * Function中这个函数的默认实现是为每个分量调用value_list()。为了提高性能，可以在派生类中重新实现以加快性能。
   *
   */
  virtual void
  vector_gradients(
    const std::vector<Point<dim>> &                            points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const;

  /**
   * 设置<tt>gradients</tt>为函数在<tt>points</tt>处的梯度，适用于所有组件。假设<tt>gradients</tt>已经有正确的大小，即与<tt>points</tt>数组的大小相同。
   * <tt>gradients</tt>的外循环是在列表中的点上，内循环是在函数的不同组成部分上。
   *
   */
  virtual void
  vector_gradient_list(
    const std::vector<Point<dim>> &                            points,
    std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients) const;

  /**
   * 计算给定分量在点<tt>p</tt>的拉普拉斯。
   *
   */
  virtual RangeNumberType
  laplacian(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * 计算在<tt>p</tt>点的所有分量的拉普拉斯，并将它们存储在<tt>值</tt>中。
   *
   */
  virtual void
  vector_laplacian(const Point<dim> &p, Vector<RangeNumberType> &values) const;

  /**
   * 计算在一组点上的一个分量的拉普拉斯。
   *
   */
  virtual void
  laplacian_list(const std::vector<Point<dim>> &points,
                 std::vector<RangeNumberType> & values,
                 const unsigned int             component = 0) const;

  /**
   * 计算一组点上的所有分量的拉普拉斯。
   *
   */
  virtual void
  vector_laplacian_list(const std::vector<Point<dim>> &       points,
                        std::vector<Vector<RangeNumberType>> &values) const;

  /**
   * 计算一个给定分量在点<tt>p</tt>的Hessian，也就是函数的梯度。
   *
   */
  virtual SymmetricTensor<2, dim, RangeNumberType>
  hessian(const Point<dim> &p, const unsigned int component = 0) const;

  /**
   * 计算点<tt>p</tt>处所有分量的Hessian，并将其存储在<tt>values</tt>中。
   *
   */
  virtual void
  vector_hessian(
    const Point<dim> &                                     p,
    std::vector<SymmetricTensor<2, dim, RangeNumberType>> &values) const;

  /**
   * 计算在一组点上的一个分量的Hessian。
   *
   */
  virtual void
  hessian_list(const std::vector<Point<dim>> &                        points,
               std::vector<SymmetricTensor<2, dim, RangeNumberType>> &values,
               const unsigned int component = 0) const;

  /**
   * 计算一组点上的所有分量的Hessians。
   *
   */
  virtual void
  vector_hessian_list(
    const std::vector<Point<dim>> &                                     points,
    std::vector<std::vector<SymmetricTensor<2, dim, RangeNumberType>>> &values)
    const;


  /**
   * 返回这个对象的内存消耗估计值，以字节为单位。
   * 这个函数是虚拟的，可以被派生类重载。
   *
   */
  virtual std::size_t
  memory_consumption() const;
};


namespace Functions
{
  /**
   * 提供一个函数，总是返回交给构造函数的常量值。
   * @ingroup functions
   *
   */
  template <int dim, typename RangeNumberType = double>
  class ConstantFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * 构造函数；将所有组件的值设置为所提供的一个。默认的组件数量是1。
     *
     */
    explicit ConstantFunction(const RangeNumberType value,
                              const unsigned int    n_components = 1);

    /**
     * 构造函数；接受一个 <tt>std::vector<RangeNumberType></tt>
     * 对象作为参数。组件的数量由<tt>values.size()</tt>决定。
     *
     */
    explicit ConstantFunction(const std::vector<RangeNumberType> &values);

    /**
     * 构造函数；接收一个<tt>Vector<RangeNumberType></tt>对象作为参数。组件的数量由<tt>values.size()</tt>决定。
     *
     */
    explicit ConstantFunction(const Vector<RangeNumberType> &values);

    /**
     * 构造函数；使用存储在[begin_ptr,
     * begin_ptr+n_components]中的任何东西来初始化一个新对象。
     *
     */
    ConstantFunction(const RangeNumberType *begin_ptr,
                     const unsigned int     n_components);

    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    vector_value(const Point<dim> &       p,
                 Vector<RangeNumberType> &return_value) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<RangeNumberType> & return_values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(
      const std::vector<Point<dim>> &       points,
      std::vector<Vector<RangeNumberType>> &return_values) const override;

    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    vector_gradient(
      const Point<dim> &                            p,
      std::vector<Tensor<1, dim, RangeNumberType>> &gradients) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &               points,
                  std::vector<Tensor<1, dim, RangeNumberType>> &gradients,
                  const unsigned int component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &                            points,
      std::vector<std::vector<Tensor<1, dim, RangeNumberType>>> &gradients)
      const override;

    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim> & point,
            const unsigned int component = 0) const override;

    virtual RangeNumberType
    laplacian(const Point<dim> & point,
              const unsigned int component = 0) const override;

    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * 存储常量函数值向量。
     *
     */
    std::vector<RangeNumberType> function_value_vector;
  };



  /**
   * 提供一个总是返回0的函数。很明显，这个函数的导数也是零。另外，如果该函数不是标量函数，它的所有分量都返回0，这可以通过向构造函数传递适当数量的分量来获得。
   * 当你想实现同质边界条件或零初始条件时，这个函数是有用的。
   * @ingroup functions
   *
   */
  template <int dim, typename RangeNumberType = double>
  class ZeroFunction : public ConstantFunction<dim, RangeNumberType>
  {
  public:
    /**
     * 构造函数。组件的数量被预设为1。
     *
     */
    explicit ZeroFunction(const unsigned int n_components = 1);
  };

  /**
   * 一个函数，其输出也是其输入。这个函数的一个可能的应用是插值或投影一个代表空间坐标的有限元场：例如，我们可以设置一个有限元场，用这个函数插值一个三角形的单元的位置（通过
   * VectorTools::interpolate()),
   * ，在拉格朗日参考框架中进行计算时很有用。
   * @ingroup functions
   *
   */
  template <int dim, typename RangeNumberType = double>
  class IdentityFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * 构造函数。组件的数量被设置为dim。
     *
     */
    IdentityFunction();

    /**
     * @copydoc   Function::value() .
     *
     */
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * @copydoc   Function::gradient() .
     *
     */
    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * @copydoc   Function::laplacian()
     *
     */
    virtual RangeNumberType
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * @copydoc   Function::hessian()   Function::hessian() 。
     *
     */
    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;
  };
} // namespace Functions

/**
 * 提供一个函数，它总是返回交给构造函数的常量值。
 * @deprecated  用 Functions::ConstantFunction 代替。
 *
 *
 */
template <int dim, typename RangeNumberType = double>
using ConstantFunction DEAL_II_DEPRECATED =
  Functions::ConstantFunction<dim, RangeNumberType>;

/**
 * 提供一个总是返回0的函数。
 * @deprecated  使用 Functions::ZeroFunction 代替。
 *
 *
 */
template <int dim, typename RangeNumberType = double>
using ZeroFunction DEAL_II_DEPRECATED =
  Functions::ZeroFunction<dim, RangeNumberType>;



/**
 * 这是一个常数向量值的函数，其中向量的一个或多个分量有一个常数值，所有其他分量都是零。 它作为 VectorTools::integrate_difference, 的权重函数特别有用，它允许只整合一个或几个向量分量，而不是整个向量值的解决方案。换句话说，它作为一个分量掩码，只选择一个分量（见 @ref GlossComponentMask "关于分量掩码的术语条目"
 * ）。参见 step-20 教程程序的详细解释和使用案例。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class ComponentSelectFunction
  : public Functions::ConstantFunction<dim, RangeNumberType>
{
public:
  /**
   * 如果只有一个分量为非零，则构造函数。参数表示选择的分量，该分量的值和向量分量的总数。
   *
   */
  ComponentSelectFunction(const unsigned int    selected,
                          const RangeNumberType value,
                          const unsigned int    n_components);

  /**
   * 构造函数。如前所述，但所选分量的值被假定为1。从本质上讲，这个函数就像一个掩码。
   *
   */
  ComponentSelectFunction(const unsigned int selected,
                          const unsigned int n_components);

  /**
   * 如果多个组件应具有非零的单位值（即这应该是一个多个组件的掩码），则构造函数。第一个参数表示一个半开的分量区间（例如
   * std::pair(0,dim)
   * 表示第一个dim分量），第二个参数是矢量分量的总数。
   *
   */
  ComponentSelectFunction(const std::pair<unsigned int, unsigned int> &selected,
                          const unsigned int n_components);


  /**
   * 用<tt>ConstantFunction  @<dim,  RangeNumberType  @></tt>
   * 对象的值代替函数值，并保持当前的选择模式。
   * 如果你想在不同的组件中有不同的值，这很有用，因为<tt>ComponentSelectFunction
   * @<dim,  RangeNumberType  @></tt>
   * 类提供的构造函数只能对所有组件有相同的值。
   * @note  我们从 @p f 开始复制基础组件值数据。所以 @p f
   * 的组件数量不能少于调用对象。
   *
   */
  virtual void
  substitute_function_value_with(
    const Functions::ConstantFunction<dim, RangeNumberType> &f);

  /**
   * 返回函数在给定点的所有组件的值。
   *
   */
  virtual void
  vector_value(const Point<dim> &       p,
               Vector<RangeNumberType> &return_value) const override;

  /**
   * 将<tt>values</tt>设为函数在<tt>点</tt>的点值，适用于所有组件。假设<tt>values</tt>已经有正确的大小，即与<tt>points</tt>数组的大小相同。
   *
   */
  virtual void
  vector_value_list(
    const std::vector<Point<dim>> &       points,
    std::vector<Vector<RangeNumberType>> &values) const override;

  /**
   * 返回这个对象的内存消耗估计值，单位是字节。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  /**
   * 所选组件的索引的半开区间。
   *
   */
  const std::pair<unsigned int, unsigned int> selected_components;
};



/**
 * 该类提供了一种转换标量函数的方法，即
 *
 * @code
 * RangeNumberType foo (const Point<dim> &);
 * @endcode
 * 转化为Function  @<dim@>.
 * 类型的对象，因为参数返回的是标量，所以结果显然是一个Function对象，对于它来说
 * <code>function.n_components == 1</code>
 * 。该类的工作方式是存储一个指向给定函数的指针，每次
 * <code>function.value(p,component)</code> 被调用时，都会调用
 * <code>foo(p)</code> 并返回相应的值。它还确保
 * <code>component</code> 实际上是零，这也是标量函数的需要。
 * 该类提供了一种简单的方法，可以将一个简单的全局函数变成具有所需的Function
 * @<dim@> 接口的东西，用于 VectorTools::interpolate_boundary_values()
 * 等操作，从而允许更简单的实验，而不必写所有的锅炉板代码，即声明一个从Function派生的类并实现
 * Function::value() 函数。在  step-53
 * 的结果部分给出了一个例子。
 * 这个类获得了额外的表达能力，因为它所接受的参数不需要是一个指向实际函数的指针。相反，它是一个函数对象，也就是说，它也可以是一个lambda函数的结果或者其他一些可以用一个参数调用的对象。例如，如果你需要一个函数对象来返回一个点的法线，你可以这样写。
 *
 * @code
 * template <int dim, typename RangeNumberType>
 * class Norm : public Function<dim, RangeNumberType>
 * {
 * public:
 * virtual RangeNumberType
 * value(const Point<dim> & p,
 *       const unsigned int component) const
 * {
 *   Assert (component == 0, ExcMessage ("This object is scalar!"));
 *   return p.norm();
 * }
 * };
 *
 * Norm<2> my_norm_object;
 * @endcode
 * 然后把 <code>my_norm_object</code>
 * 传过来，或者你可以这样写。
 *
 * @code
 * ScalarFunctionFromFunctionObject<dim, RangeNumberType> my_norm_object(
 * &Point<dim>::norm);
 * @endcode
 *
 * 同样地，为了生成一个计算到点 <code>q</code>
 * 的距离的对象，我们可以这样做。
 * @code
 * template <int dim, typename RangeNumberType>
 * class DistanceTo : public Function<dim, RangeNumberType>
 * {
 * public:
 * DistanceTo (const Point<dim> &q) : q(q) {}
 *
 * virtual RangeNumberType
 * value (const Point<dim> & p,
 *        const unsigned int component) const
 * {
 *   Assert(component == 0, ExcMessage("This object is scalar!"));
 *   return q.distance(p);
 * }
 * private:
 * const Point<dim> q;
 * };
 *
 * Point<2> q (2, 3);
 * DistanceTo<2> my_distance_object;
 * @endcode
 * 或者我们可以这样写。
 *
 * @code
 * ScalarFunctionFromFunctionObject<dim, RangeNumberType> my_distance_object(
 * [&q](const Point<dim> &p){return q.distance(p);});
 * @endcode
 * 这样写，节省的工作量是显而易见的。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class ScalarFunctionFromFunctionObject : public Function<dim, RangeNumberType>
{
public:
  /**
   * 给出一个接受Point并返回RangeNumberType值的函数对象，将其转换为一个符合Function<dim,
   * RangeNumberType>接口的对象。
   *
   */
  explicit ScalarFunctionFromFunctionObject(
    const std::function<RangeNumberType(const Point<dim> &)> &function_object);

  /**
   * 返回函数在给定点的值。返回给与构造函数的函数在这个点上产生的值。
   *
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

private:
  /**
   * 当这个类的value()或value_list()函数被调用时，我们调用的函数对象。
   *
   */
  const std::function<RangeNumberType(const Point<dim> &)> function_object;
};



/**
 * 这个类与ScalarFunctionFromFunctionObject类类似，它允许轻松地将一个函数对象转换为满足Function基类接口的东西。不同的是，这里给定的函数对象仍然是一个标量函数（即它在每个空间点都有一个单一的值），但生成的Function对象是矢量值的。矢量分量的数量是在构造函数中指定的，在构造函数中还选择了这些矢量分量中的一个，应该由传递的对象来填充。结果是一个向量Function对象，每个分量的返回值都是0，除了所选的那个分量，它返回的是构造函数第一个参数所给的值。
 *
 *
 * @note
 * 在上面的讨论中，请注意（标量）"函数对象"（即可以像
 * <code>x(p)</code> 那样调用的C++对象 <code>x</code>
 * ）和大写的（矢量值）"函数对象"（即从Function基类派生的类的对象）之间的区别。
 * 为了更加具体，让我们考虑下面的例子。
 *
 * @code
 * RangeNumberType
 * one(const Point<2> &p)
 * {
 * return 1.0;
 * }
 *
 * VectorFunctionFromScalarFunctionObject<2, RangeNumberType>
 * component_mask(&one, 1, 3);
 * @endcode
 * 这里， <code>component_mask</code>
 * 就代表了一个Function对象，对于每一个点都会返回向量
 * $(0, 1, 0)^T$ ，即一个掩码函数，例如，可以传递给
 * VectorTools::integrate_difference().
 * 这个效果也可以用ComponentSelectFunction类来实现，但显然很容易扩展到那些在其一个分量中非恒定的函数。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class VectorFunctionFromScalarFunctionObject
  : public Function<dim, RangeNumberType>
{
public:
  /**
   * 给出一个接受Point并返回RangeNumberType值的函数对象，将其转换为符合Function
   * @<dim@> 接口的对象。      @param  function_object
   * 标量函数，它将构成结果Function对象的一个组成部分。
   * @param  n_components 产生的Function对象的向量组件总数。
   * @param  selected_component
   * 应该由第一个参数填充的单个组件。
   *
   */
  VectorFunctionFromScalarFunctionObject(
    const std::function<RangeNumberType(const Point<dim> &)> &function_object,
    const unsigned int selected_component,
    const unsigned int n_components);

  /**
   * 返回函数在给定点的值。返回给与构造函数的函数在此点上产生的值。
   *
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * 返回一个矢量值函数在给定点的所有分量。
   * <tt>值</tt>应事先有正确的大小，即#n_components。
   *
   */
  virtual void
  vector_value(const Point<dim> &       p,
               Vector<RangeNumberType> &values) const override;

private:
  /**
   * 当这个类的value()或value_list()函数被调用时，我们调用的函数对象。
   *
   */
  const std::function<RangeNumberType(const Point<dim> &)> function_object;

  /**
   * 矢量分量，其值要由给定的标量函数来填充。
   *
   */
  const unsigned int selected_component;
};


/**
 * 这个类与ScalarFunctionFromFunctionObject和VectorFunctionFromFunctionObject类类似，它允许轻松地将函数对象的向量转换为满足Function基类接口的东西。
 * 不同的是，这里生成的Function对象可以是矢量值的，而且你可以指定函数的梯度。矢量分量的数量是由构造函数中的矢量大小推导出来的。
 * 为了更加具体，让我们考虑下面的例子。
 *
 *
 * @code
 * RangeNumberType
 * first_component(const Point<2> &p)
 * {
 * return 1.0;
 * }
 *
 * RangeNumberType
 * second_component(const Point<2> &p)
 * {
 * return 2.0;
 * }
 *
 * Tensor<1, 2, RangeNumberType>
 * zero_gradient(const Point<2> &) {
 * return Tensor<1, 2, RangeNumberType>();
 * }
 *
 * FunctionFromFunctionObjects<2, RangeNumberType>
 *   custom_function({&first_component, &second_component},
 *                   {&zero_gradient, &zero_gradient});
 * @endcode
 *
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class FunctionFromFunctionObjects : public Function<dim, RangeNumberType>
{
public:
  /**
   * 默认构造函数。
   * 这个构造函数不初始化内部方法。为了得到一个可用的函数，你至少需要调用set_function_values()方法。如果你还需要解的梯度，那么你还必须调用set_function_gradients()方法。
   *
   */
  explicit FunctionFromFunctionObjects(const unsigned int n_components = 1,
                                       const double       initial_time = 0);

  /**
   * 函数的构造函数，你只知道它的值。
   * 由此产生的函数将有一个等于向量大小的分量  @p values.
   * 对 FunctionFromFunctionObject::gradient()
   * 方法的调用将触发一个异常，除非你先调用set_function_gradients()方法。
   *
   */
  explicit FunctionFromFunctionObjects(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
      &          values,
    const double initial_time = 0.0);

  /**
   * 函数的构造函数，你知道它的值和梯度。
   * 由此产生的函数将具有与向量 @p values.
   * 大小相等的分量，如果 @p values 和 @p gradients
   * 的大小不匹配，将触发一个异常。
   *
   */
  FunctionFromFunctionObjects(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
      &values,
    const std::vector<
      std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>>
      &          gradients,
    const double initial_time = 0.0);


  /**
   * 返回函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
   *
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * 返回该函数在给定点的梯度。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
   *
   */
  virtual Tensor<1, dim, RangeNumberType>
  gradient(const Point<dim> & p,
           const unsigned int component = 0) const override;

  /**
   * 重置这个对象的函数值。如果 @p values
   * 参数的大小与此对象的组件数不匹配，就会抛出一个断言。
   *
   */
  void
  set_function_values(
    const std::vector<std::function<RangeNumberType(const Point<dim> &)>>
      &values);

  /**
   * 重置此对象的函数梯度。如果 @p gradients
   * 参数的大小与此对象的组件数量不匹配，则抛出一个断言。
   *
   */
  void
  set_function_gradients(
    const std::vector<
      std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>>
      &gradients);

private:
  /**
   * 实际的函数值。
   *
   */
  std::vector<std::function<RangeNumberType(const Point<dim> &)>>
    function_values;

  /**
   * 实际的函数梯度。
   *
   */
  std::vector<
    std::function<Tensor<1, dim, RangeNumberType>(const Point<dim> &)>>
    function_gradients;
};


/**
 * 这个类是作为翻译<code>Tensor<1,dim,
 * RangeNumberType></code>类型的对象所产生的值的一种手段，并将它们作为同一事物的多组件版本以矢量形式返回，以便用于例如
 * VectorTools::interpolate
 * 或许多其他采取Function对象的函数。它允许用户将所需的组件放入一个<tt>n_components</tt>长向量中，从该向量中的<tt>selected_component</tt>位置开始，并让所有其他组件为0。
 * 比如说。假设你创建了一个叫做
 *
 * @code
 * class RightHandSide : public TensorFunction<rank,dim, RangeNumberType>
 * @endcode
 * 它扩展了TensorFunction类，你有一个对象
 *
 * @code
 * RightHandSide<1,dim, RangeNumberType> rhs;
 * @endcode
 * 该类的对象，你想用 VectorTools::interpolate
 * 函数插值到你的网格上，但你用于DoFHandler对象的有限元有3份<tt>dim</tt>分量的有限元，总共有3*dim分量。为了插值到该DoFHandler上，你需要一个具有3*dim向量分量的Function类型的对象。使用这段代码可以从现有的
 * <code>rhs</code> 对象创建这样一个对象。
 *
 * @code
 * RightHandSide<1,dim, RangeNumberType> rhs;
 * VectorFunctionFromTensorFunction<dim, RangeNumberType> rhs_vector_function(
 * rhs, 0, 3*dim);
 * @endcode
 * 其中张量函数的 <code>dim</code>
 * 分量被放置到函数对象的第一个 <code>dim</code> 分量中。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim, typename RangeNumberType = double>
class VectorFunctionFromTensorFunction : public Function<dim, RangeNumberType>
{
public:
  /**
   * 给定一个TensorFunction对象，它接受一个<tt>Point</tt>并返回一个<tt>Tensor<1,dim,
   * RangeNumberType></tt>值，把它转换成一个符合Function  @<dim@>
   * 接口的对象。
   * 默认情况下，创建一个与<tt>tensor_function</tt>返回的大小相同的Vector对象，即有<tt>dim</tt>组件。
   * @param  tensor_function
   * 将形成生成的矢量函数对象的一个分量的TensorFunction。
   * @param  n_components
   * 产生的TensorFunction对象的矢量分量的总数。    @param
   * selected_component 应该由第一个参数填充的第一个分量。
   * 这应该使整个张量函数适合在<tt>n_component</tt>长度的返回向量内。
   *
   */
  explicit VectorFunctionFromTensorFunction(
    const TensorFunction<1, dim, RangeNumberType> &tensor_function,
    const unsigned int                             selected_component = 0,
    const unsigned int                             n_components       = dim);

  /**
   * 这个析构器被定义为虚拟的，以便与类的所有其他方面相吻合。
   *
   */
  virtual ~VectorFunctionFromTensorFunction() override = default;

  /**
   * 在给定的点上返回一个矢量值函数的单个分量。
   *
   */
  virtual RangeNumberType
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * 返回一个向量值函数在某一点的所有分量。
   * <tt>values</tt>应事先有正确的大小，即#n_components。
   *
   */
  virtual void
  vector_value(const Point<dim> &       p,
               Vector<RangeNumberType> &values) const override;

  /**
   * 返回一个矢量值函数在一个点的所有分量。
   * <tt>value_list</tt>的大小应与<tt>points</tt>相同，向量的每个元素将被传递给vector_value()来评估函数。
   *
   */
  virtual void
  vector_value_list(
    const std::vector<Point<dim>> &       points,
    std::vector<Vector<RangeNumberType>> &value_list) const override;

private:
  /**
   * TensorFunction对象，当这个类的vector_value()或vector_value_list()函数被调用时，我们会调用它。
   *
   */
  const TensorFunction<1, dim, RangeNumberType> &tensor_function;

  /**
   * 第一个向量分量，其值将由给定的TensorFunction填充。
   * 对于一个<tt>TensorFunction<1,dim,
   * RangeNumberType></tt>对象，这些值将被放置在selected_component到selected_component+dim-1的分量中。
   *
   */
  const unsigned int selected_component;
};


#ifndef DOXYGEN
// icc 2018 complains about an undefined reference
// if we put this in the templates.h file
//
// The destructor is pure virtual so we can't default it
// in the declaration.
template <int dim, typename RangeNumberType>
inline Function<dim, RangeNumberType>::~Function() = default;
#endif


DEAL_II_NAMESPACE_CLOSE

#endif


