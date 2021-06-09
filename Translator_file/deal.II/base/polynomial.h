//include/deal.II-translator/base/polynomial_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_polynomial_h
#define dealii_polynomial_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials   
     * @{ 
 *
 */

/**
 * 一个命名空间，与描述1d多项式空间有关的类在其中被声明。
 *
 *
 */
namespace Polynomials
{
  /**
   * 所有一维多项式的基类。一个多项式在这个类中由它的系数来表示，系数是通过构造函数或派生类来设置的。
   * 多项式的评估有两条路径。一个是基于系数，通过Horner方案进行评估，这是一个强大的通用方案。对于根在单位区间内的高阶多项式，另一种更稳定的评估方法是通过根的乘积来提供。这种形式可用于特殊多项式，如Lagrange多项式或Legendre多项式，并与相应的构造器一起使用。为了获得这种更稳定的评估形式，必须使用拉格朗日多项式形式的根的构造器。如果做了一个改变根的操作，则表示法将被切换为系数形式。
   * 这个类是TensorProductPolynomials类可能的模板参数的一个典型例子。
   *
   */
  template <typename number>
  class Polynomial : public Subscriptor
  {
  public:
    /**
     * 构造函数。多项式的系数作为参数传递，表示多项式
     * $\sum_i a[i] x^i$
     * ，即数组的第一个元素表示常数项，第二个表示线性项，以此类推。因此这个对象所代表的多项式的度数是<tt>系数</tt>数组中的元素数减去1。
     *
     */
    Polynomial(const std::vector<number> &coefficients);

    /**
     * 构造函数创建一个零度的多项式  @p n.  。
     *
     */
    Polynomial(const unsigned int n);

    /**
     * 拉格朗日多项式的构造函数和它的评估点。其思路是构造
     * $\prod_{i\neq j} \frac{x-x_i}{x_j-x_i}$
     * ，其中j是作为参数指定的评估点，支持点包含所有的点（包括x_j，内部将不被存储）。
     *
     */
    Polynomial(const std::vector<Point<1>> &lagrange_support_points,
               const unsigned int           evaluation_point);

    /**
     * 默认构造函数创建一个非法对象。
     *
     */
    Polynomial();

    /**
     * 返回这个多项式在给定点的值。
     * 这个函数对所提供的多项式的形式使用最稳定的数值评估算法。如果多项式是根的乘积形式，那么评估是基于(x)形式的乘积。
     *
     * - x_i), 而对于系数形式的多项式，则使用霍纳方案。
     *
     */
    number
    value(const number x) const;

    /**
     * 返回多项式在<tt>x</tt>点的值和导数。 <tt>values[i],
     * i=0,...,values.size()-1</tt>包括<tt>i</tt>的导数。因此，要计算的导数的数量由传递的数组的大小决定。
     * 这个函数使用Horner方案对系数形式的多项式或涉及根的项的乘积（如果使用该表示法）进行数值稳定评估。
     *
     */
    void
    value(const number x, std::vector<number> &values) const;

    /**
     * 返回多项式在<tt>x</tt>点的值和导数。 <tt>values[i],
     * i=0,...,n_derivatives</tt>包括<tt>i</tt>的导数。要计算的导数数量由
     * @p n_derivatives 决定， @p values 必须为 @p n_derivatives
     * +1个值提供足够的空间。
     * 这个函数对于所提供的多项式的形式使用最稳定的数值评估算法。如果多项式是根的乘积形式，评估是基于形式为(x
     *
     * - x_i），而Horner方案则用于系数形式的多项式。        模板类型`Number2`必须实现算术运算，如与多项式的`number`类型的加法或乘法，并且必须通过`operator=`从`number`转换过来。
     *
     */
    template <typename Number2>
    void
    value(const Number2      x,
          const unsigned int n_derivatives,
          Number2 *          values) const;

    /**
     * 多项式的程度。这是由构造函数提供的系数数所反映的程度。领先的非零系数不被单独处理。
     *
     */
    unsigned int
    degree() const;

    /**
     * 缩放多项式的底线。 给出多项式<i>p(t)</i>和缩放<i>t =
     * ax</i>，那么这个操作的结果就是多项式<i>q</i>，如<i>q(x)
     * = p(t)</i>。        该操作是在原地进行的。
     *
     */
    void
    scale(const number factor);

    /**
     * 将多项式的标点移出。 给出多项式<i>p(t)</i>和移位<i>t
     * = x +
     * a</i>，那么这个操作的结果就是多项式<i>q</i>，比如<i>q(x)
     * = p(t)</i>。
     * 模板参数允许以更高的精度计算新的系数，因为所有的计算都以<tt>number2</tt>类型进行。这可能是必要的，因为这个操作涉及大量的加法。在装有Solaris
     * 2.8的Sun Sparc Ultra上，<tt>double</tt>和<tt>long
     * double</tt>之间的差异并不明显，不过。
     * 操作是就地进行的，也就是说，现在对象的系数被改变。
     *
     */
    template <typename number2>
    void
    shift(const number2 offset);

    /**
     * 计算一个多项式的导数。
     *
     */
    Polynomial<number>
    derivative() const;

    /**
     * 计算一个多项式的基数。多项式的零阶项的系数为零。
     *
     */
    Polynomial<number>
    primitive() const;

    /**
     * 与一个标量相乘。
     *
     */
    Polynomial<number> &
    operator*=(const double s);

    /**
     * 与另一个多项式相乘。
     *
     */
    Polynomial<number> &
    operator*=(const Polynomial<number> &p);

    /**
     * 添加第二个多项式。
     *
     */
    Polynomial<number> &
    operator+=(const Polynomial<number> &p);

    /**
     * 减去第二个多项式。
     *
     */
    Polynomial<number> &
    operator-=(const Polynomial<number> &p);

    /**
     * 测试两个多项式是否相等。
     *
     */
    bool
    operator==(const Polynomial<number> &p) const;

    /**
     * 打印系数。
     *
     */
    void
    print(std::ostream &out) const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * 返回此对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const;

  protected:
    /**
     * 这个函数执行实际的缩放操作。
     *
     */
    static void
    scale(std::vector<number> &coefficients, const number factor);

    /**
     * 这个函数执行实际的移位
     *
     */
    template <typename number2>
    static void
    shift(std::vector<number> &coefficients, const number2 shift);

    /**
     * 将多项式乘以一个系数。
     *
     */
    static void
    multiply(std::vector<number> &coefficients, const number factor);

    /**
     * 将线性因子的乘积的多项式形式转化为标准形式，
     * $\sum_i a_i x^i$
     * 。删除所有与乘积形式有关的数据结构。
     *
     */
    void
    transform_into_standard_form();

    /**
     * 多项式的系数  $\sum_i a_i x^i$
     * 。这个向量由这个类的构造函数填充，并且可以由派生类传递下去。
     * 这个向量不能是常数，因为我们要允许复制多项式。
     *
     */
    std::vector<number> coefficients;

    /**
     * 存储该多项式是否为拉格朗日乘积形式，即构造为乘积
     * $(x-x_0) (x-x_1) \ldots (x-x_n)/c$ ，或者不是。
     *
     */
    bool in_lagrange_product_form;

    /**
     * 如果多项式是拉格朗日积形式，即构造为积 $(x-x_0)
     * (x-x_1) \ldots (x-x_n)/c$ ，则存储移位 $x_i$  。
     *
     */
    std::vector<number> lagrange_support_points;

    /**
     * 如果多项式是拉格朗日积的形式，即构造为积 $(x-x_0)
     * (x-x_1) \ldots (x-x_n)/c$ ，则存储权重c。
     *
     */
    number lagrange_weight;
  };


  /**
   * 类产生代表n度单项式的多项式对象，即函数  $x^n$  。
   *
   */
  template <typename number>
  class Monomial : public Polynomial<number>
  {
  public:
    /**
     * 构造函数，以单项式的度数和可选的系数作为参数。
     *
     */
    Monomial(const unsigned int n, const double coefficient = 1.);

    /**
     * 返回一个零度到<tt>度</tt>的单项式对象的向量，然后横跨到给定度数的全部多项式空间。这个函数可以用来初始化TensorProductPolynomials和PolynomialSpace类。
     *
     */
    static std::vector<Polynomial<number>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * 构造函数所需要的。
     *
     */
    static std::vector<number>
    make_vector(unsigned int n, const double coefficient);
  };


  /**
   * 拉格朗日多项式，其等距插值点在[0,1]。度数为<tt>n</tt>的多项式得到<tt>n+1</tt>的插值点。这些插值点按升序排列。这个顺序给每个插值点一个索引。
   * 一个拉格朗日多项式在其 "支持点
   * "等于1，在所有其他插值点等于0。例如，如果度数是3，支持点是1，那么这个对象所代表的多项式是立方的，它在<tt>x=1/3</tt>点的值是1，在<tt>x=0</tt>、<tt>x=2/3</tt>和<tt>x=1</tt>点是0。所有的多项式的度数都等于<tt>度数</tt>，但它们一起跨越了度数小于或等于<tt>度数</tt>的整个多项式空间。
   * 拉格朗日多项式的实现最高可达10度。
   *
   */
  class LagrangeEquidistant : public Polynomial<double>
  {
  public:
    /**
     * 构造函数。取拉格朗日多项式的度数<tt>n</tt>和支持点的索引<tt>support_point</tt>。填充基类多项式的<tt>系数</tt>。
     *
     */
    LagrangeEquidistant(const unsigned int n, const unsigned int support_point);

    /**
     * 返回一个度数为<tt>degree</tt>的多项式对象的向量，然后横跨到给定度数的全部多项式空间。多项式是通过调用这个类的构造函数生成的，其度数相同，但支持点从零到<tt>度</tt>。
     * 这个函数可以用来初始化TensorProductPolynomials和PolynomialSpace类。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * 计算基类Polynomial的<tt>系数</tt>。这个函数是<tt>静态的</tt>，允许在构造函数中调用。
     *
     */
    static void
    compute_coefficients(const unsigned int   n,
                         const unsigned int   support_point,
                         std::vector<double> &a);
  };



  /**
   * 给定一组沿实数轴的点，该函数返回用于插值这些点的所有拉格朗日多项式。多项式的数量等于点的数量，最大的度数是少一个。
   *
   */
  std::vector<Polynomial<double>>
  generate_complete_Lagrange_basis(const std::vector<Point<1>> &points);



  /**
   * 任意度的Legendre多项式。构造一个学位为<tt>p</tt>的Legendre多项式,
   * 根将由各自的点数的高斯公式和多项式的根的表示法来计算。
   * @note
   * 这类定义的多项式与通常所说的Legendre多项式在两个方面有所不同。(i)
   * 该类在参考区间  $[0,1]$
   * 上定义它们，而不是常用的区间  $[-1,1]$  。(ii)
   * 多项式的缩放方式使得它们是正交的，而不仅仅是正交的；因此，多项式的边界值不一定等于1。
   *
   */
  class Legendre : public Polynomial<double>
  {
  public:
    /**
     * 度数为<tt>p</tt>的多项式的构造函数。
     *
     */
    Legendre(const unsigned int p);

    /**
     * 返回一个零度到<tt>度</tt>的Legendre多项式对象的向量，然后横跨到给定度数的全部多项式空间。这个函数可用于初始化TensorProductPolynomials和PolynomialSpace类。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);
  };

  /**
   * <tt>[0,1]</tt>上任意程度的Lobatto多项式。
   * 这些多项式是[0,1]上的综合Legendre多项式。前两个多项式是由
   * $l_0(x) = 1-x$  和  $l_1(x) = x$
   * 给出的标准线性形状函数。对于  $i\geq2$  我们使用定义
   * $l_i(x) = \frac{1}{\Vert L_{i-1}\Vert_2}\int_0^x L_{i-1}(t)\,dt$
   * ，其中  $L_i$  表示  $i$  -次 Legendre 多项式在  $[0,1]$
   * 。Lobatto多项式 $l_0,\ldots,l_k$ 构成度数 $k$
   * 的多项式空间的完整基础。
   * 以给定的索引<tt>k</tt>调用构造函数将生成索引为<tt>k</tt>的多项式。但是只有对于
   * $k\geq 1$
   * ，索引等于多项式的度数。对于<tt>k==0</tt>也会生成一个1度的多项式。
   * 这些多项式被用于构建任意阶数的N&eacute;d&eacute;lec元素的形状函数。
   *
   */
  class Lobatto : public Polynomial<double>
  {
  public:
    /**
     * 度数为<tt>p</tt>的多项式的构造器。对于<tt>p==0</tt>有一个例外，见一般文档。
     *
     */
    Lobatto(const unsigned int p = 0);

    /**
     * 返回索引为<tt>0</tt>的多项式，直至<tt>度</tt>。
     * 对于<tt>p==0</tt>有一个例外，见一般文档。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int p);

  private:
    /**
     * 递归地计算系数。
     *
     */
    std::vector<double>
    compute_coefficients(const unsigned int p);
  };



  /**
   * 在<tt>[0,1]</tt>上任意程度的层次多项式。
   * 当构造一个度数为<tt>p</tt>的层次多项式时，系数将通过递归公式计算出来。
   * 系数被存储在一个静态数据向量中，以便下次需要时可以使用。
   * 这些分层多项式是基于Demkowicz, Oden,
   * Rachowicz和Hardy的多项式（CMAME 77 (1989) 79-112,
   * 第4节）。前两个多项式是由  $\phi_{0}(x) = 1
   *
   * - x$  和  $\phi_{1}(x) = x$  给出的标准线性形状函数。对于
   * $l \geq 2$  我们使用定义  $\phi_{l}(x) = (2x-1)^l
   *
   * - 1, l = 2,4,6,...$  和  $\phi_{l}(x) = (2x-1)^l
   *
   * - (2x-1), l = 3,5,7,...$  。这些满足递归关系  $\phi_{l}(x) =
   * (2x-1)\phi_{l-1}, l=3,5,7,...$  和  $\phi_{l}(x) = (2x-1)\phi_{l-1} +
   * \phi_{2}, l=4,6,8,...$  。
   * 自由度是顶点的值和中点的导数。目前，我们没有以任何方式对多项式进行缩放，尽管通过缩放可以实现对元素刚度矩阵更好的调节。
   * 以给定的索引<tt>p</tt>调用构造函数将产生以下结果：如果<tt>p==0</tt>，那么产生的多项式是与左边顶点相关的线性函数，如果<tt>p==1</tt>是与右边顶点相关的函数。对于更高的<tt>p</tt>值，你会得到与之前所有的正交的<tt>p</tt>度的多项式。注意，对于<tt>p==0</tt>，你因此确实<b>not</b>得到一个零度的多项式，但得到一个一度的多项式。这是为了允许生成一个完整的多项式空间的基础，只需迭代给构造函数的指数。
   * 另一方面，函数 generate_complete_basis()
   * 创建了一个给定度数的完整基。为了与多项式程度的概念相一致，如果给定的参数是零，它不会<b>not</b>返回上述的线性多项式，而是返回一个常数多项式。
   *
   */
  class Hierarchical : public Polynomial<double>
  {
  public:
    /**
     * 度数为<tt>p</tt>的多项式的构造函数。对于<tt>p==0</tt>有一个例外，见一般文档。
     *
     */
    Hierarchical(const unsigned int p);

    /**
     * 返回一个零度到<tt>degree</tt>的Hierarchical多项式对象的向量，然后横跨到给定度数的全部多项式空间。注意，如果给定的<tt>度数</tt>等于零，会有一个例外，请看这个类的一般文档。
     * 这个函数可以用来初始化TensorProductPolynomials,
     * AnisotropicPolynomials和PolynomialSpace类。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);

  private:
    /**
     * 递归地计算系数。
     *
     */
    static void
    compute_coefficients(const unsigned int p);

    /**
     * 为构造函数获取系数。
     * 这样，它就可以使用Polynomial的非标准构造函数。
     *
     */
    static const std::vector<double> &
    get_coefficients(const unsigned int p);

    /**
     * 带有已经计算过的系数的向量。对于多项式的每一个度数，我们保留一个指向系数列表的指针；我们这样做而不是保留一个向量的向量，以便简化多线程安全编程。为了避免内存泄漏，我们使用一个唯一的指针，以便在调用全局析构器时正确释放向量的内存。
     *
     */
    static std::vector<std::unique_ptr<const std::vector<double>>>
      recursive_coefficients;
  };



  /**
   * 用于Hermite插值条件的多项式。
   * 这是至少三度的多项式的集合，这样可以满足以下插值条件：多项式和它们的一阶导数在值<i>x</i>处消失。
   * ]=0和<i>x</i>=1，但<i>p</i><sub>0</sub>(0)=1，<i>p</i><sub><i>1</i></sub>(1)=1，<i>p</i>'<sub>2</sub>(0)=1，<<i>p'</i><sub>3</sub>(1)=1。
   * 对于三度，我们得到标准的四赫米特插值多项式，例如，见<a
   * href="http://en.wikipedia.org/wiki/Cubic_Hermite_spline">Wikipedia</a>。
   * 对于更高的度数，首先通过在<i>x</i>=0和<i>x</i>=1处具有消失值和导数的四度多项式进行增强，然后通过这个四阶多项式与增阶Legendre多项式的乘积进行增强。实现方法是
   * @f{align*}{
   * p_0(x) &= 2x^3-3x^2+1 \\
   * p_1(x) &=
   *
   * -2x^3+3x^2 \\
   * p_2(x) &= x^3-2x^2+x  \\
   * p_3(x) &= x^3-x^2 \\
   * p_4(x) &= 16x^2(x-1)^2 \\
   * \ldots & \ldots \\
   * p_k(x) &= x^2(x-1)^2 L_{k-4}(x)
   * @f}
   *
   *
   */
  class HermiteInterpolation : public Polynomial<double>
  {
  public:
    /**
     * 多项式的构造函数，索引为<tt>p</tt>。参见类文件中关于多项式序列的定义。
     *
     */
    HermiteInterpolation(const unsigned int p);

    /**
     * 返回指数<tt>0</tt>到<tt>p+1</tt>的多项式在一个度数不超过<tt>p</tt>的空间中。这里，<tt>p</tt>必须至少是3。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int p);
  };



  /**
   * Hermite多项式的一个变体，在插值中的条件数比HermiteInterpolation的基础更好。
   * 与适当的Hermite多项式相类似，这个基础在 $p_0$
   * 处将第一个多项式 $x=0$ 评估为1，在 $x=1$
   * 处有一个零值和零导数。同样地，最后一个多项式  $p_n$
   * 在  $x=1$  处评估为1，在  $x=0$
   * 处有零值和零导数。第二个多项式 $p_1$
   * 和倒数第二个多项式 $p_{n-1}$ 分别代表 $x=0$ 和 $x=1$
   * 处的导数自由度。  它们在两个端点 $x=0, x=1$
   * 处均为零，在另一端 $p_1'(1)=0$ 和 $p_{n-1}'(0)=0$
   * 处导数为零。与原来的Hermite多项式不同， $p_0$ 在 $x=0$
   * 处没有零导数。额外的自由度被用来使 $p_0$ 和 $p_1$
   * 正交，对于 $n=3$ ，其结果是 $p_0$ 的根在 $x=\frac{2}{7}$ ，
   * $p_n$ 的根在 $x=\frac{5}{7}$
   * ，分别。此外，这些多项式扩展到更高的度数 $n>3$
   * 是通过在单位区间内添加额外的节点来构建的，再次确保了更好的调节。这些节点被计算为
   * $\alpha=\beta=4$ 的雅可比多项式的根，它与生成函数
   * $x^2(1-x)^2$
   * 的平方是正交的，具有Hermite特性。然后，这些多项式以通常的方式构建为拉格朗日多项式，其双根在
   * $x=0$  和  $x=1$  。例如，有了 $n=4$ ，所有的 $p_0, p_1, p_3,
   * p_4$ 通过因子 $(x-0.5)$ 在 $x=0.5$
   * 处得到一个额外的根。综上所述，这个基础是由节点贡献主导的，但它不是一个节点的基础，因为第二个和倒数第二个多项式是非节点的，而且由于
   * $x=0$  和  $x=1$
   * 中存在双节点。基函数的权重设定为：所有具有单位权重的多项式之和代表常数函数1，与Lagrange多项式类似。
   * 基础只包含 <code>degree>=3</code>
   * 的Hermite信息，但对于0到2之间的度数也可以实现。对于线性情况，实现了通常的帽子函数，而对于
   * <code>degree=2</code> 的多项式是 $p_0(x)=(1-x)^2$ ，
   * $p_1(x)=2x(x-1)$ ，和 $p_2(x)=x^2$
   * ，根据3度的构造原理。这两个放松明显改善了质量矩阵的条件数（即内插），从下表可以看出这一点。
   * <table align="center" border="1"> <tr> <th>&nbsp;</th> <th
   * colspan="2">Condition number mass matrix</th> </tr> <tr> <th>degree</th>
   * <th>HermiteInterpolation</th> <th>HermiteLikeInterpolation</th> </tr>
   * <tr> <th>n=3</th> <th>1057</th> <th>17.18</th> </tr> <tr> <th>n=4</th>
   * <th>6580</th> <th>16.83</th> </tr> <tr> <th>n=5</th> <th>1.875e+04</th>
   * <th>15.99</th> </tr> <tr> <th>n=6</th> <th>6.033e+04</th> <th>16.34</th>
   * </tr> <tr> <th>n=10</th> <th>9.756e+05</th> <th>20.70</th> </tr> <tr>
   * <th>n=15</th> <th>9.431e+06</th> <th>27.91</th> </tr> <tr> <th>n=25</th>
   * <th>2.220e+08</th> <th>43.54</th> </tr> <tr> <th>n=35</th>
   * <th>2.109e+09</th> <th>59.51</th> </tr> </table>
   * 这个多项式继承了Hermite多项式的有利特性，在一个面上只有两个函数的值和/或导数不为零，这对不连续Galerkin方法是有利的，但却给出了更好的插值条件数，这改善了一些迭代方案的性能，如共轭梯度与点-Jacobi。这个多项式在FE_DGQHermite中使用。
   *
   */
  class HermiteLikeInterpolation : public Polynomial<double>
  {
  public:
    /**
     * 在设定的度数的多项式中，索引为<tt>index</tt>的多项式的构造函数
     * @p degree.  。
     *
     */
    HermiteLikeInterpolation(const unsigned int degree,
                             const unsigned int index);

    /**
     * 返回索引为<tt>0</tt>至<tt>度+1</tt>的多项式在度数为<tt>度</tt>的空间中的情况。
     *
     */
    static std::vector<Polynomial<double>>
    generate_complete_basis(const unsigned int degree);
  };



  /* 评估由参数 @p alpha,   @p beta,   @p n, 指定的雅可比多项式 $ P_n^{\alpha, \beta}(x) $ ，其中 @p n 是雅可比多项式的程度。   
*  @note  雅可比多项式不是正交的，像通常的deal.II一样定义在单位区间 $[0, 1]$ ，而不是文献中经常使用的 $[-1, +1]$ 。  @p x 是评价的点。 
*
*/
  template <typename Number>
  Number
  jacobi_polynomial_value(const unsigned int degree,
                          const int          alpha,
                          const int          beta,
                          const Number       x);


  /**
   * 计算单位区间 $[0, 1]$
   * 上给定度数的雅可比多项式的根。这些根在deal.II库中有多处使用，如Gauss-Lobatto正交公式或用于Hermite-like插值。
   * 该算法使用牛顿算法，使用切比雪夫多项式的零点作为初始猜测。这段代码已经在α和β等于0（Legendre情况）、1（Gauss-Lobatto情况）以及2的情况下进行了测试，所以在对其他数值使用时要小心，因为牛顿迭代可能会或可能不会收敛。
   *
   */
  template <typename Number>
  std::vector<Number>
  jacobi_polynomial_roots(const unsigned int degree,
                          const int          alpha,
                          const int          beta);
} // namespace Polynomials


 /** @} */ 

 /* -------------------------- inline functions --------------------- */ 

namespace Polynomials
{
  template <typename number>
  inline Polynomial<number>::Polynomial()
    : in_lagrange_product_form(false)
    , lagrange_weight(1.)
  {}



  template <typename number>
  inline unsigned int
  Polynomial<number>::degree() const
  {
    if (in_lagrange_product_form == true)
      {
        return lagrange_support_points.size();
      }
    else
      {
        Assert(coefficients.size() > 0, ExcEmptyObject());
        return coefficients.size() - 1;
      }
  }



  template <typename number>
  inline number
  Polynomial<number>::value(const number x) const
  {
    if (in_lagrange_product_form == false)
      {
        Assert(coefficients.size() > 0, ExcEmptyObject());

        // Horner scheme
        const unsigned int m     = coefficients.size();
        number             value = coefficients.back();
        for (int k = m - 2; k >= 0; --k)
          value = value * x + coefficients[k];
        return value;
      }
    else
      {
        // direct evaluation of Lagrange polynomial
        const unsigned int m     = lagrange_support_points.size();
        number             value = 1.;
        for (unsigned int j = 0; j < m; ++j)
          value *= x - lagrange_support_points[j];
        value *= lagrange_weight;
        return value;
      }
  }



  template <typename number>
  template <typename Number2>
  inline void
  Polynomial<number>::value(const Number2      x,
                            const unsigned int n_derivatives,
                            Number2 *          values) const
  {
    // evaluate Lagrange polynomial and derivatives
    if (in_lagrange_product_form == true)
      {
        // to compute the value and all derivatives of a polynomial of the
        // form (x-x_1)*(x-x_2)*...*(x-x_n), expand the derivatives like
        // automatic differentiation does.
        const unsigned int n_supp = lagrange_support_points.size();
        const number       weight = lagrange_weight;
        switch (n_derivatives)
          {
            default:
              values[0] = 1.;
              for (unsigned int d = 1; d <= n_derivatives; ++d)
                values[d] = 0.;
              for (unsigned int i = 0; i < n_supp; ++i)
                {
                  const Number2 v = x - lagrange_support_points[i];

                  // multiply by (x-x_i) and compute action on all derivatives,
                  // too (inspired from automatic differentiation: implement the
                  // product rule for the old value and the new variable 'v',
                  // i.e., expand value v and derivative one). since we reuse a
                  // value from the next lower derivative from the steps before,
                  // need to start from the highest derivative
                  for (unsigned int k = n_derivatives; k > 0; --k)
                    values[k] = (values[k] * v + values[k - 1]);
                  values[0] *= v;
                }
              // finally, multiply by the weight in the Lagrange
              // denominator. Could be done instead of setting values[0] = 1
              // above, but that gives different accumulation of round-off
              // errors (multiplication is not associative) compared to when we
              // computed the weight, and hence a basis function might not be
              // exactly one at the center point, which is nice to have. We also
              // multiply derivatives by k! to transform the product p_n =
              // p^(n)(x)/k! into the actual form of the derivative
              {
                number k_factorial = 1;
                for (unsigned int k = 0; k <= n_derivatives; ++k)
                  {
                    values[k] *= k_factorial * weight;
                    k_factorial *= static_cast<number>(k + 1);
                  }
              }
              break;

            // manually implement case 0 (values only), case 1 (value + first
            // derivative), and case 2 (up to second derivative) since they
            // might be called often. then, we can unroll the inner loop and
            // keep the temporary results as local variables to help the
            // compiler with the pointer aliasing analysis.
            case 0:
              {
                Number2 value = 1.;
                for (unsigned int i = 0; i < n_supp; ++i)
                  {
                    const Number2 v = x - lagrange_support_points[i];
                    value *= v;
                  }
                values[0] = weight * value;
                break;
              }

            case 1:
              {
                Number2 value      = 1.;
                Number2 derivative = 0.;
                for (unsigned int i = 0; i < n_supp; ++i)
                  {
                    const Number2 v = x - lagrange_support_points[i];
                    derivative      = derivative * v + value;
                    value *= v;
                  }
                values[0] = weight * value;
                values[1] = weight * derivative;
                break;
              }

            case 2:
              {
                Number2 value      = 1.;
                Number2 derivative = 0.;
                Number2 second     = 0.;
                for (unsigned int i = 0; i < n_supp; ++i)
                  {
                    const Number2 v = x - lagrange_support_points[i];
                    second          = second * v + derivative;
                    derivative      = derivative * v + value;
                    value *= v;
                  }
                values[0] = weight * value;
                values[1] = weight * derivative;
                values[2] = static_cast<number>(2) * weight * second;
                break;
              }
          }
        return;
      }

    Assert(coefficients.size() > 0, ExcEmptyObject());

    // if derivatives are needed, then do it properly by the full
    // Horner scheme
    const unsigned int   m = coefficients.size();
    std::vector<Number2> a(coefficients.size());
    std::copy(coefficients.begin(), coefficients.end(), a.begin());
    unsigned int j_factorial = 1;

    // loop over all requested derivatives. note that derivatives @p{j>m} are
    // necessarily zero, as they differentiate the polynomial more often than
    // the highest power is
    const unsigned int min_valuessize_m = std::min(n_derivatives + 1, m);
    for (unsigned int j = 0; j < min_valuessize_m; ++j)
      {
        for (int k = m - 2; k >= static_cast<int>(j); --k)
          a[k] += x * a[k + 1];
        values[j] = static_cast<number>(j_factorial) * a[j];

        j_factorial *= j + 1;
      }

    // fill higher derivatives by zero
    for (unsigned int j = min_valuessize_m; j <= n_derivatives; ++j)
      values[j] = 0.;
  }



  template <typename number>
  template <class Archive>
  inline void
  Polynomial<number>::serialize(Archive &ar, const unsigned int)
  {
    // forward to serialization function in the base class.
    ar &static_cast<Subscriptor &>(*this);
    ar &coefficients;
    ar &in_lagrange_product_form;
    ar &lagrange_support_points;
    ar &lagrange_weight;
  }



  template <typename Number>
  Number
  jacobi_polynomial_value(const unsigned int degree,
                          const int          alpha,
                          const int          beta,
                          const Number       x)
  {
    Assert(alpha >= 0 && beta >= 0,
           ExcNotImplemented("Negative alpha/beta coefficients not supported"));
    // the Jacobi polynomial is evaluated using a recursion formula.
    Number p0, p1;

    // The recursion formula is defined for the interval [-1, 1], so rescale
    // to that interval here
    const Number xeval = Number(-1) + 2. * x;

    // initial values P_0(x), P_1(x):
    p0 = 1.0;
    if (degree == 0)
      return p0;
    p1 = ((alpha + beta + 2) * xeval + (alpha - beta)) / 2;
    if (degree == 1)
      return p1;

    for (unsigned int i = 1; i < degree; ++i)
      {
        const Number v  = 2 * i + (alpha + beta);
        const Number a1 = 2 * (i + 1) * (i + (alpha + beta + 1)) * v;
        const Number a2 = (v + 1) * (alpha * alpha - beta * beta);
        const Number a3 = v * (v + 1) * (v + 2);
        const Number a4 = 2 * (i + alpha) * (i + beta) * (v + 2);

        const Number pn = ((a2 + a3 * xeval) * p1 - a4 * p0) / a1;
        p0              = p1;
        p1              = pn;
      }
    return p1;
  }



  template <typename Number>
  std::vector<Number>
  jacobi_polynomial_roots(const unsigned int degree,
                          const int          alpha,
                          const int          beta)
  {
    std::vector<Number> x(degree, 0.5);

    // compute zeros with a Newton algorithm.

    // Set tolerance. For long double we might not always get the additional
    // precision in a run time environment (e.g. with valgrind), so we must
    // limit the tolerance to double. Since we do a Newton iteration, doing
    // one more iteration after the residual has indicated convergence will be
    // enough for all number types due to the quadratic convergence of
    // Newton's method

    const Number tolerance =
      4 * std::max(static_cast<Number>(std::numeric_limits<double>::epsilon()),
                   std::numeric_limits<Number>::epsilon());

    // The following implementation follows closely the one given in the
    // appendix of the book by Karniadakis and Sherwin: Spectral/hp element
    // methods for computational fluid dynamics (Oxford University Press,
    // 2005)

    // If symmetric, we only need to compute the half of points
    const unsigned int n_points = (alpha == beta ? degree / 2 : degree);
    for (unsigned int k = 0; k < n_points; ++k)
      {
        // we take the zeros of the Chebyshev polynomial (alpha=beta=-0.5) as
        // initial values, corrected by the initial value
        Number r = 0.5 - 0.5 * std::cos(static_cast<Number>(2 * k + 1) /
                                        (2 * degree) * numbers::PI);
        if (k > 0)
          r = (r + x[k - 1]) / 2;

        unsigned int converged = numbers::invalid_unsigned_int;
        for (unsigned int it = 1; it < 1000; ++it)
          {
            Number s = 0.;
            for (unsigned int i = 0; i < k; ++i)
              s += 1. / (r - x[i]);

            // derivative of P_n^{alpha,beta}, rescaled to [0, 1]
            const Number J_x =
              (alpha + beta + degree + 1) *
              jacobi_polynomial_value(degree - 1, alpha + 1, beta + 1, r);

            // value of P_n^{alpha,beta}
            const Number f = jacobi_polynomial_value(degree, alpha, beta, r);
            const Number delta = f / (f * s - J_x);
            r += delta;
            if (converged == numbers::invalid_unsigned_int &&
                std::abs(delta) < tolerance)
              converged = it;

            // do one more iteration to ensure accuracy also for tighter
            // types than double (e.g. long double)
            if (it == converged + 1)
              break;
          }

        Assert(converged != numbers::invalid_unsigned_int,
               ExcMessage("Newton iteration for zero of Jacobi polynomial "
                          "did not converge."));

        x[k] = r;
      }

    // in case we assumed symmetry, fill up the missing values
    for (unsigned int k = n_points; k < degree; ++k)
      x[k] = 1.0 - x[degree - k - 1];

    return x;
  }

} // namespace Polynomials
DEAL_II_NAMESPACE_CLOSE

#endif


