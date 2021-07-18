//include/deal.II-translator/base/quadrature_lib_0.txt
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

#ifndef dealii_quadrature_lib_h
#define dealii_quadrature_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Quadrature */ 
 /*@{*/ 

/**
 * 用于数值积分的Gauss-Legendre系列的正交规则。
 * 这些正交规则的系数是由<a
 * href="http://en.wikipedia.org/wiki/Numerical_Recipes">Numerical
 * Recipes</a>中描述的函数计算的。
 *
 *
 */
template <int dim>
class QGauss : public Quadrature<dim>
{
public:
  /**
   * 生成一个有<tt>n</tt>个正交点（在每个空间方向）的公式，对<tt>2n-1</tt>度的多项式是精确的。
   *
   */
  QGauss(const unsigned int n);
};


/**
 * 高斯-洛巴托系列的正交规则用于数值积分。
 * 这种对高斯正交的修改也使用了两个区间的端点。对于度数为<i>2n-3</i>的多项式是精确的，这个公式有两个度数的次优。
 * 正交点是区间端点加上度数为<i>n-1</i>的Legendre多项式<i>P<sub>n-1</sub></i>的导数的根。正交权重为<i>2/(n(n-1)(P<sub>n-1</sub>(x<sub>i</sub>)<sup>2</sup>)</i>。
 *
 *
 * @note
 * 这个实现没有在数值稳定性和效率方面进行优化。它可以很容易地适应具有任意参数的Gauss-Lobatto-Jacobi-Bouzitat正交的一般情况
 * $\alpha$  ,  $\beta$  ，其中Gauss-Lobatto-Legendre正交（  $\alpha =
 * \beta = 0$  ）是一个特殊情况。 @sa
 * http://en.wikipedia.org/wiki/Handbook_of_Mathematical_Functions  @sa
 * Karniadakis, G.E. and Sherwin, S.J.: Spectral/HP element methods for
 * computational fluid dynamics. 牛津。牛津大学出版社，2005
 *
 */
template <int dim>
class QGaussLobatto : public Quadrature<dim>
{
public:
  /**
   * 生成一个具有<tt>n</tt>正交点的公式（在每个空间方向）。
   *
   */
  QGaussLobatto(const unsigned int n);
};



/**
 * 数字正交的中点规则。这个单点公式对线性多项式是精确的。
 *
 *
 */
template <int dim>
class QMidpoint : public Quadrature<dim>
{
public:
  QMidpoint();
};


/**
 * 数字正交的辛普森规则。这个有3个正交点的公式对3度的多项式是精确的。
 *
 *
 */
template <int dim>
class QSimpson : public Quadrature<dim>
{
public:
  QSimpson();
};



/**
 * 数字正交的梯形规则。这个有两个正交点的公式对线性多项式是精确的，并使用区间的端点进行1d中的函数评估，见https://en.wikipedia.org/wiki/Trapezoidal_rule
 * 。在更高的维度上，该类是通过张量积构建的，然后使用四边形或六面体的顶点进行函数评估。
 *
 *
 */
template <int dim>
class QTrapezoid : public Quadrature<dim>
{
public:
  QTrapezoid();
};


/**
 * QTrapezoid的一个别名，由于历史原因可用。这个名字已经被废弃了。
 * 该类最初被命名为QTrapez，这是一个糟糕的命名选择，因为正交公式的正确名称是
 * "梯形规则"，或者有时也被称为
 * "梯形规则"。这个错误的名字是由于它的原作者的英语水平很差，导致他们错误地从德语
 * "Trapezregel "翻译了这个名字。
 *
 *
 */
template <int dim>
using QTrapez DEAL_II_DEPRECATED = QTrapezoid<dim>;



/**
 * 数字正交公式的米尔恩规则。Milne规则是一个封闭的Newton-Cotes公式，对于5度的多项式是精确的。
 * @sa  Stoer: Einführung in die Numerische Mathematik I, p. 102
 *
 */
template <int dim>
class QMilne : public Quadrature<dim>
{
public:
  QMilne();
};


/**
 * 数值正交的韦德尔规则。Weddle规则是一个封闭的Newton-Cotes公式，对于7度的多项式是精确的。
 * @sa  Stoer: Einführung in die Numerische Mathematik I, p. 102
 *
 */
template <int dim>
class QWeddle : public Quadrature<dim>
{
public:
  QWeddle();
};



/**
 * 带对数加权函数的高斯正交的一个类。这个公式用于在区间
 * $[0,1]$ 上积分 $\ln|x|\;f(x)$ ，其中 $f$
 * 是一个没有奇点的平滑函数。正交点和权重的集合已经用<tt>数字配方</tt>得到。
 * 注意只应提供函数 $f(x)$ ，即 $\int_0^1 f(x) \ln|x| dx =
 * \sum_{i=0}^N w_i f(q_i)$  。在构造时将 @p revert 标志设置为
 * "真"，可以将权重从 $\ln|x|$ 切换到 $\ln|1-x|$  。
 * 权重和函数已经表到了12阶。
 *
 *
 */
template <int dim>
class QGaussLog : public Quadrature<dim>
{
public:
  /**
   * 产生一个具有<tt>n</tt>正交点的公式
   *
   */
  QGaussLog(const unsigned int n, const bool revert = false);

private:
  /**
   * 计算正交公式的点。
   *
   */
  static std::vector<double>
  get_quadrature_points(const unsigned int n);

  /**
   * 计算正交公式的权重。
   *
   */
  static std::vector<double>
  get_quadrature_weights(const unsigned int n);
};



/**
 * 一个用于高斯正交的类，具有任意的对数加权函数。这个公式用于在区间
 * $[0,1]$ 上对 $\ln(|x-x_0|/\alpha)\;f(x)$ 进行积分，其中 $f$
 * 是一个没有奇点的平滑函数， $x_0$ 和 $\alpha$
 * 是在构造时给出的，是奇点的位置 $x_0$
 * 和奇点中一个任意的缩放因子。 你必须确保点 $x_0$
 * 不是阶数为 $N$
 * 的高斯正交点之一，否则会出现异常，因为正交权重不能被正确计算。
 * 这个正交公式相当昂贵，因为它在内部使用了两个n阶的高斯正交公式来积分因子的非辛格部分，以及两个GaussLog正交公式来积分独立的片段
 * $[0,x_0]$  和  $[x_0,1]$
 * 。如果奇点是其中一个极值，且因子α为1，那么这个正交就与QGaussLog相同。
 * 构造函数的最后一个参数允许你以两种可能的方式之一使用这个正交规则。\f[
 * \int_0^1 g(x) dx = \int_0^1 f(x) \ln\left(\frac{|x-x_0|}{\alpha}\right) dx
 * = \sum_{i=0}^N w_i g(q_i) = \sum_{i=0}^N \bar{w}_i f(q_i) \f]
 * 提供两组权重中的哪一组，可以通过 @p
 * factor_out_singular_weight参数选择。如果它是假的（默认），那么将计算
 * $\bar{w}_i$ 权重，你应该只提供平滑函数 $f(x)$
 * ，因为奇异点被包含在正交里面。如果参数设置为
 * "true"，那么奇异点将从正交公式中剔除，你应该提供一个函数
 * $g(x)$  ，它至少应该类似于  $\ln(|x-x_0|/\alpha)$  。
 * 请注意，如果你试图将这个正交规则用于常规函数，一旦你将奇异点剔除，这个正交规则就没有价值。
 * 权重和函数已经列到了12阶。
 *
 *
 */
template <int dim>
class QGaussLogR : public Quadrature<dim>
{
public:
  /**
   * 构造函数需要四个参数：高斯公式在每段 $[0,x_0]$ 和
   * $[x_0,1]$
   * 上的顺序，奇点的实际位置，对数函数内部的比例因子和一个标志，决定奇点是留在正交公式内还是被因子化掉，以包括在积分里。
   *
   */
  QGaussLogR(const unsigned int n,
             const Point<dim> & x0                         = Point<dim>(),
             const double       alpha                      = 1,
             const bool         factor_out_singular_weight = false);

  /**
   * 移动构造器。我们不能依赖`正交'的移动构造函数，因为它不知道这个类的额外成员`分数'。
   *
   */
  QGaussLogR(QGaussLogR<dim> &&) noexcept = default;

protected:
  /**
   * 这是区间 $(0,origin)$
   * 的长度，如果两个极端中的任何一个已经被选中，则为1。
   *
   */
  const double fraction;
};


/**
 * 一个用于高斯正交的类，带有 $1/R$
 * 加权函数。这个公式可以用来在参考元素 $[0,1]^2$ 上积分
 * $1/R \ f(x)$ ，其中 $f$ 是一个没有奇点的平滑函数， $R$
 * 是点 $x$ 到顶点 $\xi$
 * 的距离，在构造时通过指定其索引给出。请注意，这个距离是在参考元素中评估的。
 * 这个正交公式是由两个QGauss正交公式得到的，将它们转换为以奇点为中心的极坐标系，然后再转换为另一个参考元素。这使得奇点可以被转换的部分雅各布系数所抵消，其中包含
 * $R$
 * 。在实践中，参考元素通过折叠与奇点相邻的一条边而被转化为一个三角形。这个变换的Jacobian包含
 * $R$
 * ，在对原始正交进行缩放之前，这个Jacobian被去除，这个过程对下一个半元素重复进行。
 * 在构建时，可以指定我们是否要去除奇点。换句话说，这个正交可以用来整合
 * $g(x) = 1/R\ f(x)$ ，或者简单地整合 $f(x)$ ， $1/R$
 * 因子已经包含在正交的权重中。
 *
 *
 */
template <int dim>
class QGaussOneOverR : public Quadrature<dim>
{
public:
  /**
   * 这个构造函数需要三个参数：高斯公式的顺序，奇点所在的参考元素的点，以及我们是将加权奇点函数包含在正交里面，还是将其留在用户函数中进行积分。
   * 传统上，正交公式包括其加权函数，最后一个参数默认设置为假。然而，有些情况下这是不可取的（例如，当你只知道你的奇点具有1/R的相同阶数，但不能准确地以这种方式书写）。
   * 换句话说，你可以用以下两种方式使用这个函数，获得相同的结果。
   * @code
   * QGaussOneOverR singular_quad(order, q_point, false);
   * // This will produce the integral of f(x)/R
   * for(unsigned int i=0; i<singular_quad.size(); ++i)
   * integral += f(singular_quad.point(i))*singular_quad.weight(i);
   *
   * // And the same here
   * QGaussOneOverR singular_quad_noR(order, q_point, true);
   *
   * // This also will produce the integral of f(x)/R, but 1/R has to
   * // be specified.
   * for(unsigned int i=0; i<singular_quad.size(); ++i) {
   * double R = (singular_quad_noR.point(i)-cell->vertex(vertex_id)).norm();
   * integral += f(singular_quad_noR.point(i))*singular_quad_noR.weight(i)/R;
   * }
   * @endcode
   *
   *
   */
  QGaussOneOverR(const unsigned int n,
                 const Point<dim> & singularity,
                 const bool         factor_out_singular_weight = false);
  /**
   * 构造函数需要三个参数：高斯公式的顺序，奇点所在顶点的索引，以及我们是将加权奇点函数包含在正交里面，还是将它留在用户函数中进行积分。请注意，这是前一个构造函数的专门版本，只对四边形的顶点起作用。
   * 传统上，正交公式包括其加权函数，最后一个参数默认设置为假。然而，有些情况下这是不可取的（例如，当你只知道你的奇点具有1/R的相同阶数，但不能准确地以这种方式书写）。
   * 换句话说，你可以用以下两种方式使用这个函数，获得相同的结果。
   * @code
   * QGaussOneOverR singular_quad(order, vertex_id, false);
   * // This will produce the integral of f(x)/R
   * for(unsigned int i=0; i<singular_quad.size(); ++i)
   * integral += f(singular_quad.point(i))*singular_quad.weight(i);
   *
   * // And the same here
   * QGaussOneOverR singular_quad_noR(order, vertex_id, true);
   *
   * // This also will produce the integral of f(x)/R, but 1/R has to
   * // be specified.
   * for(unsigned int i=0; i<singular_quad.size(); ++i) {
   * double R = (singular_quad_noR.point(i)-cell->vertex(vertex_id)).norm();
   * integral += f(singular_quad_noR.point(i))*singular_quad_noR.weight(i)/R;
   * }
   * @endcode
   *
   *
   */
  QGaussOneOverR(const unsigned int n,
                 const unsigned int vertex_index,
                 const bool         factor_out_singular_weight = false);

private:
  /**
   * 给定一个正交点和一个度数n，这个函数返回奇异正交规则的大小，考虑到该点是否在单元格内，是否在单元格的边缘，是否在单元格的角落。
   *
   */
  static unsigned int
  quad_size(const Point<dim> &singularity, const unsigned int n);
};



/**
 * 排序的正交法。给定一个任意的正交公式，该类生成一个正交公式，其中正交点根据权重排序，从对应权重小的点到对应权重大的点。这可能是必要的，例如，在整合高阶多项式时，因为在这些情况下，你可能会用非常小的数字来求和，如果要求和的数字不相近，求和就不稳定。
 *
 *
 */
template <int dim>
class QSorted : public Quadrature<dim>
{
public:
  /**
   * 构造函数接受一个任意的正交公式 @p quad
   * ，并按照升序权重对其点和权重进行排序。
   *
   */
  QSorted(const Quadrature<dim> &quad);

private:
  /**
   * 一个用于 std::sort 的规则，对点和权重进行重新排序。
   * @p a 和 @p b
   * 是权重数组的索引，结果将通过比较权重来确定。
   *
   */
  bool
  compare_weights(const unsigned int a, const unsigned int b) const;
};

/**
 * 任意阶数的Telles正交。
 * 这些正交规则的系数是使用非线性的变量变化从Gauss-Legendre正交公式开始计算的。这是用一个三次多项式，
 * $n = a x^3 + b x^2 + c x + d$
 * 来完成的，以便整合一个奇异的积分，奇异点在给定的x_0。
 * 我们从一个具有任意函数的高斯正交公式开始。然后我们应用三维变量变化。在论文中，J.C.F.Telles：A
 * Self-Adaptive Co-ordinate Transformation For Efficient Numerical Evaluation
 * of General Boundary Element Integrals. 1987年，作者在参考单元
 * $[-1, 1]$ 上应用变换。
 *
 * @f{align*}{
 * n(1) &= 1, \\ n(-1) &=
 *
 * -1, \\ \frac{dn}{dx} &= 0 \text{ at }
 * x = x_0, \\ \frac{d^2n}{dx^2} &= 0 \text{ at  } x = x_0
 * @f}
 * 我们得到
 *
 * @f{align*}{
 * a &= \frac{1}{q}, \\
 * b &=
 *
 * -3 \frac{\bar{\Gamma}}{q}, \\
 * c &= 3 \frac{\bar{\Gamma}}{q}, \\
 * d &=
 *
 * -b,
 * @f}
 * 与
 * @f{align*}{
 * \eta^{*} &= \bar{\eta}^2
 *
 * - 1, \\
 * \bar{\Gamma}  &= \sqrt[3]{\bar{\eta} \eta^{*} + |\eta^{*} | }
 *                + \sqrt[3]{ \bar{\eta} \eta^{*}
 *
 * - |\eta^{*} | }
 *                + \bar{\eta}, \\
 * q &= (\Gamma-\bar{\Gamma})^3 + \bar{\Gamma}
 *    \frac{\bar{\Gamma}^2+3}{1+3\bar{\Gamma}^2}
 * @f}
 * 由于库中假设 $[0,1]$
 * 为参考区间，我们将在实现中把这些值映射到适当的参考区间上。
 * 这种变量变化可以用来整合奇异积分。一个例子是
 * $f(x)/|x-x_0|$ 在参考区间 $[0,1]$ 上，其中 $x_0$
 * 是在构造时给出的，是奇点 $x_0$ 的位置，而 $f(x)$
 * 是一个平滑的非奇异函数。
 * 奇异正交公式是相当昂贵的，然而Telles的正交公式相对于Lachat-Watson等其他奇异积分技术来说，更容易计算。
 * 我们已经实现了 $dim = 1$ 的情况。当我们处理 $dim >1$
 * 的情况时，我们已经计算出正交公式有一个一维Telles正交公式的张量乘积，考虑到奇点的不同组成部分。
 * 高斯Legendre公式的权重和函数已被制成表格，最高可达12阶。
 *
 *
 */
template <int dim>
class QTelles : public Quadrature<dim>
{
public:
  /**
   * 一个构造函数，接受一个正交公式和一个奇点作为参数。正交公式将使用Telles规则进行映射。请确保正交规则的顺序适合于有关的奇点。
   *
   */
  QTelles(const Quadrature<1> &base_quad, const Point<dim> &singularity);
  /**
   * 上述构造函数的一个变体，将阶数 @p n
   * 和奇点的位置作为参数。将使用n阶的高斯勒格朗德正交。
   *
   */
  QTelles(const unsigned int n, const Point<dim> &singularity);
};

/**
 * 高斯-切比雪夫正交规则整合了加权乘积 $\int_{-1}^1 f(x) w(x)
 * dx$ ，其权重由以下因素给出。  $w(x) = 1/\sqrt{1-x^2}$  .
 * 节点和权重是已知的分析结果，对单项式来说，精确到阶
 * $2n-1$  ，其中 $n$
 * 是正交点的数量。在这里，我们重新调整正交公式，使其定义在区间
 * $[0,1]$  而不是  $[-1,1]$  上。所以正交公式恰恰是对积分
 * $\int_0^1 f(x) w(x) dx$  进行了加权整合。  $w(x) = 1/\sqrt{x(1-x)}$
 * . 详见。M. Abramowitz & I.A. Stegun:
 * 数学函数手册》，第2页。25.4.38
 *
 *
 */
template <int dim>
class QGaussChebyshev : public Quadrature<dim>
{
public:
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussChebyshev(const unsigned int n);
};


/**
 * 高斯-拉道-切比雪夫正交规则整合了加权乘积 $\int_{-1}^1
 * f(x) w(x) dx$ ，其权重由。  $w(x) = 1/\sqrt{1-x^2}$
 * ，附加约束条件是正交点位于区间的两个极值之一。节点和权重是已知的分析结果，对单项式来说是精确的，直到阶
 * $2n-2$  ，其中 $n$
 * 是正交点的数量。在这里，我们重新调整正交公式，使其定义在区间
 * $[0,1]$ 而不是 $[-1,1]$ 。所以正交公式恰恰是对积分
 * $\int_0^1 f(x) w(x) dx$  进行了加权整合。  $w(x) = 1/\sqrt{x(1-x)}$
 * .
 * 默认情况下，正交公式是以左端点为正交节点构建的，但正交节点可以通过变量ep强加在右端点，变量ep可以向左或向右取值。
 *
 *
 */
template <int dim>
class QGaussRadauChebyshev : public Quadrature<dim>
{
public:
  /* EndPoint用于指定单位区间的两个端点中的哪一个也作为正交点使用。 
*
*/
  enum EndPoint
  {
    /**
     * 左侧端点。
     *
     */
    left,
    /**
     * 右端点。
     *
     */
    right
  };
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussRadauChebyshev(const unsigned int n,
                       EndPoint           ep = QGaussRadauChebyshev::left);

  /**
   * 移动构造器。我们不能依靠 "正交
   * "的移动构造函数，因为它不知道这个类的额外成员`ep'。
   *
   */
  QGaussRadauChebyshev(QGaussRadauChebyshev<dim> &&) noexcept = default;

private:
  const EndPoint ep;
};

/**
 * Gauss-Lobatto-Chebyshev正交规则整合了加权乘积 $\int_{-1}^1 f(x)
 * w(x) dx$ ，其权重由以下因素给出。  $w(x) = 1/\sqrt{1-x^2}$
 * ，附加约束条件是两个正交点位于正交区间的端点。节点和权重是已知的分析结果，对单项式来说是准确的，直到阶数
 * $2n-3$ ，其中 $n$
 * 是正交点的数量。在这里，我们重新调整正交公式，使其定义在区间
 * $[0,1]$  而不是  $[-1,1]$  上。所以正交公式恰恰是对积分
 * $\int_0^1 f(x) w(x) dx$  进行了加权整合。  $w(x) = 1/\sqrt{x(1-x)}$
 * . 详见。M. Abramowitz & I.A. Stegun:
 * 数学函数手册》，第2页。25.4.40
 *
 *
 */
template <int dim>
class QGaussLobattoChebyshev : public Quadrature<dim>
{
public:
  /// Generate a formula with <tt>n</tt> quadrature points
  QGaussLobattoChebyshev(const unsigned int n);
};

/**
 * 给定一个任意的正交公式，返回一个在 $\sum_i x_i = 1$
 * 所定义的超平面之上的正交点。换句话说，它从基础公式中提取那些满足
 * $\sum_i (\mathbf x_q)_i \le 1+10^{-12}$ 的正交点。"
 * 一般来说，所得到的正交点不是很有用，除非你开始的正交点是专门为在三角形或四面体上积分而构建的。该类仅确保所产生的正交公式仅在参考单轴或其边界上有正交点。
 * 没有对权重进行转换，提及参考单轴以外的点的权重被简单地丢弃了。
 * 这个正交公式的主要用途不是用来切张量积正交的。理想情况下，你应该向这个类传递一个直接使用参考单轴中的点和权重构造的正交公式，能够在三角形或四面体上积分。
 * 对于基于四边形和六面体的有限元，QSimplex正交公式本身并不十分有用。这个类通常和其他类一起使用，比如QSplit，用几个QSimplex正交公式来修补参考元素。
 * 这样的正交公式对于整合在某些点上有奇点的函数，或者在参考元素内部沿同维度表面出现跳跃的函数是非常有用的，比如在扩展有限元方法（XFEM）中。
 *
 *
 */
template <int dim>
class QSimplex : public Quadrature<dim>
{
public:
  /**
   * 构建一个只包含在左下角参考单轴的点的正交。
   * @param[in]  quad 输入的正交点。
   *
   */
  QSimplex(const Quadrature<dim> &quad);

  /**
   * 返回该正交的仿射变换，可用于在`vertices'标识的单轴上进行积分。
   * 正交点的位置和权重都被转换，因此你可以有效地使用所得到的正交点在单线上积分。
   * 转换定义为\f[ x = v_0 + B \hat x
   * \f]，其中矩阵 $B$ 由 $B_{ij} = v[j][i]-v[0][i]$ 给出。
   * 权重以 $B$ 的行列式的绝对值为尺度，即 $J \dealcoloneq
   * |\text{det}(B)|$ 。如果 $J$
   * 为零，将返回一个空的正交。这种情况可能会发生，在二维空间，如果三个顶点是对齐的，或者在三维空间，如果四个顶点在同一个平面上。
   * @param[in]  顶点 你希望在单线上积分的顶点  @return
   * 一个正交对象，可用于在单线上积分。
   *
   */
  Quadrature<dim>
  compute_affine_transformation(
    const std::array<Point<dim>, dim + 1> &vertices) const;
};

/**
 * 一个实现从正方形到三角形的极性变换的正交，以整合参考单线原点的奇异点。该正交是通过以下极坐标变换得到的。
 * \f[ \begin{pmatrix} x \\ y \end{pmatrix} = \begin{pmatrix} \frac{\hat
 * x}{\sin(\theta)+\cos(\theta)} cos(\theta) \\ \frac{\hat
 * x}{\sin(\theta)+\cos(\theta)} sin(\theta) \end{pmatrix} \qquad \theta
 * \dealcoloneq \frac\pi 2 \hat y \f]
 *
 *
 */
class QTrianglePolar : public QSimplex<2>
{
public:
  /**
   * 构建一个QTrianglePolar正交，在径向和角度方向的公式不同。
   * @param  radial_quadrature 径向正交  @param  angular_quadrature
   * 角度正交
   *
   */
  QTrianglePolar(const Quadrature<1> &radial_quadrature,
                 const Quadrature<1> &angular_quadrature);

  /**
   * 调用另一个构造函数，用QGauss<1>(n)表示径向和角向四分仪。
   * @param  n QGauss正交的阶数
   *
   */
  QTrianglePolar(const unsigned int n);
};

/**
 * 一个实现Duffy变换的正交，从一个正方形到一个三角形，以整合参考单轴的原点的奇异性。
 * 达菲变换被定义为\f[ \begin{pmatrix} x\\ y \end{pmatrix} =
 * \begin{pmatrix} \hat x^\beta (1-\hat y)\\ \hat x^\beta \hat y \end{pmatrix}
 * \f]。 Jacobian的行列式等于 $J= \beta \hat x^{2\beta-1}$
 * 。这种变换将参考正方形 $[0,1]\times[0,1]$
 * 映射到参考单轴，通过折叠正方形的左侧并将正交点向原点挤压，然后将得到的三角形剪切到参考三角形。当
 * $\beta = 1$ 在原点有阶 $1/R$
 * 的奇异点时，这种变换显示出良好的收敛特性，但当使用更高阶的高斯规则时，可以选择不同的
 * $\beta$ 值来提高收敛性和/或准确性（见 "Generalized Duffy
 * transformation for integrating vertex singularities", S. E. Mousavi, N.
 * Sukumar, Computational Mechanics 2009）。 当 $\beta = 1$
 * ，这种变换也被称为Lachat-Watson变换。
 *
 *
 */
class QDuffy : public QSimplex<2>
{
public:
  /**
   * 构造函数，允许沿 "径向 "和 "角度
   * "方向指定不同的正交规则。
   * 由于这种正交不是基于坐标的极地变化，所以谈论径向和角度方向是不完全恰当的。然而，达菲变换的效果类似于极坐标的改变，因为所产生的正交点是相对于奇点的径向排列。
   * @param  radial_quadrature 在径向使用的基础正交  @param
   * angular_quadrature 在角度使用的基础正交  @param  beta
   * 变换中使用的指数
   *
   */
  QDuffy(const Quadrature<1> &radial_quadrature,
         const Quadrature<1> &angular_quadrature,
         const double         beta = 1.0);

  /**
   * 用QGauss<1>(n)正交公式调用上述构造函数，以获得径向和角向的正交公式。
   * @param  n QGauss正交的阶数  @param  变换中使用的β指数
   *
   */
  QDuffy(const unsigned int n, const double beta);
};

/**
 * 当单元应被分割成子区域以使用一个或多个基础正交率进行积分时，使用的正交率。
 *
 *
 */
template <int dim>
class QSplit : public Quadrature<dim>
{
public:
  /**
   * 通过将参考超立方体分割成顶点零点与 @p split_point,
   * 重合的最小数量的单纯点，并将 @p base
   * 正交的仿生变换修补在一起，构建一个正交公式。点 @p
   * split_point
   * 应该在参考元素中，如果不是这样就会出现异常。
   * 在二维空间中，如果 @p split_point
   * 与其中一个顶点重合，如果它位于其中一条边上，或者如果它在参考元素的内部，那么产生的正交公式将分别由两个、三个或四个三角正交公式组成。
   * 三维情况也是如此，如果 @p split_point
   * 与其中一个顶点重合，如果它位于其中一条边上，位于其中一个面上，或者如果它是参考元素的内部，则分别有六个、八个、十个或十二个四面体正交公式。
   * 由此产生的正交可以用于，例如，在分裂点上具有可整定奇点的函数的积分，只要你选择一个可以在参考单纯线的顶点零上积分奇点的正交作为基础。
   * 一个维度为2的例子是这样的。
   * @code
   * const unsigned int order = 5;
   * QSplit<2> quad(QTrianglePolar(order), Point<2>(.3,.4));
   * @endcode
* 由此产生的正交函数将如下所示。    @image html split_quadrature.png ""   @param  base 要使用的基础QSimplex正交函数  @param  split_point 在哪里分割超立方体？
   *
   */
  QSplit(const QSimplex<dim> &base, const Point<dim> &split_point);
};

/**
 * 单元实体的整合规则。
 * 用户指定一个数字`n_points_1D`，作为准确集成的多项式程度的指示，类似于QGauss正交对象中的点数，尽管目前的正交公式不是张量积。对于n_points_1D=1,2,3,4，给定的值被转化为二维和三维的正交点的数量。
 *
 *
 *
 *
 *
 * - 2D: 1, 3, 7, 15
 *
 *
 *
 *
 * - 3D: 1, 4, 10, 35
 * 对于一维，正交规则退化为一个
 * `dealii::QGauss<1>(n_points_1D)`. 。
 *
 *
 * @ingroup simplex
 *
 *
 */
template <int dim>
class QGaussSimplex : public QSimplex<dim>
{
public:
  /**
   * 构造函数获取一维方向上的正交点数量  @p n_points_1D.  。
   *
   */
  explicit QGaussSimplex(const unsigned int n_points_1D);
};

/**
 * 单轴实体的Witherden-Vincent规则。
 * 像QGauss一样，用户应该指定一个数字`n_points_1D`，作为确切地整合什么程度的多项式的指示（例如，对于
 * $n$ 点，该规则可以确切地整合程度为 $2 n
 *
 * - 1$ 的多项式）。n_points_1D = 1, 2, 3, 4,
 * 5的给定值导致2D和3D的正交点数量如下。
 *
 *
 *
 * - 2D: 1, 6, 7, 15, 19
 *
 *
 *
 * - 3D: 1, 8, 14, 35, 59
 * 对于一维，正交规则退化为一个
 * `dealii::QGauss<1>(n_points_1D)`. 。 这些规则与quadpy  @cite quadpy
 * 库中为Witherden-Vincent列出的规则相匹配，并在  @cite
 * witherden2015identification  中首次描述。
 *
 *
 * @ingroup simplex
 *
 */
template <int dim>
class QWitherdenVincentSimplex : public QSimplex<dim>
{
public:
  /**
   * 构造函数，获取一维方向上的正交点数量  @p n_points_1D.
   * 。
   *
   */
  explicit QWitherdenVincentSimplex(const unsigned int n_points_1D);
};

/**
 * 楔形实体的积分规则。
 *
 *
 */
template <int dim>
class QGaussWedge : public Quadrature<dim>
{
public:
  /**
   * 用户指定一个数字`n_points_1D`作为精确集成的多项式程度的指示。详情请见QGaussSimplex的注释。
   *
   */
  explicit QGaussWedge(const unsigned int n_points_1D);
};

/**
 * 金字塔实体的积分规则。
 *
 *
 */
template <int dim>
class QGaussPyramid : public Quadrature<dim>
{
public:
  /**
   * 用户指定一个数字`n_points_1D`，作为准确集成的多项式程度的指示。详情请见QGaussSimplex的注释。
   *
   */
  explicit QGaussPyramid(const unsigned int n_points_1D);
};

 /*@}*/ 

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN
template <>
QGauss<1>::QGauss(const unsigned int n);
template <>
QGaussLobatto<1>::QGaussLobatto(const unsigned int n);

template <>
std::vector<double>
QGaussLog<1>::get_quadrature_points(const unsigned int);
template <>
std::vector<double>
QGaussLog<1>::get_quadrature_weights(const unsigned int);

template <>
QMidpoint<1>::QMidpoint();
template <>
QTrapezoid<1>::QTrapezoid();
template <>
QSimpson<1>::QSimpson();
template <>
QMilne<1>::QMilne();
template <>
QWeddle<1>::QWeddle();
template <>
QGaussLog<1>::QGaussLog(const unsigned int n, const bool revert);
template <>
QGaussLogR<1>::QGaussLogR(const unsigned int n,
                          const Point<1> &   x0,
                          const double       alpha,
                          const bool         flag);
template <>
QGaussOneOverR<2>::QGaussOneOverR(const unsigned int n,
                                  const unsigned int index,
                                  const bool         flag);
template <>
QTelles<1>::QTelles(const Quadrature<1> &base_quad,
                    const Point<1> &     singularity);
#endif // DOXYGEN



DEAL_II_NAMESPACE_CLOSE
#endif


