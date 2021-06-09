//include/deal.II-translator/grid/manifold_lib_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_manifold_lib_h
#define dealii_manifold_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/grid/manifold.h>

DEAL_II_NAMESPACE_OPEN

// forward declaration
namespace internal
{
  namespace MappingQGenericImplementation
  {
    template <int, int>
    class InverseQuadraticApproximation;
  }
} // namespace internal


/**
 * 极坐标系统的歧管描述。
 * 你可以使用这个Manifold对象来描述二维或三维的任何球体、圆、超球体或超盘，既可以作为共维一的流形描述符，也可以作为共维零的流形描述符，前提是南北两极（三维）和中心（二维和三维）不在Manifold中（因为它们是极坐标变化的奇点）。
 * 两个模板参数与Triangulation<dim,
 * spacedim>中的两个模板参数的含义一致，然而这个Manifold可以用来描述薄的和厚的物体，当dim
 * <=
 * spacedim时，行为是相同的，也就是说，PolarManifold<2,3>的功能与PolarManifold<3,3>相同。
 * 这个类的工作原理是将点转换为极坐标（包括二维和三维），在该坐标系中取平均值，然后再将点转换回笛卡尔坐标。为了使这个流形能够正常工作，它不能连接到包含坐标系中心或三维中的南北极的单元。这些点是坐标转换的奇异点，围绕这些点取平均值没有任何意义。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class PolarManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  /**
   * 构造函数获取球面坐标系的中心。
   * 该类使用pull_back和push_forward机制从笛卡尔坐标系转换到球面坐标系，考虑到二维中基数Manifold的周期性，而在三维中，它取中间点，并使用周围点的平均半径沿半径投影。
   *
   */
  PolarManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * 对这个Manifold对象进行克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 从欧几里得空间拉回给定的点。将返回与该点相关的极坐标
   * @p space_point.  只在spacedim = 2时使用。
   *
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * 给定球面坐标系中的一个点，该方法返回与极坐标相关的欧几里得坐标
   * @p chart_point.  只在spacedim = 3时使用。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override;

  /**
   * 给定spacedim维欧氏空间中的一个点，该方法返回从极坐标系映射到欧氏坐标系的函数
   * $F$ 的导数。换句话说，它是一个大小为
   * $\text{spacedim}\times\text{spacedim}$ 的矩阵。
   * 这个函数被用于get_tangent_vector()函数所要求的计算中。
   * 更多信息请参考该类的一般文档。
   *
   */
  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override;

  /**
   * @copydoc   Manifold::normal_vector() .
   *
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &p) const override;

  /**
   * 球形坐标系的中心。
   *
   */
  const Point<spacedim> center;

private:
  /**
   * 辅助函数，根据dim、chartdim和spacedim，返回与该坐标系相关的周期。
   *
   */
  static Tensor<1, spacedim>
  get_periodicity();
};



/**
 * 球形空间坐标系的流形描述。
 * 你可以使用这个流形对象来描述二维或三维的任何球体、圆、超球体或超盘。该流形可以作为嵌入高维空间的球面的同维度流形描述符，也可以作为具有正体积的体的同维度零流形描述符，前提是球面空间的中心不在该域中。使用这个函数的一个例子是在描述一个超壳或超球的几何形状时，例如在使用
 * GridGenerator::hyper_ball().
 * 创建一个粗略的网格后（然而，值得一提的是，为一个盘或球生成一个好的网格是复杂的，需要增加步骤。参见
 * step-6 中的 "扩展的可能性
 * "一节，对如何构建这样的网格以及需要做什么进行了广泛的讨论）。)
 * 两个模板参数与Triangulation<dim,
 * spacedim>中的两个模板参数的含义一致，然而这个Manifold可以用来描述薄的和厚的对象，当dim
 * <=
 * spacedim时，行为是相同的，也就是说，SphericalManifold<2,3>的功能与SphericalManifold<3,3>相同。
 * 虽然PolarManifold反映了通常的极坐标概念，但它可能不适合于包含南北两极的域。
 * 例如，考虑极坐标中的一对点 $x_1=(1,\pi/3,0)$ 和
 * $x_2=(1,\pi/3,\pi)$ （位于半径为1的球体表面，高度为 $\pi/3$
 * 的平行线上）。在这种情况下，用极地坐标的直线连接这两个点，将需要绕地球一圈的长路，而不经过北极。
 * 这两个点将通过曲线连接（使用PolarManifold）。
 *
 * @f{align*}{
 * s: [0,1]  & \rightarrow &  \mathbb S^3 \\
 *         t & \mapsto     &  (1,\pi/3,0) + (0,0,t\pi)
 * @f}
 * 这条曲线不是球体上的测地线，也不是我们连接这两点的方式。更好的曲线是通过北极的那条。@f[
 * s(t) = x_1 \cos(\alpha(t)) + \kappa \times x_1 \sin(\alpha(t)) +
 * \kappa ( \kappa \cdot x_1) (1-\cos(\alpha(t))).
 * @f] 其中 $\kappa = \frac{x_1 \times x_2}{\Vert x_1 \times x_2 \Vert}$ 和 $\alpha(t) = t \cdot \arccos(x_1 \cdot x_2)$ 为 $t\in[0,1]$  。事实上，这是一个测地线，在连接球面上的点时，它是自然的选择。在上面的例子中，PolarManifold类实现了连接球体表面两点的第一种方式，而SphericalManifold则实现了第二种方式，也就是说，这个Manifold使用测地线连接点。如果通过 SphericalManifold::get_new_points() 调用涉及两个以上的点，则使用所谓的球形平均，其中最终的点通过测地线最小化到所有其他点的加权距离。
 * 特别是，这个类实现了一个Manifold，它连接空间中的任何两点，首先将它们投射到具有单位半径的球面上，然后用一个测地线将它们连接起来，最后重新调整最终半径，使之成为起始半径的加权平均值。这个Manifold与二维的PolarManifold相同，而对于三维，它返回的点在球面上的分布更加均匀，而且它对于坐标系的旋转是不变的，因此避免了PolarManifold在两极出现的问题。请注意，特别是用PolarManifold计算两极的切向量是不好定义的，而用这个类是完全可以的。
 * 由于数学上的原因，不可能只用测地曲线来构造球体的唯一地图，因此不鼓励将这个类与MappingManifold一起使用。如果你使用此
 * Manifold 来描述球体的几何结构，你应该使用 MappingQ
 * 作为底层映射，而不是 MappingManifold。
 * 这个Manifold只能用于*从中心移出一个有限半径的球的几何体。事实上，中心是这个流形的一个奇异点，如果你试图将两个点穿过中心连接起来，它们会在球面坐标上移动，避开中心。
 * 这个流形的理想几何结构是一个超壳。如果你打算在一个HyperBall上使用这个Manifold，你必须确保不把这个Manifold连接到包含中心的单元。建议将该类与TransfiniteInterpolationManifold相结合，以确保从曲线形状平滑过渡到球中心的直线坐标系。(也可参见
 * step-65  中的广泛讨论) 。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class SphericalManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * 构造函数获取球面坐标的中心。
   *
   */
  SphericalManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * 制作此Manifold对象的一个克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 给出空间中的任意两点，首先将它们投影在具有单位半径的球面上，然后用测地线连接它们并找到中间点，最后重新调整最终半径，使得到的半径是起始半径的凸组合。
   *
   */
  virtual Point<spacedim>
  get_intermediate_point(const Point<spacedim> &p1,
                         const Point<spacedim> &p2,
                         const double           w) const override;

  /**
   * 计算参数w等于零的get_intermediate_point()函数的导数。
   *
   */
  virtual Tensor<1, spacedim>
  get_tangent_vector(const Point<spacedim> &x1,
                     const Point<spacedim> &x2) const override;

  /**
   * @copydoc   Manifold::normal_vector() .
   *
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &p) const override;

  /**
   * 计算每个顶点到边界的法向量。
   *
   */
  virtual void
  get_normals_at_vertices(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
    const override;

  /**
   * 计算一组新的点，在给定的点之间进行插值  @p
   * surrounding_points。  @p weights 是一个表格，其列数与 @p
   * surrounding_points.size()相同。 @p weights 中的行数必须与 @p
   * new_points.
   * 的长度相匹配。这个函数被优化为在新点的集合上执行，通过在所有新点的循环之外收集不依赖于权重的操作。
   * 该实现不允许 @p surrounding_points 和 @p new_points
   * 指向同一个数组，所以要确保将不同的对象传入该函数。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double> &                weights,
                 ArrayView<Point<spacedim>> new_points) const override;

  /**
   * 返回球面流形上的一个点，该点相对于周围的点来说是中间点。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &vertices,
                const ArrayView<const double> &         weights) const override;

  /**
   * 球面坐标系的中心。
   *
   */
  const Point<spacedim> center;

private:
  /**
   * 返回球面流形上的一个点，该点相对于周围的点而言是中间点。这个函数使用方向的线性平均值来寻找一个估计的点。它返回一对从中心点到候选点的半径和方向。
   *
   */
  std::pair<double, Tensor<1, spacedim>>
  guess_new_point(const ArrayView<const Tensor<1, spacedim>> &directions,
                  const ArrayView<const double> &             distances,
                  const ArrayView<const double> &             weights) const;

  /**
   * 返回球面流形上的一个点，该点相对于周围的点来说是中间的。这个函数使用一个候选点作为猜测，并执行牛顿式迭代来计算正确的点。
   * 实现的主要部分使用了出版物Buss, Samuel R.和Jay P.
   * Fillmore中的观点。
   * "球面平均数和球面样条和插值的应用"。ACM Transactions on
   * Graphics (TOG) 20.2 (2001):
   * 95-126.特别是在http://math.ucsd.edu/~sbuss/ResearchWeb/spheremean/提供的实现。
   *
   */
  Point<spacedim>
  get_new_point(const ArrayView<const Tensor<1, spacedim>> &directions,
                const ArrayView<const double> &             distances,
                const ArrayView<const double> &             weights,
                const Point<spacedim> &candidate_point) const;

  /**
   * 计算一组新的点，在给定的点之间进行插值  @p
   * surrounding_points。  @p weights
   * 是一个数组视图，其条目数为 @p
   * surrounding_points.size()乘以 @p new_points.size().
   * 这个函数被优化为对新点的集合执行，通过收集所有新点上的循环之外的不依赖于权重的操作。
   * 该实现不允许 @p surrounding_points 和 @p new_points
   * 指向同一个数组，所以要确保将不同的对象传入该函数。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const ArrayView<const double> &         weights,
                 ArrayView<Point<spacedim>>              new_points) const;

  /**
   * 一个流形描述，用于二维的get_new_point。
   *
   */
  const PolarManifold<spacedim> polar_manifold;
};

/**
 * 圆柱形流形描述。
 * 在三维空间中，点使用圆柱坐标系沿<tt>x-</tt>,
 * <tt>y-</tt>或<tt>z</tt>轴进行转换（当使用该类的第一个构造函数时），或者使用一个由其轴线方向和位于轴线上的点描述的任意方向的圆柱。
 * 这个类的开发是为了与GridGenerator的 @p cylinder 或 @p
 * cylinder_shell
 * 函数结合使用。只要spacedim不等于3，这个函数就会抛出一个运行时异常。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class CylindricalManifold : public ChartManifold<dim, spacedim, 3>
{
public:
  /**
   * 构造函数。使用构造函数参数的默认值可以得到一个沿x轴的圆柱体（<tt>axis=0</tt>）。选择<tt>axis=1</tt>或<tt>axis=2</tt>，分别得到一个沿y轴或z轴的管。公差值是用来确定一个点是否在轴上的。
   *
   */
  CylindricalManifold(const unsigned int axis      = 0,
                      const double       tolerance = 1e-10);

  /**
   * 构造函数。如果用这个构造函数构造，描述的流形是一个圆柱体，其轴线指向#方向，并穿过给定的#点_在轴上。方向可以任意缩放，而给定的点可以是轴上的任何一点。公差值用于确定一个点是否在轴上。
   *
   */
  CylindricalManifold(const Tensor<1, spacedim> &direction,
                      const Point<spacedim> &    point_on_axis,
                      const double               tolerance = 1e-10);

  /**
   * 对这个Manifold对象进行克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 计算给定空间点的圆柱坐标 $(r, \phi, \lambda)$ ，其中 $r$
   * 表示与轴线的距离， $\phi$
   * 给定点与计算出的法线方向的角度，以及 $\lambda$
   * 轴向位置。
   *
   */
  virtual Point<3>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * 计算以圆柱坐标 $(r, \phi, \lambda)$
   * 给出的图表点的直角坐标，其中 $r$
   * 表示与轴线的距离， $\phi$
   * 给出的点与计算的法线方向之间的角度，以及 $\lambda$
   * 的轴向位置。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<3> &chart_point) const override;

  /**
   * 计算从圆柱坐标 $(r, \phi, \lambda)$
   * 到车轴坐标的映射的导数，其中 $r$ 表示与轴的距离，
   * $\phi$ 表示给定点与计算的法线方向之间的角度，以及
   * $\lambda$ 表示轴向位置。
   *
   */
  virtual DerivativeForm<1, 3, spacedim>
  push_forward_gradient(const Point<3> &chart_point) const override;

  /**
   * 计算CylindricalManifold上的新点。关于这个函数的详细描述，请参见基类的文档。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double> &         weights) const override;

protected:
  /**
   * 一个与法线方向正交的矢量。
   *
   */
  const Tensor<1, spacedim> normal_direction;

  /**
   * 轴的方向向量。
   *
   */
  const Tensor<1, spacedim> direction;

  /**
   * 轴上的一个任意点。
   *
   */
  const Point<spacedim> point_on_axis;

private:
  /**
   * 用于测量零距离的相对公差。
   *
   */
  double tolerance;
};

/**
 * 椭圆流形描述，源自ChartManifold。关于椭圆坐标系的更多信息可以在<a
 * href="https://en.wikipedia.org/wiki/Elliptic_coordinate_system">Wikipedia
 * </a>找到。
 * 这是基于椭圆坐标的定义 $(u,v)$ @f[
 * \left\lbrace\begin{aligned}
 * x &=  x_0 + c \cosh(u) \cos(v) \\
 * y &=  y_0 + c \sinh(u) \sin(v)
 * \end{aligned}\right.
 * @f]，其中 $(x_0,y_0)$ 是笛卡尔系统的中心坐标。
 * 目前的实现使用坐标 $(c,v)$ ，而不是 $(u,v)$
 * ，并根据给定的偏心率来固定 $u$
 * 。因此，这种坐标的选择产生了一个椭圆流形，其特点是偏心率不变。
 * $e=\frac{1}{\cosh(u)}$  ，其中 $e\in\left]0,1\right[$  。
 * 如果dim和spacedim都不同于2，这个类的构造函数将抛出一个异常。
* 这个流形可以用来产生具有椭圆曲率的超壳。作为一个例子，测试<B>elliptical_manifold_01</B>产生了以下三角结构。  @image html elliptical_hyper_shell.png
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class EllipticalManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  /**
   * 构造函数，接收流形系统的中心、主轴方向和流形偏心率。
   * 默认的主轴是<tt>x</tt>-轴。流形系统会被旋转，以使主轴与输入中指定的方向一致。
   * @param  中心 分流板的中心。    @param  major_axis_direction
   * 分流板的主轴方向。    @param  偏心率 分流板的偏心率
   * $e\in\left]0,1\right[$  。
   *
   */
  EllipticalManifold(const Point<spacedim> &    center,
                     const Tensor<1, spacedim> &major_axis_direction,
                     const double               eccentricity);

  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * @copydoc   ChartManifold::pull_back()  。
   *
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const override;

  /**
   * @copydoc   ChartManifold::push_forward()  。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const override;

  /**
   * @copydoc   ChartManifold::push_forward_gradient()
   *
   */
  virtual DerivativeForm<1, spacedim, spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const override;


protected:
  /**
   * 主轴的方向矢量。
   *
   */
  Tensor<1, spacedim> direction;
  /**
   * 流形的中心。
   *
   */
  const Point<spacedim> center;
  /**
   * 从流形体的偏心率得出的参数。
   *
   */
  const double cosh_u;
  const double sinh_u;

private:
  /**
   * @copydoc   ChartManifold::get_periodicity()  对于 $\text{dim}=2$ 和
   * $\text{spacedim}=2$
   * ，第一个坐标是非周期性的，而第二个坐标的周期性为
   * $2\pi$  。
   *
   */
  static Tensor<1, spacedim>
  get_periodicity();
};


/**
 * 从ChartManifold派生的Manifold描述，基于明确的Function<spacedim>和Function<chartdim>对象，描述push_forward（）和pull_back（）函数。
 * 你可以使用这个Manifold对象来描述任何任意的形状域，只要你能用可逆映射来表达，你可以提供正向表达式和逆向表达式。
 * 在调试模式下，会进行一个检查，以验证这些变换实际上是一个反的变换。
 *
 *
 * @ingroup manifold
 *
 */
template <int dim, int spacedim = dim, int chartdim = dim>
class FunctionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:
  /**
   * 明确的函数构造器。接受一个spacedim组件的push_forward函数，以及一个
   * @p chartdim 组件的pull_back函数。关于可选的 @p periodicity
   * 参数的含义，请参阅基类ChartManifold的文档。
   * 容许参数在调试模式下被用来实际检查这两个函数是否是一个反义词。
   * 注意：以这种方式构造的对象存储了指向push_forward和pull_back函数的指针。因此，必须保证这些函数对象只在构造的流形之后被销毁。
   *
   */
  FunctionManifold(
    const Function<chartdim> & push_forward_function,
    const Function<spacedim> & pull_back_function,
    const Tensor<1, chartdim> &periodicity = Tensor<1, chartdim>(),
    const double               tolerance   = 1e-10);

  /**
   * 和前面一样，只是这个构造函数对作为第一和第二个参数传递的Function对象拥有所有权，并在FunctionManifold对象被销毁时最终负责删除这些指针。
   * 这个构造函数很有用，因为它允许在调用构造函数的地方创建函数对象，而不需要命名和随后删除这些对象。这允许以下的成语。
   * FunctionManifold<dim>  manifold(std::make_unique<MyPushForward>(...),
   * std::make_unique<MyPullBack>(...)); 。
   *
   */
  FunctionManifold(
    std::unique_ptr<Function<chartdim>> push_forward,
    std::unique_ptr<Function<spacedim>> pull_back,
    const Tensor<1, chartdim> &         periodicity = Tensor<1, chartdim>(),
    const double                        tolerance   = 1e-10);

  /**
   * 表达式构造器。接受spacedim组件的push_forward函数和 @p
   * chartdim组件的pull_back函数的表达。关于可选的 @p
   * periodicity 参数的含义，请参阅基类ChartManifold的文档。
   * 字符串应该是FunctionParser类的默认构造函数所能读取的。你可以用最后两个可选参数指定自定义的变量表达式。如果你不这样做，就会使用默认的名称，即
   * "x,y,z"。
   * 容许参数在调试模式下用于实际检查这两个函数是否是一个反义词。
   *
   */
  FunctionManifold(
    const std::string          push_forward_expression,
    const std::string          pull_back_expression,
    const Tensor<1, chartdim> &periodicity = Tensor<1, chartdim>(),
    const typename FunctionParser<spacedim>::ConstMap =
      typename FunctionParser<spacedim>::ConstMap(),
    const std::string chart_vars =
      FunctionParser<chartdim>::default_variable_names(),
    const std::string space_vars =
      FunctionParser<spacedim>::default_variable_names(),
    const double tolerance = 1e-10,
    const double h         = 1e-8);

  /**
   * 如果需要，我们会删除我们拥有的指针。
   *
   */
  virtual ~FunctionManifold() override;

  /**
   * 对这个Manifold对象做一个克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 给定 @p chartdim
   * 坐标系中的一个点，使用push_forward_函数来计算 @p
   * chartdim空间维度中的点到 @p spacedim
   * 空间维度的push_forward。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const override;

  /**
   * 给定chartdim维度欧几里得空间中的一个点，该方法返回从sub_manifold坐标系映射到欧几里得坐标系的函数
   * $F$ 的导数。换句话说，它是一个大小为
   * $\text{spacedim}\times\text{chartdim}$ 的矩阵。
   * 这个函数被用于get_tangent_vector()函数所要求的计算中。默认实现是调用
   * FunctionManifold::push_forward_function()
   * 成员类的get_gradient()方法。如果你使用接受两个字符串表达式的构造函数来构造这个对象，那么这个方法的默认实现使用有限差分方案来计算梯度（详见AutoDerivativeFunction()类），你可以在构造时用
   * @p h 参数指定空间步长的大小。
   * 更多信息请参考该类的一般文档。
   *
   */
  virtual DerivativeForm<1, chartdim, spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const override;

  /**
   * 给定spacedim坐标系中的一个点，使用pull_back_函数来计算
   * @p spacedim 空间维度中的点对 @p chartdim 空间维度的回撤。
   *
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const override;

private:
  /**
   * FunctionParser类的常量。
   *
   */
  const typename FunctionParser<spacedim>::ConstMap const_map;

  /**
   * 指向push_forward函数的指针。
   *
   */
  SmartPointer<const Function<chartdim>,
               FunctionManifold<dim, spacedim, chartdim>>
    push_forward_function;

  /**
   * 指向pull_back函数的指针。
   *
   */
  SmartPointer<const Function<spacedim>,
               FunctionManifold<dim, spacedim, chartdim>>
    pull_back_function;

  /**
   * 相对公差。在调试模式下，我们检查在构造时提供的两个函数实际上是另一个函数的逆向。
   * 在这个检查中，这个值被用来作为相对公差。
   *
   */
  const double tolerance;

  /**
   * 检查智能指针的所有权。表示这个类是否是前两个成员变量所指向的对象的所有者。
   * 这个值是在类的构造函数中设置的。如果 @p true,
   * ，那么析构器将删除这两个指针指向的函数对象。
   *
   */
  bool owns_pointers;

  /**
   * 用来构造push_forward函数的表达式。
   *
   */
  const std::string push_forward_expression;

  /**
   * 用于构造pull_back函数的表达式。
   *
   */
  const std::string pull_back_expression;

  /**
   * 图表域中的变量名称。
   *
   */
  const std::string chart_vars;

  /**
   * 空间域中的变量名称。
   *
   */
  const std::string space_vars;

  /**
   * 内部使用的有限差分步骤。
   *
   */
  const double finite_difference_step;
};



/**
 * 三维环状体表面的流形描述。环状体被假定为在x-z平面内。参考坐标系由围绕Y轴的角度
 * $phi$ 、围绕环状体中心线的角度 $theta$ 和到中心线的距离
 * $w$ （在0和1之间）给出。 这个类的开发是为了与
 * GridGenerator::torus. 一起使用。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim>
class TorusManifold : public ChartManifold<dim, 3, 3>
{
public:
  static const int chartdim = 3;
  static const int spacedim = 3;

  /**
   * 构造函数。指定中心线的半径 @p R 和环本身的半径(  @p
   * r).  变量的含义与 GridGenerator::torus(). 中的参数相同。
   *
   */
  TorusManifold(const double R, const double r);

  /**
   * 对这个Manifold对象做一个克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, 3>>
  clone() const override;

  /**
   * 拉回操作。
   *
   */
  virtual Point<3>
  pull_back(const Point<3> &p) const override;

  /**
   * 前推操作。
   *
   */
  virtual Point<3>
  push_forward(const Point<3> &chart_point) const override;

  /**
   * 梯度。
   *
   */
  virtual DerivativeForm<1, 3, 3>
  push_forward_gradient(const Point<3> &chart_point) const override;

private:
  double r, R;
};



/**
 * 一个映射类，它将曲线边界描述扩展到计算域的内部。外侧的弯曲边界描述被假定为由另一个流形（例如圆上的极地流形）给出。扩展边界信息的机制是一个所谓的转折性插值。该类方法的使用在
 * step-65 中得到了广泛的讨论。
 * 在二维中扩展这种描述的公式，例如，在<a
 * href="https://en.wikipedia.org/wiki/Transfinite_interpolation">
 * Wikipedia</a>上有描述。 给定图表上的一个点 $(u,v)$
 * ，这个点在实空间中的图像由以下公式给出
 *
 * @f{align*}{
 * \mathbf S(u,v) &= (1-v)\mathbf c_0(u)+v \mathbf c_1(u) + (1-u)\mathbf c_2(v)
 * + u \mathbf c_3(v) \\
 * &\quad
 *
 * - \left[(1-u)(1-v) \mathbf x_0 + u(1-v) \mathbf x_1 + (1-u)v \mathbf
 * x_2 + uv \mathbf x_3 \right]
 * @f}
 * 其中 $\bf x_0, \bf x_1, \bf x_2, \bf x_3$
 * 表示限定图像空间的四个边界顶点， $\bf c_0, \bf c_1, \bf
 * c_2, \bf c_3$
 * 是描述单元格线条的四条曲线。如果弯曲的流形连接到这些线条中的任何一条，则根据
 * Manifold::get_new_point()
 * 的规定，用线条的两个端点和适当的权重进行评估。在三维中，这个公式的一般化被实施，创建一个顶点（正贡献）、线（负贡献）和面（正贡献）的加权和。
 * 这个流形通常被附加到一个粗略的网格上，然后将新的点作为边界上的描述的组合，根据点在原图坐标中的位置进行适当的加权
 * $(u,v)$
 * 。在大多数情况下，这种流形应该比只在网格的边界上设置一个弯曲的流形要好，因为随着网格的细化，它可以产生更均匀的网格分布，因为它在这个流形所连接的初始粗单元的所有子节点上从弯曲的描述转换为直线描述。这样一来，一旦网格被细化，原本包含在一个<i>coarse</i>网格层中的流形的弯曲性质将被应用到多个<i>fine</i>网格层。请注意，当只有一个单元的表面受到曲面描述时，TransfiniteInterpolationManifold的机制也被内置于MappingQGeneric类中，确保在应用曲面边界描述时，即使没有这个流形的默认情况也能获得最佳收敛率。
 * 如果没有弯曲的边界围绕着一个粗大的单元，这个类就会还原成一个平面流形描述。
 * 为了举一个使用这个类的例子，下面的代码将一个转折流形附加到一个圆上。
 *
 *
 * @code
 * PolarManifold<dim> polar_manifold;
 * TransfiniteInterpolationManifold<dim> inner_manifold;
 *
 * Triangulation<dim> triangulation;
 * GridGenerator::hyper_ball (triangulation);
 *
 * triangulation.set_all_manifold_ids(1);
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, polar_manifold);
 * inner_manifold.initialize(triangulation);
 * triangulation.set_manifold (1, inner_manifold);
 * triangulation.refine_global(4);
 * @endcode
 *
 * 在这段代码中，我们首先将所有流形的id设置为转折性插值的id，然后重新设置边界上的流形id，以确定极地流形所描述的弯曲边界。使用这段代码，我们可以得到一个非常漂亮的网格。
 * <p ALIGN="center">
 @image html circular_mesh_transfinite_interpolation.png
 * </p> 这显然比只应用于边界的极地流形要好得多。 <p
 * ALIGN="center">
 @image html circular_mesh_only_boundary_manifold.png
 * </p> 这个流形被用于一些GridGenerator函数中，包括
 * GridGenerator::channel_with_cylinder. 。 <h3>Implementation details</h3>
 * 在这个类的实现中，围绕一个粗略单元的流形被反复查询，以计算其内部的点。为了获得最佳的网格质量，这些流形应该与一个图表概念兼容。例如，使用两个顶点的权重0.25和0.75来计算两个顶点之间沿线的0.25的点，应该得到与首先计算0.5的中点，然后再次计算第一个顶点和粗略的中点相同的结果。deal.II提供的大多数流形类都是如此，如SphericalManifold或PolarManifold，但天真的实现可能会违反这一规定。如果流形的质量不够好，在网格细化时，可能会发生get_new_point()或get_new_points()方法中的图表转换产生的点在单元格外。那么这个类就会抛出一个类型为
 * Mapping::ExcTransformationFailed.
 * 的异常。在这种情况下，应该在附加这个类之前对网格进行细化，就像下面的例子那样。
 *
 *
 * @code
 * SphericalManifold<dim> spherical_manifold;
 * TransfiniteInterpolationManifold<dim> inner_manifold;
 * Triangulation<dim> triangulation;
 * GridGenerator::hyper_ball (triangulation);
 *
 * triangulation.set_all_manifold_ids(1);
 * triangulation.set_all_manifold_ids_on_boundary(0);
 * triangulation.set_manifold (0, spherical_manifold);
 * inner_manifold.initialize(triangulation);
 * triangulation.set_manifold (1, inner_manifold);
 * triangulation.refine_global(1);
 *
 * // initialize the transfinite manifold again
 * inner_manifold.initialize(triangulation);
 * triangulation.refine_global(4);
 * @endcode
 *
 *
 *
 * @note
 * 出于性能和精度的考虑，建议将转义流形应用于尽可能粗的网格。关于精度，曲面描述只能应用于从给定邻域创建的新点，在尽可能大的域上扩展曲面描述时，网格质量通常会更高。关于性能，get_new_point()方法中正确的粗单元的识别需要通过所有的粗单元，因此预计每个单一的映射操作的粗单元数量的线性复杂性，也就是说，对于整个网格的任何全局操作，至少是粗网格单元数量的二次方。因此，目前的实现只有在粗单元数量不超过几百个的情况下才是经济的。为了使更多单元的性能得到提高，我们可以通过预先识别具有轴对齐边界盒的相关单元来扩展当前的实现。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class TransfiniteInterpolationManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * 构建器。
   *
   */
  TransfiniteInterpolationManifold();

  /**
   * 解构器。
   *
   */
  virtual ~TransfiniteInterpolationManifold() override;

  /**
   * 对这个Manifold对象进行克隆。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 用一个粗略的网格来初始化流形。使用该类的前提条件是，输入的三角形是均匀细化的，并且流形后来被附加到同一三角形上。
   * 每当流形ID的分配在初始化该类的三角形上发生变化时，必须再次调用initialize()，以更新连接到粗大单元的流形ID。
   * @note
   * 在使用此对象的过程中，不得销毁用于构建流形的三角图。
   *
   */
  void
  initialize(const Triangulation<dim, spacedim> &triangulation);

  /**
   * 返回将成为新顶点的点，该点被给定的点所包围  @p
   * surrounding_points.   @p weights
   * 包含周围点的适当权重，流形根据该权重决定新点的位置。
   * 这个类的实现覆盖了基类中的方法，并通过转折插值计算新点。实现的第一步是确定周围点所处的粗略单元。然后，通过牛顿迭代将坐标转换为粗单元上的单位坐标，然后根据权重计算出新的点。最后，根据无限插值将其向前推至实空间。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double> &         weights) const override;

  /**
   * 计算一组新的点，在给定的点之间进行插值  @p
   * surrounding_points。  @p weights 是一个表，其列数与 @p
   * surrounding_points.size()相同。 @p weights 中的列数必须与 @p
   * new_points.
   * 的长度相匹配。这个类中的实现覆盖了基类中的方法，并通过一个转折性的插值计算新点。实现的第一步是确定周围点所处的粗略单元。然后，通过牛顿迭代将坐标转换为粗单元上的单位坐标，然后根据权重计算出新的点。最后，根据转折内插法将其向前推到实空间。
   * 实现不允许 @p surrounding_points 和 @p new_points
   * 指向同一个向量，所以要确保将不同的对象传入该函数。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double> &                weights,
                 ArrayView<Point<spacedim>> new_points) const override;

private:
  /**
   * 内部函数，用于识别给定的周围点所在的最合适的单元（=图表）。我们使用一种廉价的算法来识别单元格，并在实际进行相关单元格内的搜索之前按概率对单元格进行排序。这些单元是按照逆映射的Q1近似值与周围点的单元格的距离来排序的。我们期望最多有20个单元（在三维结构的网格上最多有8个候选单元，在非结构的网格上会更多一些，通常我们只得到两到三个），所以得到一个有20个条目的数组，其索引为<tt>cell->index()</tt>。
   *
   */
  std::array<unsigned int, 20>
  get_possible_cells_around_points(
    const ArrayView<const Point<spacedim>> &surrounding_points) const;

  /**
   * 最终确定正确的图表，并用周围点的回撤来填充 @p
   * chart_points。这个方法在内部调用 @p
   * get_possible_cells_around_points().
   * 返回一个迭代器到定义了图表的单元格。
   *
   */
  typename Triangulation<dim, spacedim>::cell_iterator
  compute_chart_points(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    ArrayView<Point<dim>>                   chart_points) const;

  /**
   * 在给定的粗略单元上对单位坐标进行回拉操作。
   * 这个方法目前是基于一个类似牛顿的迭代来寻找原点的。我们可以通过提供一个好的初始猜测作为第三个参数来加快迭代速度。如果没有更好的点，可以使用cell->real_to_unit_cell_affine_approximation(p)
   * @note  这个内部函数目前与 ChartManifold::pull_back()
   * 函数不兼容，因为给定的类代表一个图表图集，而不是一个单一的图表。因此，pull_back()操作只对图表的附加信息有效，这些信息由粗略网格上的单元格给出。另一种实现方式可以根据粗网格的单元来转移索引，在图表空间和图像空间之间形成1对1的关系。
   *
   */
  Point<dim>
  pull_back(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
            const Point<spacedim> &                                     p,
            const Point<dim> &initial_guess) const;

  /**
   * 前推操作。
   * @note  这个内部函数目前与 ChartManifold::push_forward()
   * 函数不兼容，因为给定的类代表一个图表图集，而不是一个单一的图表。因此，push_forward()操作只对图表的额外信息有效，这些信息由粗略网格上的单元格给出。另一种实现方式可以根据粗网格的单元来转移索引，在图表空间和图像空间之间形成1比1的关系。
   *
   */
  Point<spacedim>
  push_forward(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
               const Point<dim> &chart_point) const;

  /**
   * push_forward方法的梯度。
   * @note 这个内部函数与 ChartManifold::push_forward_gradient()
   * 函数不兼容，因为给定的类代表一个图表图集，而不是一个单一的图表。此外，这个私有函数还要求用户为这个函数的单一使用情况提供对图表点的push_forward()调用的结果，即在牛顿迭代内部，梯度是通过有限差分计算的。
   *
   */
  DerivativeForm<1, dim, spacedim>
  push_forward_gradient(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &                                          chart_point,
    const Point<spacedim> &pushed_forward_chart_point) const;

  /**
   * 底层的三角结构。
   *
   */
  const Triangulation<dim, spacedim> *triangulation;

  /**
   * 应用转折性近似的网格单元的级别，通常为0级。
   *
   */
  int level_coarse;

  /**
   * 如果周围的流形都是转折流形或者有默认的（无效的）流形ID，流形就会退化为平流形，我们可以为push_forward方法选择便宜的算法。
   *
   */
  std::vector<bool> coarse_cell_is_flat;

  /**
   * 用于计算图表空间中的新点的平坦流形，我们使用FlatManifold描述。
   *
   */
  FlatManifold<dim> chart_manifold;

  /**
   * 对每个粗网格单元来说，从实数点到图表点的逆向映射的四次方近似的向量。
   *
   */
  std::vector<internal::MappingQGenericImplementation::
                InverseQuadraticApproximation<dim, spacedim>>
    quadratic_approximation;

  /**
   * 与 Triangulation::signals::clear
   * 的连接，一旦这个类出了范围，必须重新设置。
   *
   */
  boost::signals2::connection clear_signal;
};

DEAL_II_NAMESPACE_CLOSE

#endif


