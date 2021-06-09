//include/deal.II-translator/grid/manifold_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_tria_manifold_h
#define dealii_tria_manifold_h


 /*----------------------------   manifold.h     ---------------------------*/ 

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

// forward declaration
#ifndef DOXYGEN
template <int, typename>
class Table;
#endif

/**
 * 我们在这里收集一些Manifold<dim,spacedim>类中使用的辅助函数。
 *
 *
 */
namespace Manifolds
{
  /**
   * 一个 <code>constexpr</code> 辅助函数，用于返回给定的
   * <code>MeshIteratorType</code>
   * 所指向的结构类型的默认点的数量。更多信息请参见
   * Manifolds::get_default_points_and_weights()  的文档。
   *
   */
  template <typename MeshIteratorType>
  inline constexpr std::size_t
  n_default_points_per_cell()
  {
    // Note that in C++11 a constexpr function can only have a return
    // statement, so we cannot alias the structure dimension
    return GeometryInfo<MeshIteratorType::AccessorType::structure_dimension>::
             vertices_per_cell +
           GeometryInfo<MeshIteratorType::AccessorType::structure_dimension>::
             lines_per_cell +
           GeometryInfo<MeshIteratorType::AccessorType::structure_dimension>::
             quads_per_cell +
           GeometryInfo<MeshIteratorType::AccessorType::structure_dimension>::
             hexes_per_cell -
           1; // don't count the cell itself, just the bounding objects
  }

  /**
   * 给定一个一般的网格迭代器，构建包含以下点数的正交点和权重数组。
   *
   *
   *
   *
   *
   *
   * - 如果迭代器指向一条直线，那么正交点就是该直线的两个顶点。这就产生了一个有两个点的点向量。
   *
   *
   *
   *
   *
   *
   * - 如果迭代器指向一个四边形，那么四边形的点就是顶点和直线中点。这将导致一个有八个（4+4）点的点向量。
   *
   *
   *
   *
   *
   *
   * - 如果迭代器指向一个六面体，那么正交点就是顶点、线中点和面中点。    这将产生一个有26（8+12+6）个点的向量。    这些点的正交权重的选择是相同的，并且等于正交点数量的1（如果 @p  with_interpolation是 @p false), ），或者以一种方式给更接近单元中心（在参考单元上测量）的点更高的权重。这些权重对应于无限插值中应用于线和顶点的权重（见 TransfiniteInterpolationManifold 的更详尽描述），如果  @p with_interpolation  是  @p true.  该函数主要用于构建  Manifold::get_new_point()  函数的输入参数，该函数基于由正交点和存储在返回的一对向量中的权重代表的 "周围 "点的加权平均值计算流形上的一个新点。这个函数根据 "围绕 "一个单元、面或边缘的点来创建这样一个对象，而权重的选择则适合于计算指向的对象的新 "中点"。一个需要这样做的例子是网格细化，（以2D情况为例）我们需要首先创建新的边中点，然后是新的单元格点。      @param[in]  iterator 一个网格迭代器，指向直线、四边形或六边形。    @param[in]  with_interpolation 是否从转折插值计算正交权重，如上所述。    @tparam  MeshIteratorType 对应于 Triangulation::cell_iterator （或诸如 Triangulation::active_cell_iterator 或 DoFHandler::cell_iterator) 等变体的迭代器类型，或者是诸如 <code>cell-@>face(f)</code> or <code>cell-@>line(l)</code> 等语句的结果.
   *
   */
  template <typename MeshIteratorType>
  std::pair<std::array<Point<MeshIteratorType::AccessorType::space_dimension>,
                       n_default_points_per_cell<MeshIteratorType>()>,
            std::array<double, n_default_points_per_cell<MeshIteratorType>()>>
  get_default_points_and_weights(const MeshIteratorType &iterator,
                                 const bool with_interpolation = false);
} // namespace Manifolds



/**
 * 漫游是用来描述域的边界几何以及内部几何的。因此，流形对象与单元格、面和/或边相关联，可以由用户直接操作，或者，如果用户程序没有明确这样做，则使用默认流形对象。
 * 歧管最好用微分几何的语言来理解，但它们的常见用途可以通过例子简单描述。在 @ref geometry_paper "几何学论文 "
 * 中提供了关于如何、在何处以及为何使用该类的详尽讨论。
 *
 *  <h3>Common use case: Creating a new vertex</h3>
 * 在流形最本质的使用中，流形描述被用来创建一个
 * "其他点之间的点"。例如，当一个三角形在单元格、面或边上创建一个新的顶点时，它通过以下函数调用确定新顶点的坐标。
 * @code
 *   ...
 *   Point<spacedim> new_vertex = manifold.get_new_point (points,weights);
 *   ...
 * @endcode
 * 这里， @p points 是 @p spacedim 维度上的点的集合， @p a
 * 是相应权重的集合。在这种情况下，点将是单元格、面或边的顶点，而当需要一个新的单元格、面或边的中点时，权重通常是点的数量的1。然后，派生类将以计算这个新点的位置的方式实现
 * Manifold::get_new_point()
 * 函数。在最简单的情况下，例如在FlatManifold类中，该函数只是计算给定点的算术平均数（有给定权重）。然而，其他类的做法有所不同；例如，SphericalManifold类，它用于描述形成球体（部分）的域，将确保在边界处给定边缘的两个顶点，新返回的点将位于连接这两点的大圆上，而不是选择一个位于
 * ${\mathbb R}^d$ 中两点之间的一半的点。
 *
 *
 *
 * @note
 * 与库中几乎所有其他情况不同，我们在这里将点解释为在实空间，而不是在参考单元上。
 * Manifold::get_new_point()
 * 有一个默认的实现，可以在一定程度上简化这个过程：在内部，这个函数调用
 * Manifold::get_intermediate_point()
 * 来计算成对的中间点。在内部，
 * Manifold::get_intermediate_point()
 * 在计算完给定点的凸组合后调用 Manifold::project_to_manifold()
 * 函数。 这允许派生类只在简单的情况下重载
 * Manifold::project_to_manifold()
 * 。在描述嵌入高维空间的流形时，例如球体的表面，这通常很有用。
 * 在这些情况下，所需的新点可以简单地通过所提供的点的（加权）平均值来计算，再投射到球体上。
 *
 *  <h3>Common use case: Computing tangent vectors</h3>
 * 这个类别的第二个用途是在计算域和边界上的方向。例如，我们可能需要计算一个面的法向量，以便施加无流边界条件
 * $\mathbf u \cdot \mathbf n = 0$ （见
 * VectorTools::compute_no_normal_flux_constraints()
 * 的例子）。同样地，我们在计算数值解的梯度的法向分量时可能需要法向量，以便计算误差估计器中解的梯度的跳跃（例如，见KellyErrorEstimator类）。
 * 为了实现这一点，流形类提供了一个成员函数（由派生类实现），通过
 * Manifold::get_tangent_vector() 函数计算
 * "在某一点与流形相切的矢量，在另一点的方向"。例如，在2D中，我们可以用这个函数在边界上的两条边的顶点上计算一个沿边的
 * "切向
 * "矢量，然后通过旋转90度得到法向量。在3D中，我们可以计算与边界顶点相邻的边界面的两条边的
 * "切向
 * "向量，然后取这两个向量的交积来获得边界的法向量。
 * 由于比较难以理解的原因，这些方向向量以一种非常特殊的方式被归一化，而不是具有单位规范。更多信息请参见
 * Manifold::get_tangent_vector(), 的文档以及下文。
 * 在最简单的情况下（即FlatManifold类），这些切向量只是两个给定点之间的差向量。然而，在更复杂（和更有趣）的情况下，方向可能是不同的。例如，对于SphericalManifold的情况，如果给定的两个点位于围绕原点的共同大圆上，那么切向量将与大圆相切，而不是从一个点直指另一个点。
 *
 *  <h3>A unified description</h3> 要理解这个类的作用，"真正的
 * "方法是在微分几何的框架内看它。更具体地说，微分几何从根本上是基于这样的假设：两个足够近的点通过一条
 * "最短距离 "的线相连。这条线被称为
 * "测地线"，它从连接两点的所有其他线中被挑选出来，因为如果距离是以描述流形的
 * "公制
 * "来衡量的，它就是最短的。为了举例说明，请回忆一下，平流形的测地线（在FlatManifold类中实现）只是连接两点的直线，而对于球面流形（见SphericalManifold类），距离相同的两点之间的测地线是大圆，当连接离原点不同距离的两条线时，一般是弯曲的线。
 * 在下面的讨论中，以及为了实现当前类的目的，对微分几何学来说非常基本的
 * "度量
 * "概念对我们来说不再重要了。相反，一切都可以简单地通过假设流形上连接点的测地线的存在来描述。
 * 考虑到测地线，前两节讨论的操作可以用更正式的方式来描述。从本质上讲，它们依赖于这样一个事实：我们可以假设一个测地线是由一个类似于
 * "时间 "的变量 $t$ 来参数化的，因此 $\mathbf s(t)$
 * 描述了曲线，因此 $\mathbf s(0)$ 是第一个点的位置，
 * $\mathbf s(1)$ 是第二个点的位置。此外， $\mathbf s(t)$
 * 以恒定的速度追踪测地线，在相等的时间内覆盖相等的距离（以公制来衡量）。请注意，这个参数化使用时间而不是弧长来表示沿测地线的进展。
 * 在这幅图中，计算点 $\mathbf x_1$ 和 $\mathbf x_2$
 * 之间的中点，权重为 $w_1$ 和 $w_2=1-w_1$ ，只需要计算点
 * $\mathbf s(w_1)$
 * 。计算一个新的点作为两个以上的点的加权平均，可以通过考虑成对的测地线来完成，在前两个点之间的测地上找到合适的点，然后在这个新点和第三个给定点之间的测地上找到合适的点，等等。
 * 同样，上述的 "切向 "矢量只是速度矢量， $\mathbf s'(t)$
 * ，在测地线的一个端点（即在 $t=0$ 或 $t=1$
 * ）评估。在平坦流形的情况下，测地线只是连接两点的直线，速度矢量只是这种情况下的连接矢量。另一方面，对于球面流形上的两点，测地线是一个大圆，速度矢量是与球面的切线。
 * 请注意，如果我们愿意，我们可以利用这一点来计算连接两点
 * $\mathbf x_1$ 和 $\mathbf x_2$
 * 的测地线的长度，即沿着连接它们的测地线计算 $\int_0^1
 * \|\mathbf s'(t)\| dt$
 * ，但这种操作在实践中对我们没有用。我们也可以设想使用上面的
 * "新点 "操作来计算方向向量，使用公式 $\mathbf
 * s'(0)=\lim_{w\rightarrow 0} \frac{\mathbf s(w)-\mathbf s(0)}{w}$
 * ，我们需要做的就是沿着连接 $\mathbf x_1$ 和 $\mathbf x_2$
 * 的测地线计算新点 $\mathbf s(w)$ ，权重为 $w$ 和 $1-w$
 * 。该函数的默认实现是这样做的，对一个小但有限的权重
 * $w$
 * 进行商数评估。然而，在实践中，几乎总是可以明确地计算方向向量，也就是说，不需要对极限过程进行数值近似，派生类应该这样做。
 *
 *
 *
 * @ingroup manifold
 *
 */
template <int dim, int spacedim = dim>
class Manifold : public Subscriptor
{
public:
  // explicitly check for sensible template arguments
  static_assert(dim <= spacedim,
                "The dimension <dim> of a Manifold must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");


  /**
   * 类型保持单元格面的顶点的法线信息。因此，有
   * <tt>GeometryInfo<dim>::vertices_per_face</tt>
   * 个法线向量，它们定义了顶点处边界的切线空间。请注意，存储在这个对象中的向量不需要被归一化，也不需要实际指向外部，因为人们通常只想检查正交性以定义切平面；如果一个函数要求法线被归一化，那么它必须自己这样做。
   * 由于明显的原因，这种类型在1d中是没有用的。
   *
   */
  using FaceVertexNormals =
    std::array<Tensor<1, spacedim>, GeometryInfo<dim>::vertices_per_face>;


  /**
   * 毁灭器。在这里不做任何事情，但需要声明为虚，以便使从这个类派生的类的层次结构成为可能。
   *
   */
  virtual ~Manifold() override = default;

  /**
   * 返回此流形的一个副本。
   * 每个派生类都应该以一种合理的方式实现这一操作。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const = 0;

  /**
   * @name  计算点的位置。
   *
   */
  /// @{

  /**
   * 返回两个给定点之间的一个中间点。重载这个函数可以让默认的对偶缩小法的实现get_new_point()正常工作，该方法将一个正交对象作为输入。
   * 这个函数的实现应该在流形上返回一条参数曲线，连接点`p1`和`p2`，参数`w`在区间[0,1]。特别是`get_intermediate_point(p1,
   * p2, 0.0)`应该返回`p1`，`get_intermediate_point(p1, p2,
   * 1.0)`应该返回`p2`。
   * 在其默认实现中，这个函数用`p1`和`p2`的凸形组合调用project_to_manifold()方法。用户类可以通过简单地实现project_to_manifold()方法来摆脱。
   *
   */
  virtual Point<spacedim>
  get_intermediate_point(const Point<spacedim> &p1,
                         const Point<spacedim> &p2,
                         const double           w) const;

  /**
   * 返回将成为新顶点的点，该点被给定的点所包围  @p
   * surrounding_points.   @p weights
   * 包含周围点的适当权重，流形根据该权重决定新点的位置。
   * 在其默认的实现中，它采用了成对减少点的方法，在前两个点上调用函数get_intermediate_point()，然后在产生的点和下一个点上调用，直到矢量中的所有点都被考虑在内。用户类可以通过简单地实现get_intermediate_point()函数来实现。请注意，默认情况下，get_intermediate_point()函数调用project_to_manifold()函数的参数的凸组合。对于简单的情况，你可以只实现project_to_manifold()函数。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double> &         weights) const;


  /**
   * 计算一组新的点，在给定的点 @p
   * surrounding_points之间进行插值。  @p weights
   * 是一个表格，其列数与 @p  surrounding_points.size()相同。 @p
   * weights 中的行数必须与 @p new_points. 的长度相匹配。
   * 在其默认实现中，该函数只是在 @p weights
   * 的每一行上调用get_new_point()，并将这些点写入输出数组
   * @p new_points.
   * 。然而，如果像MappingQGeneric中那样需要生成多个新点，并且流形在图表空间和物理空间之间进行昂贵转换，该函数就更加高效。对于这个函数，周围的点只需要被转换回图表稀疏空间一次，而不是每次调用get_new_point()。如果效率并不重要，你可以通过只实现get_new_point()函数来解决。
   * 该实现不允许 @p surrounding_points 和 @p new_points
   * 指向同一个数组，所以要确保将不同的对象传入该函数。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double> &                weights,
                 ArrayView<Point<spacedim>>              new_points) const;

  /**
   * 给定一个接近给定流形的点，它对其进行修改并将其投射到流形本身。
   * 该类被函数get_new_point()的默认实现所使用，并且应该被派生类所实现。默认实现只是在调用时抛出一个异常。
   * 如果你的manifold很简单，你可以只实现这个函数，默认行为应该是开箱即用的。
   *
   */
  virtual Point<spacedim>
  project_to_manifold(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const Point<spacedim> &                 candidate) const;

  /**
   * 向后兼容接口。
   * 返回应成为一条规则线的两个子线的新的中间顶点的点。在2D中，这条线是边界上的一条线，而在3D中，它是边界上的一个面的边界（因此线也在边界上）。
   * 这个函数的默认实现将其参数传递给
   * Manifolds::get_default_points_and_weights() 函数，然后调用
   * Manifold<dim,spacedim>::get_new_point()
   * 函数。用户派生类可以重载
   * Manifold<dim,spacedim>::get_new_point() 或
   * Manifold<dim,spacedim>::project_to_manifold(), ，该函数被
   * Manifold<dim,spacedim>::get_new_point(). 的默认实现所调用。
   *
   */
  virtual Point<spacedim>
  get_new_point_on_line(
    const typename Triangulation<dim, spacedim>::line_iterator &line) const;

  /**
   * 向后兼容接口。返回在三个或更多空间维度上将成为一个四边形的四个子体的边界的公共点。因此，这个函数只在至少三个维度上有用，对于更低的维度不应调用。
   * 这个函数是在限定给定 @p quad
   * 的四条线被细化后调用的，所以你可能想使用<tt>quad->line(i)->child(j)</tt>提供的信息，<tt>i=0...3</tt>,
   * <tt>j=0,1</tt>。    这个函数的默认实现将其参数传递给
   * Manifolds::get_default_points_and_weights() 函数，然后调用
   * Manifold<dim,spacedim>::get_new_point()
   * 函数。用户派生类可以重载
   * Manifold<dim,spacedim>::get_new_point() 或
   * Manifold<dim,spacedim>::project_to_manifold(), ，该函数被
   * Manifold<dim,spacedim>::get_new_point(). 的默认实现所调用。
   *
   */
  virtual Point<spacedim>
  get_new_point_on_quad(
    const typename Triangulation<dim, spacedim>::quad_iterator &quad) const;

  /**
   * 向后兼容接口。
   * 返回在三维或空间维度上将成为一个六边形的八个子节点的公共点。因此，这个函数只在至少三个维度上有用，对更低的维度不应调用。
   * 这个函数是在给定的 @p hex
   * 的所有边界对象被细化后调用的，所以你可能想使用<tt>hex->quad(i)->line(j)->child(k)</tt>提供的信息，<tt>i=0...5</tt>,
   * <tt>j=0...3</tt>, <tt>k=0,1</tt>。
   * 这个函数的默认实现将其参数传递给
   * Manifolds::get_default_points_and_weights() 函数，然后调用
   * Manifold<dim,spacedim>::get_new_point()
   * 函数。用户派生类可以重载
   * Manifold<dim,spacedim>::get_new_point() 或
   * Manifold<dim,spacedim>::project_to_manifold(), ，该函数被
   * Manifold<dim,spacedim>::get_new_point(). 的默认实现所调用。
   *
   */
  virtual Point<spacedim>
  get_new_point_on_hex(
    const typename Triangulation<dim, spacedim>::hex_iterator &hex) const;


  /**
   * 向后兼容接口。根据<tt>dim=2</tt>或<tt>dim=3</tt>这个函数调用get_new_point_on_line或get_new_point_on_quad函数。对于<tt>dim=1</tt>，它会抛出一个异常。这个封装器允许独立于维度的编程。
   *
   */
  Point<spacedim>
  get_new_point_on_face(
    const typename Triangulation<dim, spacedim>::face_iterator &face) const;


  /**
   * 向后兼容接口。 根据<tt>dim=1</tt>, <tt>dim=2</tt> 或
   * <tt>dim=3</tt> 这个函数调用get_new_point_on_line,
   * get_new_point_on_quad 或者 get_new_point_on_hex
   * 函数。这个封装器允许独立于维度的编程。
   *
   */
  Point<spacedim>
  get_new_point_on_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

  /// @}

  /**
   * @name  计算正切向量
   *
   */
  /// @{

  /**
   * 返回一个在  $\mathbf x_1$  处与连接两点  $\mathbf x_1,\mathbf
   * x_2$
   * 的测地线相切的向量。测地线是这两点之间最短的线，其中
   * "最短
   * "是通过派生类中该类的特定实现的度量来定义的。例如，在FlatManifold的情况下，两点之间的最短线只是直线，在这种情况下，切线矢量只是差值
   * $\mathbf d=\mathbf x_2-\mathbf x_1$
   * 。另一方面，对于描述嵌入高维空间的表面的流形（例如，球体的表面），那么切向量是与表面相切的，因此可能指向与连接两点的直线不同的方向。
   * 虽然切向量经常被归一化为单位长度，但这个函数返回的向量是归一化的，如本类的介绍中所述。具体来说，如果
   * $\mathbf s(t)$ 在 $\mathbf x_1 = \mathbf s(0)$ 和 $\mathbf x_2 =
   * \mathbf s(1)$
   * 两点之间追踪出测地线，那么返回的向量必须等于
   * $\mathbf s'(0)$
   * 。换句话说，从某种意义上说，返回的向量也编码了测地线的<i>length</i>，因为如果曲线
   * $\mathbf s(t)$ 在参数 $t=0$ 和 $t=1$
   * 之间连接的两点距离较远，那么它的移动速度一定会更快。
   * 这个函数的默认实现是对 $\epsilon$ 的一个小值进行逼近
   * $\mathbf s'(0) \approx \frac{\mathbf s(\epsilon)-\mathbf x_1}{\epsilon}$
   * ，而对 $\mathbf
   * s(\epsilon)$
   * 的评估是通过调用get_new_point()完成的。如果可能的话，派生类应该通过精确导数的实现来覆盖这个函数。
   * @param  x1 描述测地线的第一个点，也是要评估 "方向
   * "的那个点。    @param  x2 描述测地线的第二个点。
   * @return  一个与测地线相切的 "方向 "向量。
   *
   */
  virtual Tensor<1, spacedim>
  get_tangent_vector(const Point<spacedim> &x1,
                     const Point<spacedim> &x2) const;

  /// @}

  /**
   * @name  计算法向量
   *
   */
  /// @{

  /**
   * 如果p实际上不在表面上，而只是在附近，则尽量返回一些合理的东西，例如最接近p的表面点的法向量。（实际上p点通常不在实际表面上，而是由一些多项式映射的正交点；然而，映射的表面通常不会与实际表面重合）。
   * 这个函数只有在dim==spacedim的情况下才有意义，因为否则就没有唯一的法线矢量，而实际上是一个(spacedim-dim+1)-dimensional的矢量切线空间，这些矢量既是面的法线，又是生活在spacedim-dimensional空间中的dim-dimensional曲面的法线。例如，想想一个二维网格，它覆盖了三维空间中的一个二维表面。在这种情况下，每个面（边）都是一维的，有两个线性独立的向量都是对边的法线：一个是对边的法线，并与曲面相切（直观地说，如果曲面是局部平坦的，那就是从当前单元指向邻近的单元），另一个是以边为根，但指向垂直于曲面（这也是垂直于住在曲面中的边）。
   * 因此，因为如果spacedim大于dim，这个函数没有明显正确的语义，所以在这种情况下，这个函数会直接抛出一个错误。
   * 面的迭代器给出了这个函数应该为哪个面计算法线向量的指示。
   * 如果域的边界是由不同的非差分组成的，这就很有用了（例如，当使用FlatManifold类来近似一个完全由粗网格描述的几何体时，顶点之间有成片的（双）线性分量，但边界可能在顶点本身有一个扭结）。
   * @note
   * 在2D中，这个函数的默认实现是通过从p到构成边缘的两个顶点中的另一个顶点的切线方向来计算法向量，然后将其向外旋转（相对于边缘的坐标系）90度。在3D中，默认的实现方式更为复杂，目的是为了避免靠近其中一个顶点的点的数值舍入问题，并避免切线方向的线性依赖。
   *
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &                                     p) const;

  /**
   * 计算给定面嵌入漫游中的每个顶点到边界的法向量。并不要求法向量以某种方式被规范化。
   * 也不要求法线实际指向外部。
   * 这个函数是用来计算C1映射的数据的。默认实现是在每个顶点上调用normal_vector()。
   * 注意，当计算边界不可微调的顶点的法线向量时，你必须确保你计算的是单边极限，即关于给定面内的点的极限。
   *
   */
  virtual void
  get_normals_at_vertices(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    FaceVertexNormals &face_vertex_normals) const;

  /// @}
};


/**
 * Manifold<dim,spacedim>的特化，它代表了嵌入 @p dim
 * 维的欧氏空间中的可能的周期性欧氏空间。这个曼福尔德的主要特征是，函数
 * FlatManifold<dim,spacedim>::project_to_manifold() 是身份函数。
 *
 *
 * @ingroup manifold
 *
 *
 */
template <int dim, int spacedim = dim>
class FlatManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * 默认构造函数。可选的参数可以用来指定间隔维流形的周期性（每个方向一个周期）。周期性值为0意味着沿该方向没有周期性。默认情况下，不假定有周期性。
   * 周期性会影响到计算中间点的方式。假设两点之间的距离超过半个周期，那么应该通过跨越周期性边界来计算距离，也就是说，通过在两点之和上加一个完整的周期来计算平均值。例如，如果沿0方向我们有2*pi的周期性，那么（2*pi-eps）和（eps）的平均值就不是pi，而是2*pi（或零），因为在一个周期性流形上，这两个点的距离是2*eps，而不是（2*pi-eps）。特殊情况会被考虑在内，以确保行为总是符合预期。第三个参数在计算距离时被用作相对公差。
   * 周期性将以下列方式进行：域被认为是包含在[Point<spacedim>(),
   * periodicity)中的盒子，其中右边的极端被排除。如果这个盒子的任何一个分量的长度为零，那么在这个方向上就没有周期性。
   * 每当一个试图计算平均数的函数被调用时，如果你用来计算平均数的一个点位于周期性框外，就会产生一个异常。返回的点保证位于周期性框内，加上或减去公差*周期性.规范（）。
   *
   */
  FlatManifold(const Tensor<1, spacedim> &periodicity = Tensor<1, spacedim>(),
               const double               tolerance   = 1e-10);

  /**
   * 返回该流形的一个副本。
   *
   */
  virtual std::unique_ptr<Manifold<dim, spacedim>>
  clone() const override;

  /**
   * 让新点成为周围顶点的平均和。
   * 这种特殊的实现方式构建了周围点的加权平均，然后在内部调用函数project_to_manifold()。我们之所以这样做，是为了让懒惰的程序员只为自己的流形类实现project_to_manifold()函数，这些流形类是平流形的小（或微不足道的）扰动。只要粗略的网格是流形几何的一个体面的近似，就会出现这种情况。在这种情况下，单元格的中间点接近流形的真正中间点，投影可能就足够了。
   * 对于大多数简单的几何形状，通过从FlatManifold派生出自己的Manifold类，并只为project_to_manifold函数编写一个新的接口，可以得到合理的结果。只要在你试图细化的最粗的网格尺寸中，中间点离流形中点不是太远，也就是说，只要粗的网格尺寸足够小，你也会有很好的近似结果。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double> &         weights) const override;


  /**
   * 计算一组新的点，在给定的点之间进行插值  @p
   * surrounding_points。  @p weights 是一个列数与 @p
   * surrounding_points.size()一样多的表。 @p weights
   * 中的行数必须与 @p new_points. 的长度相匹配。
   * 对于这个特定的实现，根据 @p weights 对 @p surrounding_points
   * 的插值在笛卡尔空间中简单地执行。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double> &                weights,
                 ArrayView<Point<spacedim>> new_points) const override;

  /**
   * 投射到FlatManifold。这是对平坦的欧几里得空间的认同函数。但是请注意，这个函数可以被派生类重载，然后受益于get_new_point()函数背后的逻辑，这些逻辑通常与这个类中实现的逻辑非常相似（如果不是完全相同）。
   *
   */
  virtual Point<spacedim>
  project_to_manifold(const ArrayView<const Point<spacedim>> &points,
                      const Point<spacedim> &candidate) const override;

  /**
   * 返回一个在 $\mathbf x_1$ 处与连接两点 $\mathbf x_1,\mathbf
   * x_2$ 的测地线相切的向量。
   * 对于当前类，我们假设流形是平的，所以测地线是两点之间的直线，我们返回
   * $\mathbf x_2-\mathbf x_1$  。选择矢量的归一化，使其符合
   * Manifold::get_tangent_vector().  中描述的惯例。
   * @note
   * 如果你把这个类作为垫脚石，通过重载project_to_manifold()函数，建立一个只
   * "稍微 "偏离平面流形的流形。      @param  x1
   * 描述测地线的第一个点，也是要评估 "方向 "的那个点。
   * @param  x2 描述测地线的第二个点。    @return
   * 一个与测地线相切的 "方向 "向量。这里，这是 $\mathbf
   * x_2-\mathbf x_1$
   * ，可能被构造函数中设置的域的周期性所修改，以在必要时使用通过周期边界的各点之间的
   * "最短 "连接。
   *
   */
  virtual Tensor<1, spacedim>
  get_tangent_vector(const Point<spacedim> &x1,
                     const Point<spacedim> &x2) const override;

  /**
   * 返回p点处给定面的法向量，考虑到三维六面体单元的四边形面可能不是平面的。
   * 在这些情况下，假设该面有一个由双线性函数描述的几何形状，法向量是通过将这个双线性形式嵌入一个具有平面公制的笛卡尔空间来计算的。
   *
   */
  virtual Tensor<1, spacedim>
  normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim> &p) const override;

  /**
   * 计算给定面的每个顶点到边界的法向量，考虑到三维六面体单元的四边形面可能不是平面的。在这些情况下，假定面的几何形状是由双线性函数描述的，法向量的计算是通过将这个双线性形式嵌入一个具有平面度量的笛卡尔空间。
   *
   */
  virtual void
  get_normals_at_vertices(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
    const override;

  /**
   * 返回该芒果的周期性。
   *
   */
  const Tensor<1, spacedim> &
  get_periodicity() const;

private:
  /**
   * 这个流形的周期性。周期性会影响到计算中间点的方式。假设两点之间的距离超过半个周期，那么就应该通过跨越周期性边界来计算距离，也就是说，通过在两点之和上加一个完整的周期来计算平均值。例如，如果沿0方向我们有2*pi的周期性，那么（2*pi-eps）和（eps）的平均值就不是pi，而是2*pi（或0），因为在一个周期性流形上，这两点的距离是2*eps而不是（2*pi-eps）。
   * 沿着一个方向的周期性为0意味着没有周期性。这是对所有方向的默认值。
   *
   */
  const Tensor<1, spacedim> periodicity;

  DeclException3(ExcPeriodicBox,
                 int,
                 Point<spacedim>,
                 double,
                 << "The component number " << arg1 << " of the point [ "
                 << arg2 << " ] is not in the interval [ 0, " << arg3
                 << "), bailing out.");

  /**
   * 相对公差。这个公差用于计算双精度的距离。
   *
   */
  const double tolerance;
};


/**
 * 这个类描述了可以用图表来表达的映射关系。具体来说，这个类及其模板参数描述了一个尺寸为chartdim的图表，它是Manifold<dim,spacedim>的一部分，并被用于一个类型为Triangulation<dim,spacedim>的对象中：它特化了一个尺寸为chartdim的Manifold嵌入到尺寸为spacedim的manifold中，对于它，你有明确的pull_back（）和push_forward（）转换方式。它的使用在
 * step-53  中有详细的解释。
 * 这是一个辅助类，当你有一个从维数为chartdim的欧几里得空间到维数为spacedim的欧几里得空间的显式映射时，这个类很有用，也就是说，当你的流形
 * $\mathcal{M}$ 可以由一个映射\f[ F: \mathcal{B} \subset
 * R^{\text{chartdim}} \mapsto \mathcal{M} \subset R^{\text{spacedim}}
 * \f]（push_forward()函数）表示，并且接受反向变换\f[ F^{-1}:
 * \mathcal{M} \subset R^{\text{spacedim}} \mapsto \mathcal{B} \subset
 * R^{\text{chartdim}} \f]（pull_back()函数）。
 * ChartManifold类的get_new_point()函数是通过调用所有<tt>surrounding_points</tt>的pull_back()方法来实现的，计算它们在chartdim
 * Euclidean空间的加权平均值，然后用得到的点调用push_forward()方法，即\f[
 * \mathbf x^{\text{new}} = F(\sum_i w_i F^{-1}(\mathbf x_i)).  \f] 。
 * 派生类需要实现push_forward()和pull_back()方法。然后，所有其他映射需要的函数（除了push_forward_gradient()函数，见下文）都将由这个类提供。
 *
 *  <h3>Providing function gradients</h3>
 * 为了计算与流形相切的向量（例如，与嵌入高维空间的曲面相切，或者仅仅是
 * ${\mathbb R}^3$
 * 的三个单位向量），需要同时获得push-forward函数的<i>gradient</i>。梯度是矩阵
 * $(\nabla F)_{ij}=\partial_j F_i$ ，其中我们对 $\mathcal B$
 * 所在的平坦欧氏空间上的chartdim参考坐标进行导数。换句话说，在
 * $\mathbf x$ 点， $\nabla F(\mathbf x)$ 是一个大小为 @p spacedim 乘
 * @p chartdim. 的矩阵。 只有 ChartManifold::get_tangent_vector()
 * 函数使用前推的梯度，但所有有限元代码中只有一个子集实际需要计算切向量。因此，虽然派生类需要实现该类的抽象虚拟push_forward()和pull_back()函数，但它们不需要实现虚拟push_forward_gradient()函数。相反，该函数有一个默认的实现（因此不是抽象的，因此不强迫派生类重载它），但是默认的实现显然不能计算任何有用的东西，因此只是触发了一个异常。
 *
 *  <h3>A note on the template arguments</h3> 尺寸参数  @p chartdim,   @p
 * dim  和  @p spacedim  必须满足以下关系。
 * @code
 *    dim <= spacedim
 *    chartdim <= spacedim
 * @endcode
 * 然而， @p dim 和 @p chartdim. 之间没有先验的关系。 例如，如果你想描述嵌入三维空间的二维三角形中的一个边（一个一维物体）的映射，你可以通过线@f[
 *    F: [0,1] \rightarrow {\mathbb R}^3
 * @f]来进行参数化。 另一方面，我们没有理由不把它描述为一个映射@f[
 * F: {\mathbb R}^3 \rightarrow {\mathbb R}^3 @f]，使线 $[0,1]\times
 * \{0\}\times \{0\}$ 刚好映射到有关的边上。这里， @p chartdim
 * 是3。这可能看起来很麻烦，但只要有可能从边到回拉空间然后再返回，就能满足可逆函数
 * $F$
 * 的要求。最后，考虑到我们正在处理三维中的二维三角，我们经常会有一个从例如二维单位方块或单位圆盘到三维空间中的域的映射，而有关的边可能只是二维空间中单位域的映射的边。在这种情况下，
 * @p  chartdim是2。
 *
 *
 * @ingroup manifold
 *
 */
template <int dim, int spacedim = dim, int chartdim = dim>
class ChartManifold : public Manifold<dim, spacedim>
{
public:
  // explicitly check for sensible template arguments
  static_assert(dim <= spacedim,
                "The dimension <dim> of a ChartManifold must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");

  /**
   * 构造函数。可选参数可用于指定chartdim-dimensional流形的周期性（每个方向一个周期）。周期性值为0意味着沿该方向没有周期性。默认情况下，不假定有周期性。
   * 周期性会影响到计算中间点的方式。假设两点之间的距离超过半个周期，那么应该通过跨越周期性边界来计算距离，也就是说，然后通过在两点之和上加一个完整的周期来计算平均值。例如，如果沿0方向我们有2*pi的周期性，那么（2*pi-eps）和（eps）的平均值就不是pi，而是2*pi（或0），因为在流形上，这两个点的距离是2*eps，而不是（2*pi-eps）。
   *
   */
  ChartManifold(const Tensor<1, chartdim> &periodicity = Tensor<1, chartdim>());

  /**
   * 破坏器。在这里不做任何事情，但需要声明以使其成为虚拟的。
   *
   */
  virtual ~ChartManifold() override = default;

  /**
   * 更多信息请参考该类的一般文档和基类的文档。
   *
   */
  virtual Point<spacedim>
  get_intermediate_point(const Point<spacedim> &p1,
                         const Point<spacedim> &p2,
                         const double           w) const override;

  /**
   * 请参考这个类的一般文档和基类的文档以了解更多信息。
   *
   */
  virtual Point<spacedim>
  get_new_point(const ArrayView<const Point<spacedim>> &surrounding_points,
                const ArrayView<const double> &         weights) const override;

  /**
   * 计算一组新的点，在给定的点之间进行插值  @p
   * surrounding_points。  @p weights 是一个表，其列数与 @p
   * surrounding_points.size()相同。 @p weights 中的行数必须与 @p
   * new_points. 的长度相匹配。
   * 这个函数的实现首先通过调用pull_back()将 @p
   * surrounding_points 转换到图表空间。然后，根据给定的 @p
   * weights,
   * ，通过通常的插值在图表上计算新的点，最后通过push_forward()转换到图像空间。
   * 在pull_back()操作昂贵的情况下，与单独调用get_new_point()相比，这种实现方式对于从相同的周围点计算多个新点要高效得多。这是因为pull_back()只对周围的点调用一次，并且使用这组点对所有给定的权重进行插值。通常，pull_back()也比push_forward()昂贵，因为前者可能涉及到非三流形中的某种牛顿迭代。
   *
   */
  virtual void
  get_new_points(const ArrayView<const Point<spacedim>> &surrounding_points,
                 const Table<2, double> &                weights,
                 ArrayView<Point<spacedim>> new_points) const override;
  /**
   * 将spacedim中的给定点拉回欧几里得图表维度空间。
   * 更多信息请参考该类的一般文档。
   *
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const = 0;

  /**
   * 给定chartdim维度欧氏空间中的一个点，该方法返回嵌入spacedim欧氏空间的流形上的一个点。
   * 更多信息请参考该类的一般文档。
   *
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const = 0;

  /**
   * 给定Chartdim维欧几里得空间中的一个点，该方法返回从Chartdim维空间映射到spacedim维空间的函数
   * $F$ 的导数。换句话说，它是一个大小为
   * $\text{spacedim}\times\text{chartdim}$ 的矩阵。
   * 这个函数被用于get_tangent_vector()函数所要求的计算中。由于不是所有Manifold类接口的用户都需要调用该函数，所以当前的函数被实现了，但每当调用时都会触发一个异常。这允许派生类避免实现push_forward_gradient函数，如果用户程序中不需要这个功能。
   * 更多信息请参考该类的一般文档。
   *
   */
  virtual DerivativeForm<1, chartdim, spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

  /**
   * 返回一个在 $\mathbf x_1$ 处与连接两点 $\mathbf x_1,\mathbf
   * x_2$ 的测地线相切的向量。
   * 更详细的描述，请参见Manifold类和
   * Manifold::get_tangent_vector() 的文档。
   * 对于当前的类，我们假设这个测地线是 @p x1 和 @p x2
   * 的预像（其中预像是通过拉回位置 @p x1 和 @p x2).
   * 来计算的）的直线在push_forward()操作下的图像
   * 换句话说，如果这些预像是 $\xi_1=F^{-1}(\mathbf x_1),
   * \xi_2=F^{-1}(\mathbf x_2)$
   * ，那么测地线在预像（chartdim维欧几里得）空间是
   * @f{align*}{
   * \zeta(t) &= \xi_1 +  t (\xi_2-\xi_1)
   * \\          &= F^{-1}(\mathbf x_1) + t\left[F^{-1}(\mathbf x_2)
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -F^{-1}(\mathbf x_1)\right]
   * @f}
   * 在图像空间中，即在我们操作的空间中，这导致了曲线的出现
   * @f{align*}{
   * \mathbf s(t) &= F(\zeta(t))
   * \\          &= F(\xi_1 +  t (\xi_2-\xi_1))
   * \\          &= F\left(F^{-1}(\mathbf x_1) + t\left[F^{-1}(\mathbf x_2)
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -F^{-1}(\mathbf x_1)\right]\right).
   * @f}
   * 当前函数应该返回的是  $\mathbf s'(0)$
   * 。根据连锁规则，这等于
   * @f{align*}{
   * \mathbf s'(0) &=
   *   \frac{d}{dt}\left. F\left(F^{-1}(\mathbf x_1)
   *                      + t\left[F^{-1}(\mathbf x_2)
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -F^{-1}(\mathbf x_1)\right]\right)
   *               \right|_{t=0}
   * \\ &= \nabla_\xi F\left(F^{-1}(\mathbf x_1)\right)
   *                  \left[F^{-1}(\mathbf x_2)
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -F^{-1}(\mathbf x_1)\right].
   * @f}
   * 然后这个公式可能要稍作修改，考虑到在调用构造函数时假定的任何周期性。
   * 因此，正切向量的计算也需要实现<i>derivatives</i>
   * $\nabla_\xi F(\xi)$ 的前推映射。这里， $F^{-1}(\mathbf
   * x_2)-F^{-1}(\mathbf x_1)$ 是一个chartdim维的向量，而 $\nabla_\xi
   * F\left(F^{-1}(\mathbf
   * x_1)\right) = \nabla_\xi F\left(\xi_1\right)$
   * 是一个spacedim-times-chartdim维的矩阵。因此，正如所期望的那样，该操作的结果是一个间隔维的向量。
   * @param  x1 描述测地线的第一个点，也是 "方向
   * "要被评估的点。    @param  x2 描述测地线的第二个点。
   * @return  一个与测地线相切的 "方向 "向量。
   *
   */
  virtual Tensor<1, spacedim>
  get_tangent_vector(const Point<spacedim> &x1,
                     const Point<spacedim> &x2) const override;

  /**
   * 返回与子模子相关的周期性。
   *
   */
  const Tensor<1, chartdim> &
  get_periodicity() const;

private:
  /**
   * sub_manifold对象用于计算图表坐标系中各点的平均值。
   * 在一个理想的世界里，它的类型是FlatManifold<dim,chartdim>。然而，这将实例化dim>spacedim的情况，从而导致无效的情况。我们改用<chartdim,chartdim>，这(i)总是有效的，而且(ii)根本不重要，因为就流形功能而言，流形的第一个（dim）参数事实上是被忽略的。
   *
   */
  const FlatManifold<chartdim, chartdim> sub_manifold;
};


 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

template <>
Point<1>
Manifold<1, 1>::get_new_point_on_face(
  const Triangulation<1, 1>::face_iterator &) const;

template <>
Point<2>
Manifold<1, 2>::get_new_point_on_face(
  const Triangulation<1, 2>::face_iterator &) const;


template <>
Point<3>
Manifold<1, 3>::get_new_point_on_face(
  const Triangulation<1, 3>::face_iterator &) const;


template <>
Point<1>
Manifold<1, 1>::get_new_point_on_quad(
  const Triangulation<1, 1>::quad_iterator &) const;

template <>
Point<2>
Manifold<1, 2>::get_new_point_on_quad(
  const Triangulation<1, 2>::quad_iterator &) const;


template <>
Point<3>
Manifold<1, 3>::get_new_point_on_quad(
  const Triangulation<1, 3>::quad_iterator &) const;


template <>
Point<3>
Manifold<3, 3>::get_new_point_on_hex(
  const Triangulation<3, 3>::hex_iterator &) const;

 /*---Templated functions---*/ 

namespace Manifolds
{
  template <typename MeshIteratorType>
  std::pair<std::array<Point<MeshIteratorType::AccessorType::space_dimension>,
                       n_default_points_per_cell<MeshIteratorType>()>,
            std::array<double, n_default_points_per_cell<MeshIteratorType>()>>
  get_default_points_and_weights(const MeshIteratorType &iterator,
                                 const bool              with_interpolation)
  {
    const int dim      = MeshIteratorType::AccessorType::structure_dimension;
    const int spacedim = MeshIteratorType::AccessorType::space_dimension;
    constexpr std::size_t points_per_cell =
      n_default_points_per_cell<MeshIteratorType>();

    std::pair<std::array<Point<spacedim>, points_per_cell>,
              std::array<double, points_per_cell>>
      points_weights;


    // note that the exact weights are chosen such as to minimize the
    // distortion of the four new quads from the optimal shape; their
    // derivation and values is copied over from the
    // interpolation function in the mapping
    switch (dim)
      {
        case 1:
          Assert(points_weights.first.size() == 2, ExcInternalError());
          Assert(points_weights.second.size() == 2, ExcInternalError());
          points_weights.first[0]  = iterator->vertex(0);
          points_weights.second[0] = .5;
          points_weights.first[1]  = iterator->vertex(1);
          points_weights.second[1] = .5;
          break;
        case 2:
          Assert(points_weights.first.size() == 8, ExcInternalError());
          Assert(points_weights.second.size() == 8, ExcInternalError());

          for (unsigned int i = 0; i < 4; ++i)
            {
              points_weights.first[i] = iterator->vertex(i);
              points_weights.first[4 + i] =
                (iterator->line(i)->has_children() ?
                   iterator->line(i)->child(0)->vertex(1) :
                   iterator->line(i)->get_manifold().get_new_point_on_line(
                     iterator->line(i)));
            }

          if (with_interpolation)
            {
              std::fill(points_weights.second.begin(),
                        points_weights.second.begin() + 4,
                        -0.25);
              std::fill(points_weights.second.begin() + 4,
                        points_weights.second.end(),
                        0.5);
            }
          else
            std::fill(points_weights.second.begin(),
                      points_weights.second.end(),
                      1.0 / 8.0);
          break;
        case 3:
          {
            TriaIterator<TriaAccessor<3, 3, 3>> hex =
              static_cast<TriaIterator<TriaAccessor<3, 3, 3>>>(iterator);
            const unsigned int np = GeometryInfo<dim>::vertices_per_cell +
                                    GeometryInfo<dim>::lines_per_cell +
                                    GeometryInfo<dim>::faces_per_cell;
            Assert(points_weights.first.size() == np, ExcInternalError());
            Assert(points_weights.second.size() == np, ExcInternalError());
            auto *sp3 = reinterpret_cast<
              std::array<Point<3>, n_default_points_per_cell<decltype(hex)>()>
                *>(&points_weights.first);

            unsigned int j = 0;

            // note that the exact weights are chosen such as to minimize the
            // distortion of the eight new hexes from the optimal shape through
            // transfinite interpolation from the faces and vertices, see
            // TransfiniteInterpolationManifold for a deeper explanation of the
            // mechanisms
            if (with_interpolation)
              {
                for (unsigned int i = 0;
                     i < GeometryInfo<dim>::vertices_per_cell;
                     ++i, ++j)
                  {
                    (*sp3)[j]                = hex->vertex(i);
                    points_weights.second[j] = 1.0 / 8.0;
                  }
                for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell;
                     ++i, ++j)
                  {
                    (*sp3)[j] =
                      (hex->line(i)->has_children() ?
                         hex->line(i)->child(0)->vertex(1) :
                         hex->line(i)->get_manifold().get_new_point_on_line(
                           hex->line(i)));
                    points_weights.second[j] = -1.0 / 4.0;
                  }
                for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell;
                     ++i, ++j)
                  {
                    (*sp3)[j] =
                      (hex->quad(i)->has_children() ?
                         hex->quad(i)->isotropic_child(0)->vertex(3) :
                         hex->quad(i)->get_manifold().get_new_point_on_quad(
                           hex->quad(i)));
                    points_weights.second[j] = 1.0 / 2.0;
                  }
              }
            else
              // Overwrite the weights with 1/np if we don't want to use
              // interpolation.
              std::fill(points_weights.second.begin(),
                        points_weights.second.end(),
                        1.0 / np);
          }
          break;
        default:
          Assert(false, ExcInternalError());
          break;
      }
    return points_weights;
  }
} // namespace Manifolds

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


