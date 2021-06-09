//include/deal.II-translator/opencascade/manifold_lib_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


#ifndef dealii_occ_manifold_lib_h
#define dealii_occ_manifold_lib_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_OPENCASCADE

#  include <deal.II/grid/manifold.h>

#  include <deal.II/opencascade/utilities.h>

// opencascade needs "HAVE_CONFIG_H" to be exported...
#  define HAVE_CONFIG_H
#  include <Adaptor3d_Curve.hxx>
#  include <Adaptor3d_HCurve.hxx>
#  include <BRepAdaptor_Curve.hxx>
#  undef HAVE_CONFIG_H

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup  OpenCASCADE  
     * @{ 
 *
 */

namespace OpenCASCADE
{
  /**
   * 一个基于OpenCASCADE
   * TopoDS_Shape的Manifold对象，其中新的点首先通过与FlatManifold相同的方式对周围的点进行平均计算，然后使用OpenCASCADE工具在法线方向进行投影。
   * 这个类对你传递给它的形状不做任何假设，而Manifold的拓扑维度是由TopoDS_Shape本身推断出来的。在调试模式下，有一个理智检查，以确保周围的点（在project_to_manifold()中使用的点）确实存在于Manifold上，也就是说，在这些点上调用
   * OpenCASCADE::closest_point()
   * 会使它们不被触动。如果不是这种情况，就会抛出ExcPointNotOnManifold。
   * 例如，如果你试图使用TopoDS_Edge类型的形状在一个面上投影时，这种情况可能发生。在这种情况下，面的顶点会被折叠到边缘上，而你周围的点不会位于给定的形状上，从而引发一个异常。
   *
   */
  template <int dim, int spacedim>
  class NormalProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * 标准构造函数接受一个通用的TopoDS_Shape  @p sh,
     * 和一个用于内部计算距离的公差。
     * TopoDS_Shape可以是任意的，即一个形状、面、边的集合，或一个单一的面或边。
     *
     */
    NormalProjectionManifold(const TopoDS_Shape &sh,
                             const double        tolerance = 1e-7);

    /**
     * 克隆当前的Manifold。
     *
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * 执行对流形的实际投影。该函数在调试模式下，检查每个
     * @p surrounding_points
     * 是否在给定的TopoDS_Shape的公差范围内。如果不是这样，就会抛出一个异常。
     * 投射点是使用OpenCASCADE的正常投影算法计算的。
     *
     */
    virtual Point<spacedim>
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim> &                 candidate) const override;


  protected:
    /**
     * 内部用于投影点的拓扑形状。你可以通过调用
     * OpenCASCADE::read_IGES()
     * 函数来构建这样的形状，它将用IGES文件中包含的几何图形创建一个TopoDS_Shape。
     *
     */
    const TopoDS_Shape sh;

    /**
     * 本类用于计算距离的相对公差。
     *
     */
    const double tolerance;
  };

  /**
   * 一个基于OpenCASCADE
   * TopoDS_Shape的流形对象，其中新的点首先通过与FlatManifold相同的方式对周围的点进行平均计算，然后使用OpenCASCADE工具沿着建造时指定的方向将它们投影到流形上。
   * 这个类对你传递给它的形状不做任何假设，而流形的拓扑维度是由TopoDS_Shape本身推断出来的。在调试模式下，有一个理智检查，以确保周围的点（在project_to_manifold()中使用的点）确实存在于Manifold上，也就是说，在这些点上调用
   * OpenCASCADE::closest_point()
   * 会使它们不被触动。如果不是这样的话，就会抛出ExcPointNotOnManifold。
   * 请注意，如果要精化的三角形接近给定的TopoDS_Shape的边界，或者当你在构造时使用的方向不与形状相交时，这种类型的Manifold描述符可能无法得到结果。当这种情况发生时，会抛出一个异常。
   *
   */
  template <int dim, int spacedim>
  class DirectionalProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * 构建一个Manifold对象，该对象将沿着给定的 @p direction.
     * 在TopoDS_Shape  @p sh, 上投影点。
     *
     */
    DirectionalProjectionManifold(const TopoDS_Shape &       sh,
                                  const Tensor<1, spacedim> &direction,
                                  const double               tolerance = 1e-7);

    /**
     * 克隆当前的Manifold。
     *
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * 执行对流形的实际投影。该函数在调试模式下，检查每个
     * @p surrounding_points
     * 是否在给定的TopoDS_Shape的公差范围内。如果不是这样，就会抛出一个异常。
     * 投射点是使用OpenCASCADE的方向性投影算法计算的。
     *
     */
    virtual Point<spacedim>
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim> &                 candidate) const override;

  protected:
    /**
     * 内部用于投影点的拓扑形状。你可以通过调用
     * OpenCASCADE::read_IGES()
     * 函数来构建这样的形状，它将用IGES文件中包含的几何图形创建一个TopoDS_Shape。
     *
     */
    const TopoDS_Shape sh;

    /**
     * 用来在形状上投射新点的方向。
     *
     */
    const Tensor<1, spacedim> direction;

    /**
     * 该类用于计算距离的相对公差。
     *
     */
    const double tolerance;
  };


  /**
   * 一个基于OpenCASCADE
   * TopoDS_Shape的流形对象，其中新的点首先通过与FlatManifold相同的方式对周围的点进行平均计算，然后使用OpenCASCADE工具沿着一个方向将它们投影到流形上，这个方向是对周围点（因此是网格单元）法线的估计。
   * 网格的法线方向是特别有用的，因为它是网格中缺少节点的方向。例如，在细化一个单元的过程中，最初会在该单元的baricenter周围创建一个新的节点。这个位置在某种程度上保证了与旧单元的节点有一个统一的距离。沿着原单元的法线方向将这样的单元棒心投射到CAD表面，就可以保持与原单元的点的统一距离。当然，在网格生成阶段，没有定义dof处理程序和有限元，这样的方向必须被估计。对于存在8个周围点的情况，4个不同的三角形被确定为指定的点，这些三角形的法线被平均化以获得对单元的法线的近似值。
   * 当然，存在2个周围点的情况（即：一个单元格边缘被细化）就比较麻烦了。首先计算2个周围点的CAD表面法线的平均值，然后投射到连接周围点的那段法线上。这也是为了让新的点与周围的点有相等的距离。这个类只对CAD面进行操作，并假设你传递给它的形状至少包含一个面。如果不是这样的话，就会抛出一个异常。在调试模式下，有一个理智的检查，以确保周围的点（在project_to_manifold()中使用的点）确实存在于Manifold上，也就是说，在这些点上调用
   * OpenCASCADE::closest_point()
   * 会使它们不被触动。如果不是这样的话，就会抛出ExcPointNotOnManifold。
   * 请注意，如果要精化的三角形接近给定的TopoDS_Shape的边界，或者从周围的点估计的法线方向不与形状相交，这种类型的Manifold描述符可能无法提供结果。
   * 当这种情况发生时，会抛出一个异常。
   *
   */
  template <int dim, int spacedim>
  class NormalToMeshProjectionManifold : public FlatManifold<dim, spacedim>
  {
  public:
    /**
     * 构建一个Manifold对象，将TopoDS_Shape  @p sh,
     * 上的点沿着与网格单元近似的法线方向投影。
     *
     */
    NormalToMeshProjectionManifold(const TopoDS_Shape &sh,
                                   const double        tolerance = 1e-7);

    /**
     * 克隆当前的Manifold。
     *
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * 执行对流形的实际投影。这个函数在调试模式下，检查每个
     * @p surrounding_points
     * 是否在给定的TopoDS_Shape的公差范围内。如果不是这样，就会抛出一个异常。
     *
     */
    virtual Point<spacedim>
    project_to_manifold(
      const ArrayView<const Point<spacedim>> &surrounding_points,
      const Point<spacedim> &                 candidate) const override;

  protected:
    /**
     * 内部用于投影点的拓扑形状。你可以通过调用
     * OpenCASCADE::read_IGES()
     * 函数来构建这样的形状，该函数将用IGES文件中包含的几何图形创建一个TopoDS_Shape。
     *
     */
    const TopoDS_Shape sh;

    /**
     * 该类用于计算距离的相对公差。
     *
     */
    const double tolerance;
  };

  /**
   * 一个基于OpenCASCADE
   * TopoDS_Shape对象的Manifold对象，其拓扑维度等于1（TopoDS_Edge或TopoDS_Wire），新点位于周围点的arclength平均值。如果给定的TopoDS_Shape可以被铸造为周期性（封闭）曲线，那么这个信息将被内部用来设置基础ChartManifold类的周期性。
   * 这个类只能在TopoDS_Edge或TopoDS_Wire对象上工作，而且只有当spacedim为3时才有意义。如果你使用了一个拓扑维度为1的对象，就会产生一个异常。
   * 在调试模式下，有一个额外的理智检查，以确保周围的点确实存在于Manifold上，也就是说，在这些点上调用
   * OpenCASCADE::closest_point()
   * 会让它们不被触动。如果不是这样的话，就会抛出一个ExcPointNotOnManifold。
   *
   */
  template <int dim, int spacedim>
  class ArclengthProjectionLineManifold : public ChartManifold<dim, spacedim, 1>
  {
  public:
    /**
     * 带有TopoDS_Edge的默认构造函数。
     *
     */
    ArclengthProjectionLineManifold(const TopoDS_Shape &sh,
                                    const double        tolerance = 1e-7);

    /**
     * 克隆当前的Manifold。
     *
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * 给出实空间上的一个点，找到其arclength参数。如果该点不在构造时给定的TopoDS_Edge上，在调试模式下会抛出一个错误。
     *
     */
    virtual Point<1>
    pull_back(const Point<spacedim> &space_point) const override;

    /**
     * 给定一个arclength参数，找到它在实空间的图像。
     *
     */
    virtual Point<spacedim>
    push_forward(const Point<1> &chart_point) const override;

  protected:
    /**
     * 用于构建此物体的实际形状。
     *
     */
    const TopoDS_Shape sh;

    /**
     * 一个曲线适配器。这是在计算中使用的一个，它指向上面的右边。
     *
     */
    Handle_Adaptor3d_HCurve curve;

    /**
     * 在所有内部计算中使用的相对公差。
     *
     */
    const double tolerance;

    /**
     * 曲线的总长度。如果边缘是周期性的，这也会被用作周期。
     *
     */
    const double length;
  };

  /**
   * 使用OpenCASCADE导入的CAD的面的歧管描述。
   * @ingroup manifold
   *
   */
  template <int dim, int spacedim>
  class NURBSPatchManifold : public ChartManifold<dim, spacedim, 2>
  {
  public:
    /**
     * 构造函数接收一个OpenCASCADE TopoDS_Face  @p face
     * 和一个可选的  @p tolerance.
     * 该类使用区间OpenCASCADE变量u, v来描述流形。
     *
     */
    NURBSPatchManifold(const TopoDS_Face &face, const double tolerance = 1e-7);

    /**
     * 克隆当前的流形。
     *
     */
    virtual std::unique_ptr<Manifold<dim, spacedim>>
    clone() const override;

    /**
     * 从欧几里得空间拉回给定的点。将返回与该点相关的uv坐标
     * @p space_point.  。
     *
     */
    virtual Point<2>
    pull_back(const Point<spacedim> &space_point) const override;

    /**
     * 在uv坐标系中给定一个 @p chart_point
     * ，该方法返回与之相关的欧几里得坐标。
     *
     */
    virtual Point<spacedim>
    push_forward(const Point<2> &chart_point) const override;

    /**
     * 给定间隔维欧氏空间中的一个点，该方法返回从uv坐标系映射到欧氏坐标系的函数
     * $F$ 的导数。换句话说，它是一个大小为
     * $\text{spacedim}\times\text{chartdim}$ 的矩阵。
     * 这个函数被用于get_tangent_vector()函数所要求的计算中。
     * 更多信息请参考该类的一般文档。
     *
     */
    virtual DerivativeForm<1, 2, spacedim>
    push_forward_gradient(const Point<2> &chart_point) const override;

  protected:
    /**
     * 返回一个代表u和v的最小值和最大值的元组。准确地说，它返回（u_min,
     * u_max, v_min, v_max)
     *
     */
    std::tuple<double, double, double, double>
    get_uv_bounds() const;

    /**
     * 一个OpenCASCADE TopoDS_Face  @p face  由CAD给出。
     *
     */
    TopoDS_Face face;

    /**
     * OpenCASCADE在每次操作中用来识别点的公差。
     *
     */
    double tolerance;
  };

} // namespace OpenCASCADE

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_OPENCASCADE
#endif // dealii_occ_manifold_lib_h


