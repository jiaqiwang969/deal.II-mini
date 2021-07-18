//include/deal.II-translator/gmsh/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


#ifndef dealii_gmsh_parameters_h
#define dealii_gmsh_parameters_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GMSH

#  ifdef DEAL_II_WITH_OPENCASCADE
#    include <TopoDS_CompSolid.hxx>
#    include <TopoDS_Compound.hxx>
#    include <TopoDS_Edge.hxx>
#    include <TopoDS_Face.hxx>
#    include <TopoDS_Shape.hxx>
#    include <TopoDS_Vertex.hxx>
#  endif

#  include <deal.II/base/parameter_handler.h>

#  include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个%Gmsh相关的实用程序和类的集合。
 *
 *
 */
namespace Gmsh
{
  /**
   * 一个参数类，用于向%Gmsh可执行程序传递选项。
   *
   */
  class AdditionalParameters
  {
  public:
    /**
     * 将所有的附加参数设置为它们的默认值。
     *
     */
    AdditionalParameters(const double       characteristic_length = 1.0,
                         const std::string &output_base_name      = "");

    /**
     * 为AdditionalParameters类的每个成员调用prm.add_parameter。
     *
     */
    void
    add_parameters(ParameterHandler &prm);

    /**
     * 用于定义%Gmsh网格的特征长度。
     * %Gmsh会尽量确保每条边的大小与这个值相近。
     *
     */
    double characteristic_length = 1.0;

    /**
     * 输出文件的基本名称。
     * 如果这个名字留空，那么就会使用临时文件，当不再需要时就会删除。
     *
     */
    std::string output_base_name = "";
  };

#  ifdef DEAL_II_WITH_OPENCASCADE
  /**
   * 给出一条平滑的闭合曲线，用 %Gmsh 创建一个三角形。
   * 输入的曲线  @p boundary  应该是封闭的。
   *
   */
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(
    const TopoDS_Edge &         boundary,
    Triangulation<2, spacedim> &tria,
    const AdditionalParameters &prm = AdditionalParameters());
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


