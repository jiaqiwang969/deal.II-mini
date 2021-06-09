//include/deal.II-translator/meshworker/local_integrator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


#ifndef dealii_mesh_worker_local_integrator_h
#define dealii_mesh_worker_local_integrator_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <functional>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  // Forward declarations
#ifndef DOXYGEN
  template <int dim, int spacedim, typename number>
  class DoFInfo;
  template <int dim, int spacedim>
  class IntegrationInfo;
#endif

  /**
   * 一个局部积分器对象，它可以用来简化loop()的调用。我们没有单独提供这三个局部积分函数，而是在这个类中把它们捆绑成虚拟函数。
   * 此外，由于我们不能有一个虚拟的空函数，我们提供了标志，允许我们指示是否要在边界和内部面进行积分。这些标志在默认情况下为真，但可以由应用程序修改，以加快循环速度。
   * 如果一个函数在派生类中没有被重载，但它的使用标志是真的，那么这个函数将导致一个异常ExcPureFunction。
   * @ingroup MeshWorker
   *
   */
  template <int dim, int spacedim = dim, typename number = double>
  class LocalIntegrator : public Subscriptor
  {
  public:
    /**
     * 构造函数设置默认值，即所有集成标志为真。
     *
     */
    LocalIntegrator();

    /**
     * 构造函数将集成标志设置为指定值。
     *
     */
    LocalIntegrator(bool use_cell, bool use_boundary, bool use_face);

    /**
     * 空的虚拟析构器。
     *
     */
    virtual ~LocalIntegrator() override = default;

    /**
     * 用于在单元格上集成的虚拟函数。如果没有被派生类重载，则抛出PureFunctionCalled异常。
     *
     */
    virtual void
    cell(DoFInfo<dim, spacedim, number> &dinfo,
         IntegrationInfo<dim, spacedim> &info) const;
    /**
     * 用于边界面积分的虚拟函数。如果没有被派生类重载，则抛出PureFunctionCalled异常。
     *
     */
    virtual void
    boundary(DoFInfo<dim, spacedim, number> &dinfo,
             IntegrationInfo<dim, spacedim> &info) const;
    /**
     * 用于对内部面进行积分的虚拟函数。如果没有被派生类重载，则抛出PureFunctionCalled异常。
     *
     */
    virtual void
    face(DoFInfo<dim, spacedim, number> &dinfo1,
         DoFInfo<dim, spacedim, number> &dinfo2,
         IntegrationInfo<dim, spacedim> &info1,
         IntegrationInfo<dim, spacedim> &info2) const;

    /**
     * 指示单元格积分器cell()是否要在循环中使用的标志。默认为<tt>true</tt>。
     *
     */
    bool use_cell;

    /**
     * 表示是否在循环中使用边界积分器boundary()的标志。默认为<tt>true</tt>。
     *
     */
    bool use_boundary;

    /**
     * 表示是否在循环中使用内部面积分器face()的标志。默认为<tt>true</tt>。
     *
     */
    bool use_face;

    /**
     * 输入向量的名称。如果这个向量不是空的，它可以被应用程序用来自动选择和验证用于积分的输入向量。
     * @note
     * 这个变量目前不被库使用，但提供它是为了帮助开发应用程序。
     *
     */
    std::vector<std::string> input_vector_names;

    /**
     * 产生的结果的名称。如果这个向量是非空的，它可以被应用程序用来为输出值自动分配名称和/或验证向量的名称。
     * @note
     * 这个变量目前不被库使用，但提供它是为了帮助开发应用程序。
     *
     */
    std::vector<std::string> output_names;

    /**
     * 如果虚拟函数cell()、boundary()或face()之一在派生类中没有被重载就被调用，就会产生这个错误。可以考虑将#use_cell、#use_boundary和#use_face分别设置为假。
     * @ingroup Exceptions
     *
     */
    DeclException0(ExcPureFunction);
  };
} // namespace MeshWorker



DEAL_II_NAMESPACE_CLOSE

#endif


