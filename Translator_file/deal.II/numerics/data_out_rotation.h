//include/deal.II-translator/numerics/data_out_rotation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_data_out_rotation_h
#define dealii_data_out_rotation_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutRotationImplementation
  {
    /**
     * 一个派生类，用于DataOutFaces类中。这是一个用于WorkStream类的文档中讨论的AdditionalData那种数据结构的类。
     *
     */
    template <int dim, int spacedim>
    struct ParallelData
      : public internal::DataOutImplementation::ParallelDataBase<dim, spacedim>
    {
      ParallelData(const unsigned int               n_datasets,
                   const unsigned int               n_subdivisions,
                   const unsigned int               n_patches_per_circle,
                   const std::vector<unsigned int> &n_postprocessor_outputs,
                   const Mapping<dim, spacedim> &   mapping,
                   const std::vector<
                     std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                     &               finite_elements,
                   const UpdateFlags update_flags);

      const unsigned int n_patches_per_circle;

      std::vector<Point<spacedim>> patch_evaluation_points;
    };
  } // namespace DataOutRotationImplementation
} // namespace internal



/**
 * 这个类在全域的计算中生成输出，这些计算是利用域和解的旋转对称性完成的。特别是，如果一个三维问题的计算是围绕
 * @p z-axis 旋转对称的（即在 @p r-z-plane)
 * 中完成的，那么这个类可以用来生成原始 @p x-y-z
 * 空间中的输出。为了做到这一点，它从计算网格中的每个单元生成一个空间中的单元，其尺寸比DoFHandler对象的尺寸大。这样的输出将由六面体组成，形成一个围绕Z轴旋转对称的对象。由于大多数图形程序不能表示环状结构，角度（旋转）变量也被离散成有限数量的区间；这些区间的数量必须交给
 * @p build_patches
 * 函数。然而，需要注意的是，虽然这个函数生成了整个领域的漂亮图片，但它经常产生
 * <em> 非常 </em> 大的输出文件。
 *
 *  <h3>Interface</h3>
 * 这个类的接口是从DataOut类复制过来的。此外，它们共享共同的父类DataOut_DoFData()。关于接口的讨论以及如何通过从这个类派生出更多的类来扩展它，请参见这两个类的参考文献。
 *
 *  <h3>Details for 1d computations</h3>
 * 传递给这个类的DoFHandler对象所使用的三角测量中的一个坐标被当作径向变量，然后输出将是一个圆或一个环域。用户有责任保证径向坐标只达到非负值。
 *
 *  <h3>Details for 2d computations</h3>
 * 我们认为计算（由附加到该类的DoFHandler对象表示）发生在
 * @p r-z-plane, 中，其中 @p r 是径向变量， @p z
 * 表示围绕解决方案对称的旋转轴。输出是在 @p x-y-z
 * 空间，其中径向依赖被转换到 @p x-y
 * 平面。目前，不可能交换模拟所在平面的第一个和第二个变量的含义，即生成模拟的输出，其中第一个变量表示对称轴，第二个表示径向变量。在第一次为你的应用程序编程时，你必须考虑到这一点。
 * 用户有责任确保径向变量只达到非负值。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim, int spacedim = dim>
class DataOutRotation
  : public DataOut_DoFData<dim, dim + 1, spacedim, spacedim + 1>
{
  static_assert(dim == spacedim, "Not implemented for dim != spacedim.");

public:
  /**
   * 补丁的尺寸参数。
   *
   */
  static constexpr int patch_dim      = dim + 1;
  static constexpr int patch_spacedim = spacedim + 1;

  /**
   * 对所考虑的dof处理程序类的迭代器类型的类型定义。
   *
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
      cell_iterator;

  /**
   * 这是该类的核心功能，因为它建立了由基类的低级函数编写的补丁列表。从本质上讲，补丁是三角形和DoFHandler对象的每个单元上的数据的一些中间表示，然后可以用来以某种可视化程序可读的格式编写文件。
   * 你可以在这个类的一般文档中找到关于这个函数的使用概述。在这个类的基类DataOut_DoFData的文档中也提供了一个例子。
   * @param  n_patches_per_circle
   * 表示角度（旋转）变量将被细分为多少个间隔。
   * @param  n_subdivisions 参见 DataOut::build_patches()
   * 对该参数的详细描述。
   *
   */
  virtual void
  build_patches(const unsigned int n_patches_per_circle,
                const unsigned int n_subdivisions = 0);

  /**
   * 返回我们想要输出的第一个单元格。默认实现返回第一个 @ref GlossActive "活动单元"
   * ，但你可能想在派生类中返回其他单元。
   *
   */
  virtual cell_iterator
  first_cell();

  /**
   * 返回 @p cell
   * 之后我们想要输出的下一个单元格。如果没有更多的单元格，应返回<tt>dofs->end()</tt>。
   * 默认实现返回下一个活动单元，但你可能想在派生类中返回其他单元。请注意，默认实现假设给定的
   * @p cell 是活动的，只要 @p first_cell
   * 也被用于默认实现，就可以保证。只重载这两个函数中的一个可能不是个好主意。
   *
   */
  virtual cell_iterator
  next_cell(const cell_iterator &cell);

  /**
   * 异常情况
   *
   */
  DeclException1(ExcRadialVariableHasNegativeValues,
                 double,
                 << "You are attempting to use this class on a triangulation "
                    "in which some vertices have a negative radial coordinate "
                    "value of "
                 << arg1
                 << ". If you rotate such a triangulation around an "
                    "axis, you will get (dim+1)-dimensional meshes "
                    "that are not likely what you hoped to see.");

private:
  /**
   * 建立与第一个参数中给出的单元格相对应的所有补丁。在WorkStream中使用第二个参数作为并行调用的抓取空间，并将结果放到最后一个参数中。
   *
   */
  void
  build_one_patch(
    const cell_iterator *                                                 cell,
    internal::DataOutRotationImplementation::ParallelData<dim, spacedim> &data,
    std::vector<DataOutBase::Patch<patch_dim, patch_spacedim>> &my_patches);
};

namespace Legacy
{
  /**
   * @deprecated  使用没有DoFHandlerType模板的 dealii::DataOutRotation
   * 代替。
   *
   */
  template <int dim, typename DoFHandlerType = DoFHandler<dim>>
  using DataOutRotation DEAL_II_DEPRECATED =
    dealii::DataOutRotation<dim, DoFHandlerType::space_dimension>;
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


