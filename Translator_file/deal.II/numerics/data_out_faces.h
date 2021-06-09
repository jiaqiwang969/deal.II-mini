//include/deal.II-translator/numerics/data_out_faces_0.txt
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

#ifndef dealii_data_out_faces_h
#define dealii_data_out_faces_h


#include <deal.II/base/config.h>

#include <deal.II/numerics/data_out.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace DataOutFacesImplementation
  {
    /**
     * 一个派生类，用于DataOutFaces类中。这是一个用于WorkStream上下文的文档中讨论的AdditionalData那种数据结构的类。
     *
     */
    template <int dim, int spacedim>
    struct ParallelData
      : public internal::DataOutImplementation::ParallelDataBase<dim, spacedim>
    {
      ParallelData(const unsigned int               n_datasets,
                   const unsigned int               n_subdivisions,
                   const std::vector<unsigned int> &n_postprocessor_outputs,
                   const Mapping<dim, spacedim> &   mapping,
                   const std::vector<
                     std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
                     &               finite_elements,
                   const UpdateFlags update_flags);

      std::vector<Point<spacedim>> patch_evaluation_points;
    };
  } // namespace DataOutFacesImplementation
} // namespace internal


/**
 * 这个类从三角形的面生成输出。它可能只用于生成三角形表面的输出（这是该类的默认值），或者用于活动单元的所有面，如构造函数中指定的那样。这个类的输出是一组补丁（如类
 * DataOutBase::Patch()),
 * 所定义的，每个面都要生成输出。这些补丁可以通过基础类的功能写成几种图形数据格式。
 * <h3>Interface</h3>
 * 这个类的接口是从DataOut类复制过来的。此外，它们共享共同的父类DataOut_DoFData。关于接口的讨论，请看这两个类的参考资料。
 *
 *  <h3>Extending this class</h3>
 * 生成补丁的面的序列是以与DataOut类相同的方式生成的；关于各自接口的描述见那里。产生面片序列的函数将被用于生成输出，这些函数被称为first_face()和next_face()。
 * 由于我们需要用这些函数生成的面来初始化FEValues类型的对象，所以它们只返回面的迭代器是不够的。相反，我们需要一对单元格和面的编号，因为有限元场的值在一个面上不一定是唯一的（想想不连续的有限元，有限元场的值取决于你接近一个面的方向，因此有必要使用一对单元格和面，而不是只有一个面的迭代器）。因此，这个类定义了一个
 * @p alias ，它创建了一个类型 @p FaceDescriptor
 * ，是一对单元格迭代器和面数的缩写。函数 @p first_face 和
 * @p next_face 对该类型的对象进行操作。
 * 例如，如果你只想从边界的某些部分输出，例如由各自面的边界指示器指示的，那么扩展这个类可能是有用的。然而，我们也可以想象，我们不是从边界面生成补丁，而是从根据其他标准选择的内部面生成补丁；一种应用可能是只使用那些解的一个分量达到某个值的面，以便在这些面上显示其他解分量的值。其他的应用当然也存在，对于这些应用，作者没有足够的想象力。
 * @todo  使用实际的FEFaceValues和MeshWorker重新实现这整个类。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim, int spacedim = dim>
class DataOutFaces : public DataOut_DoFData<dim, dim - 1, spacedim, spacedim>
{
  static_assert(dim == spacedim, "Not implemented for dim != spacedim.");

public:
  /**
   * 补丁的尺寸参数。
   *
   */
  static constexpr int patch_dim      = dim - 1;
  static constexpr int patch_spacedim = spacedim;

  /**
   * 对所考虑的dof处理程序类的迭代器类型的别名。
   *
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, patch_dim, spacedim, patch_spacedim>::
      cell_iterator;

  /**
   * 决定是否写入表面网（默认）或整个线筐的构造函数。
   *
   */
  DataOutFaces(const bool surface_only = true);

  /**
   * 这是该类的核心功能，因为它建立了由基类的低级函数写入的补丁列表。从本质上讲，补丁是三角形和DoFHandler对象的每个面上的数据的一些中间表示，然后可以用来以某种可视化程序可读的格式编写文件。
   * 你可以在这个类的一般文档中找到关于这个函数的使用概述。在这个类的基类DataOut_DoFData的文档中也提供了一个例子。
   * @param  n_subdivisions 参见 DataOut::build_patches()
   * 关于这个参数的广泛描述。
   *
   */
  virtual void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * 与上述相同，只是额外的第一个参数定义了一个映射，该映射将用于生成输出。如果<tt>n_subdivisions>1</tt>，源于域的边界单元的细分斑块内部的点可以用映射来计算，也就是说，高阶映射导致通过使用更多的细分来表示弯曲的边界。
   * 即使对于非弯曲的单元，映射参数也可以用于欧拉映射（见MappingQ1Eulerian类），在这里，映射不仅用于确定单元内部的点的位置，也用于确定顶点的位置。
   * 它提供了一个机会，可以在实际进行计算的变形三角形上观察解决方案，即使网格在内部是以未变形的配置存储的，变形只是通过一个额外的向量来跟踪每个顶点的变形情况。
   * @todo  如果是具有hp-capabilities的DoFHandler， @p mapping
   * 参数应该被 hp::MappingCollection 取代。
   *
   */
  virtual void
  build_patches(const Mapping<dim, spacedim> &mapping,
                const unsigned int            n_subdivisions = 0);

  /**
   * 声明一种描述我们希望产生输出的面的方式。通常的方法当然是使用
   * <tt>DoFHandler<dim>::face_iterator</tt>,
   * 类型的对象，但由于我们必须对FEValues类型的对象描述面孔，我们只能用单元格和面孔的编号对来表示面孔。这一对在这里被别名为一个名字，最好是输入。
   *
   */
  using FaceDescriptor = typename std::pair<cell_iterator, unsigned int>;


  /**
   * 返回我们想要输出的第一个面。默认的实现是返回一个（本地拥有的）活动单元的第一个面，或者，如果在析构器中设置了
   * @p surface_only
   * 选项（如默认），则返回位于边界上的第一个此类面。
   * 如果你想使用不同的逻辑来决定哪些面应该有助于图形输出的创建，你可以在派生类中重载这个函数。
   *
   */
  virtual FaceDescriptor
  first_face();

  /**
   * 返回我们想要输出的下一个面。如果没有更多的面，<tt>dofs->end()</tt>将作为返回值的第一个组成部分返回。
   * 默认的实现是返回一个（本地拥有的）活动单元的下一个面，或者边界上的下一个这样的面（取决于是否向构造函数提供了
   * @p surface_only 选项）。
   * 这个函数逐个遍历网格中的活动单元（仅限于本地拥有的单元），然后遍历该单元的所有面。因此，内部面会被输出两次，这个特性对于不连续的Galerkin方法或使用数据后处理器可能产生的结果在单元之间不连续的情况下是很有用的）。)
   * 这个函数可以在派生类中被重载，以选择不同的面组。请注意，默认实现假设给定的
   * @p face
   * 是有效的，只要first_face()也是从默认实现中使用，就可以保证。只重载这两个函数中的一个应该谨慎行事。
   *
   */
  virtual FaceDescriptor
  next_face(const FaceDescriptor &face);

private:
  /**
   * 在表面网格和全线篮之间决定的参数。
   *
   */
  const bool surface_only;

  /**
   * 建立一个补丁。这个函数是在WorkStream上下文中调用的。
   *
   */
  void
  build_one_patch(
    const FaceDescriptor *cell_and_face,
    internal::DataOutFacesImplementation::ParallelData<dim, spacedim> &data,
    DataOutBase::Patch<patch_dim, patch_spacedim> &                    patch);
};

namespace Legacy
{
  /**
   * @deprecated  使用没有DoFHandlerType模板的 dealii::DataOutFaces
   * 代替。
   *
   */
  template <int dim, typename DoFHandlerType = DoFHandler<dim>>
  using DataOutFaces DEAL_II_DEPRECATED =
    dealii::DataOutFaces<dim, DoFHandlerType::space_dimension>;
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


