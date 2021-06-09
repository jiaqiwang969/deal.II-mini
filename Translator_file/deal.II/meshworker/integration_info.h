//include/deal.II-translator/meshworker/integration_info_0.txt
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


#ifndef dealii_mesh_worker_integration_info_h
#define dealii_mesh_worker_integration_info_h

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/block_info.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/vector_selector.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * 用于交给本地集成函数的对象的类。
   * 这个类的对象包含一个或多个类型为FEValues、FEFaceValues或FESubfaceValues的对象，用于局部整合。它们被存储在一个指向基类FEValuesBase的数组中。模板参数VectorType允许在全局系统中使用不同的数据类型。
   * 此外，这个函数包含了存储在#global_data中的正交点的有限元函数值的空间。这些向量在每个单元或面都会自动初始化。为了避免初始化未使用的向量，你可以使用initialize_selector()来按名称选择你真正想要使用的向量。
   * <h3>Integration models</h3>
   * 该类支持两种局部集成模型，对应于Assembler命名空间的文档中的数据模型。一种是由使用FESystem建议的标准模型。也就是说，这个类中有一个FEValuesBase对象，包含了整个系统的所有形状函数，并且拥有和系统一样多的组件。使用这个模型需要对所有的系统形状函数进行循环。它要求识别每个形状函数的系统组件，并选择正确的双线性形式，通常在
   * @p if 或 @p switch 语句中。
   * 第二个集成模型为系统的每个基元建立一个FEValuesBase对象。每个单元上的自由度按块重新编号，这样它们就代表了与全局系统相同的块结构。然后，执行集成的对象可以单独处理每个块，这大大改善了代码的可重复使用性。
   * @note
   * 正如DoFInfo中描述的那样，在这个类中使用initialize()之前，通过调用
   * BlockInfo::initialize_local() 来触发对本地块模型的使用。
   * @ingroup MeshWorker
   *
   */
  template <int dim, int spacedim = dim>
  class IntegrationInfo
  {
  private:
    /// vector of FEValues objects
    std::vector<std::shared_ptr<FEValuesBase<dim, spacedim>>> fevalv;

  public:
    static const unsigned int dimension       = dim;
    static const unsigned int space_dimension = spacedim;

    /**
     * 构造函数。
     *
     */
    IntegrationInfo();

    /**
     * 复制构造函数，创建一个克隆，供 WorkStream::run().
     * 使用。
     *
     */
    IntegrationInfo(const IntegrationInfo<dim, spacedim> &other);

    /**
     * 建立所有的内部结构，特别是FEValuesBase对象，并为数据向量分配空间。
     * @param  el是DoFHandler的有限元素。          @param
     * mapping是用于映射网格单元的Mapping对象。          @param
     * quadrature是在FEVALUES对象的构造器中使用的正交公式。
     * @param
     * flags是在FEVALUES对象的构造函数中使用的UpdateFlags。
     * @param
     * local_block_info是PDE系统的一个可选参数。如果它提供了合理的数据，那么单元格上的自由度将被重新排序以反映系统的块结构。
     *
     */
    template <class FEVALUES>
    void
    initialize(const FiniteElement<dim, spacedim> &            el,
               const Mapping<dim, spacedim> &                  mapping,
               const Quadrature<FEVALUES::integral_dimension> &quadrature,
               const UpdateFlags                               flags,
               const BlockInfo *local_block_info = nullptr);

    /**
     * 初始化数据向量并缓存选择器。
     *
     */
    void
    initialize_data(const std::shared_ptr<VectorDataBase<dim, spacedim>> &data);

    /**
     * 删除由initialize()创建的数据。
     *
     */
    void
    clear();

    /**
     * 返回一个对用于初始化此对象的FiniteElement的引用。
     *
     */
    const FiniteElement<dim, spacedim> &
    finite_element() const;

    /// This is true if we are assembling for multigrid
    bool multigrid;
    /// Access to finite element
    /**
     * 如果对单个元素使用initialize()，这就是正在使用的访问函数（没有BlockInfo参数）。如果应用于一个元素的向量，它会抛出一个异常。
     *
     */
    const FEValuesBase<dim, spacedim> &
    fe_values() const;

    /// Access to finite elements
    /**
     * 如果对一组元素使用了initialize()（有一个有效的BlockInfo对象），必须使用这个访问函数。
     *
     */
    const FEValuesBase<dim, spacedim> &
    fe_values(const unsigned int i) const;

    /**
     * 包含正交点中的有限元函数值的向量。
     * 每个选定的有限元函数有一个向量，每个分量有一个向量，包含每个正交点的值的向量。
     *
     */
    std::vector<std::vector<std::vector<double>>> values;

    /**
     * 包含有限元函数在正交点的导数的向量。
     * 每个选定的有限元函数有一个向量，每个分量包含一个向量，包含每个正交点的值的向量。
     *
     */
    std::vector<std::vector<std::vector<Tensor<1, spacedim>>>> gradients;

    /**
     * 包含有限元函数在正交点的二阶导数的向量。
     * 每个选定的有限元函数有一个向量，每个分量包含一个向量，包含每个正交点的值的向量。
     *
     */
    std::vector<std::vector<std::vector<Tensor<2, spacedim>>>> hessians;

    /**
     * 重新初始化内部数据结构，以便在一个单元上使用。
     *
     */
    template <typename number>
    void
    reinit(const DoFInfo<dim, spacedim, number> &i);

    /**
     * 使用#global_data中的有限元函数，填充向量#values、#gradients和#hessians。
     *
     */
    template <typename number>
    void
    fill_local_data(const DoFInfo<dim, spacedim, number> &info,
                    bool                                  split_fevalues);

    /**
     * 用于计算正交点的函数值的全局数据向量。
     *
     */
    std::shared_ptr<VectorDataBase<dim, spacedim>> global_data;

    /**
     * 这个对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * 用于初始化的（系统）元素的指针。
     *
     */
    SmartPointer<const FiniteElement<dim, spacedim>,
                 IntegrationInfo<dim, spacedim>>
      fe_pointer;

    /**
     * 使用#global_data中的有限元函数，并根据选择器的值填充向量#values、#gradients和#hessians。
     *
     */
    template <typename TYPE>
    void
    fill_local_data(std::vector<std::vector<std::vector<TYPE>>> &data,
                    VectorSelector &                             selector,
                    bool split_fevalues) const;
    /**
     * 缓存系统元素的组件数量。
     *
     */
    unsigned int n_components;
  };

  /**
   * 持有用于整合单元格和面的划痕数据的对象。  IntegrationInfoBox有三个主要用途。      <ol>   <li>  它提供了 MeshWorker::loop(), 所需的接口，即两个函数post_cell()和post_faces()，以及数据成员#cell、#boundary、#face、#subface和#neighbor。      <li>  它包含了初始化 IntegrationInfo 数据成员中的 FEValues 和 FEFaceValues 对象所需的所有信息。      <li>  它存储了关于有限元向量的信息，以及它们的数据是否应该被用来计算正交点的函数值或导数。      <li>  它对正交规则和更新标志进行有根据的猜测，因此在默认参数足够的情况下，需要编写最少的代码。    </ol>  为了允许足够的通用性，必须采取几个步骤来使用这个类。    首先，你应该考虑你是否需要AnyData对象中任何向量的值。如果是的话，将VectorSelector对象#cell_selector、#boundary_selector和#face_selector的名称和要提取的数据类型（值、梯度、Hessian）填入。    之后，你将需要考虑FEValues对象的UpdateFlags。一个好的开始是initialize_update_flags()，它查看之前填写的选择器，并添加所有需要的标志以获得选择。  其他的标志可以用add_update_flags()来设置。    最后，我们需要选择正交公式。在最简单的情况下，你可能对默认设置感到满意，即<i>n</i>-点高斯公式。如果只使用形状函数的导数（#update_values没有设置），<i>n</i>等于FiniteElement中的最高多项式度数，如果#update_values设置了，<i>n</i>是比这个度数高一个。 如果你选择使用其他大小的高斯公式，请使用initialize_gauss_quadrature()的适当值。否则，你可以直接填写#cell_quadrature、#boundary_quadrature和#face_quadrature这些变量。    为了节省时间，你可以将基类的boundary_fluxes和internal_fluxes变量设置为false，从而告诉 Meshworker::loop() 不要在这些面上循环。    这里的所有信息都是用来正确设置IntegrationInfo对象的，通常是在一个IntegrationInfoBox中。
   * @ingroup MeshWorker
   *
   */
  template <int dim, int spacedim = dim>
  class IntegrationInfoBox
  {
  public:
    /**
     * 单元的 @p info 对象的类型。
     *
     */
    using CellInfo = IntegrationInfo<dim, spacedim>;

    /**
     * 默认构造函数。
     *
     */
    IntegrationInfoBox();

    /**
     * 初始化包含的IntegrationInfo对象。
     * 在这样做之前，添加必要的更新标志以产生所需的数据，同时将未初始化的正交规则设置为高斯公式，该公式精确地整合多项式双线性形式。
     *
     */
    void
    initialize(const FiniteElement<dim, spacedim> &el,
               const Mapping<dim, spacedim> &      mapping,
               const BlockInfo *                   block_info = nullptr);

    /**
     * 初始化包含的IntegrationInfo对象。
     * 在这样做之前，添加必要的更新标志以产生所需的数据，同时将未初始化的正交规则设置为高斯公式，该公式完全集成多项式双线性形式。
     *
     */
    template <typename VectorType>
    void
    initialize(const FiniteElement<dim, spacedim> &el,
               const Mapping<dim, spacedim> &      mapping,
               const AnyData &                     data,
               const VectorType &                  dummy,
               const BlockInfo *                   block_info = nullptr);
    /**
     * 初始化包含的IntegrationInfo对象。
     * 在这样做之前，添加必要的更新标志以产生所需的数据，同时将未初始化的正交规则设置为高斯公式，该公式完全集成多项式双线性形式。
     *
     */
    template <typename VectorType>
    void
    initialize(const FiniteElement<dim, spacedim> &el,
               const Mapping<dim, spacedim> &      mapping,
               const AnyData &                     data,
               const MGLevelObject<VectorType> &   dummy,
               const BlockInfo *                   block_info = nullptr);
    /**
     * @name  FEValues设置
     *
     */
     /* 
     * @{ */ 

    /**
     * 在initialize()之前调用这个函数，以便根据选择的数据猜测需要的更新标志。
     * 在计算面通量时，我们通常可以使用原始单元的几何图形（积分权重和法向量），从而可以避免在邻近单元上更新这些值。将<tt>neighbor_geometry</tt>设置为true，以便将这些值也初始化。
     *
     */
    void
    initialize_update_flags(bool neighbor_geometry = false);

    /**
     * 增加FEValues
     * UpdateFlags，以便在所有对象（单元格、边界面和所有内部面）上进行整合。
     *
     */
    void
    add_update_flags_all(const UpdateFlags flags);

    /**
     * 添加FEValues UpdateFlags用于单元格的整合。
     *
     */
    void
    add_update_flags_cell(const UpdateFlags flags);

    /**
     * 添加FEValues UpdateFlags用于边界面的积分。
     *
     */
    void
    add_update_flags_boundary(const UpdateFlags flags);

    /**
     * 添加FEValues UpdateFlags用于内部面的积分。
     *
     */
    void
    add_update_flags_face(const UpdateFlags flags);

    /**
     * 在本程序中已经设置的更新标志的基础上，增加额外的更新标志。
     * 四个布尔标志表示是否应该为单元、边界、单元本身或相邻单元的元素间面，或它们的任何组合设置额外的标志。
     *
     */
    void
    add_update_flags(const UpdateFlags flags,
                     const bool        cell     = true,
                     const bool        boundary = true,
                     const bool        face     = true,
                     const bool        neighbor = true);

    /**
     * 为每个正交规则分配n点高斯正交值。这里，0点的大小意味着不应该对这些网格实体进行循环。
     * 如果参数<tt>force</tt>为真，那么所有的正交集都被填入新的正交规则。如果它是假的，那么只有空的规则被改变。
     *
     */
    void
    initialize_gauss_quadrature(unsigned int n_cell_points,
                                unsigned int n_boundary_points,
                                unsigned int n_face_points,
                                const bool   force = true);

    /**
     * 这个对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 用于边界单元集成的更新标志集。
     * 默认为#update_JxW_values。
     *
     */
    UpdateFlags cell_flags;
    /**
     * 用于边界面集成的更新标志集。
     * 默认为#update_JxW_values和#update_normal_vectors。
     *
     */
    UpdateFlags boundary_flags;

    /**
     * 用于内部面集成的更新标志集。
     * 默认为#update_JxW_values和#update_normal_vectors。
     *
     */
    UpdateFlags face_flags;

    /**
     * 用于内部面集成的更新标志集。
     * 默认为#update_default，因为正交权重取自其他单元。
     *
     */
    UpdateFlags neighbor_flags;

    /**
     * 在单元格上使用的正交规则。
     *
     */
    Quadrature<dim> cell_quadrature;

    /**
     * 用于边界面的正交规则。
     *
     */
    Quadrature<dim - 1> boundary_quadrature;

    /**
     * 用于内部面的正交规则。
     *
     */
    Quadrature<dim - 1> face_quadrature;
     /* @} */ 

    /**
     * @name  数据向量
     *
     */
     /* 
     * @{ */ 

    /**
     * 初始化VectorSelector对象#cell_selector、#boundary_selector和#face_selector，以节省计算工作。如果没有使用选择器，那么
     * DoFInfo::global_data
     * 中所有命名的向量的值都将在所有正交点中计算。
     * 这个函数还将把UpdateFlags添加到存储在这个类中的标志中。
     *
     */
    /**
     * 从 DoFInfo::global_data
     * 中选择应该在单元格上的正交点计算的向量。
     *
     */
    MeshWorker::VectorSelector cell_selector;

    /**
     * 从 DoFInfo::global_data
     * 中选择应该在边界面上的正交点计算的向量。
     *
     */
    MeshWorker::VectorSelector boundary_selector;

    /**
     * 从 DoFInfo::global_data
     * 中选择应该在内部面的正交点上计算的向量。
     *
     */
    MeshWorker::VectorSelector face_selector;

    std::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim>> cell_data;
    std::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim>> boundary_data;
    std::shared_ptr<MeshWorker::VectorDataBase<dim, spacedim>> face_data;
     /* @} */ 

    /**
     * @name   MeshWorker::loop() 的接口。
     *
     */
     /* 
     * @{ */ 
    /**
     * 一个回调函数，它在所有单元格的循环中被调用，在对一个单元格的操作被执行后，在面的处理之前。
     * 为了使这个函数产生这种效果，至少loop()的参数<tt>boundary_worker</tt>或<tt>face_worker</tt>应该是非零的。此外，<tt>cells_first</tt>应该为真。如果<tt>cells_first</tt>是假的，那么在对一个单元格采取任何行动之前就会调用这个函数。
     * 在这个类中是空函数，但在其他类中可以用loop()代替。
     * 参见loop()和cell_action()以了解该函数的更多使用细节。
     *
     */
    template <class DOFINFO>
    void
    post_cell(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * 一个回调函数，它在所有单元格的循环中被调用，在对一个单元格的面进行操作后，在处理该单元格本身之前（假设<tt>cells_first</tt>为假）。
     * 为了使这个函数有合理的效果，至少loop()的参数<tt>boundary_worker</tt>或<tt>face_worker</tt>应该是非零的。此外，<tt>cells_first</tt>应该是假的。
     * 而且在这个类中是空函数，但在其他类中可以用loop()代替。
     * 参见loop()和cell_action()以了解该函数的详细使用方法。
     *
     */
    template <class DOFINFO>
    void
    post_faces(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * 一个单元格的 @p info 对象。
     *
     */
    CellInfo cell;
    /**
     * 边界面的 @p info 对象。
     *
     */
    CellInfo boundary;
    /**
     * 从第一个单元看到的 @p info
     * 对象为一个规则的内部面。
     *
     */
    CellInfo face;
    /**
     * 从第一个单元看到的内部面的 @p info 对象的精炼面。
     *
     */
    CellInfo subface;
    /**
     * 从另一个单元看到的内部面的 @p info 对象。
     *
     */
    CellInfo neighbor;

     /* @} */ 
  };


  //----------------------------------------------------------------------//

  template <int dim, int sdim>
  inline IntegrationInfo<dim, sdim>::IntegrationInfo()
    : fevalv(0)
    , multigrid(false)
    , global_data(std::make_shared<VectorDataBase<dim, sdim>>())
    , n_components(numbers::invalid_unsigned_int)
  {}


  template <int dim, int sdim>
  inline IntegrationInfo<dim, sdim>::IntegrationInfo(
    const IntegrationInfo<dim, sdim> &other)
    : multigrid(other.multigrid)
    , values(other.values)
    , gradients(other.gradients)
    , hessians(other.hessians)
    , global_data(other.global_data)
    , fe_pointer(other.fe_pointer)
    , n_components(other.n_components)
  {
    fevalv.resize(other.fevalv.size());
    for (unsigned int i = 0; i < other.fevalv.size(); ++i)
      {
        const FEValuesBase<dim, sdim> &p = *other.fevalv[i];
        const FEValues<dim, sdim> *    pc =
          dynamic_cast<const FEValues<dim, sdim> *>(&p);
        const FEFaceValues<dim, sdim> *pf =
          dynamic_cast<const FEFaceValues<dim, sdim> *>(&p);
        const FESubfaceValues<dim, sdim> *ps =
          dynamic_cast<const FESubfaceValues<dim, sdim> *>(&p);

        if (pc != nullptr)
          fevalv[i] =
            std::make_shared<FEValues<dim, sdim>>(pc->get_mapping(),
                                                  pc->get_fe(),
                                                  pc->get_quadrature(),
                                                  pc->get_update_flags());
        else if (pf != nullptr)
          fevalv[i] =
            std::make_shared<FEFaceValues<dim, sdim>>(pf->get_mapping(),
                                                      pf->get_fe(),
                                                      pf->get_quadrature(),
                                                      pf->get_update_flags());
        else if (ps != nullptr)
          fevalv[i] = std::make_shared<FESubfaceValues<dim, sdim>>(
            ps->get_mapping(),
            ps->get_fe(),
            ps->get_quadrature(),
            ps->get_update_flags());
        else
          Assert(false, ExcInternalError());
      }
  }



  template <int dim, int sdim>
  template <class FEVALUES>
  inline void
  IntegrationInfo<dim, sdim>::initialize(
    const FiniteElement<dim, sdim> &                el,
    const Mapping<dim, sdim> &                      mapping,
    const Quadrature<FEVALUES::integral_dimension> &quadrature,
    const UpdateFlags                               flags,
    const BlockInfo *                               block_info)
  {
    fe_pointer = &el;
    if (block_info == nullptr || block_info->local().size() == 0)
      {
        fevalv.resize(1);
        fevalv[0] = std::make_shared<FEVALUES>(mapping, el, quadrature, flags);
      }
    else
      {
        fevalv.resize(el.n_base_elements());
        for (unsigned int i = 0; i < fevalv.size(); ++i)
          fevalv[i] = std::make_shared<FEVALUES>(mapping,
                                                 el.base_element(i),
                                                 quadrature,
                                                 flags);
      }
    n_components = el.n_components();
  }


  template <int dim, int spacedim>
  inline const FiniteElement<dim, spacedim> &
  IntegrationInfo<dim, spacedim>::finite_element() const
  {
    Assert(fe_pointer != nullptr, ExcNotInitialized());
    return *fe_pointer;
  }

  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim> &
  IntegrationInfo<dim, spacedim>::fe_values() const
  {
    AssertDimension(fevalv.size(), 1);
    return *fevalv[0];
  }


  template <int dim, int spacedim>
  inline const FEValuesBase<dim, spacedim> &
  IntegrationInfo<dim, spacedim>::fe_values(unsigned int i) const
  {
    AssertIndexRange(i, fevalv.size());
    return *fevalv[i];
  }


  template <int dim, int spacedim>
  template <typename number>
  inline void
  IntegrationInfo<dim, spacedim>::reinit(
    const DoFInfo<dim, spacedim, number> &info)
  {
    for (unsigned int i = 0; i < fevalv.size(); ++i)
      {
        FEValuesBase<dim, spacedim> &febase = *fevalv[i];
        if (info.sub_number != numbers::invalid_unsigned_int)
          {
            // This is a subface
            FESubfaceValues<dim, spacedim> &fe =
              dynamic_cast<FESubfaceValues<dim, spacedim> &>(febase);
            fe.reinit(info.cell, info.face_number, info.sub_number);
          }
        else if (info.face_number != numbers::invalid_unsigned_int)
          {
            // This is a face
            FEFaceValues<dim, spacedim> &fe =
              dynamic_cast<FEFaceValues<dim, spacedim> &>(febase);
            fe.reinit(info.cell, info.face_number);
          }
        else
          {
            // This is a cell
            FEValues<dim, spacedim> &fe =
              dynamic_cast<FEValues<dim, spacedim> &>(febase);
            fe.reinit(info.cell);
          }
      }

    const bool split_fevalues = info.block_info != nullptr;
    if (!global_data->empty())
      fill_local_data(info, split_fevalues);
  }



  //----------------------------------------------------------------------//

  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::initialize_gauss_quadrature(unsigned int cp,
                                                             unsigned int bp,
                                                             unsigned int fp,
                                                             bool         force)
  {
    if (force || cell_quadrature.size() == 0)
      cell_quadrature = QGauss<dim>(cp);
    if (force || boundary_quadrature.size() == 0)
      boundary_quadrature = QGauss<dim - 1>(bp);
    if (force || face_quadrature.size() == 0)
      face_quadrature = QGauss<dim - 1>(fp);
  }


  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::add_update_flags_all(const UpdateFlags flags)
  {
    add_update_flags(flags, true, true, true, true);
  }


  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::add_update_flags_cell(const UpdateFlags flags)
  {
    add_update_flags(flags, true, false, false, false);
  }


  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::add_update_flags_boundary(
    const UpdateFlags flags)
  {
    add_update_flags(flags, false, true, false, false);
  }


  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::add_update_flags_face(const UpdateFlags flags)
  {
    add_update_flags(flags, false, false, true, true);
  }


  template <int dim, int sdim>
  inline void
  IntegrationInfoBox<dim, sdim>::initialize(const FiniteElement<dim, sdim> &el,
                                            const Mapping<dim, sdim> &mapping,
                                            const BlockInfo *block_info)
  {
    initialize_update_flags();
    initialize_gauss_quadrature((cell_flags & update_values) ?
                                  (el.tensor_degree() + 1) :
                                  el.tensor_degree(),
                                (boundary_flags & update_values) ?
                                  (el.tensor_degree() + 1) :
                                  el.tensor_degree(),
                                (face_flags & update_values) ?
                                  (el.tensor_degree() + 1) :
                                  el.tensor_degree(),
                                false);

    cell.template initialize<FEValues<dim, sdim>>(
      el, mapping, cell_quadrature, cell_flags, block_info);
    boundary.template initialize<FEFaceValues<dim, sdim>>(
      el, mapping, boundary_quadrature, boundary_flags, block_info);
    face.template initialize<FEFaceValues<dim, sdim>>(
      el, mapping, face_quadrature, face_flags, block_info);
    subface.template initialize<FESubfaceValues<dim, sdim>>(
      el, mapping, face_quadrature, face_flags, block_info);
    neighbor.template initialize<FEFaceValues<dim, sdim>>(
      el, mapping, face_quadrature, neighbor_flags, block_info);
  }


  template <int dim, int sdim>
  template <typename VectorType>
  void
  IntegrationInfoBox<dim, sdim>::initialize(const FiniteElement<dim, sdim> &el,
                                            const Mapping<dim, sdim> &mapping,
                                            const AnyData &           data,
                                            const VectorType &,
                                            const BlockInfo *block_info)
  {
    initialize(el, mapping, block_info);
    std::shared_ptr<VectorData<VectorType, dim, sdim>> p;
    VectorDataBase<dim, sdim> *                        pp;

    p = std::make_shared<VectorData<VectorType, dim, sdim>>(cell_selector);
    // Public member function of parent class was not found without
    // explicit cast
    pp = &*p;
    pp->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p  = std::make_shared<VectorData<VectorType, dim, sdim>>(boundary_selector);
    pp = &*p;
    pp->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p  = std::make_shared<VectorData<VectorType, dim, sdim>>(face_selector);
    pp = &*p;
    pp->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }

  template <int dim, int sdim>
  template <typename VectorType>
  void
  IntegrationInfoBox<dim, sdim>::initialize(const FiniteElement<dim, sdim> &el,
                                            const Mapping<dim, sdim> &mapping,
                                            const AnyData &           data,
                                            const MGLevelObject<VectorType> &,
                                            const BlockInfo *block_info)
  {
    initialize(el, mapping, block_info);
    std::shared_ptr<MGVectorData<VectorType, dim, sdim>> p;
    VectorDataBase<dim, sdim> *                          pp;

    p = std::make_shared<MGVectorData<VectorType, dim, sdim>>(cell_selector);
    // Public member function of parent class was not found without
    // explicit cast
    pp = &*p;
    pp->initialize(data);
    cell_data = p;
    cell.initialize_data(p);

    p =
      std::make_shared<MGVectorData<VectorType, dim, sdim>>(boundary_selector);
    pp = &*p;
    pp->initialize(data);
    boundary_data = p;
    boundary.initialize_data(p);

    p  = std::make_shared<MGVectorData<VectorType, dim, sdim>>(face_selector);
    pp = &*p;
    pp->initialize(data);
    face_data = p;
    face.initialize_data(p);
    subface.initialize_data(p);
    neighbor.initialize_data(p);
  }

  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim, sdim>::post_cell(const DoFInfoBox<dim, DOFINFO> &)
  {}


  template <int dim, int sdim>
  template <class DOFINFO>
  void
  IntegrationInfoBox<dim, sdim>::post_faces(const DoFInfoBox<dim, DOFINFO> &)
  {}


} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


