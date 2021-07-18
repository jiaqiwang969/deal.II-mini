//include/deal.II-translator/meshworker/dof_info_0.txt
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


#ifndef dealii_mesh_worker_dof_info_h
#define dealii_mesh_worker_dof_info_h

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/block_info.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/vector_selector.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  // Forward declaration
#ifndef DOXYGEN
  template <int dim, class DOFINFO>
  class DoFInfoBox;
#endif

  /**
   * 一个包含网格对象的几何和自由度信息的类。
   * 这些对象中的信息通常会被某个Assembler类所使用。这也是基于网格的矩阵（通常被称为无矩阵方法）中需要的信息。
   * 除了存储在该类中的自由度信息外，它还为在LocalResults中操作的工作对象提供本地计算空间。这个基类将在每个单元上自动被重新初始化，但初始设置由用户决定，应该在调用这个类的initialize()时完成。
   * 该类在两种不同的模式下运行，与Assembler命名空间文档中讨论的数据模型相对应。
   * 本地数据模型的选择由向量 BlockInfo::local_renumbering,
   * 触发，而向量通常由 BlockInfo::initialize_local().
   * 填充。如果使用了这个函数，或者向量已经从零长度改变，那么存储在这个对象中的本地dof指数将自动重新编号以反映本地块结构。这意味着，
   * @p indices
   * 中的第一个条目将指的是系统的第一个块，然后是第二个块，依此类推。
   * BlockInfo对象是作为一个指针来存储的。因此，如果块结构发生变化，例如因为网格细化，DoFInfo类将自动使用新的结构。
   * @ingroup MeshWorker
   *
   */
  template <int dim, int spacedim = dim, typename number = double>
  class DoFInfo : public LocalResults<number>
  {
  public:
    /// The current cell
    typename Triangulation<dim, spacedim>::cell_iterator cell;

    /// The current face
    typename Triangulation<dim, spacedim>::face_iterator face;

    /**
     * 当前单元上的当前面的编号。        如果 @p info
     * 对象被初始化为一个单元格，那么这个数字就是
     * numbers::invalid_unsigned_int 。
     *
     */
    unsigned int face_number;

    /**
     * 当前面上的子面的编号 如果 @p info
     * 对象没有被初始化为子面，这个编号是
     * numbers::invalid_unsigned_int 。
     *
     */
    unsigned int sub_number;

    /**
     * 当前单元的DoF指数
     *
     */
    std::vector<types::global_dof_index> indices;

    /**
     * 当前单元上的DoF指数，按本地块组织。这个向量的大小是零，除非使用局部块。
     *
     */
    std::vector<std::vector<types::global_dof_index>> indices_by_block;

    /**
     * 设置#block_info指针的构造函数。
     *
     */
    DoFInfo(const BlockInfo &block_info);

    /**
     * 构造函数让#block_info指针为空，但设置#aux_local_indices。
     *
     */
    DoFInfo(const DoFHandler<dim, spacedim> &dof_handler);

    /**
     * 设置当前单元格并填充 @p indices. 。
     *
     */
    template <class DHCellIterator>
    void
    reinit(const DHCellIterator &c);

    /**
     * 如果#单元格发生变化，设置当前面和填充 @p indices 。
     *
     */
    template <class DHCellIterator, class DHFaceIterator>
    void
    reinit(const DHCellIterator &c,
           const DHFaceIterator &f,
           const unsigned int    face_no);

    /**
     * 设置当前的子面，如果#单元格发生变化，则填充 @p
     * indices 。
     *
     */
    template <class DHCellIterator, class DHFaceIterator>
    void
    reinit(const DHCellIterator &c,
           const DHFaceIterator &f,
           const unsigned int    face_no,
           const unsigned int    subface_no);

    /**
     * 切换到同一单元格的新面。不改变 @p indices
     * ，不重置LocalResults中的数据。
     *
     */
    template <class DHFaceIterator>
    void
    set_face(const DHFaceIterator &f, const unsigned int face_no);

    /**
     * 切换到同一单元的一个新的子面。不改变 @p indices
     * ，不重置LocalResults中的数据。
     *
     */
    template <class DHFaceIterator>
    void
    set_subface(const DHFaceIterator &f,
                const unsigned int    face_no,
                const unsigned int    subface_no);

    const BlockIndices &
    local_indices() const;


    /// The block structure of the system
    SmartPointer<const BlockInfo, DoFInfo<dim, spacedim>> block_info;

    /**
     * 该结构指的是具有水平数据的单元，而不是活动数据。
     *
     */
    bool level_cell;

  private:
    /**
     * 标准构造函数，不设置任何区块索引。不推荐使用这个构造函数，但DoFInfoBox中的数组需要它。
     *
     */
    DoFInfo();

    /// Set up local block indices
    void
    set_block_indices();

    /// Fill index vector with active indices
    template <class DHCellIterator>
    void
    get_indices(const DHCellIterator &c);

    /// Auxiliary vector
    std::vector<types::global_dof_index> indices_org;

    /**
     * 如果#block_info没有被设置，则创建一个辅助的本地BlockIndices对象。它只包含每个单元的自由度大小的单一块。
     *
     */
    BlockIndices aux_local_indices;

    friend class DoFInfoBox<dim, DoFInfo<dim, spacedim, number>>;
  };


  /**
   * 一个捆绑单元上使用的 MeshWorker::DoFInfo 对象的类。
   * @todo
   * 目前，我们为单元存储一个对象，为每个面存储两个对象。我们可以将所有与单元格本身有关的面的数据收集在一个对象中，这样可以节省一点内存和一些操作，但会牺牲一些清洁度。
   * @ingroup MeshWorker
   *
   */
  template <int dim, class DOFINFO>
  class DoFInfoBox
  {
  public:
    /**
     * 构造函数将种子复制到所有其他对象中。
     *
     */
    DoFInfoBox(const DOFINFO &seed);

    /**
     * 复制构造函数，取#cell并将其作为其他构造函数中的种子。
     *
     */
    DoFInfoBox(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * 复制赋值运算符，以另一个对象为种子。
     *
     */
    DoFInfoBox &
    operator=(const DoFInfoBox<dim, DOFINFO> &);

    /**
     * 重置所有的可用标志。
     *
     */
    void
    reset();

    /**
     * 在所有DOFINFO对象被适当填充后，使用ASSEMBLER对象将它们组装成全局数据。参见
     * MeshWorker::Assembler 中的可用类。
     *
     */
    template <class ASSEMBLER>
    void
    assemble(ASSEMBLER &ass) const;


    /**
     * 单元的数据。
     *
     */
    DOFINFO cell;
    /**
     * 里面的面的数据。
     *
     */
    DOFINFO interior[GeometryInfo<dim>::faces_per_cell];
    /**
     * 来自外部的面的数据。
     *
     */
    DOFINFO exterior[GeometryInfo<dim>::faces_per_cell];

    /**
     * 一组标志，表明内部面的数据是否可用。
     *
     */
    bool interior_face_available[GeometryInfo<dim>::faces_per_cell];

    /**
     * 一组标志，表明外部面的数据是否可用。
     *
     */
    bool exterior_face_available[GeometryInfo<dim>::faces_per_cell];

    /**
     * 一个标志，指定当前对象是否已被设置为有效单元。
     *
     */
    bool cell_valid;
  };

  //----------------------------------------------------------------------//

  template <int dim, int spacedim, typename number>
  DoFInfo<dim, spacedim, number>::DoFInfo()
    : face_number(numbers::invalid_unsigned_int)
    , sub_number(numbers::invalid_unsigned_int)
    , level_cell(false)
  {}



  template <int dim, int spacedim, typename number>
  DoFInfo<dim, spacedim, number>::DoFInfo(
    const DoFHandler<dim, spacedim> &dof_handler)
    : face_number(numbers::invalid_unsigned_int)
    , sub_number(numbers::invalid_unsigned_int)
    , level_cell(false)
  {
    std::vector<types::global_dof_index> aux(1);
    aux[0] = dof_handler.get_fe().n_dofs_per_cell();
    aux_local_indices.reinit(aux);
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator>
  inline void
  DoFInfo<dim, spacedim, number>::get_indices(const DHCellIterator &c)
  {
    indices.resize(c->get_fe().n_dofs_per_cell());
    if (block_info == nullptr || block_info->local().size() == 0)
      c->get_active_or_mg_dof_indices(indices);
    else
      {
        indices_org.resize(c->get_fe().n_dofs_per_cell());
        c->get_active_or_mg_dof_indices(indices_org);
        set_block_indices();
      }
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator>
  inline void
  DoFInfo<dim, spacedim, number>::reinit(const DHCellIterator &c)
  {
    get_indices(c);
    level_cell = c->is_level_cell();

    cell        = typename Triangulation<dim, spacedim>::cell_iterator(*c);
    face_number = numbers::invalid_unsigned_int;
    sub_number  = numbers::invalid_unsigned_int;
    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template <int dim, int spacedim, typename number>
  template <class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim, number>::set_face(const DHFaceIterator &f,
                                           const unsigned int    face_no)
  {
    face = static_cast<typename Triangulation<dim, spacedim>::face_iterator>(f);
    face_number = face_no;
    sub_number  = numbers::invalid_unsigned_int;
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim, number>::reinit(const DHCellIterator &c,
                                         const DHFaceIterator &f,
                                         const unsigned int    face_no)
  {
    if ((cell.state() != IteratorState::valid) ||
        cell != typename Triangulation<dim, spacedim>::cell_iterator(*c))
      get_indices(c);
    level_cell = c->is_level_cell();

    cell = typename Triangulation<dim, spacedim>::cell_iterator(*c);
    set_face(f, face_no);

    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template <int dim, int spacedim, typename number>
  template <class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim, number>::set_subface(const DHFaceIterator &f,
                                              const unsigned int    face_no,
                                              const unsigned int    subface_no)
  {
    face = static_cast<typename Triangulation<dim, spacedim>::face_iterator>(f);
    face_number = face_no;
    sub_number  = subface_no;
  }


  template <int dim, int spacedim, typename number>
  template <class DHCellIterator, class DHFaceIterator>
  inline void
  DoFInfo<dim, spacedim, number>::reinit(const DHCellIterator &c,
                                         const DHFaceIterator &f,
                                         const unsigned int    face_no,
                                         const unsigned int    subface_no)
  {
    if (cell.state() != IteratorState::valid ||
        cell !=
          static_cast<typename Triangulation<dim, spacedim>::cell_iterator>(c))
      get_indices(c);
    level_cell = c->is_level_cell();

    cell = static_cast<typename Triangulation<dim, spacedim>::cell_iterator>(c);
    set_subface(f, face_no, subface_no);

    if (block_info)
      LocalResults<number>::reinit(block_info->local());
    else
      LocalResults<number>::reinit(aux_local_indices);
  }


  template <int dim, int spacedim, typename number>
  inline const BlockIndices &
  DoFInfo<dim, spacedim, number>::local_indices() const
  {
    if (block_info)
      return block_info->local();
    return aux_local_indices;
  }

  //----------------------------------------------------------------------//

  template <int dim, class DOFINFO>
  inline DoFInfoBox<dim, DOFINFO>::DoFInfoBox(const DOFINFO &seed)
    : cell(seed)
    , cell_valid(true)
  {
    for (unsigned int i : GeometryInfo<dim>::face_indices())
      {
        exterior[i]                = seed;
        interior[i]                = seed;
        interior_face_available[i] = false;
        exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  inline DoFInfoBox<dim, DOFINFO>::DoFInfoBox(
    const DoFInfoBox<dim, DOFINFO> &other)
    : cell(other.cell)
    , cell_valid(other.cell_valid)
  {
    for (unsigned int i : GeometryInfo<dim>::face_indices())
      {
        exterior[i]                = other.exterior[i];
        interior[i]                = other.interior[i];
        interior_face_available[i] = false;
        exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  inline DoFInfoBox<dim, DOFINFO> &
  DoFInfoBox<dim, DOFINFO>::operator=(const DoFInfoBox<dim, DOFINFO> &other)
  {
    cell       = other.cell;
    cell_valid = other.cell_valid;
    for (unsigned int i : GeometryInfo<dim>::face_indices())
      {
        exterior[i]                = other.exterior[i];
        interior[i]                = other.interior[i];
        interior_face_available[i] = false;
        exterior_face_available[i] = false;
      }
    return *this;
  }


  template <int dim, class DOFINFO>
  inline void
  DoFInfoBox<dim, DOFINFO>::reset()
  {
    cell_valid = false;
    for (unsigned int i : GeometryInfo<dim>::face_indices())
      {
        interior_face_available[i] = false;
        exterior_face_available[i] = false;
      }
  }


  template <int dim, class DOFINFO>
  template <class ASSEMBLER>
  inline void
  DoFInfoBox<dim, DOFINFO>::assemble(ASSEMBLER &assembler) const
  {
    if (!cell_valid)
      return;

    assembler.assemble(cell);
    for (unsigned int i : GeometryInfo<dim>::face_indices())
      {
        // Only do something if data available
        if (interior_face_available[i])
          {
            // If both data
            // available, it is an
            // interior face
            if (exterior_face_available[i])
              assembler.assemble(interior[i], exterior[i]);
            else
              assembler.assemble(interior[i]);
          }
      }
  }
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


