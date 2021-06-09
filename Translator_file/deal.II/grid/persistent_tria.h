//include/deal.II-translator/grid/persistent_tria_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_persistent_tria_h
#define dealii_persistent_tria_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/grid/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类处理三角形的历史，可以在它被删除一段时间后重建它。它的主要目的是支持与时间有关的问题，即由于内存压力，人们经常删除一个三角形，然后想重建它；这个类有所有的信息，可以完全按照以前的方式重建它，包括将单元数映射到几何单元上。
 * 基本上，这是对三角剖面的一个直接替换。因为它是从三角形类派生出来的，所以它共享所有的功能，但是它重写了一些虚拟函数，也增加了一些函数。基类的主要变化是它重写了
 * @p
 * execute_coarsening_and_refinement函数，新版本首先存储所有细化和粗化标志，然后才调用基类的相应函数。存储的标志以后可以用来恢复网格，就像以前一样。其他一些函数也有小幅扩展，更多信息请看它们的文档。
 * 我们注意到，由于三角剖面是在与之前完全相同的状态下创建的，其他对象对它的处理也应该是相同的状态。这尤其适用于DoFHandler对象，它将为原始单元和重建三角剖分后的单元分配相同的自由度。因此，你也可以安全地在重建的网格上使用在原始网格上计算的数据向量。
 *
 *  <h3>Usage</h3>
 * 你可以用几乎与三角形类的对象相同的方式使用这个类的对象。为数不多的区别之一是，你只能通过给构造函数一个粗略的网格来构造这样一个对象。粗略的网格将被用来作为三角形的基础，因此，粗略的网格的寿命必须长于这个类的对象的寿命。
 * 基本上，用法是这样的。
 *
 * @code
 * Triangulation<dim> coarse_grid;
 * ...                     // initialize coarse grid
 *
 * PersistentTriangulation<dim> grid (coarse_grid);
 *
 * for (...)
 *   {
 *                         // restore grid from coarse grid
 *                         // and stored refinement flags
 *     grid.restore ();
 *     ...                 // do something with the grid
 *
 *     ...                 // flag some cells for refinement
 *                         // or coarsening
 *     grid.execute_coarsening_and_refinement ();
 *                         // actually refine grid and store
 *                         // the flags
 *
 *     ...                 // so something more with the grid
 *
 *     grid.clear ();      // delete the grid, but keep the
 *                         // refinement flags for later use
 *                         // in grid.restore() above
 *
 *     ...                 // do something where the grid
 *                         // is not needed anymore, e.g.
 *                         // working with another grid
 *   };
 * @endcode
 *
 * 请注意，最初，PersistentTriangulation对象并不构成一个三角形；只有在第一次调用
 * @p restore 之后，它才成为一个三角形。还要注意的是， @p
 * execute_coarsening_and_refinement
 * 存储了所有必要的标志，以便以后使用 @p restore
 * 函数进行重建。  Triangulation::clear()
 * 将底层三角剖面重置为原始状态，但不影响存储的细化标志，这些标志是以后重建所需的，也不触及restore()中使用的粗略网格。
 *
 *
 * @ingroup grid
 *
 *
 */
template <int dim, int spacedim = dim>
class PersistentTriangulation : public Triangulation<dim, spacedim>
{
public:
  /**
   * 使该尺寸在函数模板中可用。
   *
   */
  static const unsigned int dimension      = dim;
  static const unsigned int spacedimension = spacedim;

  /**
   * 以后从粗略的网格中建立三角测量。也可以从该网格复制平滑标志等。请注意，三角形的初始状态是空的，直到
   * @p restore_grid 被首次调用。
   * 粗略的网格必须持续到这个对象的结束，因为它将在重建网格时被使用。
   *
   */
  PersistentTriangulation(const Triangulation<dim, spacedim> &coarse_grid);

  /**
   * 复制构造函数。只有当要复制的对象所依据的三角形目前是空的，才允许进行这个操作。然而，细化标志以及粗略网格的指针会被复制。
   *
   */
  PersistentTriangulation(
    const PersistentTriangulation<dim, spacedim> &old_tria);

  /**
   * 解构器。
   *
   */
  virtual ~PersistentTriangulation() override = default;

  /**
   * 基类中相同函数的重载版本，存储细化和粗化标志，以便以后重建三角形，之后调用基类的相应函数。
   *
   */
  virtual void
  execute_coarsening_and_refinement() override;

  /**
   * 根据保存的数据恢复网格。为此，粗略的网格被复制，并使用保存的标志逐步重建网格。
   * 请注意，如果底层三角图不是空的，这个函数会导致错误，也就是说，只有在这个对象是新创建的或者之前对它调用过基类的<tt>clear()</tt>函数时才会成功。
   * 重复调用<tt>restore(unsigned
   * int)</tt>函数，在所有细化步骤中循环进行。
   *
   */
  void
  restore();

  /**
   * 差分恢复。执行 @p step_noth
   * 局部细化和粗化步骤。步骤0代表复制粗略网格。
   * 这个函数只有在三角剖面正好处于之前从<tt>step=0...step_no-1</tt>调用restore的状态下才会成功。
   *
   */
  void
  restore(const unsigned int step_no);

  /**
   * 返回细化和粗化步骤的数量。这是由 @p refine_flags
   * 矢量的大小给出的。
   *
   */
  unsigned int
  n_refinement_steps() const;

  /**
   * 重载此函数以使用 @p tria
   * 作为新的粗略网格。目前的三角形和所有存储其历史的细化和粗化标志被删除，并且底层三角形的状态被重置为空，直到下一次调用
   * @p restore_grid 。
   * 粗略的网格必须持续到这个对象的结束，因为它将在重建网格时使用。
   *
   */
  virtual void
  copy_triangulation(const Triangulation<dim, spacedim> &tria) override;

  /**
   * 抛出一个错误，因为这个函数在这个类的范围内没有用。
   *
   */
  virtual void
  create_triangulation(const std::vector<Point<spacedim>> &vertices,
                       const std::vector<CellData<dim>> &  cells,
                       const SubCellData &subcelldata) override;

  /**
   * @copydoc   Triangulation::create_triangulation()
   * Triangulation::create_triangulation() .
   * @note  还没有实现。
   *
   */
  virtual void
  create_triangulation(
    const TriangulationDescription::Description<dim, spacedim>
      &construction_data) override;

  /**
   * 基类中相应函数的重载。
   * 抛出一个错误，因为这个函数在这个类的上下文中没有用。
   *
   */
  DEAL_II_DEPRECATED
  virtual void
  create_triangulation_compatibility(
    const std::vector<Point<spacedim>> &vertices,
    const std::vector<CellData<dim>> &  cells,
    const SubCellData &                 subcelldata) override;

  /**
   * 将所有的细化和粗化标志写到ostream中  @p out.  。
   *
   */
  virtual void
  write_flags(std::ostream &out) const;

  /**
   * 读取之前由<tt>write_flags(...)</tt>写入的所有细化和粗化标志。这对于在程序结束或崩溃后重建三角剖面及其重新启动特别有用。
   *
   */
  virtual void
  read_flags(std::istream &in);

  /**
   * 清除所有标志。保留相同的粗略网格。
   *
   */
  virtual void
  clear_flags();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * 异常情况。
   *
   */
  DeclException0(ExcTriaNotEmpty);
  /**
   * 异常情况。
   *
   */
  DeclException0(ExcFlagsNotCleared);

private:
  /**
   * 此网格应作为粗略的网格使用。
   *
   */
  SmartPointer<const Triangulation<dim, spacedim>,
               PersistentTriangulation<dim, spacedim>>
    coarse_grid;

  /**
   * 持有这个时间层上不同扫描的细化和粗化标志的向量。因此，这些向量持有该网格的历史。
   *
   */
  std::vector<std::vector<bool>> refine_flags;

  /**
   * @ref refine_flags
   *
   */
  std::vector<std::vector<bool>> coarsen_flags;
};


DEAL_II_NAMESPACE_CLOSE

#endif


