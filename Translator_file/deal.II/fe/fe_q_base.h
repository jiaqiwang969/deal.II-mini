//include/deal.II-translator/fe/fe_q_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_fe_q_base_h
#define dealii_fe_q_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 该类收集了FE_Q、FE_Q_DG0和FE_Q_Bubbles中使用的基本方法。该类没有公共构造函数，因为它没有独立的功能。定义的完成是留给派生类的。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Q_Base : public FE_Poly<dim, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  FE_Q_Base(const ScalarPolynomialsBase<dim> &poly_space,
            const FiniteElementData<dim> &    fe_data,
            const std::vector<bool> &         restriction_is_additive_flags);

  /**
   * 返回从给定的有限元插值到现在的矩阵。然后矩阵的大小是
   * @p dofs_per_cell 乘以<tt>source.n_dofs_per_cell()</tt>。
   * 这些矩阵只有在源元素也是 @p FE_Q
   * 元素时才可用。否则，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;


  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。FE_Q元素家族只为同一类型的元素和FE_Nothing提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。FE_Q元素家族只为同一类型的元素和FE_Nothing提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零的函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 从精细网格空间投射到粗略网格空间。重写FiniteElement中的相应方法，实现懒人评估（在请求时初始化）。
   * 如果这个投影运算符与一个矩阵 @p P,
   * 相关联，那么这里将返回这个矩阵 @p P_i
   * 对单个子单元的限制。    矩阵 @p P 是单元格矩阵 @p
   * P_i的串联或相加，取决于#restriction_is_additive_flags。这区分了插值（连接）和标量积（求和）方面的投影。
   * 行和列指数分别与粗网格和细网格空间相关，与相关运算符的定义一致。
   * 如果投影矩阵没有在派生的有限元类中实现，这个函数将以ExcProjectionVoid中止。你可以通过调用restriction_is_implemented()或isotropic_restriction_is_implemented()函数来检查是否属于这种情况。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 网格间的嵌入矩阵。重述FiniteElement中的相应方法，实现懒人评估（查询时初始化）。
   * 从粗网格空间到细网格空间的身份运算符与一个矩阵 @p
   * P. 相关联，该矩阵 @p P_i
   * 对单个子单元的限制在此返回。    矩阵 @p P
   * 是串联的，而不是单元格矩阵 @p
   * P_i的总和。也就是说，如果同一个非零条目<tt>j,k</tt>存在于两个不同的子矩阵
   * @p P_i,
   * 中，该值在两个矩阵中应该是相同的，它只被复制到矩阵
   * @p P 中一次。
   * 行和列指数分别与细格和粗格空间相关，与相关运算符的定义一致。
   * 这些矩阵被组装多层次方法的延长矩阵的程序所使用。
   * 在使用这个矩阵阵列组装单元间的转移矩阵时，延长矩阵中的零元素被丢弃，不会填满转移矩阵。
   * 如果投影矩阵没有在派生的有限元类中实现，这个函数会以ExcEmbeddingVoid中止。你可以通过调用prolongation_is_implemented()或isotropic_prolongation_is_implemented()函数来检查是否属于这种情况。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 给出一个面的指数自然排序中的指数，返回单元格上相同自由度的指数。
   * 为了解释这个概念，考虑这样的情况：我们想知道一个面的自由度，例如作为FESystem元素的一部分，是否是原始的。不幸的是，FiniteElement类中的is_primitive()函数需要一个单元格索引，所以我们需要找到对应于当前面的索引的形状函数的单元格索引。
   * 这个函数可以做到这一点。
   * 实现这一点的代码将看起来像这样。
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   * if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   * ... do whatever
   * @endcode
   * 这个函数需要额外的参数，这些参数考虑到实际的面可以是相对于所考虑的单元格的标准排序，或者可以是翻转的，定向的，等等。
   * @param  face_dof_index
   * 一个面的自由度的索引。这个指数必须在零和每个面的自由度之间。
   * @param  face
   * 这个自由度所在的面的编号。这个数字必须介于零和
   * GeometryInfo::faces_per_cell.   @param  face_orientation
   * 描述面的方向的一部分。见  @ref GlossFaceOrientation  。
   * @param  face_flip 对脸部方向的一部分描述。参见  @ref
   * GlossFaceOrientation  。    @param  face_rotation
   * 描述脸部方向的一部分。见  @ref GlossFaceOrientation  。
   * @return
   * 这个自由度在整个单元上的自由度集合中的索引。返回值将介于0和dofs_per_cell之间。
   *
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int face_dof_index,
                     const unsigned int face,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation = false) const override;

  /**
   * 返回元素的恒定模式的列表。对于这个元素，该列表由所有组件的真参数组成。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @name  支持hp的函数  
     * @{ 
   *
   */

  /**
   * 返回这个元素是否以新的方式实现了它的悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 对于FE_Q类来说，结果总是真的（与元素的程度无关），因为它实现了hp-capability所需的完整功能集。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FE的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  //@}

  /**
   * 尝试构造一个0度的FE_Q对象
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcFEQCannotHaveDegree0,
                   "FE_Q can only be used for polynomial degrees "
                   "greater than zero. If you want an element of polynomial "
                   "degree zero, then it cannot be continuous and you "
                   "will want to use FE_DGQ<dim>(0).");

protected:
  /**
   * 仅供内部使用。其全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数中需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 执行基于一维支持点的元素初始化，即设置重新编号，初始化单位支持点，初始化约束以及限制和延长矩阵。
   *
   */
  void
  initialize(const std::vector<Point<1>> &support_points_1d);

  /**
   * 初始化悬挂节点的约束矩阵。从initialize()调用。
   *
   */
  void
  initialize_constraints(const std::vector<Point<1>> &points);

  /**
   * 初始化FiniteElement类的 @p unit_support_points 字段。
   * 从initialize()中调用。
   *
   */
  void
  initialize_unit_support_points(const std::vector<Point<1>> &points);

  /**
   * 初始化FiniteElement类的 @p unit_face_support_points
   * 字段。从initialize()中调用。
   *
   */
  void
  initialize_unit_face_support_points(const std::vector<Point<1>> &points);

  /**
   * 初始化FiniteElement类的 @p
   * adjust_quad_dof_index_for_face_orientation_table
   * 字段。从initialize()中调用。
   *
   */
  void
  initialize_quad_dof_index_permutation();

  /**
   * 向前声明一个类，我们把实现的重要部分放入其中。
   * 参见.cc文件以了解更多信息。
   *
   */
  struct Implementation;

  // Declare implementation friend.
  friend struct FE_Q_Base<dim, spacedim>::Implementation;

private:
  /**
   * 用于保护限制和嵌入矩阵的初始化的Mutex。
   *
   */
  mutable Threads::Mutex mutex;

  /**
   * 底层张量乘积空间的最高多项式度数，不需要任何丰富的内容。对于FE_Q*(p)来说，这是p。请注意，富集可能会导致度数的差异。
   *
   */
  const unsigned int q_degree;
};


 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


