//include/deal.II-translator/fe/fe_base_0.txt
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

#ifndef dealii_fe_base_h
#define dealii_fe_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/lac/block_indices.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个专门用于定义Domination枚举以及相关运算符的命名空间。
 *
 *
 */
namespace FiniteElementDomination
{
  /**
   *
   */
  enum Domination
  {
    /**
     * 当前元素占主导地位。
     *
     */
    this_element_dominates,
    /**
     * 其他元素占主导地位。
     *
     */
    other_element_dominates,
    /**
     * 两个元素都不占优势。
     *
     */
    neither_element_dominates,
    /**
     * 任何一个元素都可能占主导地位。
     *
     */
    either_element_can_dominate,
    /**
     * 没有任何要求。
     *
     */
    no_requirements
  };


  /**
   * 二进制 <code>and</code>
   * 运算符的一般化，用于比较关系。其工作方式与你想为向量定义比较关系时差不多：要么第一个向量的所有元素都比第二个向量的元素小、等于或大，要么有些是，有些不是。
   * 这个运算符基本相同：如果两个参数都是
   * <code>this_element_dominates</code> 或
   * <code>other_element_dominates</code>
   * ，那么返回值就是这个值。另一方面，如果其中一个值是
   * <code>either_element_can_dominate</code>
   * ，那么返回值就是另一个参数的值。如果其中一个参数是
   * <code>neither_element_dominates</code> ，或者两个参数是
   * <code>this_element_dominates</code> 和
   * <code>other_element_dominates</code> ，那么返回值是
   * <code>neither_element_dominates</code>  。
   *
   */
  inline Domination operator&(const Domination d1, const Domination d2);
} // namespace FiniteElementDomination

namespace internal
{
  /**
   * 用于设置FiniteElementData的内部数据结构。它为每个对象存储（包括/不包括）自由度的数量，以及它在一个单元中的第一个自由度的索引和每个面中的第一个d维对象的索引。
   * 这些信息被保存为一个向量的向量。人们可以通过：dofs_per_object_inclusive[d][i]来查询第i个d维对象的自由度的包容数。
   * 作为一个例子，数据显示的是一个四边形楔形。它由6个顶点、9条线和5个面（两个三角形和三个四边形）组成。
   * @code
   *            vertices                  lines                  faces cell
   * dpo_excl  1  1  1  1  1  1 | 1  1  1  1  1  1  1  1  1 |  0  0  1  1  1 |  0
   * dpo_incl  1  1  1  1  1  1 | 3  3  3  3  3  3  3  3  3 |  6  6  9  9  9 | 18
   * obj_index 0  1  2  3  4  5 | 6  7  8  9 10 11 12 13 14 | 15 15 15 16 17 | 18
   * @endcode
   * 由于上述表格看起来如下： 对于。
   *
   *
   *
   *
   *
   * - 一个三角形。
   * @code
   * dpo_excl  1  1  1 | 1  1  1 |  0
   * obj_index 0  1  2 | 3  4  5 |  6
   * @endcode
   *
   *
   *
   *
   * - 四边形。
   * @code
   * dpo_excl  1  1  1  1 | 1  1  1  1 |  1
   * obj_index 0  1  2  3 | 4  5  6  7 |  8
   * @endcode
   * 每个面内的第一个d维物体的索引结果为。
   * @code
   *                       vertices      lines       face
   * first_obj_index_on_face 0 0 0 0 0 | 3 3 4 4 4 | 6 6 8 8 8
   * @endcode
   *
   *
   */
  struct GenericDoFsPerObject
  {
    /**
     * 每个物体的专属自由度数。
     *
     */
    std::vector<std::vector<unsigned int>> dofs_per_object_exclusive;

    /**
     * 每个对象的自由度的包容数。
     *
     */
    std::vector<std::vector<unsigned int>> dofs_per_object_inclusive;

    /**
     * 一个对象的第一个索引。
     *
     */
    std::vector<std::vector<unsigned int>> object_index;

    /**
     * 一个物体在一个面中的第一个索引。
     *
     */
    std::vector<std::vector<unsigned int>> first_object_index_on_face;
  };
} // namespace internal

/**
 * 一个声明若干标量常量变量的类，描述了有限元实现的基本属性。这包括，例如，每个顶点、线或单元的自由度的数量；矢量分量的数量；等等。
 * 这里存储的信息是在有限元对象的初始化过程中计算出来的，并通过其构造函数传递给这个类。这个类所存储的数据是有限元类的公共接口的一部分（它派生自当前的类）。更多信息见那里。
 *
 *
 * @ingroup febase
 *
 *
 */
template <int dim>
class FiniteElementData
{
public:
  /**
   * 一个有限元可能具有的不同类型的连续性的枚举器。连续性是由包含构建的有限元空间的Sobolev空间来衡量的，也这样称呼。    请注意，某些连续性可能意味着其他连续性。例如，<i>H<sup>1</sup></i>中的函数在<i>H<sup>curl</sup></i>和<i>H<sup>div</sup></i>中也是如此。    如果你对经典意义上的连续性感兴趣，那么以下关系是成立的。      <ol>   <li>  <i>H<sup>1</sup></i>意味着该函数在单元格边界上是连续的。      <li>  <i>H<sup>2</sup></i>意味着该函数在单元格边界上是连续可微的。      <li>  <i>L<sup>2</sup></i> 表示该元素是不连续的。  由于不连续元素在网格单元之间没有拓扑耦合，代码实际上可能取决于这一属性，<i>L<sup>2</sup></i>符合性以特殊方式处理，即它是<b>not</b>由任何更高符合性所暗示的。    </ol>  为了测试一个有限元是否符合某个空间，使用 FiniteElementData<dim>::conforms().  。
   *
   */
  enum Conformity
  {
    /**
     * 表示一个系统的不兼容的连续性。
     *
     */
    unknown = 0x00,

    /**
     * 不连续的元素。见上文!
     *
     */
    L2 = 0x01,

    /**
     * 与空间的一致性
     * <i>H<sup>curl</sup></i>（矢量场的连续切向分量）。
     *
     */
    Hcurl = 0x02,

    /**
     * 与空间<i>H<sup>div</sup></i>（矢量场的连续法向分量）的符合性
     *
     */
    Hdiv = 0x04,

    /**
     * 与空间<i>H<sup>1</sup></i>的符合性（连续）。
     *
     */
    H1 = Hcurl | Hdiv,

    /**
     * 与空间<i>H<sup>2</sup></i>的符合性（连续可微）。
     *
     */
    H2 = 0x0e
  };

  /**
   * 有限元的尺寸，也就是模板参数<tt>dim</tt>。
   *
   */
  static const unsigned int dimension = dim;

private:
  /**
   * 参考单元的类型。
   *
   */
  const ReferenceCell reference_cell_kind;

  /**
   * 唯一四边形的数量。如果所有的四边形都有相同的类型，其值为1；否则等于四边形的数量。
   *
   */
  const unsigned int number_unique_quads;

  /**
   * 独特的面的数量。如果所有的面都有相同的类型，值是1；否则等于面的数量。
   *
   */
  const unsigned int number_unique_faces;

public:
  /**
   * 一个顶点上的自由度数量。
   *
   */
  const unsigned int dofs_per_vertex;

  /**
   * 一条线的自由度数；不包括线的顶点上的自由度。
   *
   */
  const unsigned int dofs_per_line;

private:
  /**
   * 四边形上的自由度数。如果所有四边形都有相同的自由度数，则数值等于dofs_per_quad。
   *
   */
  const std::vector<unsigned int> n_dofs_on_quad;

public:
  /**
   * 四边形的自由度数；不包括四边形的线和顶点的自由度。
   *
   */
  const unsigned int dofs_per_quad;

private:
  /**
   * 任何四边形上的最大自由度数。
   *
   */
  const unsigned int dofs_per_quad_max;

public:
  /**
   * 六面体的自由度数；不包括六面体的四边形、线和顶点上的自由度。
   *
   */
  const unsigned int dofs_per_hex;

  /**
   * 线上自由度的第一个索引。
   *
   */
  const unsigned int first_line_index;

private:
  /**
   * 一个四边形的第一个索引。如果所有的四边形具有相同的自由度，则只存储第一个四边形的第一个索引，因为其他四边形的索引可以简单地重新计算。
   *
   */
  const std::vector<unsigned int> first_index_of_quads;

public:
  /**
   * 一个四边形上的第一个自由度的索引。
   *
   */
  const unsigned int first_quad_index;

  /**
   * 六面体上的第一个索引。
   *
   */
  const unsigned int first_hex_index;

private:
  /**
   * 所有面的第一行的索引。
   *
   */
  const std::vector<unsigned int> first_line_index_of_faces;

public:
  /**
   * 面的数据在一行中的第一个索引。
   *
   */
  const unsigned int first_face_line_index;

private:
  /**
   * 所有面孔的第一个四边形的索引。
   *
   */
  const std::vector<unsigned int> first_quad_index_of_faces;

public:
  /**
   * 脸部数据在一个四边形上的第一个索引。
   *
   */
  const unsigned int first_face_quad_index;

private:
  /**
   * 面孔上的自由度数。如果所有的面都有相同的自由度数，那么这些值等于dofs_per_quad。
   *
   */
  const std::vector<unsigned int> n_dofs_on_face;

public:
  /**
   * 一个面的自由度数。这是构成一个面的<tt>dim-1</tt>以内的所有物体上自由度的累积数。
   *
   */
  const unsigned int dofs_per_face;

private:
  /**
   * 任何面的最大自由度数。
   *
   */
  const unsigned int dofs_per_face_max;

public:
  /**
   * 一个单元上的总自由度数。这是构成一个单元的所有尺寸到<tt>dim</tt>的对象上的自由度的累积数。
   *
   */
  const unsigned int dofs_per_cell;

  /**
   * 该有限元的矢量分量的数量，以及图像空间的维度。对于矢量值的有限元（即当这个数字大于1时），矢量分量的数量在很多情况下等于在FESystem类的帮助下粘在一起的基本元素的数量。然而，对于像Nedelec元素这样的元素，尽管我们只有一个基础元素，但这个数字还是大于1的。
   *
   */
  const unsigned int components;

  /**
   * 形状函数在单一坐标方向上的最大多项式程度。
   *
   */
  const unsigned int degree;

  /**
   * 表示这个元素所符合的空间。
   *
   */
  const Conformity conforming_space;

  /**
   * 存储一个描述复合元素每个块的尺寸的对象。对于一个不是FESystem的元素，这只包含一个长度为#dofs_per_cell的单一块。
   *
   */
  const BlockIndices block_indices_data;

  /**
   * 构造函数，计算从dofs分布到几何对象的所有必要值。
   * @param[in]  dofs_per_object
   * 一个向量，描述每个维度的几何对象的自由度数量。这个向量的大小必须是dim+1，条目0描述每个顶点的自由度数，条目1描述每条线的自由度数，等等。作为一个例子，对于2d中常见的
   * $Q_1$ 拉格朗日元素，这个向量的元素是 <code>(1,0,0)</code>
   * 。另一方面，对于3D中的 $Q_3$ 元素，它将有条目
   * <code>(1,2,4,8)</code>  。      @param[in]  n_components
   * 元素的向量分量的数量。      @param[in]  degree
   * 这个元素的任何形状函数在参考元素上的任何变量的最大多项式程度。例如，对于
   * $Q_1$
   * 元素(在任何空间维度)，这将是一个；尽管该元素具有
   * $\hat x\hat y$ (在2D)和 $\hat x\hat y\hat z$
   * (在3D)形式的形状函数，尽管是二次和三次多项式，但仍然只在每个参考变量中分别是线性的，这一点就是如此。这个变量所提供的信息通常用于确定什么是合适的正交公式。
   * @param[in]  符合性
   * 描述这个元素符合哪个Sobolev空间的一个变量。例如，
   * $Q_p$ 拉格朗日元素（由FE_Q类实现）是 $H^1$
   * 符合的，而拉维亚-托马斯元素（由FE_RaviartThomas类实现）是
   * $H_\text{div}$
   * 符合的；最后，完全不连续的元素（由FE_DGQ类实现）只有
   * $L_2$ 符合。      @param[in]  block_indices
   * 一个描述有限元的基本元素如何分组的参数。默认值是构建一个由所有
   * @p dofs_per_cell 自由度组成的单一块。这适用于所有 "原子
   * "元素（包括非原始元素），因此这些元素可以省略这个参数。另一方面，像FESystem这样的组成元素会希望在这里传递一个不同的值。
   *
   */
  FiniteElementData(const std::vector<unsigned int> &dofs_per_object,
                    const unsigned int               n_components,
                    const unsigned int               degree,
                    const Conformity                 conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * 与上述相同，但不同的是，也可以指定基础几何实体的类型。
   *
   */
  FiniteElementData(const std::vector<unsigned int> &dofs_per_object,
                    const ReferenceCell              reference_cell,
                    const unsigned int               n_components,
                    const unsigned int               degree,
                    const Conformity                 conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * 与上述相同，但不是传递一个包含每个对象的自由度的向量，而是一个GenericDoFsPerObject类型的结构。这允许二维对象有不同的自由度，这对以三角形和四边形为面的单元格特别有用。
   *
   */
  FiniteElementData(const internal::GenericDoFsPerObject &data,
                    const ReferenceCell                   reference_cell,
                    const unsigned int                    n_components,
                    const unsigned int                    degree,
                    const Conformity                      conformity = unknown,
                    const BlockIndices &block_indices = BlockIndices());

  /**
   * 返回这个元素所定义的参考单元的种类。例如，该元素的参考单元是正方形还是三角形，或更高维度的类似选择。
   *
   */
  ReferenceCell
  reference_cell() const;

  /**
   * 唯一的四边形的数量。如果所有的四边形都有相同的类型，该值为1；否则等于四边形的数量。
   *
   */
  unsigned int
  n_unique_quads() const;

  /**
   * 独特的面的数量。如果所有的面都有相同的类型，值是1；否则等于面的数量。
   *
   */
  unsigned int
  n_unique_faces() const;

  /**
   * 每个顶点的道夫数。
   *
   */
  unsigned int
  n_dofs_per_vertex() const;

  /**
   * 每行的道夫数。不包括低维物体上的道夫。
   *
   */
  unsigned int
  n_dofs_per_line() const;

  /**
   * 每个四边形的道夫数。不包括低维物体上的道夫。
   *
   */
  unsigned int
  n_dofs_per_quad(unsigned int face_no = 0) const;

  /**
   * 每个四边形的最多道夫数。不包括低维物体上的道夫。
   *
   */
  unsigned int
  max_dofs_per_quad() const;

  /**
   * 每个六面体的道夫数。不包括低维物体上的道夫。
   *
   */
  unsigned int
  n_dofs_per_hex() const;

  /**
   * 每个面的度数，累积所有低维物体的自由度。
   *
   */
  unsigned int
  n_dofs_per_face(unsigned int face_no = 0, unsigned int child = 0) const;

  /**
   * 每个面的最大度数，累积所有低维物体的自由度。
   *
   */
  unsigned int
  max_dofs_per_face() const;

  /**
   * 每个单元的自由度数，累积所有低维物体的自由度。
   *
   */
  unsigned int
  n_dofs_per_cell() const;

  /**
   * 返回每个structdim维度对象的度数。对于
   * structdim==0，该函数因此返回 dofs_per_vertex，对于
   * structdim==1
   * dofs_per_line，等等。这个函数主要用于允许一些模板技巧，这些函数应该在各种对象上工作，而不想使用与这些对象相关的不同名称（顶点、线...）。
   *
   */
  template <int structdim>
  unsigned int
  n_dofs_per_object(const unsigned int i = 0) const;

  /**
   * 组件的数量。参见 @ref GlossComponent "术语表 "
   * 以获得更多信息。
   *
   */
  unsigned int
  n_components() const;

  /**
   * 块的数量。参见 @ref GlossBlock "术语表 "
   * 以了解更多信息。
   *
   */
  unsigned int
  n_blocks() const;

  /**
   * 关于区块大小的详细信息。
   *
   */
  const BlockIndices &
  block_indices() const;

  /**
   * 形状函数在单一坐标方向上的最大多项式程度。
   * 这个函数可用于确定最佳正交规则。
   *
   */
  unsigned int
  tensor_degree() const;

  /**
   * 测试一个有限元空间是否符合某个Sobolev空间。
   * @note
   * 即使有限元空间具有比要求的更高的规则性，这个函数也会返回一个真值。
   *
   */
  bool
  conforms(const Conformity) const;

  /**
   * 比较运算符。
   *
   */
  bool
  operator==(const FiniteElementData &) const;

  /**
   * 返回一行中dof的第一个索引。
   *
   */
  unsigned int
  get_first_line_index() const;

  /**
   * 返回四边形上的第一个索引。
   *
   */
  unsigned int
  get_first_quad_index(const unsigned int quad_no = 0) const;

  /**
   * 返回六面体上的第一个索引。
   *
   */
  unsigned int
  get_first_hex_index() const;

  /**
   * 返回面的数据在一条线上的第一个索引。
   *
   */
  unsigned int
  get_first_face_line_index(const unsigned int face_no = 0) const;

  /**
   * 返回一个四面体上的脸部数据的第一个索引。
   *
   */
  unsigned int
  get_first_face_quad_index(const unsigned int face_no = 0) const;
};

namespace internal
{
  /**
   * 转换 @p dim 维度参考单元 @p reference_cell. 的
   * "每个对象的道夫 "信息的实用函数。
   *
   */
  internal::GenericDoFsPerObject
  expand(const unsigned int               dim,
         const std::vector<unsigned int> &dofs_per_object,
         const dealii::ReferenceCell      reference_cell);
} // namespace internal



// --------- inline and template functions ---------------


#ifndef DOXYGEN

namespace FiniteElementDomination
{
  inline Domination operator&(const Domination d1, const Domination d2)
  {
    // go through the entire list of possibilities. note that if we were into
    // speed, obfuscation and cared enough, we could implement this operator
    // by doing a bitwise & (and) if we gave these values to the enum values:
    // neither_element_dominates=0, this_element_dominates=1,
    // other_element_dominates=2, either_element_can_dominate=3
    // =this_element_dominates|other_element_dominates
    switch (d1)
      {
        case this_element_dominates:
          if ((d2 == this_element_dominates) ||
              (d2 == either_element_can_dominate) || (d2 == no_requirements))
            return this_element_dominates;
          else
            return neither_element_dominates;

        case other_element_dominates:
          if ((d2 == other_element_dominates) ||
              (d2 == either_element_can_dominate) || (d2 == no_requirements))
            return other_element_dominates;
          else
            return neither_element_dominates;

        case neither_element_dominates:
          return neither_element_dominates;

        case either_element_can_dominate:
          if (d2 == no_requirements)
            return either_element_can_dominate;
          else
            return d2;

        case no_requirements:
          return d2;

        default:
          // shouldn't get here
          Assert(false, ExcInternalError());
      }

    return neither_element_dominates;
  }
} // namespace FiniteElementDomination


template <int dim>
inline ReferenceCell
FiniteElementData<dim>::reference_cell() const
{
  return reference_cell_kind;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_unique_quads() const
{
  return number_unique_quads;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_unique_faces() const
{
  return number_unique_faces;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_vertex() const
{
  return dofs_per_vertex;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_line() const
{
  return dofs_per_line;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_quad(unsigned int face_no) const
{
  return n_dofs_on_quad[n_dofs_on_quad.size() == 1 ? 0 : face_no];
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::max_dofs_per_quad() const
{
  return dofs_per_quad_max;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_hex() const
{
  return dofs_per_hex;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_face(unsigned int face_no,
                                        unsigned int child_no) const
{
  (void)child_no;

  return n_dofs_on_face[n_dofs_on_face.size() == 1 ? 0 : face_no];
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::max_dofs_per_face() const
{
  return dofs_per_face_max;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_cell() const
{
  return dofs_per_cell;
}



template <int dim>
template <int structdim>
inline unsigned int
FiniteElementData<dim>::n_dofs_per_object(const unsigned int i) const
{
  switch (structdim)
    {
      case 0:
        return n_dofs_per_vertex();
      case 1:
        return n_dofs_per_line();
      case 2:
        return n_dofs_per_quad((structdim == 2 && dim == 3) ? i : 0);
      case 3:
        return n_dofs_per_hex();
      default:
        Assert(false, ExcInternalError());
    }
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_components() const
{
  return components;
}



template <int dim>
inline const BlockIndices &
FiniteElementData<dim>::block_indices() const
{
  return block_indices_data;
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::n_blocks() const
{
  return block_indices_data.size();
}



template <int dim>
inline unsigned int
FiniteElementData<dim>::tensor_degree() const
{
  return degree;
}


template <int dim>
inline bool
FiniteElementData<dim>::conforms(const Conformity space) const
{
  return ((space & conforming_space) == space);
}



template <int dim>
unsigned int
FiniteElementData<dim>::get_first_line_index() const
{
  return first_line_index;
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_quad_index(const unsigned int quad_no) const
{
  if (first_index_of_quads.size() == 1)
    return first_index_of_quads[0] + quad_no * n_dofs_per_quad(0);
  else
    return first_index_of_quads[quad_no];
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_hex_index() const
{
  return first_hex_index;
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_face_line_index(
  const unsigned int face_no) const
{
  return first_line_index_of_faces[first_line_index_of_faces.size() == 1 ?
                                     0 :
                                     face_no];
}

template <int dim>
unsigned int
FiniteElementData<dim>::get_first_face_quad_index(
  const unsigned int face_no) const
{
  return first_quad_index_of_faces[first_quad_index_of_faces.size() == 1 ?
                                     0 :
                                     face_no];
}


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif


