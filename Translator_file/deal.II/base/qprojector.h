//include/deal.II-translator/base/qprojector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_qprojector_h
#define dealii_qprojector_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup Quadrature */ 
 /*@{*/ 


/**
 * 该类是一个辅助类，以方便在单元格的面或子面上使用正交公式。它从正交对象中计算出单元格上正交点的位置，其流形的维度小于单元格的维度和面的数量。例如，给出一维的辛普森规则并使用面数为1的project_to_face()函数，返回的点将是（1,0）、（1,0.5）和（1,1）。请注意，面有一个方向，所以当投射到面3时，你会得到(0,0),
 * (0,0.5)和(0,1)，这是顺时针方向的，而面1的点是逆时针方向的。
 * 对子面的投影（即对单元格的一个面的子面），与上述情况相同。请注意一个面的子女被编号的顺序，这在二维空间中与该面的方向相吻合。
 * 第二组函数通过在<b>all</b>面和子面上投影给定的正交规则来生成一个正交公式。这在FEFaceValues和FESubfaceValues类中使用。由于我们现在在一个数组中拥有所有面和子面的正交点，我们需要有一种方法来找到这个数组中一个面或子面所对应的点和权值的起始索引。这可以通过DataSetDescriptor成员类来完成。
 * 不同的函数被分组到一个共同的类中，以避免将它们放入全局命名空间。然而，由于它们没有本地数据，所有的函数都被声明为<tt>static</tt>，无需创建这个类的对象就可以调用。
 * 对于三维的情况，你应该注意面的方向比二维的更加复杂。正交公式是在面的标准方向上投影的，而不是在六面体的内部或外部。让事情更复杂的是，在三维中我们允许面有两个方向（可以用<tt>cell->face_orientation(face)</tt>来识别），所以我们必须将正交公式投射到两个方向的面和子面。(关于不同面的方向的描述，请参考Triangulation类的文档，以及 @ref GlossFaceOrientation "面的方向的词汇表条目"
 * ，以获得更多相关信息)。DataSetDescriptor成员类用于识别每个数据集的起始位置。
 *
 *
 */
template <int dim>
class QProjector
{
public:
  /**
   * 定义一个正交的别名，作用于少一维的对象。对于单元格来说，这将是一个面的四分法。
   *
   */
  using SubQuadrature = Quadrature<dim - 1>;

  /**
   * 如果给定的正交公式用在面<tt>face_no</tt>上，计算单元格上的正交点。进一步的细节，请看这个类的一般文档。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六角形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static void
  project_to_face(const SubQuadrature &    quadrature,
                  const unsigned int       face_no,
                  std::vector<Point<dim>> &q_points);

  /**
   * 如果给定的正交公式用在面<tt>face_no</tt>上，计算单元格上的正交点。进一步的细节，请看这个类的一般文档。
   *
   */
  static void
  project_to_face(const ReferenceCell      reference_cell,
                  const SubQuadrature &    quadrature,
                  const unsigned int       face_no,
                  std::vector<Point<dim>> &q_points);

  /**
   * 计算单元格正交公式，对应于在面<tt>face_no</tt>上使用<tt>quadrature</tt>。进一步的细节，请看这个类的一般文档。
   * @note
   * 这个函数已被弃用，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_face(const SubQuadrature &quadrature, const unsigned int face_no);

  /**
   * 计算对应于在面<tt>face_no</tt>上使用<tt>quadrature</tt>的单元格正交公式。进一步的细节，请看这个类的一般文档。
   *
   */
  static Quadrature<dim>
  project_to_face(const ReferenceCell  reference_cell,
                  const SubQuadrature &quadrature,
                  const unsigned int   face_no);

  /**
   * 如果在面<tt>face_no</tt>上使用给定的正交公式，计算单元格上的正交点，子面编号<tt>subface_no</tt>对应于
   * RefineCase::Type
   * <tt>ref_case</tt>。最后一个参数只在3D中使用。
   * @note
   * 只对点进行变换。正交权重与原始规则的权重相同。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（一维）、一个四边形（二维）或一个六边形（三维）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static void
  project_to_subface(const SubQuadrature &          quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     std::vector<Point<dim>> &      q_points,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * 如果给定的正交公式用在面<tt>face_no</tt>上，子面编号<tt>subface_no</tt>对应于
   * RefineCase::Type
   * <tt>ref_case</tt>，计算单元格上的正交点。最后一个参数只在3D中使用。
   * @note
   * 只对点进行变换。正交权重与原始规则的权重相同。
   *
   */
  static void
  project_to_subface(const ReferenceCell            reference_cell,
                     const SubQuadrature &          quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     std::vector<Point<dim>> &      q_points,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * 计算对应于在面<tt>face_no</tt>的子面<tt>quadrature</tt>上使用RefinementCase<dim-1>
   * <tt>ref_case</tt>的单元正交公式。最后一个参数只在3D中使用。
   * @note
   * 只对点进行转换。正交权重与原始规则的权重相同。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_subface(const SubQuadrature &          quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * 计算对应于在面<tt>face_no</tt>的子面<tt>quadrature</tt>上使用RefinementCase<dim-1>
   * <tt>ref_case</tt>的单元格正交公式。最后一个参数只在3D中使用。
   * @note
   * 只对点进行变换。正交权重与原始规则的权重相同。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  static Quadrature<dim>
  project_to_subface(const ReferenceCell            reference_cell,
                     const SubQuadrature &          quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * 取一个面的正交公式并从中生成一个单元格的正交公式，其中给定参数的正交点被投射到所有面上。
   * 新规则的权重是原始权重的复制。
   * 因此，权重的总和不是一个，而是面的数量，也就是参考单元的表面。
   * 这特别允许我们提取对应于一个面的点的子集，并将其作为这个面的正交点，就像在FEFaceValues中做的那样。
   * @note
   * 在三维中，这个函数为每个面产生八组正交点，以应对可能不同方向的网格。
   * @note
   * 这个函数已被弃用，因为它隐含了一个假设，即单元格是一条线（1D），一个四边形（2D），或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_all_faces(const Quadrature<dim - 1> &quadrature);

  /**
   * 取一个面的正交公式集合，并从中生成一个单元格正交公式，其中给定参数的正交点被投射到所有面上。
   * 新规则的权重是原始权重的复制。
   * 因此，权重的总和不是一个，而是面的数量，也就是参考单元格的表面。
   * 这尤其允许我们提取对应于一个面的点的子集，并将其作为这个面的正交点，正如在FEFaceValues中所做的那样。
   * @note
   * 在三维中，这个函数为每个面产生八组正交点，以应对可能不同方向的网格。
   *
   */
  static Quadrature<dim>
  project_to_all_faces(const ReferenceCell             reference_cell,
                       const hp::QCollection<dim - 1> &quadrature);

  /**
   * 像上面的函数一样，在所有的面都应用相同的面正交公式。
   *
   */
  static Quadrature<dim>
  project_to_all_faces(const ReferenceCell        reference_cell,
                       const Quadrature<dim - 1> &quadrature);

  /**
   * 取一个面的正交公式并从中生成一个单元格的正交公式，其中给定参数的正交点被投射到所有子面上。
   * 和project_to_all_faces()一样，新规则的权重加起来是面的数量（不是子面），也就是参考单元格的面。
   * 这尤其允许我们提取对应于一个子面的点的子集，并将其作为这个面的正交点，就像在FESubfaceValues中所做的那样。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六角形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_all_subfaces(const SubQuadrature &quadrature);

  /**
   * 取一个面的正交公式并从中生成一个单元格的正交公式，其中给定参数的正交点被投射到所有子面上。
   * 和project_to_all_faces()一样，新规则的权重加起来是面的数量（不是子面），也就是参考单元格的面。
   * 这特别允许我们提取对应于一个子面的点的子集，并将其作为这个面的正交点，就像在FESubfaceValues中所做的那样。
   *
   */
  static Quadrature<dim>
  project_to_all_subfaces(const ReferenceCell  reference_cell,
                          const SubQuadrature &quadrature);

  /**
   * 将一个给定的正交公式投射到一个单元格的孩子身上。你可能想使用这个函数，以防你只想在一个潜在的子单元所占的区域内扩展一个积分。子单元的编号与细化单元时的编号相同。
   * 由于使用这个正交公式的积分现在只扩展到了单元格的一部分，所产生的对象的权重被
   * GeometryInfo<dim>::children_per_cell. 除以。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_child(const Quadrature<dim> &quadrature,
                   const unsigned int     child_no);

  /**
   * 将一个给定的正交公式投射到一个单元格的子单元。你可能想使用这个函数，如果你想只在一个潜在的子单元所占的区域内扩展一个积分。子单元的编号与细化单元时的编号相同。
   * 由于使用这个正交公式的积分现在只扩展到了单元格的一部分，所产生的对象的权重被
   * GeometryInfo<dim>::children_per_cell. 除以。
   *
   */
  static Quadrature<dim>
  project_to_child(const ReferenceCell    reference_cell,
                   const Quadrature<dim> &quadrature,
                   const unsigned int     child_no);

  /**
   * 将一个正交规则投射到一个单元格的所有子单元。与project_to_all_subfaces()类似，这个函数将project_to_child()生成的公式复制到所有子单元，这样权重加起来就是1，也就是整个单元的体积。
   * 子单元的编号与细化单元时的编号相同。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六角形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_all_children(const Quadrature<dim> &quadrature);

  /**
   * 投射一个正交规则到一个单元格的所有子单元。与project_to_all_subfaces()类似，这个函数为所有子单元复制project_to_child()生成的公式，这样权重加起来就是1，也就是总单元的体积。
   * 子单元的编号与细化单元时的编号相同。
   *
   */
  static Quadrature<dim>
  project_to_all_children(const ReferenceCell    reference_cell,
                          const Quadrature<dim> &quadrature);

  /**
   * 将一维规则<tt>正交</tt>投影到连接点<tt>p1</tt>和<tt>p2</tt>的直线上。
   * @note
   * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条直线（一维）、一个四边形（二维）或一个六边形（三维）。请使用该函数的另一个版本，该版本接受参考单元格类型。
   *
   */
  DEAL_II_DEPRECATED static Quadrature<dim>
  project_to_line(const Quadrature<1> &quadrature,
                  const Point<dim> &   p1,
                  const Point<dim> &   p2);

  /**
   * 将一维规则<tt>正交</tt>投影到连接点<tt>p1</tt>和<tt>p2</tt>的直线。
   *
   */
  static Quadrature<dim>
  project_to_line(const ReferenceCell  reference_cell,
                  const Quadrature<1> &quadrature,
                  const Point<dim> &   p1,
                  const Point<dim> &   p2);

  /**
   * 由于project_to_all_faces()和project_to_all_subfaces()函数将一个面的正交公式的所有投影的正交点和权重连锁到一个单元格的面或子面上，我们需要一种方法来识别一个特定面或子面的点和权重的起始索引在哪里。这个类提供了这一点：有一些静态成员函数可以生成这种类型的对象，给定面或子面的索引，然后你可以用生成的对象来代替一个整数，表示一个给定数据集的偏移。
   *
   */
  class DataSetDescriptor
  {
  public:
    /**
     * 默认构造函数。除了生成一个无效的索引外，这并没有什么作用，因为你没有给出你想要的单元、面或子面的有效描述符。
     *
     */
    DataSetDescriptor();

    /**
     * 静态函数用于生成一个单元的偏移量。由于我们每个正交对象只有一个单元，这个偏移量当然是零，但是为了与其他静态函数保持一致，我们把这个函数带在身边。
     *
     */
    static DataSetDescriptor
    cell();

    /**
     * 静态函数为一个单元格的给定面生成一个具有给定面方向、翻转和旋转的偏移对象。当然，这个函数只允许在<tt>dim>=2</tt>的情况下使用，如果空间维度等于2，则忽略面的方向、翻转和旋转。
     * 最后一个参数表示低维面的正交公式（已经投射到面的正交公式）的正交点的数量。
     * @note
     * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
     *
     */
    DEAL_II_DEPRECATED static DataSetDescriptor
    face(const unsigned int face_no,
         const bool         face_orientation,
         const bool         face_flip,
         const bool         face_rotation,
         const unsigned int n_quadrature_points);

    /**
     * 静态函数，为一个单元格的给定面生成一个具有给定面方向、翻转和旋转的偏移对象。当然，这个函数只允许在<tt>dim>=2</tt>的情况下使用，如果空间维度等于2，面的方向、翻转和旋转就被忽略了。
     * 最后一个参数表示低维面的正交公式（已经投射到面的正交公式）的正交点的数量。
     *
     */
    static DataSetDescriptor
    face(const ReferenceCell reference_cell,
         const unsigned int  face_no,
         const bool          face_orientation,
         const bool          face_flip,
         const bool          face_rotation,
         const unsigned int  n_quadrature_points);

    /**
     * 像上面的函数一样，但取一个正交集合，使每个面可能有不同数量的正交点。
     *
     */
    static DataSetDescriptor
    face(const ReferenceCell             reference_cell,
         const unsigned int              face_no,
         const bool                      face_orientation,
         const bool                      face_flip,
         const bool                      face_rotation,
         const hp::QCollection<dim - 1> &quadrature);

    /**
     * 静态函数，为一个单元格的给定子面生成一个偏移对象，并给出面的方向、翻转和旋转。当然这个函数只允许在<tt>dim>=2</tt>的情况下使用，如果空间维度等于2，面的方向、翻转和旋转将被忽略。
     * 最后但一个参数表示低维面的正交公式（被投射到面的正交公式）的正交点的数量。
     * 通过最后一个参数，各向异性的细化可以得到尊重。
     * @note
     * 这个函数已被废弃，因为它隐含了一个假设，即单元格是一条线（1D）、一个四边形（2D）或一个六边形（3D）。请使用该函数的另一个版本，该版本接受参考单元格类型。
     *
     */
    DEAL_II_DEPRECATED static DataSetDescriptor
    subface(const unsigned int               face_no,
            const unsigned int               subface_no,
            const bool                       face_orientation,
            const bool                       face_flip,
            const bool                       face_rotation,
            const unsigned int               n_quadrature_points,
            const internal::SubfaceCase<dim> ref_case =
              internal::SubfaceCase<dim>::case_isotropic);

    /**
     * 静态函数，为一个单元格的给定子面生成一个具有给定面方向、翻转和旋转的偏移对象。当然，这个函数只允许在<tt>dim>=2</tt>的情况下使用，如果空间维度等于2，面的方向、翻转和旋转将被忽略。
     * 最后但一个参数表示低维面的正交公式（被投射到面的正交公式）的正交点的数量。
     * 通过最后一个参数，各向异性的细化可以得到尊重。
     *
     */
    static DataSetDescriptor
    subface(const ReferenceCell              reference_cell,
            const unsigned int               face_no,
            const unsigned int               subface_no,
            const bool                       face_orientation,
            const bool                       face_flip,
            const bool                       face_rotation,
            const unsigned int               n_quadrature_points,
            const internal::SubfaceCase<dim> ref_case =
              internal::SubfaceCase<dim>::case_isotropic);

    /**
     * 转化为一个整数，表示该数据集的第一个元素在所有投射到面和子面上的正交公式集中的偏移。这个转换操作符允许我们使用偏移量描述符对象来代替整数偏移量。
     *
     */
    operator unsigned int() const;

  private:
    /**
     * 存储一个给定单元、面或子面的整数偏移。
     *
     */
    const unsigned int dataset_offset;

    /**
     * 这是真正的构造函数，但它是私有的，因此只对上面的静态成员函数可用。
     *
     */
    DataSetDescriptor(const unsigned int dataset_offset);
  };
};

 /*@}*/ 


// -------------------  inline and template functions ----------------



template <int dim>
inline QProjector<dim>::DataSetDescriptor::DataSetDescriptor(
  const unsigned int dataset_offset)
  : dataset_offset(dataset_offset)
{}


template <int dim>
inline QProjector<dim>::DataSetDescriptor::DataSetDescriptor()
  : dataset_offset(numbers::invalid_unsigned_int)
{}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::cell()
{
  return 0;
}



template <int dim>
inline QProjector<dim>::DataSetDescriptor::operator unsigned int() const
{
  return dataset_offset;
}



template <int dim>
Quadrature<dim> inline QProjector<dim>::project_to_all_faces(
  const Quadrature<dim - 1> &quadrature)
{
  return project_to_all_faces(ReferenceCells::get_hypercube<dim>(), quadrature);
}


template <int dim>
Quadrature<dim> inline QProjector<dim>::project_to_all_faces(
  const ReferenceCell        reference_cell,
  const Quadrature<dim - 1> &quadrature)
{
  return project_to_all_faces(reference_cell,
                              hp::QCollection<dim - 1>(quadrature));
}


 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN


template <>
void
QProjector<1>::project_to_face(const Quadrature<0> &,
                               const unsigned int,
                               std::vector<Point<1>> &);
template <>
void
QProjector<1>::project_to_face(const ReferenceCell reference_cell,
                               const Quadrature<0> &,
                               const unsigned int,
                               std::vector<Point<1>> &);
template <>
void
QProjector<2>::project_to_face(const Quadrature<1> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points);
template <>
void
QProjector<2>::project_to_face(const ReferenceCell    reference_cell,
                               const Quadrature<1> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points);
template <>
void
QProjector<3>::project_to_face(const Quadrature<2> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points);
template <>
void
QProjector<3>::project_to_face(const ReferenceCell    reference_cell,
                               const Quadrature<2> &  quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points);

template <>
Quadrature<1>
QProjector<1>::project_to_all_faces(const ReferenceCell       reference_cell,
                                    const hp::QCollection<0> &quadrature);


template <>
void
QProjector<1>::project_to_subface(const Quadrature<0> &,
                                  const unsigned int,
                                  const unsigned int,
                                  std::vector<Point<1>> &,
                                  const RefinementCase<0> &);
template <>
void
QProjector<1>::project_to_subface(const ReferenceCell reference_cell,
                                  const Quadrature<0> &,
                                  const unsigned int,
                                  const unsigned int,
                                  std::vector<Point<1>> &,
                                  const RefinementCase<0> &);
template <>
void
QProjector<2>::project_to_subface(const Quadrature<1> &  quadrature,
                                  const unsigned int     face_no,
                                  const unsigned int     subface_no,
                                  std::vector<Point<2>> &q_points,
                                  const RefinementCase<1> &);
template <>
void
QProjector<2>::project_to_subface(const ReferenceCell    reference_cell,
                                  const Quadrature<1> &  quadrature,
                                  const unsigned int     face_no,
                                  const unsigned int     subface_no,
                                  std::vector<Point<2>> &q_points,
                                  const RefinementCase<1> &);
template <>
void
QProjector<3>::project_to_subface(const Quadrature<2> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>> &  q_points,
                                  const RefinementCase<2> &face_ref_case);
template <>
void
QProjector<3>::project_to_subface(const ReferenceCell      reference_cell,
                                  const Quadrature<2> &    quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>> &  q_points,
                                  const RefinementCase<2> &face_ref_case);

template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const Quadrature<0> &quadrature);
template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const ReferenceCell  reference_cell,
                                       const Quadrature<0> &quadrature);


#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


