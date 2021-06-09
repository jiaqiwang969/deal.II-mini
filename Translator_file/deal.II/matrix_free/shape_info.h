//include/deal.II-translator/matrix_free/shape_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_shape_info_h
#define dealii_matrix_free_shape_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe.h>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 一个枚举，编码初始化时检测到的元素类型。FEEvaluation将根据给定的元素类型选择最有效的算法。
     * 在类型 ElementType::tensor_symmetric
     * 中有一个隐含的排序，即
     * ElementType::tensor_symmetric_collocation 和
     * ElementType::tensor_symmetric_hermite 也是
     * ElementType::tensor_symmetric. 的类型，同样，
     * ElementType::tensor_symmetric 类型的配置也是
     * ElementType::tensor_general. 的类型。
     * 因此，我们支持在具有这种排序的类型之间进行`<=`操作，但不能针对更高的索引类型，如
     * ElementType::truncated_tensor.
     * @ingroup matrixfree
     *
     */
    enum ElementType
    {
      /**
       * 张量积形状函数，在正交点的形状值评估对应于身份操作，不需要进行插值（拼合方法，也叫谱系评估）。例如，对于一个在Gauss-Lobatto支持点上有节点，在Gauss-Lobatto正交点上有相同阶数的积分的元素来说就是这种情况。
       *
       */
      tensor_symmetric_collocation = 0,

      /**
       * 对称张量积形状函数，满足Hermite特性，其值和一阶导数在一维的元素端点为零。
       *
       */
      tensor_symmetric_hermite = 1,

      /**
       * 通常的张量积形状函数，其形状值和正交点围绕单位区间0.5的中点对称。
       *
       */
      tensor_symmetric = 2,

      /**
       * 没有进一步特殊属性的张量积形状函数
       *
       */
      tensor_general = 3,

      /**
       * 可由截断的张量积描述的完全度数而非张量度数的多项式
       *
       */
      truncated_tensor = 4,

      /**
       * 张量积形状函数，围绕单位区间0.5的中点对称，另外根据FE_Q_DG0增加一个常数形状函数。
       *
       */
      tensor_symmetric_plus_dg0 = 5,

      /**
       * 没有张量积属性的形状函数。
       *
       */
      tensor_none = 6
    };



    /**
     * 该结构存储了张量积有限元和张量积正交公式在参考坐标下的一维截面所评估的形状函数、其梯度和Hessians。这个数据结构还包括在单元边界和子区间
     * $(0, 0.5)$ 和 $(0.5, 1)$ 上对面积分的量的评估。
     *
     */
    template <typename Number>
    struct UnivariateShapeData
    {
      /**
       * 空构造函数。设置默认配置。
       *
       */
      UnivariateShapeData();

      /**
       * 返回这个类的内存消耗，单位是字节。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 编码构造时检测到的元素的类型。FEEvaluation将根据给定的元素类型选择最有效的算法。
       *
       */
      ElementType element_type;

      /**
       * 存储所有一维正交点上评估的一维有限元的形状值。这个数组的长度为<tt>n_dofs_1d
       * n_q_points_1d</tt>，正交点是运行最快的索引。
       *
       */
      AlignedVector<Number> shape_values;

      /**
       * 存储所有一维正交点上评估的一维有限元的形状梯度。这个数组的长度为<tt>n_dofs_1d
       * n_q_points_1d</tt>，正交点是运行最快的索引。
       *
       */
      AlignedVector<Number> shape_gradients;

      /**
       * 存储所有一维正交点上评估的一维有限元的形状Hessians。这个数组的长度为<tt>n_dofs_1d
       * n_q_points_1d</tt>，正交点是运行最快的索引。
       *
       */
      AlignedVector<Number> shape_hessians;

      /**
       * 存储与正交（collocation）相关的形状函数空间的形状梯度，由FE_DGQ<1>（Quadrature<1>）给出。
       *
       */
      AlignedVector<Number> shape_gradients_collocation;

      /**
       * 存储与正交（同位）相关的形状函数空间的形状犹豫度，由FE_DGQ<1>(Quadrature<1>)给出。
       *
       */
      AlignedVector<Number> shape_hessians_collocation;

      /**
       * 以不同的格式存储形状值，即所谓的偶数方案，其中shape_values中的对称性被用于快速评估。
       *
       */
      AlignedVector<Number> shape_values_eo;

      /**
       * 以不同的格式存储形状梯度，即所谓的偶数方案，其中
       * shape_gradients 中的对称性被用于快速评估。
       *
       */
      AlignedVector<Number> shape_gradients_eo;

      /**
       * 以不同的格式存储形状二阶导数，即所谓的偶数方案，其中shape_hessians中的对称性被用于更快的评估。
       *
       */
      AlignedVector<Number> shape_hessians_eo;

      /**
       * 存储与正交（collocation）相关的形状函数空间的形状梯度，由FE_DGQ<1>（Quadrature<1>）给出。这个数组提供了偶数格式的
       * shape_gradients_collocation 字段的替代表示。
       *
       */
      AlignedVector<Number> shape_gradients_collocation_eo;

      /**
       * 存储由FE_DGQ<1>(Quadrature<1>)给出的与正交(collocation)相关的形状函数空间的形状斜率。这个数组提供了偶数格式的shape_hessians_collocation字段的替代表示。
       *
       */
      AlignedVector<Number> shape_hessians_collocation_eo;

      /**
       * 存储从正交点的数据到shape_values字段所定义的基础的反变换。正交点的数据通过其多项式插值隐含地解释，或者显式地用单独的多项式解释，如用`_collocation`字段。数组的大小等于`形状值`数组的布局，它与形状值数组相结合，使这个矩阵是形状值的伪逆。在一维正交点的数量等于基础的大小的情况下，这个数组正好是shape_values数组的逆值。这个数组的长度为<tt>n_dofs_1d
       * n_q_points_1d</tt>，正交点的索引运行速度最快。
       *
       */
      AlignedVector<Number> inverse_shape_values;

      /**
       * 存储`inverse_shape_values`字段的偶数变体。
       *
       */
      AlignedVector<Number> inverse_shape_values_eo;

      /**
       * 在一个数据结构中收集在0和1点（顶点）评估的一维形状值的所有数据。排序首先是数值，然后是梯度，然后是二阶导数。
       *
       */
      std::array<AlignedVector<Number>, 2> shape_data_on_face;

      /**
       * 在一个数据结构中收集所有在0和1点（顶点）评估的一维节点形状值（由正交规则点的拉格朗日多项式定义）的数据。
       * 这个数据结构可以用来从单元格到面的正交点进行插值。
       * @note 与shape_data_on_face相反，只有值被评估。
       *
       */
      std::array<AlignedVector<Number>, 2> quadrature_data_on_face;

      /**
       * 存储子表面上形状函数的一维值。因为有两个子表面，所以存储两个变体。
       *
       */
      std::array<AlignedVector<Number>, 2> values_within_subface;

      /**
       * 存储子表面上的形状函数的一维梯度。因为有两个子表面，所以存储两个变体。
       *
       */
      std::array<AlignedVector<Number>, 2> gradients_within_subface;

      /**
       * 存储子表面上的形状函数的一维梯度。因为有两个子表面，所以存储两个变体。
       *
       */
      std::array<AlignedVector<Number>, 2> hessians_within_subface;

      /**
       * 我们存储一份用于初始化的一维正交公式的副本。
       *
       */
      Quadrature<1> quadrature;

      /**
       * 存储该元素的度数。
       *
       */
      unsigned int fe_degree;

      /**
       * 存储每个维度的正交点的数量。
       *
       */
      unsigned int n_q_points_1d;

      /**
       * 表示基函数是否在0和1的结点，即单元格的端点。
       *
       */
      bool nodal_at_cell_boundaries;

      /**
       * 存储在所有面和方向的所有正交点上评估的有限元的形状值（不利用张量-乘积结构）。
       *
       */
      Table<3, Number> shape_values_face;

      /**
       * 存储在所有面、方向和方位的所有正交点上评估的有限元的形状梯度（不利用张量-乘积结构）。
       *
       */
      Table<4, Number> shape_gradients_face;
    };



    /**
     * 该结构存储用于评估的有限元和正交公式的张量（Kronecker）积视图。它基于描述一维成分的单一或UnivariateShapeData对象集合，再加上一些关于这些成分如何组合以及在多维情况下指数如何布置的额外信息，如分层的
     *
     * -> FE_Q的词法排序。
     * @ingroup matrixfree
     *
     */
    template <typename Number>
    struct ShapeInfo
    {
      /**
       * 编码构建时检测到的元素类型。FEEvaluation将根据给定的元素类型选择最有效的算法。
       *
       */
      ElementType element_type;

      /**
       * 空的构造函数。不做任何事情。
       *
       */
      ShapeInfo();

      /**
       * 使用reinit方法初始化数据字段的构造函数。
       *
       */
      template <int dim, int dim_q>
      ShapeInfo(const Quadrature<dim_q> & quad,
                const FiniteElement<dim> &fe,
                const unsigned int        base_element = 0);

      /**
       * 初始化数据字段。接受一个一维正交公式和一个有限元作为参数，评估一维单元上的形状函数、梯度和Hessians。
       * 这个函数假定有限元是通过张量积从一维元素派生出来的，并且零中的第2个形状函数求值为1。
       *
       */
      template <int dim, int dim_q>
      void
      reinit(const Quadrature<dim_q> & quad,
             const FiniteElement<dim> &fe_dim,
             const unsigned int        base_element = 0);

      /**
       * 返回MatrixFree支持哪些种类的元素。
       *
       */
      template <int dim, int spacedim>
      static bool
      is_supported(const FiniteElement<dim, spacedim> &fe);

      /**
       * 返回单变量形状函数的数据，它定义了关于底层有限元的矢量分量
       * @p component 的张量积形状函数的维度 @p dimension 。
       *
       */
      const UnivariateShapeData<Number> &
      get_shape_data(const unsigned int dimension = 0,
                     const unsigned int component = 0) const;

      /**
       * 返回该类的内存消耗，单位为字节。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 从deal.II的单元自由度编号重新编号为DoFHandler中底层元素的FEEvaluation方案内使用的lexicographic编号。对于矢量值元素，重新编号从第一个分量的lexicographic编号开始，然后是第二个分量的所有内容，以此类推。
       *
       */
      std::vector<unsigned int> lexicographic_numbering;

      /**
       * 存储定义底层张量积有限元的单变量形状函数的数据。
       *
       */
      std::vector<UnivariateShapeData<Number>> data;

      /**
       * 允许访问给定维度和矢量分量的单变量形状函数数据。行标识尺寸，列标识矢量分量。
       *
       */
      dealii::Table<2, UnivariateShapeData<Number> *> data_access;

      /**
       * 存储空间维度的数量。
       *
       */
      unsigned int n_dimensions;

      /**
       * 存储底层矢量值有限元的矢量分量的数量。
       *
       */
      unsigned int n_components;

      /**
       * 存储一个单元的 @p dim 维度的正交点的数量。
       *
       */
      unsigned int n_q_points;

      /**
       * 存储标量元素的每个单元的DoF数量，以 @p dim
       * 维度计算。
       *
       */
      unsigned int dofs_per_component_on_cell;

      /**
       * 存储每个面的正交点数量，以 @p dim 维度。
       *
       */
      unsigned int n_q_points_face;

      /**
       * 存储一个面的正交点的数量，以 @p dim
       * 维度，用于单轴、楔形和金字塔参考单元。
       *
       */
      std::vector<unsigned int> n_q_points_faces;

      /**
       * 存储每个面的DoF数量，以 @p dim 维度计算。
       *
       */
      unsigned int dofs_per_component_on_face;

      /**
       * 对于节点位于单元格边界的节点基函数，只涉及形状函数值的面积分（DG中一阶导数的近似值）不需要加载单元格的所有自由度，而只需要加载位于面的自由度。虽然也可以在飞行中计算这些索引，但我们选择简化代码，将间接寻址存储在这个类中。
       * 第一个表的索引贯穿一个单元的面，第二个索引贯穿面的节点自由度，使用
       * @p dofs_per_face 条目。
       * 存储在这个成员变量中的索引如下。例如，考虑一个度数为3的二维元素，其自由度按词法编号如下。
       * @code
       * 12   13   14   15
       * 8    9    10   11
       * 4    5     6    7
       * 0    1     2    3
       * @endcode
       * 第一行存储索引为0的面的指数，即数字 <code>0, 4, 8,
       * 12</code>  ，第二行存储面1的指数 <code>3, 7, 11, 15</code>
       * ，第三行存储面2的指数 <code>0, 1, 2, 3</code>
       * ，最后（第四）行存储指数 <code>12, 13, 14, 15</code>
       * 。类似地，指数被存储在3D中。(注意，由于坐标系的方向性，三维中的y面使用的指数与lexicographic数字相反。)
       * @note 这个对象只在 @p nodal_at_cell_boundaries 被评估为 @p
       * true. 的情况下被填充。
       *
       */
      dealii::Table<2, unsigned int> face_to_cell_index_nodal;

      /**
       * @p face_to_cell_index_nodal
       * 为面的值的评估提供了一个快捷方式。对于Hermite类型的基函数，我们可以更进一步，也可以使用捷径来更便宜地获得导数，其中只有两层自由度对面的导数有贡献。在lexicographic排序中，与节点值相比，各自的指数是在自由度的下一个
       * "层"。这个数组存储了数值和梯度的间接寻址。
       * 第一个表索引贯穿一个单元的面，第二个表索引贯穿该面的节点自由度和导数的对，使用
       * <code>2*dofs_per_face</code> 条目。
       * 存储在这个成员变量中的指数如下。例如，考虑一个度数为3的二维元素，其自由度按词法编号如下。
       * @code
       * 20   21   22   23   24
       * 15   16   17   18   19
       * 10   11   12   13   14
       * 5    6     7    8    9
       * 0    1     2    3    4
       * @endcode
       * 第一行存储索引为0的面的值和梯度的指数，即数字<code>0,
       * 1, 5, 6, 10, 11, 15, 16, 20,
       * 21</code>，第二行存储面1的指数<code>4, 3, 9, 8, 14, 13, 19,
       * 18, 24, 23</code>，第三行存储面2的指数【 <code>0, 5, 1, 6,
       * 2, 7, 3, 8, 4, 9</code>
       * ]，最后一行（第四行）是脸部的指数<code>20, 15, 21,
       * 16, 22, 17, 23, 18, 24,
       * 19</code>。同样地，指数也是以3D方式存储的。(注意，由于坐标系的方向性，三维中的y面使用的指数与lexicographic数字相反)。
       * @note 这个对象只在 @p element_type 被评估为 @p
       * tensor_symmetric_hermite. 的情况下被填充。
       *
       */
      dealii::Table<2, unsigned int> face_to_cell_index_hermite;

      /**
       * 对于面的度数，如果一个面不在给定元素的标准方向上，基础函数的顺序就不正确。这个数据结构被用来重新排列基函数以表示正确的顺序。
       *
       */
      dealii::Table<2, unsigned int> face_orientations;

    private:
      /**
       * 检查我们在形状值中是否有对称性。在这种情况下，也要填充shape_??_eo字段。
       *
       */
      bool
      check_1d_shapes_symmetric(
        UnivariateShapeData<Number> &univariate_shape_data);

      /**
       * 检查对称的一维基础函数是否使形状值形成一个对角线矩阵，也就是说，节点与正交点是同位的。这使得专门的算法可以在评估中节省一些操作。
       *
       */
      bool
      check_1d_shapes_collocation(
        const UnivariateShapeData<Number> &univariate_shape_data) const;
    };



    // ------------------------------------------ inline functions

    template <typename Number>
    template <int dim, int dim_q>
    inline ShapeInfo<Number>::ShapeInfo(const Quadrature<dim_q> & quad,
                                        const FiniteElement<dim> &fe_in,
                                        const unsigned int base_element_number)
      : element_type(tensor_general)
      , n_dimensions(0)
      , n_components(0)
      , n_q_points(0)
      , dofs_per_component_on_cell(0)
      , n_q_points_face(0)
      , dofs_per_component_on_face(0)
    {
      reinit(quad, fe_in, base_element_number);
    }

    template <typename Number>
    inline const UnivariateShapeData<Number> &
    ShapeInfo<Number>::get_shape_data(const unsigned int dimension,
                                      const unsigned int component) const
    {
      AssertDimension(n_dimensions, data_access.size(0));
      AssertDimension(n_components, data_access.size(1));
      AssertIndexRange(dimension, n_dimensions);
      AssertIndexRange(component, n_components);
      return *(data_access(dimension, component));
    }

  } // end of namespace MatrixFreeFunctions

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


