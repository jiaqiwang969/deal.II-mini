//include/deal.II-translator/matrix_free/mapping_info_0.txt
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


#ifndef dealii_matrix_free_mapping_info_h
#define dealii_matrix_free_mapping_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/helper_functions.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 一个枚举用于识别各种类型的单元和面。最一般的类型是我们通常在FEValues上下文中计算的，但对于许多几何形状，我们可以节省大量的存储空间。
     * @ingroup matrixfree
     *
     */
    enum GeometryType : unsigned char
    {
      /**
       * 单元或面是笛卡尔式的。
       *
       */
      cartesian = 0,

      /**
       * 该单元或面可以用仿生映射来描述。
       *
       */
      affine = 1,

      /**
       * 面是平的，即面的法线因子在所有正交点上都是相同的。这种类型不分配给单元格。
       *
       */
      flat_faces = 2,

      /**
       * 没有特殊信息可用于压缩所考虑的对象的表示。
       *
       */
      general = 3
    };



    /**
     * 定义一个结构，存储所有与来自映射的评估几何体有关的缓存数据。为了支持hp-adaptivity和压缩存储（特别是对于Jacobian、JxW值和法线），不同行的存储长度可以不同。因此，它允许在单个行的数据上进行跳跃，类似于稀疏矩阵中的压缩行存储。对于不同大小的字段，我们有两种不同的起始索引。第一类偏移是单位到实数单元转换的雅各布指数（我们存储逆雅各布）、二阶导数、JxW值和法向量的指数。我们为所有这些数据结构保留单独的数组，因为用户代码可能只访问其中的一部分。在这种情况下，一个数组将以连续的顺序访问所有条目，这使得处理器很容易预取数据。如果将所有数据放在一个数组中，就需要在访问模式上有一些跨度，这对处理器来说，预测起来要复杂得多（而且确实导致了在英特尔处理器（如BroadwellEP）上不使用的数据的预取）。
     * 第二类指数是正交点的偏移量。正交点可以比其他字段压缩得更少，因此需要更长的字段。正交点指数经常被用于其他场合，如评估右手边的情况。
     * 第三个组件是来自单元格的数据描述符，称为QuadratureDescriptor，它包含正交权重和在脸部数据的情况下如何通过正交点的排列组合。后者是支持hp-adaptivity的向量，有几个数据字段用于各个正交公式。
     * @ingroup matrixfree
     *
     */
    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    struct MappingInfoStorage
    {
      static_assert(
        std::is_same<Number, typename VectorizedArrayType::value_type>::value,
        "Type of Number and of VectorizedArrayType do not match.");

      struct QuadratureDescriptor
      {
        /**
         * 构造函数。不做任何事情。
         *
         */
        QuadratureDescriptor();

        /**
         * 设置此结构的各个成员中的长度。
         *
         */
        template <int dim_q>
        void
        initialize(const Quadrature<dim_q> &quadrature,
                   const UpdateFlags update_flags_inner_faces = update_default);

        /**
         * 设置这个结构的各个成员的长度。
         *
         */
        void
        initialize(const Quadrature<1> &quadrature_1d,
                   const UpdateFlags update_flags_inner_faces = update_default);

        /**
         * 返回以字节为单位的内存消耗。
         *
         */
        std::size_t
        memory_consumption() const;

        /**
         * 应用在给定单元或面的正交点的数量。
         *
         */
        unsigned int n_q_points;

        /**
         * 应用在给定单元格或面的原始一维正交公式。
         *
         */
        Quadrature<1> quadrature_1d;

        /**
         * 应用于给定单元格或面的正交公式。
         *
         */
        Quadrature<structdim> quadrature;

        /**
         * 按尺寸分开的正交权重，用于特定情况。
         *
         */
        std::array<AlignedVector<Number>, structdim> tensor_quadrature_weights;

        /**
         * 以给定的数字格式缓存的正交权重向量（非矢量化，因为在矢量化的情况下使用该值时，向所有车道广播的成本很低）。
         *
         */
        AlignedVector<Number> quadrature_weights;

        /**
         * 对于面的正交，如果一个面不在给定元素的标准方向上，基函数的评估就没有正确的顺序。这个数据结构用于重新排列在正交点上评估的数据，以表示正确的顺序。
         *
         */
        dealii::Table<2, unsigned int> face_orientations;
      };

      /**
       * 一个描述 @p data_storage
       * 字段的布局的类，还包括一些取决于hp-context中正交点数量的数据，如内部正交公式和对不在标准方向上的面的重新索引。
       *
       */
      std::vector<QuadratureDescriptor> descriptor;

      /**
       * 应用于给定面的正交公式的集合。
       * @note 只对面填充，因为面可能是四边形或三角形的。
       *
       */
      std::vector<dealii::hp::QCollection<structdim>> q_collection;

      /**
       * 存储索引偏移到数组 @p jxw_values,   @p jacobians,   @p
       * normal_vectors
       * 和二次导数。请注意，仿生单元的长度为1，其他单元的长度等于给定单元的正交点的数量。
       *
       */
      AlignedVector<unsigned int> data_index_offsets;

      /**
       * 正交点上的雅各布行列式（在变换为非正交的情况下乘以正交权重）的存储。
       * 以 @p data_index_offsets. 为索引。
       *
       */
      AlignedVector<VectorizedArrayType> JxW_values;

      /**
       * 存储法向量。            以 @p data_index_offsets. 为索引。
       *
       */
      AlignedVector<Tensor<1, spacedim, VectorizedArrayType>> normal_vectors;

      /**
       * 存储正交点上的协变变换，即从单位单元到实数单元的变换的逆和转置的雅各布系数。
       * 索引为  @p data_index_offsets.
       * 包含两个字段，用于从两侧访问内部面，但默认情况下（单元格积分或边界积分）只填充第2个分量，忽略第一个分量。
       *
       */
      std::array<AlignedVector<Tensor<2, spacedim, VectorizedArrayType>>, 2>
        jacobians;

      /**
       * 反雅各布变换的梯度的存储。由于对称性，只需要上对角线和对角线部分。第一个索引贯穿导数，从对角线开始，然后逐行进行，即先是
       * $\partial^2/\partial x_1 \partial x_2$ ，然后是
       * $\partial^2/\partial x_1 \partial x_3$
       * ，以此类推。第二个索引是空间坐标。
       * 索引为  @p data_index_offsets.
       * 包含两个字段，用于从两侧访问内部面，但默认情况下（单元格积分或边界积分）只填充第2个分量，忽略第一个分量。
       *
       */
      std::array<
        AlignedVector<Tensor<1,
                             spacedim *(spacedim + 1) / 2,
                             Tensor<1, spacedim, VectorizedArrayType>>>,
        2>
        jacobian_gradients;

      /**
       * 存储雅各布变换乘以法向量（这代表一个经常访问的快捷方式，因此可以获得更高的性能）。
       * 以 @p data_index_offsets. 为索引。
       *
       */
      std::array<AlignedVector<Tensor<1, spacedim, VectorizedArrayType>>, 2>
        normals_times_jacobians;

      /**
       * 将某一单元的索引偏移量以实数坐标存储到正交点阵列中。请注意，直角坐标单元比非直角坐标单元（长度为
       * @p n_q_points_1d) 或面。
       *
       */
      AlignedVector<unsigned int> quadrature_point_offsets;

      /**
       * 以实坐标存储正交点，包括直角坐标单元的压缩方案，我们不需要存储所有点的完整数据。
       * 以 @p quadrature_point_offsets. 为索引。
       *
       */
      AlignedVector<Point<spacedim, VectorizedArrayType>> quadrature_points;

      /**
       * 清除除描述符向量之外的所有数据字段。
       *
       */
      void
      clear_data_fields();

      /**
       * 返回给定数量的正交点的正交指数。如果不是在hp模式下或者没有找到索引，这个函数总是返回索引0。因此，这个函数不检查给定的度数是否真的存在。
       *
       */
      unsigned int
      quad_index_from_n_q_points(const unsigned int n_q_points) const;

      /**
       * 在给定的输出流中打印该类不同结构的内存消耗的详细摘要。
       *
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &task_info) const;

      /**
       * 返回以字节为单位的内存消耗。
       *
       */
      std::size_t
      memory_consumption() const;
    };



    /**
     * 存储所有与单元格内部相关的几何数据的类，以便在无矩阵类中使用。
     * @ingroup matrixfree
     *
     */
    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfo
    {
      /**
       * 计算给定单元和面的信息。单元由层次和层次内的索引指定（如
       * CellIterator::level() 和 CellIterator::index(),
       * 给出的，以允许不同种类的迭代器，如标准DoFHandler、多网格等）在一个固定的三角结构上。此外，还给出了一个映射和几个一维正交公式。
       *
       */
      void
      initialize(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()> &             faces,
        const std::vector<unsigned int> &active_fe_index,
        const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping,
        const std::vector<dealii::hp::QCollection<dim>> &          quad,
        const UpdateFlags update_flags_cells,
        const UpdateFlags update_flags_boundary_faces,
        const UpdateFlags update_flags_inner_faces,
        const UpdateFlags update_flags_faces_by_cells);

      /**
       * 更新给定单元和面的信息，该信息是给定的`映射`类变化的结果，保持单元、正交公式和其他未知数不变。这个调用只有在
       * MappingInfo::initialize() 之前被调用过才有效。
       *
       */
      void
      update_mapping(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const FaceInfo<VectorizedArrayType::size()> &             faces,
        const std::vector<unsigned int> &active_fe_index,
        const std::shared_ptr<dealii::hp::MappingCollection<dim>> &mapping);

      /**
       * 返回初始化时检测到的给定单元的类型。
       *
       */
      GeometryType
      get_cell_type(const unsigned int cell_chunk_no) const;

      /**
       * 清除该类中的所有数据字段。
       *
       */
      void
      clear();

      /**
       * 返回该类的内存消耗，单位为字节。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 在给定的输出流中打印出这个类的不同结构的内存消耗的详细摘要。
       *
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &task_info) const;

      /**
       * 给定的更新标志用于计算单元格上的几何图形。
       *
       */
      UpdateFlags update_flags_cells;

      /**
       * 给出的用于计算边界面的几何图形的更新标志。
       *
       */
      UpdateFlags update_flags_boundary_faces;

      /**
       * 给出的用于计算内部面的几何图形的更新标志。
       *
       */
      UpdateFlags update_flags_inner_faces;

      /**
       * 给出的用于计算以单元格为中心的循环的面的几何图形的更新标志。
       *
       */
      UpdateFlags update_flags_faces_by_cells;

      /**
       * 存储一个单元是否是笛卡尔的（单元类型0），是否有恒定的变换数据（雅各布）（单元类型1），或者是一般的（单元类型3）。类型2仅用于面，没有单元格被分配此值。
       *
       */
      std::vector<GeometryType> cell_type;

      /**
       * 存储一个面（以及与该面相邻的两个单元）是否是笛卡尔的（面类型0），是否代表一个仿生情况（面类型1），是否是一个平坦的面，法向量在整个面中是相同的（面类型2），或者是一般的（面类型3）。
       *
       */
      std::vector<GeometryType> face_type;

      /**
       * 单元的数据缓存。
       *
       */
      std::vector<MappingInfoStorage<dim, dim, Number, VectorizedArrayType>>
        cell_data;

      /**
       * 面的数据缓存。
       *
       */
      std::vector<MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>>
        face_data;

      /**
       * 面的数据缓存--与细胞相关的拓扑结构，遵循细胞类型的
       * @p cell_type 变量。
       *
       */
      std::vector<MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>>
        face_data_by_cells;

      /**
       * 指向底层 hp::MappingCollection 对象的指针。
       *
       */
      std::shared_ptr<dealii::hp::MappingCollection<dim>> mapping_collection;

      /**
       * 指向mapping_collection的第一个条目的指针。
       *
       */
      SmartPointer<const Mapping<dim>> mapping;

      /**
       * 与每个正交和活动正交索引相关的参考单元类型。
       *
       */
      std::vector<std::vector<dealii::ReferenceCell>> reference_cell_types;

      /**
       * 内部函数，用于计算映射是MappingQ且每个槽使用一个正交公式的情况下的几何图形（非hp情况）。该方法使用无矩阵框架本身的快速算子评估技术计算来自底层单元格正交点的所有数据，也就是说，它使用单元格几何的多项式描述（在第一步计算），然后根据这些信息计算所有雅各布和法向量。这种优化的方法比通过FEValues和FEFaceValues要快得多，特别是当涉及到几个不同的正交公式时，消耗的内存也更少。
       * @param  tria 用于设置的三角形  @param  cells
       * 要处理的三角形的实际单元，以级别和级别内索引的元组形式给出，在类的主初始化中使用
       * @param  faces
       * 从面到单元的连接描述，在MatrixFree类中填写。
       *
       */
      void
      compute_mapping_q(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
          &faces);

      /**
       * 计算给定单元中的信息，在初始化中调用。
       *
       */
      void
      initialize_cells(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<unsigned int> &         active_fe_index,
        const dealii::hp::MappingCollection<dim> &mapping);

      /**
       * 计算给定面的信息，在initialize中调用。
       *
       */
      void
      initialize_faces(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const std::vector<FaceToCellTopology<VectorizedArrayType::size()>>
          &                                       faces,
        const std::vector<unsigned int> &         active_fe_index,
        const dealii::hp::MappingCollection<dim> &mapping);

      /**
       * 计算给定面的信息，在初始化时调用。
       *
       */
      void
      initialize_faces_by_cells(
        const dealii::Triangulation<dim> &                        tria,
        const std::vector<std::pair<unsigned int, unsigned int>> &cells,
        const dealii::hp::MappingCollection<dim> &                mapping);

      /**
       * 帮助函数，确定在内部函数中必须设置哪些更新标志，以便按照用户的要求初始化所有数据。
       *
       */
      static UpdateFlags
      compute_update_flags(
        const UpdateFlags                                update_flags,
        const std::vector<dealii::hp::QCollection<dim>> &quad =
          std::vector<dealii::hp::QCollection<dim>>());
    };



    /**
     * 一个帮助类，用于从映射信息中提取单元或面孔数据，以便在FEEvaluationBase中使用。
     *
     */
    template <int, typename, bool, typename>
    struct MappingInfoCellsOrFaces;

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, false, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim, dim, Number, VectorizedArrayType> *
      get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
          const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.cell_data.size());
        return &mapping_info.cell_data[quad_no];
      }
    };

    template <int dim, typename Number, typename VectorizedArrayType>
    struct MappingInfoCellsOrFaces<dim, Number, true, VectorizedArrayType>
    {
      static const MappingInfoStorage<dim - 1, dim, Number, VectorizedArrayType>
        *
        get(const MappingInfo<dim, Number, VectorizedArrayType> &mapping_info,
            const unsigned int                                   quad_no)
      {
        AssertIndexRange(quad_no, mapping_info.face_data.size());
        return &mapping_info.face_data[quad_no];
      }
    };



    /**
     * 一个用于比较浮点数组的类（例如： std::vectors,
     * Tensor<1,dim>，等等）。这个类的想法是，如果两个数组在一个给定的公差范围内是相同的，就认为它们是相等的。我们在给定数组的一个
     * std::map<>
     * 内使用这个比较器类。请注意，这个比较器并不满足人们通常想要的所有数学属性（考虑例如数字a=0,
     * b=0.1,
     * c=0.2，公差为0.15；该运算器给出a<c，但既不满足a<b?也不满足b<c?）。这在这个类的用例中不是问题，但在其他情况下使用它时要小心。
     *
     */
    template <typename Number,
              typename VectorizedArrayType = VectorizedArray<Number>>
    struct FPArrayComparator
    {
      FPArrayComparator(const Number scaling);

      /**
       * 比较两个数字向量（不一定是相同长度的）。
       *
       */
      bool
      operator()(const std::vector<Number> &v1,
                 const std::vector<Number> &v2) const;

      /**
       * 比较两个向量的数组（存储为张量，以避免对齐问题）。
       *
       */
      bool
      operator()(
        const Tensor<1, VectorizedArrayType::size(), Number> &t1,
        const Tensor<1, VectorizedArrayType::size(), Number> &t2) const;

      /**
       * 比较两个向量数组的秩-1张量（存储为张量，以避免对齐问题）。
       *
       */
      template <int dim>
      bool
      operator()(
        const Tensor<1, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t1,
        const Tensor<1, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t2) const;

      /**
       * 比较向量数组的两个秩-2张量（存储为张量，以避免对齐问题）。
       *
       */
      template <int dim>
      bool
      operator()(
        const Tensor<2, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t1,
        const Tensor<2, dim, Tensor<1, VectorizedArrayType::size(), Number>>
          &t2) const;

      /**
       * 比较两个数组的张量。
       *
       */
      template <int dim>
      bool
      operator()(const std::array<Tensor<2, dim, Number>, dim + 1> &t1,
                 const std::array<Tensor<2, dim, Number>, dim + 1> &t2) const;

      Number tolerance;
    };



     /* ------------------- inline functions ----------------------------- */ 

    template <int structdim,
              int spacedim,
              typename Number,
              typename VectorizedArrayType>
    inline unsigned int
    MappingInfoStorage<structdim, spacedim, Number, VectorizedArrayType>::
      quad_index_from_n_q_points(const unsigned int n_q_points) const
    {
      for (unsigned int i = 0; i < descriptor.size(); ++i)
        if (n_q_points == descriptor[i].n_q_points)
          return i;
      return 0;
    }



    template <int dim, typename Number, typename VectorizedArrayType>
    inline GeometryType
    MappingInfo<dim, Number, VectorizedArrayType>::get_cell_type(
      const unsigned int cell_no) const
    {
      AssertIndexRange(cell_no, cell_type.size());
      return cell_type[cell_no];
    }

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


