//include/deal.II-translator/Matrix_free/face_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_face_info_h
#define dealii_matrix_free_face_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/table.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 数据类型，用于建立用于面积分矢量化的批次信息。面的批次的设置是独立于单元的，因此，我们必须存储与单元索引的关系，以便访问自由度。
     * 内部面由两个相邻的单元存储，我们将其标记为面的
     * "内部 "和 "外部
     * "侧。存储在MappingInfo中的法线向量只存储一次，是
     * "内部 "面的单元的外侧法线，而 "外部
     * "面的符号则相反。
     * 这个数据字段被存储为所有参与计算的面的一个向量。为了避免在内存表示中出现空隙，四个
     * "char
     * "变量被放在彼此旁边，这与大多数架构上的无符号整数所占的大小相同。
     *
     */
    template <int vectorization_width>
    struct FaceToCellTopology
    {
      /**
       * 当前面批中的面的指数，与面的逻辑 "内部
       * "侧的单元格的数字相比，该面的逻辑 "内部 "侧与
       * FEEvaluation::get_normal_vector(). 的方向对齐。
       *
       */
      std::array<unsigned int, vectorization_width> cells_interior;

      /**
       * 当前面批中的面的指数与面的逻辑 "外部
       * "侧的单元格数量相比，该面与
       * FEEvaluation::get_normal_vector().
       * 的方向相反。注意，内部和外部面的区分是纯逻辑的，仅指法线方向。在问题的实际离散化过程中，离散化通常需要确保内部和外部面得到正确的处理，例如上风通量。
       * 对于边界面，数字被设置为 `numbers::invalid_unsigned_int`.
       * 。
       *
       */
      std::array<unsigned int, vectorization_width> cells_exterior;

      /**
       * 在面的 "外部 "的单元格内，0和
       * GeometryInfo::faces_per_cell 之间的面的索引。
       * 对于一个边界面，这个数据字段存储了边界ID。
       *
       */
      types::boundary_id exterior_face_no;

      /**
       * 在面的 "内部 "一侧的单元格中，介于0和
       * GeometryInfo::faces_per_cell 之间的面的索引。
       *
       */
      unsigned char interior_face_no;

      /**
       * 对于自适应细化网格，面的外侧的单元格可能比内侧的细化程度低。这个指数表示外侧的可能的子面指数。
       *
       */
      unsigned char subface_index;

      /**
       * 在三维中，与一个面相邻的两个单元中的一个可能使用与标准方向不同的方向（也称为面朝向、面翻转和面旋转）。这个变量在第一个自由位中存储了本批面的面朝向、面翻转和面旋转（用于内部或外部的一个面）的值。如果内部单元有非标准方向，则第四位为一。
       * @note
       * 与库中其他地方相比，面的方向位（第一位）是翻转的。
       *
       */
      unsigned char face_orientation;

      /**
       * 给定面的参考单元类型。0代表直线或四边形，1代表三角形。
       *
       */
      unsigned char face_type;

      /**
       * 返回当前数据结构的内存消耗。
       *
       */
      std::size_t
      memory_consumption() const
      {
        return sizeof(*this);
      }
    };



    /**
     * 保存面和单元之间的连接性的数据结构。
     *
     */
    template <int vectorization_width>
    struct FaceInfo
    {
      /**
       * 清除所有数据字段，使其处于类似于调用默认构造函数后的状态。
       *
       */
      void
      clear()
      {
        faces = std::vector<FaceToCellTopology<vectorization_width>>();
        cell_and_face_to_plain_faces.reinit(TableIndices<3>(0, 0, 0));
        cell_and_face_boundary_id.reinit(TableIndices<3>(0, 0, 0));
      }

      /**
       * 返回当前数据结构的内存消耗。
       *
       */
      std::size_t
      memory_consumption() const
      {
        return sizeof(faces) +
               cell_and_face_to_plain_faces.memory_consumption() +
               cell_and_face_boundary_id.memory_consumption();
      }

      /**
       * 内部面的矢量存储，链接到矢量单元存储中的两个单元。
       *
       */
      std::vector<FaceToCellTopology<vectorization_width>> faces;

      /**
       * 该表将宏观单元编号、单元内的面的索引和向量化的单元批内的索引三者转化为
       * @p faces 数组内的索引。
       *
       */
      ::dealii::Table<3, unsigned int> cell_and_face_to_plain_faces;

      /**
       * 使用与cell_and_face_to_plain_faces数据结构相同的索引，以矢量化格式存储面的边界ID
       *
       */
      ::dealii::Table<3, types::boundary_id> cell_and_face_boundary_id;
    };
  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


