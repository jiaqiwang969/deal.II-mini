//include/deal.II-translator/grid/tria_levels_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_tria_levels_h
#define dealii_tria_levels_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_objects.h>

#include <boost/serialization/utility.hpp>

#include <cstdint>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * 存储属于多级层次结构中某一层次的所有信息。
     * 在TriaLevel中，存储所有不依赖于维度的单元数据，例如，存储单元的细化标志的字段（单元的实际内容在其他地方声明），等等。对于非面向水平的数据，也请参见TriaObjects。
     * 还有一个字段，可能适合在这里，即材料数据（对于单元）或边界指标（对于面），但由于我们需要一个线或四边形的边界信息或材料数据，我们将它们与线和四边形一起存储，而不是与普通数据一起存储。
     * 同样地，在3D中，我们需要线和四边形的边界指标（如果相邻的两个面有不同的边界指标，我们需要知道如何细化一条线），以及单元格的材料数据。
     *
     */
    class TriaLevel
    {
    public:
      /**
       * 构造函数。              @param  三角形的dim尺寸。
       *
       */
      TriaLevel(const unsigned int dim)
        : dim(dim)
        , cells(dim)
      {}

      /**
       * 默认构造函数（Boost需要）。
       *
       */
      TriaLevel()
        : dim(numbers::invalid_unsigned_int)
        , cells(numbers::invalid_unsigned_int)
      {}

      /**
       * 三角形的尺寸。
       *
       */
      unsigned int dim;

      /**
       * @p RefinementCase<dim>::Type 单元格是否需要精化的标志
       * (RefinementCase<dim>::no_refinement).
       * 单元格的含义是与维度有关的，因此这个向量的长度也取决于维度：在一维中，这个向量的长度等于
       * @p lines 向量的长度，在二维中等于 @p quads
       * 向量的长度，等等。
       *
       */
      std::vector<std::uint8_t> refine_flags;

      /**
       * 与上面的含义相同，但指定了一个单元是否必须被粗化。
       *
       */
      std::vector<bool> coarsen_flags;


      /**
       * 一个整数，对于每一个活动单元，它存储了多少个活动单元。对于非活动单元，该值未使用，并设置为无效值。
       *
       */
      std::vector<unsigned int> active_cell_indices;

      /**
       * 每个活动单元的全局单元索引。
       *
       */
      std::vector<types::global_cell_index> global_active_cell_indices;

      /**
       * 给定级别上每个单元的全局单元索引。
       *
       */
      std::vector<types::global_cell_index> global_level_cell_indices;

      /**
       * 级别和单元格邻居的指数。惯例是，索引为 @p i
       * 的单元格的邻居存储在 $i*(2*real\_space\_dimension)$
       * 之后的字段中，例如，在一个空间维度中，单元格0的邻居存储在<tt>neighbors[0]</tt>和<tt>neighbors[1]</tt>，单元格1的邻居存储在<tt>neighbors[2]</tt>和<tt>neighbors[3]</tt>，以此类推。
       * 在邻居中，<tt>neighbors[i].first</tt>是级别，而<tt>neighbors[i].second</tt>是邻居的索引。
       * 如果一个邻居不存在（单元格在边界），<tt>level=index=-1</tt>被设置。
       * <em>  公约。 </em>  一个单元的 @p ith
       * 邻居是与该单元的 @p ith 面（二维的 @p Line ，三维的
       * @p Quad ）共享的那个。
       * 一个单元格的邻居最多具有与该单元格相同的水平，也就是说，它可能是也可能不是精炼的。
       * 在一个维度上，一个邻居可能具有小于或等于这个单元的水平的任何水平。如果它有相同的层次，它可以被精炼任意次数，但邻居指针仍然指向同一层次的单元，而邻居的子单元的邻居可以指向这个单元或其子单元。
       * 在二维和更多维度上，邻居要么在同一层次上并被细化（在这种情况下，它的子代有指向本单元或其直接子代的邻居指针），要么在同一层次上未被细化，或下一层次（在这种情况下，它的邻居指针指向本单元的母单元）。
       *
       */
      std::vector<std::pair<int, int>> neighbors;

      /**
       * 每个单元有一个整数，用来存储它属于哪个子域。这个字段最常被用于并行计算，它表示哪个处理器应该工作在/拥有给定子域编号的单元。
       * 这个号码只在活动单元上使用。
       *
       */
      std::vector<types::subdomain_id> subdomain_ids;

      /**
       * 在每一层用于并行多网格的子域ID。
       * 与subdomain_id相反，一旦网格被划分到多网格层次的下层，这个数字也会用于非活动单元上。
       *
       */
      std::vector<types::subdomain_id> level_subdomain_ids;

      /**
       * 每一对连续的单元都有一个整数，用于存储它们的父级单元的索引。
       * (我们为每对单元存储一次这一信息，因为每一次细化，无论是各向同性还是各向异性，在任何空间维度上，总是以2的倍数创建子单元，所以没有必要为每一个单元存储父单元的索引。)
       *
       */
      std::vector<int> parents;

      /**
       * 每个单元有一个bool来表示法线的方向 true:
       * 使用来自顶点的方向 false: 恢复方向。参见  @ref
       * GlossDirectionFlag  。            这仅用于codim==1的网格。
       *
       */
      std::vector<bool> direction_flags;

      /**
       * 包含线的数据和相关函数的对象
       *
       */
      TriaObjects cells;

      /**
       * 对于边，我们强制执行一个标准惯例，即相对的边应该是平行的。现在，这在大多数情况下是可以执行的，我们有代码确保如果一个网格允许这种情况发生，我们就有这个约定。我们也知道，总是有可能让相对的面具有平行的法向量。(关于这两点，请参见GridReordering类的文档中提到的ACM
       * Transactions on Mathematical Software中Agelek, Anderson, Bangerth,
       * Barth的论文)。
       * 问题是，我们原本还有一个条件，即0、2和4号面的法线指向单元格内，而其他面的法线指向外面。事实证明，这并不总是可能的。实际上，我们必须存储每个单元格的每个面的法线向量是否遵循这一惯例。如果是这样，那么这个变量就存储一个
       * @p true 值，否则就是 @p false 值。
       * 实际上，这个字段有 <code>6*n_cells</code>
       * 个元素，是每个单元格的六个面的倍数。
       * @note  只有dim=3时才需要。
       *
       */
      std::vector<unsigned char> face_orientations;

      /**
       * 每个单元格的参考单元格类型。
       * @note  仅用于dim=2和dim=3。
       *
       */
      std::vector<dealii::ReferenceCell> reference_cell;

      /**
       * 单元顶点索引的缓存（`structdim ==
       * dim`），以便更快速地检索这些频繁访问的数量。
       * 为了简化寻址，这些信息以一个单元（四边形/六面体）可能的最大顶点数为索引。
       *
       */
      std::vector<unsigned int> cell_vertex_indices_cache;

      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };


    template <class Archive>
    void
    TriaLevel::serialize(Archive &ar, const unsigned int)
    {
      ar &dim;

      ar &refine_flags &coarsen_flags;

      // do not serialize `active_cell_indices` and `vertex_indices_cache`
      // here. instead of storing them to the stream and re-reading them again
      // later, we just rebuild them in Triangulation::load()

      ar &neighbors;
      ar &subdomain_ids;
      ar &level_subdomain_ids;
      ar &parents;
      ar &direction_flags;
      ar &cells;

      if (dim == 3)
        ar &face_orientations;

      if (dim == 2 || dim == 3)
        ar &reference_cell;
    }

  } // namespace TriangulationImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


