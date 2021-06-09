//include/deal.II-translator/grid/tria_description_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_grid_construction_utilities_h
#define dealii_grid_construction_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/la_parallel_vector.h>


DEAL_II_NAMESPACE_OPEN

 /*------------------------------------------------------------------------*/ 

/**
 * CellData类（以及相关的SubCellData类）用于在通过
 * Triangulation::create_triangulation().
 * 创建三角图时提供全面但最低限度的单元描述。
 * 具体来说，每个CellData对象
 *
 * - 描述三角测量中的一个单元
 *
 * 当作为SubCellData类的成员时，这个结构也被用来表示面和边的数据。在这种情况下，对象的模板参数 @p structdim 将小于三角形的尺寸 @p dim 。如果是这样，那么#vertices数组代表了传递给 Triangulation::create_triangulation(). 的一个单元格的一个面或边的顶点的索引。此外，对于面来说，材料id没有任何意义， @p material_id 字段被重新用来存储一个 @p boundary_id 来代替指定面或边属于边界的哪一部分（见 @ref GlossBoundaryIndicator  "边界id的词汇条"
 * ）。 一个显示该类如何使用的例子是在  step-14  的
 * <code>create_coarse_grid()</code>
 * 函数中。在实现GridGenerator命名空间的功能时，还有很多用例。
 *
 *
 * @ingroup grid
 *
 *
 */
template <int structdim>
struct CellData
{
  /**
   * 这个单元的顶点的索引。这些指数对应于传递给
   * Triangulation::create_triangulation().
   * 的顶点位置向量中的条目。
   *
   */
  std::vector<unsigned int> vertices;

  /**
   * 此单元的材料或边界指标。
   * 这个字段是一个联合体，存储<i>either</i>边界或材料ID，取决于当前对象是用来描述一个单元（在CellData对象的矢量中）还是一个面或边（作为SubCellData对象的一部分）。
   *
   */
  union
  {
    /**
     * 被描述的单元格的材质ID。关于如何使用这个字段的例子，请参见CellData类的文档。
     * 这个变量只有在当前对象被用来描述一个单元时才能使用，即如果
     * @p structdim 等于一个三角形的尺寸 @p dim 。
     *
     */
    types::material_id material_id;

    /**
     * 被描述的面或边的边界ID。关于如何使用这个字段的例子，请参见CellData类的文档。
     * 这个变量只能在当前对象用于描述一个面或边的情况下使用，也就是说，如果
     * @p structdim 小于三角形的尺寸 @p dim
     * 。在这种情况下，这个变量所属的CellData对象将是SubCellData对象的一部分。
     *
     */
    types::boundary_id boundary_id;
  };

  /**
   * 该对象的歧管标识符。这个标识符应该被用来识别这个对象所属的流形，这个对象将从流形中收集关于细化后如何添加点的信息。
   * 关于如何使用这个字段的例子，请参见CellData类的文档。
   *
   */
  types::manifold_id manifold_id;

  /**
   * 默认构造函数。将成员变量设置为以下值。
   *
   *
   *
   *
   *
   *
   * -顶点指数为无效值
   *
   *
   *
   * - 边界或材料ID为零（边界或材料ID的默认值）。
   *
   *
   *
   *
   *
   *
   * - 歧管ID为 numbers::flat_manifold_id 。
   *
   */
  CellData(
    const unsigned int n_vertices = GeometryInfo<structdim>::vertices_per_cell);

  /**
   * 比较运算符。
   *
   */
  bool
  operator==(const CellData<structdim> &other) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  static_assert(structdim > 0,
                "The class CellData can only be used for structdim>0.");
};



/**
 * SubCellData类在通过 Triangulation::create_triangulation().
 * 创建三角形时用于描述网格边界的面和边的信息，它包含描述边界边和边界四边形的成员变量。
 * 该类没有模板参数，既可用于描述2D的边界边（在这种情况下，
 * @p boundary_quads
 * 成员变量的内容被忽略），也可用于描述3D的边界边和面（在这种情况下，
 * @p boundary_lines 和 @p boundary_quads
 * 成员都可以使用）。它也被用作1d中
 * Triangulation::create_triangulation()
 * 的参数，其中当前类型的对象的内容被简单地忽略。
 * @p boundary_lines 和 @p boundary_quads
 * 向量中的每个条目需要对应于由传递给
 * Triangulation::create_triangulation().
 * 的CellData对象向量所描述的单元格的边或四边，即每个条目中存储的顶点指数需要对应于具有相同顶点指数集的三角形的边或面，且顺序相同。对于这些边界边或四边形，可以设置
 * CellData::boundary_id 和 CellData::manifold_id. 中的一个或两个。
 * 也有一些用例，人们可能想设置<i>interior</i>边或面的流形ID。这种由其顶点指数识别的面也可能出现在
 * @p boundary_lines 和 @p boundary_quads
 * 向量中（尽管这些成员变量的名称不同）。然而，这显然不允许设置边界ID（因为该对象实际上不是边界的一部分）。因此，为了有效，内部边缘或面的
 * CellData::boundary_id 需要等于 numbers::internal_face_boundary_id.  。
 *
 *
 * @ingroup grid
 *
 *
 */
struct SubCellData
{
  /**
   * 一个CellData<1>对象的向量，描述2D或3D三角形的边缘的边界和流形信息。
   * 这个向量不能用于创建1d三角形。
   *
   */
  std::vector<CellData<1>> boundary_lines;

  /**
   * 一个CellData<2>对象的向量，描述三维三角形的四边形的边界和流形信息。
   * 这个向量不能用于创建1d或2d三角形。
   *
   */
  std::vector<CellData<2>> boundary_quads;

  /**
   * 判断上述可能不会在给定维度中使用的成员变量是否真的为空。换句话说，当
   * @p dim 等于1时，此函数返回 @p boundary_lines 和 @p
   * boundary_quads 是否都是空向量，当 @p dim 等于2时， @p
   * boundary_quads 向量是否为空。
   *
   */
  bool
  check_consistency(const unsigned int dim) const;
};


template <int structdim>
template <class Archive>
void
CellData<structdim>::serialize(Archive &ar, const unsigned int  /*version*/ )
{
  ar &vertices;
  ar &material_id;
  ar &boundary_id;
  ar &manifold_id;
}

/**
 * 一个专门用于结构描述的命名空间，可以在
 * Triangulation::create_triangulation(). 中使用。
 *
 *
 */
namespace TriangulationDescription
{
  /**
   * Triangulations的配置标志。  设置可以用位法OR来组合。
   *
   */
  enum Settings
  {
    /**
     * 默认设置，其他选项被禁用。
     *
     */
    default_setting = 0x0,
    /**
     * 这个标志需要被设置以使用几何多网格功能。这个选项需要额外的计算和通信。
     *
     */
    construct_multigrid_hierarchy = 0x1
  };

  /**
   * 每个本地相关单元所需的信息，存储在Description中，并在构建三角结构时使用。该结构存储了单元ID、子域ID和水平子域ID，以及与manifold_id和boundary_id相关的信息。
   * @note  与 dealii::CellData,
   * 类似，该结构也存储单元格的信息。然而，与
   * dealii::CellData,
   * 不同的是，它还存储一个唯一的id、分区信息以及与单元格面和边相关的信息。
   *
   */
  template <int dim>
  struct CellData
  {
    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int  /*version*/ );

    /**
     * 比较运算符。
     *
     */
    bool
    operator==(const CellData<dim> &other) const;

    /**
     * 单元的唯一CellID。
     *
     */
    CellId::binary_type id;

    /**
     * 单元的subdomain_id。
     *
     */
    types::subdomain_id subdomain_id;

    /**
     * 该单元格的level_subdomain_id。
     *
     */
    types::subdomain_id level_subdomain_id;

    /**
     * 单元的Manifold id。
     *
     */
    types::manifold_id manifold_id;

    /**
     * 该单元的所有行的Multifold id。
     * @note  仅用于2D和3D。
     *
     */
    std::array<types::manifold_id, GeometryInfo<dim>::lines_per_cell>
      manifold_line_ids;

    /**
     * 单元中所有面的四边形的集合体ID。
     * @note  仅用于三维。
     *
     */
    std::array<types::manifold_id,
               dim == 1 ? 1 : GeometryInfo<3>::quads_per_cell>
      manifold_quad_ids;

    /**
     * 单元的所有非内部面的面号和边界ID的列表。
     *
     */
    std::vector<std::pair<unsigned int, types::boundary_id>> boundary_ids;
  };

  /**
   * 在 Triangulation::create_triangulation(). 中使用的数据。
   *
   */
  template <int dim, int spacedim>
  struct Description
  {
    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int  /*version*/ );

    /**
     * 比较运算符。
     *
     */
    bool
    operator==(const Description<dim, spacedim> &other) const;

    /**
     * 本地相关的粗网格三角测量的单元格。
     *
     */
    std::vector<dealii::CellData<dim>> coarse_cells;

    /**
     * 本地相关的粗网格三角结构的顶点。
     *
     */
    std::vector<Point<spacedim>> coarse_cell_vertices;

    /**
     * 为每个本地相关的粗网格单元提供相应的全局 @ref
     * GlossCoarseCellId 的列表。
     *
     */
    std::vector<types::coarse_cell_id> coarse_cell_index_to_coarse_cell_id;

    /**
     * cell_infos[i]包含第i层的每个本地相关单元的CellData。
     *
     */
    std::vector<std::vector<CellData<dim>>> cell_infos;

    /**
     * 用于创建此结构的MPI通信器。它将与Triangulation内部的通信器进行比较，如果不匹配则抛出断言。
     * @note 请注意这是必要的，因为 parallel::TriangulationBase
     * 中的通信器是常量，在构造函数被调用后不能被改变。
     *
     */
    MPI_Comm comm;

    /**
     * 在构建三角形时要使用的属性。
     *
     */
    Settings settings;

    /**
     * 网格平滑的类型。
     *
     */
    typename Triangulation<dim, spacedim>::MeshSmoothing smoothing;
  };


  /**
   * 一个用于 TriangulationDescription::Description
   * 实用函数的命名空间。
   * @ingroup TriangulationDescription
   *
   */
  namespace Utilities
  {
    /**
     * 从一个给定的分区三角形`tria`和一个指定的过程中构造
     * TriangulationDescription::Description 。
     * 输入的三角图可以是类型为 dealii::Triangulation
     * 的串行三角图，它已被着色（子域_id和/或level_subdomain_id已被设置）或类型为
     * dealii::parallel::distributed::Triangulation,
     * 的分布式三角图，其中分区被采用而不被改变。
     * @param  tria 分布式输入三角图。      @param  comm
     * 要使用的MPI_Communicator。在
     * dealii::parallel::distributed::Triangulation,
     * 的情况下，通信器必须匹配。      @param  settings
     * 参见设置枚举器的描述。      @param  my_rank_in 构造
     * 指定等级的描述（仅适用于已被
     * GridToold::partition_triangulation()).
     * 等函数分割的序列三角图）  @return
     * 用于设置三角图的描述。
     * @note
     * 如果在设置中设置了construct_multigrid_hierarchy，则必须用limit_level_difference_at_vertices设置源三角形。
     *
     */
    template <int dim, int spacedim = dim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const dealii::Triangulation<dim, spacedim> &tria,
      const MPI_Comm &                            comm,
      const TriangulationDescription::Settings    settings =
        TriangulationDescription::Settings::default_setting,
      const unsigned int my_rank_in = numbers::invalid_unsigned_int);

    /**
     * 与上述函数类似，但活动单元的所有者由单元向量提供（另见
     * parallel::TriangulationBase::global_active_cell_index_partitioner() 和
     * CellAccessor::global_active_cell_index()).
     * 该函数允许重新划分分布式三角测量对象。
     * @note  从矢量中提取通讯器  @p partition.  。
     * @note  三角测量 @p tria 可以在 @p partitioner.
     * 的通信器的子通信器上设置，所有不属于该子通信器的进程需要用特殊目的通信器MPI_COMM_NULL设置本地三角测量。
     * @note 目前没有构建多网格层次，因为 @p partition
     * 只描述了活动层次的划分。
     *
     */
    template <int dim, int spacedim>
    Description<dim, spacedim>
    create_description_from_triangulation(
      const Triangulation<dim, spacedim> &              tria,
      const LinearAlgebra::distributed::Vector<double> &partition);


    /**
     * 构建一个 TriangulationDescription::Description.
     * 与上面的函数不同，这个函数也负责创建一个序列三角形，并负责其分区（通过调用提供的
     * `std::function`
     * 对象）。在内部，只有选定的进程（每一个n-th/每一个大小为group_size的组的根）才会为其组中的所有进程创建一个串行三角形和
     * TriangulationDescription::Description ，这是被通报的。
     * @note
     * 合理的组大小是一个NUMA域的大小或一个计算节点的大小。
     * @param  serial_grid_generator 一个创建串行三角形的函数。
     * @param  serial_grid_partitioner
     * 一个可以分割串行三角图的函数，即设置活动单元的sudomain_ids。
     * 该函数的第一个参数是一个串行三角图，第二个参数是MPI通信器，第三个参数是组的大小。
     * @param  comm MPI communicator。      @param  group_size
     * 每个组的大小。      @param  smoothing 网格平滑类型。
     * @param  setting 参见设置枚举器的描述。      @return
     * 用于设置三角测量的描述。
     * @note  如果在设置中设置了construct_multigrid_hierarchy， @p
     * smoothing
     * 参数会被扩展为limit_level_difference_at_vertices标志。
     *
     */
    template <int dim, int spacedim = dim>
    Description<dim, spacedim>
    create_description_from_triangulation_in_groups(
      const std::function<void(dealii::Triangulation<dim, spacedim> &)>
        &                                            serial_grid_generator,
      const std::function<void(dealii::Triangulation<dim, spacedim> &,
                               const MPI_Comm &,
                               const unsigned int)> &serial_grid_partitioner,
      const MPI_Comm &                               comm,
      const int                                      group_size = 1,
      const typename Triangulation<dim, spacedim>::MeshSmoothing smoothing =
        dealii::Triangulation<dim, spacedim>::none,
      const TriangulationDescription::Settings setting =
        TriangulationDescription::Settings::default_setting);

  } // namespace Utilities



  template <int dim>
  template <class Archive>
  void
  CellData<dim>::serialize(Archive &ar, const unsigned int  /*version*/ )
  {
    ar &id;
    ar &subdomain_id;
    ar &level_subdomain_id;
    ar &manifold_id;
    if (dim >= 2)
      ar &manifold_line_ids;
    if (dim >= 3)
      ar &manifold_quad_ids;
    ar &boundary_ids;
  }


  template <int dim, int spacedim>
  template <class Archive>
  void
  Description<dim, spacedim>::serialize(Archive &ar,
                                        const unsigned int  /*version*/ )
  {
    ar &coarse_cells;
    ar &coarse_cell_vertices;
    ar &coarse_cell_index_to_coarse_cell_id;
    ar &cell_infos;
    ar &settings;
    ar &smoothing;
  }



  template <int dim>
  bool
  CellData<dim>::operator==(const CellData<dim> &other) const
  {
    if (this->id != other.id)
      return false;
    if (this->subdomain_id != other.subdomain_id)
      return false;
    if (this->level_subdomain_id != other.level_subdomain_id)
      return false;
    if (this->manifold_id != other.manifold_id)
      return false;
    if (dim >= 2 && this->manifold_line_ids != other.manifold_line_ids)
      return false;
    if (dim >= 3 && this->manifold_quad_ids != other.manifold_quad_ids)
      return false;
    if (this->boundary_ids != other.boundary_ids)
      return false;

    return true;
  }



  template <int dim, int spacedim>
  bool
  Description<dim, spacedim>::
  operator==(const Description<dim, spacedim> &other) const
  {
    if (this->coarse_cells != other.coarse_cells)
      return false;
    if (this->coarse_cell_vertices != other.coarse_cell_vertices)
      return false;
    if (this->coarse_cell_index_to_coarse_cell_id !=
        other.coarse_cell_index_to_coarse_cell_id)
      return false;
    if (this->cell_infos != other.cell_infos)
      return false;
    if (this->settings != other.settings)
      return false;
    if (this->smoothing != other.smoothing)
      return false;

    return true;
  }
} // namespace TriangulationDescription


DEAL_II_NAMESPACE_CLOSE

#endif


