//include/deal.II-translator/grid/grid_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#ifndef dealii_grid_tools_h
#  define dealii_grid_tools_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/bounding_box.h>
#  include <deal.II/base/geometry_info.h>
#  include <deal.II/base/std_cxx17/optional.h>

#  include <deal.II/boost_adaptors/bounding_box.h>

#  include <deal.II/distributed/shared_tria.h>

#  include <deal.II/dofs/dof_handler.h>

#  include <deal.II/fe/fe_values.h>
#  include <deal.II/fe/mapping.h>
#  include <deal.II/fe/mapping_q1.h>

#  include <deal.II/grid/manifold.h>
#  include <deal.II/grid/tria.h>
#  include <deal.II/grid/tria_accessor.h>
#  include <deal.II/grid/tria_iterator.h>

#  include <deal.II/hp/dof_handler.h>

#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/la_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/sparsity_tools.h>
#  include <deal.II/lac/trilinos_vector.h>

#  include <deal.II/numerics/rtree.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <boost/archive/binary_iarchive.hpp>
#  include <boost/archive/binary_oarchive.hpp>
#  include <boost/geometry/index/rtree.hpp>
#  include <boost/random/mersenne_twister.hpp>
#  include <boost/serialization/array.hpp>
#  include <boost/serialization/vector.hpp>

#  ifdef DEAL_II_WITH_ZLIB
#    include <boost/iostreams/device/back_inserter.hpp>
#    include <boost/iostreams/filter/gzip.hpp>
#    include <boost/iostreams/filtering_stream.hpp>
#    include <boost/iostreams/stream.hpp>
#  endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#  include <bitset>
#  include <list>
#  include <set>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
namespace parallel
{
  namespace distributed
  {
    template <int, int>
    class Triangulation;
  }
} // namespace parallel

namespace hp
{
  template <int, int>
  class MappingCollection;
}

class SparsityPattern;
#  endif

namespace internal
{
  template <int dim, int spacedim, class MeshType>
  class ActiveCellIterator
  {
  public:
#  ifndef _MSC_VER
    using type = typename MeshType::active_cell_iterator;
#  else
    using type = TriaActiveIterator<dealii::CellAccessor<dim, spacedim>>;
#  endif
  };

#  ifdef _MSC_VER
  template <int dim, int spacedim>
  class ActiveCellIterator<dim, spacedim, dealii::DoFHandler<dim, spacedim>>
  {
  public:
    using type =
      TriaActiveIterator<dealii::DoFCellAccessor<dim, spacedim, false>>;
  };
#  endif
} // namespace internal

/**
 * 这个命名空间是工作在三角形上的算法的集合，比如移动或旋转三角形，但也可以找到包含给定点的单元。更多信息请参见各个函数的描述。
 *
 *
 * @ingroup grid
 *
 *
 */
namespace GridTools
{
  template <int dim, int spacedim>
  class Cache;

  /**
   * @name  关于网格和单元的信息
   *
   */
   /*@{*/ 

  /**
   * 返回一个三角结构的直径。直径的计算只使用顶点，也就是说，如果由于高阶映射的原因，直径应该大于边界顶点之间的最大距离，那么这个函数将不会捕捉到这个。
   *
   */
  template <int dim, int spacedim>
  double
  diameter(const Triangulation<dim, spacedim> &tria);

  /**
   * 计算三角形的体积（即二维度量）。我们使用积分
   * $\sum_K \int_K 1
   * \; dx$ 来计算这个度量，其中 $K$
   * 是给定三角结构的单元。该积分通过正交来逼近，为此我们需要映射参数。
   * 如果三角形是一个嵌入到高维空间中的二维空间，那么返回的值就是二维度量。例如，对于三维空间中的二维三角形，返回的值是所描述的曲面的面积。这显然是有意义的，因为如果dim
   * @<  spacedim，那么dim-dimensional triangulation的spacedim-dimensional
   * measure将总是零。    这个函数也适用于
   * parallel::distributed::Triangulation,
   * 类型的对象，在这种情况下，这个函数是一个集体操作。
   * @param  tria 三角化。    @param  mapping
   * 一个可选的参数，用于表示在描述单元格是由直线还是曲线面所包围时应该使用的映射。默认情况下使用
   * $Q_1$ 映射，它对应于单元格的直线边界。    @return
   * 三角形所描述的域的dim-dimensional度量，如上所述。
   *
   */
  template <int dim, int spacedim>
  double
  volume(const Triangulation<dim, spacedim> &tria,
         const Mapping<dim, spacedim> &      mapping =
           (ReferenceCells::get_hypercube<dim>()
              .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 返回一个三角形的最小活动单元的直径近似值。参见
   * step-24 中关于这个函数的使用实例。
   * 请注意，即使你传递了一个非琐碎的映射，返回的值也只是使用三角结构的顶点信息来计算的，可能是通过映射来转换的。虽然这在大多数情况下是准确的，但当三角结构包含非常扭曲的单元时，它可能无法给出正确的结果。
   *
   */
  template <int dim, int spacedim>
  double
  minimal_cell_diameter(
    const Triangulation<dim, spacedim> &triangulation,
    const Mapping<dim, spacedim> &      mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 返回一个三角形的最大活动单元的直径近似值。
   * 请注意，即使你向这个函数传递了一个非琐碎的映射，返回的值也只是使用三角结构的顶点信息来计算的，可能会被映射转化。虽然这在大多数情况下是准确的，但当三角结构包含非常扭曲的单元时，它可能无法给出正确的结果。
   *
   */
  template <int dim, int spacedim>
  double
  maximal_cell_diameter(
    const Triangulation<dim, spacedim> &triangulation,
    const Mapping<dim, spacedim> &      mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 给出一个顶点列表（通常使用 Triangulation::get_vertices)
   * 获得）作为第一个参数，以及一个表征单个单元的顶点索引列表作为第二个参数，返回这个单元的度量（面积、体积）。如果这是一个真实的单元格，那么你可以使用
   * <code>cell-@>measure()</code>
   * 得到同样的结果，但这个函数也适用于不存在的单元格，只是你通过从列表中命名它的顶点来编造它。
   * @deprecated
   * 使用更通用的函数，该函数需要一个ArrayView来代替。
   *
   */
  template <int dim>
  DEAL_II_DEPRECATED double
  cell_measure(
    const std::vector<Point<dim>> &all_vertices,
    const unsigned int (&vertex_indices)[GeometryInfo<dim>::vertices_per_cell]);

  /**
   * 给出一个顶点列表（通常使用 Triangulation::get_vertices())
   * 获得）作为第一个参数，以及一个表征单个单元的顶点索引列表作为第二个参数，返回该单元的度量（面积、体积）。如果这是一个真实的单元格，那么你可以用
   * <code>cell-@>measure()</code>
   * 得到同样的结果，但这个函数也适用于不存在的单元格，只是你通过从列表中命名它的顶点来编造它。
   * 参数 @p vertex_indices 被期望有
   * GeometryInfo<dim>::vertices_per_cell 项。一个 std::vector
   * 可以隐式转换为ArrayView，所以它可以直接传递。更多信息请参见ArrayView类。
   * @note  这个函数只对二维度为零的对象实现。
   *
   */
  template <int dim>
  double
  cell_measure(const std::vector<Point<dim>> &      all_vertices,
               const ArrayView<const unsigned int> &vertex_indices);

  /**
   * 这个函数通过对代表四边形或六面体单元的 $p_\text{real}
   * = A p_\text{unit} + b $
   * 顶点进行最小二乘法拟合，计算出从单位坐标到实坐标形式的仿生近似图，其维度为`spacedim`。结果是以矩阵<i>A</i>作为第一个参数，向量<i>b</i>描述平面到原点的距离，以一对方式返回。
   * 对于任何有效的网格单元，其几何形状不是退化的，这个操作的结果是唯一的仿生映射，即使在双/三线或高阶映射的实际转换可能是单数的情况下。如果从单元到实际单元的转换确实是仿生的，例如在一维或二维/三维的笛卡尔和仿生（平行四边形）网格中，其结果是精确的。
   * 这种近似是函数
   * TriaAccessor::real_to_unit_cell_affine_approximation() 功能的基础。
   * 对于单元格的精确变换，使用
   * Mapping::transform_real_to_unit_cell(). 。
   *
   */
  template <int dim, int spacedim>
  std::pair<DerivativeForm<1, dim, spacedim>, Tensor<1, spacedim>>
  affine_cell_approximation(const ArrayView<const Point<spacedim>> &vertices);

  /**
   * 计算所有本地拥有的活动单元的长宽比度量，并填充一个每个单元有一个条目的向量，给定一个
   * @p triangulation 和 @p mapping.
   * ，返回的向量的大小等于活动单元的数量。对于非本地拥有的单元，该向量包含零。单元的长宽比定义为Jacobian的最大奇异值与最小奇异值之比，取通过
   * @p quadrature.  指定的正交规则的所有正交点的最大值。
   * 例如，对于尺寸为  $a$  和  $b$
   * 的2D矩形元素的特殊情况（  $a \geq b$
   * ），这个函数返回通常的长宽比定义  $a/b$
   * 。上述使用奇异值的定义是对任意变形元素的一种概括。这个函数旨在用于
   * $d=2,3$ 空间维度，但它也可以用于 $d=1$ 返回值为1。
   * @note
   * 颠倒的元素不会抛出一个异常。相反，在倒置元素的情况下，一个inf的值被写入向量。
   * @note
   * 确保使用足够的正交点，以便在变形元素的情况下精确计算纵横比。
   * @note
   * 在并行计算中，返回值将有n_active_cells的长度，但长宽比只计算本地拥有的、分别放置在索引
   * CellAccessor::active_cell_index(),
   * 的单元。所有其他的值都被设置为0。
   *
   */
  template <int dim>
  Vector<double>
  compute_aspect_ratio_of_cells(const Mapping<dim> &      mapping,
                                const Triangulation<dim> &triangulation,
                                const Quadrature<dim> &   quadrature);

  /**
   * 通过取所有单元的最大值来计算最大长宽比。
   * @note
   * 当与支持MPI的Triangulation并行运行时，这是一个集体调用，其返回值是所有处理器的最大值。
   *
   */
  template <int dim>
  double
  compute_maximum_aspect_ratio(const Mapping<dim> &      mapping,
                               const Triangulation<dim> &triangulation,
                               const Quadrature<dim> &   quadrature);

  /**
   * 计算包含整个三角形的最小的盒子。
   * 如果输入的三角形是一个 `parallel::distributed::Triangulation`,
   * ，那么每个处理器将计算一个包围所有本地拥有的、幽灵的和人工的单元的包围盒。在一个没有弯曲边界的域的情况下，这些边界盒在处理器之间都是一致的，因为人工和鬼魂单元所占据的区域的联合等于其他处理器拥有的单元所占据的区域的联合。
   * 然而，如果域有弯曲的边界，情况就不再是这样了。
   * 返回的边界盒可能适合于当前的处理器，但与其他处理器上计算的边界盒不同。
   *
   */
  template <int dim, int spacedim>
  BoundingBox<spacedim>
  compute_bounding_box(const Triangulation<dim, spacedim> &triangulation);

  /**
   * 返回几何对象 @p object 上最接近给定点 @p trial_point.
   * 的点。例如，如果 @p object
   * 是一条一维线或边，那么返回的点将是连接顶点的几何线上的一个点，因为与该对象相关的流形看到它（即，如果该几何线生活在高维空间，它可能是弯曲的）。如果迭代器指向高维空间中的四边形，那么返回的点位于相关流形所看到的四边形顶点的凸壳内。
   * @note
   * 这种投影通常不是很好解决，因为对象上可能有多个点使距离最小化。这个函数中使用的算法是稳健的（而且输出保证在给定的
   * @p object)
   * 上，但如果物体具有高曲率，可能只提供几个正确的数字。如果你的流形支持它，那么专门的函数
   * Manifold::project_to_manifold() 可能表现得更好。
   *
   */
  template <typename Iterator>
  Point<Iterator::AccessorType::space_dimension>
  project_to_object(
    const Iterator &                                      object,
    const Point<Iterator::AccessorType::space_dimension> &trial_point);

  /**
   * 返回定义三角网格的粗略网格的数组。这个函数是
   * Triangulation::create_triangulation().
   * 的逆函数，返回值是一个包含顶点向量、单元向量和SubCellData结构的元组。后者包含关于面和线的额外信息。
   * 这个函数在需要解构三角图或以某种方式操作顶点编号的情况下非常有用：一个例子是
   * GridGenerator::merge_triangulations().  。
   *
   */
  template <int dim, int spacedim>
  std::
    tuple<std::vector<Point<spacedim>>, std::vector<CellData<dim>>, SubCellData>
    get_coarse_mesh_description(const Triangulation<dim, spacedim> &tria);

   /*@}*/ 
  /**
   * @name  支持创建网格的函数
   *
   */
   /*@{*/ 

  /**
   * 删除不被任何单元格引用的顶点。这个函数被所有
   * <tt>GridIn::read_*</tt>
   * 函数调用，以消除输入文件中列出的但不被输入文件中的单元所使用的顶点。虽然这些顶点从一开始就不应该出现在输入文件中，但有时会出现，最常见的是当一些单元格被手工删除而不想更新顶点列表时，因为它们可能很冗长。
   * 这个函数被所有 <tt>GridIn::read_*</tt>
   * 函数调用，因为三角形类要求它们只用已使用的顶点来调用。
   * 因为顶点是由该类逐字复制的，所以我们必须事先消除未使用的顶点。
   * 在维度为1的情况下没有实现。
   *
   */
  template <int dim, int spacedim>
  void
  delete_unused_vertices(std::vector<Point<spacedim>> &vertices,
                         std::vector<CellData<dim>> &  cells,
                         SubCellData &                 subcelldata);

  /**
   * 移除重复的顶点，例如由于输入了结构化的网格而导致的。如果这些顶点没有被移除，这些顶点所包围的面会成为边界的一部分，即使它们在网格的内部。
   * 这个函数被一些 <tt>GridIn::read_*</tt>
   * 函数所调用。只有索引在 @p considered_vertices
   * 中的顶点才会被测试是否相等。这加快了算法的速度，对于最坏的超立方体几何形状
   * $O(N^{3/2})$ 的二维和 $O(N^{5/3})$
   * 的三维，算法是相当慢的。
   * 然而，如果你希望考虑所有顶点，只需传递一个空矢量。在这种情况下，该函数会将所有顶点填入
   * @p considered_vertices 。
   * 如果两个顶点在每个坐标方向上的差异小于 @p tol.
   * ，则认为它们是相等的。
   *
   */
  template <int dim, int spacedim>
  void
  delete_duplicated_vertices(std::vector<Point<spacedim>> &all_vertices,
                             std::vector<CellData<dim>> &  cells,
                             SubCellData &                 subcelldata,
                             std::vector<unsigned int> &   considered_vertices,
                             const double                  tol = 1e-12);

  /**
   * 由网格生成器生成的网格可能有一个单元的方向，这个方向是deal.II要求的方向的倒数。
   * 在2D和3D中，这个函数检查所有单元是否有负的或正的量度/体积。在前一种情况下，所有的单元格都是倒置的。在1d中，它没有任何作用。
   * 当所有单元格中只有一个子集的体积为负时，单元格的反转也可能起作用。然而，由负向和正向单元混合组成的网格很可能被打破。因此，如果单元格的方向不一致，就会抛出一个异常。
   * @note  这个函数应该在  GridTools::consistently_order_cells().
   * @param  all_vertices 网格的顶点前调用。    @param  cells
   * 描述网格拓扑结构的CellData对象的阵列。
   *
   */
  template <int dim, int spacedim>
  void
  invert_all_negative_measure_cells(
    const std::vector<Point<spacedim>> &all_vertices,
    std::vector<CellData<dim>> &        cells);

  /**
   * 给出一个描述网格的CellData对象的向量，重新排列其顶点，使所有线条的方向一致。    关于方向的期望和这个函数的讨论可以在 @ref reordering "重新排序模块 "
   * 中找到。      @param  cells
   * 描述网格拓扑结构的CellData对象的数组。
   *
   */
  template <int dim>
  void
  consistently_order_cells(std::vector<CellData<dim>> &cells);

   /*@}*/ 
  /**
   * @name  旋转、拉伸和其他变换网格的方法
   *
   */
   /*@{*/ 

  /**
   * 通过对所有顶点应用作为第一个参数提供的函数对象来变换给定的三角网格的顶点。
   * 作为参数给出的变换被用于变换每个顶点。
   * 其各自的类型必须提供类似函数的语法，即谓词要么是一个具有<tt>operator()</tt>类型的对象，要么是一个指向函数的指针。在这两种情况下，参数和返回值都必须是<tt>Point
   * @<spacedim@></tt>.  类型。
   * @note
   * 与该函数一起使用的有意义的变换应该有一个具有正行列式的雅各布系数。例如，旋转、剪切、拉伸或缩放都满足这一点（尽管没有要求使用的变换实际上是线性的，因为所有这些例子都是）。另一方面，反射或反转有一个雅各布式的负行列式。目前的函数没有办法断定雅各布的正行列式，但是如果你碰巧使用了这样的变换，其结果将是一个单元格体积为负的三角结构。
   * @note  如果你使用的是 parallel::distributed::Triangulation
   * ，即使你的 "全局
   * "网格没有悬空节点，你的局部三角结构也会有悬空节点。如果你调用当前的函数，这将导致悬空节点在幽灵单元中的错误定位问题。所有本地拥有的单元的顶点将是正确的，但一些幽灵单元的顶点可能不是。这意味着像KellyErrorEstimator这样的计算可能会给出错误的答案。
   * @note
   * 这个函数一般来说与附加在三角形上的流形不兼容。例如，为了在网格转换后细化网格（使用流形），你必须确保原始流形对转换后的几何体仍然有效。这在一般情况下是不成立的，在这种情况下，有必要清除流形，并为转换后的几何体附加一个新的流形。
   * 如果你想根据附加到三角形上的原始流形描述进行细化，你应该先进行细化，随后停用所有流形，最后再调用transform()函数。其结果是一个具有正确转换顶点的三角形，但在其他方面是直边元素。建议采用以下程序
   * @code
   * ...
   * triangulation.refine_global(n_refinements);
   * triangulation.reset_all_manifolds();
   * Transformation<dim> transformation;
   * GridTools::transform(transformation, triangulation);
   * ...
   * @endcode
   * 这个函数在  step-38  的 "扩展的可能性
   * "部分中使用。它也在  step-49  和  step-53  中使用。
   *
   */
  template <int dim, typename Transformation, int spacedim>
  void
  transform(const Transformation &        transformation,
            Triangulation<dim, spacedim> &triangulation);

  /**
   * 通过给定的移位矢量对三角形的每个顶点进行移位。这个函数使用了上面的transform()函数，所以那里所说的对三角形的要求也适用于这个函数。
   *
   */
  template <int dim, int spacedim>
  void
  shift(const Tensor<1, spacedim> &   shift_vector,
        Triangulation<dim, spacedim> &triangulation);


  /**
   * 将给定的二维三角结构的所有顶点围绕坐标系的原点逆时针旋转给定的角度（用弧度而不是度数给定）。这个函数使用了上面的transform()函数，所以那里所说的对三角形的要求也适用于这个函数。
   * @note 这个函数只支持dim=2的情况。
   *
   */
  template <int dim>
  void
  rotate(const double angle, Triangulation<dim> &triangulation);

  /**
   * 将给定 @p triangulation
   * 的所有顶点以逆时针方向围绕给定索引的轴旋转。否则就像上面的函数一样。
   * @param[in]  angle 以弧度为单位，将三角函数旋转的角度。
   * @param[in]  axis
   * 围绕坐标轴的索引，保持该坐标的固定（0=x轴，1=y轴，2=z轴）。
   * @param[in,out]  triangulation 要旋转的三角测量对象。
   * @note 对dim=1、2和3实施。
   *
   */
  template <int dim>
  void
  rotate(const double           angle,
         const unsigned int     axis,
         Triangulation<dim, 3> &triangulation);

  /**
   * 将给定的三角剖分平滑地转换到一个不同的域，通常，三角剖分边界的每个顶点都被映射到
   * @p new_points 地图中的相应点。    方向 $u_d(\mathbf x)$
   * 的未知位移场 $d$ 由受规定约束的最小化问题\f[ \min\,
   * \int \frac{1}{2} c(\mathbf x) \mathbf \nabla u_d(\mathbf x) \cdot \mathbf
   * \nabla u_d(\mathbf x) \,\rm d x
   * \f]获得。最小化器是通过解决位移场的dim分量的拉普拉斯方程得到的，该位移场将当前域映射成由
   * @p new_points
   * 描述的域。使用线性有限元，每个方向上有四个高斯正交点。因此，
   * @p new_points 中指定的顶点位置与 @p tria
   * 中的当前值之间的差值代表该位移场在域边界的规定值，或者更准确地说，在
   * @p new_points
   * 提供的所有位置（可能是在边界的一部分，甚至在域的内部）。然后，该函数在每个无约束的顶点上评估这个位移场，并使用它将映射的顶点放在位移场定位的地方。因为拉普拉斯方程的解是平滑的，这保证了从旧域到新域的平滑映射。
   * @param[in]  new_points
   * 要放置现有顶点的子集的位置。通常，这将是一个从边界上所有节点的顶点指数到其新位置的映射，从而完全指定了映射域的几何形状。然而，如果有必要，它也可以包括内部的点，而且它不需要包括所有的边界顶点（尽管你会失去对映射域的确切形状的控制）。
   * @param[in,out]  tria
   * Triangulation对象。这个对象被就地改变，即之前的顶点位置被覆盖。
   * @param[in]  coefficient 拉普拉斯问题的可选系数。
   * 较大的值使单元格不容易变形（有效地增加其刚度）。该系数是在三角形的旧的、未变形的配置的坐标系中作为输入进行评估的，也就是说，在应用变换之前。
   * 如果提供这个函数，只有在所有的系数都是正数的情况下才能期望得到合理的结果。
   * @param[in]  solve_for_absolute_positions 如果设置为
   * <code>true</code>
   * ，则最小化问题是针对最终顶点位置而非其位移制定的。这两个公式对于同质问题是等价的（默认值为
   * @p coefficient),
   * ），但在其他情况下，它们会导致非常不同的网格运动。由于在大多数情况下，我们会在位移公式中使用一个非恒定系数，这个参数的默认值是
   * <code>false</code>  。
   * @note 这个功能目前还没有在1d情况下实现。
   *
   */
  template <int dim>
  void
  laplace_transform(const std::map<unsigned int, Point<dim>> &new_points,
                    Triangulation<dim> &                      tria,
                    const Function<dim, double> *coefficient = nullptr,
                    const bool solve_for_absolute_positions  = false);

  /**
   * 返回一个  std::map  具有位于边界的所有顶点的面
   * @param[in]  Triangulation 对象。
   *
   */
  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  get_all_vertices_at_boundary(const Triangulation<dim, spacedim> &tria);

  /**
   * 按给定的系数缩放整个三角图。为了保持三角形的方向，该因子必须是正的。
   * 这个函数使用了上面的transform()函数，所以对三角形的要求也适用于这个函数。
   *
   */
  template <int dim, int spacedim>
  void
  scale(const double                  scaling_factor,
        Triangulation<dim, spacedim> &triangulation);

  /**
   * 通过随机移动网格中的所有顶点来扭曲给定的三角结构。
   * 每个顶点的移动方向是随机的，而移位矢量的长度为 @p
   * factor 乘以与此顶点相邻的活动边的最小长度。  注意，
   * @p factor 显然应该远低于<tt>0.5</tt>。    如果 @p keep_boundary
   * 被设置为 @p true
   * （这是默认的），那么边界顶点不会被移动。      @p seed
   * 用于随机引擎的初始化。它的默认值用以前版本的deal.II中的相同状态来初始化引擎。
   *
   */
  template <int dim, int spacedim>
  void
  distort_random(
    const double                  factor,
    Triangulation<dim, spacedim> &triangulation,
    const bool                    keep_boundary = true,
    const unsigned int            seed = boost::random::mt19937::default_seed);

  /**
   从网格中移除悬挂的节点。如果 @p isotropic 参数设置为 @p false （默认值），该函数会检测到有悬空节点的单元格，并在去除悬空节点的方向上细化邻域。  如果 @p isotropic 参数设置为 @p true, ，则在每个方向上进行邻居细化。  为了去除所有悬空节点，这个过程必须重复进行：这可能需要大量的迭代。  为了避免这种情况，我们提供了一个最大的迭代数（ @p max_iterations) ）。    考虑以下网格。    @image html remove_hanging_nodes-hanging.png   @p isotropic  ==  @p false  将返回。    @image html remove_hanging_nodes-aniso.png   @p isotropic  ==  @p true  会返回。    @image html remove_hanging_nodes-isotro.png   @param[in,out]  tria Triangulation来细化。      @param[in]  isotropic 如果为真，则在每个方向上细化单元，否则（默认值）则在删除悬挂节点的方向上细化单元。      @param[in]  max_iterations 在每一步中，只有最接近悬挂节点的单元被精炼。该代码可能需要大量的迭代来移除所有悬空节点。  @p max_iterations  是允许的最大迭代数。如果  @p max_iterations  ==  numbers::invalid_unsigned_int  这个函数继续精炼，直到没有悬空节点。
   * @note  在并行代码的情况下，该函数应与
   * GridGenerator::flatten_triangulation. 相结合。
   *
   */
  template <int dim, int spacedim>
  void
  remove_hanging_nodes(Triangulation<dim, spacedim> &tria,
                       const bool                    isotropic      = false,
                       const unsigned int            max_iterations = 100);

  /**
   各向异性地细化一个网格，使得到的网格由尺寸间最大比率小于 @p max_ratio. 的单元组成，这个过程需要一个可能不会终止的算法。因此，可以通过 @p max_iterations 参数设置一个最大的迭代次数。    从这样的一个单元格开始。    @image html remove_anisotropy-coarse.png  这个函数将返回。    @image html remove_anisotropy-refined.png   @param[in,out]  tria 要细化的三角图。      @param[in]  max_ratio 每个单元的尺寸之间的比率允许的最大值。      @param[in]  max_iterations 允许的最大迭代次数。
   * @note  如果是并行代码，该函数应与
   * GridGenerator::flatten_triangulation 和 GridTools::remove_hanging_nodes.
   * 相结合。
   *
   */
  template <int dim, int spacedim>
  void
  remove_anisotropy(Triangulation<dim, spacedim> &tria,
                    const double                  max_ratio      = 1.6180339887,
                    const unsigned int            max_iterations = 5);

  /**
   * 分析网格的边界单元，如果发现有一个单元位于角的位置（边界上有昏暗的相邻面），并且其二维角度分数超过
   * @p limit_angle_fraction,
   * ，则全局细化一次，并且用角不再违反给定角度分数的子单元替换该单元的子单元。
   * 如果不存在边界上有两个相邻面的边界单元，那么三角形就不会被触动。如果我们确实有边界上有dim相邻面的单元，那么将根据参数
   * @p limit_angle_fraction.
   * 检查dim-dimensional实体角和dim*pi/2之间的分数，如果它更高，网格将被细化一次，并且违规单元的子单元将被替换成一些尊重极限的单元。
   * 在这个过程之后，三角形被压平，所有的Manifold对象都被恢复到原来的三角形中。
   * 以下的网格就是一个例子，它是将一个SphericalManifold附加到使用
   * GridGenerator::hyper_cube: 生成的网格上而得到的。
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * tria.refine_global(4);
   * @endcode
   * <p ALIGN="center">
   @image html regularize_mesh_01.png
   * </p>
   * 四个原本是正方形的角的单元在计算过程中会给你带来一些麻烦，因为从参考单元到这些单元的变换的贾可宾会归零，影响有限元估计的误差常数。
   * 那些单元格的角非常接近180度，也就是说，角分值非常接近1。
   * 同样的代码，加入对regularize_corner_cells的调用。
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * GridTools::regularize_corner_cells(tria);
   * tria.refine_global(2);
   * @endcode
   * 生成的网格在Mapping的jacobian方面有更好的表现。      <p
   * ALIGN="center">
   @image html regularize_mesh_02.png
   * </p>  这个网格与 GridGenerator::hyper_ball.
   * 得到的网格非常相似。然而，使用
   * GridTools::regularize_corner_cells
   * 可以自由选择何时应用正则化，也就是说，原则上可以先细化几次，然后再调用regularize_corner_cells函数。
   * @code
   * const SphericalManifold<dim> m0;
   * Triangulation<dim> tria;
   * GridGenerator::hyper_cube(tria,-1,1);
   * tria.set_all_manifold_ids_on_boundary(0);
   * tria.set_manifold(0, m0);
   * tria.refine_global(2);
   * GridTools::regularize_corner_cells(tria);
   * tria.refine_global(1);
   * @endcode
   * 这就产生了下面的网格。      <p ALIGN="center">
   @image html regularize_mesh_03.png
   * </p>  这个函数目前只在dim = 2的情况下实现，如果在dim =
   * 3的情况下调用，会产生一个异常。      @param[in,out]  tria
   * 正则化的三角图。      @param[in]  limit_angle_fraction
   * 网格中的角元素允许的最大角度或实体角度的比率。
   *
   */
  template <int dim, int spacedim>
  void
  regularize_corner_cells(Triangulation<dim, spacedim> &tria,
                          const double limit_angle_fraction = .75);

   /*@}*/ 
  /**
   * @name  寻找三角形的单元和顶点
   *
   */
   /*@{*/ 

  /**
   * 给定一个三角形的  @p cache  和一个  @p points,
   * 的列表，在  @p points,  的每个元素上调用
   * find_active_cell_around_point() 并返回  @p cells,  参考位置  @p
   * qpoints,  和一个从局部到全局索引的映射  @p maps  到  @p
   * points  。      @param[in]  缓存 三角形的 GridTools::Cache  .
   * @param[in]  点 点的向量。    @param[in]  cell_hint (可选)
   * 可能包含  @p points.   @return  的第一个点的单元格迭代器
   * 包含以下信息。
   *
   *
   *
   *
   *
   * -  @p cells  : 所有包含至少一个 @p points  的单元格的向量。
   *
   *
   *
   *
   *
   *
   * -  @p qpoints  : 一个点的向量。  @p qpoints[i]  包含落在单元格内的所有点的参考位置  @p cells[i]  。
   *
   *
   *
   *
   *
   * -  @p indices  : 一个整数向量，包含  @p qpoints  中的局部编号，和  @p points  中的全局索引之间的映射。    如果 @p points[a] 和 @p points[b] 是落在 @p cells[c], 中的唯一两个点，那么 @p qpoints[c][0] 和 @p qpoints[c][1] 是 @p points[a] 和 @p points[b] 在 @p cells[c], 和 @p indices[c][0] 的参考位置。 ]=a， @p indices[c][1] =b。函数 Mapping::transform_unit_to_real(qpoints[c][0]) 返回 @p points[a]. 。算法建立了一个 @p points 的rtree，对它们进行空间排序，然后尝试调用find_active_cell_around_point（）。
   * @note  这个函数没有在二维一的情况下实现（<tt>spacedim !=
   * dim</tt>）。
   * @note  如果一个点在网格内没有找到，或者位于一个
   * parallel::TriangulationBase,
   * 的人工单元内，这个点会被默默地忽略掉。如果你想推断哪些点的搜索失败了，请使用函数compute_point_locations_try_all()，该函数也会返回一个索引的向量，表示搜索失败的点。
   * @note  这个函数的实际返回类型，即上面提到的 @p
   * return_type, 的类型是
   * @code
   * std::tuple<
   * std::vector<
   *   typename Triangulation<dim, spacedim>::active_cell_iterator>,
   * std::vector<std::vector<Point<dim>>>,
   * std::vector<std::vector<unsigned int>>>
   * @endcode
   * 在线文档中对该类型进行了缩写，以提高本页面的可读性。
   * @note  这个函数通过利用
   * GridTools::Cache::get_cell_bounding_boxes_rtree(),
   * 来优化搜索，该函数要么返回一个缓存的rtree，要么建立并存储一个。如果该函数只在少数几个点上调用一次，建立一个rtree可能会妨碍性能。
   *
   */
  template <int dim, int spacedim>
#  ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>>
#  else
  return_type
#  endif
  compute_point_locations(
    const Cache<dim, spacedim> &        cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint =
        typename Triangulation<dim, spacedim>::active_cell_iterator());

  /**
   * 这个函数与 GridTools::compute_point_locations(),
   * 类似，但是compute_point_locations()默默地忽略了所有find_active_cell_around_point()失败的点，这个函数也返回一个包含find_active_cell_around_point()失败的点的索引的向量。
   * @return  一个包含四个元素的元组；前三个元素在
   * GridTools::compute_point_locations().  中有记录。  @p return_type
   * 的最后一个元素包含了既没有在网格内发现也没有位于人工单元内的点的索引。
   * @p return_type  等于以下元组类型。
   * @code
   * std::tuple<
   *   std::vector<
   *      typename Triangulation<dim,spacedim>::active_cell_iterator>,
   *   std::vector<std::vector<Point<dim>>>,
   *   std::vector<std::vector<unsigned int>>,
   *   std::vector<unsigned int>
   * >
   * @endcode
   *
   * @note  这个函数在二维一的情况下没有实现（<tt>spacedim !=
   * dim</tt>）。
   * @note  这个函数通过使用
   * GridTools::Cache::get_cell_bounding_boxes_rtree(),
   * 来优化搜索，该函数要么返回一个缓存的rtree，要么建立并存储一个。如果该函数只在少数几个点上调用一次，建立一个rtree可能会妨碍性能。
   * 更详细的文档见 GridTools::compute_point_locations(). 。
   *
   */
  template <int dim, int spacedim>
#  ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<unsigned int>>
#  else
  return_type
#  endif
  compute_point_locations_try_all(
    const Cache<dim, spacedim> &        cache,
    const std::vector<Point<spacedim>> &points,
    const typename Triangulation<dim, spacedim>::active_cell_iterator
      &cell_hint =
        typename Triangulation<dim, spacedim>::active_cell_iterator());

  /**
   * 给出一个 @p cache 和每个进程的 @p local_points
   * 列表，找到位于网格本地拥有部分的点并计算它们的正交规则。
   * 分布式计算点位置是一个类似于
   * GridTools::compute_point_locations 的函数，但对
   * parallel::TriangulationBase
   * 对象起作用，而且，与它的序列版本不同，也适用于分布式三角测量（见
   * parallel::distributed::Triangulation).   @param[in] 缓存一个
   * GridTools::Cache 对象 @param[in]
   * local_points是当前进程拥有的点的阵列。
   * 每个进程可以有一个不同的点阵列，这个阵列可以是空的，不包含在三角形的本地拥有的部分内
   * @param[in]  global_bboxes
   * 一个边界盒的向量；它描述了每个进程的网格的本地拥有部分。global_boxes[rk]中包含了描述网格的哪一部分被等级为rk的进程局部拥有的界线盒。局部描述可以从
   * GridTools::compute_mesh_predicate_bounding_box;
   * 中获得，然后全局描述可以通过
   * GridTools::exchange_local_bounding_boxes 或 Utilities::MPI::all_gather
   * 获得  @param[in]  容差
   * 以单元格坐标计算的容差。根据问题的不同，可能需要调整公差，以便能够确定一个单元。浮点运算意味着，一般来说，一个点不会完全位于一个顶点、边缘或面。在任何一种情况下，都无法预测这个函数会返回哪个与顶点或边/面相邻的单元。
   * 因此，调用这个函数的算法需要考虑到，返回的单元格将只包含近似的点。
   * @return  一个包含正交信息的元组 输出元组的元素是。
   *
   *
   *
   *
   *
   * - cells : 所有包含至少一个点的单元格的向量。
   *
   *
   *
   *
   *
   *
   * - qpoints : 一个点的向量；包含在 @p qpoints[i] 中的所有点的参考位置，这些点位于单元格 @p cells[i] 中。
   *
   *
   *
   *
   *
   *
   * - maps : 一个整数向量，包含qpoints中的编号（元组的前一个元素）与拥有这些点的过程的局部点向量之间的映射。
   *
   *
   *
   *
   *
   *
   * - points : 一个点的向量。  @p points[i][j]  是实空间中对应于  @p qpoints[i][j]  的点。注意 @p points 是位于网格本地所有部分的点；因此这些可以是 @p local_points 的副本或从其他进程收到的点，即其他进程的本地_点
   *
   *
   *
   *
   *
   *
   * - owners : 一个向量的向量； @p owners[i][j]  包含拥有point[i][j]的进程的等级（元组的前一个元素）。    该函数使用三角形的mpi通信器：由于这个原因，如果三角形不是从 parallel::TriangulationBase 派生的，它会抛出一个断言错误。    在一个串行执行中，元组的前三个元素与  GridTools::compute_point_locations  中的相同。    注意：这个函数是一个集体操作。
   * @note  这个函数的实际返回类型，即上面提到的  @p
   * return_type,  的类型是
   * @code
   * std::tuple<
   * std::vector<
   *   typename Triangulation<dim, spacedim>::active_cell_iterator>,
   * std::vector<std::vector<Point<dim>>>,
   * std::vector<std::vector<unsigned int>>,
   * std::vector<std::vector<Point<spacedim>>>,
   * std::vector<std::vector<unsigned int>>>
   * @endcode
   * 在线文档中对该类型进行了缩写，以提高本页面的可读性。
   *
   */
  template <int dim, int spacedim>
#  ifndef DOXYGEN
  std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>,
    std::vector<std::vector<Point<dim>>>,
    std::vector<std::vector<unsigned int>>,
    std::vector<std::vector<Point<spacedim>>>,
    std::vector<std::vector<unsigned int>>>
#  else
  return_type
#  endif
  distributed_compute_point_locations(
    const GridTools::Cache<dim, spacedim> &                cache,
    const std::vector<Point<spacedim>> &                   local_points,
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const double                                           tolerance = 1e-10);

  namespace internal
  {
    /**
     * 由 GridTools::internal::distributed_compute_point_locations().
     * 返回的数据结构 它提供信息以执行
     * GridTools::distributed_compute_point_locations() 并在
     * Utilities::MPI::RemotePointEvaluation::reinit().
     * 内设置通信模式。
     * @note  字段的名称是在考虑到
     * Utilities::MPI::RemotePointEvaluation
     * 的情况下选择的。在这里，数量在指定的任意定位点（甚至在MPI宇宙中的远程进程上）被逐个计算，这些值被发送给请求进程，请求进程接收结果，并根据点的情况诉诸结果。
     *
     */
    template <int dim, int spacedim>
    struct DistributedComputePointLocationsInternal
    {
      /**
       * 发送/评估方的每个点的信息。该元组的元素如下。0）单元水平和索引，1）拥有进程的等级，2）拥有进程的本地索引，3）参考位置，4）实际位置，5）发送缓冲区内的包络索引。
       *
       */
      std::vector<std::tuple<std::pair<int, int>,
                             unsigned int,
                             unsigned int,
                             Point<dim>,
                             Point<spacedim>,
                             unsigned int>>
        send_components;

      /**
       * 要发送的等级。
       *
       */
      std::vector<unsigned int> send_ranks;

      /**
       * 发送缓冲区内的范围指针，将被发送至send_ranks指定的等级。发送缓冲区的大小由send_ptrs.back()给出。
       *
       */
      std::vector<unsigned int> send_ptrs;

      /**
       * 每个收到的数据值的信息。该元组的元素如下。0）发送者的等级，1）本地索引，2）枚举索引。
       * @note  向量按照1）、0）、2）进行排序。)
       * @note
       * 每个点都可能有多个数据值与之相关。如果一个点与一个由多个单元共享的几何实体（例如，顶点）重合，就可能是这种情况。
       *
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        recv_components;

      /**
       * 从哪里接收数据的等级。
       *
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * 接收缓冲区内的范围指针，由recv_ranks指定的等级来填充。接收缓冲区的大小由recv_ptrs.back()给出。
       *
       */
      std::vector<unsigned int> recv_ptrs;
    };

    /**
     * 一个填充DistributedComputePointLocationsInternal的函数。
     * 如果输入参数 @p perform_handshake 被设置为false，只有
     * GridTools::internal::distributed_compute_point_locations()
     * 需要的字段被填充。    如果输入参数被设置为
     * "true"，则会设置额外的数据结构，以便能够在
     * Utilities::MPI::RemotePointEvaluation::reinit().
     * 中设置通信模式。
     *
     */
    template <int dim, int spacedim>
    DistributedComputePointLocationsInternal<dim, spacedim>
    distributed_compute_point_locations(
      const GridTools::Cache<dim, spacedim> &                cache,
      const std::vector<Point<spacedim>> &                   points,
      const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
      const double                                           tolerance,
      const bool                                             perform_handshake,
      const bool enforce_unique_mapping = false);

  } // namespace internal

  /**
   * 返回一个地图`顶点索引
   *
   * ->
   * Point<spacedim>`，包含给定的`容器'的使用顶点。返回的地图的键（即上面一对的第一个元素）是三角形中的全局索引，而每一对的值是对应顶点的物理位置。使用的顶点是通过在所有单元中循环得到的，并通过（可选）`mapping`参数查询每个单元的顶点在哪里。
   * 在序列Triangulation对象和 parallel::shared::Triangulation
   * 对象中，返回的地图大小等于 Triangulation::n_used_vertices()
   * （而不是 Triangulation::n_vertices()).  注意，在
   * parallel::distributed::Triangulation
   * 对象中，只返回本地拥有的单元和幽灵单元中的顶点，因为所有其他顶点的真实位置可能不知道（例如，对于使用MappingQEulerian的分布计算）。
   * 如果你使用默认的`mapping'，返回的地图满足以下等价关系。
   * @code
   * const auto used_vertices = extract_used_vertices(tria);
   * auto all_vertices = tria.get_vertices();
   *
   * for(const auto &id_and_v : used_vertices)
   * all_vertices[id_and_v.first] == id_and_v.second; // true
   * @endcode
   * 注意，对于改变顶点位置的映射，如MappingQEulerian，则不满足上述规定。      @ref ConceptMeshType  "MeshType概念"
   * 。    @param  container 要提取顶点的容器。    @param  mapping
   * 用来计算点的位置的映射。
   *
   */
  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim>>
  extract_used_vertices(
    const Triangulation<dim, spacedim> &container,
    const Mapping<dim, spacedim> &      mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 在作为第一个参数传递的顶点映射中，查找并返回离给定点最近的顶点的索引。
   * @param  vertices 索引->顶点的地图，如
   * GridTools::extract_used_vertices().   @param  p 目标点。    @return
   * 最接近目标点`p`的顶点的索引。
   *
   */
  template <int spacedim>
  unsigned int
  find_closest_vertex(const std::map<unsigned int, Point<spacedim>> &vertices,
                      const Point<spacedim> &                        p);

  /**
   * 查找并返回给定网格中最接近给定点的已使用顶点（或标记顶点）的索引。    这个函数使用存储在三角结构中的顶点位置。这通常是足够的，除非你使用一个移动顶点的Mapping（例如，MappingQEulerian）。在这种情况下，你应该用相同的名字和额外的Mapping参数来调用这个函数。      @param  mesh 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型的变量。    @param  p
   * 我们想找到最接近的顶点的点。    @param  marked_vertices
   * 表示 @p mesh
   * 的哪些顶点将在搜索中被视为潜在的最近顶点的一个布尔数组。当收到一个非空的
   * @p marked_vertices, 时，该函数将只在 @p marked_vertices
   * 中搜索最接近的顶点。  这个数组的大小应该等于
   * Triangulation::n_vertices()
   * 对给定网格的三角结构返回的值（而不是
   * Triangulation::n_used_vertices()). 返回的值）  @return
   * 找到的最接近顶点的索引。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex(const MeshType<dim, spacedim> &mesh,
                      const Point<spacedim> &        p,
                      const std::vector<bool> &      marked_vertices = {});

  /**
   * 在给定的网格中找到并返回最接近给定点的已使用顶点（或标记顶点）的索引。使用给定的映射来计算顶点的实际位置。    如果Mapping不修改网格顶点的位置（例如，MappingQEulerian），那么这个函数等同于同名的函数，并且没有`mapping`参数。      @param  mapping 用于计算顶点位置的映射  @param  mesh 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型的变量。    @param  p
   * 我们想找到最近的顶点的点。    @param  marked_vertices 表示
   * @p mesh
   * 的哪些顶点将在搜索中被视为潜在的最近顶点的一个布尔数组。当收到一个非空的
   * @p marked_vertices, 时，该函数将只在 @p marked_vertices
   * 中搜索最接近的顶点。  这个数组的大小应该等于
   * Triangulation::n_vertices()
   * 对给定网格的三角结构返回的值（而不是
   * Triangulation::n_used_vertices()). 返回的值）  @return
   * 找到的最接近顶点的索引。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex(const Mapping<dim, spacedim> & mapping,
                      const MeshType<dim, spacedim> &mesh,
                      const Point<spacedim> &        p,
                      const std::vector<bool> &      marked_vertices = {});


  /**
   * 找到并返回一个迭代器的向量，这些迭代器围绕着给定顶点的索引  @p vertex_index.  对于局部细化网格，顶点本身可能不是返回的所有相邻单元的顶点。然而，它将始终是一个单元的顶点，或者是位于面或边上的一个悬挂节点。      @param  容器 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型的变量。    @param  vertex_index
   * 我们试图找到相邻单元的顶点的索引。    @return
   * 与给定顶点相邻的单元格的一个向量。
   * @note
   * 目前还不完全清楚该函数是否对各向异性的细化网格做出正确的处理。它需要对这种情况进行检查。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  std::vector<typename MeshType<dim, spacedim>::active_cell_iterator>
#  else
  std::vector<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type>
#  endif
  find_cells_adjacent_to_vertex(const MeshType<dim, spacedim> &container,
                                const unsigned int             vertex_index);

  /**
   * 查找围绕给定点的活动非人工单元  @p p.  返回类型是一对活动单元的迭代器以及该点的单元坐标。    这个函数使用的算法是首先寻找最接近给定点的顶点，见 GridTools::find_closest_vertex().  其次，在网格中找到这个顶点的所有相邻单元，见 GridTools::find_cells_adjacent_to_vertex().  最后，对于每个单元，函数测试点是否在里面。这个检查是使用给定的 @p mapping 参数来确定单元的边界是直的还是弯的。    如果一个点位于两个或多个单元的边界上，那么该算法将试图确定细化程度最高的那个单元。    如果请求的点不在本地拥有的单元或幽灵单元中，那么这个函数将返回（无效的）MeshType<dim,  spacedim>::end()  迭代器。这种情况可以类似于各种 `std::find()` 和 `std::lower_bound()` 函数的处理方式。      @param  映射 用于确定给定点是否在给定单元内的映射。    @param  mesh 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型的变量。    @param  p
   * 我们想找到周围单元的点。    @param  marked_vertices
   * 一个`bool'数组，表示顶点数组中的某个条目是否应该被视为（而其他的必须被忽略）可能是离指定点最近的顶点。在指定一个非默认的
   * @p marked_vertices, 时，find_closest_vertex()只会在 @p
   * marked_vertices 中搜索最近的顶点。
   * 这个数组的大小应该等于三角形的n_vertices()（而不是n_used_vertices()）。使用
   * @p marked_vertices
   * 的动机是为了减少顶点的搜索空间，如果人们对感兴趣的点可能接近的顶点集合有先验的了解。
   * @param  容差
   * 以单元格坐标为单位的容差。根据问题的不同，可能有必要调整公差，以便能够识别一个单元。浮点运算意味着，一般来说，一个点不会完全位于一个顶点、边缘或面。在任何一种情况下，都无法预测这个函数会返回哪个与顶点或边/面相邻的单元。
   * 因此，调用这个函数的算法需要考虑到，返回的单元格将只包含点的近似值。
   * @return
   * 一对进入网格的迭代器，指向周围的单元格，以及该点的单元格坐标。由于数字上的舍入，这个局部位置可能位于实际单元格之外。因此，这个函数返回的点应该被投影到单元格上，使用
   * GeometryInfo::project_to_unit_cell().
   * 这不是由算法自动执行的。返回的单元格可以是本地拥有的单元格或幽灵单元格（但不是人造单元格）。即使给定的点是本地拥有的单元格的一个顶点，返回的单元格也可能是一个幽灵单元。
   * 背后的原因是，这是保证所有参与平行三角形计算的处理器都同意哪个单元包含一个点的唯一方法。例如，如果两个处理器聚集在一个顶点，并且用这个顶点调用该函数，那么一个处理器将返回一个本地拥有的单元，另一个则返回一个幽灵单元。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#  else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#  endif
  find_active_cell_around_point(const Mapping<dim, spacedim> & mapping,
                                const MeshType<dim, spacedim> &mesh,
                                const Point<spacedim> &        p,
                                const std::vector<bool> &marked_vertices = {},
                                const double             tolerance = 1.e-10);

  /**
   * 上述函数的一个版本，假定边界是直的，因此只是用MappingQ1作为映射参数调用上述函数。
   * @return 一个进入网格的迭代器，指向周围的单元。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  typename MeshType<dim, spacedim>::active_cell_iterator
#  else
  typename dealii::internal::
    ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type
#  endif
  find_active_cell_around_point(const MeshType<dim, spacedim> &mesh,
                                const Point<spacedim> &        p,
                                const std::vector<bool> &marked_vertices = {},
                                const double             tolerance = 1.e-10);

  /**
   * 另一个版本，我们在一个给定的单元上使用该映射，该映射对应于该单元的活动有限元索引。
   * 这显然只对hp-problems有用，因为所有其他DoF处理程序的活动有限元索引总是零。
   *
   */
  template <int dim, int spacedim>
  std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           mesh,
    const Point<spacedim> &                     p,
    const double                                tolerance = 1.e-10);

  /**
   * 在一个点周围寻找一个活跃的非人工单元在计算成本上可能是非常昂贵的。这个函数旨在通过使用空间树来加速几何体的搜索，提供上述函数的快速版本。
   * @param  cache 包含三角形空间树信息的对象，见
   * GridTools::Cache.   @param  p 我们要为其寻找周围的单元。
   * @param  cell_hint
   * 给出几何搜索的提示，如果有关于该点可能位于哪个单元的先验知识，这将是有益的。一个典型的用例是，这个搜索必须针对一个相互靠近的点的阵列，并且前一个点的相邻单元是阵列中下一个点的良好提示。
   * @param  marked_vertices 见上文。    @param  tolerance 见上文。
   * 下面的代码示例显示了如何使用这个函数。
   * @code
   * GridTools::Cache<dim, dim> cache(triangulation, mapping);
   * auto cell_hint = typename Triangulation<dim, dim>::active_cell_iterator();
   * std::vector<bool> marked_vertices = {};
   * double tolerance = 1.e-10;
   *
   * std::vector<Point<dim>> points; // a vector of many points
   * ...
   *
   * for(auto p : points)
   * {
   * auto cell_and_ref_point = GridTools::find_active_cell_around_point(
   *   cache, p, cell_hint, marked_vertices, tolerance);
   *
   * if (cell_and_ref_point.first != triangulation.end())
   *   {
   *    // use current cell as hint for the next point
   *    cell_hint = cell_and_ref_point.first;
   *    // do something with cell_and_ref_point
   *    ...
   * }
   * else
   *  {
   *     // The function did not find a locally owned or ghost cell in which
   *     // the point is located. We ought to handle this somehow here.
   *  }
   * ...
   * }
   * @endcode
   *
   *
   */
  template <int dim, int spacedim>
  std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator,
            Point<dim>>
  find_active_cell_around_point(
    const Cache<dim, spacedim> &cache,
    const Point<spacedim> &     p,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &
                             cell_hint = typename Triangulation<dim, spacedim>::active_cell_iterator(),
    const std::vector<bool> &marked_vertices = {},
    const double             tolerance       = 1.e-10);

  /**
   * 前一个函数的一个版本，利用顶点和单元格之间已经存在的映射（使用函数
   * GridTools::vertex_to_cell_map()),
   * 构建一个顶点_到单元格_中心的映射（通过
   * GridTools::vertex_to_cell_centers_directions()),
   * 获得，也可以选择从三角结构的使用顶点构建的RTree。
   * @note  所有这些结构都可以从一个 GridTools::Cache
   * 对象中查询到。但是请注意，在这种情况下，MeshType必须是Triangulation，所以在这种情况下，直接调用上面的函数，参数为`cache'可能更合适。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim>>
#  else
  std::pair<typename dealii::internal::
              ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
            Point<dim>>
#  endif
  find_active_cell_around_point(
    const Mapping<dim, spacedim> & mapping,
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim> &        p,
    const std::vector<
      std::set<typename MeshType<dim, spacedim>::active_cell_iterator>>
      &                                                  vertex_to_cell_map,
    const std::vector<std::vector<Tensor<1, spacedim>>> &vertex_to_cell_centers,
    const typename MeshType<dim, spacedim>::active_cell_iterator &cell_hint =
      typename MeshType<dim, spacedim>::active_cell_iterator(),
    const std::vector<bool> &                              marked_vertices = {},
    const RTree<std::pair<Point<spacedim>, unsigned int>> &used_vertices_rtree =
      RTree<std::pair<Point<spacedim>, unsigned int>>{},
    const double tolerance = 1.e-10,
    const RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      *relevant_cell_bounding_boxes_rtree = nullptr);

  /**
   * 与上面的函数相比，这个函数以单位坐标的方式识别一个给定的容忍度`tolerance`的点周围所有活跃的非人工单元。给定一个参考坐标为参数
   * @p first_cell,
   * 的第一个单元，例如通过上面的一个函数得到的，所有相应的具有单位坐标点的邻近单元也被确定。
   * 这个函数对不连续函数空间很有用，例如，对于给定的点`p`位于一个顶点、边缘或面的情况，几个单元可能持有独立的解的值，在用户代码中以某种方式组合。
   * 这个函数的使用方法如下
   * @code
   * auto first_pair = GridTools::find_active_cell_around_point(...);
   * auto all_cells  = GridTools::find_all_active_cells_around_point(
   * 			   mapping, mesh, p, tolerance, first_pair);
   * @endcode
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#  else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#  endif
  find_all_active_cells_around_point(
    const Mapping<dim, spacedim> & mapping,
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim> &        p,
    const double                   tolerance,
    const std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                    Point<dim>> &  first_cell);

  /**
   * 前一个函数的变体，在内部调用其中一个函数find_active_cell_around_point()来获得第一个单元，随后通过调用上面的函数find_all_active_cells_around_point()来增加所有其他活跃的非人工单元。
   *
   */
  template <int dim, template <int, int> class MeshType, int spacedim>
#  ifndef _MSC_VER
  std::vector<std::pair<typename MeshType<dim, spacedim>::active_cell_iterator,
                        Point<dim>>>
#  else
  std::vector<std::pair<
    typename dealii::internal::
      ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim>>::type,
    Point<dim>>>
#  endif
  find_all_active_cells_around_point(
    const Mapping<dim, spacedim> & mapping,
    const MeshType<dim, spacedim> &mesh,
    const Point<spacedim> &        p,
    const double                   tolerance       = 1e-10,
    const std::vector<bool> &      marked_vertices = {});

  /**
   * 返回给定单元格的所有活跃的后代的列表。例如，如果当前单元格曾经被精炼过，但是它的子代没有任何进一步的精炼，那么返回的列表将包含它的所有子代。    如果当前单元格已经被激活，那么返回的列表是空的（因为该单元格没有可能被激活的子代）。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。    @param  cell
   * 指向Mesh的一个单元的迭代器。    @return
   * 给定单元格的活动子孙列表。
   * @note
   * 因为在C++中MeshType模板参数不能从函数调用中推导出来，所以你必须在函数名称后面指定它，例如
   * @code
   * GridTools::get_active_child_cells<DoFHandler<dim> > (cell)
   * @endcode
   *
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  get_active_child_cells(const typename MeshType::cell_iterator &cell);

  /**
   * 提取给定单元格 @p cell 周围的活动单元，并在向量 @p active_neighbors. 中返回这些邻居，这些邻居具体是指单元格的<i>face</i>邻居，如果该邻居被进一步细化，则是其与该面交界的活动子女。另一方面，返回的邻居不包括那些位于，例如，与一个顶点对角线相对但本身不是面的邻居的单元。在3D中，它也不包括与当前单元格的一条边相邻，但不是面的邻居的单元格）。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。    @param[in]  cell
   * 一个指向Mesh中某一单元的迭代器。    @param[out]
   * active_neighbors 指向给定单元格的活跃子孙的列表。
   * @note
   * 因为在C++中，MeshType模板参数不能从函数调用中推导出来，你必须在函数名称后指定它，例如
   * @code
   * GridTools::get_active_neighbors<DoFHandler<dim>>(cell, active_neighbors)
   * @endcode
   *
   *
   */
  template <class MeshType>
  void
  get_active_neighbors(
    const typename MeshType::active_cell_iterator &       cell,
    std::vector<typename MeshType::active_cell_iterator> &active_neighbors);

  /**
   * 提取并返回 @p mesh
   * 中的子域（活动单元的集合）周围的活动单元层（即那些与子域共享一组顶点但不属于子域的单元）。在这里，"子域
   * "恰好包括 @p 谓词返回 @p true. 的所有单元。
   * 一个自定义谓词的例子是检查一个给定的材料id
   * @code
   * template <int dim>
   * bool
   * pred_mat_id(const typename Triangulation<dim>::active_cell_iterator & cell)
   * {
   * return cell->material_id() ==  1;
   * }
   * @endcode
   * 然后我们可以通过以下调用提取这个材料周围的细胞层。
   * @code
   * GridTools::compute_active_cell_halo_layer(tria, pred_mat_id<dim>);
   * @endcode
   * 经常有用的谓词可以在命名空间IteratorFilters中找到。例如，可以提取所有具有给定材料ID的细胞周围的细胞层。
   * @code
   * GridTools::compute_active_cell_halo_layer(
   * tria, IteratorFilters::MaterialIdEqualTo(1, true));
   * @endcode
   * 或者在具有hp-capabilities的DoFHandler的一组活动FE指数的所有单元周围提取一层单元。
   * @code
   * GridTools::compute_active_cell_halo_layer(
   * hp_dof_handler, IteratorFilters::ActiveFEIndexEqualTo({1,2}, true));
   * @endcode
   * 注意，在最后两个例子中，我们确保谓词只对本地拥有的单元返回真。这意味着光环层将不包含任何人工单元。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。    @param[in]  mesh
   * 一个网格（即Triangulation或DoFHandler类型的对象）。
   * @param[in]  谓词
   * 一个函数（或带有operator()的类型对象），定义要提取晕层的子域。它是一个接收活动单元并返回一个布尔值的函数。
   * @return
   * 一个与所指子域至少有一个共同顶点的活动单元的列表。
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_active_cell_halo_layer(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &predicate);


  /**
   * 提取并返回 @p mesh
   * 指定层次上的子域（单元格集合）周围的单元格层（即该层次上与子域共享一组共同顶点但不属于子域的那些单元格）。在这里，"子域
   * "恰好由 @p predicate 返回 @p true. 的所有单元组成。
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::cell_iterator>
  compute_cell_halo_layer_on_level(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::cell_iterator &)>
      &                predicate,
    const unsigned int level);


  /**
   * 提取并返回幽灵单元，这些单元是所有本地拥有的单元周围的活动单元层。这与 parallel::shared::Triangulation 最为相关，它将返回一个处理器上所有幽灵单元的子集，但对于 parallel::distributed::Triangulation 来说，这将返回所有的幽灵单元。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。    @param[in]  mesh
   * 一个网格（即Triangulation或DoFHandler类型的对象）。
   * @return  一个幽灵单元的列表。
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_ghost_cell_halo_layer(const MeshType &mesh);

  /**
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_active_cell_layer_within_distance(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &          predicate,
    const double layer_thickness);

  /**
   * 提取并返回一组幽灵单元，这些单元在所有本地拥有的单元周围的 @p layer_thickness 内。  这与 parallel::shared::Triangulation 最相关，它将返回一个进程中所有幽灵单元的子集，但对于 parallel::distributed::Triangulation 这将返回所有的幽灵单元。  对于 parallel::shared::Triangulation 类来说，所有不属于当前处理器的单元格都可以被认为是幽灵单元格；特别是，它们不仅仅是在本地拥有的单元格周围形成一个单层。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。    @param  mesh
   * 一个网格（即Triangulation或DoFHandler类型的对象）。
   * @param  layer_thickness
   * 指定函数从本地拥有的单元中搜索活动单元的几何距离。
   * @return  在给定的几何距离 @p
   * layer_thickness内的鬼魂单元子集与当前进程的本地拥有的单元。
   * 参见compute_ghost_cell_halo_layer() 和
   * compute_active_cell_layer_within_distance()。
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_ghost_cell_layer_within_distance(const MeshType &mesh,
                                           const double    layer_thickness);

  /**
   * 计算并返回一个通过左下角和右上角的一对点定义的包围盒，该包围盒围绕着
   * @p mesh. 的一个子域。这里，"子域 "恰好由 @p predicate
   * 返回 @p true. 的所有活动单元组成。 关于 @p predicate
   * 如何工作的描述，见compute_active_cell_halo_layer（）。
   * @note  这个函数是在BoundingBox类被发明之前写的。
   * 因此，它返回一对点，而不是人们期望的BoundingBox对象。然而，BoundingBox有一个从点对转换的构造函数，所以这个函数的结果仍然可以被分配给一个BoundingBox对象。
   *
   */
  template <class MeshType>
  std::pair<Point<MeshType::space_dimension>, Point<MeshType::space_dimension>>
  compute_bounding_box(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &predicate);

  /**
   * 计算一个边界框的集合，使所有给定谓词为真的活动单元都完全被包围在至少一个边界框中。请注意，这个包围只保证包含所有这些活动单元，但它不一定是精确的，也就是说，它可以包括比它们的联合体更大的区域。
   * 对于一个给定的细化级别中包含 @p predicate
   * 为真的活动单元的每个单元，该函数创建一个 @p predicate
   * 为真的子单元的边界盒。    这导致了对 @p predicate
   * 为真的所有活动单元的覆盖；参数 @p allow_merge 和 @p
   * max_boxes
   * 用于减少计算成本的单元数量，覆盖更大的n维体积。
   * 控制该算法的参数是。
   *
   *
   *
   *
   *
   *
   * -  @p predicate  : 要包围的单元格的属性，例如  IteratorFilters::LocallyOwnedCell  。   该谓词仅在活动单元格上进行测试。
   *
   *
   *
   *
   *
   * -  @p refinement_level  : 它定义了创建初始边界盒的级别。细化应该被设置为粗略的细化级别。如果 @p refinement_level 高于三角形的层数，将为每个活动单元创建一个比 @p refinement_level; 更粗的包围盒，将产生一个异常。
   *
   *
   *
   *
   *
   *
   * -  @p allow_merge  : 这个标志允许盒子合并，默认为假。该算法的成本为O(N^2)，其中N是由细化级别创建的边界盒的数量；由于这个原因，如果该标志被设置为真，请确保明智地选择一个足够粗的  @p refinement_level.
   *
   *
   *
   *
   *
   *
   * -  @p max_boxes  : 要计算的边界盒的最大数量。如果创建了更多的盒子，那么小的盒子就会与相邻的盒子合并。默认情况下，在合并了可以表示为一个的盒子后，不再合并更多的盒子。详见 BoundingBox::get_neighbor_type （）函数。   注意只有相邻的单元格会被合并（见边界盒类中的 @p get_neighbor_type 函数）：如果边界盒的目标数量max_boxes不能通过合并相邻的单元格来达到，则会抛出一个异常。 下面的图片描述了一个算法的例子， @p refinement_level  = 2,  @p allow_merge  = true and  @p max_boxes  = 1。带有属性谓词的单元格是红色的，包围盒的区域略带橙色。    @image html bounding_box_predicate.png
   *
   *
   *
   *
   *
   * - 1.在黑色中我们可以看到当前级别的单元格。
   *
   *
   *
   *
   *
   * - 2. 对于每个包含红色区域的单元格，都会创建一个边界框：默认情况下，这些框会被返回。
   *
   *
   *
   *
   * 因为 @p allow_merge  =
   * true，所以在不改变封面的情况下减少了包围盒的数量。
   * 如果 @p max_boxes
   * 被保留为默认值或大于1，这两个盒子将被返回。
   *
   *
   *
   *
   *
   *
   * 因为 @p max_boxes
   * =1，最小的边界盒被合并到较大的边界盒。
   * 注意，明智地选择参数是很重要的。例如， @p allow_merge
   * =false和 @p refinement_level
   * =1会返回非常相同的边界框，但计算成本只有一小部分。
   * 这个函数没有考虑到单元格的曲率，因此它不适合处理弯曲的几何图形：映射被假定为线性。
   *
   */
  template <class MeshType>
  std::vector<BoundingBox<MeshType::space_dimension>>
  compute_mesh_predicate_bounding_box(
    const MeshType &mesh,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &                predicate,
    const unsigned int refinement_level = 0,
    const bool         allow_merge      = false,
    const unsigned int max_boxes        = numbers::invalid_unsigned_int);

  /**
   * 给定一个点阵列，使用使用
   * GridTools::compute_mesh_predicate_bounding_box
   * 获得的全局边界盒描述来猜测，对于每个点，哪个进程可能拥有它。
   * @param[in]  global_bboxes
   * 描述每个进程拥有属性的网格部分的边界盒的矢量。
   * @param[in]  points 要测试的点的阵列。      @return
   * 一个包含以下信息的元组。
   *
   *
   *
   *
   *
   *
   * - 一个以进程的等级为标志的向量。对于每个等级，它包含一个它可能拥有的点的索引的向量。
   *
   *
   *
   *
   *
   * - 从 <code>unsigned int</code> 中的点的索引 @p points 到所有者的等级的地图。
   *
   *
   *
   *
   *
   *
   * - 从 @p points 中的点的索引 <code>unsigned int</code> 到被猜测的所有者的行列的地图。
   * @note  这个函数的实际返回类型，即上面提到的 @p
   * return_type, 的类型是
   * @code
   * std::tuple<std::vector<std::vector<unsigned int>>,
   *          std::map< unsigned int, unsigned int>,
   *          std::map< unsigned int, std::vector<unsigned int>>>
   * @endcode
   * 在线文档中对该类型进行了缩写，以提高本页面的可读性。
   *
   */
  template <int spacedim>
#  ifndef DOXYGEN
  std::tuple<std::vector<std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#  else
  return_type
#  endif
  guess_point_owner(
    const std::vector<std::vector<BoundingBox<spacedim>>> &global_bboxes,
    const std::vector<Point<spacedim>> &                   points);


  /**
   * 给定一个覆盖的rtree（见 GridTools::Cache::get_covering_rtree()),
   * ）和一个点的数组，找到一个进程的超集，这个超集可以单独拥有包含这些点的单元。
   * 进一步的细节见 GridTools::guess_point_owner;
   * 这里只报告不同的输入/输出类型。      @param[in]
   * covering_rtree
   * RTRee，它使我们能够识别并行计算中哪些进程可能拥有围绕给定点的单元。
   * @param[in]  points 要考虑的点的一个向量。      @return
   * 一个包含以下信息的元组。
   *
   *
   *
   *
   *
   *
   * - 一个以处理器等级为索引的地图。对于每个等级，它包含一个它可能拥有的点的索引向量。
   *
   *
   *
   *
   *
   *
   * - 从 <code>unsigned int</code> 中的点的索引 @p points 到所有者的等级的地图；这些是找到单一可能所有者的点。
   *
   *
   *
   *
   *
   *
   * - 从 <code>unsigned int</code> 中的点的索引 @p points 到猜测的所有者行列的地图；这些是发现有多个可能的所有者的点。
   * @note  这个函数的实际返回类型，即上面提到的 @p
   * return_type, 的类型是
   * @code
   * std::tuple<std::map<unsigned int, std::vector<unsigned int>>,
   *          std::map<unsigned int, unsigned int>,
   *          std::map<unsigned int, std::vector<unsigned int>>>
   * @endcode
   * 在线文档中对该类型进行了缩写，以提高本页面的可读性。
   *
   */
  template <int spacedim>
#  ifndef DOXYGEN
  std::tuple<std::map<unsigned int, std::vector<unsigned int>>,
             std::map<unsigned int, unsigned int>,
             std::map<unsigned int, std::vector<unsigned int>>>
#  else
  return_type
#  endif
  guess_point_owner(
    const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &covering_rtree,
    const std::vector<Point<spacedim>> &                         points);


  /**
   * 返回所有顶点的相邻单元。如果一个顶点也是一个悬空的节点，也会返回相关的粗略单元。顶点是按顶点索引排序的。这是由函数
   * <code>cell-@>vertex_index()</code>
   * 返回的数字。注意，只使用由
   * Triangulation<dim,spacedim>::get_used_vertices()
   * 返回的数组中标记的索引。
   *
   */
  template <int dim, int spacedim>
  std::vector<
    std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
  vertex_to_cell_map(const Triangulation<dim, spacedim> &triangulation);

  /**
   * 为 GridTools::vertex_to_cell_map()
   * 输出的每个顶点-单元组合返回一个归一化张量的向量（期望作为此函数的输入参数）。每个张量代表一个从顶点到各自单元中心的几何向量。
   * 如果输入向量的大小不等于三角形的顶点数量，将抛出一个断言。
   * result[v][c]是顶点索引v的单位张量，表示第c个单元的中心相对于顶点v的方向。
   *
   */
  template <int dim, int spacedim>
  std::vector<std::vector<Tensor<1, spacedim>>>
  vertex_to_cell_centers_directions(
    const Triangulation<dim, spacedim> &mesh,
    const std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      &vertex_to_cells);


  /**
   * 返回最接近给定位置的单元格 @p cell 的局部顶点索引  @p
   * position.  顶点的位置从（可选） @p mapping
   * 参数中提取，以保证在底层映射修改顶点位置时返回正确答案。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  find_closest_vertex_of_cell(
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const Point<spacedim> &                                            position,
    const Mapping<dim, spacedim> &                                     mapping =
      (ReferenceCells::get_hypercube<dim>()
         .template get_default_linear_mapping<dim, spacedim>()));

  /**
   * 为与本地拥有的活动单元相关的每个顶点和悬挂节点计算一个全局唯一的索引。作为本地拥有的单元格的悬挂节点的幽灵单元格的顶点有一个全局索引。
   * 然而，不<i>touch</i>一个活动单元的其他顶点在这个处理器上没有全局索引。
   * 地图的键是顶点的本地索引，值是全局索引。这些索引在细化或粗化后需要重新计算，并且可能是不同的。
   *
   */
  template <int dim, int spacedim>
  std::map<unsigned int, types::global_vertex_index>
  compute_local_to_global_vertex_index_map(
    const parallel::distributed::Triangulation<dim, spacedim> &triangulation);

  /**
   * 返回一个 @p cell.
   * 的每个坐标方向上的外延之间的比率中的最高值
   * 此外，返回相对于最高伸长率的尺寸。      @param[in]
   * cell一个指向单元格的迭代器。      @return  一个
   * std::pair<unsigned  int, double>，这样 @p first
   * 值是最高伸长率的尺寸， @p second 值是 @p cell.
   * 尺寸中的比率。
   *
   */
  template <int dim, int spacedim>
  std::pair<unsigned int, double>
  get_longest_direction(
    typename Triangulation<dim, spacedim>::active_cell_iterator cell);

   /*@}*/ 
  /**
   * @name  三角形的分区和子域
   *
   */
   /*@{*/ 

  /**
   * 产生一个稀疏模式，其中非零条目表示两个单元格通过一个共同的面连接。稀疏模式的对角线条目也被设定。
   * 行和列指的是使用单元格迭代器按自然顺序遍历的单元格。
   *
   */
  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern &            connectivity);

  /**
   * 产生一个稀疏度模式，其中非零条目表示两个单元格通过一个共同的顶点连接。稀疏模式的对角线条目也被设定。
   * 行和列指的是使用单元格迭代器按自然顺序遍历的单元格。
   *
   */
  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells(
    const Triangulation<dim, spacedim> &triangulation,
    DynamicSparsityPattern &            connectivity);

  /**
   * 为一个给定的水平网格产生一个稀疏模式，其中非零条目表示两个单元通过一个共同的顶点连接。稀疏模式的对角线条目也被设置。
   * 行和列指的是使用单元格迭代器按自然顺序遍历的单元格。
   *
   */
  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells_on_level(
    const Triangulation<dim, spacedim> &triangulation,
    const unsigned int                  level,
    DynamicSparsityPattern &            connectivity);

  /**
   * 使用图形分割器来分割构成整个域的活动单元。调用此函数后，所有活动单元的子域id的值将在0和
   * @p n_partitions-1. 之间，你可以通过使用<tt>cell-
   * @>subdomain_id()</tt>.
   * 使用第三个参数来选择METIS或ZOLTAN提供的分区算法。METIS是默认的分区器。
   * 如果deal.II没有与ZOLTAN或METIS一起安装，当选择相应的分区方法时，这个函数将产生一个错误，除非
   * @p n_partitions 是一个。
   * 即，你可以写一个程序，使其在单处理器单分区的情况下运行，而不安装软件包，只有在需要多分区时才需要安装。
   * @note 如果 @p cell_weight 信号已被附加到 @p triangulation,
   * ，那么这将被使用并传递给分区器。
   *
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          Triangulation<dim, spacedim> &   triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * 这个函数执行的操作与上面的函数相同，只是它考虑到了一组特定的
   * @p cell_weights,
   * ，它允许分区器平衡图形，同时考虑到每个单元所花费的计算努力。
   * @note  如果 @p cell_weights
   * 向量为空，则不考虑加权。如果不是，那么这个向量的大小必须等于三角形中有效单元的数量。
   *
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          Triangulation<dim, spacedim> &   triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * 这个函数与前一个函数的作用相同，即使用分区算法将一个三角形划分为由
   * <code>cell-@>subdomain_id()</code> 标志确定的若干子域。
   * 与前一个函数不同的是第二个参数，一个代表单元格之间连接模式的稀疏模式。
   * 虽然上面的函数通过考虑哪些单元彼此相邻而直接从三角图中建立，但这个函数可以采用更精细的连接图。稀疏模式的大小需要是
   * $N\times N$ ，其中 $N$
   * 是三角形中活动单元的数量。如果稀疏模式在位置
   * $(i,j)$ 处包含一个条目，那么这意味着单元格 $i$ 和 $j$
   * （按照主动单元格迭代器遍历的顺序）将被视为连接；然后分区算法将尝试以这样的方式划分域：（i）子域的大小大致相等，以及（ii）最小数量的连接被破坏。
   * 这个函数主要适用于单元格之间存在仅在三角形中不存在的连接的情况（否则前面的函数将是更简单的用法）。这种连接可能包括域的边界的某些部分通过对称边界条件或积分进行耦合（例如，域中裂缝两边的摩擦接触），或者如果使用的数值方案不仅连接紧邻的单元，而且连接更大的邻近单元（例如，在求解积分方程时）。
   * 此外，在默认的稀疏模式不完全足够的情况下，这个函数可能是有用的。这种情况可能会发生，因为默认情况下只是考虑面的邻居，而不是由边或顶点连接的相邻单元。虽然在使用连续有限元时，后者夫妇在邻接图中通常仍是紧密相连的，分区算法在这种情况下通常不会切断重要的连接。然而，如果网格中存在许多单元的顶点（分别比2D和3D中常见的4个或6个多得多）聚集在一起，那么就会有大量的单元跨顶点连接，但在仅使用面邻关系构建的连接图中却有几度的距离。在这样的情况下，分区算法有时可能会做出错误的决定，你可能想建立自己的连接图。
   * @note 如果 @p cell_weight 信号已被附加到 @p triangulation,
   * ，那么这将被使用并传递给分区器。
   *
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int            n_partitions,
                          const SparsityPattern &       cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * 这个函数执行的操作与上面的函数相同，只是它考虑到了一组特定的
   * @p cell_weights,
   * ，它允许分区器平衡图形，同时考虑到每个单元上所花费的计算努力。
   * @note  如果 @p cell_weights
   * 向量为空，则不考虑加权。如果不是，那么这个向量的大小必须等于三角形中有效单元的数量。
   *
   */
  template <int dim, int spacedim>
  void
  partition_triangulation(const unsigned int               n_partitions,
                          const std::vector<unsigned int> &cell_weights,
                          const SparsityPattern &       cell_connection_graph,
                          Triangulation<dim, spacedim> &triangulation,
                          const SparsityTools::Partitioner partitioner =
                            SparsityTools::Partitioner::metis);

  /**
   * 如果标志 @p group_siblings 被设置为
   * "true"（此函数的默认行为），则使用与p4est库中相同的分区方案生成构成整个域的活动单元的分区。
   * 调用此函数后，所有活动单元的子域id的值将在0和 @p
   * n_partitions-1. 之间。你可以通过使用<tt>cell-
   * @>subdomain_id()</tt>. 来访问一个单元的子域id。
   * @note  如果标志 @p group_siblings
   * 被设置为false，一个单元的子域可能会被放在不同的处理器上，即使它们都处于活动状态，这是p4est的一个假设。通过放宽这一点，我们可以创建拥有单个单元的分区（也适用于精炼网格）。
   *
   */
  template <int dim, int spacedim>
  void
  partition_triangulation_zorder(const unsigned int            n_partitions,
                                 Triangulation<dim, spacedim> &triangulation,
                                 const bool group_siblings = true);

  /**
   * 通过使用 "最年轻的孩子
   * "规则分配级别子域id来划分多网格层次结构的单元，也就是说，层次结构中的每个单元都由在森林中拥有其最左边孩子的处理器拥有，活跃的单元具有相同的子域id和级别子域id。你可以通过使用<tt>cell-
   * @>level_subdomain_id()</tt>.
   * 注意：这个函数假定活动单元已经被分区。
   *
   */
  template <int dim, int spacedim>
  void
  partition_multigrid_levels(Triangulation<dim, spacedim> &triangulation);

  /**
   * 该函数允许询问由CellId对象识别的单元格的所属子域，该对象在当前进程中不一定存在。
   * @note 这个函数还没有为 parallel::fullydistributed::Triangulation.
   * 实现。
   *
   */
  template <int dim, int spacedim>
  std::vector<types::subdomain_id>
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            const std::vector<CellId> &         cell_ids);

  /**
   * 对于每个活动单元，在输出数组中返回它属于哪个子域（由<tt>cell->subdomain_id()</tt>函数给出）。在调用此函数时，输出数组应该已经有了合适的大小。
   * 这个函数返回每个单元格与一个子域的关联。如果你要寻找每个
   * @em DoF与一个子域的关联，请使用
   * <tt>DoFTools::get_subdomain_association</tt> 函数。
   *
   */
  template <int dim, int spacedim>
  void
  get_subdomain_association(const Triangulation<dim, spacedim> &triangulation,
                            std::vector<types::subdomain_id> &  subdomain);

  /**
   * 计算有多少个单元与给定的 @p subdomain 索引唯一相关。
   * 如果没有具有给定 @p
   * 子域索引的单元格，该函数可能返回0。这种情况可能发生，例如，如果你试图将一个粗略的网格划分为更多的分区（每个处理器一个），而不是网格中的单元。
   * 这个函数返回与一个子域相关的单元数。
   * 如果你正在寻找 @em DoF与这个子域的关联，请使用
   * <tt>DoFTools::count_dofs_with_subdomain_association</tt> 函数。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  count_cells_with_subdomain_association(
    const Triangulation<dim, spacedim> &triangulation,
    const types::subdomain_id           subdomain);

  /**
   * 对于一个三角形，返回一个掩码，代表哪些顶点被当前进程 "拥有"
   * ，就像我们谈论本地拥有的单元或自由度一样（见 @ref
   * GlossLocallyOwnedCell 和 @ref GlossLocallyOwnedDof  ）。
   * 为了这个函数的目的，我们对本地拥有的顶点定义如下：一个顶点是由与该顶点相邻的所有单元的所有者中具有最小的子域id（相当于该处理器的MPI等级）的那个处理器所拥有。换句话说，位于三角形分区内部的顶点由这个分区的所有者拥有；对于位于两个或多个分区之间边界的顶点，所有者是所有相邻子域中拥有最小子域id的处理器。
   * 对于顺序三角计算（相对于，例如
   * parallel::distributed::Triangulation),
   * 每个用户顶点当然是由当前处理器拥有的，即函数返回
   * Triangulation::get_used_vertices().
   * 对于并行三角计算，返回的掩码是
   * Triangulation::get_used_vertices() 返回的一个子集。      @param
   * triangulation
   * 该函数评估哪些顶点是本地拥有的三角结构。    @return
   * 顶点的子集，如上所述。返回的数组长度等于Triangulation.n_vertices()，因此，可能大于
   * Triangulation::n_used_vertices(). 。
   *
   */
  template <int dim, int spacedim>
  std::vector<bool>
  get_locally_owned_vertices(const Triangulation<dim, spacedim> &triangulation);

   /*@}*/ 
  /**
   * @name  比较不同的网格
   *
   */
   /*@{*/ 

  /**
   * 给出两个基于相同粗略网格的网格（即Triangulation或DoFHandler类型的对象），这个函数找出一组在这两个网格之间匹配的单元，其中最多只有一个网格在这个单元上更精细。换句话说，它找到了两个网格共同的最小单元，并且这些单元一起完全覆盖了该领域。    这个函数很有用，例如，在随时间变化的或非线性的应用中，我们必须对一个网格（例如，前一个时间步骤或非线性迭代的网格）上定义的解决方案与另一个网格（下一个时间步骤，下一个非线性迭代）的形状函数进行积分。例如，如果新的网格更细，那么就必须在粗的网格（Mesh_1）上获得解决方案，并将其内插到Mesh_2的相应单元中。反之，如果新的网格更粗，我们就必须用细网格形状函数的线性组合来表达粗网格的形状函数。无论哪种情况，我们都需要循环计算两个三角形中共同的最细的单元。这个函数返回一个与两个网格中的单元匹配的迭代器对的列表，可以用来达到这个目的。    请注意，这些迭代器的列表不一定是有序的，也不一定与作为参数的一个或两个网格中的单元格被遍历的顺序相一致。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。
   * @note  这个函数只能与 parallel::distributed::Triangulation
   * 一起使用，当两个网格都使用相同的三角法时，因为在分布式三角法中，并非所有的单元都存储在本地，所以产生的列表可能不会覆盖整个域。
   *
   */
  template <typename MeshType>
  std::list<std::pair<typename MeshType::cell_iterator,
                      typename MeshType::cell_iterator>>
  get_finest_common_cells(const MeshType &mesh_1, const MeshType &mesh_2);

  /**
   * 如果两个三角形是基于相同的粗略网格，则返回true。
   * 这是通过检查它们在最粗层次上是否有相同数量的单元来确定的，然后再检查它们是否有相同的顶点。
   * 这两个网格可能有不同的细化历史，超出了粗略的网格。
   *
   */
  template <int dim, int spacedim>
  bool
  have_same_coarse_mesh(const Triangulation<dim, spacedim> &mesh_1,
                        const Triangulation<dim, spacedim> &mesh_2);

  /**
   * 与上面的函数相同，但对DoFHandler类型的参数工作。  提供这个函数是为了允许对所有类型的代表三角形的容器或建立在三角形上的类调用have_same_coarse_mesh。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。
   *
   */
  template <typename MeshType>
  bool
  have_same_coarse_mesh(const MeshType &mesh_1, const MeshType &mesh_2);

   /*@}*/ 
  /**
   * @name  处理扭曲的单元格
   *
   */
   /*@{*/ 

  /**
   * 给出一个三角网格和一个单元格的列表，这些单元格的子节点由于网格细化而变得扭曲，尝试通过移动中心节点来修复这些单元格。    该函数返回一个子节点变形的单元格列表，这些单元格由于某种原因不能被修复。因此，返回的列表是输入参数的一个子集。    关于扭曲的单元格的概念的定义，见 @ref GlossDistorted "词汇表条目"
   * 。  传递给当前函数的第一个参数通常是
   * Triangulation::execute_coarsening_and_refinement 函数抛出的异常。
   *
   */
  template <int dim, int spacedim>
  typename Triangulation<dim, spacedim>::DistortedCellList
  fix_up_distorted_child_cells(
    const typename Triangulation<dim, spacedim>::DistortedCellList
      &                           distorted_cells,
    Triangulation<dim, spacedim> &triangulation);



   /*@}*/ 
  /**
   * @name  提取和创建单元格斑块 这些函数提取和创建围绕单个单元格的单元格斑块，并从中创建三角图。
   *
   */
   /*@{*/ 


  /**
   * 该函数返回给定活动单元的所有活动邻居单元的列表。 这里，邻居被定义为与给定的单元至少有一部分共同的面，但不是边缘（在3D）或顶点邻居（在2D和3D）。    返回列表中的第一个元素是作为参数提供的单元格。  其余的是邻居。该函数在给定的单元格的所有面上循环，并检查该面是否不在域的边界上。然后，如果邻居单元没有任何子单元（也就是说，它与当前单元处于相同的细化水平，或者更粗），那么这个邻居单元将被添加到单元列表中。否则，如果邻接单元是细化的，因此有孩子，那么这个函数就会在当前面的所有子面中循环，将这些子面后面的邻接单元添加到要返回的列表中。      @tparam  MeshType 一个满足 @ref ConceptMeshType  "MeshType概念 "
   * 要求的类型。  在C++中，编译器不能从函数调用中确定
   * <code>MeshType</code>
   * 。你需要把它作为一个明确的模板参数跟在函数名后面指定。
   * @param[in]  cell 指向网格中某一单元的迭代器。    @return
   * 构成给定单元周围补丁的活动单元的列表
   * @note
   * 补丁通常用于定义误差估计器，需要解决网格中每个单元周围的补丁上的局部问题。这也需要操作与补丁的单元相关的自由度。为此，在命名空间DoFTools中有更多的函数在处理补丁。
   * @note
   * 在并行分布式计算的背景下，只有在本地拥有的单元上调用这个函数才有意义。这是因为本地拥有的单元的邻居要么是本地拥有的单元，要么是幽灵单元。对于这两种情况，我们知道这些单元实际上是完整的、平行的三角形的真实单元。我们还可以查询这些单元的自由度。
   *
   */
  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  get_patch_around_cell(const typename MeshType::active_cell_iterator &cell);


  /**
   * 这个函数接收一个活动单元的向量（以下命名为 @p
   * patch_cells）作为输入参数，并返回一个它们的父单元的向量，其细化程度为最粗的公共水平。换句话说，找到那个生活在同一细化水平的细胞集合，使输入向量中的所有细胞都是该集合中的细胞的子代，或者本身就在该集合中。
   * @tparam  容器 在C++中，编译器不能从函数调用中确定
   * <code>Container</code>
   * 的类型。你需要把它作为一个明确的模板参数在函数名后面指定。这个类型必须满足网状容器的要求（见
   * @ref ConceptMeshType  ）。      @param[in]  patch_cells
   * 一个活动单元的向量，本函数为其找到最粗的公共层的父单元。这个单元格向量通常是调用函数的结果
   * GridTools::get_patch_around_cell().   @return
   * 具有输入单元格的最粗共同细化水平的单元格列表。
   *
   */
  template <class Container>
  std::vector<typename Container::cell_iterator>
  get_cells_at_coarsest_common_level(
    const std::vector<typename Container::active_cell_iterator> &patch_cells);

  /**
   *
   */
  template <class Container>
  void
  build_triangulation_from_patch(
    const std::vector<typename Container::active_cell_iterator> &patch,
    Triangulation<Container::dimension, Container::space_dimension>
      &local_triangulation,
    std::map<
      typename Triangulation<Container::dimension,
                             Container::space_dimension>::active_cell_iterator,
      typename Container::active_cell_iterator> &patch_to_global_tria_map);

  /**
   * 这个函数通过DoFHandler定义的自由度运行，并为每个自由度构建一个active_cell_iterators的向量，代表该自由度下相关基元的支持单元。这个函数最初是为实现局部投影而设计的，例如Clement插值，结合其他局部修补函数，如
   * GridTools::build_triangulation_from_patch.
   * DoFHandler的建立在Triangulation之上或
   * parallel:distributed::Triangulation ，都得到支持和适当处理。
   * 其结果是代表与自由度相关的基元支持的单元补丁。
   * 例如，使用FE_Q有限元，我们得到了接触自由度的标准单元补丁，然后添加其他单元来处理可能的悬挂节点约束。
   * 使用FE_DGQ有限元，自由度在逻辑上被认为是单元的
   * "内部"，所以补丁将只由自由度所在的单个单元组成。
   * @param[in]  dof_handler DoFHandler可以建立在三角或
   * parallel::distributed::Triangulation
   * 有限元上，其自由度在逻辑上与顶点、直线、四边形或六边形相关。
   * @return 从局部相关单元上的自由度的global_dof_index到包含
   * DoFHandler::active_cell_iterators
   * 在该自由度的基函数支持中的单元的向量的映射。
   *
   */
  template <int dim, int spacedim>
  std::map<
    types::global_dof_index,
    std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>>
  get_dof_to_support_patch_map(DoFHandler<dim, spacedim> &dof_handler);


   /*@}*/ 

  /**
   * @name  处理周期性域的问题
   *
   */
   /*@{*/ 

  /**
   * 提供所有必要信息的数据类型，以创建周期性约束和相对于两个
   * "周期性 "单元面的周期性p4est森林。
   *
   */
  template <typename CellIterator>
  struct PeriodicFacePair
  {
    /**
     * 与两个'周期性'面相关的单元格。
     *
     */
    CellIterator cell[2];

    /**
     * 两个 "周期性 "面的局部面指数（相对于指定单元）。
     *
     */
    unsigned int face_idx[2];

    /**
     * 在orthogonal_equality()和 DoFTools::make_periodicity_constraints()
     * 中描述的第一个面相对于第二个面的相对方向（并存储为比特集）。
     *
     */
    std::bitset<3> orientation;

    /**
     * 一个 @p dim   $\times$   @p dim 旋转矩阵，描述了在约束到第二个面的DoF之前，应该如何修改第一个面的矢量值DoF。        旋转矩阵在 DoFTools::make_periodicity_constraints() 中使用，对有限元空间的参数 @p first_vector_components 中列出的所有矢量值块进行旋转。更多细节见 DoFTools::make_periodicity_constraints() 和词汇表 @ref GlossPeriodicConstraints "关于周期性条件的词汇表条目"
     * 。
     *
     */
    FullMatrix<double> matrix;
  };


  /**
   * 对面的正交平等测试。      @p face1 和 @p face2
   * 被认为是相等的，如果其顶点之间可以通过正交平等关系实现一对一的匹配。
   * 这里，两个顶点<tt>v_1</tt>和<tt>v_2</tt>被认为是相等的，如果
   * $M\cdot v_1 + offset
   *
   * - v_2$ 与单位方向的单位向量 @p direction. 平行，如果参数
   * @p matrix 是对spacedim x spacedim矩阵的引用， $M$ 被设置为 @p
   * matrix, ，否则 $M$ 为身份矩阵。    如果匹配成功， @p
   * face1 相对于 @p face2 的_相对方向被返回到比特集 @p
   * orientation, 中，其中
   * @code
   * orientation[0]
   *
   * -> face_orientation
   * orientation[1]
   *
   * -> face_flip
   * orientation[2]
   *
   * -> face_rotation
   * @endcode
   * 在2D中，<tt>face_orientation</tt>总是<tt>true</tt>，<tt>face_rotation</tt>总是<tt>false</tt>，而face_flip具有<tt>line_flip</tt>的含义。更确切地说，在3D中。
   * <tt>face_orientation</tt>: <tt>真</tt>如果 @p face1 和 @p face2
   * 有相同的方向。否则， @p face1 的顶点指数与 @p face2
   * 的顶点指数以下列方式匹配。
   * @code
   * face1:           face2:
   *
   * 1
   *
   * - 3            2
   *
   * - 3
   * |   |    <-->    |   |
   * 0
   *
   * - 2            0
   *
   * - 1
   * @endcode
   * <tt>face_flip</tt>: <tt>真</tt>如果匹配的顶点旋转了180度。
   * @code
   * face1:           face2:
   *
   * 1
   *
   * - 0            2
   *
   * - 3
   * |   |    <-->    |   |
   * 3
   *
   * - 2            0
   *
   * - 1
   * @endcode
   * <tt>face_rotation</tt>:
   * <tt>真</tt>如果匹配的顶点逆时针旋转90度。
   * @code
   * face1:           face2:
   *
   * 0
   *
   * - 2            2
   *
   * - 3
   * |   |    <-->    |   |
   * 1
   *
   * - 3            0
   *
   * - 1
   * @endcode
   * 以及任何的组合... 关于该主题的更多信息可以在 @ref GlossFaceOrientation "词汇表 "
   * 文章中找到。
   *
   */
  template <typename FaceIterator>
  bool orthogonal_equality(
    std::bitset<3> &                                              orientation,
    const FaceIterator &                                          face1,
    const FaceIterator &                                          face2,
    const int                                                     direction,
    const Tensor<1, FaceIterator::AccessorType::space_dimension> &offset =
      Tensor<1, FaceIterator::AccessorType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * 与上述函数相同，但不返回实际方向
   *
   */
  template <typename FaceIterator>
  bool
  orthogonal_equality(
    const FaceIterator &                                          face1,
    const FaceIterator &                                          face2,
    const int                                                     direction,
    const Tensor<1, FaceIterator::AccessorType::space_dimension> &offset =
      Tensor<1, FaceIterator::AccessorType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * 这个函数将在给定的 @p mesh
   * （一个三角形或DoFHandler）的最粗的网格层次上收集周期性的面对，并将它们添加到矢量
   * @p matched_pairs 中，而不改变原来的内容。    定义 "第一
   * "边界为所有边界面，其边界ID为 @p  b_id1，"第二
   * "边界由所有属于 @p  b_id2的面组成。
   * 这个函数试图在orthogonal_equality()的帮助下将所有属于第一个边界的面与属于第二个边界的面进行匹配。
   * 在PeriodicFacePair中返回的比特集编码了第一个面相对于第二个面的_相对_方向，更多细节请参见orthogonal_equality()的文档。
   * @p direction
   * 指的是周期性被强制执行的空间方向。当匹配周期性面时，这个向量分量被忽略。
   * @p offset 是一个与面相切的矢量，当试图将 "第一
   * "边界的顶点与 "第二
   * "边界的相应顶点匹配时，该矢量将被添加到 "第一
   * "边界的顶点位置。这可以用来实现诸如 $u(0,y)=u(1,y+1)$
   * 等条件。    可以选择指定一个 $dim\times dim$ 旋转 @p matrix
   * ，描述在约束到第二个面的DoF之前，应该如何修改第一个面的矢量值DoF。
   * @p matrix 在两个地方使用。首先， @p matrix
   * 将被提供给orthogonal_equality()并用于匹配面。如果
   * $\text{matrix}\cdot v_1 + \text{offset}
   *
   * @note  创建的 std::vector 可以在
   * DoFTools::make_periodicity_constraints() 和
   * parallel::distributed::Triangulation::add_periodicity()
   * 中使用，以代数方式强制实现周期性。
   * @note  因为元素将被添加到 @p matched_pairs
   * 中（而现有的条目将被保留），所以可以用不同的边界ID多次调用这个函数来生成一个具有所有周期对的向量。
   * @note
   * 由于周期性面对是在最粗的网格层次上找到的，因此有必要确保最粗层次的面有正确的边界指标设置。一般来说，这意味着在进行任何全局或局部网格细化之前，必须首先在粗大的网格上设置所有的边界指标。
   *
   */
  template <typename MeshType>
  void
  collect_periodic_faces(
    const MeshType &         mesh,
    const types::boundary_id b_id1,
    const types::boundary_id b_id2,
    const int                direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
      &                                         matched_pairs,
    const Tensor<1, MeshType::space_dimension> &offset =
      dealii::Tensor<1, MeshType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());


  /**
   * @note 这个版本的collect_periodic_faces()将不会在单元格不在 @ref GlossFaceOrientation "标准方向 "
   * 的网格上工作。
   *
   */
  template <typename MeshType>
  void
  collect_periodic_faces(
    const MeshType &         mesh,
    const types::boundary_id b_id,
    const int                direction,
    std::vector<PeriodicFacePair<typename MeshType::cell_iterator>>
      &                                                 matched_pairs,
    const dealii::Tensor<1, MeshType::space_dimension> &offset =
      dealii::Tensor<1, MeshType::space_dimension>(),
    const FullMatrix<double> &matrix = FullMatrix<double>());

   /*@}*/ 
  /**
   * @name  处理边界和流形ID的问题
   *
   */
   /*@{*/ 

  /**
   * 将边界ID复制到边界上的面和边的流形ID。新的三角测量对象的默认manifold_id是
   * numbers::flat_manifold_id.
   * 这个函数将边界面和边的边界_id复制到相同面和边的manifold_id上，允许用户改变边界_id并将其用于边界条件，而不考虑几何形状，这将使用manifold_id来创建新点。只有活动单元会被迭代。当你的三角网格上只有一个活动层时，通常会调用这个函数。然后，网格细化会将这些指标继承到子单元、面和边上。
   * 可选参数 @p reset_boundary_ids,
   * 表示该函数在将边界面和边的值复制到manifold_id后，是否应将其重置为默认值0。默认情况下，boundary_ids不做任何改动。
   * @ingroup manifold   @relatesalso boundary
   *
   */
  template <int dim, int spacedim>
  void
  copy_boundary_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool reset_boundary_ids = false);

  /**
   * 将给定的边界id映射到边界上的面和边上的给定流形id。
   * 该函数将参数 @p src_boundary_ids
   * 中存在的边界面和边的边界id复制到 @p dst_manifold_ids,
   * 中相同面和边的相应流形id。    如果可选的参数 @p
   * reset_boundary_ids 非空， @p src_boundary_ids,
   * 中的每个边界id将被替换为 @p reset_boundary_ids.
   * 中的相应边界id。
   * 如果输入向量的大小不匹配，将抛出异常。如果 @p
   * src_boundary_ids
   * 中指出的边界ID不存在于三角形中，在这个过程中会被简单地忽略。
   * @ingroup manifold   @relatesalso boundary
   *
   */
  template <int dim, int spacedim>
  void
  map_boundary_to_manifold_ids(
    const std::vector<types::boundary_id> &src_boundary_ids,
    const std::vector<types::manifold_id> &dst_manifold_ids,
    Triangulation<dim, spacedim> &         tria,
    const std::vector<types::boundary_id> &reset_boundary_ids = {});

  /**
   * 将材料id复制到流形id。新Triangulation对象的默认manifold_id为
   * numbers::flat_manifold_id.
   * 当细化发生时，Triangulation会询问在何处将新点定位到底层流形。
   * 当从支持的输入格式中读取Triangulation时，可以存储在文件中的典型信息是边界面的边界条件（我们将其存储在面的边界_id中）、单元的材料类型（我们将其存储在单元的材料_id中）以及在某些情况下单元的子域ID（我们将其存储在单元的子域_id中）。
   * 如果您将这些网格之一读入
   * Triangulation，您可能仍然希望使用输入文件中指定的
   * material_id 作为 manifold_id
   * 描述。在这种情况下，您可以将一个Manifold对象与内部单元格相关联，该对象将被Triangulation用于查询Manifold对象的新点。该函数对活动单元进行迭代，并将
   * material_ids 复制到 manifold_ids 中。    可选参数 @p
   * compute_face_ids,
   * 表示该函数是否也应设置面的manifold_ids（包括内部面和边界上的面）。如果设置为
   * "true"，那么每个面的manifold_id将等于周围manifold_ids的最小值，确保为三角结构的每个面选择一个唯一的manifold
   * id。默认情况下，不计算面的manifold_id。
   * @ingroup manifold
   *
   */
  template <int dim, int spacedim>
  void
  copy_material_to_manifold_id(Triangulation<dim, spacedim> &tria,
                               const bool compute_face_ids = false);

  /**
   * 将与三角形 @p tria
   * 的单元格相关的流形指标传播到其同维度的一和二对象。
   * 这个函数将面和边（包括内部和边界）的 @p manifold_id
   * 设置为 @p disambiguation_function
   * 方法返回的值，该方法与共享相同面或边的单元格的流形指标集合一起调用。
   * 默认情况下， @p disambiguation_function
   * 在集合尺寸大于1时（即无法根据相邻单元格的流形指标来决定一个面或边应该有什么流形指标时）返回
   * numbers::flat_manifold_id
   * ，在集合尺寸为1时（即所有相邻单元格和面有相同的流形指标时）返回集合中包含的流形指标。
   * 参数 @p overwrite_only_flat_manifold_ids
   * 允许您指定当一个面或一个边已经具有不同于
   * numbers::flat_manifold_id. 的流形指标时该如何处理。
   * 如果标志是 @p true, ，该边或面将保持其原始流形指标。
   * 如果是 @p false, ，那么这些面和边的流形指标也将根据
   * @p disambiguation_function. 的返回值进行设置。
   *
   */
  template <int dim, int spacedim>
  void
  assign_co_dimensional_manifold_indicators(
    Triangulation<dim, spacedim> &            tria,
    const std::function<types::manifold_id(
      const std::set<types::manifold_id> &)> &disambiguation_function =
      [](const std::set<types::manifold_id> &manifold_ids) {
        if (manifold_ids.size() == 1)
          return *manifold_ids.begin();
        else
          return numbers::flat_manifold_id;
      },
    bool overwrite_only_flat_manifold_ids = true);
   /*@}*/ 

  /**
   * 将函数对象提供的 @p DataType 类型的任意数据从本地拥有的单元交换到其他处理器上的幽灵单元。    在这个调用之后，你通常会从每个幽灵单元上收到 @p unpack 的数据，因为它是由拥有处理器上的 @p pack 提供的。  你是否收到某个鬼魂单元上的 @p unpack 的信息，取决于 @p pack 函数是否决定需要发送什么。它使用 std_cxx17::optional 机制来做：如果 std_cxx17::optional 函数的 @p pack 返回对象是空的，那么这意味着它所调用的本地所有单元不需要发送数据。在这种情况下， @p unpack 也不会在接收方与之对应的幽灵单元上被调用。另一方面，如果 std_cxx17::optional 对象不是空的，那么存储在它里面的数据将被发送到接收方，并通过它调用 @p unpack 函数。      @tparam  DataType 要通信的数据的类型。在许多情况下，这个数据类型不能被编译器推断出来，例如，如果你为这个函数的第二个和第三个参数提供lambda函数。在这种情况下，你必须明确指定 @p DataType 作为函数调用的一个模板参数。    @tparam  MeshType  @p mesh.   @param  mesh 类型的变量，满足 @ref ConceptMeshType  "MeshType概念 "
   * 的要求。    @param  pack
   * 对每个本地拥有的、在其他地方是幽灵单元的单元，将被调用的函数。如上所述，该函数可以返回一个类型为
   * @p DataType
   * 的常规数据对象，表示应该发送数据，或者返回一个空的
   * <code>std_cxx17::optional@<DataType@></code>
   * ，表示这个单元不需要发送任何东西。    @param  unpack
   * 对于每个发送了数据的ghost单元，即发送方的 @p pack
   * 函数返回一个非空的 std_cxx17::optional
   * 对象，将调用该函数。    然后， @p unpack
   * 函数会被拥有该单元的处理器所发送的数据调用。
   * @param  cell_filter
   * 只有在这个过滤函数返回值为`true'的单元才会被通信。在默认情况下，该函数对所有单元都返回true，因此，所有相关的单元都被传送。
   * <h4> An example </h4>
   * 下面是一个例子，显示了这个函数在具体环境中的使用。它取自于确保
   * @p active_fe_index
   * （一个无符号整数）从本地拥有的单元（人们可以在具有hp-capabilities的DoFHandler对象中设置它）传送到其他处理器上相应的ghost单元的代码，以确保人们也能在这些处理器上查询到正确的值。
   * @code
   * using active_cell_iterator =
   * typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator;
   * auto pack = [] (const active_cell_iterator &cell)
   *
   * -> unsigned int
   *           {
   *             return cell->active_fe_index();
   *           };
   *
   * auto unpack = [] (const active_cell_iterator &cell,
   *                 const unsigned int active_fe_index)
   *
   * -> void
   *             {
   *               cell->set_active_fe_index(active_fe_index);
   *             };
   *
   * GridTools::exchange_cell_data_to_ghosts<
   * unsigned int, dealii::DoFHandler<dim,spacedim>> (dof_handler,
   *                                                  pack,
   *                                                  unpack);
   * @endcode
   * 你会注意到 @p pack
   * lambda函数返回一个`无符号的int`，而不是
   * `std_cxx17::optional<unsigned
   * int>`。前者会自动转换为后者，意味着数据将总是被传送到另一个处理器。
   * (实际上， @p unpack
   * 函数需要更复杂一些，因为它不允许在幽灵单元上调用
   * DoFAccessor::set_active_fe_index() 。相反， @p unpack
   * 函数直接访问内部数据结构。但你会明白的
   *
   * - 该代码也可以通过与上述类似的调用来交换材料ID、用户索引、边界指示器或任何种类的其他数据）。)
   *
   */
  template <typename DataType, typename MeshType>
  void
  exchange_cell_data_to_ghosts(
    const MeshType &                                     mesh,
    const std::function<std_cxx17::optional<DataType>(
      const typename MeshType::active_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::active_cell_iterator &,
                             const DataType &)> &        unpack,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &cell_filter =
        always_return<typename MeshType::active_cell_iterator, bool>{true});

  /**
   * 交换由函数对象提供的 @p DataType
   * 类型的任意数据，从本地拥有的水平单元到其他进程上的幽灵水平单元。
   * 除了 exchange_cell_data_to_ghosts()
   * 的参数外，这个函数允许提供一个  @p cell_filter
   * 函数，它可以用来只交流有标记的单元。在默认情况下，所有相关单元都会被通信。
   *
   */
  template <typename DataType, typename MeshType>
  void
  exchange_cell_data_to_level_ghosts(
    const MeshType &                                    mesh,
    const std::function<std_cxx17::optional<DataType>(
      const typename MeshType::level_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::level_cell_iterator &,
                             const DataType &)> &       unpack,
    const std::function<bool(const typename MeshType::level_cell_iterator &)> &
      cell_filter = always_return<typename MeshType::level_cell_iterator, bool>{
        true});

  /* 与MPI通信器的所有处理器交换边界框的向量  @p local_bboxes.  这个函数的目的是交换边界框，描述用函数  GridTools::compute_mesh_predicate_bounding_box  获得的分布式三角形中本地拥有的单元。    输出向量的大小是MPI通信器的进程数：其第i个条目包含第i个进程的向量 @p local_bboxes 。 
*
*/
  template <int spacedim>
  std::vector<std::vector<BoundingBox<spacedim>>>
  exchange_local_bounding_boxes(
    const std::vector<BoundingBox<spacedim>> &local_bboxes,
    const MPI_Comm &                          mpi_communicator);

  /**
   * 在这个集体操作中，每个进程提供一个边界框的向量和一个通信器。
   * 所有这些向量被收集在每个进程中，组织在一个搜索树中，然后返回。
   * 我们的想法是，边界框的向量描述了每个进程上计算的相关属性，这也可能对其他进程有用。一个例子是，如果输入的边界盒向量对应于一个
   * @ref GlossLocallyOwnedCell
   * 对象的局部拥有的网格分区的覆盖（见
   * parallel::distributed::Triangulation
   * ）。虽然这些边界盒可能与其他进程的边界盒重叠，但如果试图找出这一点的进程有一个其他进程的边界盒列表，那么找到哪个进程拥有包围给定点的单元格就容易得多。
   * 返回的搜索树对象是一个带有打包算法的r-tree，由boost库提供。更多信息见https://www.boost.org/doc/libs/1_67_0/libs/geometry/doc/html/geometry/spatial_indexes/introduction.html。
   * 在返回的树中，每个节点都包含一对元素：第一个是一个边界框，第二个是其本地描述包含边界框的进程的等级。
   * @note  这个函数是一个集体操作。
   *
   */
  template <int spacedim>
  RTree<std::pair<BoundingBox<spacedim>, unsigned int>>
  build_global_description_tree(
    const std::vector<BoundingBox<spacedim>> &local_description,
    const MPI_Comm &                          mpi_communicator);

  /**
   * 对于一个给定的三角形，收集所有由于周期性而重合的本地相关顶点。
   * 重合的顶点被放入一个组，例如。[1, 25,
   * 51]，用其中的一个任意元素来标记，例如。"1".
   * 所有重合顶点都将标签存储到它的组中，这样它们就可以快速访问该组中的所有重合顶点：例如。51
   *
   * -> "1"
   *
   * --> [1, 25, 51]  @param[in]  tria Triangulation.     @param[out]
   * 重合顶点组（coinciding_vertex_groups）
   * 用其中的任意元素标注的等价类（重合顶点）的映射。
   * 不重合的顶点被忽略。    @param[out]
   * vertex_to_coinciding_vertex_group
   * 一个顶点到一组重合顶点的标签的映射。不包含在这个向量中的顶点不与任何其他顶点重合。
   *
   */
  template <int dim, int spacedim>
  void
  collect_coinciding_vertices(
    const Triangulation<dim, spacedim> &               tria,
    std::map<unsigned int, std::vector<unsigned int>> &coinciding_vertex_groups,
    std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group);

  /**
   * 返回一个地图，对于每个顶点，列出其子域与该顶点相邻的所有进程。
   * @param[in]  tria Triangulation。
   *
   */
  template <int dim, int spacedim>
  std::map<unsigned int, std::set<dealii::types::subdomain_id>>
  compute_vertices_with_ghost_neighbors(
    const Triangulation<dim, spacedim> &tria);

  /**
   * 一种结构，允许将 @p T
   * 类型的单元数据从一个处理器传输到另一个处理器。它相当于一个打包的缓冲区，存储了一个CellId的向量和一个
   * @p T.
   * 类型的向量。这个类通过提供保存/加载函数，能够把CellId的向量和
   * @p T 类型的相关数据打包成一个流，从而促进了传输。
   * 类型 @p T 被假定为可被 <code>boost::serialization</code>
   * 序列化（例如 <code>unsigned int</code> or
   * <code>std::vector@<double@></code>  ）。
   *
   */
  template <int dim, typename T>
  struct CellDataTransferBuffer
  {
    /**
     * 一个向量，用于存储要转移的单元格的ID。
     *
     */
    std::vector<CellId> cell_ids;

    /**
     * 一个要传输的单元格数据的向量。
     *
     */
    std::vector<T> data;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)
     * @pre
     * 将此对象的数据写入一个流，以便进行序列化。用户有责任保持
     * @p data 的大小与 @p cell_ids 的大小相等。
     *
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以便进行序列化。
     * 扔掉之前的内容。
     *
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#  ifdef DOXYGEN
    /**
     * 为了序列化的目的，使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中。
     *
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int version);
#  else
    // This macro defines the serialize() method that is compatible with
    // the templated save() and load() method that have been implemented.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#  endif
  };



  /**
   * 行进方程（2D）和行进立方程算法的实现，用于创建数据结构（Point和CellData的向量），在标量场的等高线/轮廓上创建线性/双线性表面网格。
   * 为了提高等值线/轮廓的近似度和产生的线性表面网格，可以增加细分的数量，使算法不是在一个单元上运行，而是在子单元上运行，顶点值由单元值插值而成。
   * @note  得到的网格将包含二维的线和三维的三角形。
   * @note
   * 生成的网格将不是高质量的，因为如果网格在靠近顶点处被切割，可能会包含直径非常小的单元。
   *
   */
  template <int dim, typename VectorType>
  class MarchingCubeAlgorithm
  {
  public:
    /**
     * 矢量的数值类型。
     *
     */
    using value_type = typename VectorType::value_type;

    /**
     * 构造函数。
     *
     */
    MarchingCubeAlgorithm(const Mapping<dim, dim> &      mapping,
                          const FiniteElement<dim, dim> &fe,
                          const unsigned int             n_subdivisions = 1,
                          const double                   tolerance = 1e-10);

    /**
     * 处理所有本地拥有的单元格并为所有被切割的单元格填充
     * @p vertices 和 @p cells 。
     *
     */
    void
    process(const DoFHandler<dim> &         background_dof_handler,
            const VectorType &              ls_vector,
            const double                    iso_level,
            std::vector<Point<dim>> &       vertices,
            std::vector<CellData<dim - 1>> &cells) const;

    /**
     * 处理所提供的单元格并为所有被切割的单元格填充 @p
     * vertices 和 @p cells 。
     * @note 如果单元格没有被切割，产生的向量为空。
     *
     */
    void
    process_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 const VectorType &              ls_vector,
                 const double                    iso_level,
                 std::vector<Point<dim>> &       vertices,
                 std::vector<CellData<dim - 1>> &cells) const;

  private:
    /**
     * 内部函数创建一个具有n_subdivisions+1同等位置正交点的正交规则。
     *
     */
    static Quadrature<dim>
    create_quadrature_rule(const unsigned int n_subdivisions);

    /**
     * 处理一个单元。
     *
     */
    void
    process_cell(std::vector<value_type> &       ls_values,
                 const std::vector<Point<dim>> & points,
                 const double                    iso_level,
                 std::vector<Point<dim>> &       vertices,
                 std::vector<CellData<dim - 1>> &cells) const;

    /**
     * 处理一个子单元（2D）。
     * @note
     * 带有鞍状点的子单元被忽略。在这种情况下，请增加子单元的数量。
     *
     */
    void
    process_sub_cell(const std::vector<value_type> & ls_values,
                     const std::vector<Point<2>> &   points,
                     const std::vector<unsigned int> mask,
                     const double                    iso_level,
                     std::vector<Point<2>> &         vertices,
                     std::vector<CellData<1>> &      cells) const;

    /**
     * 处理一个子单元（3D）。
     *
     */
    void
    process_sub_cell(const std::vector<value_type> & ls_values,
                     const std::vector<Point<3>> &   points,
                     const std::vector<unsigned int> mask,
                     const double                    iso_level,
                     std::vector<Point<3>> &         vertices,
                     std::vector<CellData<2>> &      cells) const;

    /**
     * 每个单元在每个方向上被细分的数量，以提高近似度。
     *
     */
    const unsigned int n_subdivisions;

    /**
     * 绝对公差，指定顶点和切割点之间的最小距离，从而使一条线被认为是被切割的。
     *
     */
    const double tolerance;

    /**
     * 内部使用的FEValues，并以正确的细分数设置正交规则。
     *
     */
    mutable FEValues<dim> fe_values;
  };



  /**
   * @name  异常情况
   *
   */
   /*@{*/ 

  /**
   * 例外情况
   *
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 int,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcNonExistentSubdomain,
                 int,
                 << "The subdomain id " << arg1
                 << " has no cells associated with it.");
  /**
   * 异常情况
   *
   */
  DeclException0(ExcTriangulationHasBeenRefined);

  /**
   * 异常情况
   *
   */
  DeclException1(ExcScalingFactorNotPositive,
                 double,
                 << "The scaling factor must be positive, but it is " << arg1
                 << ".");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcVertexNotUsed,
                 unsigned int,
                 << "The given vertex with index " << arg1
                 << " is not used in the given triangulation.");

   /*@}*/ 

}  /*namespace GridTools*/ 


/**
 * 当一个网格的边缘不能被定向时，就会抛出一个异常。
 *
 *
 * @note
 * 为了向后兼容旧的GridReordering类，这个异常不在GridTools命名空间。
 *
 *
 * @ingroup Exceptions
 *
 */
DeclExceptionMsg(ExcMeshNotOrientable,
                 "The edges of the mesh are not consistently orientable.");



 /* ----------------- Template function --------------- */ 

#  ifndef DOXYGEN

namespace GridTools
{
  template <int dim>
  double
  cell_measure(
    const std::vector<Point<dim>> &all_vertices,
    const unsigned int (&indices)[GeometryInfo<dim>::vertices_per_cell])
  {
    // We forward call to the ArrayView version:
    const ArrayView<const unsigned int> view(
      indices, GeometryInfo<dim>::vertices_per_cell);
    return cell_measure(all_vertices, view);
  }



  // This specialization is defined here so that the general template in the
  // source file doesn't need to have further 1D overloads for the internal
  // functions it calls.
  template <>
  inline Triangulation<1, 1>::DistortedCellList
  fix_up_distorted_child_cells(const Triangulation<1, 1>::DistortedCellList &,
                               Triangulation<1, 1> &)
  {
    return {};
  }



  template <int dim, typename Predicate, int spacedim>
  void
  transform(const Predicate &             predicate,
            Triangulation<dim, spacedim> &triangulation)
  {
    std::vector<bool> treated_vertices(triangulation.n_vertices(), false);

    // loop over all active cells, and
    // transform those vertices that
    // have not yet been touched. note
    // that we get to all vertices in
    // the triangulation by only
    // visiting the active cells.
    typename Triangulation<dim, spacedim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (const unsigned int v : cell->vertex_indices())
        if (treated_vertices[cell->vertex_index(v)] == false)
          {
            // transform this vertex
            cell->vertex(v) = predicate(cell->vertex(v));
            // and mark it as treated
            treated_vertices[cell->vertex_index(v)] = true;
          };


    // now fix any vertices on hanging nodes so that we don't create any holes
    if (dim == 2)
      {
        typename Triangulation<dim, spacedim>::active_cell_iterator
          cell = triangulation.begin_active(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                Assert(cell->reference_cell() ==
                         ReferenceCells::get_hypercube<dim>(),
                       ExcNotImplemented());

                // this line has children
                cell->face(face)->child(0)->vertex(1) =
                  (cell->face(face)->vertex(0) + cell->face(face)->vertex(1)) /
                  2;
              }
      }
    else if (dim == 3)
      {
        typename Triangulation<dim, spacedim>::active_cell_iterator
          cell = triangulation.begin_active(),
          endc = triangulation.end();
        for (; cell != endc; ++cell)
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->has_children() &&
                !cell->face(face)->at_boundary())
              {
                Assert(cell->reference_cell() ==
                         ReferenceCells::get_hypercube<dim>(),
                       ExcNotImplemented());

                // this face has hanging nodes
                cell->face(face)->child(0)->vertex(1) =
                  (cell->face(face)->vertex(0) + cell->face(face)->vertex(1)) /
                  2.0;
                cell->face(face)->child(0)->vertex(2) =
                  (cell->face(face)->vertex(0) + cell->face(face)->vertex(2)) /
                  2.0;
                cell->face(face)->child(1)->vertex(3) =
                  (cell->face(face)->vertex(1) + cell->face(face)->vertex(3)) /
                  2.0;
                cell->face(face)->child(2)->vertex(3) =
                  (cell->face(face)->vertex(2) + cell->face(face)->vertex(3)) /
                  2.0;

                // center of the face
                cell->face(face)->child(0)->vertex(3) =
                  (cell->face(face)->vertex(0) + cell->face(face)->vertex(1) +
                   cell->face(face)->vertex(2) + cell->face(face)->vertex(3)) /
                  4.0;
              }
      }

    // Make sure FEValues notices that the mesh has changed
    triangulation.signals.mesh_movement();
  }



  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  get_active_child_cells(const typename MeshType::cell_iterator &cell)
  {
    std::vector<typename MeshType::active_cell_iterator> child_cells;

    if (cell->has_children())
      {
        for (unsigned int child = 0; child < cell->n_children(); ++child)
          if (cell->child(child)->has_children())
            {
              const std::vector<typename MeshType::active_cell_iterator>
                children = get_active_child_cells<MeshType>(cell->child(child));
              child_cells.insert(child_cells.end(),
                                 children.begin(),
                                 children.end());
            }
          else
            child_cells.push_back(cell->child(child));
      }

    return child_cells;
  }



  template <class MeshType>
  void
  get_active_neighbors(
    const typename MeshType::active_cell_iterator &       cell,
    std::vector<typename MeshType::active_cell_iterator> &active_neighbors)
  {
    active_neighbors.clear();
    for (const unsigned int n : cell->face_indices())
      if (!cell->at_boundary(n))
        {
          if (MeshType::dimension == 1)
            {
              // check children of neighbor. note
              // that in 1d children of the neighbor
              // may be further refined. In 1d the
              // case is simple since we know what
              // children bound to the present cell
              typename MeshType::cell_iterator neighbor_child =
                cell->neighbor(n);
              if (!neighbor_child->is_active())
                {
                  while (neighbor_child->has_children())
                    neighbor_child = neighbor_child->child(n == 0 ? 1 : 0);

                  Assert(neighbor_child->neighbor(n == 0 ? 1 : 0) == cell,
                         ExcInternalError());
                }
              active_neighbors.push_back(neighbor_child);
            }
          else
            {
              if (cell->face(n)->has_children())
                // this neighbor has children. find
                // out which border to the present
                // cell
                for (unsigned int c = 0;
                     c < cell->face(n)->n_active_descendants();
                     ++c)
                  active_neighbors.push_back(
                    cell->neighbor_child_on_subface(n, c));
              else
                {
                  // the neighbor must be active
                  // himself
                  Assert(cell->neighbor(n)->is_active(), ExcInternalError());
                  active_neighbors.push_back(cell->neighbor(n));
                }
            }
        }
  }



  namespace internal
  {
    namespace ProjectToObject
    {
      /**
       * 方法 GridTools::project_to_object
       * 需要沿着单线的表面求导数。一般来说，这些导数不能用有限差分来近似，而是用df/dx_i形式的特殊差分。
       *
       * - df/dx_j  <em> 可以 </em> 被逼近。这 <code>struct</code> 只是存储了由模版近似的两个导数（在上面的例子中 <code>i</code> and <code>j</code> ）。
       *
       */
      struct CrossDerivative
      {
        const unsigned int direction_0;
        const unsigned int direction_1;

        CrossDerivative(const unsigned int d0, const unsigned int d1);
      };

      inline CrossDerivative::CrossDerivative(const unsigned int d0,
                                              const unsigned int d1)
        : direction_0(d0)
        , direction_1(d1)
      {}



      /**
       * 用两点居中的方案对一阶导数进行标准的二阶逼近。这在下面的一维牛顿方法中使用。
       *
       */
      template <typename F>
      inline auto
      centered_first_difference(const double center,
                                const double step,
                                const F &f) -> decltype(f(center) - f(center))
      {
        return (f(center + step) - f(center - step)) / (2.0 * step);
      }



      /**
       * 二阶导数的标准二阶近似，采用三点为中心的方案。这在下面的一维牛顿方法中使用。
       *
       */
      template <typename F>
      inline auto
      centered_second_difference(const double center,
                                 const double step,
                                 const F &f) -> decltype(f(center) - f(center))
      {
        return (f(center + step) - 2.0 * f(center) + f(center - step)) /
               (step * step);
      }



      /**
       * 导数的四阶近似 df/dx_i
       *
       * - df/dx_j 其中 <code>i</code> and <code>j</code> 是由 @p cross_derivative指定。导数近似在 @p center ，步长为 @p step ，函数为 @p f.  。
       *
       */
      template <int structdim, typename F>
      inline auto
      cross_stencil(
        const CrossDerivative cross_derivative,
        const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &center,
        const double                                                 step,
        const F &f) -> decltype(f(center) - f(center))
      {
        Tensor<1, GeometryInfo<structdim>::vertices_per_cell> simplex_vector;
        simplex_vector[cross_derivative.direction_0] = 0.5 * step;
        simplex_vector[cross_derivative.direction_1] = -0.5 * step;
        return (-4.0 * f(center) - 1.0 * f(center + simplex_vector) -
                1.0 / 3.0 * f(center - simplex_vector) +
                16.0 / 3.0 * f(center + 0.5 * simplex_vector)) /
               step;
      }



      /**
       * 在 GridTools::project_to_object
       * 中使用的优化算法本质上是一种梯度下降法。这个函数计算目标函数梯度的条目；更多信息请参见
       * GridTools::project_to_object 里面的注释描述。
       *
       */
      template <int spacedim, int structdim, typename F>
      inline double
      gradient_entry(
        const unsigned int     row_n,
        const unsigned int     dependent_direction,
        const Point<spacedim> &p0,
        const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &center,
        const double                                                 step,
        const F &                                                    f)
      {
        Assert(row_n < GeometryInfo<structdim>::vertices_per_cell &&
                 dependent_direction <
                   GeometryInfo<structdim>::vertices_per_cell,
               ExcMessage("This function assumes that the last weight is a "
                          "dependent variable (and hence we cannot take its "
                          "derivative directly)."));
        Assert(row_n != dependent_direction,
               ExcMessage(
                 "We cannot differentiate with respect to the variable "
                 "that is assumed to be dependent."));

        const Point<spacedim>     manifold_point = f(center);
        const Tensor<1, spacedim> stencil_value  = cross_stencil<structdim>(
          {row_n, dependent_direction}, center, step, f);
        double entry = 0.0;
        for (unsigned int dim_n = 0; dim_n < spacedim; ++dim_n)
          entry +=
            -2.0 * (p0[dim_n] - manifold_point[dim_n]) * stencil_value[dim_n];
        return entry;
      }

      /**
       * 投射到一个d-线性对象上。这比project_to_object中的一般算法更精确，但只适用于由线性、双线性或三线性映射描述的几何图形。
       *
       */
      template <typename Iterator, int spacedim, int structdim>
      Point<spacedim>
      project_to_d_linear_object(const Iterator &       object,
                                 const Point<spacedim> &trial_point)
      {
        // let's look at this for simplicity for a quad (structdim==2) in a
        // space with spacedim>2 (notate trial_point by y): all points on the
        // surface are given by
        //   x(\xi) = sum_i v_i phi_x(\xi)
        // where v_i are the vertices of the quad, and \xi=(\xi_1,\xi_2) are the
        // reference coordinates of the quad. so what we are trying to do is
        // find a point x on the surface that is closest to the point y. there
        // are different ways to solve this problem, but in the end it's a
        // nonlinear problem and we have to find reference coordinates \xi so
        // that J(\xi) = 1/2 || x(\xi)-y ||^2 is minimal. x(\xi) is a function
        // that is structdim-linear in \xi, so J(\xi) is a polynomial of degree
        // 2*structdim that we'd like to minimize. unless structdim==1, we'll
        // have to use a Newton method to find the answer. This leads to the
        // following formulation of Newton steps:
        //
        // Given \xi_k, find \delta\xi_k so that
        //   H_k \delta\xi_k = - F_k
        // where H_k is an approximation to the second derivatives of J at
        // \xi_k, and F_k is the first derivative of J.  We'll iterate this a
        // number of times until the right hand side is small enough. As a
        // stopping criterion, we terminate if ||\delta\xi||<eps.
        //
        // As for the Hessian, the best choice would be
        //   H_k = J''(\xi_k)
        // but we'll opt for the simpler Gauss-Newton form
        //   H_k = A^T A
        // i.e.
        //   (H_k)_{nm} = \sum_{i,j} v_i*v_j *
        //                   \partial_n phi_i *
        //                   \partial_m phi_j
        // we start at xi=(0.5, 0.5).
        Point<structdim> xi;
        for (unsigned int d = 0; d < structdim; ++d)
          xi[d] = 0.5;

        Point<spacedim> x_k;
        for (const unsigned int i : GeometryInfo<structdim>::vertex_indices())
          x_k += object->vertex(i) *
                 GeometryInfo<structdim>::d_linear_shape_function(xi, i);

        do
          {
            Tensor<1, structdim> F_k;
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              F_k +=
                (x_k - trial_point) * object->vertex(i) *
                GeometryInfo<structdim>::d_linear_shape_function_gradient(xi,
                                                                          i);

            Tensor<2, structdim> H_k;
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              for (const unsigned int j :
                   GeometryInfo<structdim>::vertex_indices())
                {
                  Tensor<2, structdim> tmp = outer_product(
                    GeometryInfo<structdim>::d_linear_shape_function_gradient(
                      xi, i),
                    GeometryInfo<structdim>::d_linear_shape_function_gradient(
                      xi, j));
                  H_k += (object->vertex(i) * object->vertex(j)) * tmp;
                }

            const Tensor<1, structdim> delta_xi = -invert(H_k) * F_k;
            xi += delta_xi;

            x_k = Point<spacedim>();
            for (const unsigned int i :
                 GeometryInfo<structdim>::vertex_indices())
              x_k += object->vertex(i) *
                     GeometryInfo<structdim>::d_linear_shape_function(xi, i);

            if (delta_xi.norm() < 1e-7)
              break;
          }
        while (true);

        return x_k;
      }
    } // namespace ProjectToObject
  }   // namespace internal



  namespace internal
  {
    // We hit an internal compiler error in ICC 15 if we define this as a lambda
    // inside the project_to_object function below.
    template <int structdim>
    inline bool
    weights_are_ok(
      const Tensor<1, GeometryInfo<structdim>::vertices_per_cell> &v)
    {
      // clang has trouble figuring out structdim here, so define it
      // again:
      static const std::size_t n_vertices_per_cell =
        Tensor<1, GeometryInfo<structdim>::vertices_per_cell>::
          n_independent_components;
      std::array<double, n_vertices_per_cell> copied_weights;
      for (unsigned int i = 0; i < n_vertices_per_cell; ++i)
        {
          copied_weights[i] = v[i];
          if (v[i] < 0.0 || v[i] > 1.0)
            return false;
        }

      // check the sum: try to avoid some roundoff errors by summing in order
      std::sort(copied_weights.begin(), copied_weights.end());
      const double sum =
        std::accumulate(copied_weights.begin(), copied_weights.end(), 0.0);
      return std::abs(sum - 1.0) < 1e-10; // same tolerance used in manifold.cc
    }
  } // namespace internal

  template <typename Iterator>
  Point<Iterator::AccessorType::space_dimension>
  project_to_object(
    const Iterator &                                      object,
    const Point<Iterator::AccessorType::space_dimension> &trial_point)
  {
    const int spacedim  = Iterator::AccessorType::space_dimension;
    const int structdim = Iterator::AccessorType::structure_dimension;

    Point<spacedim> projected_point = trial_point;

    if (structdim >= spacedim)
      return projected_point;
    else if (structdim == 1 || structdim == 2)
      {
        using namespace internal::ProjectToObject;
        // Try to use the special flat algorithm for quads (this is better
        // than the general algorithm in 3D). This does not take into account
        // whether projected_point is outside the quad, but we optimize along
        // lines below anyway:
        const int                      dim = Iterator::AccessorType::dimension;
        const Manifold<dim, spacedim> &manifold = object->get_manifold();
        if (structdim == 2 && dynamic_cast<const FlatManifold<dim, spacedim> *>(
                                &manifold) != nullptr)
          {
            projected_point =
              project_to_d_linear_object<Iterator, spacedim, structdim>(
                object, trial_point);
          }
        else
          {
            // We want to find a point on the convex hull (defined by the
            // vertices of the object and the manifold description) that is
            // relatively close to the trial point. This has a few issues:
            //
            // 1. For a general convex hull we are not guaranteed that a unique
            //    minimum exists.
            // 2. The independent variables in the optimization process are the
            //    weights given to Manifold::get_new_point, which must sum to 1,
            //    so we cannot use standard finite differences to approximate a
            //    gradient.
            //
            // There is not much we can do about 1., but for 2. we can derive
            // finite difference stencils that work on a structdim-dimensional
            // simplex and rewrite the optimization problem to use those
            // instead. Consider the structdim 2 case and let
            //
            // F(c0, c1, c2, c3) = Manifold::get_new_point(vertices, {c0, c1,
            // c2, c3})
            //
            // where {c0, c1, c2, c3} are the weights for the four vertices on
            // the quadrilateral. We seek to minimize the Euclidean distance
            // between F(...) and trial_point. We can solve for c3 in terms of
            // the other weights and get, for one coordinate direction
            //
            // d/dc0 ((x0 - F(c0, c1, c2, 1 - c0 - c1 - c2))^2)
            //      = -2(x0 - F(...)) (d/dc0 F(...) - d/dc3 F(...))
            //
            // where we substitute back in for c3 after taking the
            // derivative. We can compute a stencil for the cross derivative
            // d/dc0 - d/dc3: this is exactly what cross_stencil approximates
            // (and gradient_entry computes the sum over the independent
            // variables). Below, we somewhat arbitrarily pick the last
            // component as the dependent one.
            //
            // Since we can now calculate derivatives of the objective
            // function we can use gradient descent to minimize it.
            //
            // Of course, this is much simpler in the structdim = 1 case (we
            // could rewrite the projection as a 1D optimization problem), but
            // to reduce the potential for bugs we use the same code in both
            // cases.
            const double step_size = object->diameter() / 64.0;

            constexpr unsigned int n_vertices_per_cell =
              GeometryInfo<structdim>::vertices_per_cell;

            std::array<Point<spacedim>, n_vertices_per_cell> vertices;
            for (unsigned int vertex_n = 0; vertex_n < n_vertices_per_cell;
                 ++vertex_n)
              vertices[vertex_n] = object->vertex(vertex_n);

            auto get_point_from_weights =
              [&](const Tensor<1, n_vertices_per_cell> &weights)
              -> Point<spacedim> {
              return object->get_manifold().get_new_point(
                make_array_view(vertices.begin(), vertices.end()),
                make_array_view(weights.begin_raw(), weights.end_raw()));
            };

            // pick the initial weights as (normalized) inverse distances from
            // the trial point:
            Tensor<1, n_vertices_per_cell> guess_weights;
            double                         guess_weights_sum = 0.0;
            for (unsigned int vertex_n = 0; vertex_n < n_vertices_per_cell;
                 ++vertex_n)
              {
                const double distance =
                  vertices[vertex_n].distance(trial_point);
                if (distance == 0.0)
                  {
                    guess_weights           = 0.0;
                    guess_weights[vertex_n] = 1.0;
                    guess_weights_sum       = 1.0;
                    break;
                  }
                else
                  {
                    guess_weights[vertex_n] = 1.0 / distance;
                    guess_weights_sum += guess_weights[vertex_n];
                  }
              }
            guess_weights /= guess_weights_sum;
            Assert(internal::weights_are_ok<structdim>(guess_weights),
                   ExcInternalError());

            // The optimization algorithm consists of two parts:
            //
            // 1. An outer loop where we apply the gradient descent algorithm.
            // 2. An inner loop where we do a line search to find the optimal
            //    length of the step one should take in the gradient direction.
            //
            for (unsigned int outer_n = 0; outer_n < 40; ++outer_n)
              {
                const unsigned int dependent_direction =
                  n_vertices_per_cell - 1;
                Tensor<1, n_vertices_per_cell> current_gradient;
                for (unsigned int row_n = 0; row_n < n_vertices_per_cell;
                     ++row_n)
                  {
                    if (row_n != dependent_direction)
                      {
                        current_gradient[row_n] =
                          gradient_entry<spacedim, structdim>(
                            row_n,
                            dependent_direction,
                            trial_point,
                            guess_weights,
                            step_size,
                            get_point_from_weights);

                        current_gradient[dependent_direction] -=
                          current_gradient[row_n];
                      }
                  }

                // We need to travel in the -gradient direction, as noted
                // above, but we may not want to take a full step in that
                // direction; instead, guess that we will go -0.5*gradient and
                // do quasi-Newton iteration to pick the best multiplier. The
                // goal is to find a scalar alpha such that
                //
                // F(x - alpha g)
                //
                // is minimized, where g is the gradient and F is the
                // objective function. To find the optimal value we find roots
                // of the derivative of the objective function with respect to
                // alpha by Newton iteration, where we approximate the first
                // and second derivatives of F(x - alpha g) with centered
                // finite differences.
                double gradient_weight = -0.5;
                auto   gradient_weight_objective_function =
                  [&](const double gradient_weight_guess) -> double {
                  return (trial_point -
                          get_point_from_weights(guess_weights +
                                                 gradient_weight_guess *
                                                   current_gradient))
                    .norm_square();
                };

                for (unsigned int inner_n = 0; inner_n < 10; ++inner_n)
                  {
                    const double update_numerator = centered_first_difference(
                      gradient_weight,
                      step_size,
                      gradient_weight_objective_function);
                    const double update_denominator =
                      centered_second_difference(
                        gradient_weight,
                        step_size,
                        gradient_weight_objective_function);

                    // avoid division by zero. Note that we limit the gradient
                    // weight below
                    if (std::abs(update_denominator) == 0.0)
                      break;
                    gradient_weight =
                      gradient_weight - update_numerator / update_denominator;

                    // Put a fairly lenient bound on the largest possible
                    // gradient (things tend to be locally flat, so the gradient
                    // itself is usually small)
                    if (std::abs(gradient_weight) > 10)
                      {
                        gradient_weight = -10.0;
                        break;
                      }
                  }

                // It only makes sense to take convex combinations with weights
                // between zero and one. If the update takes us outside of this
                // region then rescale the update to stay within the region and
                // try again
                Tensor<1, n_vertices_per_cell> tentative_weights =
                  guess_weights + gradient_weight * current_gradient;

                double new_gradient_weight = gradient_weight;
                for (unsigned int iteration_count = 0; iteration_count < 40;
                     ++iteration_count)
                  {
                    if (internal::weights_are_ok<structdim>(tentative_weights))
                      break;

                    for (unsigned int i = 0; i < n_vertices_per_cell; ++i)
                      {
                        if (tentative_weights[i] < 0.0)
                          {
                            tentative_weights -=
                              (tentative_weights[i] / current_gradient[i]) *
                              current_gradient;
                          }
                        if (tentative_weights[i] < 0.0 ||
                            1.0 < tentative_weights[i])
                          {
                            new_gradient_weight /= 2.0;
                            tentative_weights =
                              guess_weights +
                              new_gradient_weight * current_gradient;
                          }
                      }
                  }

                // the update might still send us outside the valid region, so
                // check again and quit if the update is still not valid
                if (!internal::weights_are_ok<structdim>(tentative_weights))
                  break;

                // if we cannot get closer by traveling in the gradient
                // direction then quit
                if (get_point_from_weights(tentative_weights)
                      .distance(trial_point) <
                    get_point_from_weights(guess_weights).distance(trial_point))
                  guess_weights = tentative_weights;
                else
                  break;
                Assert(internal::weights_are_ok<structdim>(guess_weights),
                       ExcInternalError());
              }
            Assert(internal::weights_are_ok<structdim>(guess_weights),
                   ExcInternalError());
            projected_point = get_point_from_weights(guess_weights);
          }

        // if structdim == 2 and the optimal point is not on the interior then
        // we may be able to get a more accurate result by projecting onto the
        // lines.
        if (structdim == 2)
          {
            std::array<Point<spacedim>, GeometryInfo<structdim>::lines_per_cell>
              line_projections;
            for (unsigned int line_n = 0;
                 line_n < GeometryInfo<structdim>::lines_per_cell;
                 ++line_n)
              {
                line_projections[line_n] =
                  project_to_object(object->line(line_n), trial_point);
              }
            std::sort(line_projections.begin(),
                      line_projections.end(),
                      [&](const Point<spacedim> &a, const Point<spacedim> &b) {
                        return a.distance(trial_point) <
                               b.distance(trial_point);
                      });
            if (line_projections[0].distance(trial_point) <
                projected_point.distance(trial_point))
              projected_point = line_projections[0];
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
        return projected_point;
      }

    return projected_point;
  }



  template <int dim, typename T>
  template <class Archive>
  void
  CellDataTransferBuffer<dim, T>::save(Archive &ar,
                                       const unsigned int  /*version*/ ) const
  {
    Assert(cell_ids.size() == data.size(),
           ExcDimensionMismatch(cell_ids.size(), data.size()));
    // archive the cellids in an efficient binary format
    const std::size_t n_cells = cell_ids.size();
    ar &              n_cells;
    for (const auto &id : cell_ids)
      {
        CellId::binary_type binary_cell_id = id.template to_binary<dim>();
        ar &                binary_cell_id;
      }

    ar &data;
  }



  template <int dim, typename T>
  template <class Archive>
  void
  CellDataTransferBuffer<dim, T>::load(Archive &ar,
                                       const unsigned int  /*version*/ )
  {
    std::size_t n_cells;
    ar &        n_cells;
    cell_ids.clear();
    cell_ids.reserve(n_cells);
    for (unsigned int c = 0; c < n_cells; ++c)
      {
        CellId::binary_type value;
        ar &                value;
        cell_ids.emplace_back(value);
      }
    ar &data;
  }


  namespace internal
  {
    template <typename DataType,
              typename MeshType,
              typename MeshCellIteratorType>
    inline void
    exchange_cell_data(
      const MeshType &mesh,
      const std::function<
        std_cxx17::optional<DataType>(const MeshCellIteratorType &)> &pack,
      const std::function<void(const MeshCellIteratorType &, const DataType &)>
        &                                                         unpack,
      const std::function<bool(const MeshCellIteratorType &)> &   cell_filter,
      const std::function<void(
        const std::function<void(const MeshCellIteratorType &,
                                 const types::subdomain_id)> &)> &process_cells,
      const std::function<std::set<types::subdomain_id>(
        const parallel::TriangulationBase<MeshType::dimension,
                                          MeshType::space_dimension> &)>
        &compute_ghost_owners)
    {
#    ifndef DEAL_II_WITH_MPI
      (void)mesh;
      (void)pack;
      (void)unpack;
      (void)cell_filter;
      (void)process_cells;
      (void)compute_ghost_owners;
      Assert(false, ExcNeedsMPI());
#    else
      constexpr int dim      = MeshType::dimension;
      constexpr int spacedim = MeshType::space_dimension;
      auto          tria =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
          &mesh.get_triangulation());
      Assert(
        tria != nullptr,
        ExcMessage(
          "The function exchange_cell_data_to_ghosts() only works with parallel triangulations."));

      if (const auto tria = dynamic_cast<
            const parallel::shared::Triangulation<dim, spacedim> *>(
            &mesh.get_triangulation()))
        {
          Assert(
            tria->with_artificial_cells(),
            ExcMessage(
              "The functions GridTools::exchange_cell_data_to_ghosts() and "
              "GridTools::exchange_cell_data_to_level_ghosts() can only "
              "operate on a single layer ghost cells. However, you have "
              "given a Triangulation object of type "
              "parallel::shared::Triangulation without artificial cells "
              "resulting in arbitrary numbers of ghost layers."));
        }

      // build list of cells to request for each neighbor
      std::set<dealii::types::subdomain_id> ghost_owners =
        compute_ghost_owners(*tria);
      std::map<dealii::types::subdomain_id,
               std::vector<typename CellId::binary_type>>
        neighbor_cell_list;

      for (const auto ghost_owner : ghost_owners)
        neighbor_cell_list[ghost_owner] = {};

      process_cells([&](const auto &cell, const auto key) {
        if (cell_filter(cell))
          neighbor_cell_list[key].emplace_back(
            cell->id().template to_binary<spacedim>());
      });

      Assert(ghost_owners.size() == neighbor_cell_list.size(),
             ExcInternalError());


      // Before sending & receiving, make sure we protect this section with
      // a mutex:
      static Utilities::MPI::CollectiveMutex      mutex;
      Utilities::MPI::CollectiveMutex::ScopedLock lock(
        mutex, tria->get_communicator());

      const int mpi_tag =
        Utilities::MPI::internal::Tags::exchange_cell_data_request;
      const int mpi_tag_reply =
        Utilities::MPI::internal::Tags::exchange_cell_data_reply;

      // send our requests:
      std::vector<MPI_Request> requests(ghost_owners.size());
      {
        unsigned int idx = 0;
        for (const auto &it : neighbor_cell_list)
          {
            // send the data about the relevant cells
            const int ierr = MPI_Isend(it.second.data(),
                                       it.second.size() * sizeof(it.second[0]),
                                       MPI_BYTE,
                                       it.first,
                                       mpi_tag,
                                       tria->get_communicator(),
                                       &requests[idx]);
            AssertThrowMPI(ierr);
            ++idx;
          }
      }

      using DestinationToBufferMap =
        std::map<dealii::types::subdomain_id,
                 GridTools::CellDataTransferBuffer<dim, DataType>>;
      DestinationToBufferMap destination_to_data_buffer_map;

      // receive requests and reply with the ghost indices
      std::vector<std::vector<typename CellId::binary_type>> cell_data_to_send(
        ghost_owners.size());
      std::vector<std::vector<types::global_dof_index>>
                                     send_dof_numbers_and_indices(ghost_owners.size());
      std::vector<MPI_Request>       reply_requests(ghost_owners.size());
      std::vector<std::vector<char>> sendbuffers(ghost_owners.size());

      for (unsigned int idx = 0; idx < ghost_owners.size(); ++idx)
        {
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               mpi_tag,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int len;
          ierr = MPI_Get_count(&status, MPI_BYTE, &len);
          AssertThrowMPI(ierr);
          Assert(len % sizeof(cell_data_to_send[idx][0]) == 0,
                 ExcInternalError());

          const unsigned int n_cells =
            len / sizeof(typename CellId::binary_type);
          cell_data_to_send[idx].resize(n_cells);

          ierr = MPI_Recv(cell_data_to_send[idx].data(),
                          len,
                          MPI_BYTE,
                          status.MPI_SOURCE,
                          status.MPI_TAG,
                          tria->get_communicator(),
                          &status);
          AssertThrowMPI(ierr);

          // store data for each cell
          for (unsigned int c = 0; c < static_cast<unsigned int>(n_cells); ++c)
            {
              const auto cell =
                tria->create_cell_iterator(CellId(cell_data_to_send[idx][c]));

              MeshCellIteratorType                mesh_it(tria,
                                           cell->level(),
                                           cell->index(),
                                           &mesh);
              const std_cxx17::optional<DataType> data = pack(mesh_it);

              if (data)
                {
                  typename DestinationToBufferMap::iterator p =
                    destination_to_data_buffer_map
                      .insert(std::make_pair(
                        idx,
                        GridTools::CellDataTransferBuffer<dim, DataType>()))
                      .first;

                  p->second.cell_ids.emplace_back(cell->id());
                  p->second.data.emplace_back(*data);
                }
            }

          // send reply
          GridTools::CellDataTransferBuffer<dim, DataType> &data =
            destination_to_data_buffer_map[idx];

          sendbuffers[idx] =
            Utilities::pack(data,  /*enable_compression*/  false);
          ierr = MPI_Isend(sendbuffers[idx].data(),
                           sendbuffers[idx].size(),
                           MPI_BYTE,
                           status.MPI_SOURCE,
                           mpi_tag_reply,
                           tria->get_communicator(),
                           &reply_requests[idx]);
          AssertThrowMPI(ierr);
        }

      // finally receive the replies
      std::vector<char> receive;
      for (unsigned int idx = 0; idx < ghost_owners.size(); ++idx)
        {
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               mpi_tag_reply,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int len;
          ierr = MPI_Get_count(&status, MPI_BYTE, &len);
          AssertThrowMPI(ierr);

          receive.resize(len);

          char *ptr = receive.data();
          ierr      = MPI_Recv(ptr,
                          len,
                          MPI_BYTE,
                          status.MPI_SOURCE,
                          status.MPI_TAG,
                          tria->get_communicator(),
                          &status);
          AssertThrowMPI(ierr);

          auto cellinfo =
            Utilities::unpack<CellDataTransferBuffer<dim, DataType>>(
              receive,  /*enable_compression*/  false);

          DataType *data = cellinfo.data.data();
          for (unsigned int c = 0; c < cellinfo.cell_ids.size(); ++c, ++data)
            {
              const typename Triangulation<dim, spacedim>::cell_iterator
                tria_cell = tria->create_cell_iterator(cellinfo.cell_ids[c]);

              MeshCellIteratorType cell(tria,
                                        tria_cell->level(),
                                        tria_cell->index(),
                                        &mesh);

              unpack(cell, *data);
            }
        }

      // make sure that all communication is finished
      // when we leave this function.
      if (requests.size() > 0)
        {
          const int ierr =
            MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
      if (reply_requests.size() > 0)
        {
          const int ierr = MPI_Waitall(reply_requests.size(),
                                       reply_requests.data(),
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }


#    endif // DEAL_II_WITH_MPI
    }

  } // namespace internal

  template <typename DataType, typename MeshType>
  inline void
  exchange_cell_data_to_ghosts(
    const MeshType &                                     mesh,
    const std::function<std_cxx17::optional<DataType>(
      const typename MeshType::active_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::active_cell_iterator &,
                             const DataType &)> &        unpack,
    const std::function<bool(const typename MeshType::active_cell_iterator &)>
      &cell_filter)
  {
#    ifndef DEAL_II_WITH_MPI
    (void)mesh;
    (void)pack;
    (void)unpack;
    (void)cell_filter;
    Assert(false, ExcNeedsMPI());
#    else
    internal::exchange_cell_data<DataType,
                                 MeshType,
                                 typename MeshType::active_cell_iterator>(
      mesh,
      pack,
      unpack,
      cell_filter,
      [&](const auto &process) {
        for (const auto &cell : mesh.active_cell_iterators())
          if (cell->is_ghost())
            process(cell, cell->subdomain_id());
      },
      [](const auto &tria) { return tria.ghost_owners(); });
#    endif
  }



  template <typename DataType, typename MeshType>
  inline void
  exchange_cell_data_to_level_ghosts(
    const MeshType &                                    mesh,
    const std::function<std_cxx17::optional<DataType>(
      const typename MeshType::level_cell_iterator &)> &pack,
    const std::function<void(const typename MeshType::level_cell_iterator &,
                             const DataType &)> &       unpack,
    const std::function<bool(const typename MeshType::level_cell_iterator &)>
      &cell_filter)
  {
#    ifndef DEAL_II_WITH_MPI
    (void)mesh;
    (void)pack;
    (void)unpack;
    (void)cell_filter;
    Assert(false, ExcNeedsMPI());
#    else
    internal::exchange_cell_data<DataType,
                                 MeshType,
                                 typename MeshType::level_cell_iterator>(
      mesh,
      pack,
      unpack,
      cell_filter,
      [&](const auto &process) {
        for (const auto &cell : mesh.cell_iterators())
          if (cell->level_subdomain_id() !=
                dealii::numbers::artificial_subdomain_id &&
              !cell->is_locally_owned_on_level())
            process(cell, cell->level_subdomain_id());
      },
      [](const auto &tria) { return tria.level_ghost_owners(); });
#    endif
  }
} // namespace GridTools

#  endif

DEAL_II_NAMESPACE_CLOSE

 /*----------------------------   grid_tools.h     ---------------------------*/ 
 /* end of #ifndef dealii_grid_tools_h */ 
#endif
 /*----------------------------   grid_tools.h     ---------------------------*/ 


