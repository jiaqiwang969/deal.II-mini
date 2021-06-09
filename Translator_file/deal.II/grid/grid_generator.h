//include/deal.II-translator/grid/grid_generator_0.txt
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

#ifndef dealii_grid_generator_h
#define dealii_grid_generator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria.h>

#include <array>
#include <map>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个命名空间提供了一个函数集合，用于为一些基本的几何体生成三角图。
 * 如果域是弯曲的，应该按照适当的Manifold描述进行细化的每个域部分将收到不同的 @ref GlossManifoldIndicator "流形指标"
 * ，并且正确的Manifold描述符将被附加到三角图中。请注意，如果你以后对三角测量进行转换，你必须确保将正确的新Manifold附加到三角测量中。
 *
 *
 * @ingroup grid
 *
 *
 */
namespace GridGenerator
{
  /**
   * @name  为基本几何形状创建网格
   *
   */
  ///@{

  /**
   * 用一个恰好由一个单元组成的超立方体（一维的线，二维的方，等等）初始化给定的三角结构。超立方体的体积是目前维数中的张量乘积区间 $[left,right]^{\text{dim}}$ ，其中极限作为参数给出。它们默认为零和一，然后产生单位超立方体。    如果参数 @p colorize 是假的，那么2d和3d的所有边界指标都被设置为零（默认边界指标）。如果它是真，那么边界就会像hyper_rectangle()中那样被 @ref GlossColorization "着色"
   * 。在1d中，指标总是被着色的，见hyper_rectangle()。
   * @image html hyper_cubes.png  如果 @p dim  <  @p spacedim,
   * 这将在第一个 @p dim 坐标方向创建一个 @p dim
   * 维的对象，嵌入到 @p spacedim
   * 维的空间，其余的条目设置为零。例如，一个<tt>Triangulation
   * @<2,3@></tt>
   * 将是xy平面上的一个正方形，z=0。对于由多个单元组成的粗略网格，也请参见subdivided_hyper_cube()。如果需要在不同的序数方向上有不同的长度，参见hyper_rectangle()。
   * @pre
   * 当调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim, int spacedim>
  void
  hyper_cube(Triangulation<dim, spacedim> &tria,
             const double                  left     = 0.,
             const double                  right    = 1.,
             const bool                    colorize = false);

  /**
   创建一个 $d$  -<a href="https://en.wikipedia.org/wiki/Simplex">simplex</a>（即2D中的三角形，或3D中的四面体），其角为 $d+1$ 。由于deal.II不支持三角形和四面体单元，输入参数所描述的单面体通过增加边缘、面和单面体中点被细分为四边形和六面体，从而得到一个由 $d+1$ 四边形或六面体单元组成的网格。     @p vertices  参数包含一个向量，其中包含所有定义单线角的d+1个顶点。它们必须以这样的顺序给出，即从第一个顶点到其他每个顶点的向量形成一个右手系统。    在二维和三维中生成的网格是。      @image html simplex_2d.png   @image html simplex_3d.png   @param  tria 要创建的三角结构。在调用此函数时，它需要是空的。      @param  顶点 单纯线的dim+1角。
   * @note  为<tt>三角法 @<2,2@></tt>,  <tt>三角法 @<3,3@></tt>.
   * 实现。
   *
   */
  template <int dim>
  void
  simplex(Triangulation<dim, dim> &      tria,
          const std::vector<Point<dim>> &vertices);

  /* 创建一个具有所提供的参考单元形状的单一单元的（粗略）网格。这是上面hyper_cube()和simplex()函数的一个概括。 
*
*/
  template <int dim, int spacedim>
  void
  reference_cell(Triangulation<dim, spacedim> &tria,
                 const ReferenceCell &         reference_cell);


  /**
   * 与hyper_cube()相同，但不同的是，不仅创建一个单元，而且每个坐标方向被细分为
   * @p repetitions
   * 个单元。因此，填充给定体积的单元的数量是<tt>repetitions<sup>dim</sup></tt>。
   * 如果 @p dim < @p spacedim, ，这将在第一个 @p dim
   * 坐标方向创建一个 @p dim 维度的对象，嵌入到 @p spacedim
   * 维度的空间，其余条目设置为零。例如， @<2,3@></tt>
   * 三角形将是xy平面上的一个正方形，z=0。  @pre
   * 在调用此函数时，作为参数传递的三角形需要为空。
   * @param  tria
   * 要创建的三角结构。调用此函数时，它必须为空。
   * @param  repetitions
   * 一个无符号整数，表示在每个方向上生成的单元数。
   * @param  left 用于创建超立方体的区间的下限。      @param
   * 右 用于创建超立方体的区间的上界。      @param  colorize
   * 如果设置为 "true"，则指定不同的边界ID。
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_cube(Triangulation<dim, spacedim> &tria,
                        const unsigned int            repetitions,
                        const double                  left     = 0.,
                        const double                  right    = 1.,
                        const bool                    colorize = false);

  /**
   *
   */
  template <int dim, int spacedim>
  void
  hyper_rectangle(Triangulation<dim, spacedim> &tria,
                  const Point<dim> &            p1,
                  const Point<dim> &            p2,
                  const bool                    colorize = false);

  /**
   * 从两个斜对角的角点 @p p1 和 @p p2. 创建一个坐标平行的砖块，坐标方向 @p i 的单元格数量由整数<tt>repetitions[i]<tt>给出。    为了得到与域的长宽比不同的单元，使用不同的细分数，由不同坐标方向的 @p repetitions, 给出。如果 @p colorize 的标志是 <code>true</code> ，则分配曲面的 @p boundary_ids ，这样 @p x-direction 中较低的是0，较高的是1（左边和右边的垂直面）。 @p y-direction 中表面的指标是2和3， @p z 的指标是4和5。 此外，材料ID是根据单元格的中心所在的八度空间来分配的：对于任何坐标方向<i>x<sub>i</sub></i>来说，在右半面都会增加2<sup>i</sup>（见 @ref GlossColorization "关于着色的词汇条"
   * ）。  例如，中心点（1,-1,1）产生一个材料id
   * 5（这意味着在2d中只有材料id
   * 0,1,2,3被分配，与重复的数量无关）。    请注意， @p
   * colorize
   * 标志在1d中被忽略，并被假定为始终为真。这意味着边界指标在左边是0，在右边是1。
   * 详见 step-15 。    如果 @p dim < @p spacedim, ，这将在第一个
   * @p dim 坐标方向创建一个 @p dim 维度的对象，嵌入到 @p
   * spacedim
   * 维度空间，其余条目设置为零。例如，一个<tt>三角法
   * @<2,3@></tt>
   * 将是xy平面上的一个矩形，z=0，由两个相对的角 @p p1 和
   * @p p2. 定义。
   * @note  关于这个函数的使用实例，见 step-28 教程程序。
   * @param  tria
   * 要创建的三角结构。在调用此函数时，它需要是空的。
   * @param  Repetitions  @p dim
   * 正值的向量，表示在该方向上生成的单元的数量。
   * @param  p1 第一个角点。      @param  p2 第二个角点与 @p p1.
   * 相对  @param  colorize
   * 如果设置为真，则指定不同的边界id。与hyper_rectangle()函数适用的注释相同。
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle(Triangulation<dim, spacedim> &   tria,
                             const std::vector<unsigned int> &repetitions,
                             const Point<dim> &               p1,
                             const Point<dim> &               p2,
                             const bool                       colorize = false);

  /**
   * 和前面的函数一样。然而，这里的第二个参数并不表示每个坐标方向上的细分数量，而是表示每个坐标方向上的步长序列。因此，域将在坐标方向
   * <code>i</code> 上被细分为 <code>step_sizes[i].size()</code>
   * 个单元，宽度为 <code>step_sizes[i][j]</code> for the
   * <code>j</code>  个单元。
   * 因此，这个函数适合于生成分级网格，单元集中在某些区域，而不是像前一个函数生成的均匀细分的网格。
   * 步长必须加到由点 @p p1 和 @p p2. 指定的超矩形的尺寸。
   *
   */
  template <int dim>
  void
  subdivided_hyper_rectangle(Triangulation<dim> &                    tria,
                             const std::vector<std::vector<double>> &step_sizes,
                             const Point<dim> &                      p_1,
                             const Point<dim> &                      p_2,
                             const bool colorize = false);

  /**
   * 就像之前的函数一样，但有以下变化： @p
   * 的material_id参数是一个dim-dimensional数组，对于每个单元，它指示应该设置哪个material_id。此外，这是主要的新功能，如果一个单元格的
   * material_id 是 <tt>(unsigned
   * char)(-1)</tt>，那么这个单元格就会从三角结构中删除，也就是说，域中会有一个空白。
   * @note 如果你需要大量的孔，你可以考虑cheese()。
   *
   */
  template <int dim>
  void
  subdivided_hyper_rectangle(Triangulation<dim> &                    tria,
                             const std::vector<std::vector<double>> &spacing,
                             const Point<dim> &                      p,
                             const Table<dim, types::material_id> &material_id,
                             const bool colorize = false);

  /**
   \brief 带矩形孔的矩形域 域本身是矩形的，非常像由subdivided_hyper_rectangle()生成的。参数 <code>holes</code> 指定了域在每个坐标方向上应该有多少个方孔。该方向上的网格单元总数是这个数字的两倍加一。    一个方向上的孔的数量必须至少是一个。    一个有二乘三孔的例子是  @image html cheese_2d.png  如果  @p dim  <  @p spacedim,  这将在第一个  @p dim  坐标方向创建一个  @p dim  维度的对象，嵌入到  @p spacedim  维度的空间，其余条目设置为零。      @param  tria 要创建的三角结构。在调用此函数时，它需要是空的。      @param  holes 在每个dim方向上的孔的正数。
   *
   */
  template <int dim, int spacedim>
  void
  cheese(Triangulation<dim, spacedim> &   tria,
         const std::vector<unsigned int> &holes);

  /**
   \brief 带有一个（偏移）圆柱形孔的矩形板。    生成一个带有（偏移）圆柱形孔的矩形板。该几何体由2个区域组成。  第一个是一个正方形区域，长度为 @p outer_radius ，半径为 @p inner_radius 的孔。  这个区域的单元格将有流形id为 @p tfi_manifold_id 的TransfiniteInterpolationManifold连接到它们。此外，洞的边界面将与PolarManifold（二维）或CylindricalManifold（三维）相关。这个区域的中心可以通过 @p center 规定，即孔的轴线将位于 @p center 处。  第二个区域描述了散装材料的剩余部分。它通过填充参数 @p pad_bottom,   @p padding_top,   @p padding_left 和 @p padding_right. 指定，这个区域的所有单元将有一个FlatManifold连接到它们。  板块的最终宽度将是<code>padding_left + 2*outer_radius + padding_right</code>，而其长度为<code>padding_top + 2*outer_radius + padding_bottom</code>。    下面是非对称网格（经过一次全局细化，根据流形ID着色），分别是2D和3D的。    \htmlonly <style>div.image img[src="plate_with_a_hole.png"]{width:25%;}</style> \endhtmlonly  @image html plate_with_a_hole.png  \htmlonly <style>div.image img[src="plate_with_a_hole_3D.png"]{width:25%;}</style> \endhtmlonly  @image html plate_with_a_hole_3D.png  ] 在3D中，三角形将在Z方向上被挤压，总高度为 @p L 使用 @p n_slices 片（最小为2）。
   * 如果 @p colorize 标志是 <code>true</code> ，则边界面的边界_id被分配为，在x方向上的低位是0，高位是1（见 @ref GlossColorization "关于着色的词汇条"
   * ）。
   * y方向的面的指标是2和3，z方向的是5和6。孔洞边界的指标是4。
   * @p tria
   * 是要创建的三角形。在调用这个函数时，它需要是空的。
   *
   */
  template <int dim>
  void
  plate_with_a_hole(Triangulation<dim> &     tria,
                    const double             inner_radius      = 0.4,
                    const double             outer_radius      = 1.,
                    const double             pad_bottom        = 2.,
                    const double             pad_top           = 2.,
                    const double             pad_left          = 1.,
                    const double             pad_right         = 1.,
                    const Point<dim> &       center            = Point<dim>(),
                    const types::manifold_id polar_manifold_id = 0,
                    const types::manifold_id tfi_manifold_id   = 1,
                    const double             L                 = 1.,
                    const unsigned int       n_slices          = 2,
                    const bool               colorize          = false);

  /**
   * 生成一个由通道与圆柱体组成的网格。这是纳维-斯托克斯求解器的一个常见基准。该几何体包括一个尺寸为 $[0, 2.2] \times [0, 0.41] \times [0, 0.41] $ 的通道（其中 $z$ 尺寸在二维中被省略），该通道有一个圆柱体，平行于 $z$ 轴，直径为 $0.1$ ，中心为 $(0.2, 0.2, 0)$ 。该通道有三个不同的区域。    <ol>   <li>  如果 @p n_shells 大于零，那么就有那么多以圆柱体为中心的壳， </li>   <li> 是壳和三角形其他部分之间的混合区域，以及 </li>   <li> 是由笛卡尔单元组成的散装区域。 </li>   </ol>  由于圆柱体略微偏离通道中心，这种几何形状导致了中等雷诺数下的涡流脱落。下面是二维的网格（没有额外的全局细化）： @image html channel_with_cylinder_2d.png 和三维的网格： @image html channel_with_cylinder_3d.png 由此产生的三角形使用了三个流形：一个PolarManifold（二维）或CylindricalManifold（三维），流形ID为 $0$ ，一个TransfiniteInterpolationManifold，流形ID为 $1$ ，其他地方为FlatManifold。关于这一主题的更多信息，请参见 @ref GlossManifoldIndicator "流形指标的词汇表条目"。  圆柱体和周围贝壳上的单元面的流形指标为 $0$ ，而与贝壳相邻的单元体（如果不存在，则为圆柱体）的流形指标为 $1$  。换句话说：这个网格使用TransfiniteInterpolationManifold来平滑地从壳（用 GridGenerator::concentric_hyper_shells) 生成）过渡到体块区域。所有其他单元体和面的流形ID为 numbers::flat_manifold_id 并使用FlatManifold。所有id为 numbers::flat_manifold_id 的单元都是与坐标轴对齐的矩形棱镜。    下图显示了两次全局细化后的部分二维网格（使用该函数的所有默认参数）。流形标识 $0$ 的单元为橙色（极坐标流形标识），流形标识 $1$ 的单元为黄色（无限插值流形标识），流形标识 numbers::flat_manifold_id 的单元为青色：  @image html channel_with_cylinder_2d_manifolds.png   @param  tria 要创建的三角形。调用此函数时必须为空。      @param  shell_region_width 围绕圆柱体的壳层的宽度。  这个值应该在  $0$  和  $0.05$  之间；默认值是  $0.03$  。      @param  n_shells 在壳层中使用的壳的数量。      @param  skewness 控制壳与圆柱体的接近程度的参数：见 GridGenerator::concentric_hyper_shells. 中给出的数学定义  @param  colorize 如果设置为true，则分配不同的边界ID。关于边界指示器的更多信息见 @ref GlossBoundaryIndicator "本词汇表条目"。  左边的边界（在  $x = 0$  处）被分配一个  $0$  的id，右边的边界（在  $x = 2.2$  处）被分配一个  $1$  的id，圆柱体的边界被分配一个  $2$  的id，而通道壁被分配一个  $3$  的id 。    更多信息请参见原始论文。
   * @code{.bib}
   * @inbook{schafer1996,
   * author    = {Sch{\"a}fer, M. and Turek, S. and Durst, F. and Krause, E.
   *            and Rannacher, R.},
   * title     = {Benchmark Computations of Laminar Flow Around a Cylinder},
   * bookTitle = {Flow Simulation with High-Performance Computers II: DFG
   *            Priority Research Programme Results 1993--1995},
   * year      = {1996},
   * publisher = {Vieweg+Teubner Verlag},
   * address   = {Wiesbaden},
   * pages     = {547--566},
   * isbn      = {978-3-322-89849-4},
   * doi       = {10.1007/978-3-322-89849-4_39},
   * url       = {https://doi.org/10.1007/978-3-322-89849-4_39}
   * }
   * @endcode
   *
   */
  template <int dim>
  void
  channel_with_cylinder(Triangulation<dim> &tria,
                        const double        shell_region_width = 0.03,
                        const unsigned int  n_shells           = 2,
                        const double        skewness           = 2.0,
                        const bool          colorize           = false);

  /**
   * 一个一般的  @p dim
   *
   * - 浸入 @p dim 的一般单元（如果dim为1，则为线段，如果 @p dim 为2，则为四边形，如果 @p dim 为3，则为六面体）。
   *
   * - 膨胀空间。用户有责任按照正确的顺序提供顶点（参见GeometryInfo类的文档），因为顶点是按照给出的顺序存储的。确保单元格的体积为正值也很重要。    如果参数 @p colorize 是假的，那么2d和3d的所有边界指标都被设置为零。如果是true，边界就会像hyper_rectangle()中那样被着色（见 @ref GlossColorization "关于着色的词汇表条目"）。  在1d中，指标总是被着色的，见hyper_rectangle()。      @param  tria 将被创建的三角形  @param  vertices 单元的2^dim顶点  @param  colorize 如果为真，设置不同的边界id。
   *
   */
  template <int dim, int spacedim>
  void
  general_cell(Triangulation<dim, spacedim> &      tria,
               const std::vector<Point<spacedim>> &vertices,
               const bool                          colorize = false);

  /**
   * 一个平行四边形。第一个角点是原点。接下来的 @p dim
   * 顶点是第二个参数中给出的顶点，最后一个顶点将是连接原点和这些点的两个向量之和。着色的方式与hyper_rectangle()中的方式相同。
   * @note  这个函数只在2d中实现。      @param  tria 要创建的三角结构。在调用此函数时，它需要是空的。      @param  corners 平行四边形的第二个和第三个顶点。      @param  colorize 如果为真，则指定不同的边界ID。参见  @ref GlossColorization  "关于着色的词汇表条目"
   * ）。
   *
   */
  template <int dim>
  void
  parallelogram(Triangulation<dim> &tria,
                const Point<dim> (&corners)[dim],
                const bool colorize = false);

  /**
   * 一个平行四边形。第一个角点是原点。 @p dim
   * 相邻的点是描述平行四边形相对于原点的边缘的向量。额外的点是这些凹陷向量的总和。着色是根据hyper_rectangle()进行的。
   * @note
   * 这个函数默默地将单元格上的顶点重新排序为lexicographic排序（见
   * <code>GridReordering::reorder_grid</code>  ）。
   * 换句话说，如果顶点的重新排序确实发生了，
   * <code>corners</code>
   * 的数组中的顶点排序将不再是指同一个三角形。      @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  parallelepiped(Triangulation<dim> &tria,
                 const Point<dim> (&corners)[dim],
                 const bool colorize = false);

  /**
   * 一个细分的平行四边形。第一个角点是原点。 @p
   * 相邻的dim点是描述平行四边形相对于原点的边缘的向量。其他的点是这些dim向量的总和。变量
   * @p n_subdivisions 指定了每个 @p dim
   * 方向上的细分数量。着色是根据hyper_rectangle()来完成的。
   * @pre  调用此函数时，作为参数传递的三角图需要为空。
   *
   */
  template <int dim>
  void
  subdivided_parallelepiped(Triangulation<dim> &tria,
                            const unsigned int  n_subdivisions,
                            const Point<dim> (&corners)[dim],
                            const bool colorize = false);

  /**
   * 一个细分的平行四边形，即与上述相同，但每个 @p dim
   * 方向上的细分数量可能不同。
   * 着色是根据hyper_rectangle()进行的。      @pre
   * 调用此函数时，作为参数传递的三角图需要为空。
   *
   */
  template <int dim>
  void
  subdivided_parallelepiped(Triangulation<dim> &tria,
#ifndef _MSC_VER
                            const unsigned int (&n_subdivisions)[dim],
#else
                            const unsigned int *n_subdivisions,
#endif
                            const Point<dim> (&corners)[dim],
                            const bool colorize = false);

  /**
   * @note  对 @p dim 和 @p spacedim. 的所有组合实施。
   * @note
   * 你可能需要帮助编译器，在调用此函数时明确指定两个模板参数。
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_parallelepiped(Triangulation<dim, spacedim> &              tria,
                            const Point<spacedim> &                     origin,
                            const std::array<Tensor<1, spacedim>, dim> &edges,
                            const std::vector<unsigned int> &subdivisions = {},
                            const bool                       colorize = false);

  /**
   * 超立方体，周围有一层超立方体。参数 @p left 和 @p right 给出了所有坐标方向上的内超立方体的下限和上限。   @p thickness 标记了层单元的大小。    如果标志 @p colorize 被设置，外部单元根据以下方案获得材料ID：在（+/-）x方向1/2、y方向4/8、z方向16/32的内立方体上延伸。在角落和边缘（3D），使用一个比特OR操作来获得这些值，（也见 @ref GlossColorization "关于着色的术语条目"
   * ）。    目前只有2d和3d版本可用。      @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  enclosed_hyper_cube(Triangulation<dim> &tria,
                      const double        left      = 0.,
                      const double        right     = 1.,
                      const double        thickness = 1.,
                      const bool          colorize  = false);

  /**
   * @note
   * 由于这可能是用户通常最早考虑的创建具有曲线边界的网格的函数之一，让我们也来评论一下经常令人困惑的一个方面。也就是说，人们看到的并不总是实际发生的情况。具体来说，如果你使用默认选项，用
   * GridOut::write_vtk()
   * 这样的函数输出粗略的网格，那么人们一般不会看到边界上的曲面。
   * 这是因为大多数文件格式默认只存储顶点位置，隐含的理解是单元由这些顶点组成，并以直边为界。同时，这个函数将SphericalManifold对象附加到边界面的事实意味着，至少在内部*，边缘确实是弯曲的。如果你想看到它们，你需要确保你用来输出网格的函数实际上是将边界面绘制成曲线，而不是仅由两个端点的位置来描述的直线。例如，如果你在
   * GridOutFlags::Gnuplot 结构中设置相应的标志，
   * GridOut::write_gnuplot()
   * 就可以做到这一点。然而，你是否真的在曲线单元上进行计算*，这是一个完全独立的考虑。在典型的有限元计算中，我们必须计算积分，这些积分是通过使用参考单元的映射将实际单元转换回来计算的。使用什么样的映射决定了这些内部计算的单元的形状。例如，使用广泛使用的
   * $Q_1$ 映射（隐含在 step-6
   * 中使用），积分总是发生在假定只有顶点位置描述的直线边界的单元上。换句话说，如果使用这样的映射，那么域的单元就真的有直的边缘，而不管这些边缘的流形描述如何，也不管生成输出时的标志是什么。综上所述，有必要区分三件事。(i)附加在网格中物体上的流形描述；(ii)用于集成的映射；以及(iii)用于输出网格图形信息的风格。所有这些都可以或多或少地独立选择，你所看到的可视化的东西不一定就是正在发生的。
   * @pre  调用此函数时，作为参数传递的三角形需要为空。
   *
   */
  template <int dim>
  void
  hyper_ball(Triangulation<dim> &tria,
             const Point<dim> &  center = Point<dim>(),
             const double        radius = 1.,
             const bool attach_spherical_manifold_on_boundary_cells = false);

  /**
   * 这是一个替代hyper_ball的方法，在2D中使用12个单元，在3D中使用32个单元，这样可以更好地平衡外部弯曲边界周围的单元和内部的单元的大小。网格是基于
   * GridGenerator::quarter_hyper_ball()
   * 所使用的单元，并进行适当的复制和旋转以填充整个球。
   * 下面的图片显示了二维（左）和三维的网格结果：
   * <table align="center" class="doxtable"> <tr> <td> \htmlonly
   * <style>div.image img[src="hyper_ball_balanced_2d.png"]{width:40%}</style>
   * \endhtmlonly
         @image html hyper_ball_balanced_2d.png
   * </td> <td> \htmlonly <style>div.image
   * img[src="hyper_ball_balanced_3d.png"]{width:40%}</style> \endhtmlonly
         @image html hyper_ball_balanced_3d.png
   * </td> </tr> </table>
   * 默认情况下，manifold_id在边界面设置为0，在边界单元设置为1，
   * numbers::flat_manifold_id  在中心单元和内部面设置为1。
   * @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  hyper_ball_balanced(Triangulation<dim> &tria,
                      const Point<dim> &  center = Point<dim>(),
                      const double        radius = 1.);

  /**
   * 生成一个二维网格，该网格由单元格和一个移位  $s =
   * (1,0)$
   * 的副本组成。根据所传递的标志，右边或左边的正方形将被旋转
   * $\pi/2$
   * 。这样就可以产生一个网格，其中一个正方形可能包含一个与另一个正方形的相邻边的切向（因此也是相反的法向）的边。
   * 从实用的角度来看，这种网格并没有太大用处。出于调试的目的，它可以用来检查矢量或张量有限元的方向问题。
   * @note  如果 <code>rotate_left_square==rotate_right_square</code>
   * 网格的方向一致。      @param[out]  tria 输入的三角剖面。
   * @param[in]  rotate_left_square  <code>true</code>
   * 如果左边的正方形被旋转  $\pi/2$  .     @param[in]
   * rotate_right_square  <code>true</code>  如果右边的正方形被
   * $\pi/2$  旋转。
   *
   */
  void non_standard_orientation_mesh(Triangulation<2> &tria,
                                     const bool        rotate_left_square,
                                     const bool        rotate_right_square);

  /**
   * 生成一个由单位立方体和一个移位  $s = (1,0,0)$
   * 的副本组成的三维网格。根据所传递的标志，右边或左边的立方体（当看到正向的(x,z)平面时）包含一个不是标准方向的面和/或被
   * $\pi/2$  ,  $\pi$  或  $3/2\pi$  旋转过的面。
   * 从实用的角度来看，这个网格没有太大用处。出于调试的目的，它可以用来检查矢量或张量的有限元的方向问题。
   * @param[out]  tria 输入的三角形网格。    @param[in]
   * face_orientation  <code>true</code>  如果该面不是标准方向。
   * @param[in]  face_flip  <code>true</code>  如果面被旋转+180度
   * @param[in]  face_rotation  <code>true</code>
   * 如果面被旋转（另外）+90度  @param[in]  manipulate_left_cube
   * <code>true</code>
   * 如果左侧立方体要被重新排序。如果`false`，则是右方的立方体。
   *
   */
  void non_standard_orientation_mesh(Triangulation<3> &tria,
                                     const bool        face_orientation,
                                     const bool        face_flip,
                                     const bool        face_rotation,
                                     const bool        manipulate_left_cube);


  /**
   * 创建一个超球体，即在 @p spacedim
   * 维度上的球的表面。这个函数只存在于dim+1=spacedim的2和3空间维度。(要创建一个球的网格，请使用
   * GridGenerator::hyper_ball().)
   * 默认情况下，三角形的所有流形id被设置为零，并且一个SphericalManifold被附加到网格上。
   * 下面的图片是用以下方法生成的。
   * @code
   * Triangulation<2,3>   triangulation;
   * GridGenerator::hyper_sphere(triangulation);
   * triangulation.refine_global(3);
   * @endcode
   * 参见 @ref manifold "流形的文件模块"
   * ，以了解更多细节。      @image html sphere.png   @image html
   * sphere_section.png   @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int spacedim>
  void hyper_sphere(Triangulation<spacedim - 1, spacedim> &tria,
                    const Point<spacedim> &center = Point<spacedim>(),
                    const double           radius = 1.);

  /**
   * 这个函数产生一个与相对于 @p center,
   * 的正交的超球体，它包含了2d的三个元素和3d的四个元素。网格的内部点的选择是为了平衡内部点周围的单元中从参考坐标到实坐标的映射的最小单值，这相当于一个高的网格质量。
   * 最终三角化的边界指标是：曲线边界为0，切割面为1。弯曲边界的流形ID被设置为0，并且一个SphericalManifold被附加到它。
   * 由此产生的二维和三维网格看起来如下。    <table
   * align="center" class="doxtable"> <tr> <td> \htmlonly <style>div.image
   * img[src="quarter_hyper_ball_2d.png"]{width:50%}</style> \endhtmlonly
         @image html quarter_hyper_ball_2d.png
   * </td> <td> \htmlonly <style>div.image
   * img[src="quarter_hyper_ball_3d.png"]{width:46%}</style> \endhtmlonly
         @image html quarter_hyper_ball_3d.png
   * </td> </tr> </table>   @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  quarter_hyper_ball(Triangulation<dim> &tria,
                     const Point<dim> &  center = Point<dim>(),
                     const double        radius = 1.);

  /**
   * 这个函数在 @p center,
   * 周围产生一个半超球，它在2d中包含4个元素，在3d中包含6个。切割面垂直于<i>x</i>轴。
   * 最终三角化的边界指标为：弯曲边界为0，切割面为1。弯曲边界的流形ID被设置为0，并且一个SphericalManifold被连接到它。
   * @pre
   * 调用此函数时，作为参数传递的三角剖面需要为空。
   *
   */
  template <int dim>
  void
  half_hyper_ball(Triangulation<dim> &tria,
                  const Point<dim> &  center = Point<dim>(),
                  const double        radius = 1.);

  /**
   * 创建一个 @p dim 维圆柱体，其中 $x$
   * -轴作为圆柱体的轴。在这个函数中，圆柱体被定义为一个（
   * @p dim  ）。
   *
   * -1）维的圆盘，其给定的 @p radius, 沿着圆柱体的轴线（即第一坐标方向）挤压。因此，在三维空间中，圆柱体从`x=-半长`延伸到`x=+半长`，它在 @p yz-plane 的投影是一个半径为 @p 的圆。在二维空间中，圆柱体是一个从`x=-半长`到`x=+半长`，从`y=-半径`到`y=半径`的矩形。    边界按照以下方案着色：0代表圆柱体的外壳，1代表左手面，2代表右手面（见 @ref GlossColorization "关于着色的术语条目"
   * ）。
   * 圆柱体外壳的流形ID被设置为0，并将一个CylindricalManifold连接到它。
   * @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  cylinder(Triangulation<dim> &tria,
           const double        radius      = 1.,
           const double        half_length = 1.);


  /**
   * 创建一个 @p dim 维圆柱体，其中 $x$
   * -轴作为圆柱体的轴。在本函数中，圆柱体被定义为一个（
   * @p dim  ）。
   *
   * -1）的维度盘，其给定的 @p radius, 沿着圆柱体的轴线（即第一坐标方向）挤压。因此，在三维空间中，圆柱体从`x=-半长`延伸到`x=+半长`，它在 @p yz-plane 的投影是一个半径为 @p 的圆。在二维空间中，圆柱体是一个从`x=-半长`到`x=+半长`，从`y=-半径`到`y=半径`的矩形。这个函数只在dim==3时实现。    边界按照以下方案着色：0代表圆柱体的外壳，1代表左手面，2代表右手面（见 @ref GlossColorization "关于着色的词汇表条目"
   * ）。
   * 圆柱体的流形ID被设置为0，并将一个CylindricalManifold连接到它。
   * @image html subdivided_cylinder_3D.png   @param  tria
   * 要创建的三角剖面。在调用此函数时，它需要是空的。
   * @param  x_subdivisions
   * 一个正整数，表示在x方向上生成的单元的数量。默认圆柱体的x_repetitions=2。
   * @param  radius 用于挤出圆柱体的YZ平面上的圆的半径。
   * @param  half_length 圆柱体在x方向的半长。
   *
   */
  template <int dim>
  void
  subdivided_cylinder(Triangulation<dim> &tria,
                      const unsigned int  x_subdivisions,
                      const double        radius      = 1.,
                      const double        half_length = 1.);


  /**
   * 围绕x轴创建一个切割锥体。 该圆锥体从<tt>x=-半长</tt>延伸到<tt>x=半长</tt>，其在 @p yz-plane 的投影是在<tt>x=半长</tt>处半径为 @p radius_0 的圆和在<tt>x=+半长</tt>处半径为 @p radius_1 的圆。在这两者之间，半径是线性递减的。    在二维空间中。圆锥体是一个梯形，从<tt>x=-半长</tt>到<tt>x=+半长</tt>，从<tt>y=-半径_0</tt>到<tt>y=半径_0</tt>在<tt>x=-半长</tt>，从<tt>y=-半径_1</tt>到<tt>y=半径_1</tt>在<tt>x=+半长>。  在这之间，<tt>y</tt>的范围是线性递减的。    边界按照以下方案着色：0代表圆锥体，1代表左手面，2代表右手面（见 @ref GlossColorization "关于着色的术语条目"
   * ）。  边界指标和流形指标都被设定。
   * 在三维空间中，船体的流形指标被设置为零，并且一个CylindricalManifold被附加到它上面。
   * 下面是两次网格细化后的二维和三维网格。      @image
   * html truncated_cone_2d.png   @image html truncated_cone_3d.png   @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  truncated_cone(Triangulation<dim> &tria,
                 const double        radius_0    = 1.0,
                 const double        radius_1    = 0.5,
                 const double        half_length = 1.0);

  /**
   * \brief 一个中心单元，每个表面都有堆积的单元突出。    每个正方形网格单元都是笛卡尔的，在每个坐标方向都有一个大小。零号单元的中心是原点。      @param  tria 一个三角形对象，必须是空的。      @param  sizes 一个维度为 GeometryInfo<dim>::faces_per_cell 的整数向量，其含义如下：十字形的腿在中心单元格的面上堆叠，按照deal.II单元格的通常顺序，即首先是 $-x$ ，然后是 $x$ ，然后是 $-y$ 等等。在 <code>sizes</code> 中的相应条目命名了堆积在这个面上的单元格的数量。所有的数字都可能是零，因此L型和T型域是这个域的特化。      @param  colorize_cells 如果着色被启用，那么一个单元格的材料ID就对应于它所在的腿。中心单元的id为0，然后各腿从1开始编号（参见 @ref GlossColorization "关于着色的术语条目"
   * ）。    二维和三维的例子是。      @image html
   * hyper_cross_2d.png   @image html hyper_cross_3d.png 。
   *
   */
  template <int dim, int spacedim>
  void
  hyper_cross(Triangulation<dim, spacedim> &   tria,
              const std::vector<unsigned int> &sizes,
              const bool                       colorize_cells = false);

  /**
   * 用一个恰好由<tt>2^dim-1</tt>单元组成的超L（2D或3D）初始化给定的三角结构。它产生的超立方体的区间[<i>left,right</i>]没有每个坐标的区间[<i>(left+right)/2,right</i>]做出来的超立方体。因为该域是大约最简单的一个有重心（即非凸）角的域，许多偏微分方程的解在这个角有奇点。也就是说，在拐角处，梯度或高阶导数（取决于所选择的边界条件）并不保持有界。因此，当溶液缺乏规律性时，这个域经常被用来测试方案的收敛性。    如果 @p colorize 的标志是 <code>true</code> ，则曲面的 @p boundary_ids 被分配为左边的边界为0，其他的则按逆时针升序分配（见 @ref GlossColorization "关于着色的词汇条"
   * ）。 @p 着色选项只在二维空间工作。
   * 这个函数将在二维中创建经典的L形，在三维中看起来就像下面这样。
   * @image html hyper_l.png
   * @note  3d域也经常被称为 "Fichera角"，这是以Gaetano
   * Fichera（1922-1996）命名的，他首次计算了域的最低特征函数的角奇异指数的近似值。
   * 这个函数存在于所有空间维度的三角形中，但如果在1d中调用，则会产生错误。
   * @pre  调用此函数时，作为参数传递的三角形需要为空。
   *
   */
  template <int dim>
  void
  hyper_L(Triangulation<dim> &tria,
          const double        left     = -1.,
          const double        right    = 1.,
          const bool          colorize = false);

  /**
   在二维或三维中用广义的细分超L来初始化给定的三角结构。    这个函数产生一个细分的超矩形，其尺寸由 @p bottom_left 和 @p top_right, 给出，每个方向上的细分数量由向量 @p repetitions, 给出，并去除一定数量的单元格，由向量 @p n_cells_to_remove. 给出。注意， @p n_cells_to_remove 包含整数，意味着其条目可以是正数和负数。正数表示在 "正 "的方向切割细胞，例如在X方向从左到右，在Y方向从下到上，在Z方向从前到后。负数表示以相反的方向切割单元，如从右到左，从上到下，从后到前。    这个网格的演示可以在  step-75  中找到。    这个函数可以用来生成一个面向后方的网格，这是一个对流体动力学基准问题有用的领域。  第一张图片是一个3D的后向台阶，通过去除z方向的所有单元，以及正x和y方向的2个单元而生成。    @image html subdivided_hyper_L_3d.png  而在二维中，我们可以剪去负x方向的1个单元，负y方向的2个单元。    @image html subdivided_hyper_L_2d.png
   * @note
   * 这个函数被声明存在于所有空间维度的三角形中，但如果在一维中调用则会产生错误。
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_L(Triangulation<dim, spacedim> &   tria,
                     const std::vector<unsigned int> &repetitions,
                     const Point<dim> &               bottom_left,
                     const Point<dim> &               top_right,
                     const std::vector<int> &         n_cells_to_remove);

  /**
   * 用一个带狭缝的超立方体初始化给定的三角测量。在每个坐标方向上，超立方体从 @p left 延伸到 @p right. 。 在2D中，分割在垂直方向上从<tt>x=(left+right)/2, y=left</tt>到广场中心<tt>x=y=(left+right)/2</tt>。    在3D中，2D域只是在<i>z</i>方向延伸，这样一个平面将矩形的下半部分一分为二。 这个函数被声明存在于所有空间维度的三角形中，但如果在1d中调用，会抛出一个错误。    如果 @p colorize 被设置为 @p true, ，形成狭缝的面将分别被标记为边界id 1和2（见 @ref GlossColorization "关于着色的词汇条"
   * ）。      @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  hyper_cube_slit(Triangulation<dim> &tria,
                  const double        left     = 0.,
                  const double        right    = 1.,
                  const bool          colorize = false);

  /**
   * 产生一个超壳，即围绕<tt>中心</tt>的两个球体之间的区域，给定<tt>内半径</tt>和<tt>外半径</tt>。数字<tt>n_cells</tt>表示所产生的三角形的单元数，即有多少单元构成环（在2D）或壳（在3D）。    如果标志 @p colorize 是 <code>true</code> ，那么外边界将有指标1，而内边界的id为0。在三维中，这同时适用于这些边界的面和边。如果标志是 @p false, ，则两者的指标都是0（见 @ref GlossColorization "关于着色的词汇条"）。    所有流形的id都被设置为零，并且一个SphericalManifold被附加到三角形的每个单元和面。    在2d中，这个初始三角形的元素数量<tt>n_cells</tt>可以任意选择。如果初始单元的数量为零（默认情况下），那么它将被自适应地计算，从而使生成的元素具有最小的长宽比。    在3D中，只有某些数字是允许的。    <ul>   <li>  6（或默认的0）用于基于六面体的曲面（即在内球面上的6个面板沿径向挤压形成6个单元）， <li>  12用于菱形十二面体， <li>  ] 24个为基于六面体的表面，在方位角方向上精炼一次，但不在径向方向上精炼， <li>  48个为菱形十二面体，在方位角方向上精炼一次，但不在径向方向上精炼， <li>  96个为菱形十二面体精炼一次。这个选择可以追溯到Manifold类实现之前的旧版本deal.II：今天这个选择等同于执行一次全局细化后的菱形十二面体。    <li>   $192\times 2^m$ 类的数字与 $m\geq 0$ 的整数。这种选择与24和48单元的情况类似，但提供了方位角方向的额外细化与径向方向的单层细化相结合。基本网格是6个或12个单元的版本，分别取决于权力中的 $m$ 是奇数还是偶数。    </ul> 24、48和 $2^m 192$ 单元的版本在壳薄而径向长度应与周向长度更相似的情况下是有用的。    下面是12和96单元的三维网格图。      @image html hypershell3d-12.png   @image html hypershell3d-96.png 。
   * @note
   * 这个函数被声明存在于所有空间维度的三角形中，但如果在1d中调用则会抛出一个错误。
   * @pre  调用此函数时，作为参数传递的三角形需要为空。
   *
   */
  template <int dim>
  void
  hyper_shell(Triangulation<dim> &tria,
              const Point<dim> &  center,
              const double        inner_radius,
              const double        outer_radius,
              const unsigned int  n_cells  = 0,
              bool                colorize = false);

  /**
   产生一个偏心的超壳，即以两个不同中心点为中心的两个球体之间的区域。我们必须指定<tt>inner_center</tt>和<tt>outer_center</tt>，并给出<tt>inner_radius</tt>和<tt>outer_radius</tt>。数字<tt>n_cells</tt>表示所产生的三角形的单元数，也就是说，有多少单元构成环（在2D）或壳（在3D）。    默认情况下，外边界的指标为1，而内边界的指标为0。在三维中，这适用于这些边界的面和边。    一个SphericalManifold连接到外边界，ID为1，而另一个SphericalManifold连接到内边界，ID为0。一个TransfiniteInterpolationManifold连接到三角形的所有其他单元和面，ID为2。这里，元素的数量<tt>n_cells</tt>与 GridGenerator::hyper_shell. 中的含义相同。      @image html eccentric_hyper_shell_2D.png   @image html eccentric_hyper_shell_3D.png 。
   * @note
   * 因为它使用了超壳的定义，这个函数被声明为存在于所有空间维度的三角形中，但如果在1d中调用，会抛出一个错误。
   * @pre  调用此函数时，作为参数传递的三角形需要为空。
   *
   */
  template <int dim>
  void
  eccentric_hyper_shell(Triangulation<dim> &triangulation,
                        const Point<dim> &  inner_center,
                        const Point<dim> &  outer_center,
                        const double        inner_radius,
                        const double        outer_radius,
                        const unsigned int  n_cells);

  /**
   * 产生一个半超壳，即两个空间维度中两个圆之间的空间和三维中两个球体之间的区域，对于这个初始三角剖面，具有给定的内外半径和给定的元素数。 然而，与前一个函数相反，它并不产生整个壳，而只是产生它的一半，即第一分量被限制为非负值的那部分。这个函数的目的是为了能够计算具有旋转对称性的解决方案，在这种情况下，2D的半壳代表3D的壳。    如果2d中初始单元 @p n_cells 的数量为零（默认值），那么它将被自适应地计算，从而使得到的元素具有最小的长宽比。这个参数在3D中被忽略，在3D中粗略的网格总是有5个单元。    如果colorize设置为 <code>true</code> ，内部、外部和边界 $x=0$ 的部分分别得到指标0、1和2。此外，在2d中，边界指标3被赋予X轴以下的垂直边缘。否则，如果colorize设置为 <code>false</code> ，所有指标都被设置为0（见 @ref GlossColorization "关于着色的词汇表条目"
   * ）。
   * 所有流形ID都被设置为0，并将SphericalManifold附加到三角剖面上。
   * @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  half_hyper_shell(Triangulation<dim> &tria,
                   const Point<dim> &  center,
                   const double        inner_radius,
                   const double        outer_radius,
                   const unsigned int  n_cells  = 0,
                   const bool          colorize = false);


  /**
   * 产生一个域，该域是具有给定内半径和外半径的超壳之间的交点，即两个空间维度中两个圆之间的空间和三维中两个球体之间的区域，以及正象限（在2d中）或八角形（在三维）。在2D中，这确实是全环形的四分之一，而在3D中，这个函数是一个错误的名称，因为在那里，域不是四分之一，而是全壳的八分之一。    如果初始单元的数量为零（默认情况下），那么它是自适应计算的，这样得到的元素具有2d中最小的长宽比。    如果 @p colorize 被设置为 <code>true</code> ，内、外、左、右边界分别得到2d中的指标0、1、2和3。在3D中，指标2位于面 $x=0$ ，3位于 $y=0$ ，4位于 $z=0$ （见 @ref GlossColorization "关于着色的词汇条"
   * ）。
   * 所有流形的id都被设置为零，并将一个SphericalManifold附加到三角结构上。
   * @pre
   * 调用此函数时，作为参数传递的三角结构需要为空。
   *
   */
  template <int dim>
  void
  quarter_hyper_shell(Triangulation<dim> &tria,
                      const Point<dim> &  center,
                      const double        inner_radius,
                      const double        outer_radius,
                      const unsigned int  n_cells  = 0,
                      const bool          colorize = false);

  /**
   * 产生一个域，该域是3D中两个圆柱体之间的空间，具有给定的长度、内半径和外半径以及给定的元素数量。圆柱体外壳围绕
   * $z$ 轴建立，两个面位于 $z = 0$ 和 $z = $   @p length.  如果
   * @p n_radial_cells
   * 为零（默认值），那么它将被自适应地计算，使生成的元素具有最小的纵横比。对于
   * @p n_axial_cells. 也是如此。
*  @note  虽然这个函数被声明为模板，但它在一维和二维中没有意义。同时请记住，这个对象的旋转和定位与圆柱体（）创建的对象不同。    所有流形的id都被设置为零，并且一个CylindricalManifold被附加到三角形中。      @pre  调用此函数时，作为参数传递的三角剖面需要为空。      @image html cylinder_shell.png  在这张图片中，显示了一个长度为2，内半径为0.5，外半径为1的圆柱体外壳。使用了n_radial_cells和n_axial_cells的默认参数，进行了一次全局细化。
   *
   */
  template <int dim>
  void
  cylinder_shell(Triangulation<dim> &tria,
                 const double        length,
                 const double        inner_radius,
                 const double        outer_radius,
                 const unsigned int  n_radial_cells = 0,
                 const unsigned int  n_axial_cells  = 0);

  /**
   生成环状体的体积或表面网格。环状体的轴是 $y$  -轴，而环状体的平面是 $x$  -  $z$  平面。    如果 @p dim 是3，网格将是环状体的体积，使用相当于圆的极坐标的网格，横截面上有5个单元。这个函数为所有边界面附加了一个环形Manifold，其流形ID为1；为内部单元及其所有面附加了一个CylindricalManifold，其流形ID为2（代表极坐标中的平面状态）；为表面的环形Manifold和中心的环形Manifold之间的单元附加了一个TransfiniteInterpolationManifold，其流形ID为0。 ]为3，在 $z=0$ 处切开域，6个环形单元， $R=2$ 和 $r=0.5$ 没有任何全局细化的情况下，在此给出。      @image html torus_manifold_ids.png  在这张图片中，浅灰色的阴影代表了无限插值的流形id 0，它被应用于在域边界上的环形形状和内缘之间平滑地添加新的点，在那里围绕y轴的圆柱形描述被规定。具有流形id 2的内缘显示为红色阴影。    如果 @p dim 为2，网格将描述环形的表面，这个函数将TorusManifold附加到所有的单元和面（这些单元和面的流形id标记为0）。      @param  tria 要填充的三角结构。      @param  R 圆的半径，它构成了包含单元格环的环形的中间线。必须大于  @p r.   @param  r 环状体的内半径。      @param  n_cells_toroidal 可选参数，用于设置环形方向的细胞层数。默认为6个细胞层。      @param  phi 可选参数，用于生成角度为  $0 < \varphi <= 2 \pi$  的开放环形。默认值是  $2 \pi$  ，在这种情况下会生成一个封闭的环形体。如果环状体是开放的，环状体将在垂直于环状体中心线的两个平面上被切割。  这两个平面的中心位于  $(x_1, y_1, z_1) = (R, 0, 0)$  和  $(x_2, y_2, z_2) = (R \cos(\varphi), 0, R \sin(\varphi))$  。
   * @note 为Triangulation<2,3>和Triangulation<3,3>实现。
   *
   */
  template <int dim, int spacedim>
  void
  torus(Triangulation<dim, spacedim> &tria,
        const double                  R,
        const double                  r,
        const unsigned int            n_cells_toroidal = 6,
        const double                  phi              = 2.0 * numbers::PI);

  /**
   * 这个函数在<i>xy</i>平面上产生一个正方形，中间有一个圆柱形的孔。正方形和圆形都以原点为中心。在三维空间中，这个几何体沿 $z$ 方向被挤压到区间 $[0,L]$ 。    内边界的流形ID为 $0$ ，边界ID为 $6$  。这个函数将一个PolarManifold或CylindricalManifold分别附着在2d和3d的内部边界上。其他面的边界ID为 $0, 1, 2, 3, 4$ ，或 $5$ ，按照2d或3d中面的标准顺序给出。      @image html cubes_hole.png  它在2d和3d中实现，并接受以下参数。      @param  triangulation 要填充的三角结构。    @param  inner_radius 内部孔的半径。    @param  outer_radius 正方形边长的一半。    @param  L 在 @p z-direction 中的扩展（只在3D中使用）。    @param  repetitions 沿着 @p z-direction.   @param  colorize 是否给不同的面分配不同的边界指标（见 @ref GlossColorization  "关于着色的词汇表条目"
   * ）。
   * 颜色是以词法排序给出的，平坦的面（2D中的0到3，3D中的0到5）加上弯曲的孔（2D中的4，3D中的6）。如果
   * @p colorize
   * 被设置为false，那么平面就会得到数字0，孔就会得到数字1。
   *
   */
  template <int dim>
  void
  hyper_cube_with_cylindrical_hole(Triangulation<dim> &triangulation,
                                   const double        inner_radius = .25,
                                   const double        outer_radius = .5,
                                   const double        L            = .5,
                                   const unsigned int  repetitions  = 1,
                                   const bool          colorize     = false);

  /**
   * 生成一个由同心壳组成的网格。这个函数和 GridGenerator::hyper_shell() 的主要区别是，这个函数允许不均匀的间隔（在径向方向） @ref GlossCoarseMesh "粗级单元"
   * 。    参数 @p center,  @p inner_radius, 和 @p outer_radius
   * 的行为与 GridGenerator::hyper_shell. 的前三个参数相同 @p
   * n_shells
   * 给出了要使用的壳的总数（即径向上的单元数）。 $k$
   * 第1个壳的外半径由@f[ r = r_{\mathrm{inner}} + (r_\mathrm{outer}
   *
   * - r_\mathrm{inner}) \frac{1
   *
   * - \tanh(\mathrm{skewness}(1
   *
   * - k/\mathrm{n\_shells}))} {\tanh(\mathrm{skewness})} @f]给出，其中
   * @p skewness 是控制径向方向上壳间距的参数： @p skewness
   * 接近零的值对应于均匀的间距，而 @p skewness
   * 的较大值（如 $2$ 或 $3$
   * ）对应于偏重于内半径的壳体。      @p n_cells_per_shell 与
   * GridGenerator::hyper_shell:
   * 相同，在2D中默认选择为0，将导致每个壳有8个单元（在3D中为12个）。3d中唯一有效的值是6（默认值）、12和96个单元：更多信息见
   * GridGenerator::hyper_shell 的文档。    如果 @p colorize 是
   * <code>true</code> ，那么合并后的壳的外部边界的边界ID为
   * $1$ ，内部边界的边界ID为 $0$  。
   * 例子。以下代码（关于如何可视化GNUPLOT输出的说明，请参见
   * step-10 )
   * @code
   * #include <deal.II/fe/mapping_q_generic.h>
   *
   * #include <deal.II/grid/grid_generator.h>
   * #include <deal.II/grid/grid_out.h>
   * #include <deal.II/grid/tria.h>
   *
   * #include <fstream>
   *
   * int main()
   * {
   * using namespace dealii;
   *
   * Triangulation<2> triangulation;
   * GridGenerator::concentric_hyper_shells(triangulation,
   *                                        Point<2>(),
   *                                        1.0,
   *                                        2.0,
   *                                        5u,
   *                                        2.0);
   *
   * GridOut grid_out;
   * GridOutFlags::Gnuplot gnuplot_flags(false, 10, true);
   * grid_out.set_flags(gnuplot_flags);
   *
   * const MappingQGeneric<2> mapping(3);
   * std::ofstream out("out.gpl");
   * grid_out.write_gnuplot(triangulation, out, &mapping);
   * }
   * @endcode
* 产生以下输出。      @image html concentric_hyper_shells_2d.svg
   *
   */
  template <int dim>
  void
  concentric_hyper_shells(Triangulation<dim> &triangulation,
                          const Point<dim> &  center,
                          const double        inner_radius      = 0.125,
                          const double        outer_radius      = 0.25,
                          const unsigned int  n_shells          = 1,
                          const double        skewness          = 0.1,
                          const unsigned int  n_cells_per_shell = 0,
                          const bool          colorize          = false);

  /**
   * 在3D中产生一个细胞环，将其切开、扭曲并再次粘在一起。这就产生了一种莫比乌斯环。
   * @param  tria 要处理的三角结构。    @param  n_cells
   * 循环中的单元数。必须大于4。  @param  n_rotations
   * 在胶合循环之前要进行的旋转次数（ $\pi/2$  每次）。
   * @param  R
   * 圆的半径，它构成了包含细胞环的环状体的中间线。必须大于
   * @p r.   @param  r 圆柱体弯曲在一起作为环的半径。
   *
   */
  void moebius(Triangulation<3, 3> &tria,
               const unsigned int   n_cells,
               const unsigned int   n_rotations,
               const double         R,
               const double         r);

  /**
   * 调用其他GridGenerator函数之一，从字符串 @p
   * grid_generator_function_name,
   * 中解析要调用的函数名称，从字符串 @p
   * grid_generator_function_arguments.
   * 中解析函数的参数，提供参数的字符串被传递给函数
   * Patterns::Tools::Convert<TupleTyple>::to_value(),
   * ，这里的`TupleType`是一个元组，包含*所有**GridGenerator函数的参数，包括所有可选参数。
   * 这个函数的一个使用例子是由。
   * @code
   * GridGenerator::generate_from_name_and_arguments(
   * tria,
   * "hyper_ball",
   * "0.0, 0.0 : 1 : false");
   * @endcode
   * 这里，冒号分隔了函数参数，逗号分隔了一个Point<2>参数的坐标。
   * 根据`TupleType`的算数，函数的参数可以用不同的分隔符来分隔（关于如何进行转换的细节，请参见
   * Patterns::Tuple
   * 的文档）。如果使用了错误的格式，会抛出一个异常，并将预期的格式作为错误信息输出。
   * 所有的GridGenerator函数都被支持。如果你发现有一些缺失，请在GitHub上开一个问题。
   * @param  tria 要处理的三角形  @param  grid_generator_function_name
   * 要调用的函数的名称  @param  grid_generator_function_arguments
   * 函数的参数，格式为可转换元组的字符串
   *
   */
  template <int dim, int spacedim>
  void
  generate_from_name_and_arguments(
    Triangulation<dim, spacedim> &tria,
    const std::string &           grid_generator_function_name,
    const std::string &           grid_generator_function_arguments);
  ///@}

  /**
   * @name  从其他网格中创建网格
   *
   */
  ///@{

  /**
   * 给出作为前两个参数指定的两个三角形，创建包含两个三角形的单元的三角形，并将其存储在第三个参数中。之前 @p result 的内容将被删除。  两个输入三角图中的一个也可以是 @p result 三角图。    如果几何体可以由较简单的部分组成，而这些部分存在生成 @ref GlossCoarseMesh "粗略网格 "
   * 的函数，那么这个函数最常被用来为更复杂的几何体组成网格。例如，
   * step-35 中使用的通道网格原则上可以用
   * GridGenerator::hyper_cube_with_cylindrical_hole
   * 函数创建的网格和几个矩形来创建，并使用当前函数将它们合并。矩形必须向右平移，这个任务可以用
   * GridTools::shift
   * 函数来完成（其他转换单个网格构件的工具有
   * GridTools::transform,   GridTools::rotate,  和  GridTools::scale).
   * 相距小于  @p duplicated_vertex_tolerance
   * 的顶点会被合并在一起。通常有必要将这个值设置为在某种程度上取决于输入三角形的东西。一个合理的选择是使用输入网格的所有相邻顶点之间的最小距离除以某个常数。
   * @code
   * auto min_line_length = [](const Triangulation<dim> &tria)
   *
   * -> double
   * {
   * double length = std::numeric_limits<double>::max();
   * for (const auto &cell : tria.active_cell_iterators())
   *   for (unsigned int n = 0; n < GeometryInfo<dim>::lines_per_cell; ++n)
   *     length = std::min(length, (cell->line(n)->vertex(0)
   *
   * -
   *                                cell->line(n)->vertex(1)).norm());
   * return length;
   * };
   *
   * const double tolerance = std::min(min_line_length(triangulation_1),
   *                                 min_line_length(triangulation_2)) / 2.0;
   * @endcode
   * 这将合并任何比输入网格上任何一对顶点更接近的顶点。
   * @note  两个输入三角形必须是 @ref GlossCoarseMesh "粗略网格"
   * ，即不能有任何细化单元。
   * @note
   * 该函数将两个输入三角形的单元的材料ID复制到输出三角形中。如果
   * @p copy_manifold_ids 被设置为 @p true,
   * ，流形id将被复制。边界指标永远不会被复制。换句话说，如果两个粗略的网格除了默认的边界指示器之外还有其他的东西，那么你将不得不在输出的三角图中重新手工设置边界指示器。
   * @note  该函数不会将任何流形附加到 @p result,
   * ，也不会设置任何流形ID。特别是，附加到两个输入三角形的流形将在
   * @p result 三角形中丢失。
   * @note
   * 当两个网格都来自同一个粗网格时，对细化网格的相关操作，见
   * GridGenerator::create_union_triangulation().  。
   *
   */
  template <int dim, int spacedim>
  void
  merge_triangulations(const Triangulation<dim, spacedim> &triangulation_1,
                       const Triangulation<dim, spacedim> &triangulation_2,
                       Triangulation<dim, spacedim> &      result,
                       const double duplicated_vertex_tolerance = 1.0e-12,
                       const bool   copy_manifold_ids           = false);

  /**
   * 与上述相同，但允许一次合并两个以上的三角形。
   * 下面给出一个如何使用这个函数的例子。
   * @code
   * Triangulation<2> tria_1, tria_2, tria_3;
   * // initialize tria_1, tria_2 and tria_3
   * ...
   * Triangulation<2> merged_triangulation;
   * GridGenerator::merge_triangulations({&tria_1, &tria_2, &tria_3},
   *                                     merged_triangulation,
   *                                     1.0e-10,
   *                                     false);
   * @endcode
   *
   */
  template <int dim, int spacedim>
  void
  merge_triangulations(
    const std::vector<const Triangulation<dim, spacedim> *> &triangulations,
    Triangulation<dim, spacedim> &                           result,
    const double duplicated_vertex_tolerance = 1.0e-12,
    const bool   copy_manifold_ids           = false);

  /**
   * \brief Replicate a given triangulation in multiple coordinate axes
   * @param  input The triangulation which will be replicated along the
   * coordinate axes.       @param  extents
   * 一个带有<tt>dim</tt>项的向量，指定沿每个坐标轴应存在多少个三角形的副本。
   * @param  result
   * 要创建的三角剖面。调用此函数时，它需要为空。
   * 这个函数创建一个新的三角形，等于一个<tt>dim</tt>-dimensional
   * array of copies of  @p input.   @p input
   * 的副本是通过沿坐标轴平移 @p input
   * 创建的。面的边界ID（但不是3D中的线）和所有流形ID都被复制，但流形对象没有被复制，因为大多数流形对象在Triangulation被平移后不能正确工作。
   * 为了了解这一点，请考虑以下代码。
   * @code
   * Triangulation<2> input;
   * GridGenerator::hyper_cube_with_cylindrical_hole(input);
   * Triangulation<2> output;
   * GridGenerator::replicate_triangulation(input, {3, 2}, output);
   * @endcode
*结果是 @image html replicated_tria_2d.png  而且，类似地，在三维中。
   * @code
   * Triangulation<3> input;
   * GridGenerator::hyper_cross(1, 1, 1, 2, 1, 2);
   * Triangulation<3> output;
   * GridGenerator::replicate_triangulation(input, {3, 2, 1}, output);
   * @endcode
*结果是 @image html replicated_tria_3d.png  。
   * @note  这个函数根据 @p input. 的BoundingBox确定 @p input
   * 副本的间距。 如果 @p input
   * 的边界面没有与坐标轴对齐，那么这些副本可能不会有共同的面；也就是说，这个函数是针对边界面沿坐标轴对齐的简单几何图形。
   *
   */
  template <int dim, int spacedim = dim>
  void
  replicate_triangulation(const Triangulation<dim, spacedim> &input,
                          const std::vector<unsigned int> &   extents,
                          Triangulation<dim, spacedim> &      result);

  /**
   * 给出作为前两个参数指定的两个三角形，创建包含两个三角形中最细单元的三角形，并将其存储在第三个参数中。之前
   * @p result 的内容将被删除。
   * @note  这个函数的目的是创建一个自适应细化的三角剖面，该三角剖面包含了从<i>same</i> @ref GlossCoarseMesh "粗网格 "
   * 中通过自适应细化得到的两个输入三角剖面。当人们在同一领域的单独细化网格上求解一个耦合问题的两个变量时（例如，因为这些变量在不同的地方有边界层），但随后需要计算涉及这两个变量的东西，或者希望将结果输出到一个文件中，有时需要这样的操作。在这两种情况下，为了不丢失信息，这两个解不能分别内插到其他网格上，因为那可能比计算变量的网格更粗。相反，我们需要有一个至少和两个初始网格一样精细的网格。这个函数可以计算出这样的网格。
   * @note  如果你想创建一个由另外两个 @ref GlossCoarseMesh "粗略网格 "
   * 合并而成的网格，例如，为了从较简单的几何体的网格组成一个复杂的几何体，那么这个函数就不适合你。相反，可以考虑
   * GridGenerator::merge_triangulations().  。
   * @note  这个函数假设  @p triangulation_1  和  @p  triangulation_2
   * 具有相同的流形描述。输出的三角形 @p has
   * 与这两个三角形的流形ID相同。
   * @note  两个源条件都需要完全在本地提供。
   * 换言之，它们不能是 parallel::distributed::Triangulation.
   * 类型的对象。
   *
   */
  template <int dim, int spacedim>
  void
  create_union_triangulation(
    const Triangulation<dim, spacedim> &triangulation_1,
    const Triangulation<dim, spacedim> &triangulation_2,
    Triangulation<dim, spacedim> &      result);

  /**
   * 这个函数创建了一个由第一个参数中的相同单元组成的三角形，除了第二个参数中列出的那些单元。该函数的目的是从现有的三角形描述的几何体中生成<i>subtractively</i>的几何体。一个典型的案例是一个有矩形孔的2维领域。这可以通过首先对整个域进行网格划分，然后用这个函数来摆脱位于孔洞的单元来实现。这个特殊用例的演示是
   * step-27  的一部分。同样的，你可以通过从
   * GridGenerator::hyper_cube(),
   * 开始细化一次，然后在第二个参数中加入单个单元，来创建
   * GridGenerator::hyper_L() 产生的网格。      @param[in]
   * input_triangulation
   * 作为模板的原始三角图，新的三角图将由此产生。
   * @param[in]  cells_to_remove
   * 作为第一个参数提供的三角剖面中应被移除的单元格的列表（即不应该出现在结果中。
   * @param[out]  结果 由 @p input_triangulation,
   * 中的相同单元组成的三角形，但 @p cells_to_remove.
   * 中列出的单元除外。
   * @note  与大多数GridGenerator函数不同，该函数不会为 @p
   * result, 附加任何流形，也不会设置任何流形ID。      @pre
   * 因为我们不能从头开始创建包含自适应细化单元的三角形，所以输入的三角形需要将其所有单元放在同一水平线上。通常情况下，这实际上是最粗糙的层次，但它被允许在一个已经被细化<i>globally</i>多次的三角图中通过。在这种情况下，输出的三角形将是一个只有一个层次的网格，由输入的活动单元减去第二个参数中的单元组成。然而，输入的三角形必须没有经过<i>adaptively</i>的精炼。
   *
   */
  template <int dim, int spacedim>
  void
  create_triangulation_with_removed_cells(
    const Triangulation<dim, spacedim> &input_triangulation,
    const std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>
      &                           cells_to_remove,
    Triangulation<dim, spacedim> &result);

  /**
   * 沿 @p input 方向从 $z = 0$ 到 $z =
   * \text{height}$ 挤压三角图 @p input ，并将其存储在 @p result. 中。 <em> 片 </em> ，或垂直于 $z = 0$ 平面的细胞层的数量将是 @p n_slices 片（最小是2）。 @p input 面的边界指标将被分配给 $z$ 方向的相应侧壁。底部和顶部得到下两个自由边界指标：即如果 @p input 的边界ID为 $0$ 、 $1$ 和 $42$ ，那么 $z = 0$ 的边界ID将是 $43$ ，而 $z = \text{height}$ 的边界ID将是 $44$  。    该函数默认情况下不复制流形的id。原因是在没有更多信息的情况下，没有办法在所生成的三角形的线条上设置流形ID：例如，如果 @p input 中两个流形ID不同的面在一个共享顶点相遇，那么在 </em> 中创建的与 $z$ 轴平行并通过该点的线条，就没有 <em> 先验的理由选择一个流形ID或另一个。如果 @p copy_manifold_ids 是 <code>true</code> ，那么该函数通过挑选 <em> 中首先出现的 </em> 在 @p manifold_priorities. 中的那个来设置线条流形的ID。 例如：如果 @p  ] 的 manifold_priorities 是  <code>{0, 42, numbers::flat_manifold_id}</code> ，并且所考虑的线与 manifold ids 为  <code>0</code> and <code>42</code>  的面相邻，那么该线的 manifold id 为  <code>0</code>  。正确的排序几乎总是 <ol>   <li>  设置在边界上的流形id， </li>   <li>  描述三角形中大多数单元的流形id（例如， numbers::flat_manifold_id),  和  </li>   <li>  任何对应于 TransfiniteInterpolationManifold 流形的流形id。 </li>   </ol>  特别是，由于TransfiniteInterpolationManifold在周围流形之间进行插值，其流形id通常不应设置在与不同流形id的单元格相邻的线或面上。 @p manifold_priorities 的默认值遵循这个排序（其中每个类别按升序排序）。    <ol>   <li>  与非TransfiniteInterpolationManifold的流形相关的流形ID，以及  </li>   <li>  与任何TransfiniteInterpolationManifold对象相关的流形ID。 </li>   </ol>  注意 numbers::flat_manifold_id （如果它是 @p 输入的流形ID）将始终是第一类中的最后一项。      @pre  2d输入三角形 @p input 必须是 @ref GlossCoarseMesh "粗略网格"，即不能有任何细化单元。
   * @note 由于 @p input 和 @p output 具有不同的空间尺寸，无论
   * @p copy_manifold_ids.
   * 的值如何，此函数都不会复制流形对象。
   *
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const unsigned int                     n_slices,
    const double                           height,
    Triangulation<3, 3> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});

  /**
   * extrude_triangulation()的重载，允许独立于维度的代码被编译。这个函数在调用时抛出一个错误，因为extrude_triangulation()只实现了将dim=2挤压到dim=3的三角关系。
   *
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const unsigned int                     n_slices,
    const double                           height,
    Triangulation<2, 2> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});


  /**
   * 前一个函数的重载。取一个被挤出的二维三角板。与前一个函数不同的是，这个函数采用高度和片数进行均匀挤压，这个函数采用Z轴值 @p slice_coordinates ，在这里将进行切片处理。 @p input 的面的边界指标将被分配给z方向上的相应侧壁。底部和顶部得到下两个自由边界指标。      @pre  2d输入三角形 @p input 必须是 @ref GlossCoarseMesh "粗略网格"
   * ，即不能有任何细化单元。
   * @note 由于 @p input 和 @p output
   * 具有不同的空间尺寸，因此该函数不会复制流形对象（也不会设置任何流形ID）。
   *
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const std::vector<double> &            slice_coordinates,
    Triangulation<3, 3> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});

  /**
   * extrude_triangulation()的重载，允许独立于维度的代码被编译。这个函数在调用时抛出一个错误，因为extrude_triangulation()只实现了将dim=2挤压到dim=3的三角结构。
   *
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const std::vector<double> &            slice_coordinates,
    Triangulation<2, 2> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});



  /**
   * 给定一个输入三角形 @p in_tria,
   * ，该函数生成一个新的平面三角形 @p out_tria
   * ，其中包含一个具有输入三角形的所有活动单元的单一层次。如果
   * @p spacedim1 和 @p spacedim2
   * 不同，只有顶点的最小间隔分量被复制过来。这对于从三角图<2,2>中创建三角图<2,3>，或者通过忽略顶点的Z分量将三角图<2,3>投影到三角图<2,2>中是很有用的。
   * 不对顶点进行内部检查，假设顶点在目标 @p spacedim2
   * 维空间中具有拓扑意义。如果不是这样，你在以后使用三角法时就会遇到问题。
   * 所有关于单元格流形ID和材料ID的信息都会从一个三角剖面复制到另一个三角剖面，只有边界流形ID和边界ID会从
   * @p in_tria 的面复制到 @p out_tria.
   * 的面，如果你需要指定内部面的流形ID，它们必须在三角剖面创建后手动指定。
   * 如果输入的Triangulation是 parallel::distributed::Triangulation,
   * 类型，以及输入的Triangulation包含悬空节点时，该函数将失败。
   * @param[in]  in_tria 新的平面三角图的基本输入。
   * @param[out]  out_tria 由in_tria构建的所需的扁平化三角图。
   * @note  由于 @p input 和 @p output
   * 的空间尺寸不同，此函数不会复制流形对象：您必须将新的流形对象附加到
   * @p out_tria. 。
   *
   */
  template <int dim, int spacedim1, int spacedim2>
  void
  flatten_triangulation(const Triangulation<dim, spacedim1> &in_tria,
                        Triangulation<dim, spacedim2> &      out_tria);

  /**
   将仅由超立方体单元（四边形、六面体）组成的三角剖面转换为仅由简单单元（三角形、四面体）组成的三角剖面。    作为一个例子，下面的图片显示了一组八分之一球体的三个六面体是如何被细分为四面体的，以及如何将曲面考虑在内。颜色表示边界指标是如何继承的。    @image html "convert_hypercube_to_simplex_mesh_visualization_octant.png"  一般来说，每个2D的四边形被细分为八个三角形，每个3D的六面体被细分为24个四面体，如图所示。    @image html "convert_hypercube_to_simplex_mesh_visualization.png"  材料ID和边界ID在转换时将被继承。      @param  in_tria 含有六面体元素的三角结构。    @param  out_tria 转换后的包含tet元素的三角图。
   * @note
   * 此函数不复制流形对象：您必须将现有的流形对象从 @p
   * in_tria 复制到 @p out_tria, ，例如，使用以下代码。
   * @code
   * for (const auto i : in_tria.get_manifold_ids())
   * if (i != numbers::flat_manifold_id)
   *   out_tria.set_manifold(i, in_tria.get_manifold(i));
   * @endcode
   *
   *
   */
  template <int dim, int spacedim>
  void
  convert_hypercube_to_simplex_mesh(const Triangulation<dim, spacedim> &in_tria,
                                    Triangulation<dim, spacedim> &out_tria);

  /**
   * 上述函数对一维的专门化：简单地复制三角函数。
   *
   */
  template <int spacedim>
  void
  convert_hypercube_to_simplex_mesh(const Triangulation<1, spacedim> &in_tria,
                                    Triangulation<1, spacedim> &      out_tria);


  /**
   * 名称空间Airfoil包含类和函数，以便为Joukowski或NACA机翼周围的（流）场创建C型网格。
   * @ingroup GridGenerator
   *
   */
  namespace Airfoil
  {
    /**
     * AdditionalData收集了用函数生成机翼三角图所需的所有设置
     * Airfoil::create_triangulation().  。
     *
     */
    struct AdditionalData
    {
      /**
       * 机翼的类型："NACA "或
       * "Joukowksi"，在NACA和Joukowski机翼中选择机翼的几何形状。
       *
       */
      std::string airfoil_type;

      /**
       * 定义机翼形状的NACA序列号。
       * @note  目前支持长度为4的序列号。
       * 维基百科（https://en.wikipedia.org/wiki/NACA_airfoil）对NACA序列号进行了很好的概述。
       *
       */
      std::string naca_id;

      /**
       * 茹科夫斯基圆的中心。
       * @note 中心在X轴上导致了对称的机翼。
       *
       */
      Point<2, double> joukowski_center;

      /**
       * 机翼的弦长，即从前缘到后缘的距离。
       *
       */
      double airfoil_length;

      /**
       * 从机翼弦线到网格上边界的垂直距离，即网格总高度的一半。
       *
       */
      double height;

      /**
       * 从机翼后缘到外流边界的网格长度。
       *
       */
      double length_b2;

      /**
       * 定义粗网格倾斜度HG的因素
       * 图中显示了具有两种不同倾斜度的上层粗网格
       *
       *
       *
       *
       *
       *
       * - 倾斜系数 = 0
       *
       *
       *
       *
       *
       * -> 面临HG
       *
       *
       *
       *
       *
       *
       * - 倾斜系数 = 0.5
       *
       * - >面HG' G'点的坐标由内插后的倾斜系数定义 G'(0) = G(0) + 倾斜系数 (K(0)
       *
       * - G(0))，倾斜因子在[0,1]中。o-----G--G'--K / | | / | / o | / o----o H-------o
       *
       */
      double incline_factor;

      /**
       * 偏置函数：f(x) = tanh(bx) /
       * tanh(x)，x在[0,1]中，导致接近x=1的数值被压缩。
       *
       */
      double bias_factor;

      /**
       * 全局细化的数量。
       *
       */
      unsigned int refinements;

      /**
       * 沿着机翼的左块细分的数量。
       *
       */
      unsigned int n_subdivision_x_0;

      /**
       * 沿着机翼在中间区块的细分数量。
       *
       */
      unsigned int n_subdivision_x_1;

      /**
       * 机翼右侧区块的细分数量。
       *
       */
      unsigned int n_subdivision_x_2;

      /**
       * 翼面轮廓法线的细分数。
       *
       */
      unsigned int n_subdivision_y;

      /**
       * 用于提高机翼几何形状的近似度的系数，在将提供的无距离机翼点插值到等距离机翼点时发生。当生成由等距机翼点组成的所需矢量时，它是在无液体和机翼点之间内插的。
       * 增加所提供的非液体机翼点，可以更好地接近机翼的几何形状。参数
       * "airfoil_sampling_factor
       * "由此定义了所提供的非对称点与所要求的对称点的关系。
       *
       */
      unsigned int airfoil_sampling_factor;

      /**
       * 构建器。
       *
       */
      AdditionalData();

      /**
       * 该函数添加ParameterHandler条目。              @param[in]  prm
       * Parameter handler.
       *
       */
      void
      add_parameters(ParameterHandler &prm);
    };

    /**
     * 用一个机翼周围的流场初始化给定的三角形，即一个接近Joukowski或NACA（4位数）机翼的C型网格。
     * 用户可以通过为结构AdditionalData提供输入参数来指定机翼的几何形状和网格设置。
     * 因此，用户可以在不同类型的Joukowski或NACA机翼中选择不同的弦长、远场尺寸和网格密度。
     * @note
     * 该函数创建一个细化的网格（全局细化的数量可以由用户指定）。没有附加流形。最终网格中的顶点会被这个函数移动到正确的位置。
*  @note  这个函数目前只适用于二维，但当然也可以用  GridGenerator::extrude().   @param[out]  tria 要建立的三角结构。调用此函数时，它必须是空的。      @param[in]  additional_data 网格的配置。        <style>div.image img[src="https://www.dealii.org/images/grids/airfoils_naca_joukowski.png"]{width:50%;}</style>\endhtmlonly  @image html https://www.dealii.org/images/grids/airfoils_naca_joukowski.png
     *
     */
    template <int dim>
    void
    create_triangulation(
      Triangulation<dim, dim> &tria,
      const AdditionalData &   additional_data = AdditionalData());



    /**
     * 与上述相同，但在远场的上、下两面施加周期性边界条件。
     * @note  这个函数目前只对二维实现。          @param[out]
     * tria
     * 要创建的三角结构。在调用此函数时，它需要是空的。
     * @param[out]  periodic_faces 上下水平边界的周期面。
     * @param[in]  additional_data 网格的配置。
     *
     */
    template <int dim>
    void
    create_triangulation(
      Triangulation<dim, dim> &                            tria,
      std::vector<GridTools::PeriodicFacePair<
        typename Triangulation<dim, dim>::cell_iterator>> &periodic_faces,
      const AdditionalData &additional_data = AdditionalData());

  } // namespace Airfoil

  /**
   * 从两个斜对角的角点 @p p1 和 @p p2.
   * 创建一个坐标平行的砖块，坐标方向 @p i
   * 的顶点数量由<tt>repetitions[i]+1</tt>决定。
   * @note
   * 这个函数在内部将4/8个顶点连接到四边形/六面体单元，并将这些单元细分为2/5个三角形/四面体单元。
   * @note  目前，这个函数只对`dim==spacedim`起作用。
   * @ingroup simplex
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle_with_simplices(
    Triangulation<dim, spacedim> &   tria,
    const std::vector<unsigned int> &repetitions,
    const Point<dim> &               p1,
    const Point<dim> &               p2,
    const bool                       colorize = false);

  /**
   * 用一个超立方体（二维的正方形和三维的立方体）初始化给定的三角形，该超立方体在每个方向上由
   * @p repetitions 个单元组成。
   * 超立方体的体积是目前维数中的张量乘积区间
   * $[left,right]^{\text{dim}}$
   * ，其中极限作为参数给出。它们默认为零和一，然后产生单位超立方体。
   * @note
   * 这个函数在内部将4/8个顶点连接到四边形/六面体单元，并将其细分为2/5个三角形/四面体单元。
   * @ingroup simplex
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_cube_with_simplices(Triangulation<dim, spacedim> &tria,
                                       const unsigned int repetitions,
                                       const double       p1       = 0.0,
                                       const double       p2       = 1.0,
                                       const bool         colorize = false);

  ///@}

  /**
   * @name  创建低维网格 从高维网格的一部分创建。
   *
   */
  ///@{

#ifdef _MSC_VER
  // Microsoft's VC++ has a bug where it doesn't want to recognize that
  // an implementation (definition) of the extract_boundary_mesh function
  // matches a declaration. This can apparently only be avoided by
  // doing some contortion with the return type using the following
  // intermediate type. This is only used when using MS VC++ and uses
  // the direct way of doing it otherwise
  template <template <int, int> class MeshType, int dim, int spacedim>
  struct ExtractBoundaryMesh
  {
    using return_type =
      std::map<typename MeshType<dim - 1, spacedim>::cell_iterator,
               typename MeshType<dim, spacedim>::face_iterator>;
  };
#endif

  /**
   * 这个函数实现了边界子网格的提取。
   * 给定一个<dim,spacedim>三角网格（"体网格"），该函数提取其边界的一个子集（"面网格"）。
   * 要提取的边界是由一个 boundary_ids 列表指定的。
   * 如果没有指定，整个边界将被提取出来。该函数在
   * step-38  中使用。
   * 该函数还建立了一个映射，将表面网格上的单元格与体积网格上的相应面连接起来。这个映射是该函数的返回值。
   * @note
   * 上面概述的算法假定所有在较高细化水平上的面总是具有与它们的父面完全相同的边界指标。因此，我们可以从粗层次的面开始，在此基础上建立曲面网格。将这个函数扩展到将细化级面的边界指标复制到其相应的表面网格单元并不是很困难，例如，在弯曲边界的情况下适应不同的几何描述（但目前还没有实现）。
   * @note 由于 @p volume_mesh 和 @p surface_mesh
   * 有不同的空间维度，没有流形对象被这个函数复制：你必须将新的流形对象附加到
   * @p surface_mesh. 。
   *
   */
  template <template <int, int> class MeshType, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename MeshType<dim - 1, spacedim>::cell_iterator,
           typename MeshType<dim, spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<MeshType, dim, spacedim>::return_type
#endif
  extract_boundary_mesh(const MeshType<dim, spacedim> &     volume_mesh,
                        MeshType<dim - 1, spacedim> &       surface_mesh,
                        const std::set<types::boundary_id> &boundary_ids =
                          std::set<types::boundary_id>());

  ///@}


  /**
   * @name  异常情况
   *
   */
  ///@{


  /**
   * 异常情况
   *
   */
  DeclException0(ExcInvalidRadii);
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidRepetitions,
                 int,
                 << "The number of repetitions " << arg1 << " must be >=1.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidRepetitionsDimension,
                 int,
                 << "The vector of repetitions  must have " << arg1
                 << " elements.");

  /**
   * 输入方向不正确的异常情况。
   *
   */
  DeclExceptionMsg(ExcInvalidInputOrientation,
                   "The input to this function is oriented in a way that will"
                   " cause all cells to have negative measure.");
  ///@}

#ifndef DOXYGEN
  // These functions are only implemented with specializations; declare them
  // here
  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<1> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<2> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<3> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void channel_with_cylinder(Triangulation<1> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);

  template <>
  void channel_with_cylinder(Triangulation<2> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);

  template <>
  void channel_with_cylinder(Triangulation<3> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);
#endif
} // namespace GridGenerator



DEAL_II_NAMESPACE_CLOSE

#endif


