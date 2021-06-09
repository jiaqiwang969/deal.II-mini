//include/deal.II-translator/grid/grid_out_0.txt
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

#ifndef dealii_grid_out_h
#define dealii_grid_out_h



#include <deal.II/base/config.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class ParameterHandler;
template <int dim, int spacedim>
class Triangulation;
template <int dim, int spacedim>
class Mapping;
#endif


/**
 * 在这个命名空间中，我们定义了几个结构，用来描述可以给网格输出例程的标志，以修改写进文件的网格的默认装备。更多细节请参见不同的子类和GridOut类的文档。
 *
 *
 * @ingroup output
 *
 *
 */
namespace GridOutFlags
{
  /**
   * 用于OpenDX格式的网格输出的标志。
   * @ingroup output
   *
   */
  struct DX
  {
    /**
     * 写入单元格。
     *
     */
    bool write_cells;

    /**
     * 写入面。
     *
     */
    bool write_faces;

    /**
     * 写出有直径的场。
     *
     */
    bool write_diameter;

    /**
     * 写入有面积/体积的字段。
     *
     */
    bool write_measure;

    /**
     * 写所有的面，包括内部面。如果<tt>false</tt>，只写边界面。
     *
     */
    bool write_all_faces;

    /**
     * 构造函数。
     *
     */
    DX(const bool write_cells     = true,
       const bool write_faces     = false,
       const bool write_diameter  = false,
       const bool write_measure   = false,
       const bool write_all_faces = true);

    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };

  /**
   * 描述MSH格式的输出细节的标志。
   * @ingroup output
   *
   */
  struct Msh
  {
    /**
     * 写入网格时，如果边界面的指示器不是默认的边界指示器（即零），则明确写入边界面。
     * 这是必要的，如果你以后想重新读取网格，并希望得到三角形边界的不同部分的相同边界指示。
     * 如果你只想写出三角图来查看或打印，则没有必要。
     * 默认情况下。  @p false.
     *
     */
    bool write_faces;
    /**
     * 写入网格时，如果边界线的指示器不是默认的边界指示器（即零），则明确写入边界线。
     * 这是必要的，如果你以后想重新读取网格，并希望得到三角形边界的不同部分的相同边界指示。
     * 如果你只想写三角剖面图来查看或打印，则没有必要。
     * 只有在<tt>dim==3</tt>时才使用，其他情况下都忽略。
     * 默认值。  @p false.
     *
     */
    bool write_lines;

    /**
     * 构造函数。
     *
     */
    Msh(const bool write_faces = false, const bool write_lines = false);
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   * 描述UCD格式输出细节的标志。
   * @ingroup output
   *
   */
  struct Ucd
  {
    /**
     * 在文件的开头写一个注释，说明创建日期和其他一些数据。
     * 虽然UCD格式（和AVS程序）支持这一点，但其他一些程序对此感到困惑，所以默认情况是不写序言。然而，可以用这个标志写一个序言。
     * 默认值。  <code>false</code>  .
     *
     */
    bool write_preamble;

    /**
     * 写入网格时，如果边界面的指示器不是默认的边界指示器（即零），则明确写入边界面。
     * 这是必要的，如果你以后想重新读取网格，并希望得到三角形边界的不同部分的相同边界指示。
     * 如果你只想写出三角图来查看或打印，则没有必要。
     * 默认情况下。  @p false.
     *
     */
    bool write_faces;

    /**
     * 写入网格时，如果边界线的指示器不是默认的边界指示器（即零），就明确写入边界线。
     * 这是必要的，如果你以后想重新读取网格，并希望得到三角形边界的不同部分的相同边界指示。
     * 如果你只想写出三角图来查看或打印它，则没有必要。
     * 如果<tt>dim!=3</tt>，这个指令会被忽略。
     * 默认情况下。  @p false.
     *
     */
    bool write_lines;

    /**
     * 构造函数。
     *
     */
    Ucd(const bool write_preamble = false,
        const bool write_faces    = false,
        const bool write_lines    = false);

    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   * 描述GNUPLOT格式输出细节的标志。
   * @ingroup output
   *
   */
  struct Gnuplot
  {
    /**
     * 将每个单元格的编号写进输出文件，然后再开始写它所组成的行，作为注释。如果你想找出网格的细节，例如你知道编号的单元格的位置，这可能很有用。然而，它大大增加了输出的大小。
     * 默认值。  @p false.
     *
     */
    bool write_cell_numbers;

    /**
     * 点的数量， <em> 不包括 </em>
     * 顶点，用于绘制曲线。由于GNUPLOT只能绘制直线，将这个数字设置为大于0的值（对于精细的网格来说，4或5通常就足够了）会使绘图看起来是弯曲的，尽管它不是。
     *
     */
    unsigned int n_extra_curved_line_points;

    /**
     * 表示是否要用<tt>n_extra_curved_line_points</tt>线段绘制内部线的布尔值。
     *
     */
    bool curved_inner_cells;

    /**
     * 标志。如果为真，那么在写<tt>spacedim =
     * 3</tt>输出时，在边界面上写<tt>2*n_extra_curved_line_points</tt>额外的线。
     * 当<tt>spacedim =
     * 2</tt>时，设置这个选项没有影响，因为在这种情况下，边界面是线，输出额外的线是没有意义的。
     * @note  对于<tt>dim =
     * 2</tt>的情况，这个选项还没有实现。然而，为了向后兼容，这将不会引发运行时错误。
     *
     */
    bool write_additional_boundary_lines;

    /**
     * 构造函数。
     *
     */
    Gnuplot(const bool         write_cell_number               = false,
            const unsigned int n_extra_curved_line_points      = 2,
            const bool         curved_inner_cells              = false,
            const bool         write_additional_boundary_lines = true);

    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };

  /**
   * 描述封装后的postscript的输出细节的标志。
   * 在这个结构中，列出了所有尺寸通用的标志。只针对一个空间维度的标志被列在派生类中。
   * 默认情况下，图片的尺寸是按比例缩放的，因此宽度等于300单位。
   * @ingroup output
   *
   */
  struct EpsFlagsBase
  {
    /**
     * 枚举表示是否应该进行缩放，使给定的 @p size
     * 等于结果图片的宽度或高度。
     *
     */
    enum SizeType
    {
      /**
       * 按宽度缩放。
       *
       */
      width,
      /**
       * 随高度缩放。
       *
       */
      height
    };

    /**
     * 见上文。默认是 @p width. 。
     *
     */
    SizeType size_type;

    /**
     * 输出的宽度或高度以postscript单位给出，这通常是由奇怪的单位1/72英寸给出。这是高还是宽由标志
     * @p size_type. 指定 默认为300。
     *
     */
    unsigned int size;

    /**
     * 一行的宽度，以postscript为单位。默认为0.5。
     *
     */
    double line_width;

    /**
     * 设置为 @p
     * user_flag 的线条是否应该用不同的颜色（红色）绘制？
     * 参见 @ref GlossUserFlags ，了解有关用户标志的信息。
     *
     */
    bool color_lines_on_user_flag;

    /**
     * 边界面上的点的数量，除了该面的顶点外，还被绘制出来。
     * 这个数字只在所使用的映射不是简单的标准 $Q_1$
     * 映射（即MappingQGeneric(1)类型的对象）时使用，该映射可能将单元格的边缘描述为弯曲，然后将使用线段进行近似，中间点的数量由当前变量描述。
     *
     */
    unsigned int n_boundary_face_points;

    /**
     * 线条是否应该根据其细化程度来着色？这将覆盖所有级别的color_lines_on_user_flag，除了级别0。
     * 颜色是：0级：黑色，其他级别：从蓝色到红色的彩虹色。
     *
     */
    bool color_lines_level;

    /**
     * 构造函数。
     *
     */
    EpsFlagsBase(const SizeType     size_type                = width,
                 const unsigned int size                     = 300,
                 const double       line_width               = 0.5,
                 const bool         color_lines_on_user_flag = false,
                 const unsigned int n_boundary_face_points   = 2,
                 const bool         color_lines_level        = false);
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   * 描述所有尺寸的封装postscript的输出细节的标志，下面没有明确地专门说明。在基类中列出了一些对所有维度通用的标志。
   * 这个类实际上并不存在，我们只在这里声明一般的模板，并在下面声明明确的专门化。
   * @ingroup output
   *
   */
  template <int dim>
  struct Eps
  {};

  /**
   * 专门针对一个空间维度的网格输出的标志。
   * @ingroup output
   *
   */
  template <>
  struct Eps<1> : public EpsFlagsBase
  {
    /**
     * 构造函数。
     *
     */
    Eps(const SizeType     size_type                = width,
        const unsigned int size                     = 300,
        const double       line_width               = 0.5,
        const bool         color_lines_on_user_flag = false,
        const unsigned int n_boundary_face_points   = 2);
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   * 专门针对二维空间的网格输出的标志。
   * @ingroup output
   *
   */
  template <>
  struct Eps<2> : public EpsFlagsBase
  {
    /**
     * 如果这个标志被设置，那么我们就把单元格的编号放到每个单元格的中间。默认值是不这样做。
     * 写入的单元格编号的格式是<tt>level.index</tt>，或者干脆是
     * @p index, ，取决于以下标志的值。
     *
     */
    bool write_cell_numbers;
    /**
     * 如果单元格编号应被写入，使用上述标志，那么这个标志的值决定了格式是<tt>level.index</tt>，还是简单的
     * @p index.  如果 @p true, ，则采取第一种格式。默认是 @p
     * true.  如果 @p write_cell_numbers 是 @p false.
     * ，该标志显然没有作用。
     *
     */
    bool write_cell_number_level;

    /**
     * 顶点编号可以被写入顶点。这是由以下标志控制的。默认是
     * @p false. 。
     *
     */
    bool write_vertex_numbers;

    /**
     * 构造函数。
     *
     */
    Eps(const SizeType     size_type                = width,
        const unsigned int size                     = 300,
        const double       line_width               = 0.5,
        const bool         color_lines_on_user_flag = false,
        const unsigned int n_boundary_face_points   = 2,
        const bool         write_cell_numbers       = false,
        const bool         write_cell_number_level  = true,
        const bool         write_vertex_numbers     = false,
        const bool         color_lines_level        = false);
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };

  /**
   * 三维空间网格输出的特定标志。
   * @ingroup output
   *
   */
  template <>
  struct Eps<3> : public EpsFlagsBase
  {
    /**
     * 线条原点观测器对Z轴的角度，单位是度。
     * 默认为Gnuplot默认的60。
     *
     */
    double azimut_angle;

    /**
     * 从上面看时，观察者的位置投射到x-y平面上围绕z轴旋转的角度，为正数。
     * 单位是度，零等于高于或低于负y轴的位置。
     * 默认值是Gnuplot默认的30。
     *
     */
    double turn_angle;

    /**
     * 构造函数。
     *
     */
    Eps(const SizeType     size_type                = width,
        const unsigned int size                     = 300,
        const double       line_width               = 0.5,
        const bool         color_lines_on_user_flag = false,
        const unsigned int n_boundary_face_points   = 2,
        const double       azimut_angle             = 60,
        const double       turn_angle               = 30);
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };

  /**
   * 用于XFig输出的标志。
   * @ingroup output
   *
   */
  struct XFig
  {
    /**
     * 绘制边界线。默认为true。
     *
     */
    bool draw_boundary;

    /**
     * 一个枚举，用于决定哪个字段用于给单元格着色。
     *
     */
    enum Coloring
    {
      /// Convert the material id into the cell color
      material_id,
      /// Convert the level into the cell color
      level_number,
      /// Convert the global subdomain id into the cell color
      subdomain_id,
      /// Convert the level subdomain id into the cell color
      level_subdomain_id
    } color_by;

    /**
     * 代码级别到深度。默认为真。如果是假的，颜色取决于材料或边界ID。
     * 如果这个值为真，对象的深度为900级。
     *
     */
    bool level_depth;

    /**
     * 弯曲边界的附加点。默认是没有。
     *
     */
    unsigned int n_boundary_face_points;

    /**
     * 图形的缩放。默认是单位长度为一英寸。
     *
     */
    Point<2> scaling;

    /**
     * 图形的偏移。在缩放之前，坐标被移到这个值上。默认是每个方向上的零。
     *
     */
    Point<2> offset;

    /**
     * 填充单元格的样式。默认是实心填充（20）。这个值会被转发到XFig的折线对象的相应字段<tt>fill_style</tt>中而不被改变。
     *
     */
    int fill_style;

    /**
     * 绘制多边形边界线的样式。默认为实心（0），并转发到XFig。
     *
     */
    int line_style;

    /**
     * 多边形边界线的厚度。默认为1。将其设置为0，以避免非常细的网格出现边界线。
     *
     */
    int line_thickness;

    /**
     * 在边界上画线的样式。默认为实体（0）。
     *
     */
    int boundary_style;

    /**
     * 边界线的厚度。默认为3。
     *
     */
    int boundary_thickness;

    /**
     * 构造函数。
     *
     */
    XFig();
    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   控制SVG输出的标志。    下图是SVG标志的放大图，说明了SVG标志能够产生的效果。这些确切的设置与测试中使用的设置相同  <code>tests/grid/grid_out_svg_02.cc</code>  ，只是增加了标志  <code>svg_flags.label_boundary_id = true;</code>  。      @image html svg_flags.png
   * @ingroup output
   *
   */
  struct Svg
  {
    /**
     * 绘图的高度，以SVG为单位，如果为0，则由宽度计算。默认为1000。
     *
     */
    unsigned int height;

    /**
     * 绘图的宽度。如果为零，则由高度自动计算（默认）。
     *
     */
    unsigned int width;

    /**
     * 单元间的线的厚度。
     *
     */
    unsigned int line_thickness;
    /**
     * 边界处线条的厚度。
     *
     */
    unsigned int boundary_line_thickness;

    /**
     * 绘图区域周围的边距。
     *
     */
    bool margin;

    /**
     * 一个描述所有可能的背景风格的 "enum"。
     *
     */
    enum Background
    {
      /**
       * 使用SVG的透明值。
       *
       */
      transparent,

      /**
       * 使用白色背景。
       *
       */
      white,

      /**
       * 使用从白色（顶部）到钢蓝色（底部）的梯度，并添加日期和时间以及deal.II的标志。自动画出一个空白。
       *
       */
      dealii
    };

    /**
     * 用于网格背景的样式。
     *
     */
    Background background;

    // View angles for the perspective view of the grid; Default is 0, 0 (top
    // view).

    /**
     * 从???测量的方位角，单位是度。默认为0。
     *
     */
    int azimuth_angle;

    /**
     * 从xy平面垂直方向上的角度。默认为0。
     *
     */
    int polar_angle;

    /**
     * 单元的着色。
     *
     */
    enum Coloring
    {
      /// No cell coloring
      none,
      /// Convert the material id into the cell color (default)
      material_id,
      /// Convert the level number into the cell color
      level_number,
      /// Convert the subdomain id into the cell color
      subdomain_id,
      /// Convert the level subdomain id into the cell color
      level_subdomain_id
    };

    Coloring coloring;

    /// Interpret the level number of the cells as altitude over the x-y-plane
    /// (useful in the perspective view).
    bool convert_level_number_to_height;

    /**
     * 决定层间垂直距离的因素（默认=0.3。
     *
     */
    float level_height_factor;

    /**
     * 用于单元格注释的字体缩放。默认为1。
     *
     */
    float cell_font_scaling;
    /**
     * 在每个单元格中写入级别编号。默认为false。
     *
     */
    bool label_level_number;

    /**
     * 在每个单元格中写入单元格索引。默认为false。
     *
     */
    bool label_cell_index;

    /**
     * 写入每个单元格的材料ID。默认为false。
     *
     */
    bool label_material_id;

    /**
     * 写入每个单元格的子域ID。默认为false。
     *
     */
    bool label_subdomain_id;

    /**
     * 写入每个单元格的子域ID。默认为false。
     *
     */
    bool label_level_subdomain_id;

    /**
     * 在相应的边界边上的圆圈中写入每个边界面的边界ID。默认为false。
     * 注意：根据图像查看器的选择，边界ID标签可能不会出现在圆的中心。
     *
     */
    bool label_boundary_id;

    /**
     * 在绘制的网格旁边画一个颜色条，与选择的单元格的颜色有关。
     *
     */
    bool draw_colorbar;

    /**
     * 在绘制的网格旁边画一个图例，解释单元格的标签。
     *
     */
    bool draw_legend;

    /**
     * 构造函数。
     *
     */
    Svg(const unsigned int line_thickness                 = 2,
        const unsigned int boundary_line_thickness        = 4,
        const bool         margin                         = true,
        const Background   background                     = white,
        const int          azimuth_angle                  = 0,
        const int          polar_angle                    = 0,
        const Coloring     coloring                       = level_number,
        const bool         convert_level_number_to_height = false,
        const bool         label_level_number             = false,
        const bool         label_cell_index               = false,
        const bool         label_material_id              = false,
        const bool         label_subdomain_id             = false,
        const bool         draw_colorbar                  = false,
        const bool         draw_legend                    = false,
        const bool         label_boundary_id              = false);
  };

  /**
   * 用于MathGL格式的网格输出的标志。
   * @ingroup output
   *
   */
  struct MathGL
  {
    /**
     * 构造函数。
     *
     */
    MathGL();

    /**
     * 在图形周围画一个包围盒。
     *
     */
    bool draw_bounding_box;

    /**
     * 在ParameterHandler中声明参数。
     *
     */
    static void
    declare_parameters(ParameterHandler &param);

    /**
     * 解析ParameterHandler的参数。
     *
     */
    void
    parse_parameters(ParameterHandler &param);
  };


  /**
   * 用于Vtk格式的网格输出的标志。这些标志与
   * DataOutBase::VtkFlags. 中声明的标志相同。
   * @ingroup output
   *
   */
  struct Vtk : public DataOutBase::VtkFlags
  {
    /**
     * 默认构造函数。
     *
     */
    Vtk(const bool output_cells         = true,
        const bool output_faces         = true,
        const bool output_edges         = true,
        const bool output_only_relevant = true)
      : output_cells(output_cells)
      , output_faces(output_faces)
      , output_edges(output_edges)
      , output_only_relevant(output_only_relevant)
    {}

    /**
     * 输出单元。
     *
     */
    bool output_cells;

    /**
     * 输出面。
     *
     */
    bool output_faces;

    /**
     * 输出共面/边。
     *
     */
    bool output_edges;

    /**
     * 只输出与默认设置不同的面孔/共同面孔（例如边界_id）。
     *
     */
    bool output_only_relevant;
  };


  /**
   * 用于Vtu格式的网格输出的标志。这些标志与
   * DataOutBase::VtkFlags,
   * 中声明的标志相同，但增加了一个标志，决定你是否要在vtu文件（实际上是一个xml文件）中增加一个包含整个三角测量序列化的条目。
   * @ingroup output
   *
   */
  struct Vtu : public DataOutBase::VtkFlags
  {
    Vtu(const bool serialize_triangulation = false)
      : serialize_triangulation(serialize_triangulation)
    {}

    /**
     * 将序列化的三角图也添加到vtu文件中。
     *
     */
    bool serialize_triangulation;
  };
} // namespace GridOutFlags



/**
 * 这个类提供了一种方法，可以将三角图以不同的格式输出到文件中。参见枚举
 * GridOut::OutputFormat 中的格式列表和相应的输出函数名称。
 * 用法很简单：你可以使用直接的形式
 *
 * @code
 * std::ofstream output_file("some_filename");
 * GridOut().write_gnuplot (tria, output_file);
 * @endcode
 * 如果你知道你想要哪种格式，或者你想让格式成为一个运行时参数，你可以编写
 *
 * @code
 * GridOut::OutputFormat grid_format =
 *   GridOut::parse_output_format(get_format_name_from_somewhere());
 * std::ofstream output_file("some_filename"
 *                           + GridOut::default_suffix(output_format));
 * GridOut().write (tria, output_file, output_format);
 * @endcode
 * 函数<tt>get_output_format_names()</tt>以ParameterHandler类可理解的字符串提供了一个可能的输出格式名称的列表。
 * 注意，在这里，我们创建了一个未命名的GridOut类型的对象，并调用了它的一个<tt>write_*</tt>函数。这看起来好像真的可以让各自的函数
 * @p static.
 * 这不是为了让参数以一种与允许在运行时通过通用 @p write
 * 函数选择正确输出格式的方案相兼容的方式传递给不同的输出函数。
 * 为了解释这一点，考虑到每个函数都有一个或多个额外的参数，给出输出的细节，例如3D网格的观众位置，线的厚度等等。虽然这将允许每个输出函数具有它所需要的任何灵活性，但它不允许我们使用通用函数
 * @p write
 * ，该函数被赋予一个决定输出格式的参数，因为给它一个支持每一种输出格式的参数列表是不切实际的，然后它可以将其传递给各自的输出函数。
 * 相反，我们选择让这个类GridOut的每个对象为每个支持的输出格式拥有一组参数。这些参数被收集在GridOutFlags命名空间中声明的结构
 * GridOutFlags::Eps(),   GridOutFlags::Gnuplot(),
 * 等中，你可以像这样设置你的首选标志。
 *
 * @code
 * GridOut grid_out;
 * GridOutFlags::Ucd ucd_flags;
 * ...    // set some fields in ucd_flags
 * grid_out.set_flags (ucd_flags);
 * ...
 * ...    // write some file with data_out
 * @endcode
 * 各个输出函数然后使用如此设置的标志。默认情况下，它们被设置为合理的值，如上面和不同标志结构的文档中所述。重置标志可以通过调用<tt>grid_out.set_flags
 * (GridOutFlags::Ucd());</tt>,
 * 来完成，因为每个标志结构的默认构造函数都将参数设置为初始值。
 * 这种方法的好处是，可以根据你的需要改变一个或多个输出格式的标志，然后再使用通用
 * @p write
 * 函数；然后调用的实际输出函数将使用之前设置的标志。
 * 注意，一些描述不同输出格式标志的结构是空的，因为各自的格式不支持任何标志。反正结构和
 * @p set_flags
 * 函数是提供的。还需要注意的是，有些结构在这个类所支持的尺寸之间可能有所不同；然后它们有一个模板参数，像往常一样。
 *
 *
 * @ingroup grid
 *
 * @ingroup output
 *
 */
class GridOut
{
public:
  /**
   * 为每个不同的输出格式声明一个名称。这些名称被通用输出函数write()用来确定实际的输出格式。
   *
   */
  enum OutputFormat
  {
    /// Do nothing in write()
    none,
    /// write() calls write_dx()
    dx,
    /// write() calls write_gnuplot()
    gnuplot,
    /// write() calls write_eps()
    eps,
    /// write() calls write_ucd()
    ucd,
    /// write() calls write_xfig()
    xfig,
    /// write() calls write_msh()
    msh,
    /// write() calls write_svg()
    svg,
    /// write() calls write_mathgl()
    mathgl,
    /// write() calls write_vtk()
    vtk,
    /// write() calls write_vtu()
    vtu
  };

  /**
   * 构造函数。
   *
   */
  GridOut();

  /**
   * 以OpenDX格式写入三角图。
   * 单元或面与它们的水平面和它们的材料ID或边界指示器一起被写入。
   * 在二维度为1的情况下不执行。
   *
   */
  template <int dim, int spacedim>
  void
  write_dx(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;

  /**
   * 用GNUPLOT格式写三角剖面图。
   * 在GNUPLOT格式中，每个单元被写成其限定线的序列。除了线的端点坐标之外，单元格的水平和材料也被附加到每一行的输出中。因此，如果你让GNUPLOT把一个2D网格画成一个3D图，你会看到更精细的单元与更不精细的单元相对应。
   * 另外，如果你在三维网格中画一个切面，你可以在与切面正交的方向挤出细化水平。同样的方法也可以用在材料ID上，材料ID被绘制在层次之后。
   * 这个功能的一个更有用的应用是：如果你使用GNUPLOT命令（这里是指2d网格
   * @verbatim
   * splot [:][:][2.5:3.5] "grid_file.gnuplot"
   * @endverbatim
   * 那么整个x-和y-范围将被绘制出来，也就是整个网格，但只有那些z值在2.5和3.5之间的线条。由于z值被选择为一个单元格所属的级别，这就导致在这个例子中只绘制那些属于第3级的单元格。这样一来，就很容易产生不同级别网格的图。
   * @p mapping
   * 是一个指针，用于边界上的单元格转换的映射。如果为零，则使用标准的Q1映射。
   * 控制输出的额外标志的名称和值可以在 GridOutFlags::Gnuplot
   * 类的文档中找到，其中还描述了一些关于一维情况的注意事项。
   *
   */
  template <int dim, int spacedim>
  void
  write_gnuplot(const Triangulation<dim, spacedim> &tria,
                std::ostream &                      out,
                const Mapping<dim, spacedim> *      mapping = nullptr) const;

  /**
   * 以msh格式写入三角剖面图。
   * Msh是%Gmsh使用的格式，它在%Gmsh用户指南中有所描述。除了通常只输出网格外，你可以通过额外的标志（见下文，以及
   * GridOutFlags::Msh()
   * 类的文档）决定是否应将边界指示器为非零的边界面明确写入文件。这很有用，如果你以后想重新读取网格，因为<tt>deal.II</tt>默认将边界指示器设置为零；因此，为了获得和以前一样的三角测量，你必须明确指定具有不同边界指示器的面，这是由这个标志完成的。
   * 其他控制输出的标志的名称和值可以在 GridOutFlags::Msh()
   * 类的文档中找到。    在维度为1的情况下也能工作。
   *
   */
  template <int dim, int spacedim>
  void
  write_msh(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;

#ifdef DEAL_II_GMSH_WITH_API
  /**
   * 用gmsh API支持的任何格式来写三角剖面。    Gmsh
   * API允许通过其C++
   * API以多种格式写入其输出。这个函数将三角图对象翻译成gmsh的实体集合，并调用
   * gmsh::write()
   * 方法，将文件名作为参数传递。该方法为每一对独特的非默认流形id和边界id生成一个不同的实体，并为每一个独特的组合写一个gmsh物理组，允许你使用以字符串为参数的
   * GridIn::read_msh() 方法读回三角图。
   * 特别是，所有具有非默认边界ID或非默认流形ID的单元对象都被归入一个独特的物理标签，其名称包含边界和流形指标。名称是用
   * Patterns::Tools::to_value() 应用于 `std::map<std::string,
   * int>`构建的，其中键是`材料ID`、`边界ID`或`流形ID`，即一个材料ID为1、流形ID为3的单元将被分组在一个物理标签（其编号未被指定），命名为`材料ID:1,
   * 流形ID:3`。
   * 例如，用超球网格精炼一次后调用该方法，将导致输出文件中定义以下物理标签。
   * @code
   * MeshFormat
   * 4.1 0 8
   * \$EndMeshFormat
   * \$PhysicalNames
   * 3
   * 1 1 "ManifoldID:0"
   * 1 2 "BoundaryID:-1, ManifoldID:1"
   * 2 3 "ManifoldID:1"
   * \$EndPhysicalNames
   * \$Entities
   * ...
   * @endcode
   * 特殊的边界ID`-1`是用来表示内部边界的。只要有必要指定一个非平面流形的id，就必须指定内部边界。
   *
   */
  template <int dim, int spacedim>
  void
  write_msh(const Triangulation<dim, spacedim> &tria,
            const std::string &                 filename) const;
#endif

  /**
   * 以ucd格式写入三角图。
   * UCD（非结构化单元数据）是AVS和其他一些程序使用的格式。它在AVS的开发者指南中有所描述。除了通常只输出网格之外，你可以通过额外的标志（见下文，以及
   * GridOutFlags::Ucd()
   * 类的文档）决定是否将边界指示器为非零的边界面明确写入文件。这很有用，如果你以后想重新读取网格，因为<tt>deal.II</tt>默认将边界指示器设置为零；因此，为了获得和以前一样的三角测量，你必须明确指定具有不同边界指示器的面，这是由这个标志完成的。
   * 其他控制输出的标志的名称和值可以在 GridOutFlags::Ucd()
   * 类的文档中找到。    也适用于维度为1的情况。
   *
   */
  template <int dim, int spacedim>
  void
  write_ucd(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;

  /**
   * 以封装的postscript格式编写三角图。
   * 在这种格式中，三角图的每一行都是单独写的。我们对图片进行缩放，使X值或Y值的范围在0和固定大小之间。另一个轴的比例是相同的因素。哪一个轴被用来计算比例，以及它应适合的盒子的大小，由输出标志决定（见下文，以及
   * GridOutFlags::Eps() 类的文档）。
   * 界限盒的四边都接近三角形，没有额外的框架。线宽默认选择为0.5，但可以改变。线宽要与图片的延伸部分相比较，其默认值为300。
   * 标志 @p color_lines_on_user_flag 允许在设置了 @p
   * user_flag的情况下绘制红色的线条。在输出文件的序言中，黑色和红色被定义为
   * @p b 和 @p r ，可以根据需要在那里进行更改。      @p
   * mapping
   * 是一个指向用于边界处单元格转换的映射的指针。如果为零，则使用标准的Q1映射。
   * 控制输出的额外标志的名称和值可以在 GridOutFlags::Eps()
   * 类的文档中找到。特别是对于三维网格的观点在这里是很重要的。
   * 对于一维的情况没有实现。
   *
   */
  template <int dim, int spacedim>
  void
  write_eps(const Triangulation<dim, spacedim> &tria,
            std::ostream &                      out,
            const Mapping<dim, spacedim> *      mapping = nullptr) const;

  /**
   * 写二维XFig-file。
   * 这个函数将所有的网格单元写成多边形，也可以选择写成边界线。几个参数可以通过XFigFlags控制对象来调整。
   * 如果层次被编码为深度，完整的网格层次就会在其父辈之前绘制出细小的单元。这样，在xfig中可以通过选择层次来开启或关闭层次。
   * 多边形的深度为900级或900+ @p material_id, ，取决于标志 @p
   * level_depth. 相应地，边界边缘的深度为800级或800+ @p
   * boundary_id. 因此，边界边缘总是在单元的前面。
   * 在二维度为1的情况下不执行。
   *
   */
  template <int dim, int spacedim>
  void
  write_xfig(const Triangulation<dim, spacedim> &tria,
             std::ostream &                      out,
             const Mapping<dim, spacedim> *      mapping = nullptr) const;

  /**
   * 以SVG格式写入三角剖面图。    SVG（Scalable Vector
   * Graphics）是一种基于XML的矢量图像格式，由万维网联盟（W3C）开发和维护。这个函数符合2011年8月16日发布的最新规范SVG
   * 1.1。
   * 三角形的单元被写成多边形，在三角形的边界处有附加线。为了使单元格的某个属性可视化，例如它们的级别或材料ID，还可以对单元格进行着色。可以画一个色条来编码所选择的颜色。
   * 此外，还可以添加一个单元格标签，显示级别索引等。事实上，通过使用set_flags()和一个适当生成的
   * GridOutFlags::Svg,
   * 类型的对象，这个函数的可视化方式和内容的许多方面可以被定制。
   * @note
   * 这个函数目前只在两个空间维度的二维网格中实现。
   *
   */
  void
  write_svg(const Triangulation<2, 2> &tria, std::ostream &out) const;

  /**
   * 对所有其他维度和空间维度声明与上述相同的函数。这个函数目前没有实现，声明存在只是为了支持独立维度编程。
   *
   */
  template <int dim, int spacedim>
  void
  write_svg(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;


  /**
   * 以MathGL脚本格式编写三角剖面图。要解释这个文件，需要MathGL>=2.0.0的版本。
   * 为了在图形环境中掌握所产生的MathGL脚本，需要一个解释器。建议从
   * <code>mglview</code>, which is bundled with MathGL. <code>mglview</code>
   * 开始，它可以在图形窗口中解释和显示中小型MathGL脚本，并能转换为其他格式，如EPS、PNG、JPG、SVG，以及查看/显示动画。还可以进行一些小的编辑，如修改照明或alpha通道。
   * @note 未对二维码一的情况进行实现。
   *
   */
  template <int dim, int spacedim>
  void
  write_mathgl(const Triangulation<dim, spacedim> &tria,
               std::ostream &                      out) const;

  /**
   * 以VTK格式写入三角图。该函数写入一个UNSTRUCTURED_GRID文件，该文件包含以下VTK单元类型。
   * VTK_HEXAHEDRON, VTK_QUAD, 和 VTK_LINE, 取决于模板尺寸。
   * 在三个维度中，该函数写入一个文件，其中包括
   *
   *
   *
   *
   * - VTK_HEXAHEDRON单元格类型，包含三角测量的单元格信息
   *
   *
   *
   *
   * - VTK_QUAD单元类型，包含所有具有非零边界ID的边界面，以及所有具有非平面流形ID的面
   *
   *
   *
   *
   *
   * - VTK_LINE单元格类型，包含所有边界id为非零的边界边，以及所有边界id为非平坦的流形边 在二维方面。
   *
   *
   *
   * - VTK_QUAD单元格类型，包含Triangulation的单元格信息
   *
   *
   *
   *
   *
   * - VTK_LINE单元格类型，包含所有具有非零边界ID的边界面，以及所有具有非平面流形ID的面 在一维中
   *
   *
   *
   * - VTK_LINE单元类型，包含三角测量的单元信息 输出文件将包含两个CELL_DATA部分，`MaterialID`和`ManifoldID`，记录每个VTK单元类型的材料或边界ID，以及流形。请参阅[VTK文件格式](http://www.vtk.org/VTK/img/file-formats.pdf)文件，了解生成输出的解释。    可以使用配套的 GridIn::read_vtk 函数来读取用这种方法生成的VTK文件。
   *
   */
  template <int dim, int spacedim>
  void
  write_vtk(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;

  /**
   * 以VTU格式写入三角图。
   * 由于这个函数将数据写入输出流的方式，产生的输出文件对应于网格的忠实表示，即所有的单元都是可见的，以便进行可视化。一般来说，数据的格式不允许通过GridIn类再次读入这个文件。这是因为网格的每一个顶点都被复制了多次，就像相邻的单元一样。换句话说，每个单元都有自己独立的顶点集合，这些顶点与其他单元的顶点在同一位置，但被单独编号。
   * 为了创建一个可以用GridIn类读取的文件，必须将标志
   * GridOutFlags::Vtu::serialize_triangulation
   * 设置为真。在这种情况下，生成的vtu文件将在xml部分包含三角图，一般的vtu阅读器会忽略它。
   *
   */
  template <int dim, int spacedim>
  void
  write_vtu(const Triangulation<dim, spacedim> &tria, std::ostream &out) const;

  /**
   为每个处理器编写VTU格式的三角图，并为VisIt或Paraview的可视化添加一个.pvtu文件，将VTU文件的集合描述为同一模拟的所有部分。输出的形式是<tt>filename_without_extension.proc000*.vtu</tt>，其中是0,1,...,n_proc-1和<tt>filename_without_extension.pvtu</tt>。输入<tt>view_levels</tt>可以设置为true，以查看多层次方法的每一个层次。输入<tt>include_artificial</tt>可以设置为true，以查看每个处理器的人工单元。每个.vtu和.pvtu文件都会有subdomain、level_subdomain、level和proc_writing等属性。level值可用于将图像分离成多级方法中每一级的网格视图，proc_writing值可用于将图像分离成每个处理器的自有单元和幽灵单元。这是通过对每个值应用paraview中的 "标量翘曲 "过滤器来完成的。在打开输入<tt>view_levels</tt>设置为true的网格的.pvtu文件后，选择 "warp by scalar "过滤器。在 "Scalars "输入中选择<tt>proc_writing</tt>，在 "Normal "输入中输入1 0 0，然后点击应用。接下来再次选择 "标量翘曲 "过滤器。在 "Scalars "输入中选择<tt>level</tt>，在 "Normal "输入中输入0 1 0，然后点击Apply。这将给你提供以下图片。    @image html write_mesh_vtu_levels.png  如果<tt>view_levels</tt>保持为false，从而只给出活动层的网格，就足以将图像分成由不同处理器写入的视图。这在下面的图片中显示，其中<tt>include_artificial</tt>输入被设置为true。    @image html write_mesh_vtu_active.png 注意：根据网格的大小，你可能需要增加 "Scale Factor "的输入，以使每一块不重叠。
   *
   */
  template <int dim, int spacedim>
  void
  write_mesh_per_processor_as_vtu(const Triangulation<dim, spacedim> &tria,
                                  const std::string &filename_without_extension,
                                  const bool         view_levels = false,
                                  const bool include_artificial  = false) const;

  /**
   * 根据给定的数据格式，将网格写到 @p out
   * 。这个函数只是调用相应的<tt>write_*</tt>函数。
   *
   */
  template <int dim, int spacedim>
  void
  write(const Triangulation<dim, spacedim> &tria,
        std::ostream &                      out,
        const OutputFormat                  output_format,
        const Mapping<dim, spacedim> *      mapping = nullptr) const;

  /**
   * 以ParameterHandler设置的默认格式写入网格。
   *
   */
  template <int dim, int spacedim>
  void
  write(const Triangulation<dim, spacedim> &tria,
        std::ostream &                      out,
        const Mapping<dim, spacedim> *      mapping = nullptr) const;

  /**
   * 设置DX输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::DX &flags);

  /**
   * 设置%Gmsh输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::Msh &flags);

  /**
   * 设置UCD输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::Ucd &flags);

  /**
   * 设置GNUPLOT输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::Gnuplot &flags);

  /**
   * 为EPS输出的一维三角图设置标志
   *
   */
  void
  set_flags(const GridOutFlags::Eps<1> &flags);

  /**
   * 为EPS输出的二维三角图设置标志
   *
   */
  void
  set_flags(const GridOutFlags::Eps<2> &flags);

  /**
   * 为EPS输出的三维三角图设置标志
   *
   */
  void
  set_flags(const GridOutFlags::Eps<3> &flags);

  /**
   * 为EPS输出的三维三角图设置标志
   *
   */
  void
  set_flags(const GridOutFlags::XFig &flags);

  /**
   * 为SVG输出设置标志
   *
   */
  void
  set_flags(const GridOutFlags::Svg &flags);

  /**
   * 设置MathGL输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::MathGL &flags);

  /**
   * 设置VTK输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::Vtk &flags);

  /**
   * 设置VTU输出的标志
   *
   */
  void
  set_flags(const GridOutFlags::Vtu &flags);

  /**
   * 提供一个函数，可以告诉我们一个给定的输出格式通常有哪个后缀。例如，它定义了以下的映射关系。    <ul>   <li>   @p OpenDX:  <tt>.dx</tt>  <li>   @p gnuplot:  <tt>.gnuplot</tt>  <li>   @p ucd:  <tt>.inp</tt>  <li>   @p eps:  <tt>.eps</tt>。    </ul>  对所有实现的格式都提供了类似的映射。    由于这个函数不需要这个对象的数据，它是静态的，因此可以在不创建这个类的对象的情况下调用。
   *
   */
  static std::string
  default_suffix(const OutputFormat output_format);

  /**
   * 通过ParameterHandler选择的默认输出格式的默认后缀。
   *
   */
  std::string
  default_suffix() const;

  /**
   * 返回对应于给定字符串的 @p OutputFormat
   * 值。如果该字符串与任何已知的格式不匹配，就会抛出一个异常。
   * 由于这个函数不需要这个对象的数据，它是静态的，因此可以在不创建这个类的对象的情况下调用。它的主要目的是允许程序使用任何已实现的输出格式，而不需要在每次实现新格式时扩展程序的分析器。
   * 要获得目前可用的格式名称的列表，例如，将其交给ParameterHandler类，请使用函数get_output_format_names()。
   *
   */
  static OutputFormat
  parse_output_format(const std::string &format_name);

  /**
   * 返回一个已实现的输出格式的列表。不同的名称由垂直条形符号（<tt>`|'</tt>）分开，如ParameterHandler类所使用的。
   *
   */
  static std::string
  get_output_format_names();

  /**
   * 在ParameterHandler中声明参数。
   *
   */
  static void
  declare_parameters(ParameterHandler &param);

  /**
   * 解析ParameterHandler的参数。
   *
   */
  void
  parse_parameters(ParameterHandler &param);

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况
   *
   */
  DeclException0(ExcInvalidState);

private:
  /**
   * 默认的输出格式，由ParameterHandler设置。
   *
   */
  OutputFormat default_format;

  /**
   * OpenDX输出的标志。
   *
   */
  GridOutFlags::DX dx_flags;

  /**
   * 用于%Gmsh输出的标志。可以通过使用set_flags(const
   * GridOutFlags::Msh&)  函数来改变。
   *
   */
  GridOutFlags::Msh msh_flags;

  /**
   * 用于UCD输出的标志。可以通过使用set_flags(const
   * GridOutFlags::Ucd&) 函数来改变。
   *
   */
  GridOutFlags::Ucd ucd_flags;

  /**
   * 输出GNUPLOT数据时使用的标志。可以通过使用set_flags(const
   * GridOutFlags::Gnuplot&) 函数来改变。
   *
   */
  GridOutFlags::Gnuplot gnuplot_flags;

  /**
   * 在输出一个空间维度的EPS数据时使用的标志。可以通过使用set_flags(const
   * GridOutFlags::Eps<1>&) 函数来改变。
   *
   */
  GridOutFlags::Eps<1> eps_flags_1;

  /**
   * 在两个空间维度上输出EPS数据时使用的标志。可以通过使用
   * @p set_flags 函数来改变。
   *
   */
  GridOutFlags::Eps<2> eps_flags_2;

  /**
   * 在三维空间中输出EPS数据时使用的标志。可以通过使用
   * @p set_flags 函数来改变。
   *
   */
  GridOutFlags::Eps<3> eps_flags_3;

  /**
   * 用于XFig输出的标志。
   *
   */
  GridOutFlags::XFig xfig_flags;

  /**
   * 用于Svg输出的标志。
   *
   */
  GridOutFlags::Svg svg_flags;

  /**
   * 用于MathGL输出的标志。
   *
   */
  GridOutFlags::MathGL mathgl_flags;

  /**
   * 用于VTK输出的标志。
   *
   */
  GridOutFlags::Vtk vtk_flags;

  /**
   * 用于VTU输出的标志。
   *
   */
  GridOutFlags::Vtu vtu_flags;

  /**
   * 将面的网格信息写到  @p out.
   * 只有那些在边界上的面和边界指示器不等于0的面才会被打印出来，因为后者是边界面的默认值。
   * 由于在%Gmsh中，几何元素是连续编号的，这个函数需要一个参数
   * @p next_element_index
   * 提供下一个几何元素的编号。这个索引的数值应该比之前用来写入几何元素的索引多一个
   * @p out.   @returns  下一个不用的几何元素索引。      @warning
   * @p next_element_index
   * 应该（至少）比当前已经写入的三角形元素（线、单元、面）的数量多一个
   * @p out.
   * %如果有重复的索引，Gmsh将不能正确加载保存的文件。
   * 不幸的是，这个函数不能包含在常规的 @p
   * write_msh函数中，因为它需要对<tt>dim==1</tt>的情况进行特殊处理，在这种情况下，面迭代器是<tt>void*</tt>的，缺乏被调用的成员函数。我们实际上不会调用这些函数，但编译器在编译<tt>dim==1</tt>的函数时还是会抱怨。运气不好。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  write_msh_faces(const Triangulation<dim, spacedim> &tria,
                  const unsigned int                  next_element_index,
                  std::ostream &                      out) const;

  /**
   * 声明上述函数对1d的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_faces(const Triangulation<1, 1> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对1d, 2sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_faces(const Triangulation<1, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  unsigned int
  write_msh_faces(const Triangulation<1, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;



  /**
   * 将线的网格信息写到 @p out.
   * 中，只打印那些在边界上的线和边界指示器不等于0的线，因为后者是边界面的默认值。
   * 由于在%Gmsh中，几何元素是连续编号的，这个函数需要一个参数
   * @p next_element_index
   * 提供下一个几何元素编号。这个索引的数值应该比之前用来写入几何元素的索引多一个
   * @p out.   @returns  下一个未使用的几何元素索引。
   * @warning   @p next_element_index
   * 应该（至少）比当前已经写入的三角形元素（线、单元、面）的数量多一个
   * @p out.
   * %如果有重复的索引，Gmsh将不能正确加载保存的文件。
   * 不幸的是，这个函数不能包含在常规的 @p
   * write_msh函数中，因为它需要对<tt>dim==1</tt>和<tt>dim==2</tt>的情况进行特殊处理，在这种情况下，边缘迭代器是<tt>void*</tt>的，缺乏被调用的成员函数。我们实际上不会调用这些函数，但编译器在编译<tt>dim==1/2</tt>的函数时还是会抱怨。运气不好。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  write_msh_lines(const Triangulation<dim, spacedim> &tria,
                  const unsigned int                  next_element_index,
                  std::ostream &                      out) const;

  /**
   * 声明上述函数对1d的特殊化。什么也没做。
   *
   */
  unsigned int
  write_msh_lines(const Triangulation<1, 1> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;

  /**
   * 声明上述函数对1d, 2sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_lines(const Triangulation<1, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;

  /**
   * 声明上述函数对1d, 3sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_lines(const Triangulation<1, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对2d的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_lines(const Triangulation<2, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对2d, 3sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_msh_lines(const Triangulation<2, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;

  /**
   * 将面的网格信息写到 @p out.
   * 中，只打印那些在边界上的面和边界指示符不等于0的面，因为后者是边界面的默认值。
   * 由于（在UCD格式中）几何元素是连续编号的，这个函数需要一个参数
   * @p next_element_index
   * 提供下一个几何元素的编号。这个索引的数值应该比之前用来写几何元素的索引多一个，以
   * @p out.   @returns  下一个未使用的几何元素索引。
   * @warning   @p next_element_index
   * 应该（至少）比当前已经写入的三角测量元素（线、单元、面）的数量多一个
   * @p out.
   * 如果有重复的索引，可视化程序可能无法正确加载保存的文件。
   * 不幸的是，这个函数不能包含在常规的 @p
   * write_ucd函数中，因为它需要对<tt>dim==1</tt>的情况进行特殊处理，在这种情况下，面迭代器是<tt>void*</tt>的，缺乏被调用的成员函数。我们实际上不会调用这些函数，但编译器在编译<tt>dim==1</tt>的函数时还是会抱怨。运气不好。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  write_ucd_faces(const Triangulation<dim, spacedim> &tria,
                  const unsigned int                  next_element_index,
                  std::ostream &                      out) const;

  /**
   * 声明上述函数对1d的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_faces(const Triangulation<1, 1> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;

  /**
   * 声明上述函数对1d, 2sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_faces(const Triangulation<1, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  unsigned int
  write_ucd_faces(const Triangulation<1, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;


  /**
   * 将线的网格信息写到 @p out.
   * 中，只打印那些在边界上的线和边界指示符不等于0的线，因为后者是边界线的默认值。
   * 由于（在UCD格式中）几何元素是连续编号的，这个函数需要一个参数
   * @p next_element_index
   * 提供下一个几何元素的编号。这个索引的数值应该比之前用来写几何元素的索引多一个，以
   * @p out.   @returns  下一个未使用的几何元素索引。
   * @warning   @p next_element_index
   * 应该（至少）比当前已经写入的三角测量元素（线、单元、面）的数量多一个
   * @p out.
   * 如果有重复的索引，可视化程序可能无法正确加载保存的文件。
   * 不幸的是，这个函数不能包含在常规的 @p
   * write_ucd函数中，因为它需要对<tt>dim==1/2</tt>的情况进行特殊处理，在这种情况下，边缘迭代器是<tt>void*</tt>的，缺乏被调用的成员函数。我们实际上不会调用这些函数，但编译器在编译<tt>dim==1/2</tt>的函数时还是会抱怨。运气不好。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  write_ucd_lines(const Triangulation<dim, spacedim> &tria,
                  const unsigned int                  next_element_index,
                  std::ostream &                      out) const;

  /**
   * 声明上述函数对1d的特殊化。什么也没做。
   *
   */
  unsigned int
  write_ucd_lines(const Triangulation<1, 1> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对1d, 2sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_lines(const Triangulation<1, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对1d, 3sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_lines(const Triangulation<1, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;


  /**
   * 声明上述函数对2d的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_lines(const Triangulation<2, 2> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;
  /**
   * 声明上述函数对2d, 3sd的特殊化。什么也不做。
   *
   */
  unsigned int
  write_ucd_lines(const Triangulation<2, 3> &tria,
                  const unsigned int         next_element_index,
                  std::ostream &             out) const;


  /**
   * 返回三角形中边界指标不等于0的面的数量。只有这些面在<tt>write_*</tt>函数中被明确打印出来；所有指标为
   * numbers::internal_face_boundary_id
   * 的面是内部面，指标为零的边界面被认为是默认的。
   * 这个函数在一个维度上总是返回一个空列表。
   * 这个函数的原因与write_ucd_faces()的原因相同。更多信息见那里。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  n_boundary_faces(const Triangulation<dim, spacedim> &tria) const;

  /**
   * 声明上述函数对1d的特殊化。简单地返回0。
   *
   */
  unsigned int
  n_boundary_faces(const Triangulation<1, 1> &tria) const;

  /**
   * 声明上述函数对1d, 2sd的特殊化。只需返回0。
   *
   */
  unsigned int
  n_boundary_faces(const Triangulation<1, 2> &tria) const;

  /**
   * 声明上述函数对1d, 3sd的特殊化。只需返回0。
   *
   */
  unsigned int
  n_boundary_faces(const Triangulation<1, 3> &tria) const;

  /**
   * 返回三角形中边界指标不等于0的线的数量。只有这些线在<tt>write_*</tt>函数中被明确打印出来；所有指标为
   * numbers::internal_face_boundary_id
   * 的线都是内部线，边界上的面的指标值为0，被视为默认。
   * 这个函数在一维和二维中总是返回一个空列表。
   * 这个函数的原因与write_ucd_faces()的原因相同。更多信息见那里。
   *
   */
  template <int dim, int spacedim>
  unsigned int
  n_boundary_lines(const Triangulation<dim, spacedim> &tria) const;

  /**
   * 声明上述函数对1d的特殊化。简单地返回0。
   *
   */
  unsigned int
  n_boundary_lines(const Triangulation<1, 1> &tria) const;

  /**
   * 声明上述函数对1d, 2sd的特殊化。只需返回0。
   *
   */
  unsigned int
  n_boundary_lines(const Triangulation<1, 2> &tria) const;

  /**
   * 声明上述函数对1d, 3sd的特殊化。只需返回0。
   *
   */
  unsigned int
  n_boundary_lines(const Triangulation<1, 3> &tria) const;

  /**
   * 声明上述函数对2d的特殊化。简单地返回0。
   *
   */
  unsigned int
  n_boundary_lines(const Triangulation<2, 2> &tria) const;
  /**
   * 声明上述函数对2d, 3sd的特殊化。只需返回0。
   *
   */
  unsigned int
  n_boundary_lines(const Triangulation<2, 3> &tria) const;
};



DEAL_II_NAMESPACE_CLOSE

#endif


