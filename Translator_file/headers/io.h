//include/deal.II-translator/A-headers/io_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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


/**
 *    @defgroup IO Input/Output
 * 这个模块收集了用于读写网格和数据的类。有两个子模块分别用于这些操作。
 *
 */

/**
 *    @defgroup input Input
 * deal.II可以读取许多不同格式的网格。然而，它们都被限制在所谓的
 * "粗略的网格
 * "上，即没有细化层次，特别是没有悬挂节点的网格。GridIn类详细描述了支持哪些格式。
 * 此外，deal.II可以使用DataOutReader读取中间图形格式。这种格式被用作与模拟相关的数据之间的中间步骤，由DataOutBase类（或通过更多的派生类在
 * `ref输出模块'中描述）写入。DataOutReader类将这些数据读回来，然后它可以被转换为可视化程序支持的任何一种数据格式。
 * 最后，ParameterHandler和MultipleParameterLoop类（以及相关的Patterns命名空间）被用来处理描述程序运行时参数的参数文件，这些参数不想在程序源中硬编码。
 *
 *  <h3>The PathSearch class</h3>
 * PathSearch类是输入处理中的一个辅助类。它被用来在一个目录列表中寻找一个文件，就像unix系统在
 * <code>PATH</code>
 * 环境变量中列出的目录中寻找可执行文件一样。
 *
 * @ingroup IO
 *
 */

/**
 *    @defgroup output Graphical output
 * deal.II生成三种类型的输出：它可以将三角图/网格写成几种网格阅读器（包括deal.II本身的阅读器）可以理解的格式，它还可以创建用于数据可视化的输出。最后，它能以图形格式输出矩阵。
 *
 *  <h3>Visualization of data</h3>
 * deal.II通过DataOutBase类支持大量流行的可视化格式，如OpenDX、gmv或gnuplot程序所使用的格式。支持的格式的完整列表列在DataOutBase类的文档中。
 * DataOutBase类只负责在一些不同的可视化格式中实际写入一些中间格式。这个中间格式是由直接或间接从DataOutBase派生的类生成的。例如，DataOut类最常被用来从一个三角形、一个DoFHandler对象（将一个特定的有限元类与三角形联系起来）和一个或多个数据向量生成这种中间格式。DataOut类从每个单元创建中间数据，随后由DataOutBase类以某种最终格式写入。几乎所有的例子程序，从
 * step-3 开始，都使用这种方法来生成输出。
 * DataOutFaces类是另一种从模拟数据创建中间格式的方法。然而，它不是从三角形的每个单元创建可视化数据，而是只为位于表面上的单元的所有面创建信息（尽管该类有一种方法可以覆盖选择哪些面应该被生成输出）。虽然这在2D中可能不是特别有趣（这些面只是线段），但在3D中，如果人们真正想知道的是域的形状或曲面上的一个变量的值，这往往是有帮助的。使用DataOutFaces类可以节省生成和存储所有内部单元的数据的工作，这对大型三维模拟来说是非常昂贵的。
 * 第三个类，DataOutRotation类，允许采取一个二维模拟，并通过围绕给定的轴旋转二维领域来生成三维数据。这对于使用旋转对称性的模拟的可视化非常有用，例如，一个圆柱体。
 * 最后，DataOutStack类允许在时空域中对与时间有关的模拟数据进行可视化：它收集每个时间步骤的结果，并在最后将所有这些信息一次性输出为一个时空文件。
 *
 *  <h3>Grid output</h3>
 * Meshes，没有任何与之相关的数据向量，也可以用多种格式写入。这是通过GridOut类完成的，该类的文档中列出了支持的格式。
 * 一些教程程序，特别是 step-1 、 step-6 、 step-9 、 step-10 、
 * step-12  b和 step-14 演示了GridOut类的使用。
 *
 *  <h3>Matrix output</h3>
 * 通过MatrixOut类，deal.II还可以以彩色或天际线图的形式给出矩阵的图形可视化。MatrixOut类使用DataOutBase进行输出。因此，矩阵可以用后一类支持的所有格式进行可视化。
 *
 * @ingroup IO
 *
 */

/**
 *    @defgroup textoutput Textual output
 * 除了提供图形输出格式的类（见 @ref
 * 输出模块），deal.II还有一些类，以多种方式促进文本输出。它们被收集在本模块中。更多细节请参见这些类的文档。
 *
 * @ingroup IO
 *
 */


