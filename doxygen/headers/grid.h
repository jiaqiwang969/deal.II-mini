//include/deal.II-translator/A-headers/grid_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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
 *    @defgroup grid Grids and Triangulations
 * 这个模块将与网格的拓扑结构和几何形状有关的函数和类归类。一个网格可以被认为是一个单元的集合；如果网格被细化了（可能是以一种自适应的方式），那么这个集合就会被分组为一个细化等级的层次结构。除了单元之外，构成三角形的几何对象是单元的面（在三维中是单元的边）以及单元的顶点。请注意，我们有些滥用<i>triangulation</i>这个词，因为deal.II只实现了由线性、四边形和六面体单元组成的三角形，三角形和四面体不被支持。
 * 这个单元的集合由Triangulation类和派生类（如 parallel::distributed::Triangulation
 * 和 parallel::shared::Triangulation.
 * ）管理，它在内存中保存相关的数据，并提供接口来查询。你想在单元格上做的大多数事情都是在所有单元格上循环进行的。为此，Triangulation类提供了迭代器的概念（见
 * @ref Iterators
 * ）：虽然实现方式不同，但它们的行为类似于单元格或面的指针，可以查询到单元格的几何属性以及相邻单元格或单元格的面等信息。
 * 值得注意的是，Triangulation类只存储几何（即顶点和单元的位置）和网格的拓扑结构（即哪些单元是其他单元的邻居，等等）。它与网格上可能定义的有限元或自由度没有关系。这些功能由DoFHandler类（见 @ref
 * dofs
 * 模块）执行，该类获得有限元空间的描述，并分配和管理顶点、面或单元的自由度，如有限元类所描述的。这种分离使得多个DoFHandler类可以同时在同一个网格上工作。
 * 在整个计划中，deal.II中的三角计算与库中的各种其他部分交互。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid",
 * fillcolor="deepskyblue"]; fe         [label="Finite elements",    URL="\ref
 * feall"]; mapping    [label="Mapping",          URL="\ref mapping"];
 * quadrature [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs"]; fevalues   [label="FEValues",
 * URL="\ref feaccess"]; systems    [label="Linear systems",   URL="\ref
 * LAC"]; solvers    [label="Linear solvers",   URL="\ref Solvers"]; output
 * [label="Graphical output", URL="\ref output"]; manifold
 * [label="Manifold",         URL="\ref manifold"];
 *
 * tria
 *
 * -> dh              [color="black",style="solid"]; fe
 *
 * -> dh                [color="black",style="solid"]; fe
 *
 * -> fevalues          [color="black",style="solid"]; mapping
 *
 * -> fevalues     [color="black",style="solid"]; quadrature
 *
 * -> fevalues  [color="black",style="solid"]; dh
 *
 * -> systems           [color="black",style="solid"]; fevalues
 *
 * -> systems     [color="black",style="solid"]; systems
 *
 * -> solvers      [color="black",style="solid"]; solvers
 *
 * -> output       [color="black",style="solid"]; manifold
 *
 * -> tria        [color="black",style="solid"]; manifold
 *
 * -> mapping     [color="black",style="solid"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
 *
 * node [fontname="FreeSans",fontsize=12, shape=record,height=0.2,width=0.4,
 * color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
 * edge [color="gray55", weight=1];
 *
 * opencascade [label="OpenCASCADE"]; opencascade
 *
 * -> manifold [dir="none"];
 *
 *
 * node [fontname="FreeSans",fontsize=12, shape=ellipse,height=0.2,width=0.4,
 * color="gray55", fontcolor="gray55", fillcolor="white", style="filled"];
 * edge [color="gray55", weight=1];
 *
 * gmsh        [label="gmsh", URL="\ref Gmsh"]; gmsh
 *
 * -> tria       [dir="none"]; } @enddot <h3>Grid generation</h3>
 * 有三种方法来创建网格。   <ul>   <li>  由GridGenerator类创建；  <li>  从文件读取；  <li>  手工创建。   </ul>
 * 对于第一种情况，GridGenerator类提供的函数可以自动生成最简单和最常见的几何图形。例如，矩形（或砖形）几何体以及圆形、球形或圆柱形都可以用这个类中的函数生成。大多数的教程程序都使用这种机制。
 * 其次，可以使用GridIn类从输入文件中读入一些不同格式的网格。使用这个类，可以读取几万或十万个单元的网格，尽管这并不推荐：自适应有限元方法的威力只有在初始网格尽可能粗的情况下才能发挥出来，而且还有空间进行一系列的自适应细化步骤。如果初始网格已经太细了，那么在自适应网格细化能够发挥很大作用之前，就会耗尽内存或计算时间。尽管如此，GridIn类可以用于复杂的几何形状，或者用于与其他程序进行比较或交互，这些程序在网格上进行计算，然后通过这个类进行交换。
 * step-5 教程程序展示了如何使用GridIn类。
 * 第三种方式是手工创建网格，通过建立一个数据结构来描述三角形的顶点和单元。这种方法在中等复杂度的情况下非常有用，即无需借助网格生成器就可以手工建立网格，但该领域不属于GridIn类已经支持的领域。在这个方法中，所构建的数据结构被交给三角化类的create_triangulation()函数。
 * step-14 的教程程序展示了如何做到这一点。
 *
 *  <h3>Grid output</h3>
 * 网格可以被写入一些不同格式的输出文件中。如果这涉及到在这个网格上获得的仿真结果，那么这将使用DataOut类（在 @ref
 * output
 * 模块中有更详细的描述）来完成。另一方面，如果只需要将网格的几何结构和拓扑结构写入文件，GridOut类可以为你做到这一点。
 *
 *  <h3>Tool classes</h3>
 * GridTool类提供了各种作用于网格的功能。例如，这包括移动节点，拉伸或旋转整个三角形，计算域的直径，或将其细分为大致相同大小的块，以便进行并行计算。
 * GridRefinement类实现了一系列的网格细化算法，基于给定的细化指标给其成员函数。
 *
 *  <h3>Internal classes</h3>
 * 除了上述内容外，本模块中还有相当数量的类只用于网格处理的内部数据结构中。它们一般都在内部命名空间，而不是为了在应用程序代码中使用。
 *
 * @author  Wolfgang Bangerth，1998-2006年
 *
 */


