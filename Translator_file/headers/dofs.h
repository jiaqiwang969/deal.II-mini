//include/deal.II-translator/A-headers/dofs_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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
 *    @defgroup dofs Degrees of Freedom
 * 本模块将与处理自由度有关的类和命名空间分组。该组的中心类是DoFHandler类：它建立在三角形和有限元类之上，并根据有限元对象所描述的有限元空间的要求在三角形的每个单元上分配自由度。DoFHandler类还有其他变体，如
 * hp::DoFHandler ，对更特殊的情况做类似的事情。
 * DoFHandler对象与FiniteElement（或 hp::FECollection 中的
 * hp::DoFHandler)
 * ）类型的对象一起使用，以列举该特定有限元的三角结构中存在的所有自由度。因此，网格、有限元和DoF处理程序对象的组合可以被认为是提供了一个<i>basis</i>的有限元空间：网格提供了定义基函数的位置；有限元描述了存在哪些种类的基函数；DoF处理程序对象提供了基的枚举，也就是说，它提供了空间的具体结构，因此我们可以通过系数向量来描述这个有限维空间的函数。
 * DoFHandlers扩展了Triangulation对象（以及 @ref
 * 网格模块中的其他类），因为它们也提供了迭代器，在所有单元格、面或其他构成三角形的几何对象上运行。这些迭代器是从三角形迭代器派生出来的，因此提供了同样的功能，但它们也提供了额外的功能。例如，它们允许查询与当前单元相关的自由度的索引。请注意，DoFHandler类来自Triangulation<i>not
 * derived</i>，尽管它们使用Triangulation对象；原因是可以有一个以上的DoFHandler对象对同一个Triangulation对象工作。
 * 除了DoF处理程序类之外，这个模块还拥有一些在应用程序中不常用的辅助类，以及三个与DoFHandler类的数据结构没有直接联系的类。其中第一个是AffineConstraints类，用于存储和处理与悬挂节点相关的约束。其次，DoFRenumbering命名空间提供了可以重新排序自由度的函数；在它的函数中，有在下游方向排序自由度的函数，例如，有以使相关矩阵的带宽最小化的方式排序自由度。最后，DoFTools命名空间提供了各种处理自由度的算法。
 * 从总体上看，这个模块的各个部分与库的其他部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall"]; mapping
 * [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature"]; dh
 * [label="DoFHandler",       URL="\ref dofs", fillcolor="deepskyblue"];
 * fevalues   [label="FEValues",         URL="\ref feaccess"]; systems
 * [label="Linear systems",   URL="\ref LAC"]; solvers    [label="Linear
 * solvers",   URL="\ref Solvers"]; output     [label="Graphical output",
 * URL="\ref output"]; manifold   [label="Manifold",         URL="\ref
 * manifold"];
 *
 * { rank=same mapping
 *
 * -> quadrature [dir="none", color="transparent"]; quadrature
 *
 * -> fe      [dir="none", color="transparent"]; fe
 *
 * -> tria            [dir="none", color="transparent"]; }
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
 * -> mapping     [color="black",style="solid"]; } @enddot
 *
 */


