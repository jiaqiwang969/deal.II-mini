//include/deal.II-translator/A-headers/quadrature_0.txt
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
 *    @defgroup Quadrature Quadrature formulas
 * 本模块包含基类正交以及由deal.II提供的正交公式。正交公式提供了两个基本数据：单元格[0,1]^d上的正交点的位置，以及每个正交点的权重。
 * 由于deal.II使用四边形和六面体，几乎所有的正交公式都是作为定义在单位区间[0,1]上的一维正交公式的张量产物产生的，这使得它们在高维情况下的定义几乎是微不足道的。然而，通过QAnisotropic类，该库也允许各向异性的张量产品（一个坐标方向上的正交点比另一个方向上的多），以及定义不是张量产品的正交公式。
 * 从总体上看，这个模块的类与库中的其他各种部分相互作用。
 * @dot digraph G { graph[rankdir="TB",bgcolor="transparent"];
 *
 * node [fontname="FreeSans",fontsize=15, shape=box,height=0.2,width=0.4,
 * color="black", fillcolor="white", style="filled"]; edge [color="black",
 * weight=10];
 *
 * tria       [label="Triangulation",    URL="\ref grid"]; fe
 * [label="Finite elements",    URL="\ref feall"]; mapping
 * [label="Mapping",          URL="\ref mapping"]; quadrature
 * [label="Quadrature",       URL="\ref Quadrature", fillcolor="deepskyblue"];
 * dh         [label="DoFHandler",       URL="\ref dofs"]; fevalues
 * [label="FEValues",         URL="\ref feaccess"]; systems    [label="Linear
 * systems",   URL="\ref LAC"]; solvers    [label="Linear solvers",
 * URL="\ref Solvers"]; output     [label="Graphical output", URL="\ref
 * output"]; manifold   [label="Manifold",         URL="\ref manifold"];
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
 * } @enddot   <h3>Use</h3>
 * 正交公式除其他用途外，还可用于积分矩阵条目和右手边矢量的分量。为此，定义在单元格上的正交点必须被映射到实数单元格上的相应位置，并且权重必须乘以雅各布的行列式。这一步是由派生自Mapping基类的类来完成的，尽管这通常是隐藏的，因为如果没有提供特定的映射，库的许多部分会退回到使用MappingQ1类型的对象。
 * 下一步是评估形状函数和它们在这些位置的梯度。虽然从FiniteElement基类派生出来的类提供了对单元格上形状函数的描述，但在正交点上的实际评估以及将其与从映射中获得的信息结合起来的工作是由FEValues类及其关联者完成的。因此，从本质上讲，FEValues类是对有限元空间（由FiniteElement类定义）的视图，在正交点（由正交类提供）上进行评估，并映射到真实空间（而不是单元空间）的单元内的位置（由Mapping类提供映射）。
 * FEValues类作为副产品，提供了映射到实数单元的正交点的位置，也可用于其他用途。例如，这可以用来在这些点上评估一个右手边的函数。
 *
 *  <h3>QIterated</h3>
 * QIterated类用于从现有的正交公式中构造一个迭代的正交公式，从而在不增加阶数的情况下提高公式的精度。例如，通过对点在0和1、权重为1/2和1/2的梯形规则进行两次迭代，我们可以得到一个点在0、1/2和1、权重分别为1/4、1/2和1/4的正交公式。这个公式是通过将正交公式分别投射到子区间[0,1/2]和[1/2,1]上，然后将左边区间的右端点与右边区间的左端点合并得到的。以同样的方式，所有的一维正交公式都可以被迭代。高维迭代公式是作为一维迭代公式的张量积产生的。
 *
 *  <h3>QAnisotropic</h3>
 * 高维的通常正交公式产生的张量产品在每个方向上都是相等的，而QAnisotropic类产生的张量产品在每个方向上可能是不同的公式。
 *
 *  <h3>QProjector</h3>
 * QProjector类本身实际上不是正交规则，但它提供了在高维单元表面计算正交公式的函数。
 * 本模块中的所有其他类实际上实现了不同顺序和其他特征的正交规则。
 *
 *  <h3>QuadratureSelector</h3>
 * 该类用于根据标识正交公式的字符串生成一个正交对象。这在希望在输入文件中指定某个正交公式，而不是在程序中硬编码的情况下很有用。
 *
 */


