//include/deal.II-translator/A-headers/geometry_and_primitives_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2021 by the deal.II authors
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
 *    @defgroup geomprimitives Geometric and other primitives
 * 本组包含一些作为几何基元或其他数学对象的基元的类。例如，Tensor
 * @<rank,dim@>  类提供了等级为 <code>rank</code> in <code>dim</code>
 * 的空间维度的张量。同样地，SymmetricTensor提供了对称的张量。
 * 在几何学上，点类是deal.II库中所有几何描述的基础。它表示
 * <code>dim</code>
 * 维空间中的一个几何点。我们可以把一个点看作是一个坐标为
 * <code>dim</code>
 * 的矢量，它连接着原点和那个特定的点；因此，点类是从秩1的张量（即矢量）派生出来的，但与任意张量相比，点具有空间中点的特殊内涵，因此具有一些额外的属性。
 * 在deal.II中，网格是由线段、四边形或六面体（取决于空间维度）构建的。GeometryInfo类用于描述这些基本对象在单位空间中的属性（即对于单位线、单位方和单位立方）。它提供了静态数据成员，表示每个单元的顶点数量、每个面的线，或者顶点的位置。这种抽象允许编写的应用程序大多独立于实际的空间维度：所有顶点的循环只是从0到
 * GeometryInfo<dim>::vertices_per_cell
 * ，而不是从0到4（在2D）或0到8（在3D）。这样一来，程序在2D和3D中都是正确的，人们可以通过重新编译在不同的空间维度上运行程序，而不必改变代码的很大一部分。这些与维度无关的编程技术在前几个教程程序中得到了广泛的讨论，并在整个交易中得到应用。
 *
 */


