//include/deal.II-translator/A-headers/integrators_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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
 *    @defgroup Integrators Integrators
 * 一个命名空间和函数的集合，简化了有限元空间上的表格和双线性表格的编码。这里收集了两个不同目的的函数：MeshWorker中对有限元网格的抽象积分，以及LocalIntegrators中对具体问题的单元和面项的积分的实际实现。
 *
 *
 * @note
 * 关于编码惯例、类之间的关系以及实现细节的文档可以在本模块的命名空间文档中找到。
 * <h3>Integration on finite element meshes</h3>
 * 当我们在有限元空间上积分一个函数或一个函数时，积分循环的结构总是相同的。我们有3到5个嵌套循环，从外到内。   <ol>   <li>  在所有单元上循环  <li>  可选择在所有面上循环以计算通量  <li>  在单元/面的所有正交点上循环  <li>  可选择在所有试验函数上循环以计算形式  <li>  可选择在所有试验函数上循环以计算双线性形式  </ol>
 * 这些循环自然分为两类，即计算单元和面的贡献（循环3到5），以及对网格对象的外循环，通常被称为
 * <em>  装配  </em>  。
 * 在deal.II中对外循环的支持可以在命名空间MeshWorker中找到（见那里的文档）。为了支持单元和面的贡献（从现在开始称为局部贡献），deal.II提供了FEValuesBase和其派生类。虽然外循环是通用的（数据类型除外），但局部贡献的计算却取决于问题。因此，这里不可能有通用算法。不过，我们可以为此目的定义一个函数的通用接口，并提供一个本地积分器库供应用中使用。这些集合在命名空间LocalIntegrators中。
 *
 */


