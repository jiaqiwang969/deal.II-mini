//include/deal.II-translator/A-headers/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2014 by the deal.II authors
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
 *    @defgroup utilities Utility functions and classes
 * 这个模块简单地收集了一些函数和类，为那些通常与有限元程序没有太大关系，但恰好也需要的任务提供了通用工具。
 *
 */


/**
 *    @defgroup data Data storage primitives
 * 这里有几个简单的类，有助于存储和查看数据。例如，Table模板不仅允许使用对象的数组（为此可能要使用
 * std::vector
 * 类），而且还允许使用任意对象的二维（矩形）表，以及高阶类似物，直到具有（目前）七个索引的表。
 * 类似地，VectorSlice函数是一个基元，它接收任何具有类似于矢量的接口的东西（例如deal.II
 * Vector或 std::vector
 * 类），并对其呈现一个视图，就好像它本身是一个矢量。
 *
 * @ingroup utilities
 *
 */


