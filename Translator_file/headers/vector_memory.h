//include/deal.II-translator/A-headers/vector_memory_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
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
 *    @defgroup VMemory Vector memory management
 * 这个模块将一些类分组，用来避免在迭代过程中反复分配和解配向量。这些方法都使用基类VectorMemory的一个对象来获取它们的辅助向量。
 * 关于这个话题的一些讨论可以在 step-20
 * 中关于InverseMatrix类的讨论中找到。
 *
 * @ingroup LAC
 *
 */


