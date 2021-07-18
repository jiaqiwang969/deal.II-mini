//include/deal.II-translator/fe/mapping_q1_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_q1_h
#define dealii_mapping_q1_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mapping */ 
 /*@{*/ 


/**
 * 实现从参考单元到一般四边形/六面体的 $d$  -线性映射。
 * 该类实现的映射将参考（单位）单元映射到具有 $d$
 * 维度直线的一般网格单元。然而，请注意，在三维中，一般的、三线映射的单元格的<i>faces</i>可能是弯曲的，即使边缘不是）。这是用于多面体域的标准映射。它也是整个deal.II中用于许多函数的映射，这些函数有两种变体，一种允许明确传递映射参数，另一种则简单地返回到这里声明的MappingQ1类。(或者，事实上，回到一个MappingQGeneric(1)类型的对象，它完全实现了这个类的功能。)
 * 这个映射的形状函数与多项式1度的有限元FE_Q相同。
 *
 *
 * @note
 * 这个类实际上不过是调用MappingQGeneric的不同名称，参数为多项式度数为1。
 *
 *
 */
template <int dim, int spacedim = dim>
class MappingQ1 : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * 默认构造函数。
   *
   */
  MappingQ1();

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;
};



/**
 * 库中的许多地方默认使用（双，三）线性映射，除非用户明确提供不同的映射来使用。在这些情况下，被调用的函数必须创建一个
 * $Q_1$
 * 映射对象，即一个MappingQGeneric(1)类型的对象。这是很昂贵的。在受影响的函数中创建这样的对象作为静态对象也是很昂贵的，因为静态对象在程序的整个生命周期中都不会被销毁，尽管它们只需要在代码第一次运行某个特定函数时被创建一次。
 * 为了避免在整个库的这些上下文中创建（静态或动态）
 * $Q_1$ 映射对象，这个类定义了一个静态 $Q_1$
 * 映射对象。然后这个对象可以在所有需要这种对象的地方使用。
 *
 *
 * @note
 * 应该避免使用这个对象，因为它只适用于网格完全由四边形或六面体组成的情况。
 * 请使用
 * `ReferenceCells::get_hypercube<dim>().get_default_linear_mapping()`
 * 代替。
 *
 *
 */
template <int dim, int spacedim = dim>
struct StaticMappingQ1
{
  /**
   * 本类文档中讨论的静态 $Q_1$ 映射对象。
   *
   */
  static MappingQGeneric<dim, spacedim> mapping;
};


 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


