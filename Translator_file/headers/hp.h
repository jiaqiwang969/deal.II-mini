//include/deal.II-translator/A-headers/hp_0.txt
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
 *    @defgroup hp hp-finite element support
 * 与hp-finite元素有关的类和函数。 step-27
 * 教程中概述了如何使用该命名空间中的类。 step-46
 * 中给出了一个稍显奇特的应用。
 * hp-命名空间实现了在deal.II中用于hp-framework的算法和数据结构。 @ref hp_paper  "hp-paper "
 * 中给出了关于这些算法如何工作以及使用何种数据结构的细节概述。
 *
 */

/**
 *    @defgroup hpcollection hp-Collections 
 * 在hp-finite element方法的实现中，每个单元可能有不同的有限元与之关联。为了处理这个问题，
 * hp::DoFHandler
 * 必须有一整套与之相关的有限元类。这个概念由
 * hp::FECollection
 * 类代表。这种类型的对象作为容器，容纳了一整套的有限元对象。我们不在每个单元上存储指向有限元对象的指针，而只为每个单元存储一个索引，该索引标识了该单元应该使用的集合中的有限元对象。与给定单元相关的DoFHandler对象可以根据单元使用的有限元为每个单元分配自由度。
 * 在单元上积分项时也会出现类似的情况：人们可能希望对不同的有限元使用不同的正交公式。例如，在我们使用Q1元素的单元上，QGauss(2)对象（即每个空间方向上有两个点的正交公式）可能就足够了，但在另一个使用Q3元素的单元上，这将导致积分不足，我们应该使用QGauss(4)公式。就像上面一样，存在一个类
 * hp::QCollection ，作为正交公式的集合。
 * 最后，人们可能希望对具有不同阶数的单元的边界逼近使用不同阶数的有限元。
 * hp::MappingCollection 类允许这样做。 所有这三个类，
 * hp::FECollection,   hp::QCollection,  和  hp::MappingCollection
 * 类，实现了与  <code>std::vector</code>
 * 非常相似的接口。它们有函数  <code>push_back()</code>
 * 来添加有限元、正交公式或映射到集合。他们有一个
 * <code>operator[] (unsigned int)</code>
 * 函数，允许检索对集合中某个给定元素的引用。他们还有一个
 * <code>size()</code>
 * 函数可以返回集合中元素的数量。一些类，特别是持有有限元对象的类，也实现了其他特定的功能。
 * 相似性超出了接口的范围。当向集合中添加一个元素时，所有的类都会创建一个参数的副本。这样就可以把一个临时对象传递给添加元素的函数。例如，下面的工作。
 * @verbatim
 * FECollection<dim> fe_collection;
 * for (unsigned int degree=1; degree<5; ++degree)
 *   fe_collection.push_back (FE_Q<dim>(degree));
 * @endverbatim
 * 这样一来，人们就可以把多项式度数为1到4的元素添加到集合中。没有必要保留所添加的对象：集合会对其进行复制，它不仅存储了一个指向给定有限元对象的指针。这个观察同样适用于其他集合类。
 * 习惯上，在一个hp-finite
 * element程序中，人们保留了具有相同数量元素的有限元和正交公式集合，一个集合中的每个元素都与另一个集合中的元素相匹配。这不是必须的，但它常常使编码变得简单得多。如果使用映射的集合，对
 * hp::MappingCollection 对象也是如此。 每当在hp-finite
 * element程序中考虑p-adaptivity时，需要建立一个有限元的层次结构，以确定细化的后续有限元和粗化的前面的有限元。通常，这种层次结构考虑了有限元空间的嵌套方式：例如，
 * $Q_1$ 元素描述了 $Q_2$ 元素的一个子空间，因此进行 $p$
 * 细化通常意味着使用更大（更精确）的有限元空间。换句话说，有限元的层次结构是通过考虑集合中的一些元素是其他元素的子空间还是超空间来建立的。
 * 默认情况下，我们假设有限元是根据其多项式程度以升序存储的。如果元素的顺序不同，需要通过
 * hp::FECollection::set_hierarchy()
 * 成员函数向集合提供相应的层次结构。
 *
 * @ingroup hp
 *
 */


/**
 * 一个用于实现hp-finite元素特定算法和数据结构的命名空间。
 *
 * @ingroup hp
 *
 */
namespace hp
{
}


