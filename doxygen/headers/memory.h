//include/deal.II-translator/A-headers/memory_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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
 *    @defgroup memory Memory handling
 * 本组有一些用于内存处理的基本类和命名空间。Subscriptor和SmartPointer类用于计数内存处理，也就是说，每当SmartPointer被设置为指向一个对象时，它就会增加该对象中的一个计数器；当指针被设置为指向其他地方时，它又会减少。这样一来，人们总是知道一个对象还有多少用户。虽然这本身很少有用，但如果一个对象被销毁，而某个地方的指针仍然指向它，它就会被用来产生一个异常，因为在以后的时间里通过该指针的任何访问都会导致访问无效的内存区域。
 * 与此相反，MemoryConsumption命名空间提供的函数可以用来确定对象的内存消耗。对于一些简单的类，比如标准库的容器，它直接确定它们需要多少内存（或者至少给出一个估计值）。对于deal.II类，它使用大多数类有的
 * <code>memory_consumption</code> 成员函数。
 *
 * @ingroup utilities
 *
 */


