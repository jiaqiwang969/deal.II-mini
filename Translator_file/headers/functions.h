//include/deal.II-translator/A-headers/functions_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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
 *    @defgroup functions Functions
 * 函数在deal.II中被用于不同的地方，例如用来描述边界条件、方程中的系数、强制项或精确解。由于方程的封闭式表达式通常很难作为函数参数传递，deal.II使用函数基类来描述这些对象。基本上，这个基类的接口要求派生类实现返回一个或一列特定位置的函数值的能力，以及可能的（如果需要）函数的梯度或二次导数。有了这个，函数对象就可以被像
 * VectorTools::interpolate,  VectorTools::project_boundary_values,
 * 这样的算法和其他函数所使用。
 * 有些函数是反复需要的，因此已经在deal.II中提供。这包括一个具有常数值的函数；一个在任何地方都为零的函数，或者一个只有一个向量分量具有特定值而所有其他分量为零的向量值函数。在函数命名空间中还定义了一些更专门的函数。
 *
 *  <h3>Time dependent functions</h3>
 * 对于时间相关的计算，边界条件和/或右手边的函数也可能随时间变化。由于在一个给定的时间步长，人们通常只对函数的空间依赖性感兴趣，如果必须向所有使用函数对象的方法传递一个时间变量的值，那就很尴尬了。例如，
 * VectorTools::interpolate_boundary_values
 * 函数将不得不接受一个时间参数，当它想查询边界函数在特定时间步长的值时，它可以使用这个参数。然而，如果我们考虑的是一个静止的问题，它也必须这样做，因为在这个问题上没有类似时间变量的东西。
 * 为了规避这个问题，函数对象总是只被认为是空间函数。然而，Function类是由FunctionTime基类派生出来的，如果有必要的话，它可以存储一个时间变量的值。这样，人们可以定义一个作为空间函数的函数对象，但在内部可以通过引用一个特定的时间来实现。在上面的例子中，在把函数对象交给
 * VectorTools::interpolate_boundary_values
 * 方法之前，人们会把它的时间设置为现在的时间步长。
 *
 *  <h3>Tensor-valued functions</h3>
 * 函数类是最常用的，但有时人们需要一个函数，其值是张量，而不是标量。TensorFunction模板可以为你做到这一点。除了返回类型外，该接口与函数类的接口基本相同。
 *
 */


