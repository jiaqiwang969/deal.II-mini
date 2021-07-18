//include/deal.II-translator/A-headers/simplex_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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
 *    @defgroup simplex Simplex support (experimental)
 * @brief This module describes the experimental simplex support in deal.II.
 * deal.II中的Simplex和混合网格仍然是实验性的，也就是正在进行的工作。该库的大部分内容已经被移植到能够操作这种网格上。然而，仍有许多函数需要被概括。通过查看
 * "test/simplex
 * "文件夹中的测试，你可以对移植的功能有一个很好的概述。在下文中，我们提供了两个非常基本的例子来入门，并提供一些实现细节。
 *
 *
 * @section  simplex_reference_example_simplex 示例：simplex网格
 * 下面的代码显示了如何处理单纯x网格。
*  @include step_3_simplex.cc
 *
 *
 * @section  simplex_reference_example_mixed 示例：混合网格
 * 下面的代码显示了如何处理混合网格的工作。
*  @include step_3_mixed.cc
 *
 *
 * @section  simplex_reference_cells 参考单元
 * 在二维中，我们提供三角形和四边形，在三维中的可能方向如下。
 * <div class="twocolumn" style="width: 100%"> <div class="parent"> <div
 * class="img" align="center">
       @image html reference_cells_0.png
 * </div>   <div class="text" align="center"> 2D: triangle and quadrilateral
 * </div>  </div>  <div class="parent"> <div class="img" align="center">
       @image html reference_cells_1.png
 * </div>   <div class="text" align="center"> Possible orientations of
 * triangles and quadrilaterals in 3D </div>  </div> </div>
 * 在三维中，四面体、金字塔、楔形和六面体都可以使用。
 * <div class="parent"> <div class="img" align="center">
       @image html reference_cells_2.png
 * </div>   <div class="text" align="center"> 3D: Tetrahedron </div>  </div>
 * </div> <div class="parent"> <div class="img" align="center">
       @image html reference_cells_3.png
 * </div>   <div class="text" align="center"> 3D: Pyramid </div>  </div <div
 * class="parent"> <div class="img" align="center">
       @image html reference_cells_4.png
 * </div>   <div class="text" align="center"> 3D: Wedge </div>  </div <div
 * class="parent"> <div class="img" align="center">
       @image html reference_cells_5.png
 * </div>   <div class="text" align="center"> 3D: Hexahedron </div>  </div>
 * </div>
 * 一个三维参考单元的每个表面由二维参考单元组成。枚举其顶点和线的编号的文件在右列给出。
 *
 *
 */


