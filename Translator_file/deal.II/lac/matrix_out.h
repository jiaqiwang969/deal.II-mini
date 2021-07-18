//include/deal.II-translator/lac/matrix_out_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2021 by the deal.II authors
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

#ifndef dealii_matrix_out_h
#  define dealii_matrix_out_h

#  include <deal.II/base/config.h>

#  include <deal.II/base/data_out_base.h>

#  include <deal.II/lac/block_sparse_matrix.h>
#  include <deal.II/lac/sparse_matrix.h>

#  ifdef DEAL_II_WITH_TRILINOS
#    include <deal.II/lac/trilinos_block_sparse_matrix.h>
#    include <deal.II/lac/trilinos_sparse_matrix.h>
#  endif


DEAL_II_NAMESPACE_OPEN

/**
 * 使用基类的通用格式独立输出例程以图形形式输出一个矩阵。矩阵被转换为一个二维领域上的补丁列表，其中高度由矩阵的元素给出。基准类的函数可以将这个矩阵的
 * "山形表示法
 * "写成各种图形输出格式。矩阵输出的坐标是：像往常一样，列从零开始随着X轴的增加而运行，而行则从零开始进入负Y轴。注意，由于一些内部限制，这个类一次只能输出一个矩阵，也就是说，它不能利用基类的多数据集功能。
 * 这个类的一个典型用法如下。
 *
 * @code
 * FullMatrix<double> M;
 * // fill matrix M with some values
 * ...
 *
 * // now write out M:
 * MatrixOut matrix_out;
 * std::ofstream out ("M.gnuplot");
 * matrix_out.build_patches (M, "M");
 * matrix_out.write_gnuplot (out);
 * @endcode
 * 当然，你也可以选择不同的图形输出格式。另外，这个类支持任何矩阵，不仅仅是FullMatrix类型的，只要它满足一些要求，用这个类的成员函数说明。
 * 通过build_patches()函数生成的补丁可以通过给它一个持有某些标志的对象来修改。关于这些标志的描述，请参见Options类成员的文档。
 *
 *
 *
 * @ingroup output
 *
 */
class MatrixOut : public DataOutInterface<2, 2>
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 持有各种变量的类，这些变量用于修改MatrixOut类的输出。
   *
   */
  struct Options
  {
    /**
     * 如果 @p true,
     * 只显示矩阵项的绝对值，而不是包括符号的真实值。默认值为
     * @p false. 。
     *
     */
    bool show_absolute_values;

    /**
     * 如果大于1，则不显示矩阵的每个元素，而是显示若干条目的平均值。输出补丁的数量也相应变小。这个标志决定了每个显示的区块应该有多大（以行/列为单位）。例如，如果它是两个，那么总是四个条目被整理成一个。
     * 默认值是1。
     *
     */
    unsigned int block_size;

    /**
     * 如果为真，绘制不连续的斑块，每个条目一个。
     *
     */
    bool discontinuous;

    /**
     * 默认构造函数。将此结构的所有元素设置为其默认值。
     *
     */
    Options(const bool         show_absolute_values = false,
            const unsigned int block_size           = 1,
            const bool         discontinuous        = false);
  };

  /**
   * 解构器。为了使它成为虚拟的而声明的。
   *
   */
  virtual ~MatrixOut() override = default;

  /**
   * 从给定的矩阵生成一个补丁列表，并在写入文件时使用给定的字符串作为数据集的名称。一旦补丁被建立，你可以使用基类的函数将数据写入文件，使用支持的输出格式之一。
   * 你可以给一个持有各种选项的结构。更多信息请看这个结构的字段描述。
   * 注意，这个函数要求我们能够提取矩阵的元素，这可以通过内部命名空间中声明的get_element()函数来完成。通过添加特殊化，你可以将这个类扩展到目前不支持的其他矩阵类。此外，我们需要能够提取矩阵的大小，为此我们假设矩阵类型提供了成员函数<tt>m()</tt>和<tt>n()</tt>，它们分别返回行和列的数量。
   *
   */
  template <class Matrix>
  void
  build_patches(const Matrix &     matrix,
                const std::string &name,
                const Options      options = Options(false, 1, false));

private:
  /**
   * 缩写 dealii::DataOutBase::Patch 类的有点冗长的名字。
   *
   */
  using Patch = DataOutBase::Patch<2, 2>;

  /**
   * 这是一个补丁列表，每次build_patches()被调用时都会创建一个补丁。这些补丁在基类的输出例程中使用。
   *
   */
  std::vector<Patch> patches;

  /**
   * 要写入的矩阵的名称。
   *
   */
  std::string name;

  /**
   * %函数，基类的函数通过这个函数知道他们应该写什么补丁到文件中。
   *
   */
  virtual const std::vector<Patch> &
  get_patches() const override;

  /**
   * 虚拟函数，基类的输出函数通过它获得数据集的名称。
   *
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * 获取网格点<tt>(i,j)</tt>处的矩阵值。根据给定的标志，这可能意味着不同的事情，例如，如果只应显示绝对值，则取矩阵条目的绝对值。如果块的大小大于1，那么就取几个矩阵项的平均值。
   *
   */
  template <class Matrix>
  static double
  get_gridpoint_value(const Matrix &  matrix,
                      const size_type i,
                      const size_type j,
                      const Options & options);
};


 /* ---------------------- Template and inline functions ------------- */ 


namespace internal
{
  namespace MatrixOutImplementation
  {
    /**
     * 返回稀疏矩阵中具有给定索引的元素。
     *
     */
    template <typename number>
    double
    get_element(const dealii::SparseMatrix<number> &matrix,
                const types::global_dof_index       i,
                const types::global_dof_index       j)
    {
      return matrix.el(i, j);
    }



    /**
     * 返回一个块状稀疏矩阵的给定指数的元素。
     *
     */
    template <typename number>
    double
    get_element(const dealii::BlockSparseMatrix<number> &matrix,
                const types::global_dof_index            i,
                const types::global_dof_index            j)
    {
      return matrix.el(i, j);
    }


#  ifdef DEAL_II_WITH_TRILINOS
    /**
     * 以给定的指数返回一个特里诺斯稀疏矩阵的元素。
     *
     */
    inline double
    get_element(const TrilinosWrappers::SparseMatrix &matrix,
                const types::global_dof_index         i,
                const types::global_dof_index         j)
    {
      return matrix.el(i, j);
    }



    /**
     * 以给定的指数返回一个特里诺斯块状稀疏矩阵的元素。
     *
     */
    inline double
    get_element(const TrilinosWrappers::BlockSparseMatrix &matrix,
                const types::global_dof_index              i,
                const types::global_dof_index              j)
    {
      return matrix.el(i, j);
    }
#  endif


#  ifdef DEAL_II_WITH_PETSC
    // no need to do anything: PETSc matrix objects do not distinguish
    // between operator() and el(i,j), so we can safely access elements
    // through the generic function below
#  endif


    /**
     * 从上面没有声明此函数特殊化的任何矩阵类型中，返回具有给定指数的元素。这将在矩阵上调用<tt>operator()</tt>。
     *
     */
    template <class Matrix>
    double
    get_element(const Matrix &                matrix,
                const types::global_dof_index i,
                const types::global_dof_index j)
    {
      return matrix(i, j);
    }
  } // namespace MatrixOutImplementation
} // namespace internal



template <class Matrix>
inline double
MatrixOut::get_gridpoint_value(const Matrix &  matrix,
                               const size_type i,
                               const size_type j,
                               const Options & options)
{
  // special case if block size is
  // one since we then don't need all
  // that loop overhead
  if (options.block_size == 1)
    {
      if (options.show_absolute_values == true)
        return std::fabs(
          internal::MatrixOutImplementation::get_element(matrix, i, j));
      else
        return internal::MatrixOutImplementation::get_element(matrix, i, j);
    }

  // if blocksize greater than one,
  // then compute average of elements
  double    average    = 0;
  size_type n_elements = 0;
  for (size_type row = i * options.block_size;
       row <
       std::min(size_type(matrix.m()), size_type((i + 1) * options.block_size));
       ++row)
    for (size_type col = j * options.block_size;
         col < std::min(size_type(matrix.m()),
                        size_type((j + 1) * options.block_size));
         ++col, ++n_elements)
      if (options.show_absolute_values == true)
        average += std::fabs(
          internal::MatrixOutImplementation::get_element(matrix, row, col));
      else
        average +=
          internal::MatrixOutImplementation::get_element(matrix, row, col);
  average /= n_elements;
  return average;
}



template <class Matrix>
void
MatrixOut::build_patches(const Matrix &     matrix,
                         const std::string &name,
                         const Options      options)
{
  size_type gridpoints_x = (matrix.n() / options.block_size +
                            (matrix.n() % options.block_size != 0 ? 1 : 0)),
            gridpoints_y = (matrix.m() / options.block_size +
                            (matrix.m() % options.block_size != 0 ? 1 : 0));

  // If continuous, the number of
  // plotted patches is matrix size-1
  if (!options.discontinuous)
    {
      --gridpoints_x;
      --gridpoints_y;
    }

  // first clear old data and set it
  // to virgin state
  patches.clear();
  patches.resize((gridpoints_x) * (gridpoints_y));

  // now build the patches
  size_type index = 0;
  for (size_type i = 0; i < gridpoints_y; ++i)
    for (size_type j = 0; j < gridpoints_x; ++j, ++index)
      {
        // within each patch, order the points in such a way that if some
        // graphical output program (such as gnuplot) plots the quadrilaterals
        // as two triangles, then the diagonal of the quadrilateral which cuts
        // it into the two printed triangles is parallel to the diagonal of the
        // matrix, rather than perpendicular to it. this has the advantage that,
        // for example, the unit matrix is plotted as a straight rim, rather
        // than as a series of bumps and valleys along the diagonal
        patches[index].vertices[0](0) = j;
        patches[index].vertices[0](1) = -static_cast<signed int>(i);
        patches[index].vertices[1](0) = j;
        patches[index].vertices[1](1) = -static_cast<signed int>(i + 1);
        patches[index].vertices[2](0) = j + 1;
        patches[index].vertices[2](1) = -static_cast<signed int>(i);
        patches[index].vertices[3](0) = j + 1;
        patches[index].vertices[3](1) = -static_cast<signed int>(i + 1);
        // next scale all the patch
        // coordinates by the block
        // size, to get original
        // coordinates
        for (auto &vertex : patches[index].vertices)
          vertex *= options.block_size;

        patches[index].n_subdivisions = 1;

        patches[index].data.reinit(1, 4);
        if (options.discontinuous)
          {
            patches[index].data(0, 0) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 1) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 2) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 3) =
              get_gridpoint_value(matrix, i, j, options);
          }
        else
          {
            patches[index].data(0, 0) =
              get_gridpoint_value(matrix, i, j, options);
            patches[index].data(0, 1) =
              get_gridpoint_value(matrix, i + 1, j, options);
            patches[index].data(0, 2) =
              get_gridpoint_value(matrix, i, j + 1, options);
            patches[index].data(0, 3) =
              get_gridpoint_value(matrix, i + 1, j + 1, options);
          }
      };

  // finally set the name
  this->name = name;
}



 /*----------------------------   matrix_out.h     ---------------------------*/ 

DEAL_II_NAMESPACE_CLOSE

#endif
 /*----------------------------   matrix_out.h     ---------------------------*/ 


