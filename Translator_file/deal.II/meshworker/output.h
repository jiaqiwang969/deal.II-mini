//include/deal.II-translator/meshworker/output_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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


#ifndef dealii_mesh_worker_output_h
#define dealii_mesh_worker_output_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/meshworker/dof_info.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    /**
     * 一个类，它不是组装成一个矩阵或向量，而是在一个单元格上输出结果到gnuplot补丁。
     * 这个汇编器期望LocalResults包含用
     * LocalResults::quadrature_value().
     * 设置的正交值，当它被初始化为单一（！）空间方向的正交点数量和要显示的数据字段数量时，它会自动初始化LocalResults。为了适应数据点的坐标，本地结果中的数据字段数量将被增加dim。
     * 虽然空间坐标的数据槽被自动分配，但这些坐标并没有被输入。这取决于用户在每个点的第一个dim数据项中输入坐标。这增加了输出转换后的坐标或甚至完全不同的东西的灵活性。
     * @note 在目前的实现中，只有单元格数据可以被写入。
     *
     */
    class GnuplotPatch
    {
    public:
      /**
       * 构造函数。
       *
       */
      GnuplotPatch();

      /**
       * 初始化用于写入<i>n</i>数据向量。点数是指张量积公式中单个方向的正交点的数量。它必须与用于创建补丁的实际正交点中的数量一致。产生的数据向量总数为<tt>n+dim</tt>，第一个dim应该是点的空间坐标。尽管如此，还是要由用户来设置这些值，以满足其需要。
       *
       */
      void
      initialize(const unsigned int n_points, const unsigned int n_vectors);

      /**
       * 设置数据被写入的流#os。如果没有用这个函数选择流，数据就会进入
       * @p deallog. 。
       *
       */
      void
      initialize_stream(std::ostream &stream);

      /**
       * 初始化DoFInfo对象中的本地数据，以后用于组装。
       * 如果 <code>!face</code> ， @p info
       * 对象指的是一个单元，否则指的是一个内部或边界面。
       *
       */
      template <int dim>
      void
      initialize_info(DoFInfo<dim> &info, bool face);

      /**
       * 将补丁写到输出流中。
       *
       */
      template <int dim>
      void
      assemble(const DoFInfo<dim> &info);

      /**
       * @warning  还没有实现
       *
       */
      template <int dim>
      void
      assemble(const DoFInfo<dim> &info1, const DoFInfo<dim> &info2);

    private:
      /**
       * 如果调用了initialize_stream()，则将对象T写到流#os，如果没有设置指针，则写到
       * @p deallog 。
       *
       */
      template <typename T>
      void
      write(const T &t) const;

      /**
       * 如果调用了initialize_stream，则向流#os写一个行结束标记，如果没有设置指针，则向
       * @p deallog 写。
       *
       */
      void
      write_endl() const;

      /**
       * 每个点的输出组件的数量。
       *
       */
      unsigned int n_vectors;
      /**
       * 一个方向上的点的数量。
       *
       */
      unsigned int n_points;

      /**
       * 输出要写入的流。由initialize_stream()设置。
       *
       */
      std::ostream *os;
    };

    //----------------------------------------------------------------------//

    template <typename T>
    inline void
    GnuplotPatch::write(const T &d) const
    {
      if (os == nullptr)
        deallog << d;
      else
        (*os) << d;
    }


    inline void
    GnuplotPatch::write_endl() const
    {
      if (os == nullptr)
        deallog << std::endl;
      else
        (*os) << std::endl;
    }


    inline GnuplotPatch::GnuplotPatch()
      : n_vectors(numbers::invalid_unsigned_int)
      , n_points(numbers::invalid_unsigned_int)
      , os(nullptr)
    {}


    inline void
    GnuplotPatch::initialize(const unsigned int np, const unsigned int nv)
    {
      n_vectors = nv;
      n_points  = np;
    }


    inline void
    GnuplotPatch::initialize_stream(std::ostream &stream)
    {
      os = &stream;
    }


    template <int dim>
    inline void
    GnuplotPatch::initialize_info(DoFInfo<dim> &info, bool face)
    {
      if (face)
        info.initialize_quadrature(Utilities::fixed_power<dim - 1>(n_points),
                                   n_vectors + dim);
      else
        info.initialize_quadrature(Utilities::fixed_power<dim>(n_points),
                                   n_vectors + dim);
    }


    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim> &info)
    {
      const unsigned int np = info.n_quadrature_points();
      const unsigned int nv = info.n_quadrature_values();
      const unsigned int patch_dim =
        (info.face_number == numbers::invalid_unsigned_int) ? dim : (dim - 1);
      const unsigned int row_length = n_points;
      // If patches are 1D, end the
      // patch after a row, else end
      // it after a square
      const unsigned int row_length2 =
        (patch_dim == 1) ? row_length : (row_length * row_length);

      //      AssertDimension(np, Utilities::fixed_power<dim>(n_points));
      AssertDimension(nv, n_vectors + dim);


      for (unsigned int k = 0; k < np; ++k)
        {
          if (k % row_length == 0)
            write_endl();
          if (k % row_length2 == 0)
            write_endl();

          for (unsigned int i = 0; i < nv; ++i)
            {
              write(info.quadrature_value(k, i));
              write('\t');
            }
          write_endl();
        }
    }


    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim> &info1, const DoFInfo<dim> &info2)
    {
      assemble(info1);
      assemble(info2);
    }
  } // namespace Assembler
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif


