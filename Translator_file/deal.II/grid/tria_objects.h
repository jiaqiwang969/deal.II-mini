//include/deal.II-translator/grid/tria_objects_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2021 by the deal.II authors
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

#ifndef dealii_tria_objects_h
#define dealii_tria_objects_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <class Accessor>
class TriaRawIterator;
template <int, int, int>
class TriaAccessor;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    /**
     * 属于三角测量的几何对象的信息的一般模板，即线、四边形、六面体......。
     * 除了对象的向量外，还包括其他信息，即表示子代、使用状态、用户标志、材料编码的向量。
     * 这些类的对象包括在TriaLevel和TriaFaces类中。
     *
     */
    class TriaObjects
    {
    public:
      /**
       * 构造函数重置了一些数据。
       *
       */
      TriaObjects();

      /**
       * 特定维度的构造函数。
       *
       */
      TriaObjects(const unsigned int structdim);

      unsigned int structdim;

      /**
       * 属于这个层面的对象的向量。对象的索引等于这个容器中的索引。
       *
       */
      std::vector<int> cells;

      /**
       * 返回该类存储的几何对象的数量。
       *
       */
      unsigned int
      n_objects() const;

      /**
       * 返回一个关于绑定当前对象所存储的 @p
       * 第1个索引对象的对象的索引的视图。例如，如果当前对象存储的是单元格，那么这个函数返回相当于一个包含绑定
       * @p index-th 单元格的面的索引的数组。
       *
       */
      ArrayView<int>
      get_bounding_object_indices(const unsigned int index);

      /**
       * 对象的偶数子女的索引。因为当对象被细化时，所有的子代都是同时被创建的，它们至少是成对地被追加到列表中，彼此之间的关系。因此，我们只存储偶数子女的索引，非偶数子女紧随其后。
       * 如果一个对象没有子女。
       *
       * - 被存储在这个列表中。如果一个对象没有子嗣，那么它就被称为活动对象。函数 TriaAccessorBase::has_children() 对此进行测试。
       *
       */
      std::vector<int> children;

      /**
       * 存储每个单元格被细化的情况。这个向量可能会被vector<vector<bool>
       * > (dim, vector<bool> (n_cells))取代，这样更节省内存。
       *
       */
      std::vector<std::uint8_t> refinement_cases;

      /**
       * 存储对象是否被用于 @p cells 向量的向量。
       * 由于在 @p vector,
       * 中很难删除元素，当一个元素不再需要时（例如，在取消定义后），它不会被从列表中删除，而是将相应的
       * @p used 标志设置为 @p false.  。
       *
       */
      std::vector<bool> used;

      /**
       * 为用户数据提供一个字段，每个对象有一个比特。这个字段通常用于当一个操作在所有单元上运行，并且需要另一个单元（例如邻居）是否已经被处理的信息。
       * 你可以使用 Triangulation::clear_user_flags().
       * 清除所有使用的标志。
       *
       */
      std::vector<bool> user_flags;


      /**
       * 我们使用这个联盟来存储边界和材料数据。因为这里实际上只需要这两个中的一个，所以我们使用一个联合。
       *
       */
      struct BoundaryOrMaterialId
      {
        union
        {
          types::boundary_id boundary_id;
          types::material_id material_id;
        };


        /**
         * 默认构造函数。
         *
         */
        BoundaryOrMaterialId();

        /**
         * 返回这类对象的大小。
         *
         */
        static std::size_t
        memory_consumption();

        /**
         * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
         *
         */
        template <class Archive>
        void
        serialize(Archive &ar, const unsigned int version);
      };

      /**
       * 存储边界和材料数据。例如，在一个维度中，这个字段存储线的材料ID，这是一个介于0和
       * numbers::invalid_material_id-1. 之间的数字。
       * 在一个以上的维度中，线没有材料ID，但它们可能在边界上；然后，我们在这个字段中存储边界指示器，它表示这条线属于边界的哪一部分，以及这一部分的边界条件是什么。边界指示符也是一个介于零和
       * numbers::internal_face_boundary_id-1; 之间的数字，id
       * numbers::internal_face_boundary_id
       * 是为内部的线保留的，可以用来检查一条线是否在边界上，否则，如果你不知道它属于哪个单元格，就不可能。
       *
       */
      std::vector<BoundaryOrMaterialId> boundary_or_material_id;

      /**
       * 存储分流板ID。这个字段存储每个对象的流形ID，它是0到
       * numbers::flat_manifold_id-1. 之间的一个数字。
       *
       */
      std::vector<types::manifold_id> manifold_id;

      /**
       * 返回一个通往单个对象的下一个空闲槽的迭代器。这个函数只被3D中的
       * Triangulation::execute_refinement() 使用。              @warning
       * 有趣的是，这个函数没有用于一维或二维三角计算，在这里，似乎细化函数的作者坚持要重新实现其内容。
       * @todo  这个函数没有被实例化用于codim-one的情况
       *
       */
      template <int structdim, int dim, int spacedim>
      dealii::TriaRawIterator<dealii::TriaAccessor<structdim, dim, spacedim>>
      next_free_single_object(const Triangulation<dim, spacedim> &tria);

      /**
       * 返回一个通往一对对象的下一个空闲槽的迭代器。这个函数只被3D中的
       * Triangulation::execute_refinement() 使用。              @warning
       * 有趣的是，这个函数并没有用于一维或二维三角计算，在这里，似乎细化函数的作者坚持重新实现其内容。
       * @todo  这个函数没有被实例化用于codim-one的情况
       *
       */
      template <int structdim, int dim, int spacedim>
      dealii::TriaRawIterator<dealii::TriaAccessor<structdim, dim, spacedim>>
      next_free_pair_object(const Triangulation<dim, spacedim> &tria);

      /**
       * 返回一个迭代器到一对六角的下一个空闲槽。只对
       * <code>G=Hexahedron</code> 实现。
       *
       */
      template <int dim, int spacedim>
      typename Triangulation<dim, spacedim>::raw_hex_iterator
      next_free_hex(const Triangulation<dim, spacedim> &tria,
                    const unsigned int                  level);

      /**
       * 访问用户指针。
       *
       */
      void *&
      user_pointer(const unsigned int i);

      /**
       * 对用户指针的只读访问。
       *
       */
      const void *
      user_pointer(const unsigned int i) const;

      /**
       * 对用户索引的访问。
       *
       */
      unsigned int &
      user_index(const unsigned int i);

      /**
       * 对用户指针的只读访问。
       *
       */
      unsigned int
      user_index(const unsigned int i) const;

      /**
       * 将用户数据重置为零。
       *
       */
      void
      clear_user_data(const unsigned int i);

      /**
       * 清除所有的用户指针或索引，并重置它们的类型，这样，下一次访问可能是或者。
       *
       */
      void
      clear_user_data();

      /**
       * 清除所有的用户标志。
       *
       */
      void
      clear_user_flags();

      /**
       * 确定这个对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);

      /**
       * Triangulation对象可以访问用户指针或用户索引。你试图做的是在使用其中一个之后再访问另一个。
       * @ingroup Exceptions
       *
       */
      DeclException0(ExcPointerIndexClash);

      /**
       * next_free_single_*函数的计数器
       *
       */
      unsigned int next_free_single;

      /**
       * 下一个免费配对函数的计数器
       *
       */
      unsigned int next_free_pair;

      /**
       * 用于next_free_single_*函数的Bool标志
       *
       */
      bool reverse_order_next_free_single;

      /**
       * 存储用户指针或用户索引的数据类型。
       *
       */
      struct UserData
      {
        union
        {
          /// The entry used as user
          /// pointer.
          void *p;
          /// The entry used as user
          /// index.
          unsigned int i;
        };

        /**
         * 默认的构造函数。
         *
         */
        UserData()
        {
          p = nullptr;
        }

        /**
         * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
         *
         */
        template <class Archive>
        void
        serialize(Archive &ar, const unsigned int version);
      };

      /**
       * 描述用户数据的可能类型的枚举。
       *
       */
      enum UserDataType
      {
        /// No userdata used yet.
        data_unknown,
        /// UserData contains pointers.
        data_pointer,
        /// UserData contains indices.
        data_index
      };


      /**
       * 指针，不被库使用，但可以被用户访问和设置，以处理一行/四行/等的本地数据。
       *
       */
      std::vector<UserData> user_data;

      /**
       * 为了避免用户指针和索引之间的混淆，这个枚举由第一个访问这两者的函数设置，随后的访问将不被允许改变所访问的数据类型。
       *
       */
      mutable UserDataType user_data_type;
    };


    //----------------------------------------------------------------------//

    inline unsigned int
    TriaObjects::n_objects() const
    {
      // assume that each cell has the same number of faces
      const unsigned int faces_per_cell = 2 * this->structdim;
      return cells.size() / faces_per_cell;
    }



    inline ArrayView<int>
    TriaObjects::get_bounding_object_indices(const unsigned int index)
    {
      // assume that each cell has the same number of faces
      const unsigned int faces_per_cell = 2 * this->structdim;
      return ArrayView<int>(cells.data() + index * faces_per_cell,
                            faces_per_cell);
    }



    inline TriaObjects::BoundaryOrMaterialId::BoundaryOrMaterialId()
    {
      material_id = numbers::invalid_material_id;
    }



    inline std::size_t
    TriaObjects::BoundaryOrMaterialId::memory_consumption()
    {
      return sizeof(BoundaryOrMaterialId);
    }



    template <class Archive>
    void
    TriaObjects::BoundaryOrMaterialId::serialize(Archive &ar,
                                                 const unsigned int  /*version*/ )
    {
      // serialize this
      // structure by
      // writing and
      // reading the larger
      // of the two values,
      // in order to make
      // sure we get all
      // bits
      if (sizeof(material_id) > sizeof(boundary_id))
        ar &material_id;
      else
        ar &boundary_id;
    }


    inline void *&
    TriaObjects::user_pointer(const unsigned int i)
    {
      Assert(user_data_type == data_unknown || user_data_type == data_pointer,
             ExcPointerIndexClash());
      user_data_type = data_pointer;

      AssertIndexRange(i, user_data.size());
      return user_data[i].p;
    }


    inline const void *
    TriaObjects::user_pointer(const unsigned int i) const
    {
      Assert(user_data_type == data_unknown || user_data_type == data_pointer,
             ExcPointerIndexClash());
      user_data_type = data_pointer;

      AssertIndexRange(i, user_data.size());
      return user_data[i].p;
    }


    inline unsigned int &
    TriaObjects::user_index(const unsigned int i)
    {
      Assert(user_data_type == data_unknown || user_data_type == data_index,
             ExcPointerIndexClash());
      user_data_type = data_index;

      AssertIndexRange(i, user_data.size());
      return user_data[i].i;
    }


    inline void
    TriaObjects::clear_user_data(const unsigned int i)
    {
      AssertIndexRange(i, user_data.size());
      user_data[i].i = 0;
    }


    inline TriaObjects::TriaObjects()
      : structdim(numbers::invalid_unsigned_int)
      , next_free_single(numbers::invalid_unsigned_int)
      , next_free_pair(numbers::invalid_unsigned_int)
      , reverse_order_next_free_single(false)
      , user_data_type(data_unknown)
    {}


    inline TriaObjects::TriaObjects(const unsigned int structdim)
      : structdim(structdim)
      , next_free_single(numbers::invalid_unsigned_int)
      , next_free_pair(numbers::invalid_unsigned_int)
      , reverse_order_next_free_single(false)
      , user_data_type(data_unknown)
    {}


    inline unsigned int
    TriaObjects::user_index(const unsigned int i) const
    {
      Assert(user_data_type == data_unknown || user_data_type == data_index,
             ExcPointerIndexClash());
      user_data_type = data_index;

      AssertIndexRange(i, user_data.size());
      return user_data[i].i;
    }


    inline void
    TriaObjects::clear_user_data()
    {
      user_data_type = data_unknown;
      for (auto &data : user_data)
        data.p = nullptr;
    }


    inline void
    TriaObjects::clear_user_flags()
    {
      user_flags.assign(user_flags.size(), false);
    }


    template <class Archive>
    void
    TriaObjects::UserData::serialize(Archive &ar, const unsigned int)
    {
      // serialize this as an integer
      ar &i;
    }



    template <class Archive>
    void
    TriaObjects::serialize(Archive &ar, const unsigned int)
    {
      ar &structdim;
      ar &cells &children;
      ar &       refinement_cases;
      ar &       used;
      ar &       user_flags;
      ar &       boundary_or_material_id;
      ar &       manifold_id;
      ar &next_free_single &next_free_pair &reverse_order_next_free_single;
      ar &user_data &user_data_type;
    }


    //----------------------------------------------------------------------//

    template <int structdim_, int dim, int spacedim>
    dealii::TriaRawIterator<dealii::TriaAccessor<structdim_, dim, spacedim>>
    TriaObjects::next_free_single_object(
      const Triangulation<dim, spacedim> &tria)
    {
      // TODO: Think of a way to ensure that we are using the correct
      // triangulation, i.e. the one containing *this.

      AssertDimension(structdim_, this->structdim);

      int pos = next_free_single, last = used.size() - 1;
      if (!reverse_order_next_free_single)
        {
          // first sweep forward, only use really single slots, do not use
          // pair slots
          for (; pos < last; ++pos)
            if (!used[pos])
              if (used[++pos])
                {
                  // this was a single slot
                  pos -= 1;
                  break;
                }
          if (pos >= last)
            {
              reverse_order_next_free_single = true;
              next_free_single               = used.size() - 1;
              pos                            = used.size() - 1;
            }
          else
            next_free_single = pos + 1;
        }

      if (reverse_order_next_free_single)
        {
          // second sweep, use all slots, even
          // in pairs
          for (; pos >= 0; --pos)
            if (!used[pos])
              break;
          if (pos > 0)
            next_free_single = pos - 1;
          else
            // no valid single object anymore
            return dealii::TriaRawIterator<
              dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, -1, -1);
        }

      return dealii::TriaRawIterator<
        dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, 0, pos);
    }



    template <int structdim_, int dim, int spacedim>
    dealii::TriaRawIterator<dealii::TriaAccessor<structdim_, dim, spacedim>>
    TriaObjects::next_free_pair_object(const Triangulation<dim, spacedim> &tria)
    {
      // TODO: Think of a way to ensure that we are using the correct
      // triangulation, i.e. the one containing *this.

      AssertDimension(structdim_, this->structdim);

      int pos = next_free_pair, last = used.size() - 1;
      for (; pos < last; ++pos)
        if (!used[pos])
          if (!used[++pos])
            {
              // this was a pair slot
              pos -= 1;
              break;
            }
      if (pos >= last)
        // no free slot
        return dealii::TriaRawIterator<
          dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, -1, -1);
      else
        next_free_pair = pos + 2;

      return dealii::TriaRawIterator<
        dealii::TriaAccessor<structdim_, dim, spacedim>>(&tria, 0, pos);
    }
  } // namespace TriangulationImplementation
} // namespace internal



DEAL_II_NAMESPACE_CLOSE

#endif


