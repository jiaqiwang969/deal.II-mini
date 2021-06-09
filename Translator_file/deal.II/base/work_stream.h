//include/deal.II-translator/base/work_stream_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_work_stream_h
#  define dealii_work_stream_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/graph_coloring.h>
#  include <deal.II/base/iterator_range.h>
#  include <deal.II/base/multithread_info.h>
#  include <deal.II/base/parallel.h>
#  include <deal.II/base/template_constraints.h>
#  include <deal.II/base/thread_local_storage.h>
#  include <deal.II/base/thread_management.h>

#  ifdef DEAL_II_WITH_TBB
#    include <tbb/pipeline.h>
#  endif

#  include <functional>
#  include <iterator>
#  include <memory>
#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN



/**
 * 一个命名空间，其主要的模板函数支持运行多个线程，每个线程在给定的对象范围内的一个子集上操作。该类使用英特尔线程构件（TBB）来平衡各个子范围在可用线程上的负载。关于这个类的原理的长篇讨论，见 @ref threads "多处理器的并行计算 "
 * 模块。它在教程中首先用于 step-9 ，然后在 step-13 、
 * step-14 、 step-32 和其他地方再次使用。
 * 该类建立在以下前提之上：人们经常有一些工作需要在一连串的对象上完成；一个典型的例子是将单元格的贡献集合到一个系统矩阵或右手边。在许多这样的例子中，部分工作可以完全独立和并行地完成，可能使用一台具有共享内存的机器上的几个处理器核心。然而，这项工作的其他部分可能需要同步进行，并按顺序完成。在组装矩阵的例子中，局部贡献的计算可以完全并行完成，但将局部贡献复制到全局矩阵中需要一些注意。首先，几个线程不能同时写入，而是需要使用mutex来同步写入；其次，我们希望本地贡献加入全局矩阵的顺序总是相同的，因为浮点加法不是换元的，以不同的顺序将本地贡献加入全局矩阵会导致微妙的不同结果，这可能会影响迭代求解器的迭代次数，以及随机的解的舍弃误差。因此，我们要确保每次只有一个线程写入全局矩阵，而且结果是以稳定和可重复的顺序复制的。
 * 这个类实现了这种工作模式的框架。它处理一个由迭代器范围给定的对象流，在所有这些对象上并行地运行一个工作函数，然后将每个对象传递给一个后处理函数，该函数按顺序运行，并完全按照对象在输入迭代器范围中出现的顺序获得对象。所有的同步工作都不会暴露给这个类的用户。
 * 在内部，给这个类的run()函数的范围被分割成一连串的
 * "项目"，然后根据一些%的内部算法分配到可用线程的数量上。一个项目是我们要操作的迭代器范围中的一个元素；例如，为了组装矩阵或评估错误指标，一个项目可以是一个单元。TBB库决定了创建多少个线程（通常与处理器核心一样多），但在任何特定时间可能处于活动状态的项的数量由构造函数的参数指定。它应该大于或等于处理器内核的数量
 *
 * - 默认是当前系统中核心数量的四倍。
 * 每当一个工作线程处于空闲状态或预计将成为空闲状态时，TBB就会根据请求创建项目。然后，它被移交给一个工作者函数，通常是一个主类的成员函数。这些工作函数在一些线程上并行运行，而且不能保证它们被要求以任何特定的顺序处理项目，特别是不一定以项目从迭代器范围生成的顺序。
 * 通常情况下，工作者函数需要额外的数据，例如FEValues对象、输入数据向量等，其中有些数据不能在线程之间共享。为此，run()函数需要另一个模板参数ScratchData，它指定了一种类型的对象，这些对象与每个项目一起存储，线程可以将其作为私有数据使用而不必与其他线程共享。run()函数需要一个额外的参数，其中有一个ScratchData类型的对象，该对象将被复制为传递给每个工作者函数的参数。
 * 此外，工作者函数将其结果存储在模板类型CopyData的对象中。然后，这些被移交给一个单独的函数，称为copier，它可能使用存储的结果，将它们转移到永久存储中。例如，它可以将工作函数计算的矩阵的局部贡献的结果复制到全局矩阵中。然而，与工作函数不同的是，在任何给定的时间内，只有一个复制器的实例在运行；因此，它可以安全地将本地贡献复制到全局矩阵中，而不需要使用突变器或类似手段锁定全局对象。此外，它保证复制器与CopyData对象的运行顺序与相关项目的创建顺序相同；因此，即使工作线程可能以非指定的顺序计算结果，复制器也总是以项目创建时的相同顺序接收结果。
 * 一旦一个项目被复制器处理，它就会被删除，在其计算中使用的ScratchData和CopyData对象被认为是未使用的，并且可以在这个线程或另一个线程上的工作者函数的下一次调用中重新使用。然而，WorkStream函数没有试图将这些对象重置为任何一种原始状态
 *
 * - 工作者应该假定它得到的CopyData对象有先前的内容，并以任何看起来合适的方式首先清除它，然后再把内容放进它，以后可以由复制者再次处理。
 * ScratchData和CopyData中的成员变量可以独立于这些数据结构副本的其他并发使用而被访问。因此，将与ScratchData和CopyData相关的辅助数据结构的大小调整到每个单元上的不同长度是完全可以的。例如，与
 * LocalIntegrators::L2::weighted_mass_matrix()
 * 一起用于组装局部矩阵的持有每个正交点密度的向量可以被调整为具有hp-capabilities的DoFHandlers中当前单元的相应正交点的数量。同样，CopyData中的局部刚度矩阵也可以根据当前单元上的局部DoF的数量来调整大小。
 *
 *
 * @note
 * 对于单元格和面的积分，使用比当前函数更具体的方法往往是有用的（它不在乎迭代器是在单元格、向量元素还是其他类型的范围上）。
 * MeshWorker::mesh_loop()
 * 函数是一个特别适合于积分的接口的实现。
 *
 *
 * @note
 * 这个命名空间中的函数只有在deal.II配置时选择了多线程模式时才会真正并行工作。否则，它们只是按顺序对每个项目工作。
 *
 *
 * @ingroup threads
 *
 */
namespace WorkStream
{
  /**
   * 嵌套命名空间包含工作流算法的各种实现。
   *
   */
  namespace internal
  {
#  ifdef DEAL_II_WITH_TBB
    /**
     * 一个命名空间，用于实现WorkStream模式和函数的细节。这个命名空间持有处理Turcksin、Kronbichler和Bangerth的论文中描述的第二个实现的类（见 @ref
     * workstream_paper  ）。
     * 在这里，没有提供着色，所以复制是使用TBB过滤器顺序进行的。
     * 尽管这个实现比那篇论文中讨论的第三个实现要慢，我们还是需要保留它，原因有二。(i)用户可能不会给我们一个图形着色，(ii)我们想用这个实现来处理那些太小的颜色。
     *
     */
    namespace tbb_no_coloring
    {
      /**
       * 一个可以从一系列迭代器中创建一个项目序列的类。
       *
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class IteratorRangeToItemStream : public tbb::filter
      {
      public:
        /**
         * 一个数据类型，我们用它来识别要处理的项目。这是一个结构，在WorkStream实现的不同部分之间传递，以确定流水线的各个阶段需要做什么。
         *
         */
        struct ItemType
        {
          /**
           * 一个结构，它包含一个指向抓取数据对象的指针，以及一个指示该对象当前是否正在使用的标志。
           *
           */
          struct ScratchDataObject
          {
            std::unique_ptr<ScratchData> scratch_data;
            bool                         currently_in_use;

            /**
             * 默认构造函数。
             *
             */
            ScratchDataObject()
              : currently_in_use(false)
            {}

            ScratchDataObject(ScratchData *p, const bool in_use)
              : scratch_data(p)
              , currently_in_use(in_use)
            {}

            // Provide a copy constructor that actually doesn't copy the
            // internal state. This makes handling ScratchAndCopyDataObjects
            // easier to handle with STL containers.
            ScratchDataObject(const ScratchDataObject &)
              : currently_in_use(false)
            {}

            ScratchDataObject(ScratchDataObject &&o) noexcept = default;
          };


          /**
           * 类型化为一个从头数据对象的列表。这个列表的原理是在使用这些列表的变量中提供的。
           *
           */
          using ScratchDataList = std::list<ScratchDataObject>;

          /**
           * 一个需要被处理的迭代器的列表。只有前面的n_个项目是相关的。
           *
           */
          std::vector<Iterator> work_items;

          /**
           * 管道的Worker部分为每个工作项目填充的CopyData对象。同样，只有前n_items元素是我们所关心的。
           *
           */
          std::vector<CopyData> copy_datas;

          /**
           * Work_items数组所标识的、Worker和Copier管道阶段需要工作的项目数量。这个变量的最大值将是chunk_size。
           *
           */
          unsigned int n_items;

          /**
           * 指向线程局部变量的指针，识别该线程将使用的抓取数据对象。这个类的最初实现是使用线程局部变量，每个线程只提供一个从头开始的对象。这并不奏效，因为工作者函数可能会自己启动任务，然后调用
           * Threads::TaskGroup::join_all()
           * 或类似的函数，TBB调度器可能会用它来在当前线程上运行其他东西
           *
           *
           *
           *
           *
           *
           * - 例如工作者函数的另一个实例。          因此，如果我们只为每个线程提供一个抓取对象，就会有两个工作者函数的实例使用同一个抓取对象。解决方案是为每个线程提供一个从头开始的对象列表，同时提供一个标志，表明这个从头开始的对象目前是否被使用。如果一个线程需要一个抓取对象，它就会浏览这个列表，直到找到一个未使用的对象，或者，如果没有，就自己创建一个。请注意，我们不需要为这个过程使用同步原语，因为这些列表是线程本地的，只要我们在访问列表之间没有屈服点，就能保证只有一个线程访问它们。                    存放在这些列表中的指向scratch对象的指针必须是这样的：当线程本地对象被销毁时，它们会在所有线程上被删除。这可以通过使用unique_ptr来实现。                    请注意，当一个工作者需要创建一个从头开始的对象时，它会使用sample_scratch_data来分配它，以便从中复制。这样做的好处是第一次接触初始化，也就是说，从头数据对象的内存是由后来使用它的同一个线程分配和初始化的。
           *
           */
          Threads::ThreadLocalStorage<ScratchDataList> *scratch_data;

          /**
           * 指针指向一个样本的抓取数据对象，用于初始化为每个单独线程创建的抓取数据对象。
           *
           */
          const ScratchData *sample_scratch_data;

          /**
           * 如果缓冲区被使用，则标志为真；如果缓冲区可以使用，则标志为假。
           *
           */
          bool currently_in_use;


          /**
           * 默认构造函数。初始化所有本身没有默认构造函数的东西。
           *
           */
          ItemType()
            : n_items(0)
            , scratch_data(nullptr)
            , sample_scratch_data(nullptr)
            , currently_in_use(false)
          {}
        };


        /**
         * 构造函数。取一个迭代器范围，可以容纳项目的缓冲区的大小，以及将传递给每个工作者和复制者函数调用的附加数据对象样本。
         *
         */
        IteratorRangeToItemStream(const Iterator &   begin,
                                  const Iterator &   end,
                                  const unsigned int buffer_size,
                                  const unsigned int chunk_size,
                                  const ScratchData &sample_scratch_data,
                                  const CopyData &   sample_copy_data)
          : tbb::filter( /*is_serial=*/ true)
          , remaining_iterator_range(begin, end)
          , item_buffer(buffer_size)
          , sample_scratch_data(sample_scratch_data)
          , chunk_size(chunk_size)
        {
          // initialize the elements of the ring buffer
          for (unsigned int element = 0; element < item_buffer.size();
               ++element)
            {
              Assert(item_buffer[element].n_items == 0, ExcInternalError());

              item_buffer[element].work_items.resize(
                chunk_size, remaining_iterator_range.second);
              item_buffer[element].scratch_data        = &thread_local_scratch;
              item_buffer[element].sample_scratch_data = &sample_scratch_data;
              item_buffer[element].copy_datas.resize(chunk_size,
                                                     sample_copy_data);
              item_buffer[element].currently_in_use = false;
            }
        }


        /**
         * 创建一个项目并返回它的指针。
         *
         */
        virtual void *
        operator()(void *) override
        {
          // find first unused item. we know that there must be one
          // because we have set the maximal number of tokens in flight
          // and have set the ring buffer to have exactly this size. so
          // if this function is called, we know that less than the
          // maximal number of items in currently in flight
          //
          // note that we need not lock access to this array since
          // the current stage is run sequentially and we can therefore
          // enter the following block only once at any given time.
          // thus, there can be no race condition between checking that
          // a flag is false and setting it to true. (there may be
          // another thread where we release items and set 'false'
          // flags to 'true', but that too does not produce any
          // problems)
          ItemType *current_item = nullptr;
          for (unsigned int i = 0; i < item_buffer.size(); ++i)
            if (item_buffer[i].currently_in_use == false)
              {
                item_buffer[i].currently_in_use = true;
                current_item                    = &item_buffer[i];
                break;
              }
          Assert(current_item != nullptr,
                 ExcMessage("This can't be. There must be a free item!"));

          // initialize the next item. it may
          // consist of at most chunk_size
          // elements
          current_item->n_items = 0;
          while ((remaining_iterator_range.first !=
                  remaining_iterator_range.second) &&
                 (current_item->n_items < chunk_size))
            {
              current_item->work_items[current_item->n_items] =
                remaining_iterator_range.first;

              ++remaining_iterator_range.first;
              ++current_item->n_items;
            }

          if (current_item->n_items == 0)
            // there were no items
            // left. terminate the pipeline
            return nullptr;
          else
            return current_item;
        }

      private:
        /**
         * 仍待处理的迭代器的区间。这个范围会随着时间的推移而缩小。
         *
         */
        std::pair<Iterator, Iterator> remaining_iterator_range;

        /**
         * 一个缓冲区，它将存储项目。
         *
         */
        std::vector<ItemType> item_buffer;

        /**
         * 指向线程局部变量的指针，用于识别该线程将使用的抓取数据对象。这个类的最初实现是使用线程局部变量，每个线程只提供一个抓取对象。这并不奏效，因为工作者函数可能会自己启动任务，然后调用
         * Threads::TaskGroup::join_all()
         * 或类似的函数，TBB调度器可能会用它来在当前线程上运行其他东西
         *
         * --例如，工作函数的另一个实例。因此，如果我们只为每个线程提供一个划痕对象，就会有两个工作函数的实例使用同一个划痕对象。解决办法是为每个线程提供一个从头开始的对象列表，同时提供一个标志，表明这个从头开始的对象目前是否被使用。如果一个线程需要一个抓取对象，它就会浏览这个列表，直到找到一个未使用的对象，或者，如果没有，就自己创建一个。请注意，我们不需要为这个过程使用同步原语，因为这些列表是线程本地的，只要我们在访问列表之间没有屈服点，就能保证只有一个线程访问它们。
         * 存放在这些列表中的指向scratch对象的指针必须是这样的：当线程本地对象被销毁时，它们会在所有线程上被删除。这可以通过使用unique_ptr来实现。
         * 请注意，当一个工作者需要创建一个从头开始的对象时，它会使用sample_scratch_data来分配它，以便从中复制。这样做的好处是第一次接触初始化，也就是说，划痕数据对象的内存是由后来使用它的同一个线程分配和初始化的。
         *
         */
        Threads::ThreadLocalStorage<typename ItemType::ScratchDataList>
          thread_local_scratch;

        /**
         * 对抓取数据样本的引用，该样本将被用于初始化每个工作任务所使用的抓取数据对象的线程本地指针。
         *
         */
        const ScratchData &sample_scratch_data;

        /**
         * 每个线程应按顺序工作的迭代器范围的元素数；一个大的数字可以确保每个线程在下一个任务切换发生之前得到大量的工作，而一个小的数字则更有利于负载平衡。
         *
         */
        const unsigned int chunk_size;
      };



      /**
       * 一个管理在若干并行线程上调用工作者函数的类。请注意，用TBB的说法，它是一个可以并行运行的过滤器。
       *
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class TBBWorker : public tbb::filter
      {
      public:
        /**
         * 构造函数。取一个对我们要操作的对象的引用，以及一个指向将进行装配的函数的指针。
         *
         */
        TBBWorker(
          const std::function<void(const Iterator &, ScratchData &, CopyData &)>
            &  worker,
          bool copier_exist = true)
          : tbb::filter( /* is_serial= */  false)
          , worker(worker)
          , copier_exist(copier_exist)
        {}


        /**
         * 对一个项目进行操作。
         *
         */
        void *
        operator()(void *item) override
        {
          // first unpack the current item
          using ItemType =
            typename IteratorRangeToItemStream<Iterator,
                                               ScratchData,
                                               CopyData>::ItemType;

          ItemType *current_item = static_cast<ItemType *>(item);

          // we need to find an unused scratch data object in the list that
          // corresponds to the current thread and then mark it as used. if
          // we can't find one, create one
          //
          // as discussed in the discussion of the documentation of the
          // IteratorRangeToItemStream::scratch_data variable, there is no
          // need to synchronize access to this variable using a mutex
          // as long as we have no yield-point in between. this means that
          // we can't take an iterator into the list now and expect it to
          // still be valid after calling the worker, but we at least do
          // not have to lock the following section
          ScratchData *scratch_data = nullptr;
          {
            typename ItemType::ScratchDataList &scratch_data_list =
              current_item->scratch_data->get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ItemType::ScratchDataList::iterator p =
                   scratch_data_list.begin();
                 p != scratch_data_list.end();
                 ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data        = p->scratch_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no object was found, create one and mark it as used
            if (scratch_data == nullptr)
              {
                scratch_data =
                  new ScratchData(*current_item->sample_scratch_data);

                typename ItemType::ScratchDataList::value_type
                  new_scratch_object(scratch_data, true);
                scratch_data_list.push_back(std::move(new_scratch_object));
              }
          }

          // then call the worker function on each element of the chunk we were
          // given. since these worker functions are called on separate threads,
          // nothing good can happen if they throw an exception and we are best
          // off catching it and showing an error message
          for (unsigned int i = 0; i < current_item->n_items; ++i)
            {
              try
                {
                  if (worker)
                    worker(current_item->work_items[i],
                           *scratch_data,
                           current_item->copy_datas[i]);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception(exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            typename ItemType::ScratchDataList &scratch_data_list =
              current_item->scratch_data->get();

            for (typename ItemType::ScratchDataList::iterator p =
                   scratch_data_list.begin();
                 p != scratch_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }

          // if there is no copier, mark current item as usable again
          if (copier_exist == false)
            current_item->currently_in_use = false;


          // then return the original pointer
          // to the now modified object
          return item;
        }


      private:
        /**
         * 指向对单元格序列进行装配的函数的指针。
         *
         */
        const std::function<void(const Iterator &, ScratchData &, CopyData &)>
          worker;

        /**
         * 如果复制器阶段存在，该标志为真。如果它不存在，工作者必须释放缓冲区。否则，复制者将会做这件事。
         *
         */
        bool copier_exist;
      };



      /**
       * 一个管理调用复制器函数的类。请注意，在TBB符号中，它是一个按顺序运行的过滤器，确保所有项目按照它们被创建的相同顺序被复制。
       *
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class TBBCopier : public tbb::filter
      {
      public:
        /**
         * 构造函数。接受一个对我们要操作的对象的引用，以及一个指向函数的指针，该函数将完成从附加数据对象到全局矩阵的复制或类似的复制。
         *
         */
        TBBCopier(const std::function<void(const CopyData &)> &copier)
          : tbb::filter( /*is_serial=*/ true)
          , copier(copier)
        {}


        /**
         * 在一个单项上工作。
         *
         */
        void *
        operator()(void *item) override
        {
          // first unpack the current item
          using ItemType =
            typename IteratorRangeToItemStream<Iterator,
                                               ScratchData,
                                               CopyData>::ItemType;

          ItemType *current_item = static_cast<ItemType *>(item);

          // initiate copying data. for the same reasons as in the worker class
          // above, catch exceptions rather than letting it propagate into
          // unknown territories
          for (unsigned int i = 0; i < current_item->n_items; ++i)
            {
              try
                {
                  if (copier)
                    copier(current_item->copy_datas[i]);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception(exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception();
                }
            }

          // mark current item as usable again
          current_item->currently_in_use = false;


          // return an invalid item since we are at the end of the
          // pipeline
          return nullptr;
        }


      private:
        /**
         * 指向进行数据复制的函数的指针。
         *
         */
        const std::function<void(const CopyData &)> copier;
      };

      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const Iterator &                         begin,
          const typename identity<Iterator>::type &end,
          Worker                                   worker,
          Copier                                   copier,
          const ScratchData &                      sample_scratch_data,
          const CopyData &                         sample_copy_data,
          const unsigned int                       queue_length,
          const unsigned int                       chunk_size)
      {
        // create the three stages of the pipeline
        IteratorRangeToItemStream<Iterator, ScratchData, CopyData>
          iterator_range_to_item_stream(begin,
                                        end,
                                        queue_length,
                                        chunk_size,
                                        sample_scratch_data,
                                        sample_copy_data);

        TBBWorker<Iterator, ScratchData, CopyData> worker_filter(worker);
        TBBCopier<Iterator, ScratchData, CopyData> copier_filter(copier);

        // now create a pipeline from these stages
        tbb::pipeline assembly_line;
        assembly_line.add_filter(iterator_range_to_item_stream);
        assembly_line.add_filter(worker_filter);
        assembly_line.add_filter(copier_filter);

        // and run it
        assembly_line.run(queue_length);

        assembly_line.clear();
      }

    }    // namespace tbb_no_coloring
#  endif // DEAL_II_WITH_TBB


    /**
     * 一个不使用多线程的参考实现，如果我们没有多线程支持，或者用户要求按顺序运行的话，就可以使用。如果我们只有一个单线程，这比使用TBB或taskflow更有效率。
     *
     */
    namespace sequential
    {
      /**
       * 没有颜色的顺序版本。
       *
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const Iterator &                         begin,
          const typename identity<Iterator>::type &end,
          Worker                                   worker,
          Copier                                   copier,
          const ScratchData &                      sample_scratch_data,
          const CopyData &                         sample_copy_data)
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data; // NOLINT

        // Optimization: Check if the functions are not the zero function. To
        // check zero-ness, create a C++ function out of it:
        const bool have_worker =
          (static_cast<const std::function<
             void(const Iterator &, ScratchData &, CopyData &)> &>(worker)) !=
          nullptr;
        const bool have_copier =
          (static_cast<const std::function<void(const CopyData &)> &>(
            copier)) != nullptr;

        // Finally loop over all items and perform the necessary work:
        for (Iterator i = begin; i != end; ++i)
          {
            if (have_worker)
              worker(i, scratch_data, copy_data);
            if (have_copier)
              copier(copy_data);
          }
      }



      /**
       * 有颜色的顺序版本
       *
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const std::vector<std::vector<Iterator>> &colored_iterators,
          Worker                                    worker,
          Copier                                    copier,
          const ScratchData &                       sample_scratch_data,
          const CopyData &                          sample_copy_data)
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data; // NOLINT

        // Optimization: Check if the functions are not the zero function. To
        // check zero-ness, create a C++ function out of it:
        const bool have_worker =
          (static_cast<const std::function<
             void(const Iterator &, ScratchData &, CopyData &)> &>(worker)) !=
          nullptr;
        const bool have_copier =
          (static_cast<const std::function<void(const CopyData &)> &>(
            copier)) != nullptr;

        // Finally loop over all items and perform the necessary work:
        for (unsigned int color = 0; color < colored_iterators.size(); ++color)
          if (colored_iterators[color].size() > 0)
            for (auto &it : colored_iterators[color])
              {
                if (have_worker)
                  worker(it, scratch_data, copy_data);
                if (have_copier)
                  copier(copy_data);
              }
      }

    } // namespace sequential



#  ifdef DEAL_II_WITH_TBB
    /**
     * 一个用于实现WorkStream模式和功能细节的命名空间。这个命名空间拥有处理Turcksin、Kronbichler和Bangerth的论文中描述的第三个实现的类（见 @ref
     * workstream_paper  ）。
     *
     */
    namespace tbb_colored
    {
      /**
       * 一个结构，包含一个指向从头开始和复制数据对象的指针，以及一个指示该对象当前是否正在使用的标志。
       *
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      struct ScratchAndCopyDataObjects
      {
        std::unique_ptr<ScratchData> scratch_data;
        std::unique_ptr<CopyData>    copy_data;
        bool                         currently_in_use;

        /**
         * 默认构造函数。
         *
         */
        ScratchAndCopyDataObjects()
          : currently_in_use(false)
        {}

        ScratchAndCopyDataObjects(std::unique_ptr<ScratchData> &&p,
                                  std::unique_ptr<CopyData> &&   q,
                                  const bool                     in_use)
          : scratch_data(std::move(p))
          , copy_data(std::move(q))
          , currently_in_use(in_use)
        {}

        // Provide a copy constructor that actually doesn't copy the
        // internal state. This makes handling ScratchAndCopyDataObjects
        // easier to handle with STL containers.
        ScratchAndCopyDataObjects(const ScratchAndCopyDataObjects &)
          : currently_in_use(false)
        {}
      };



      /**
       * 一个管理调用worker和copyer函数的类。与其他实现不同的是，使用parallel_for而不是流水线。
       *
       */
      template <typename Iterator, typename ScratchData, typename CopyData>
      class WorkerAndCopier
      {
      public:
        /**
         * 构造函数。
         *
         */
        WorkerAndCopier(
          const std::function<void(const Iterator &, ScratchData &, CopyData &)>
            &                                          worker,
          const std::function<void(const CopyData &)> &copier,
          const ScratchData &                          sample_scratch_data,
          const CopyData &                             sample_copy_data)
          : worker(worker)
          , copier(copier)
          , sample_scratch_data(sample_scratch_data)
          , sample_copy_data(sample_copy_data)
        {}


        /**
         * 在由两个参数表示的项目范围内调用工作器和复制器函数的函数。
         *
         */
        void
        operator()(const tbb::blocked_range<
                   typename std::vector<Iterator>::const_iterator> &range)
        {
          // we need to find an unused scratch and corresponding copy
          // data object in the list that corresponds to the current
          // thread and then mark it as used. If we can't find one,
          // create one as discussed in the discussion of the documentation
          // of the IteratorRangeToItemStream::scratch_data variable,
          // there is no need to synchronize access to this variable
          // using a mutex as long as we have no yield-point in between.
          // This means that we can't take an iterator into the list
          // now and expect it to still be valid after calling the worker,
          // but we at least do not have to lock the following section.
          ScratchData *scratch_data = nullptr;
          CopyData *   copy_data    = nullptr;
          {
            ScratchAndCopyDataList &scratch_and_copy_data_list = data.get();

            // see if there is an unused object. if so, grab it and mark
            // it as used
            for (typename ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->currently_in_use == false)
                {
                  scratch_data        = p->scratch_data.get();
                  copy_data           = p->copy_data.get();
                  p->currently_in_use = true;
                  break;
                }

            // if no element in the list was found, create one and mark it as
            // used
            if (scratch_data == nullptr)
              {
                Assert(copy_data == nullptr, ExcInternalError());

                scratch_and_copy_data_list.emplace_back(
                  std::make_unique<ScratchData>(sample_scratch_data),
                  std::make_unique<CopyData>(sample_copy_data),
                  true);
                scratch_data =
                  scratch_and_copy_data_list.back().scratch_data.get();
                copy_data = scratch_and_copy_data_list.back().copy_data.get();
              }
          }

          // then call the worker and copier functions on each
          // element of the chunk we were given.
          for (typename std::vector<Iterator>::const_iterator p = range.begin();
               p != range.end();
               ++p)
            {
              try
                {
                  if (worker)
                    worker(*p, *scratch_data, *copy_data);
                  if (copier)
                    copier(*copy_data);
                }
              catch (const std::exception &exc)
                {
                  Threads::internal::handle_std_exception(exc);
                }
              catch (...)
                {
                  Threads::internal::handle_unknown_exception();
                }
            }

          // finally mark the scratch object as unused again. as above, there
          // is no need to lock anything here since the object we work on
          // is thread-local
          {
            ScratchAndCopyDataList &scratch_and_copy_data_list = data.get();

            for (typename ScratchAndCopyDataList::iterator p =
                   scratch_and_copy_data_list.begin();
                 p != scratch_and_copy_data_list.end();
                 ++p)
              if (p->scratch_data.get() == scratch_data)
                {
                  Assert(p->currently_in_use == true, ExcInternalError());
                  p->currently_in_use = false;
                }
          }
        }

      private:
        using ScratchAndCopyDataObjects = typename tbb_colored::
          ScratchAndCopyDataObjects<Iterator, ScratchData, CopyData>;

        /**
         * 对抓取数据对象的列表的类型化定义。这个列表的原理是在使用这些列表的变量中提供的。
         *
         */
        using ScratchAndCopyDataList = std::list<ScratchAndCopyDataObjects>;

        Threads::ThreadLocalStorage<ScratchAndCopyDataList> data;

        /**
         * 指向在单元格序列上进行组装的函数的指针。
         *
         */
        const std::function<void(const Iterator &, ScratchData &, CopyData &)>
          worker;

        /**
         * 指向从本地贡献复制到全局对象的函数的指针。
         *
         */
        const std::function<void(const CopyData &)> copier;

        /**
         * 当我们需要的时候，对样本刮擦和复制数据的引用。
         *
         */
        const ScratchData &sample_scratch_data;
        const CopyData &   sample_copy_data;
      };

      /**
       * 使用TBB的彩色运行函数。
       *
       */
      template <typename Worker,
                typename Copier,
                typename Iterator,
                typename ScratchData,
                typename CopyData>
      void
      run(const std::vector<std::vector<Iterator>> &colored_iterators,
          Worker                                    worker,
          Copier                                    copier,
          const ScratchData &                       sample_scratch_data,
          const CopyData &                          sample_copy_data,
          const unsigned int                        chunk_size)
      {
        // loop over the various colors of what we're given
        for (unsigned int color = 0; color < colored_iterators.size(); ++color)
          if (colored_iterators[color].size() > 0)
            {
              using WorkerAndCopier = internal::tbb_colored::
                WorkerAndCopier<Iterator, ScratchData, CopyData>;

              using RangeType = typename std::vector<Iterator>::const_iterator;

              WorkerAndCopier worker_and_copier(worker,
                                                copier,
                                                sample_scratch_data,
                                                sample_copy_data);

              parallel::internal::parallel_for(
                colored_iterators[color].begin(),
                colored_iterators[color].end(),
                [&worker_and_copier](
                  const tbb::blocked_range<
                    typename std::vector<Iterator>::const_iterator> &range) {
                  worker_and_copier(range);
                },
                chunk_size);
            }
      }

    }    // namespace tbb_colored
#  endif // DEAL_II_WITH_TBB


  } // namespace internal



  /**
   * @note
   * 如果你的数据对象很大，或者它们的构造函数很昂贵，记住<tt>queue_length</tt>拷贝的<tt>ScratchData</tt>对象和<tt>queue_length*chunk_size</tt>拷贝的<tt>CopyData</tt>对象是有帮助的。
   * @note  在复制器不做任何事情的情况下，传递
   * `std::function<void(const  CopyData &)>()`作为 @p copier
   * 以确保内部使用更有效的算法。然而，重要的是要认识到，上面创建的空函数对象并不*
   * 与具有空主体的lambda函数不同，`[](const CopyData &){}`。
   *
   * - 从这个函数的角度来看，没有办法识别作为复制者提供的λ函数在其主体中是否做了什么，所以需要复制。另一方面，一个默认构造的 `std::function` 对象可以*被识别，然后被用来选择一个更有效的算法。
   *
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const std::vector<std::vector<Iterator>> &colored_iterators,
      Worker                                    worker,
      Copier                                    copier,
      const ScratchData &                       sample_scratch_data,
      const CopyData &                          sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8);


  /**
   * 这是WorkStream概念的两个主要功能之一，做本命名空间介绍中描述的工作。它对应于Turcksin、Kronbichler和Bangerth的论文中的实现2（见 @ref
   * workstream_paper ）。
   * 这个函数可以用于worker和copyer对象，这些对象要么是指向非成员函数的指针，要么是允许用operator()调用的对象，例如lambda函数或由
   * std::bind. 创建的对象。
   * 如果copyer是一个空函数，它在管道中会被忽略。然而，一个具有空主体的lambda函数并不*等同于一个空的
   * `std::function` 对象，因此，不会被忽略。    作为 @p end
   * 传递的参数必须可以转换为与 @p
   * 开始时相同的类型，但本身不一定是相同的类型。这允许编写类似<code>WorkStream().run(dof_handler.begin_active(),
   * dof_handler.end(), ...</code>的代码，其中第一个是
   * DoFHandler::active_cell_iterator 类型，而第二个是
   * DoFHandler::raw_cell_iterator. 类型。
   * 两个数据类型<tt>ScratchData</tt>和<tt>CopyData</tt>需要有一个工作拷贝构造器。<tt>ScratchData</tt>只在<tt>worker</tt>函数中使用，而<tt>CopyData</tt>是由<tt>worker</tt>传递给<tt>copier</tt>的对象。
   * @p queue_length
   * 参数表示在任何给定时间内可以活用的项目数量。每个项目由输入流的
   * @p chunk_size
   * 个元素组成，这些元素将被worker和copier函数在同一个线程上一个接一个地处理。
   * @note
   * 如果你的数据对象很大，或者它们的构造函数很昂贵，记住<tt>queue_length</tt>拷贝的<tt>ScratchData</tt>对象和<tt>queue_length*chunk_size</tt>拷贝的<tt>CopyData</tt>对象是有帮助的。
   * @note  在拷贝器不做任何事情的情况下，传递
   * `std::function<void(const  CopyData &)>()`作为 @p copier
   * 以确保内部使用更高效的算法。然而，重要的是要认识到，上面创建的空函数对象并不*
   * 与具有空主体的lambda函数不同，`[](const CopyData &){}`。
   *
   * - 从这个函数的角度来看，没有办法识别作为复制者提供的λ函数在其主体中是否做了什么，所以需要复制。另一方面，一个默认构造的 `std::function` 对象可以*被识别，然后被用来选择一个更有效的算法。
   *
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator &                         begin,
      const typename identity<Iterator>::type &end,
      Worker                                   worker,
      Copier                                   copier,
      const ScratchData &                      sample_scratch_data,
      const CopyData &                         sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    Assert(queue_length > 0,
           ExcMessage("The queue length must be at least one, and preferably "
                      "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert(chunk_size > 0, ExcMessage("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode

    // If no work then skip. (only use operator!= for iterators since we may
    // not have an equality comparison operator)
    if (!(begin != end))
      return;

    if (MultithreadInfo::n_threads() > 1)
      {
#  ifdef DEAL_II_WITH_TBB
        if (static_cast<const std::function<void(const CopyData &)> &>(copier))
          {
            // If we have a copier, run the algorithm:
            internal::tbb_no_coloring::run(begin,
                                           end,
                                           worker,
                                           copier,
                                           sample_scratch_data,
                                           sample_copy_data,
                                           queue_length,
                                           chunk_size);
          }
        else
          {
            // There is no copier function. in this case, we have an
            // embarrassingly parallel problem where we can
            // essentially apply parallel_for. because parallel_for
            // requires subdividing the range for which operator- is
            // necessary between iterators, it is often inefficient to
            // apply it directly to cell ranges and similar iterator
            // types for which operator- is expensive or, in fact,
            // nonexistent. rather, in that case, we simply copy the
            // iterators into a large array and use operator- on
            // iterators to this array of iterators.
            //
            // instead of duplicating code, this is essentially the
            // same situation we have in the colored implementation below, so we
            // just defer to that place
            std::vector<std::vector<Iterator>> all_iterators(1);
            for (Iterator p = begin; p != end; ++p)
              all_iterators[0].push_back(p);

            run(all_iterators,
                worker,
                copier,
                sample_scratch_data,
                sample_copy_data,
                queue_length,
                chunk_size);
          }

        // exit this function to not run the sequential version below:
        return;
#  endif
      }

    // no TBB installed or we are requested to run sequentially:
    internal::sequential::run(
      begin, end, worker, copier, sample_scratch_data, sample_copy_data);
  }



  /**
   * 与上面的函数相同，但用于迭代器范围和C风格的数组。
   * 一个满足迭代器范围要求的类定义了
   * `IteratorRangeType::begin()` 和 `IteratorRangeType::end()`,
   * 两个函数，这两个函数都返回到构成范围边界的元素的迭代器。
   *
   */
  template <typename Worker,
            typename Copier,
            typename IteratorRangeType,
            typename ScratchData,
            typename CopyData,
            typename = typename std::enable_if<
              has_begin_and_end<IteratorRangeType>::value>::type>
  void
  run(IteratorRangeType  iterator_range,
      Worker             worker,
      Copier             copier,
      const ScratchData &sample_scratch_data,
      const CopyData &   sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(iterator_range.begin(),
        iterator_range.end(),
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  /**
   * 与上面的函数相同，但针对deal.II的IteratorRange。
   *
   */
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const IteratorRange<Iterator> &iterator_range,
      Worker                         worker,
      Copier                         copier,
      const ScratchData &            sample_scratch_data,
      const CopyData &               sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(iterator_range.begin(),
        iterator_range.end(),
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const std::vector<std::vector<Iterator>> &colored_iterators,
      Worker                                    worker,
      Copier                                    copier,
      const ScratchData &                       sample_scratch_data,
      const CopyData &                          sample_copy_data,
      const unsigned int                        queue_length,
      const unsigned int                        chunk_size)
  {
    Assert(queue_length > 0,
           ExcMessage("The queue length must be at least one, and preferably "
                      "larger than the number of processors on this system."));
    (void)queue_length; // removes -Wunused-parameter warning in optimized mode
    Assert(chunk_size > 0, ExcMessage("The chunk_size must be at least one."));
    (void)chunk_size; // removes -Wunused-parameter warning in optimized mode


    if (MultithreadInfo::n_threads() > 1)
      {
#  ifdef DEAL_II_WITH_TBB
        internal::tbb_colored::run(colored_iterators,
                                   worker,
                                   copier,
                                   sample_scratch_data,
                                   sample_copy_data,
                                   chunk_size);

        // exit this function to not run the sequential version below:
        return;
#  endif
      }

    // run all colors sequentially:
    {
      internal::sequential::run(colored_iterators,
                                worker,
                                copier,
                                sample_scratch_data,
                                sample_copy_data);
    }
  }



  /**
   * 这是WorkStream概念的两个主要函数之一的变体，做本命名空间介绍中描述的工作。
   * 它对应于Turcksin、Kronbichler和Bangerth的论文中的实现2（见
   * @ref workstream_paper  ）。
   * 这是可以用于作为类的成员函数的工作者和复制者函数的函数。如果复制器是一个空函数，那么它在管道中会被忽略。
   * 作为 @p end 传递的参数必须可以转换为与 @p
   * 开始相同的类型，但本身不一定是同一类型。这允许编写类似<code>WorkStream().run(dof_handler.begin_active(),
   * dof_handler.end(), ...</code>的代码，其中第一个是
   * DoFHandler::active_cell_iterator 类型，而第二个是
   * DoFHandler::raw_cell_iterator. 类型。  @p queue_length
   * 参数表示在任何特定时间可以直播的项目数量。每个项目由输入流的
   * @p chunk_size
   * 个元素组成，这些元素将由工作者和复制者函数在同一线程上一个接一个地处理。
   * @note
   * 如果你的数据对象很大，或者它们的构造函数很昂贵，记住<tt>queue_length</tt>拷贝的<tt>ScratchData</tt>对象和<tt>queue_length*chunk_size</tt>拷贝的<tt>CopyData</tt>对象是有帮助的。
   * @note  在拷贝器不做任何事情的情况下，传递
   * `std::function<void(const  CopyData &)>()`作为 @p copier
   * 以确保内部使用更有效的算法。然而，重要的是要认识到，上面创建的空函数对象并不*
   * 与具有空主体的lambda函数不同，`[](const CopyData &){}`。
   *
   * - 从这个函数的角度来看，没有办法识别作为复制者提供的λ函数在其主体中是否做了什么，所以需要复制。另一方面，一个默认构造的 `std::function` 对象可以*被识别，然后被用来选择一个更有效的算法。
   *
   */
  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator &                         begin,
      const typename identity<Iterator>::type &end,
      MainClass &                              main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData &   sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // forward to the other function
    run(begin,
        end,
        [&main_object, worker](const Iterator &iterator,
                               ScratchData &   scratch_data,
                               CopyData &      copy_data) {
          (main_object.*worker)(iterator, scratch_data, copy_data);
        },
        [&main_object, copier](const CopyData &copy_data) {
          (main_object.*copier)(copy_data);
        },
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }


  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const IteratorOverIterators<Iterator> &                         begin,
      const IteratorOverIterators<typename identity<Iterator>::type> &end,
      MainClass &main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData &   sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // forward to the other function
    run(begin,
        end,
        [&main_object, worker](const Iterator &iterator,
                               ScratchData &   scratch_data,
                               CopyData &      copy_data) {
          (main_object.*worker)(iterator, scratch_data, copy_data);
        },
        [&main_object, copier](const CopyData &copy_data) {
          (main_object.*copier)(copy_data);
        },
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  /**
   * 与上面的函数相同，但用于迭代器范围和C风格的数组。
   * 一个满足迭代器范围要求的类定义了
   * `IteratorRangeType::begin()` 和 `IteratorRangeType::end()`,
   * 两个函数，这两个函数都返回到构成范围边界的元素的迭代器。
   *
   */
  template <typename MainClass,
            typename IteratorRangeType,
            typename ScratchData,
            typename CopyData,
            typename = typename std::enable_if<
              has_begin_and_end<IteratorRangeType>::value>::type>
  void
  run(IteratorRangeType iterator_range,
      MainClass &       main_object,
      void (MainClass::*worker)(
        const typename identity<IteratorRangeType>::type::iterator &,
        ScratchData &,
        CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData &   sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(std::begin(iterator_range),
        std::end(iterator_range),
        main_object,
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }



  /**
   * 与上面的函数相同，但针对deal.II的IteratorRange。
   *
   */
  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(IteratorRange<Iterator> iterator_range,
      MainClass &             main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData &   sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // Call the function above
    run(std::begin(iterator_range),
        std::end(iterator_range),
        main_object,
        worker,
        copier,
        sample_scratch_data,
        sample_copy_data,
        queue_length,
        chunk_size);
  }

} // namespace WorkStream



DEAL_II_NAMESPACE_CLOSE



//----------------------------   work_stream.h     ---------------------------
// end of #ifndef dealii_work_stream_h
#endif
//----------------------------   work_stream.h     ---------------------------


