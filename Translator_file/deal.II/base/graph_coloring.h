//include/deal.II-translator/base/graph_coloring_0.txt

// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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

#ifndef dealii_graph_coloring_h
#  define dealii_graph_coloring_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/thread_management.h>

#  include <deal.II/lac/sparsity_tools.h>

#  include <functional>
#  include <set>
#  include <unordered_map>
#  include <unordered_set>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * 一个包含可以给图形着色的函数的命名空间。
 *
 *
 */
namespace GraphColoring
{
  namespace internal
  {
    /**
     * 给定两组假设为排序的索引，确定它们是否会有一个非空的交集。实际的交集是不计算的。
     * @param  indices1 一组指数，假定已排序。      @param
     * indices2 一组指数，假定已排序。      @return
     * 两组指数是否有一个非空的交集。
     *
     */
    inline bool
    have_nonempty_intersection(
      const std::vector<types::global_dof_index> &indices1,
      const std::vector<types::global_dof_index> &indices2)
    {
      // we assume that both arrays are sorted, so we can walk
      // them in lockstep and see if we encounter an index that's
      // in both arrays. once we reach the end of either array,
      // we know that there is no intersection
      std::vector<types::global_dof_index>::const_iterator p = indices1.begin(),
                                                           q = indices2.begin();
      while ((p != indices1.end()) && (q != indices2.end()))
        {
          if (*p < *q)
            ++p;
          else if (*p > *q)
            ++q;
          else
            // conflict found!
            return true;
        }

      // no conflict found!
      return false;
    }


    /**
     * 使用Cuthill-McKee算法（广度优先搜索算法）的简化版本创建给定范围的迭代器分区。该函数创建的分区包含迭代器的
     * "区域"，其中第一个分区包含第一个迭代器，第二个区域包含所有与第一个区域的单一元素有冲突的迭代器，第三个区域包含那些与第二个区域的迭代器有冲突并且之前没有被分配到一个区域的迭代器，等等。如果迭代器代表单元，那么这就会产生像洋葱壳一样的分区，围绕着最开始的单元。请注意，每个区的元素可能与同一区的其他元素冲突。
     * 两个迭代器是否冲突的问题由一个用户提供的函数决定。这个函数的含义在
     * GraphColoring::make_graph_coloring() 函数的文档中讨论。
     * @param[in]  begin 寻求分区的迭代器范围的第一个元素。
     * @param[in]  end 迭代器范围结束后的元素。      @param[in]
     * get_conflict_indices
     * 一个用户定义的函数对象，返回一组描述代表冲突的指标。更深入的讨论见上文。
     * @return  一组迭代器的集合（为了提高效率，集合由
     * std::vector
     * 表示）。最外层集合的每个元素都对应于指向处于同一分区（即同一区域）的对象的迭代器。
     *
     */
    template <typename Iterator>
    std::vector<std::vector<Iterator>>
    create_partitioning(
      const Iterator &                         begin,
      const typename identity<Iterator>::type &end,
      const std::function<std::vector<types::global_dof_index>(
        const Iterator &)> &                   get_conflict_indices)
    {
      // Number of iterators.
      unsigned int n_iterators = 0;

      // Create a map from conflict indices to iterators
      std::unordered_map<types::global_dof_index, std::vector<Iterator>>
        indices_to_iterators;
      for (Iterator it = begin; it != end; ++it)
        {
          const std::vector<types::global_dof_index> conflict_indices =
            get_conflict_indices(it);
          const unsigned int n_conflict_indices = conflict_indices.size();
          for (unsigned int i = 0; i < n_conflict_indices; ++i)
            indices_to_iterators[conflict_indices[i]].push_back(it);
          ++n_iterators;
        }

      // create the very first zone which contains only the first
      // iterator. then create the other zones. keep track of all the
      // iterators that have already been assigned to a zone
      std::vector<std::vector<Iterator>> zones(1,
                                               std::vector<Iterator>(1, begin));
      std::set<Iterator>                 used_it;
      used_it.insert(begin);
      while (used_it.size() != n_iterators)
        {
          // loop over the elements of the previous zone. for each element of
          // the previous zone, get the conflict indices and from there get
          // those iterators that are conflicting with the current element
          typename std::vector<Iterator>::iterator previous_zone_it(
            zones.back().begin());
          typename std::vector<Iterator>::iterator previous_zone_end(
            zones.back().end());
          std::vector<Iterator> new_zone;
          for (; previous_zone_it != previous_zone_end; ++previous_zone_it)
            {
              const std::vector<types::global_dof_index> conflict_indices =
                get_conflict_indices(*previous_zone_it);

              const unsigned int n_conflict_indices(conflict_indices.size());
              for (unsigned int i = 0; i < n_conflict_indices; ++i)
                {
                  const std::vector<Iterator> &conflicting_elements =
                    indices_to_iterators[conflict_indices[i]];
                  for (unsigned int j = 0; j < conflicting_elements.size(); ++j)
                    {
                      // check that the iterator conflicting with the current
                      // one is not associated to a zone yet and if so, assign
                      // it to the current zone. mark it as used
                      //
                      // we can shortcut this test if the conflicting iterator
                      // is the current iterator
                      if ((conflicting_elements[j] != *previous_zone_it) &&
                          (used_it.count(conflicting_elements[j]) == 0))
                        {
                          new_zone.push_back(conflicting_elements[j]);
                          used_it.insert(conflicting_elements[j]);
                        }
                    }
                }
            }

          // If there are iterators in the new zone, then the zone is added to
          // the partition. Otherwise, the graph is disconnected and we need to
          // find an iterator on the other part of the graph. start the whole
          // process again with the first iterator that hasn't been assigned to
          // a zone yet
          if (new_zone.size() != 0)
            zones.push_back(new_zone);
          else
            for (Iterator it = begin; it != end; ++it)
              if (used_it.count(it) == 0)
                {
                  zones.push_back(std::vector<Iterator>(1, it));
                  used_it.insert(it);
                  break;
                }
        }

      return zones;
    }



    /**
     * 这个函数使用DSATUR（Degree
     * SATURation）来给一个集合的元素着色。DSATUR的工作原理如下。
     *
     *
     *
     *
     *
     * - 按度数递减的顺序排列顶点。
     *
     *
     * - 用颜色1给最大程度的顶点上色。
     *
     *
     *
     * - 选择一个具有最大饱和度的顶点。如果存在平等，则选择无色子图中的任何最大程度的顶点。
     *
     *
     *
     *
     *
     * - 用最少的可能（最低编号）的颜色给所选的顶点上色。
     *
     *
     *
     *
     *
     * - 如果所有的顶点都被着色，则停止。否则，返回到3。  @param[in]  partition 应该被着色的迭代器的集合。      @param[in]  get_conflict_indices 一个用户定义的函数对象，返回一组描述代表冲突的指标。更深入的讨论见上文。      @param[out]  partition_coloring 一组迭代器的集合（其中集合由 std::vector 表示，以提高效率）。最外层集合的每个元素都对应于指向处于同一分区（具有相同颜色）的对象的迭代器，因此不会冲突。不同集合中的元素可能会发生冲突。
     *
     */
    template <typename Iterator>
    void
    make_dsatur_coloring(
      std::vector<Iterator> &             partition,
      const std::function<std::vector<types::global_dof_index>(
        const Iterator &)> &              get_conflict_indices,
      std::vector<std::vector<Iterator>> &partition_coloring)
    {
      partition_coloring.clear();

      // Number of zones composing the partitioning.
      const unsigned int        partition_size(partition.size());
      std::vector<unsigned int> sorted_vertices(partition_size);
      std::vector<int>          degrees(partition_size);
      std::vector<std::vector<types::global_dof_index>> conflict_indices(
        partition_size);
      std::vector<std::vector<unsigned int>> graph(partition_size);

      // Get the conflict indices associated to each iterator. The
      // conflict_indices have to be sorted so we can more easily find conflicts
      // later on
      for (unsigned int i = 0; i < partition_size; ++i)
        {
          conflict_indices[i] = get_conflict_indices(partition[i]);
          std::sort(conflict_indices[i].begin(), conflict_indices[i].end());
        }

      // Compute the degree of each vertex of the graph using the
      // intersection of the conflict indices.
      for (unsigned int i = 0; i < partition_size; ++i)
        for (unsigned int j = i + 1; j < partition_size; ++j)
          // If the two iterators share indices then we increase the degree of
          // the vertices and create an ''edge'' in the graph.
          if (have_nonempty_intersection(conflict_indices[i],
                                         conflict_indices[j]))
            {
              ++degrees[i];
              ++degrees[j];
              graph[i].push_back(j);
              graph[j].push_back(i);
            }

      // Sort the vertices by decreasing degree.
      std::vector<int>::iterator degrees_it;
      for (unsigned int i = 0; i < partition_size; ++i)
        {
          // Find the largest element.
          degrees_it         = std::max_element(degrees.begin(), degrees.end());
          sorted_vertices[i] = degrees_it - degrees.begin();
          // Put the largest element to -1 so it cannot be chosen again.
          *degrees_it = -1;
        }

      // Color the graph.
      std::vector<std::unordered_set<unsigned int>> colors_used;
      for (unsigned int i = 0; i < partition_size; ++i)
        {
          const unsigned int current_vertex(sorted_vertices[i]);
          bool               new_color(true);
          // Try to use an existing color, i.e., try to find a color which is
          // not associated to one of the vertices linked to current_vertex.
          // Loop over the color.
          for (unsigned int j = 0; j < partition_coloring.size(); ++j)
            {
              // Loop on the vertices linked to current_vertex. If one vertex
              // linked to current_vertex is already using the color j, this
              // color cannot be used anymore.
              bool unused_color(true);
              for (const auto adjacent_vertex : graph[current_vertex])
                if (colors_used[j].count(adjacent_vertex) == 1)
                  {
                    unused_color = false;
                    break;
                  }
              if (unused_color)
                {
                  partition_coloring[j].push_back(partition[current_vertex]);
                  colors_used[j].insert(current_vertex);
                  new_color = false;
                  break;
                }
            }
          // Add a new color.
          if (new_color)
            {
              partition_coloring.push_back(
                std::vector<Iterator>(1, partition[current_vertex]));
              std::unordered_set<unsigned int> tmp;
              tmp.insert(current_vertex);
              colors_used.push_back(tmp);
            }
        }
    }



    /**
     * 给定一个分区着色图，即一组分区（partition），每个分区都有颜色，为整个迭代器集产生一个组合着色。这是可能的，因为一个偶数（或奇数）区的任何颜色都不会与任何其他偶数（或奇数）区的任何颜色冲突。因此，我们可以将所有偶数区和奇数区的颜色结合起来。这个函数试图创建元素数量相似的颜色。
     *
     */
    template <typename Iterator>
    std::vector<std::vector<Iterator>>
    gather_colors(
      const std::vector<std::vector<std::vector<Iterator>>> &partition_coloring)
    {
      std::vector<std::vector<Iterator>> coloring;

      // Count the number of iterators in each color.
      const unsigned int partition_size(partition_coloring.size());
      std::vector<std::vector<unsigned int>> colors_counter(partition_size);
      for (unsigned int i = 0; i < partition_size; ++i)
        {
          const unsigned int n_colors(partition_coloring[i].size());
          colors_counter[i].resize(n_colors);
          for (unsigned int j = 0; j < n_colors; ++j)
            colors_counter[i][j] = partition_coloring[i][j].size();
        }

      // Find the partition with the largest number of colors for the even
      // partition.
      unsigned int       i_color(0);
      unsigned int       max_even_n_colors(0);
      const unsigned int colors_size(colors_counter.size());
      for (unsigned int i = 0; i < colors_size; i += 2)
        {
          if (max_even_n_colors < colors_counter[i].size())
            {
              max_even_n_colors = colors_counter[i].size();
              i_color           = i;
            }
        }
      coloring.resize(max_even_n_colors);
      for (unsigned int j = 0; j < colors_counter[i_color].size(); ++j)
        coloring[j] = partition_coloring[i_color][j];

      for (unsigned int i = 0; i < partition_size; i += 2)
        {
          if (i != i_color)
            {
              std::unordered_set<unsigned int> used_k;
              for (unsigned int j = 0; j < colors_counter[i].size(); ++j)
                {
                  // Find the color in the current partition with the largest
                  // number of iterators.
                  std::vector<unsigned int>::iterator it;
                  it = std::max_element(colors_counter[i].begin(),
                                        colors_counter[i].end());
                  unsigned int min_iterators(static_cast<unsigned int>(-1));
                  unsigned int pos(0);
                  // Find the color of coloring with the least number of colors
                  // among the colors that have not been used yet.
                  for (unsigned int k = 0; k < max_even_n_colors; ++k)
                    if (used_k.count(k) == 0)
                      if (colors_counter[i_color][k] < min_iterators)
                        {
                          min_iterators = colors_counter[i_color][k];
                          pos           = k;
                        }
                  colors_counter[i_color][pos] += *it;
                  // Concatenate the current color with the existing coloring.
                  coloring[pos].insert(
                    coloring[pos].end(),
                    partition_coloring[i][it - colors_counter[i].begin()]
                      .begin(),
                    partition_coloring[i][it - colors_counter[i].begin()]
                      .end());
                  used_k.insert(pos);
                  // Put the number of iterators to the current color to zero.
                  *it = 0;
                }
            }
        }

      // If there is more than one partition, do the same thing that we did for
      // the even partitions to the odd partitions
      if (partition_size > 1)
        {
          unsigned int max_odd_n_colors(0);
          for (unsigned int i = 1; i < partition_size; i += 2)
            {
              if (max_odd_n_colors < colors_counter[i].size())
                {
                  max_odd_n_colors = colors_counter[i].size();
                  i_color          = i;
                }
            }
          coloring.resize(max_even_n_colors + max_odd_n_colors);
          for (unsigned int j = 0; j < colors_counter[i_color].size(); ++j)
            coloring[max_even_n_colors + j] = partition_coloring[i_color][j];

          for (unsigned int i = 1; i < partition_size; i += 2)
            {
              if (i != i_color)
                {
                  std::unordered_set<unsigned int> used_k;
                  for (unsigned int j = 0; j < colors_counter[i].size(); ++j)
                    {
                      // Find the color in the current partition with the
                      // largest number of iterators.
                      std::vector<unsigned int>::iterator it;
                      it = std::max_element(colors_counter[i].begin(),
                                            colors_counter[i].end());
                      unsigned int min_iterators(static_cast<unsigned int>(-1));
                      unsigned int pos(0);
                      // Find the color of coloring with the least number of
                      // colors among the colors that have not been used yet.
                      for (unsigned int k = 0; k < max_odd_n_colors; ++k)
                        if (used_k.count(k) == 0)
                          if (colors_counter[i_color][k] < min_iterators)
                            {
                              min_iterators = colors_counter[i_color][k];
                              pos           = k;
                            }
                      colors_counter[i_color][pos] += *it;
                      // Concatenate the current color with the existing
                      // coloring.
                      coloring[max_even_n_colors + pos].insert(
                        coloring[max_even_n_colors + pos].end(),
                        partition_coloring[i][it - colors_counter[i].begin()]
                          .begin(),
                        partition_coloring[i][it - colors_counter[i].begin()]
                          .end());
                      used_k.insert(pos);
                      // Put the number of iterators to the current color to
                      // zero.
                      *it = 0;
                    }
                }
            }
        }

      return coloring;
    }
  } // namespace internal


  /**
   * 对给定范围的迭代器创建一个分区，这样指向冲突对象的迭代器将被放入不同的分区，其中两个对象是否冲突的问题由用户提供的函数决定。
   * 这个函数也可以看作是一个图的着色：迭代器所指向的每个对象被认为是一个节点，每两个冲突的节点之间有一条边。然后，图形着色算法为每个节点分配一个颜色，使由一条边连接的两个节点不具有相同的颜色。
   * 这个函数的一个典型的用例是在并行组装一个矩阵。在这里，人们希望同时在不同的单元上集合局部贡献（这种操作是纯粹的局部操作，因此不需要同步），但随后我们需要将这些局部贡献添加到全局矩阵中。一般来说，如果单元共享自由度，来自不同单元的贡献可能是对同一矩阵项的贡献，因此，除非我们想冒竞赛条件的风险，否则不能在同一时间发生（见http://en.wikipedia.org/wiki/Race_condition）。因此，我们称这两个单元为冲突单元，我们只能允许来自不冲突的单元的并行操作。换句话说，如果矩阵条目集（例如以行为特征的）有一个非空的交集，那么两个单元就处于冲突之中。
   * 在这种情况下，计算冲突图需要调用一个确定两个迭代器（或它们代表的两个对象）是否冲突的函数，并为每一对迭代器调用，即
   * $\frac 12 N (N-1)$
   * 次。这在一般情况下是太昂贵了。一个更好的方法是要求一个用户定义的函数，为每个被调用的迭代器返回一个表征冲突的某种指标集；如果两个迭代器的冲突指标集有一个非空的交集，那么它们就是冲突的。在组装矩阵的例子中，冲突指标集将包含指向的单元上所有自由度的指数（在连续Galerkin方法的情况下），或者当前单元和与当前单元面相邻的所有单元上自由度指数的联合（在非连续Galerkin方法的情况下，因为在那里计算面积分，耦合由共同面连接的自由度
   *
   * - 见 step-12 ）。)
   * @note
   * 由作为第三个参数传递的用户定义函数返回的冲突集需要准确地描述<i>all</i>自由度，对于这些自由度，有什么东西被写进矩阵或右手边。换句话说，如果写入是通过
   * AffineConstraints::copy_local_to_global(),
   * 这样的函数发生的，那么冲突指标集实际上不仅要包含当前单元格上的自由度，还要包含它们通过约束条件（如悬挂节点）连接的自由度。
   * 在其他情况下，冲突指标集可能代表完全不同的东西
   *
   * - 由这个函数的调用者来描述两个迭代器冲突的含义。鉴于此，计算冲突图边可以比用 ${\cal O}(N^2)$ 操作大大便宜。    在任何情况下，该函数的结果将是，其冲突指标集有重叠的迭代器不会被分配到相同的颜色。
   * @note
   * 这个函数中使用的算法在Turcksin、Kronbichler和Bangerth的论文中有所描述，见
   * @ref workstream_paper  。      @param[in]  开始
   * 寻求着色的迭代器范围的第一个元素。    @param[in]  结束
   * 迭代器范围结束后的元素。    @param[in]  get_conflict_indices
   * 一个用户定义的函数对象，返回一组代表冲突的描述性指标。
   * 更深入的讨论见上文。    @return
   * 一组迭代器的集合（为了提高效率，集合由 std::vector
   * 表示）。最外层集合的每个元素都对应于指向处于同一分区（具有相同颜色）的对象的迭代器，因此不会冲突。不同集合中的元素可能会发生冲突。
   *
   */
  template <typename Iterator>
  std::vector<std::vector<Iterator>>
  make_graph_coloring(
    const Iterator &                               begin,
    const typename identity<Iterator>::type &      end,
    const std::function<std::vector<types::global_dof_index>(
      const typename identity<Iterator>::type &)> &get_conflict_indices)
  {
    Assert(begin != end,
           ExcMessage(
             "GraphColoring is not prepared to deal with empty ranges!"));

    // Create the partitioning.
    std::vector<std::vector<Iterator>> partitioning =
      internal::create_partitioning(begin, end, get_conflict_indices);

    // Color the iterators within each partition.
    // Run the coloring algorithm on each zone in parallel
    const unsigned int partitioning_size(partitioning.size());
    std::vector<std::vector<std::vector<Iterator>>> partition_coloring(
      partitioning_size);

    Threads::TaskGroup<> tasks;
    for (unsigned int i = 0; i < partitioning_size; ++i)
      tasks += Threads::new_task(&internal::make_dsatur_coloring<Iterator>,
                                 partitioning[i],
                                 get_conflict_indices,
                                 partition_coloring[i]);
    tasks.join_all();

    // Gather the colors together.
    return internal::gather_colors(partition_coloring);
  }

  /**
   * GraphColoring::color_sparsity_pattern, 是
   * SparsityTools::color_sparsity_pattern,
   * 的一个包装函数，是使用SparsityPattern表示的图连接进行着色的另一种方法。
   * 进一步的细节，请参考 SparsityTools::color_sparsity_pattern. 。
   *
   */
  unsigned int
  color_sparsity_pattern(const SparsityPattern &    sparsity_pattern,
                         std::vector<unsigned int> &color_indices);

} // namespace GraphColoring

DEAL_II_NAMESPACE_CLOSE


//----------------------------   graph_coloring.h ---------------------------
// end of #ifndef dealii_graph_coloring_h
#endif
//----------------------------   graph_coloring.h ---------------------------


