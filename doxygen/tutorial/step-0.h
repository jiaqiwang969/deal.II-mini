/**
@page step_0 The step-0 tutorial program
@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Aboutthetutorial"> About the tutorial </a>
        <li><a href="#Videolecturesontutorialprograms"> Video lectures on tutorial programs </a>
        <li><a href="#Whatthisprogramdoes"> What this program does </a>
        <li><a href="#Aboutscientificcomputingingeneral"> About scientific computing in general </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Creatingthefirstmesh">Creating the first mesh</a>
        <li><a href="#Creatingthesecondmesh">Creating the second mesh</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibleextensions"> Possible extensions </a>
      <ul>
        <li><a href="#Differentadaptiverefinementstrategies"> Different adaptive refinement strategies </a>
        <li><a href="#Differentgeometries"> Different geometries </a>
        <li><a href="#Commentsaboutprogramminganddebugging"> Comments about programming and debugging </a>
        <li><a href="#Moreaboutgraphicaloutput"> More about graphical output </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="Intro"></a>
<a name="Introduction"></a><h1>Introduction</h1>


<a name="Aboutthetutorial"></a><h3> About the tutorial </h3>


Since this is the first tutorial program, let us comment first on how
this tutorial and the rest of the deal.II documentation is supposed to
work. The documentation for deal.II comes essentially at three
different levels:
- The tutorial: This is a collection of programs that shows how
  deal.II is used in practice. It doesn't typically discuss individual
  functions at the level of individual arguments, but rather wants to
  give the big picture of how things work together. In other words, it
  discusses "concepts": what are the building blocks of deal.II and
  how are they used together in finite element programs.
- The manual: This is the documentation of every single class and
  every single (member) function in deal.II. You get there if, for
  example, you click on the "Main page" or "Classes" tab at the top of
  this page. This is the place where you would look up what the second
  argument of Triangulation::create_triangulation_compatibility means,
  to give just one slightly obscure example. You need this level of
  documentation for when you know what you want to do, but forgot how
  exactly the function was named, what its arguments are, or what it
  returns. Note that you also get into the manual whenever you read
  through the tutorial and click on any of the class or function
  names, i.e. the tutorial contains a great many links into the manual
  for whenever you need a more detailed description of a function or
  class. On the other hand, the manual is not a good place to learn
  deal.II since it gives you a microscopic view of things without
  telling you how a function might fit into the bigger picture.
- Modules: These are groups of classes and functions that work
  together or have related functionality. If you click on the
  "Modules" tab at the top of this page, you end up on a page that
  lists a number of such groups. Each module discusses the underlying
  principles of these classes; for example, the @ref Sparsity module
  talks about all sorts of different issues related to storing
  sparsity patterns of matrices. This is documentation at an
  intermediate level: they give you an overview of what's there in a
  particular area. For example when you wonder what finite element
  classes exist, you would take a look at the @ref fe module. The
  modules are, of course, also cross-linked to the manual (and, at
  times, to the tutorial); if you click on a class name, say on
  Triangulation, would will also at the very top right under the class
  name get a link to the modules this class is a member of if you want
  to learn more about its context.

Let's come back to the tutorial, since you are looking at the first program
(or "step") of it. Each tutorial program is subdivided into the following
sections:
<ol>
  <li> <b>Introduction:</b> This is a discussion of what the program
       does, including the mathematical model, and
       what programming techniques are new compared to previous
       tutorial programs.
  <li> <b>The commented program:</b> An extensively documented listing of the
       source code. Here, we often document individual lines, or
       blocks of code, and discuss what they do, how they do it, and
       why. The comments frequently reference the introduction,
       i.e. you have to understand <i>what</i> the program wants to achieve
       (a goal discussed in the introduction) before you can
       understand <i>how</i> it intends to get there.
  <li> <b>Results:</b> The output of the program, with comments and
       interpretation. This section also frequently has a subsection
       that gives suggestions on how to extend the program in various
       direction; in the earlier programs, this is intended to give
       you directions for little experiments designed to make your
       familiar with deal.II, while in later programs it is more about
       how to use more advanced numerical techniques.
  <li> <b>The plain program:</b> The source code stripped of
       all comments. This is useful if you want to see the "big
       picture" of the code, since the commented version of the
       program has so much text in between that it is often difficult
       to see the entire code of a single function on the screen at
       once.
</ol>

The tutorials are not only meant to be static documentation, but you
should play with them. To this end, go to the
<code>examples/step-1</code> directory (or whatever the number of the
tutorial is that you're interested in) and type
@code
  cmake .
  make
  make run
@endcode
The first command sets up the files that describe which include files this
tutorial program depends on, how to compile it and how to run it. This command
should find the installed deal.II libraries as well that were generated when
you compiled and installed everything as described in the
<a href="../../readme.html" target="body">README</a> file.
If this command should fail to find the deal.II library, then you need to
provide the path to the installation using the command
@code
  cmake -DDEAL_II_DIR=/path/to/installed/deal.II .
@endcode
instead.

The second of the commands above compiles the sources into an executable, while the
last one executes it (strictly speaking, <code>make run</code> will also
compile the code if the executable doesn't exist yet, so you could
have skipped the second command if you wanted). This is all that's
needed to run the code and produce the output that is discussed in the
"Results" section of the tutorial programs. This sequence needs to be repeated
in all of the tutorial directories you want to play with.

When learning the library, you need to play with it and see what
happens. To this end, open the <code>examples/step-1/step-1.cc</code>
source file with your favorite editor and modify it in some way, save it and
run it as above. A few suggestions for possibly modifications are given at the
end of the results section of this program, where we also provide a few links
to other useful pieces of information.


<a name="Videolecturesontutorialprograms"></a><h3> Video lectures on tutorial programs </h3>


This and several of the other tutorial programs are also discussed and
demonstrated in <a
href="http://www.math.colostate.edu/~bangerth/videos.html">Wolfgang
Bangerth's video lectures</a> on deal.II and computational science. In
particular, you can see the steps he executes to run this and other
programs, and you will get a much better idea of the tools that can be
used to work with deal.II. In particular, lectures 2 and 4 give an overview of
deal.II and of the building blocks of any finite element code.
(@dealiiVideoLectureSeeAlso{2,4})

If you are not yet familiar with using Linux and running things on the
command line, you may be interested in watching lectures 2.9 and 2.91.
(@dealiiVideoLectureSeeAlso{2.9,2.91}) These give overviews over the command
line and on what happens when compiling programs, respectively.

Note that deal.II is actively developed, and in the course of this
development we occasionally rename or deprecate functions or classes
that are still referenced in these video lectures.  For
example, the step-1 code shown in video lecture 5 uses a class
HyperShellBoundary which was replaced with SphericalManifold class
later on. Additionally, as of deal.II version 9.0, GridGenerator::hyper_shell()
now automatically attaches a SphericalManifold to the Triangulation. Otherwise
the rest of the lecture material is relevant.

<a name="Whatthisprogramdoes"></a><h3> What this program does </h3>


Let's come back to step-1, the current program.
In this first example, we don't actually do very much, but show two
techniques: what is the syntax to generate triangulation objects, and
some elements of simple loops over all cells. We create two grids, one
which is a regularly refined square (not very exciting, but a common
starting grid for some problems), and one more geometric attempt: a
ring-shaped domain, which is refined towards the inner edge. Through
this, you will get to know three things every finite element program
will have to have somewhere: An object of type Triangulation for the
mesh; a call to the GridGenerator functions to generate a mesh; and
loops over all cells that involve iterators (iterators are a
generalization of pointers and are frequently used in the C++ standard
library; in the context of deal.II, the @ref Iterators module talks
about them).

The program is otherwise small enough that it doesn't need a whole lot
of introduction.

@dealiiVideoLecture{5,6}


<a name="Aboutscientificcomputingingeneral"></a><h3> About scientific computing in general </h3>


If you are reading through this tutorial program, chances are that you are
interested in continuing to use deal.II for your own projects. Thus, you are
about to embark on an exercise in programming using a large-scale scientific
computing library. Unless you are already an experienced user of large-scale
programming methods, this may be new territory for you &mdash; with all the
new rules that go along with it such as the fact that you will have to deal
with code written by others, that you may have to think about documenting your
own code because you may not remember what exactly it is doing a year down the
road (or because others will be using it as well), or coming up with ways to
test that your program is doing the right thing. None of this is something
that we typically train mathematicians, engineers, or scientists in but that
is important when you start writing software of more than a few hundred
lines. Remember: Producing software is not the same as just writing code.

To make your life easier on this journey let us point to some resources that
are worthwhile browsing through before you start any large-scale programming:

- The <a href="https://github.com/dealii/dealii/wiki/Frequently-Asked-Questions">
  deal.II FAQ</a> has a good number of answers to questions about
  particular aspects of deal.II, but also to more general questions such as "How
  do I debug scientific computing codes?" or "Can I train myself to write code
  that has fewer bugs?".

- You will benefit from becoming a better programmer. An excellent
  resource to this end is the book
  [Code Complete](https://en.wikipedia.org/wiki/Code_Complete)
  by Steve McConnell @cite CodeComplete . It's already
  a few years old, with the last edition published in 2004, but it has
  lost none of its appeal as a guide to good programming practices,
  and some of the principal developers use it as a group reading
  project with every generation of their research group members.

- The <a href="http://software-carpentry.org/">Software Carpentry project</a>
  that provides introductions to many topics that are important to dealing
  with software, such as version control, make files, testing, etc. It is
  specifically written for scientists and engineers, not for computer
  scientists, and has a focus on short, practical lessons.

- The <a href="https://bssw.io/">Better Scientific Software
  project</a> has a lot of resources (and interesting blog posts) that
  cover many aspects of writing scientific software.

- The <a href="https://ideas-productivity.org/">IDEAS
  project</a> also has resources on software development, in
  particular for parallel computing. In the "Events" section on
  that site are recorded tutorials and webinars that cover many
  interesting topics.

- An article on <a href="http://arxiv.org/abs/1210.0530">Best
  Practices for Scientific Computing</a> that gives an introduction to
  many of the ways by which you can make sure you are an efficient
  programmer writing programs that work.

As a general recommendation: If you expect to spend more than a few days
writing software in the future, do yourself the favor of learning tools that
can make your life more productive, in particular debuggers and integrated
development environments. (@dealiiVideoLectureSeeAlso{7,8,8.01,25})
You will find that you will get the time spent
learning these tools back severalfold soon by being more productive!
Several of the video lectures referenced above show how to use tools
such as integrated development environments or debuggers.
 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * The most fundamental class in the library is the Triangulation class, which
 * is declared here:
 * #include <deal.II/grid/tria.h>
 * We need the following two includes for loops over cells and/or faces:
 * #include <deal.II/grid/tria_accessor.h>
 * #include <deal.II/grid/tria_iterator.h>
 * Here are some functions to generate standard grids:
 * #include <deal.II/grid/grid_generator.h>
 * Output of grids in various graphics formats:
 * #include <deal.II/grid/grid_out.h>
 * 

 * 
 * 
 * @code
 * #include <list>
 * #include <deal.II/mystl/list.h>
 * 
 * @endcode
 * 
 * This is needed for C++ output:
 * 
 * @code
 * #include <iostream>
 * #include <fstream>
 * @endcode
 * 
 * And this for the declarations of the `std::sqrt` and `std::fabs` functions:
 * 
 * @code
 * #include <cmath>
 * 
 * @endcode
 * 
 * 省略书写std命名空间
 * using namespace std;
 * 

 * 
 * The final step in importing deal.II is this: All deal.II functions and
 * classes are in a namespace <code>dealii</code>, to make sure they don't
 * clash with symbols from other libraries you may want to use in conjunction
 * with deal.II. One could use these functions and classes by prefixing every
 * use of these names by <code>dealii::</code>, but that would quickly become
 * cumbersome and annoying. Rather, we simply import the entire deal.II
 * namespace for general use:
 * using namespace dealii;
 * 

 * 
 * 
 * <a name="Creatingthefirstmesh"></a> 
 * <h3>Creating the first mesh</h3>
 * 

 * 
 * In the following, first function, we simply use the unit square as domain
 * and produce a globally refined grid from it.
 * 
 * @code
 * void first_grid()
 * {
 * @endcode
 * 
 * The first thing to do is to define an object for a triangulation of a
 * two-dimensional domain:
 * Triangulation<2> triangulation;
 * Here and in many following cases, the string "<2>" after a class name
 * indicates that this is an object that shall work in two space
 * dimensions. Likewise, there are versions of the triangulation class that
 * are working in one ("<1>") and three ("<3>") space dimensions. The way
 * this works is through some template magic that we will investigate in
 * some more detail in later example programs; there, we will also see how
 * to write programs in an essentially dimension independent way.
 * 

 * 
 * Next, we want to fill the triangulation with a single cell for a square
 * domain. The triangulation is the refined four times, to yield $4^4=256$
 * cells in total:
 * GridGenerator::hyper_cube(triangulation);
 * triangulation.refine_global(4);
 * 

 * 
 * Now we want to write a graphical representation of the mesh to an output
 * file. The GridOut class of deal.II can do that in a number of different
 * output formats; here, we choose scalable vector graphics (SVG) format
 * that you can visualize using the web browser of your choice:
 * std::ofstream out("grid-1.svg");
 * GridOut       grid_out;
 * grid_out.write_svg(triangulation, out);
 * 
 * @code
 *   std::cout << "TEST PASS" << std::endl;
 * }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Creatingthesecondmesh"></a> 
 * <h3>Creating the second mesh</h3>
 * 

 * 
 * The grid in the following, second function is slightly more complicated in
 * that we use a ring domain and refine the result once globally.
 * void second_grid()
 * {
 * We start again by defining an object for a triangulation of a
 * two-dimensional domain:
 * Triangulation<2> triangulation;
 * 

 * 
 * We then fill it with a ring domain. The center of the ring shall be the
 * point (1,0), and inner and outer radius shall be 0.5 and 1. The number of
 * circumferential cells could be adjusted automatically by this function,
 * but we choose to set it explicitly to 10 as the last argument:
 * const Point<2> center(1, 0);
 * const double   inner_radius = 0.5, outer_radius = 1.0;
 * GridGenerator::hyper_shell(
 * triangulation, center, inner_radius, outer_radius, 10);
 * By default, the triangulation assumes that all boundaries are straight
 * lines, and all cells are bi-linear quads or tri-linear hexes, and that
 * they are defined by the cells of the coarse grid (which we just
 * created). Unless we do something special, when new points need to be
 * introduced the domain is assumed to be delineated by the straight
 * lines of the coarse mesh, and new points will simply be in the middle
 * of the surrounding ones. Here, however, we know that the domain is
 * curved, and we would like to have the Triangulation place new points
 * according to the underlying geometry. Fortunately, some good soul
 * implemented an object which describes a spherical domain, of which the
 * ring is a section; it only needs the center of the ring and
 * automatically figures out how to instruct the Triangulation where to
 * place the new points. The way this works in deal.II is that you tag
 * parts of the triangulation you want to be curved with a number that is
 * usually referred to as "manifold indicator" and then tell the
 * triangulation to use a particular "manifold object" for all places
 * with this manifold indicator. How exactly this works is not important
 * at this point (you can read up on it in step-53 and @ref manifold).
 * The functions in GridGenerator handle this for us in most
 * circumstances: they attach the correct manifold to a domain so that
 * when the triangulation is refined new cells are placed in the correct
 * places. In the present case GridGenerator::hyper_shell attaches a
 * SphericalManifold to all cells: this causes cells to be refined with
 * calculations in spherical coordinates (so new cells have edges that
 * are either radial or lie along concentric circles around the origin).
 * 

 * 
 * By default (i.e., for a Triangulation created by hand or without a
 * call to a GridGenerator function like GridGenerator::hyper_shell or
 * GridGenerator::hyper_ball), all cells and faces of the Triangulation
 * have their manifold_id set to numbers::flat_manifold_id, which is
 * the default if you want a manifold that produces straight edges, but
 * you can change this number for individual cells and faces. In that
 * case, the curved manifold thus associated with number zero will not
 * apply to those parts with a non-zero manifold indicator, but other
 * manifold description objects can be associated with those non-zero
 * indicators. If no manifold description is associated with a particular
 * manifold indicator, a manifold that produces straight edges is
 * implied. (Manifold indicators are a slightly complicated topic; if
 * you're confused about what exactly is happening here, you may want to
 * look at the @ref GlossManifoldIndicator "glossary entry on this
 * topic".) Since the default chosen by GridGenerator::hyper_shell is
 * reasonable we leave things alone.
 * 

 * 
 * In order to demonstrate how to write a loop over all cells, we will
 * refine the grid in five steps towards the inner circle of the domain:
 * for (unsigned int step = 0; step < 5; ++step)
 * {
 * Next, we need to loop over the active cells of the triangulation. You
 * can think of a triangulation as a collection of cells. If it were an
 * array, you would just get a pointer that you increment from one
 * element to the next using the operator `++`. The cells of a
 * triangulation aren't stored as a simple array, but the concept of an
 * <i>iterator</i> generalizes how pointers work to arbitrary collections
 * of objects (see <a href=
 * "http://en.wikipedia.org/wiki/Iterator#C.2B.2B">wikipedia</a> for more
 * information). Typically, any container type in C++ will return an
 * iterator pointing to the start of the collection with a method called
 * `begin`, and an iterator point to 1 past the end of the collection with
 * a method called `end`. We can increment an iterator `it` with the
 * operator `++it`, dereference it to get the underlying data with `*it`,
 * and check to see if we're done by comparing `it != collection.end()`.
 * 

 * 
 * The second important piece is that we only need the active cells.
 * Active cells are those that are not further refined, and the only
 * ones that can be marked for further refinement. deal.II provides
 * iterator categories that allow us to iterate over <i>all</i> cells
 * (including the parent cells of active ones) or only over the active
 * cells. Because we want the latter, we need to call the method
 * Triangulation::active_cell_iterators().
 * 

 * 
 * Putting all of this together, we can loop over all the active cells of
 * a triangulation with
 * <div class=CodeFragmentInTutorialComment>
 * @code{.cpp}
 *     for (auto it = triangulation.active_cell_iterators().begin();
 *          it != triangulation.active_cell_iterators().end();
 *          ++it)
 *       {
 *         auto cell = *it;
 *         // Then a miracle occurs...
 *       }
 * @endcode
 * </div>
 * In the initializer of this loop, we've used the `auto` keyword for the
 * type of the iterator `it`. The `auto` keyword means that the type of
 * the object being declared will be inferred from the context. This
 * keyword is useful when the actual type names are long or possibly even
 * redundant. If you're unsure of what the type is and want to look up
 * what operations the result supports, you can go to the documentation
 * for the method Triangulation::active_cell_iterators(). In this case,
 * the type of `it` is `Triangulation::active_cell_iterator`.
 * 

 * 
 * While the `auto` keyword can save us from having to type out long names
 * of data types, we still have to type a lot of redundant declarations
 * about the start and end iterator and how to increment it. Instead of
 * doing that, we'll use
 * <a href="http://en.cppreference.com/w/cpp/language/range-for">range-
 * based for loops</a>, which wrap up all of the syntax shown above into a
 * much shorter form:
 * for (auto &cell : triangulation.active_cell_iterators())
 * {
 * @note See @ref Iterators for more information about the iterator
 * classes used in deal.II, and @ref CPP11 for more information about
 * range-based for loops and the `auto` keyword.
 * 

 * 
 * Next, we loop over all vertices of the cells. For that purpose
 * we query an iterator over the vertex indices (in 2d, this is an
 * array that contains the elements `{0,1,2,3}`, but since
 * `cell->vertex_indices()` knows the dimension the cell lives in, the
 * array so returned is correct in all dimensions and this enables
 * this code to be correct whether we run it in 2d or 3d, i.e., it
 * enables "dimension-independent programming" -- a big part of what
 * we will discuss in step-4).
 * for (const auto v : cell->vertex_indices())
 * {
 * If this cell is at the inner boundary, then at least one of its
 * vertices must sit on the inner ring and therefore have a radial
 * distance from the center of exactly 0.5, up to floating point
 * accuracy. So we compute this distance, and if we find a vertex
 * with this property, we flag this cell for later refinement. We
 * can then also break the loop over all vertices and move on to
 * the next cell.
 * 

 * 
 * Because the distance from the center is computed as a floating
 * point number, we have to expect that whatever we compute is
 * only accurate to within
 * [round-off](https://en.wikipedia.org/wiki/Round-off_error). As
 * a consequence, we can never expect to compare the distance
 * with the inner radius by equality: A statement such as
 * `if (distance_from_center == inner_radius)` will fail
 * unless we get exceptionally lucky. Rather, we need to do this
 * comparison with a certain tolerance, and the usual way to do
 * this is to write it as `if (std::abs(distance_from_center -
 * inner_radius) <= tolerance)`
 * where `tolerance` is some small number larger
 * than round-off. The question is how to choose it: We could just
 * pick, say, `1e-10`, but this is only appropriate if the objects
 * we compare are of size one. If we had created a mesh with cells
 * of size `1e+10`, then `1e-10` would be far lower than round-off
 * and, as before, the comparison will only succeed if we get
 * exceptionally lucky. Rather, it is almost always useful to make
 * the tolerance *relative* to a typical "scale" of the objects
 * being compared. Here, the "scale" would be the inner radius, or
 * maybe the diameter of cells. We choose the former and set the
 * tolerance equal to $10^{-6}$ times the inner radius of the
 * annulus.
 * const double distance_from_center =
 * center.distance(cell->vertex(v));
 * 

 * 
 * if (std::fabs(distance_from_center - inner_radius) <=
 * 1e-6 * inner_radius)
 * {
 * cell->set_refine_flag();
 * break;
 * }
 * }
 * }
 * 

 * 
 * Now that we have marked all the cells that we want refined, we let
 * the triangulation actually do this refinement. The function that does
 * so owes its long name to the fact that one can also mark cells for
 * coarsening, and the function does coarsening and refinement all at
 * once:
 * triangulation.execute_coarsening_and_refinement();
 * }
 * 

 * 
 * 

 * 
 * Finally, after these five iterations of refinement, we want to again
 * write the resulting mesh to a file, again in SVG format. This works just
 * as above:
 * std::ofstream out("grid-2.svg");
 * GridOut       grid_out;
 * grid_out.write_svg(triangulation, out);
 * 

 * 
 * std::cout << "Grid written to grid-2.svg" << std::endl;
 * }
 * 

 * 
 * 

 * 
 * 

 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * Finally, the main function. There isn't much to do here, only to call the
 * two subfunctions, which produce the two grids.
 * 
 * @code
 * int main(int argc, const char *argv[])
 * {
 *   first_grid();
 *   return 0;
 * }
 * @endcode
<a name="Results"></a><h1>Results</h1>


Running the program produces graphics of two grids (grid-1.svg and grid-2.svg).
You can open these with most every web browser -- in the simplest case,
just open the current directory in your file system explorer and click
on the file. If you like working on the command line, you call your
web browser with the file: `firefox grid-1.svg`, `google-chrome grid-1.svg`,
or whatever the name of your browser is. If you do this, the two meshes
should look like this:

<table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-1-r9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-2-r9.2.png" alt="">
    </td>
  </tr>
</table>

The left one, well, is not very exciting. The right one is &mdash; at least
&mdash; unconventional. The pictures color-code the "refinement level" of each
cell: How many times did a coarse mesh cell have to be subdivided to obtain
the given cell. In the left image, this is boring since the mesh was
refined globally a number of times, i.e., <i>every</i> cell was
refined the same number of times.

(While the second mesh is entirely artificial and made-up, and
certainly not very practical in applications, to everyone's surprise it
has found its way into the literature: see @cite Mu05. Apparently it is
good for some things at least.)


<a name="Possibleextensions"></a><h3> Possible extensions </h3>


<a name="Differentadaptiverefinementstrategies"></a><h4> Different adaptive refinement strategies </h4>


This program obviously does not have a whole lot of functionality, but
in particular the <code>second_grid</code> function has a bunch of
places where you can play with it. For example, you could modify the
criterion by which we decide which cells to refine. An example would
be to change the condition to this:
@code
      for (auto &cell: triangulation.active_cell_iterators())
        if (cell->center()[1] > 0)
          cell->set_refine_flag ();
@endcode
This would refine all cells for which the $y$-coordinate of the cell's
center is greater than zero (the <code>TriaAccessor::center</code>
function that we call by dereferencing the <code>cell</code> iterator
returns a Point<2> object; subscripting <code>[0]</code> would give
the $x$-coordinate, subscripting <code>[1]</code> the
$y$-coordinate). By looking at the functions that TriaAccessor
provides, you can also use more complicated criteria for refinement.

In general, what you can do with operations of the form
`cell->something()` is a bit difficult to find in the documentation
because `cell` is not a pointer but an iterator. The functions you can
call on a cell can be found in the documentation of the classes
`TriaAccessor` (which has functions that can also be called on faces
of cells or, more generally, all sorts of geometric objects that
appear in a triangulation), and `CellAccessor` (which adds a few
functions that are specific to *cells*).

A more thorough description of the whole iterator concept can be found
in the @ref Iterators documentation module.


<a name="Differentgeometries"></a><h4> Different geometries </h4>


Another possibility would be to generate meshes of entirely different
geometries altogether. While for complex geometries there is no way around
using meshes obtained from mesh generators, there is a good number of
geometries for which deal.II can create meshes using the functions in the
GridGenerator namespace. Many of these geometries (such as the one used in this
example program) contain cells with curved faces: put another way, we expect the
new vertices placed on the boundary to lie along a circle. deal.II handles complex
geometries with the Manifold class (and classes inheriting from it); in particular,
the functions in GridGenerator corresponding to non-Cartesian grids (such as
GridGenerator::hyper_shell or GridGenerator::truncated_cone) attach a Manifold
object to the part of the triangulation that should be curved (SphericalManifold
and CylindricalManifold, respectively) and use another manifold on the parts that
should be flat (FlatManifold). See the documentation
of Manifold or the @ref manifold "manifold module" for descriptions of the design
philosophy and interfaces of these classes. Take a look at what they provide and
see how they could be used in a program like this.

We also discuss a variety of other ways to create and manipulate meshes (and
describe the process of attaching Manifolds) in step-49.


<a name="Commentsaboutprogramminganddebugging"></a><h4> Comments about programming and debugging </h4>


We close with a comment about modifying or writing programs with deal.II in
general. When you start working with tutorial programs or your own
applications, you will find that mistakes happen: your program will contain
code that either aborts the program right away or bugs that simply lead to
wrong results. In either case, you will find it extremely helpful to know how
to work with a debugger: you may get by for a while by just putting debug
output into your program, compiling it, and running it, but ultimately finding
bugs with a debugger is much faster, much more convenient, and more reliable
because you don't have to recompile the program all the time and because you
can inspect the values of variables and how they change.

Rather than postponing learning how to use a debugger till you really can't
see any other way to find a bug, here's the one piece of
advice we will provide in this program: learn how to use a debugger as soon as
possible. It will be time well invested.
(@dealiiVideoLectureSeeAlso{25}) The deal.II Frequently Asked
Questions (FAQ) page linked to from the top-level <a
href="http://www.dealii.org/">deal.II webpage</a> also provides a good number
of hints on debugging deal.II programs.


<a name="Moreaboutgraphicaloutput"></a><h4> More about graphical output </h4>


It is often useful to include meshes into your theses or publications.
For this, it may not be very useful to color-code the cells by
refinement level, and to print the cell number onto each cell. But
it doesn't have to be that way -- the GridOut class allows setting flags
for each possible output format (see the classes in the GridOutFlags
namespace) that control how exactly a mesh is plotted. You can of
course also choose other output file formats such as VTK or VTU; this
is particularly useful for 3d meshes where a 2d format such as SVG
is not particular useful because it fixes a particular viewpoint onto
the 3d object. As a consequence, you might want to explore other
options in the GridOut class.
 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-0.cc"
*/
