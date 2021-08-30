// B.6. The Poisson Equation on a Distributed Grid
// File: grid-distributed-poisson.cc
// Equation: Poisson equation
// Discretization: First-order Lagrange finite elements
// Grid: Unstructured simplicial grid, implemented with UGGrid
// Solver: CG solver with Jacobi preconditioner
// Distributed: yes
// Discussed in: Chapter 6.7

#include <config.h>

#include <vector>
#include <map>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class Matrix>
void assembleElementStiffnessMatrix(const LocalView &localView,
                                    Matrix &elementMatrix)
{
    using Element = typename LocalView::Element;
    constexpr int dim = Element::dimension;
    auto element = localView.element();
    auto geometry = element.geometry();

    // Get set of shape functions for this element
    const auto &localFiniteElement = localView.tree().finiteElement();

    // Set all matrix entries to zero
    elementMatrix.setSize(localView.size(), localView.size());
    elementMatrix = 0; // Fill the entire matrix with zeros

    // Get a quadrature rule
    int order = 2 * (dim * localFiniteElement.localBasis().order() - 1);
    const auto &quadRule = QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (const auto &quadPoint : quadRule)
    {
        // Position of the current quadrature point in the reference element
        const auto quadPos = quadPoint.position();

        // The transposed inverse Jacobian of the map from the reference element
        // to the grid element
        const auto jacobian = geometry.jacobianInverseTransposed(quadPos);

        // The multiplicative factor in the integral transformation formula
        const auto integrationElement = geometry.integrationElement(quadPos);

        // The gradients of the shape functions on the reference element
        std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
        localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                         referenceGradients);

        // Compute the shape function gradients on the grid element
        std::vector<FieldVector<double, dim>> gradients(referenceGradients.size());
        for (size_t i = 0; i < gradients.size(); i++)
            jacobian.mv(referenceGradients[i][0], gradients[i]);

        // Compute the actual matrix entries
        for (size_t p = 0; p < elementMatrix.N(); p++)
        {
            auto localRow = localView.tree().localIndex(p);
            for (size_t q = 0; q < elementMatrix.M(); q++)
            {
                auto localCol = localView.tree().localIndex(q);
                elementMatrix[localRow][localCol] += (gradients[p] * gradients[q]) * quadPoint.weight() * integrationElement;
            }
        }
    }
}

// Compute the source term for a single element
template <class LocalView>
void assembleElementVolumeTerm(
    const LocalView &localView,
    BlockVector<double> &localB,
    const std::function<double(FieldVector<double, LocalView::Element::dimension>)> volumeTerm)
{
    using Element = typename LocalView::Element;
    auto element = localView.element();
    constexpr int dim = Element::dimension;

    // Set of shape functions for a single element
    const auto &localFiniteElement = localView.tree().finiteElement();

    // Set all entries to zero
    localB.resize(localFiniteElement.size());
    localB = 0;

    // A quadrature rule
    int order = dim;
    const auto &quadRule = QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (const auto &quadPoint : quadRule)
    {
        // Position of the current quadrature point in the reference element
        const FieldVector<double, dim> &quadPos = quadPoint.position();

        // The multiplicative factor in the integral transformation formula
        const double integrationElement = element.geometry().integrationElement(quadPos);

        double functionValue = volumeTerm(element.geometry().global(quadPos));

        // Evaluate all shape function values at this point
        std::vector<FieldVector<double, 1>> shapeFunctionValues;
        localFiniteElement.localBasis().evaluateFunction(quadPos, shapeFunctionValues);

        // Actually compute the vector entries
        for (size_t p = 0; p < localB.size(); p++)
        {
            auto localIndex = localView.tree().localIndex(p);
            localB[localIndex] += shapeFunctionValues[p] * functionValue * quadPoint.weight() * integrationElement;
        }
    }
}

// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis &basis, MatrixIndexSet &nb)
{
    nb.resize(basis.size(), basis.size());

    auto gridView = basis.gridView();

    // A loop over all elements of the grid
    auto localView = basis.localView();

    for (const auto &element : elements(gridView))
    {
        localView.bind(element);

        for (size_t i = 0; i < localView.size(); i++)
        {
            // The global index of the i-th vertex of the element
            auto row = localView.index(i);

            for (size_t j = 0; j < localView.size(); j++)
            {
                // The global index of the j-th vertex of the element
                auto col = localView.index(j);
                nb.add(row, col);
            }
        }
    }
}

// Assemble the Laplace stiffness matrix on the given grid view
template <class Basis>
void assemblePoissonProblem(const Basis &basis,
                            BCRSMatrix<double> &matrix,
                            BlockVector<double> &b,
                            const std::function<
                                double(FieldVector<double,
                                                   Basis::GridView::dimension>)>
                                volumeTerm)
{
    auto gridView = basis.gridView();

    // MatrixIndexSets store the occupation pattern of a sparse matrix.
    // They are not particularly efficient, but simple to use.
    MatrixIndexSet occupationPattern;
    getOccupationPattern(basis, occupationPattern);

    // ... and give it the occupation pattern we want.
    occupationPattern.exportIdx(matrix);

    // Set all entries to zero
    matrix = 0;

    // Set b to correct length
    b.resize(basis.dimension());

    // Set all entries to zero
    b = 0;

    // A loop over all elements of the grid
    auto localView = basis.localView();

    for (const auto &element : elements(gridView, Partitions::interior))
    {
        // Now let’s get the element stiffness matrix
        // A dense matrix is used for the element stiffness matrix
        localView.bind(element);

        Matrix<double> elementMatrix;
        assembleElementStiffnessMatrix(localView, elementMatrix);

        for (size_t p = 0; p < elementMatrix.N(); p++)
        {
            // The global index of the p-th degree of freedom of the element
            auto row = localView.index(p);

            for (size_t q = 0; q < elementMatrix.M(); q++)
            {
                // The global index of the q-th degree of freedom of the element
                auto col = localView.index(q);
                matrix[row][col] += elementMatrix[p][q];
            }
        }

        // Now get the local contribution to the right-hand side vector
        BlockVector<double> localB;
        assembleElementVolumeTerm(localView, localB, volumeTerm);

        for (size_t p = 0; p < localB.size(); p++)
        {
            // The global index of the p-th vertex of the element
            auto row = localView.index(p);
            b[row] += localB[p];
        }
    }
}

// { lb_data_handle_begin }
template <class Grid, class AssociativeContainer>
struct LBVertexDataHandle
    : public CommDataHandleIF<LBVertexDataHandle<Grid, AssociativeContainer>,
                              typename AssociativeContainer::mapped_type>
{
    LBVertexDataHandle(const std::shared_ptr<Grid> &grid,
                       AssociativeContainer &dataContainer)
        : idSet_(grid->localIdSet()), dataContainer_(dataContainer)
    {
    }

    bool contains(int dim, int codim) const
    {
        assert(dim == Grid::dimension);
        return (codim == dim); // Only vertices have data
    }

    bool fixedSize(int dim, int codim) const
    {
        return true; // All vertices carry the same number of data items
    }

    template <class Entity>
    size_t size(const Entity &entity) const
    {
        return 1; // One data item per vertex
    }

    template <class MessageBuffer, class Entity>
    void gather(MessageBuffer &buffer, const Entity &entity) const
    {
        auto id = idSet_.id(entity);
        buffer.write(dataContainer_[id]);
    }

    template <class MessageBuffer, class Entity>
    void scatter(MessageBuffer &buffer, const Entity &entity, size_t n)
    {
        assert(n == 1); // This data handle implementation
        // transfers only one data item.
        auto id = idSet_.id(entity);
        buffer.read(dataContainer_[id]);
    }

private:
    const typename Grid::LocalIdSet &idSet_;
    AssociativeContainer &dataContainer_;
};
// { lb_data_handle_end }

// A DataHandle class to communicate and add vertex data
// { comm_data_handle_begin }
template <class GridView, class Vector>
struct VertexDataUpdate
    : public Dune::CommDataHandleIF<VertexDataUpdate<GridView, Vector>,
                                    typename Vector::value_type>
{
    using DataType = typename Vector::value_type;

    // Constructor
    VertexDataUpdate(const GridView &gridView,
                     const Vector &userDataSend,
                     Vector &userDataReceive)
        : gridView_(gridView),
          userDataSend_(userDataSend),
          userDataReceive_(userDataReceive)
    {
    }

    // True if data for this codim should be communicated
    bool contains(int dim, int codim) const
    {
        return (codim == dim); // Only vertices have data
    }

    // True if data size per entity of given codim is constant
    bool fixedSize(int dim, int codim) const
    {
        return true; // All vertices carry the same number of data items
    }

    // How many objects of type DataType have to be sent for a given entity
    template <class Entity>
    size_t size(const Entity &e) const
    {
        return 1; // One data item per vertex
    }

    // Pack user data into message buffer
    template <class MessageBuffer, class Entity>
    void gather(MessageBuffer &buffer, const Entity &entity) const
    {
        auto index = gridView_.indexSet().index(entity);
        buffer.write(userDataSend_[index]);
    }

    // Unpack user data from message buffer
    template <class MessageBuffer, class Entity>
    void scatter(MessageBuffer &buffer, const Entity &entity, size_t n)
    {
        assert(n == 1);
        DataType x;
        buffer.read(x);

        userDataReceive_[gridView_.indexSet().index(entity)] += x;
    }

private:
    const GridView gridView_;
    const Vector &userDataSend_;
    Vector &userDataReceive_;
};
// { comm_data_handle_end }

// { main_begin }
int main(int argc, char *argv[])
{
    // Set up MPI
    const MPIHelper &mpiHelper = MPIHelper::instance(argc, argv);
    // { mpihelper_end }

    //////////////////////////////////
    // Generate the grid
    //////////////////////////////////

    // { create_grid_begin }
    constexpr int dim = 2;
    using Grid = UGGrid<dim>;
    using GridView = Grid::LeafGridView;

    std::shared_ptr<Grid> grid = GmshReader<Grid>::read("l-shape-refined.msh");
    auto gridView = grid->leafGridView();
    // { create_grid_end }

    // { sample_initial_iterate_begin }
    std::vector<double> dataVector;

    if (mpiHelper.rank() == 0)
    {
        // The initial iterate as a function
        auto initialIterate = [](auto x)
        { return std::min(x[0], x[1]); };

        // Sample on the grid vertices
        dataVector.resize(gridView.size(dim));
        for (const auto &vertex : vertices(gridView,
                                           Dune::Partitions::interiorBorder))
        {
            auto index = gridView.indexSet().index(vertex);
            dataVector[index] = initialIterate(vertex.geometry().corner(0));
        }
    }
    // { sample_initial_iterate_end }

    // { data_into_map_begin }
    // Copy vertex data into associative container
    using PersistentContainer = std::map<Grid::LocalIdSet::IdType, double>;
    PersistentContainer persistentContainer;
    const auto &idSet = grid->localIdSet();

    for (const auto &vertex : vertices(gridView))
        persistentContainer[idSet.id(vertex)] = dataVector[gridView.indexSet().index(vertex)];
    // { data_into_map_end }

    // { load_balancing_begin }
    // Distribute the grid and the data
    LBVertexDataHandle<Grid, PersistentContainer>
        dataHandle(grid, persistentContainer);
    grid->loadBalance(dataHandle);
    // { load_balancing_end }

    // { data_from_map_begin }
    // Get gridView again after load-balancing, to make sure it is up-to-date
    gridView = grid->leafGridView();

    // Copy data back into the array
    dataVector.resize(gridView.size(dim));

    for (const auto &vertex : vertices(gridView))
        dataVector[gridView.indexSet().index(vertex)] = persistentContainer[idSet.id(vertex)];
    // { data_from_map_end }

    /////////////////////////////////////////////////////////
    // Stiffness matrix and right hand side vector
    /////////////////////////////////////////////////////////

    // { create_matrix_vector_begin }
    using Matrix = BCRSMatrix<double>;
    using Vector = BlockVector<double>;

    Matrix stiffnessMatrix;
    Vector b;

    auto sourceTerm = [](const FieldVector<double, dim> &x)
    { return -5.0; };

    // Assemble the Poisson system in a first-order Lagrange space
    Functions::LagrangeBasis<GridView, 1> basis(gridView);
    assemblePoissonProblem(basis, stiffnessMatrix, b, sourceTerm);
    // { call_assembler_end }

    // Obtain a consistent representation of the matrix diagonal
    // { make_consistent_diagonal_begin }
    Vector diagonal(basis.size());
    for (std::size_t i = 0; i < basis.size(); ++i)
        diagonal[i] = stiffnessMatrix[i][i];

    auto consistentDiagonal = diagonal;
    VertexDataUpdate<GridView, Vector> matrixDataHandle(gridView,
                                                        diagonal,
                                                        consistentDiagonal);

    gridView.communicate(matrixDataHandle,
                         InteriorBorder_InteriorBorder_Interface,
                         ForwardCommunication);
    // { make_consistent_diagonal_end }

    // Determine Dirichlet degrees of freedom by marking all degrees
    // of freedom whose Lagrange nodes comply with a given predicate.
    // { dirichlet_marking_begin }
    auto dirichletPredicate = [](auto p)
    {
        return p[0] < 1e-8 || p[1] < 1e-8 || (p[0] > 0.4999 && p[1] > 0.4999);
    };

    // Interpolating the predicate will mark
    // all desired Dirichlet degrees of freedom
    std::vector<bool> dirichletNodes;
    Functions::interpolate(basis, dirichletNodes, dirichletPredicate);
    // { dirichlet_marking_end }

    //////////////////////////////////////////////////////////////
    // Modify Dirichlet matrix rows and load vector entries
    //////////////////////////////////////////////////////////////

    // { dirichlet_modification_begin }
    // Loop over the matrix rows
    for (size_t i = 0; i < stiffnessMatrix.N(); i++)
    {
        if (dirichletNodes[i])
        {
            auto cIt = stiffnessMatrix[i].begin();
            auto cEndIt = stiffnessMatrix[i].end();
            // Loop over nonzero matrix entries in current row
            for (; cIt != cEndIt; ++cIt)
                *cIt = (i == cIt.index()) ? 1.0 : 0.0;

            // Modify corresponding load vector entry
            b[i] = dataVector[i];
        }
    }
    // { dirichlet_modification_end }

    ///////////////////////////
    // Compute solution
    ///////////////////////////

    // { algebraic_solving_preprocess_begin }
    // Set the initial iterate
    Vector x(basis.size());
    std::copy(dataVector.begin(), dataVector.end(), x.begin());

    // Solver parameters
    double reduction = 1e-3; // Desired residual reduction factor
    int maxIterations = 50;  // Maximum number of iterations
    // { algebraic_solving_preprocess_end }

    // Solve!
    // { algebraic_solving_begin }
    auto r = b;
    stiffnessMatrix.mmv(x, r); // r -= Ax
    // { initial_residual_end }

    // { make_residual_consistent_begin }
    // Construct consistent representation of the data in r
    auto rConsistent = r;

    VertexDataUpdate<GridView, Vector> vertexUpdateHandle(gridView,
                                                          r,
                                                          rConsistent);

    gridView.communicate(vertexUpdateHandle,
                         InteriorBorder_InteriorBorder_Interface,
                         ForwardCommunication);
    // { make_residual_consistent_end }

    // { global_initial_residual_begin }
    double defect0 = r.dot(rConsistent); // Norm on the local process
    defect0 = grid->comm().sum(defect0);
    defect0 = sqrt(defect0);
    // { global_initial_residual_end }

    // { output_header_begin }
    if (mpiHelper.rank() == 0)
    {
        std::cout << " Iteration Defect Rate" << std::endl;
        ;
        std::cout << " 0" << std::setw(16) << defect0 << std::endl;
    }
    // { output_header_end }

    // { initial_direction_begin }
    // Construct initial search direction in variable d by applying
    // the Jacobi preconditioner to the residual in r.
    Vector d(r.size());
    for (std::size_t i = 0; i < stiffnessMatrix.N(); ++i)
    {
        d[i] = 0;
        if (std::abs(consistentDiagonal[i]) > 1e-5) // Degree of freedom
            // is not on ghost vertex
            d[i] = rConsistent[i] / consistentDiagonal[i];
    }
    // { initial_direction_end }

    // { orthogonalization_begin }
    double rho = d.dot(r);
    rho = grid->comm().sum(rho);
    // { orthogonalization_end }

    // { loop_and_alpha_begin }
    // Current residual norm
    double defect = defect0;

    for (int k = 0; k < maxIterations; ++k)
    {
        // Step length in search direction d
        Vector tmp(d.size());
        stiffnessMatrix.mv(d, tmp);     // tmp=Ad^k
        double alphaDenom = d.dot(tmp); // Scalar product
        alphaDenom = grid->comm().sum(alphaDenom);
        double alpha = rho / alphaDenom;
        // { loop_and_alpha_end }

        // { update_iterate_begin }
        x.axpy(alpha, d); // Update iterate
        // { update_iterate_end }
        // { update_residual_begin }
        r.axpy(-alpha, tmp); // Update residual
        // { update_residual_end }

        // Convergence test
        // { check_convergence_begin }
        // Compute residual norm again
        rConsistent = r;
        gridView.communicate(vertexUpdateHandle,
                             InteriorBorder_InteriorBorder_Interface,
                             ForwardCommunication);

        auto residualNorm = r.dot(rConsistent);
        residualNorm = grid->comm().sum(residualNorm);
        residualNorm = sqrt(residualNorm);

        if (mpiHelper.rank() == 0)
        {
            std::cout << std::setw(5) << k + 1 << " ";
            std::cout << std::setw(16) << residualNorm << " ";
            // Estimated convergence rate
            std::cout << std::setw(16) << residualNorm / defect << std::endl;
        }

        defect = residualNorm; // Update norm

        if (defect < defect0 * reduction) // Convergence check
            break;
        // { check_convergence_end }

        // Determine new search direction
        // { compute_new_direction_begin }
        // Precondition the residual
        Vector preconditionedResidual(d.size());
        for (std::size_t i = 0; i < stiffnessMatrix.N(); i++)
        {
            preconditionedResidual[i] = 0;
            if (std::abs(consistentDiagonal[i]) > 1e-5) // Degree of freedom
                // not on ghost vertex
                preconditionedResidual[i] = rConsistent[i] / consistentDiagonal[i];
        }

        double rhoNext = preconditionedResidual.dot(r);
        rhoNext = grid->comm().sum(rhoNext);
        double beta = rhoNext / rho;

        // Compute new search direction
        d *= beta;
        d += preconditionedResidual;
        rho = rhoNext; // Remember rho for the next iterate
    }
    // { algebraic_solving_end }

    // Output result
    // { vtk_output_begin }
    // For visualization: Write the rank number for each element
    MultipleCodimMultipleGeomTypeMapper<GridView>
        elementMapper(gridView, mcmgElementLayout());
    std::vector<int> ranks(elementMapper.size());
    for (const auto &element : elements(gridView))
        ranks[elementMapper.index(element)] = mpiHelper.rank();

    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(x, "solution");
    vtkWriter.addCellData(ranks, "ranks");
    vtkWriter.write("grid-distributed-poisson-result");
    // { vtk_output_end }
}