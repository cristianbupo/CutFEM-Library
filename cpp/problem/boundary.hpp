#ifndef CUTFEM_PROBLEM_BOUNDARY_HPP
#define CUTFEM_PROBLEM_BOUNDARY_HPP

#include <algorithm>
#include <limits>
#include <map>
#include <span>
#include <utility>
#include <vector>

#include "../common/point.hpp"
#include "../FESpace/FESpace.hpp"

// ------------------------------
// Per-DOF metadata on the boundary
// ------------------------------
template <int D>
class DofData {
    using v_t = typename typeRd<D>::Rd;

public:
    DofData() = default;
    DofData(int k, int dom, int ci, v_t P) : k(k), domain(dom), ci(ci), P(P) {}

    DofData(DofData&&)            = default;
    DofData& operator=(DofData&&) = default;
    // No copy ctor on purpose (map stores/moves fine, and we iterate by const&)

    int k;          // index in BACKGROUND mesh
    int domain{0};  // active-mesh domain id
    int ci;         // component index
    v_t P;          // evaluation point (global)
};

// ------------------------------
// BoundaryDirichlet: builds dof set on (cut) boundaries and applies strong BCs
// ------------------------------
template <typename M>
class BoundaryDirichlet {
    using Rd         = typename M::Rd;
    static const int D = Rd::d;

public:
    using mesh_t     = M;
    using elt_t      = typename mesh_t::Element;
    using space_t    = GFESpace<mesh_t>;
    using cutspace_t = CutFESpace<mesh_t>;
    using fct_t      = FunFEM<mesh_t>;
    using dof_data_t = DofData<D>;

    // Construct from CutFEM space on the *outer fitted* boundary of the active mesh
    BoundaryDirichlet(const cutspace_t& Vh, std::vector<int> lab);

    // Construct from CutFEM space at the outer boundary of the active mesh (domain-wise).
    // Works for nodal FE; uses topological adjacency in the active mesh.
    BoundaryDirichlet(const cutspace_t& Vh, const int domain = 0);

    // Specialization for barycentric cut meshes is provided in boundary.cpp:
    // BoundaryDirichlet(const cutspace_t& Vh, const BarycentricActiveMesh2& active_mesh, const int domain = 0);
    BoundaryDirichlet(const cutspace_t& Vh,
                      const BarycentricActiveMesh2& active_mesh,
                      const int domain = 0);

    // Standard fitted FEM boundary by labels
    BoundaryDirichlet(const space_t& Vh, std::vector<int> lab);

    // --- Apply strong BCs (Dirichlet elimination) ---

    // Row elimination: zero each constrained row I, set A(I,I)=1.
    // Optional dof_start for when constrained unknowns live in a later block.
    void apply_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                             size_t dof_start = 0);

    // Column elimination + RHS shift for inhomogeneous g:
    // for r!=I: b[r] -= A(r,I)*g, and set A(r,I)=0. Keeps A(I,I)=1.
    void finalize_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                                std::span<double> b,
                                size_t dof_start = 0);

    // Write RHS values g into b(I) from a FE function
    void apply(std::span<double> b, const fct_t& f);

    // Write RHS values g into b(I) from a plain function f(P, component)
    template <typename Fct>
    void apply(std::span<double> b, const Fct& f);

    // Constant value
    void apply(std::span<double> b, double val);

    // Public so you can inspect/debug
    std::map<int, dof_data_t> boundary_dofs;
};

// ------------------------------
// Implementations
// ------------------------------

template <typename M>
BoundaryDirichlet<M>::BoundaryDirichlet(const typename BoundaryDirichlet<M>::cutspace_t& Vh,
                                        std::vector<int> lab)
{
    const auto& Th    = Vh.Th;          // background mesh
    const auto& cutTh = Vh.get_mesh();  // active mesh (fitted boundary elements provided)

    // iterate over fitted boundary elements of the ACTIVE mesh
    for (int idx_be = cutTh.first_boundary_element();
         idx_be < cutTh.last_boundary_element();
         idx_be += cutTh.next_boundary_element())
    {
        int idx_bdry_face = -1; // local face index in the background element
        const int kb = cutTh.Th.BoundaryElement(idx_be, idx_bdry_face); // bg elem index

        // which active-space elements correspond to this background element?
        std::vector<int> idxK = cutTh.idxAllElementFromBackMesh(kb, -1);
        if (idxK.empty()) continue;

        const int k = idxK[0];             // we only accept *uncut* bg faces (one active element)
        const auto& FK(Vh[k]);
        const int domain = FK.get_domain();

        // label filter
        const auto& BE = cutTh.be(idx_be);
        if (std::find(lab.begin(), lab.end(), BE.lab) == lab.end()) continue;

        // skip if the active element is actually cut; we want outer, fitted boundary
        if (cutTh.isCut(k, 0)) continue;

        const auto& T = cutTh.Th[kb];

        std::vector<Rd> dof_point(FK.tfe->NbPtforInterpolation);
        FK.tfe->global_dofs(T, dof_point);

        const int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

        for (int ic = 0; ic < Vh.N; ++ic) {
            for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {
                const int id_item = static_cast<int>(df_loc);
                bool is_on_border = false;

                // vertex dof
                if (id_item < T.nv) {
                    for (int i = 0; i < elt_t::nva; ++i) {
                        const int i_e = elt_t::nvedge.at(idx_bdry_face).at(i);
                        if (i_e == id_item) { is_on_border = true; break; }
                    }
                }
                // edge dof
                else if (id_item < T.nv + T.ne * ndof_edge_per_component) {
                    const int id_face = (id_item - T.nv) / ndof_edge_per_component;
                    if (id_face == idx_bdry_face) is_on_border = true;
                }

                if (is_on_border) {
                    const Rd P = dof_point.at(df_loc);
                    const size_t df_glob = FK.loc2glb(df);
                    boundary_dofs.emplace(static_cast<int>(df_glob),
                                          dof_data_t(kb, domain, ic, P));
                }
            }
        }
    }
}

// NOTE: works for nodal FE; uses "no neighbor" as boundary indicator in active mesh
template <typename M>
BoundaryDirichlet<M>::BoundaryDirichlet(const typename BoundaryDirichlet<M>::cutspace_t& Vh,
                                        const int /*domain_arg*/)
{
    const auto& Th    = Vh.Th;         // background
    const auto& cutTh = Vh.get_mesh(); // active

    for (int k = cutTh.first_element(); k < cutTh.last_element(); k += cutTh.next_element()) {
        if (!cutTh.isCut(k, 0)) continue; // only cut elements can touch the outer boundary

        const int kb = cutTh.idxElementInBackMesh(k);
        const auto& FK(Vh[k]);
        const auto& T = Th[kb];
        const int dom = FK.get_domain();

        for (int ifac = 0; ifac < M::Element::nea; ++ifac) {
            int jfac = ifac;
            const int kn = cutTh.ElementAdj(k, jfac);

            // boundary face => no neighbor in the active mesh
            if (kn != -1) continue;

            std::vector<Rd> dof_points(FK.tfe->NbPtforInterpolation);
            FK.tfe->global_dofs(T, dof_points);

            const int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

            for (int ic = 0; ic < Vh.N; ++ic) {
                for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {
                    const int id_item = static_cast<int>(df_loc);
                    bool is_on_border = false;

                    if (id_item < T.nv) {
                        for (int i = 0; i < elt_t::nva; ++i) {
                            const int i_e = elt_t::nvedge.at(ifac).at(i);
                            if (i_e == id_item) { is_on_border = true; break; }
                        }
                    } else if (id_item < T.nv + T.ne * ndof_edge_per_component) {
                        const int id_face = (id_item - T.nv) / ndof_edge_per_component;
                        if (id_face == ifac) is_on_border = true;
                    }

                    if (is_on_border) {
                        const Rd P = dof_points.at(df_loc);
                        const size_t df_glob = FK.loc2glb(df);
                        boundary_dofs.emplace(static_cast<int>(df_glob),
                                              dof_data_t(kb, dom, ic, P));
                    }
                }
            }
        }
    }
}

// Decl only; your Mesh2/BarycentricActiveMesh2 specialization stays in boundary.cpp
// template <typename M>
// BoundaryDirichlet<M>::BoundaryDirichlet(const typename BoundaryDirichlet<M>::cutspace_t&,
//                                         const BarycentricActiveMesh2&,
//                                         const int) {}

// Standard fitted FEM boundary via labels
template <typename M>
BoundaryDirichlet<M>::BoundaryDirichlet(const typename BoundaryDirichlet<M>::space_t& Vh,
                                        std::vector<int> lab)
{
    const mesh_t& Th = Vh.Th;

    for (int k = Th.first_boundary_element(); k < Th.last_boundary_element(); k += Th.next_boundary_element()) {
        const auto& BE  = Th.be(k);
        if (std::find(lab.begin(), lab.end(), BE.lab) == lab.end()) continue;

        auto [elt_idx, face_idx] = Th.getBoundaryElement(k);
        const auto& T  = Th[elt_idx];
        const auto& FK = Vh[elt_idx];
        const int dom  = FK.get_domain();

        std::vector<Rd> dof_point(FK.tfe->NbPtforInterpolation);
        FK.tfe->global_dofs(T, dof_point);

        const int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

        for (int ic = 0; ic < Vh.N; ++ic) {
            for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {
                const int id_item = static_cast<int>(df_loc);
                bool is_on_border = false;

                if (id_item < T.nv) {
                    for (int i = 0; i < elt_t::nva; ++i) {
                        const int i_e = elt_t::nvedge.at(face_idx).at(i);
                        if (i_e == id_item) { is_on_border = true; break; }
                    }
                } else if (id_item < T.nv + T.ne * ndof_edge_per_component) {
                    const int id_face = (id_item - T.nv) / ndof_edge_per_component;
                    if (id_face == face_idx) is_on_border = true;
                }

                if (is_on_border) {
                    const Rd P  = dof_point.at(df_loc);
                    const size_t df_glob = FK.loc2glb(df);
                    boundary_dofs.emplace(static_cast<int>(df_glob),
                                          dof_data_t(elt_idx, dom, ic, P));
                }
            }
        }
    }
}

// -------- Row elimination: zero row I, set A(I,I)=1
template <typename M>
void BoundaryDirichlet<M>::apply_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                                               size_t dof_start)
{
    using Key = std::pair<int,int>;

    for (const auto& kv : boundary_dofs) {
        const int I = static_cast<int>(kv.first + dof_start);

        auto row_begin = A.lower_bound(Key{I, std::numeric_limits<int>::min()});
        auto row_end   = A.lower_bound(Key{I+1, std::numeric_limits<int>::min()});
        A.erase(row_begin, row_end);

        A[Key{I, I}] = 1.0;
    }
}

// -------- Column elimination + RHS shift for inhomogeneous g
template <typename M>
void BoundaryDirichlet<M>::finalize_inhomogeneous(std::map<std::pair<int,int>, double>& A,
                                                  std::span<double> b,
                                                  size_t dof_start)
{
    for (const auto& kv : boundary_dofs) {
        const int I = static_cast<int>(kv.first + dof_start);
        const double g = b[static_cast<size_t>(I)];

        for (auto it = A.begin(); it != A.end(); ) {
            const int r = it->first.first;
            const int c = it->first.second;
            if (c == I && r != I) {
                b[static_cast<size_t>(r)] -= it->second * g; // RHS shift
                it = A.erase(it);                            // zero A(r,I)
            } else {
                ++it;
            }
        }
        // A(I,I)=1 already set; row I already zeroed.
    }

}

// -------- Write RHS values from FE function on background mesh
template <typename M>
void BoundaryDirichlet<M>::apply(std::span<double> b, const fct_t& f)
{
    for (const auto& kv : boundary_dofs) {
        const int df = kv.first;
        const auto& d = kv.second;
        b[static_cast<size_t>(df)] =
            f.evalOnBackMesh(d.k, d.domain, d.P, d.ci, op_id);
    }
}

// -------- Write RHS values from raw function f(P, component)
template <typename M>
template <typename Fct>
void BoundaryDirichlet<M>::apply(std::span<double> b, const Fct& f)
{
    for (const auto& kv : boundary_dofs) {
        const int df = kv.first;
        const auto& d = kv.second;
        b[static_cast<size_t>(df)] = f(d.P, d.ci);
    }
}

template <typename M>
void BoundaryDirichlet<M>::apply(std::span<double> b, double val)
{
    for (const auto& kv : boundary_dofs) {
        b[static_cast<size_t>(kv.first)] = val;
    }
}

#endif // CUTFEM_PROBLEM_BOUNDARY_HPP
