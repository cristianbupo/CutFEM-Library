#include "boundary.hpp"

//#include "Mesh2dn.hpp"

// template <>
// BoundaryDirichlet<Mesh2>::BoundaryDirichlet(const cutspace_t &Vh, const BarycentricActiveMesh2& active_mesh, const int domain) {
//     const auto &Th_background = Vh.Th;  // background mesh
    
//     for (int km = 0; km < active_mesh.active_macro_elements.size(); ++km) {

//         if (!active_mesh.is_macro_cut(km))
//             continue;

//         for (int ifac = 0; ifac < elt_t::nea; ++ifac) { 
//             int km_n = active_mesh.macro_adjacent(km, ifac);
//             // Only want faces on the boundary
//             if (km_n != -1)
//                 continue;

//             const auto& micro_elements = active_mesh.active_macro_elements[km];

//             for (int k_micro : micro_elements) {

//                 for (int ifac_micro = 0; ifac_micro < elt_t::nea; ++ifac_micro) {
//                     int jfac_micro = ifac_micro;
//                     int kn_micro = active_mesh.ElementAdj(k_micro, jfac_micro);

//                     if (kn_micro != -1)
//                         continue;
                    
//                     const int kb = active_mesh.idxElementInBackMesh(k_micro);
//                     const auto &FK(Vh[k_micro]);  // Get the finite element for the current element in the finite element space
//                     const auto &T(Th_background[kb]);
//                     const int domain = FK.get_domain();

//                     std::vector<Rd> dof_points(FK.tfe->NbPtforInterpolation);
//                     FK.tfe->global_dofs(T, dof_points);

//                     int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;
                
//                     for (int ic = 0; ic < Vh.N; ++ic) {
//                         for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

//                             // here need a fct that give the item from the dof for P3 for example
//                             int id_item       = df_loc;
//                             bool is_on_border = false;

//                             // case 1: dof is on node. E.g. (0, 1, 2) if triangle
//                             if (id_item < T.nv) {
//                                 for (int i = 0; i < elt_t::nva; ++i) {

//                                     int i_e = elt_t::nvedge.at(ifac_micro).at(i);
//                                     if (i_e == id_item) {
//                                         is_on_border = true;
//                                         break;
//                                     }
//                                 }
//                             }
//                             // case 2: dof is on an edge (or face)
//                             else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

//                                 int id_face = (id_item - T.nv) / ndof_edge_per_component;

//                                 if (id_face == ifac_micro) {
//                                     is_on_border = true;
//                                 }
//                             }
//                             if (is_on_border) {

//                                 Rd P                   = dof_points.at(df_loc);
//                                 std::cout << "Boundary point P = " << P << "\n";
//                                 size_t df_glob         = FK.loc2glb(df);
//                                 boundary_dofs[df_glob] = DofData<D>(kb, domain, ic, P);
//                             }
//                         }
//                     }
//                 }                
//             }
//         }
    
//         std::cout << std::endl;
//     }
// }


template <>
BoundaryDirichlet<Mesh2>::BoundaryDirichlet(const cutspace_t &Vh, const BarycentricActiveMesh2& active_mesh, const int domain) {
    const auto &Th_background = Vh.Th;  // background mesh
    
    for (int km = 0; km < active_mesh.active_macro_elements.size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!active_mesh.is_macro_cut(km))
            continue;

        const auto& micro_elements = active_mesh.active_macro_elements[km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < elt_t::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = active_mesh.ElementAdj(k_micro, jfac);

                if (kn_micro != -1) 
                    continue;

                // const int kn_macro = active_mesh.inverse_active_macro_map.at(kn_micro);
                // std::cout << "kn_micro = " << kn_micro << "\n";
                
                // assert(!active_mesh.is_macro_cut(kn_macro) || !active_mesh.is_macro_interior(kn_macro));

                const int kb = active_mesh.idxElementInBackMesh(k_micro);
                const auto &FK(Vh[k_micro]);  // Get the finite element for the current element in the finite element space
                const auto &T(Th_background[kb]);
                const int domain = FK.get_domain();

                std::vector<Rd> dof_points(FK.tfe->NbPtforInterpolation);
                FK.tfe->global_dofs(T, dof_points);

                int ndof_edge_per_component = FK.tfe->ndfonEdge / Vh.N;

                // std::cout << "Macro element " << km << "\n";
            
                for (int ic = 0; ic < Vh.N; ++ic) {
                    for (size_t df_loc = 0, df = FK.dfcbegin(ic); df < FK.dfcend(ic); ++df, ++df_loc) {

                        // here need a fct that give the item from the dof for P3 for example
                        int id_item       = df_loc;
                        bool is_on_border = false;

                        // case 1: dof is on node. E.g. (0, 1, 2) if triangle
                        if (id_item < T.nv) {
                            for (int i = 0; i < elt_t::nva; ++i) {

                                int i_e = elt_t::nvedge.at(ifac).at(i);
                                if (i_e == id_item) {
                                    is_on_border = true;
                                    break;
                                }
                            }
                        }
                        // case 2: dof is on an edge (or face)
                        else if (id_item < T.nv + T.ne * ndof_edge_per_component) { // df if on an edge

                            int id_face = (id_item - T.nv) / ndof_edge_per_component;

                            if (id_face == ifac) {
                                is_on_border = true;
                            }
                        }
                        if (is_on_border) {

                            Rd P                   = dof_points.at(df_loc);
                            // std::cout << "Boundary point P = " << P << "\n";
                            size_t df_glob         = FK.loc2glb(df);
                            boundary_dofs[df_glob] = DofData<D>(kb, domain, ic, P);
                        }
                    }
                }
            }                
        }
    
        // std::cout << std::endl;
    }
}
