#include "time_interface.hpp"
#include "cut_mesh.hpp"



BarycentricActiveMesh::BarycentricActiveMesh(const BarycentricMesh2 &th)
    : ActiveMesh<Mesh2>(static_cast<const Mesh2 &>(th)) {}
    
void BarycentricActiveMesh::truncate(const Interface<Mesh2> &interface, int sign_domain_remove) {
    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);
    not_in_active_mesh_.resize(21); // handles up to 21 domains

    // Resize time-level for stationary case
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(1); // 1 quadrature point in time

    // Prepare containers
    // for (int d = 0; d < dom_size; ++d) {
    //     idx_in_background_mesh_[d].resize(0);
    //     int nt_max = idx_from_background_mesh_[d].size();
    //     idx_in_background_mesh_[d].reserve(nt_max);
    // }
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].clear();
        this->idx_from_background_mesh_[d].clear();
    }


    std::vector<int> nt(dom_size, 0);
    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    for (std::size_t macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (std::size_t i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            if (signK.sign() != sign_domain_remove || interface.isCut(sub_k)) {
                keep_macro = true;
                break;
            }
        }

        if (!keep_macro)
            continue;

        // Keep all three subelements
        for (std::size_t i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            for (int d = 0; d < dom_size; ++d) {
                int local_id = nt[d];

                idx_in_background_mesh_[d].push_back(sub_k);
                idx_from_background_mesh_[d][sub_k] = local_id;

                if (interface.isCut(sub_k)) {
                    interface_id_[0][{d, local_id}].emplace_back(&interface, -sign_domain_remove);
                } else if (signK.sign() == sign_domain_remove) {
                    // std::cout << "sub_k = " << sub_k << "\n";
                    not_in_active_mesh_[d][0][local_id] = true; //! to exclude elements completely outside of domain from isInactive
                }

                ++nt[d];
            }
        }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
}


void BarycentricActiveMesh::truncate(const TimeInterface<Mesh2> &interface, int sign_domain_remove) {
    
    int n_tid = interface.size();
    assert(n_tid < interface_id_.size());
    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);
    not_in_active_mesh_.resize(21); // handles up to 21 domains

    // Resize time-level for stationary case
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(nb_quadrature_time_);

    // Prepare containers
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].clear();
        this->idx_from_background_mesh_[d].clear();
    }

    std::vector<int> nt(dom_size, 0);
    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    for (std::size_t macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        // Loop through subelements
        for (std::size_t i = 0; i < 3; ++i) {
            int sgn;    // hold sign of the current element

            // Loop over all time quadrature points
            for (size_t it = 0; it < interface.size() - 1; ++it) {
                int sub_k = sub_elems[i];
                const auto signKi  = interface(it)->get_SignElement(sub_k);
                const auto signKii = interface(it+1)->get_SignElement(sub_k);

                sgn = signKi.sign();
                const bool Ki_cut = interface(it)->isCut(sub_k);
                const bool Kii_cut = interface(it+1)->isCut(sub_k);
                
                if (signKi.sign() != sign_domain_remove || Ki_cut || Kii_cut || signKi.sign()*signKii.sign() <= 0) {
                    keep_macro = true;
                    break;
                }
            }
        }

        if (!keep_macro)
            continue;

        // Keep all three subelements
        for (std::size_t i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];

            for (int d = 0; d < dom_size; ++d) {
                int local_id = nt[d];

                idx_in_background_mesh_[d].push_back(sub_k);
                idx_from_background_mesh_[d][sub_k] = local_id;

                for (int it = 0; it < n_tid; ++it) {
                    const auto signK = interface(it)->get_SignElement(sub_k);
                    const bool is_cut = interface(it)->isCut(sub_k);

                    if (is_cut) {
                        interface_id_[it][{d, local_id}].emplace_back(interface[it], -sign_domain_remove);
                    } else if (signK.sign() == sign_domain_remove) {
                        not_in_active_mesh_[d][it][local_id] = true;
                    }
                }

                ++nt[d];
            }
        }

    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
}


void BarycentricActiveMesh::createSurfaceMesh(const Interface<Mesh2> &interface) {
    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);
    not_in_active_mesh_.resize(21); // handles up to 21 domains

    // Resize time-level for stationary case
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(1); // 1 quadrature point in time

    // Prepare containers
    // for (int d = 0; d < dom_size; ++d) {
    //     idx_in_background_mesh_[d].resize(0);
    //     int nt_max = idx_from_background_mesh_[d].size();
    //     idx_in_background_mesh_[d].reserve(nt_max);
    // }
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].clear();
        this->idx_from_background_mesh_[d].clear();
    }


    std::vector<int> nt(dom_size, 0);
    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    for (std::size_t macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (std::size_t i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            if (interface.isCut(sub_k)) {
                keep_macro = true;
                break;
            }
        }

        if (!keep_macro)
            continue;

        // Keep all three subelements
        for (std::size_t i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            for (int d = 0; d < dom_size; ++d) {
                int local_id = nt[d];

                idx_in_background_mesh_[d].push_back(sub_k);
                idx_from_background_mesh_[d][sub_k] = local_id;

                if (interface.isCut(sub_k)) {
                    interface_id_[0][{d, local_id}].emplace_back(&interface, 0);
                } else {
                    // std::cout << "sub_k = " << sub_k << "\n";
                    not_in_active_mesh_[d][0][local_id] = true; //! to exclude elements completely outside of domain from isInactive
                }

                ++nt[d];
            }
        }
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
}


// void BarycentricActiveMesh::truncate(const Interface<Mesh2> &interface, int sign_domain_remove) {
//     int dom_size = this->get_nb_domain();
//     idx_element_domain.resize(0);
//     not_in_active_mesh_.resize(21);

//     for (int i = 0; i < 21; ++i)
//         not_in_active_mesh_[i].resize(1);   // 1 quadrature point in time when stationary


//     {
//         // Iterate through number of remaining subdomains
//         for (int d = 0; d < dom_size; ++d) {
//             idx_in_background_mesh_[d].resize(0);
//             // Compute number of elements in subdomain d
//             int nt_max = idx_from_background_mesh_[d].size();
//             // Reserve memory for these elements
//             idx_in_background_mesh_[d].reserve(nt_max);
//         }
//     }

//     std::vector<int> nt(dom_size, 0);
//     const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

//     for (int d = 0; d < dom_size; ++d) {
        
//         for (auto it_k = idx_from_background_mesh_[d].begin(); it_k != idx_from_background_mesh_[d].end();) {
            
//             int kb = it_k->first;  // background mesh element index
//             int k  = it_k->second; // active mesh element index

//             // Get interface segment
//             auto it_gamma = interface_id_[0].find(std::make_pair(d, k));
            
//             // Get the sign of the level set function in the element kb
//             const auto signK = interface.get_SignElement(kb);
//             bool is_cut = interface.isCut(kb);

//             // Skip elements fully in the remove domain *unless* someone in the macro is cut
//             if ((signK.sign() == sign_domain_remove) && !is_cut) {

//                 const auto &elements_in_macro = Th_bary.macro_elements[Th_bary.get_macro_element(kb)];
//                 bool is_macro_cut = false;

//                 std::cout << "Element kb = " << kb << " is NOT cut and OUTSIDE domain\n";
//                 // Check subelements
//                 for (std::size_t sub_k : elements_in_macro) {
//                     if (sub_k == kb) continue;  // already checked in "is_cut"
//                         std::cout << "Neighbor sub_k = " << sub_k << "\n";
//                     if (interface.isCut(sub_k)) {
//                         is_macro_cut = true;
//                         not_in_active_mesh_[d][0][kb] = false;   //! added to try to make isInactive method exclude elements that are entirely outside the domain but in the active mesh

//                         std::cout << "Neighbor sub_k = " << sub_k << " IS CUT\n";
//                         break; // no need to check further
//                     }
//                 }

//                 // Exclude macro if none of its elements are cut
//                 if (!is_macro_cut) {
//                     // std::cout << "it_k = " << it_k->first << ", " << it_k->second << "\n";
//                     std::cout << "Macro element belonging to kb = " << kb << " is NOT active\n";
//                     it_k = idx_from_background_mesh_[d].erase(it_k);
                    
//                     continue; // skip to next
//                 }
//             }

//             // Save and erase old interfaces
//             int nb_interface = (it_gamma == interface_id_[0].end()) ? 0 : it_gamma->second.size();
//             std::vector<const Interface<Mesh2> *> old_interface(nb_interface);
//             std::vector<int> ss(nb_interface);
//             for (int i = 0; i < nb_interface; ++i)
//                 old_interface[i] = it_gamma->second[i].first;
//             for (int i = 0; i < nb_interface; ++i)
//                 ss[i] = it_gamma->second[i].second;
//             if (it_gamma != interface_id_[0].end()) {
//                 auto ittt = interface_id_[0].erase(it_gamma);
//             }

//             // Set new indices and put back interfaces
//             idx_in_background_mesh_[d].push_back(kb);
//             it_k->second = nt[d];
//             for (int i = 0; i < nb_interface; ++i) {
//                 interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(old_interface[i], ss[i]));
//             }

//             // Is cut so need to add interface and sign
//             // if (signK.cut()) {    //! Was like this before
//             if (is_cut) {
//                 interface_id_[0][std::make_pair(d, nt[d])].push_back(std::make_pair(&interface, -sign_domain_remove));
//             }
//             nt[d]++;
//             it_k++;
//         }
//     }

//     idx_element_domain.push_back(0);
//     for (int d = 0; d < dom_size; ++d) {
//         idx_in_background_mesh_[d].resize(nt[d]);
//         idx_in_background_mesh_[d].shrink_to_fit();
//         int sum_nt = idx_element_domain[d] + nt[d];
//         idx_element_domain.push_back(sum_nt);
//     }
// }
    
