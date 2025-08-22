#include "time_interface.hpp"
#include "cut_mesh.hpp"



BarycentricActiveMesh2::BarycentricActiveMesh2(const BarycentricMesh2 &th)
    : ActiveMesh<Mesh2>(static_cast<const Mesh2 &>(th)) {}
    
// --- Stationary --------------------------------------------------------------
void BarycentricActiveMesh2::truncate(const Interface<Mesh2>& interface,
                                      int sign_domain_remove)
{
    // one time slice
    nb_quadrature_time_ = 1;

    const int dom_size = this->get_nb_domain();
    assert(dom_size == 1);
    const auto& Th_bary = static_cast<const BarycentricMesh2&>(this->Th);

    // reset containers
    idx_element_domain.clear();
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].clear();
        this->idx_from_background_mesh_[d].clear();
    }

    // interface + inactive flags
    interface_id_.assign(nb_quadrature_time_, {}); // [0] cleared
    not_in_active_mesh_.assign(21, {});
    for (int i = 0; i < 21; ++i) not_in_active_mesh_[i].assign(nb_quadrature_time_, {});

    // macro bookkeeping
    std::vector<int> nt(dom_size, 0);        // number of active elements in each domain
    int nt_macros = 0; // number of active macro elements 

    // --- decide kept macros ---
    for (int macro_k = 0; macro_k < (int)Th_bary.macro_elements.size(); ++macro_k) {
        const auto& sub = Th_bary.macro_elements[macro_k]; // 3 background subelements

        bool keep_macro = false;
        for (int i = 0; i < 3 && !keep_macro; ++i) {
            const int sub_k = sub[i];
            const auto signK = interface.get_SignElement(sub_k);
            if (interface.isCut(sub_k) || signK.sign() != sign_domain_remove)
                keep_macro = true;
        }
        if (!keep_macro) continue;

        macro_in_background_mesh.push_back(macro_k);
        macro_in_active_mesh[macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // keep all three subelements (into every kept domain; adapt if you keep a subset)
        for (int i = 0; i < 3; ++i) {
            const int sub_k = sub[i];
            const bool is_cut = interface.isCut(sub_k);
            const int sign   = interface.get_SignElement(sub_k).sign();

            for (int d = 0; d < dom_size; ++d) {
                const int local_id = nt[d];     // index in active mesh
                this->idx_in_background_mesh_[d].push_back(sub_k);
                this->idx_from_background_mesh_[d][sub_k] = local_id;
                sub_elems_active[i] = local_id;
                inverse_active_macro_map.push_back(nt_macros);

                if (is_cut) {
                    interface_id_[0][{d, local_id}].emplace_back(&interface, -sign_domain_remove);
                } else if (sign == sign_domain_remove) {
                    // mark "outside" for stationary t=0
                    not_in_active_mesh_[d][0][local_id] = true;
                }
                ++nt[d];
            }
        }
        active_macro_elements.push_back(sub_elems_active);
        nt_macros++;
    }

    nb_active_macros = nt_macros;
    // finalize prefix sums
    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }

}




// --- Time-dependent ----------------------------------------------------------
void BarycentricActiveMesh2::truncate(const TimeInterface<Mesh2>& interface,
                                      int sign_domain_remove)
{
    const int n_tid = interface.size();
    nb_quadrature_time_ = n_tid;

    const int dom_size = this->get_nb_domain();
    assert(dom_size == 1);
    const auto& Th_bary = static_cast<const BarycentricMesh2&>(this->Th);

    // reset containers
    idx_element_domain.clear();
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].clear();
        this->idx_from_background_mesh_[d].clear();
    }

    interface_id_.assign(nb_quadrature_time_, {}); // time slices
    not_in_active_mesh_.assign(21, {});
    for (int i = 0; i < 21; ++i) not_in_active_mesh_[i].assign(nb_quadrature_time_, {});

    std::vector<int> nt(dom_size, 0);
    int nt_macros = 0; // number of active macro elements 

    // --- decide kept macros (exists t where inside OR cut OR sign flips) ---
    for (int macro_k = 0; macro_k < (int)Th_bary.macro_elements.size(); ++macro_k) {
        const auto& sub = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (int i = 0; i < 3 && !keep_macro; ++i) {
            const int sub_k = sub[i];

            for (int it = 0; it + 1 < n_tid; ++it) {
                const int s_i  = interface(it)->get_SignElement(sub_k).sign();
                const int s_ip = interface(it+1)->get_SignElement(sub_k).sign();
                const bool cut_i  = interface(it)->isCut(sub_k);
                const bool cut_ip = interface(it+1)->isCut(sub_k);

                if (cut_i || cut_ip || s_i != sign_domain_remove || s_i * s_ip <= 0) {
                    keep_macro = true;
                    break;
                }
            }
        }
        if (!keep_macro) continue;

        macro_in_background_mesh.push_back(macro_k);
        macro_in_active_mesh[macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // keep all three subelements
        for (int i = 0; i < 3; ++i) {
            const int sub_k = sub[i];

            for (int d = 0; d < dom_size; ++d) {
                const int local_id = nt[d];
                sub_elems_active[i] = local_id;
                this->idx_in_background_mesh_[d].push_back(sub_k);
                this->idx_from_background_mesh_[d][sub_k] = local_id;
                inverse_active_macro_map.push_back(nt_macros);

                // per-time bookkeeping
                for (int it = 0; it < n_tid; ++it) {
                    const bool is_cut = interface(it)->isCut(sub_k);
                    const int sign    = interface(it)->get_SignElement(sub_k).sign();

                    if (is_cut) {
                        interface_id_[it][{d, local_id}].emplace_back(interface[it], -sign_domain_remove);
                    } else if (sign == sign_domain_remove) {
                        not_in_active_mesh_[d][it][local_id] = true; // outside at this t
                    }
                }

                ++nt[d];
            }
        }
        active_macro_elements.push_back(sub_elems_active);
        nt_macros++;
    }

    // finalize prefix sums
    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
    

}





void BarycentricActiveMesh2::createSurfaceMesh(const Interface<Mesh2> &interface) {
    int dom_size = this->get_nb_domain();
    assert(dom_size == 1);
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
    int nt_macros = 0; // number of active macro elements 

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    for (int macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool keep_macro = false;
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            if (interface.isCut(sub_k)) {
                keep_macro = true;
                break;
            }
        }

        if (!keep_macro)
            continue;

        macro_in_background_mesh.push_back(macro_k);
        macro_in_active_mesh[macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // Keep all three subelements
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];
            const auto signK = interface.get_SignElement(sub_k);

            for (int d = 0; d < dom_size; ++d) {
                int local_id = nt[d];
                sub_elems_active[i] = local_id;
                idx_in_background_mesh_[d].push_back(sub_k);
                idx_from_background_mesh_[d][sub_k] = local_id;
                inverse_active_macro_map.push_back(nt_macros);


                if (interface.isCut(sub_k)) {
                    interface_id_[0][{d, local_id}].emplace_back(&interface, 0);
                } else {
                    // std::cout << "sub_k = " << sub_k << "\n";
                    not_in_active_mesh_[d][0][local_id] = true; //! to exclude elements completely outside of domain from isInactive
                }

                ++nt[d];
            }
        }
        active_macro_elements.push_back(sub_elems_active);
        nt_macros++;
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
}

void BarycentricActiveMesh2::createSurfaceMesh(const TimeInterface<Mesh2> &interface) {
    int n_tid = interface.size();
    nb_quadrature_time_ = n_tid;

    not_in_active_mesh_.resize(21);
    for (int i = 0; i < 21; ++i)
        not_in_active_mesh_[i].resize(n_tid);

    int dom_size = this->get_nb_domain();
    idx_element_domain.resize(0);

    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].clear();
        idx_from_background_mesh_[d].clear();
    }

    std::vector<int> nt(dom_size, 0);
    int nt_macros = 0; // number of active macro elements 

    const auto &Th_bary = static_cast<const BarycentricMesh2 &>(this->Th);

    for (int macro_k = 0; macro_k < Th_bary.macro_elements.size(); ++macro_k) {
        const auto &sub_elems = Th_bary.macro_elements[macro_k];

        bool macro_is_active = false;

        // Check all subelements and time quadrature points  
        for (int i = 0; i < 3 && !macro_is_active; ++i) {
            int sub_k = sub_elems[i];
            for (int it = 0; it < n_tid - 1; ++it) {
                const auto signKi  = interface(it)->get_SignElement(sub_k);
                const auto signKii = interface(it + 1)->get_SignElement(sub_k);
                bool Ki_cut  = interface(it)->isCut(sub_k);
                bool Kii_cut = interface(it + 1)->isCut(sub_k);

                if (Ki_cut || Kii_cut || signKi.sign() * signKii.sign() <= 0) {
                    macro_is_active = true;
                    break;
                }
            }
        }

        if (!macro_is_active)
            continue;
        
        macro_in_background_mesh.push_back(macro_k);
        macro_in_active_mesh[macro_k] = nt_macros;
        
        std::array<int, 3> sub_elems_active;

        // Keep all subelements if macro is active
        for (int i = 0; i < 3; ++i) {
            int sub_k = sub_elems[i];

            for (int d = 0; d < dom_size; ++d) {
                int local_id = nt[d];
                sub_elems_active[i] = local_id;

                idx_in_background_mesh_[d].push_back(sub_k);
                idx_from_background_mesh_[d][sub_k] = local_id;
                inverse_active_macro_map.push_back(nt_macros);


                for (int t = 0; t < n_tid; ++t) {
                    bool is_cut = interface(t)->isCut(sub_k);

                    if (is_cut) {
                        interface_id_[t][{d, local_id}].emplace_back(interface[t], 0);
                    } else {
                        not_in_active_mesh_[d][t][local_id] = true;
                    }
                }

                ++nt[d];
            }
        }
        active_macro_elements.push_back(sub_elems_active);
        nt_macros++;
    }

    idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        idx_in_background_mesh_[d].shrink_to_fit();
        idx_element_domain.push_back(idx_element_domain.back() + nt[d]);
    }
}

int BarycentricActiveMesh2::get_macro_in_background_mesh(int k_active) const {
    assert(0 <= k_active <= nb_active_macros);
    return macro_in_background_mesh[k_active];
}
int BarycentricActiveMesh2::get_macro_in_active_mesh(int k_bg) const {
    auto it = macro_in_active_mesh.find(k_bg);
    if (it != macro_in_active_mesh.end()) {
        return it->second;
    }
    return -1; // not found
}


bool BarycentricActiveMesh2::is_macro_cut(int macro_k, int t /*=0*/) const {
    assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    const auto& sub = active_macro_elements[macro_k]; // 3 active micro ids
    for (int k_micro : sub)
        if (this->isCut(k_micro, t)) return true;
    return false;
}

// "interior micro" := !cut && !inactive
bool BarycentricActiveMesh2::is_macro_interior(int macro_k, int t /*=0*/) const {
    assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    const auto& sub = active_macro_elements[macro_k];
    for (int k_micro : sub)
        if (this->isCut(k_micro, t) || this->isInactive(k_micro, t)) return false;
    return true;
}

// "macro has no active micros" at time t
bool BarycentricActiveMesh2::is_macro_inactive(int macro_k, int t /*=0*/) const {
    assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    const auto& sub = active_macro_elements[macro_k];
    for (int k_micro : sub)
        if (!this->isInactive(k_micro, t)) return false;
    return true;
}

// Stabilize if macro is cut OR inactive at ANY time instance
bool BarycentricActiveMesh2::stabilize_macro(int macro_k) const {
    if (!(0 <= macro_k && macro_k < (int)active_macro_elements.size())) {
        std::cout << "macro_k = " << macro_k << " index out of range";
        assert(0);
    }
    // assert(0 <= macro_k && macro_k < (int)active_macro_elements.size());
    for (int t = 0; t < nb_quadrature_time_; ++t)
        if (is_macro_cut(macro_k, t) || is_macro_inactive(macro_k, t)) return true;
    return false;
}







// void BarycentricActiveMesh2::truncate(const Interface<Mesh2> &interface, int sign_domain_remove) {
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
//                 for (int sub_k : elements_in_macro) {
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
    
