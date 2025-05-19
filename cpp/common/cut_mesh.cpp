#include "time_interface.hpp"
#include "cut_mesh.hpp"



BarycentricActiveMesh::BarycentricActiveMesh(const BarycentricMesh2 &th)
    : ActiveMesh<Mesh2>(static_cast<const Mesh2 &>(th)) {}
    
void BarycentricActiveMesh::truncate(const Interface<Mesh2> &interface, int sign_domain_remove) {
    int dom_size = this->get_nb_domain();
    this->idx_element_domain.clear();

    // Clear background ↔ active mesh index structures
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

            for (int d = 0; d < dom_size; ++d) {
                this->idx_in_background_mesh_[d].push_back(sub_k);
                this->idx_from_background_mesh_[d][sub_k] = nt[d];

                if (interface.isCut(sub_k)) {
                    this->interface_id_[0][{d, nt[d]}].emplace_back(&interface, -sign_domain_remove);
                }

                ++nt[d];
            }
        }
    }

    this->idx_element_domain.push_back(0);
    for (int d = 0; d < dom_size; ++d) {
        this->idx_in_background_mesh_[d].shrink_to_fit();
        this->idx_element_domain.push_back(this->idx_element_domain.back() + nt[d]);
    }
}
    
