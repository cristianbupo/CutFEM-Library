#include "baseProblem.hpp" 
#include "baseCutProblem.hpp" 
// #include <cassert>
// #include <cstddef>              // size_t

// If BaseCutFEM lives in a namespace, wrap the specialization with it:

#include "baseCutProblem.hpp"

// specialization for BarycentricMesh2
// template<>
// void BaseCutFEM<Mesh2>::addPatchStabilization(const itemVFlist_t &VF,
//                                           const BarycentricActiveMesh2 &Th)
// {
//     using Element = typename BarycentricActiveMesh2::Element;

    
//     assert(!VF.isRHS());

//     size_t num_stab_faces = 0;

//     // Helper: a “micro is exterior” test. Adjust to your flags if needed.
//     auto is_micro_exterior = [&](int kmicro) -> bool {
//         // If you have a dedicated API, prefer that:
//         // return Th.is_micro_exterior(kmicro);
//         // Fallback using your usual flags:
//         return Th.isInactive(kmicro, /*level*/0) && !Th.isCut(kmicro, /*level*/0);
//     };

//     for (int kept_macro_id = 0; kept_macro_id < (int)Th.macro_elements_amk_.size(); ++kept_macro_id) {

//         if (!Th.is_macro_cut(kept_macro_id, 0))
//             continue;

//         const auto &micro_ids = Th.macro_elements_amk_[kept_macro_id]; // 3 active element IDs

//         for (int r = 0; r < 3; ++r) {
//             int k = micro_ids[r];
//             if (k < 0) continue; // was not kept in active mesh

//             for (int ifac = 0; ifac < Element::nea; ++ifac) {
//                 int jfac = ifac;
//                 int kn   = Th.ElementAdj(k, jfac);

//                 // skip edges processed from other side
//                 // if (kn >= 0 && k > kn) continue;
//                 if (kn < k)
//                     continue;

//                 // check exclusion rule: skip if exactly one is outside
//                 bool k_out  = Th.isInactive(k, 0) && !Th.isCut(k, 0);
//                 bool kn_out = (kn >= 0) ? (Th.isInactive(kn, 0) && !Th.isCut(kn, 0)) : false;

//                 if (k_out ^ kn_out) continue;
//                 std::cout << "Here\n";
//                 BaseFEM<Mesh2>::addPatchContribution(VF, k, kn, nullptr, 0, 1.0);
//                 ++num_stab_faces;
//             }
//         }

//         this->addLocalContribution();
//     }


//     std::cout << "Number of stabilized faces (macro patches): " << num_stab_faces << "\n";
// }


template<>
void BaseCutFEM<Mesh2>::addPatchStabilization(const itemVFlist_t &VF,
                                          const BarycentricActiveMesh2 &Th)
{
    using Element = typename BarycentricActiveMesh2::Element;

    
    assert(!VF.isRHS());

    size_t num_stab_faces = 0;
    for (int km = 0; km < Th.active_macro_elements.size(); ++km) {

        if (!Th.is_macro_cut(km))
            continue;
        
        const auto& micro_elements = Th.active_macro_elements[km];

        // Loop over all micro elements in the macro element
        for (int k_micro : micro_elements) {

            // Loop over all the faces of the micro element
            for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

                int jfac = ifac;
                int kn   = Th.ElementAdj(k_micro, jfac);
                if (kn == -1)   
                    continue;

                int kn_macro = Th.inverse_active_macro_map[kn];
                // Apply "lower index only" logic only for cut elements, otherwise we might exclude the edges going to interior elements
                if ((kn < k_micro) && Th.is_macro_cut(kn_macro)) 
                    continue;
            
                BaseFEM<Mesh2>::addPatchContribution(VF, k_micro, kn, nullptr, 0, 1.);
                num_stab_faces++;
            }
        }
        
        this->addLocalContribution();
    }
    // std::cout << "Number of STABILIZED faces: " << num_stab_faces << "\n";
}


// template <>
// void BaseCutFEM<Mesh2>::addPatchStabilization(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th, const TimeSlab &In) {

//     int number_of_quadrature_points = this->get_nb_quad_point_time();

//     // Loop through time quadrature points
//     for (int itq = 0; itq < number_of_quadrature_points; ++itq) {
//         assert(!VF.isRHS());

//         // Compute contribution from time basis functions
//         auto tq    = this->get_quadrature_time(itq);
//         double tid = In.map(tq);
//         KNMK<double> basisFunTime(In.NbDoF(), 1, op_dz + 1);
//         RNMK_ bf_time(this->databf_time_, In.NbDoF(), 1, op_dz);
//         In.BF(tq.x, bf_time); // compute time basic funtions
//         double cst_time = tq.a * In.get_measure();

//         // Loop through active macro elements
//         for (int km = 0; km < Th.active_macro_elements.size(); ++km) {

//             // Exclude elements whose edges do not need stabilization
//             if (!Th.stabilize_macro(km))
//                 continue;

//             const auto& micro_elements = Th.active_macro_elements[km];

//             // Loop over all micro elements in the macro element
//             for (int k_micro : micro_elements) {

//                 // Loop through the element's edges
//                 for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

//                     int jfac = ifac;
//                     int kn   = Th.ElementAdj(k_micro, jfac); // get neighbor micro element's index
//                     int kn_macro = Th.inverse_active_macro_map[kn];
//                     // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
//                     if ((kn < k_micro) && (Th.stabilize_macro(kn_macro, itq)))
//                         continue;

//                     std::pair<int, int> e1 = std::make_pair(k_micro, ifac);  // (element index, edge index) current element
//                     std::pair<int, int> e2 = std::make_pair(kn, jfac); // (element index, edge index) neighbor element

//                     // Add patch contribution
//                     // BaseFEM<M>::addFaceContribution(VF, e1, e2, &In, itq, cst_time);
//                     BaseFEM<Mesh2>::addPatchContribution(VF, k_micro, kn, &In, itq, cst_time);
//                 }
//                 this->addLocalContribution();
//             }
//         }
//     }
// }