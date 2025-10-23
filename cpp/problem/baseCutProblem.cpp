#include "baseProblem.hpp" 
#include "baseCutProblem.hpp" 
// #include <cassert>
// #include <cstddef>              // size_t

// If BaseCutFEM lives in a namespace, wrap the specialization with it:

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

//     int num_stab_faces = 0;

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

//             std::cout << "km = " << km << " is stabilized\n";

//             const auto& micro_elements = Th.active_macro_elements[km];

//             // Loop over all micro elements in the macro element
//             for (int k_micro : micro_elements) {

//                 // Loop through the element's edges
//                 for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

//                     int jfac = ifac;
//                     int kn   = Th.ElementAdj(k_micro, jfac); // get neighbor micro element's index

//                     // std::cout << "kn = " << kn << "\n";
//                     int kn_macro = Th.inverse_active_macro_map[kn];
//                     // std::cout << "kn_macro = " << kn_macro << "\n";
//                     if (kn == -1)
//                         continue;
//                     // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
//                     if ((kn < k_micro) && (Th.stabilize_macro(kn_macro)))
//                         continue;

//                     // std::cout << "Stabilizing face between active micro element " << k_micro << " with coordinates (" << Th[k_micro][0] << ", " << Th[k_micro][1] << ", " << Th[k_micro][2] << ") " <<  " and " << kn << " with coordinates (" << Th[kn][0] << ", " << Th[kn][1] << ", " << Th[kn][2] << ") " <<  "\n";
//                     num_stab_faces++;
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
//     // std::cout << "Number of stabilized faces: " << num_stab_faces / number_of_quadrature_points << "\n";
// }


// Stabilize edges corresponding to macro elements that are entirely OUTSIDE of the domain
// template <>
// void BaseCutFEM<Mesh2>::addPatchStabilizationExterior(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th, const TimeSlab &In) {

//     int number_of_quadrature_points = this->get_nb_quad_point_time();

//     int num_stab_faces = 0;

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

//             // Exclude elements who are not entirely outside of the active mesh
//             if (!Th.is_macro_exterior(km))
//                 continue;

//             std::cout << "km = " << km << " is exterior\n";

//             const auto& micro_elements = Th.active_macro_elements[km];

//             // Loop over all micro elements in the macro element
//             for (int k_micro : micro_elements) {

//                 // Loop through the element's edges
//                 for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

//                     int jfac = ifac;
//                     int kn   = Th.ElementAdj(k_micro, jfac); // get neighbor micro element's index

//                     // std::cout << "kn = " << kn << "\n";
//                     int kn_macro = Th.inverse_active_macro_map[kn];
//                     // std::cout << "kn_macro = " << kn_macro << "\n";
//                     if (kn == -1)
//                         continue;
//                     // By skipping neighbors with smaller indices, we avoid adding contribution to the same edge twice
//                     //if ((kn < k_micro) && (Th.stabilize_macro(kn_macro)))
//                     if ((kn < k_micro))
//                         continue;

//                     // std::cout << "Stabilizing face between active micro element " << k_micro << " with coordinates (" << Th[k_micro][0] << ", " << Th[k_micro][1] << ", " << Th[k_micro][2] << ") " <<  " and " << kn << " with coordinates (" << Th[kn][0] << ", " << Th[kn][1] << ", " << Th[kn][2] << ") " <<  "\n";
//                     num_stab_faces++;
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
//     // std::cout << "Number of stabilized faces: " << num_stab_faces / number_of_quadrature_points << "\n";
// }




template <> void BaseCutFEM<Mesh2>::addBilinearInner(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);

    for (int km = 0; km < Th.active_macro_elements.size(); ++km) {
        
        if (Th.is_macro_cut(km))
            continue;

        assert(Th.is_macro_interior(km));
        
        const auto& micro_elements = Th.active_macro_elements[km];

        // Loop over all micro elements in the macro element
        for (int k_micro : micro_elements) 
            BaseFEM<Mesh2>::addElementContribution(VF, k_micro, nullptr, 0, 1.);
        
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addLinearInner(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());
    progress bar(" Add Bilinear CutMesh", Th.last_element(), globalVariable::verbose);
    for (int km = 0; km < Th.active_macro_elements.size(); ++km) {
        
        if (Th.is_macro_cut(km))
            continue;

        assert(Th.is_macro_interior(km));
        
        const auto& micro_elements = Th.active_macro_elements[km];

        // Loop over all micro elements in the macro element
        for (int k_micro : micro_elements) 
            BaseFEM<Mesh2>::addElementContribution(VF, k_micro, nullptr, 0, 1.);
        
        this->addLocalContribution();
    }
    bar.end();
}

template <> void BaseCutFEM<Mesh2>::addBilinearInnerBorder(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());

    for (int km = 0; km < Th.active_macro_elements.size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements[km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

            // a neighboring element can be 1) outside of the domain, 2) another cut element, 
            // 3) inside of the domain but in the same macro element, or 4) inside of the domain and in another macro element
                //! THIS WILL BE PROBLEMATIC IF kn_micro == -1
                if ((Th.inverse_active_macro_map[k_micro] == Th.inverse_active_macro_map[kn_micro]) || (kn_micro == -1) || Th.isCut(kn_micro, 0))
                    continue;
                
                assert(!Th.isCut(kn_micro, 0));
                assert(Th.is_macro_interior(Th.inverse_active_macro_map[kn_micro]));

                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addInnerBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

        // Find the face that neighbors an interior element
        // for (int iface = 0; iface < Element::nea; ++iface) {
        //     int kn_macro = Th.macro_adjacent(km, iface);
            
        //     if ((kn_macro == -1) || Th.is_macro_cut(kn_macro))
        //         continue;
            
        //     std::cout << "Integrating on iface = " << iface << ", kn_macro = " << kn_macro << "\n";


        //     for (int k_micro : micro_elements) 
        //         BaseFEM<Mesh2>::addInnerBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
        // }
    //     
    //     if (!Th.isCut(k, 0)) {
    //         continue;
    //     } 
        
    //     // Find the face that neighbors an interior element
    //     for (int ifac = 0; ifac < Element::nea; ++ifac) { 
    //         int jfac = ifac;
    //         int kn   = Th.ElementAdj(k, jfac);

    //         // an neighboring element can be 1) outside of the domain, 2) another cut element or 3) inside of the domain 
    //         if ((kn == -1) || Th.isCut(kn, 0))
    //             continue;
            
    //         assert(!Th.isCut(kn, 0) && Th.isCut(k, 0));

    //         BaseFEM<Mesh2>::addInnerBorderContribution(VF, k, ifac, nullptr, 0, 1.);
            
    //     }
        
    //     this->addLocalContribution();
    }

}


template<>
void BaseCutFEM<Mesh2>::addBilinearOuterBorder(const itemVFlist_t& VF,
                                               const BarycentricActiveMesh2& Th)
{
    using Element = typename BarycentricActiveMesh2::Element;
    assert(!VF.isRHS());

    for (int km = 0; km < (int)Th.active_macro_elements.size(); ++km) {
        if (!Th.is_macro_cut(km)) continue;

        const auto& micro = Th.active_macro_elements[km];
        for (int k_micro : micro) {
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                const int kn_micro = Th.ElementAdj(k_micro, jfac);
                if (kn_micro != -1) continue;         // only outer boundary

                BaseFEM<Mesh2>::addOuterBorderContribution(VF, k_micro, ifac,
                                                           /*user*/nullptr, /*w*/0, /*coef*/1.0);
            }
        }
    }
}


template <> void BaseCutFEM<Mesh2>::addLinearOuterBorder(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(VF.isRHS());

    for (int km = 0; km < Th.active_macro_elements.size(); ++km) {

        // Loop only over cut elements for efficiency
        if (!Th.is_macro_cut(km))
            continue;

        const auto& micro_elements = Th.active_macro_elements[km];

        for (int k_micro : micro_elements) {
            
            for (int ifac = 0; ifac < Element::nea; ++ifac) {
                int jfac = ifac;
                int kn_micro = Th.ElementAdj(k_micro, jfac);

                // if ((kn_micro != -1) ||(Th.inverse_active_macro_map[k_micro] == Th.inverse_active_macro_map[kn_micro]))
                //     continue;
                if (kn_micro != -1)
                    continue;
                
                // std::cout << "Integrating on edge between micro elements " << k_micro << " (macro element " << Th.inverse_active_macro_map[k_micro] << ") and " << kn_micro << " (macro element " << Th.inverse_active_macro_map[kn_micro] << ")\n";
                BaseFEM<Mesh2>::addOuterBorderContribution(VF, k_micro, ifac, nullptr, 0, 1.);
            
            }
        
        }

    }

}







// template <> void BaseCutFEM<Mesh2>::addElementStabilization(const itemVFlist_t &VF, const Interface<Mesh2> &interface, const CutMesh &Th) {
//     assert(!VF.isRHS());
//     progress bar("Add Bilinear Mesh", Th.last_element(), globalVariable::verbose);
// #pragma omp parallel for num_threads(this->get_num_threads())
//     for (int k = Th.first_element(); k < Th.last_element(); k += Th.next_element()) {
//         bar += Th.next_element();
        
//         // Compute parameter coonected to the mesh.
//         // on can take the one from the first test function
//         const FESpace &Vh(*VF[0].fespaceV);
//         const FElement &FK(Vh[k]);
//         const Element &K(FK.T);
//         double meas = K.measure();
//         double h    = K.get_h();
//         int domain  = FK.get_domain();
//         int kb      = Vh.idxElementInBackMesh(k);
//         int iam     = omp_get_thread_num();

//         // GET THE QUADRATURE RULE
//         const QF &qf(this->get_quadrature_formular_K());
    
//         double tid = 0.;

//         // LOOP OVER THE VARIATIONAL FORMULATION ITEMS
//         for (int l = 0; l < VF.size(); ++l) {
//             if (!VF[l].on(domain))
//                 continue; 

//             // FINTE ELEMENT SPACES && ELEMENTS
//             const FESpace &Vhv(VF.get_spaceV(l));
//             const FESpace &Vhu(VF.get_spaceU(l));
//             const FElement &FKv(Vhv[k]);
//             const FElement &FKu(Vhu[k]);
//             this->initIndex(FKu, FKv); //, iam);
//             // BF MEMORY MANAGEMENT -
//             bool same   = (&Vhu == &Vhv);
//             int lastop  = getLastop(VF[l].du, VF[l].dv);
//             // RNMK_ fv(this->databf_,FKv.NbDoF(),FKv.N,lastop); //  the value for
//             // basic fonction RNMK_ fu(this->databf_+ (same
//             // ?0:FKv.NbDoF()*FKv.N*lastop) ,FKu.NbDoF(),FKu.N,lastop); //  the value
//             // for basic fonction
//             long offset = iam * this->offset_bf_;
//             RNMK_ fv(this->databf_ + offset, FKv.NbDoF(), FKv.N,
//                     lastop); //  the value for basic fonction
//             RNMK_ fu(this->databf_ + offset + (same ? 0 : FKv.NbDoF() * FKv.N * lastop), FKu.NbDoF(), FKu.N, lastop);
//             What_d Fop = Fwhatd(lastop);

//             // COMPUTE COEFFICIENT
//             double coef = VF[l].computeCoefElement(h, meas, meas, meas, domain);

//             // LOOP OVER QUADRATURE IN SPACE
//             for (int ipq = 0; ipq < qf.getNbrOfQuads(); ++ipq) {
//                 typename QF::QuadraturePoint ip(qf[ipq]);
//                 const Rd mip = K.mapToPhysicalElement(ip);
//                 double Cint  = meas * ip.getWeight() * cst_time;

//                 // EVALUATE THE BASIS FUNCTIONS
//                 FKv.BF(Fop, ip, fv);
//                 if (!same)
//                     FKu.BF(Fop, ip, fu);
//                 //   VF[l].applyFunNL(fu,fv);

//                 // FIND AND COMPUTE ALL THE COEFFICENTS AND PARAMETERS
//                 Cint *= VF[l].evaluateFunctionOnBackgroundMesh(kb, domain, mip, tid);
//                 Cint *= coef * VF[l].c;

//                 if (In) {
//                     if (VF.isRHS())
//                         this->addToRHS(VF[l], *In, FKv, fv, Cint);
//                     else
//                         this->addToMatrix(VF[l], *In, FKu, FKv, fu, fv, Cint);
//                 } else {
//                     if (VF.isRHS())
//                         this->addToRHS(VF[l], FKv, fv, Cint);
//                     else
//                         this->addToMatrix(VF[l], FKu, FKv, fu, fv, Cint); //, iam);
//                 }
//             }
//         }

//         this->addLocalContribution();
//     }
//     bar.end();
// }