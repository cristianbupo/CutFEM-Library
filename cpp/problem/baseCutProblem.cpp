/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#include "baseCutProblem.hpp"

template <>
void BaseCutFEM<BarycentricMesh2>::addPatchStabilization(const itemVFlist_t &VF, const BarycentricActiveMesh2 &Th) {
    assert(!VF.isRHS());
    progress bar(" Add Patch Stabilization BarycentricActiveMesh", Th.last_element(), globalVariable::verbose);

    size_t num_stab_faces = 0;

    // Loop over macro elements
    for (std::size_t macro_k = 0; macro_k < Th.macro_elements.size(); ++macro_k) {
        bar += Th.next_element();

        if (!Th.isCut(k, 0) && !Th.isInactive(k, 0))
            continue;
        for (int ifac = 0; ifac < Element::nea; ++ifac) { // loop over the edges / faces

            int jfac = ifac;
            int kn   = Th.ElementAdj(k, jfac);
            // ONLY INNER EDGE && LOWER INDEX TAKE CARE OF THE INTEGRATION
            if (kn < k)
                continue;
        

            // std::pair<int, int> e1 = std::make_pair(k, ifac);
            // std::pair<int, int> e2 = std::make_pair(kn, jfac);
            BaseFEM<M>::addPatchContribution(VF, k, kn, nullptr, 0, 1.);
            num_stab_faces++;
        }
        this->addLocalContribution();
    }
    bar.end();
    std::cout << "Number of stabilized faces: " << num_stab_faces << "\n";
}


// template class CutFEM<BarycentricMesh2>;
// template class FEM<BarycentricMesh2>;

