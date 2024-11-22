/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

/** \file SeepageOps.hpp
 * \example SeepageOps.hpp
 */

namespace SeepageOps {

template <int DIM>
struct OpDomainRhsHydrostaticStress
    : public AssemblyDomainEleOp { // changed opfaceele to AssemblyDomainEleOp
public:
  OpDomainRhsHydrostaticStress(std::string field_name1,
                               boost::shared_ptr<VectorDouble> p_ptr)
      : AssemblyDomainEleOp(field_name1, field_name1, DomainEleOp::OPROW),
        pPtr(p_ptr) {}

  MoFEMErrorCode iNtegrate(DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;

    const int nb_dofs = data.getIndices().size();

    if (nb_dofs) {
      
      auto &nf = AssemblyDomainEleOp::locF; //access to local vector
      // get element measure
      const double measure = getMeasure();

      // get number of integration points
      const int nb_integration_points = getGaussPts().size2();
      // get integration weights
      auto t_w = getFTensor0IntegrationWeight();

      // get base function
      auto t_base_diff = data.getFTensor1DiffN<DIM>();

      //constexpr double g_acceleration = 9.81;

      FTensor::Index<'i', DIM> i;

      auto t_p = getFTensor0FromVec(*pPtr);
      for (int gg = 0; gg != nb_integration_points; gg++) { //loop over Gauss integration points
        auto t_nf = getFTensor1FromPtr<DIM>(&nf[0]);

        const double a = t_w * measure * t_p;

        for (int rr = 0; rr != nb_dofs / DIM; rr++) {  //loop over dofs (rr = alpha)
          t_nf(i) -= t_base_diff(i) * a;  //loop over dimensions

          // move to the next base function
          ++t_base_diff; // moves the pointer to the next shape function
          ++t_nf;
        }

        // move to the weight of the next integration point
        ++t_w;
        ++t_p;
      }

    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> pPtr;
};

} // namespace SeepageOps