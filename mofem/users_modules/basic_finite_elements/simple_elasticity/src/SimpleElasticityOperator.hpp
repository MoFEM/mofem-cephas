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

#ifndef __ELASTICITYNEWOPERATOR_HPP__
#define __ELASTICITYNEWOPERATOR_HPP__

#include"../../poisson/src/PoissonOperators.hpp"

namespace ElasticitySimpleExample {


  struct VolRule {
    int operator()(int,int,int p) const {
      return 2*(p-1);
    }
  };

  struct FaceRule {
    int operator()(int p_row,int p_col,int p_data) const {
      return 2*p_data+1;
    }
  };


  struct CreateFiniteElements{

    CreateFiniteElements(MoFEM::Interface &m_field):
      mField(m_field) {
    }


    /**
     * \brief Create finite element to calculate matrix and vectors
     */
    MoFEMErrorCode createFEToAssmbleMatrixAndVector(
      boost::function<double (const double,const double,const double)> f_u,
      boost::function<double (const double,const double,const double)> f_source,
      boost::shared_ptr<ForcesAndSourcesCore>& domain_lhs_fe,
      boost::shared_ptr<ForcesAndSourcesCore>& boundary_lhs_fe,
      boost::shared_ptr<ForcesAndSourcesCore>& domain_rhs_fe,
      boost::shared_ptr<ForcesAndSourcesCore>& boundary_rhs_fe,
      bool trans = true
    ) const {
      MoFEMFunctionBegin;

      // Create elements element instances
      domain_lhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      boundary_lhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(new FaceElementForcesAndSourcesCore(mField));
      domain_rhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(new VolumeElementForcesAndSourcesCore(mField));
      boundary_rhs_fe = boost::shared_ptr<ForcesAndSourcesCore>(new FaceElementForcesAndSourcesCore(mField));

      // Set integration rule to elements instances
      domain_lhs_fe->getRuleHook = VolRule();
      domain_rhs_fe->getRuleHook = VolRule();
      boundary_lhs_fe->getRuleHook = FaceRule();
      boundary_rhs_fe->getRuleHook = FaceRule();
     
      MoFEMFunctionReturn(0);
    }


  private:
    
    MoFEM::Interface &mField;
    

  };


}

#endif //__SIMPLEELASTICITYOPERATOR_HPP__
