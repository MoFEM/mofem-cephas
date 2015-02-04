/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

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

struct CohesiveInterfaceElement {

  struct BlockData {
    double h;
    double E0;
    double ft;
    double Gf;
    Range tRisms; 
  }; 
  map<int,BlockData> setOfBlocks; 

  struct CommonData {
    ublas::matrix<double> gAp;
  };
  CommonData commonData;

  CohesiveInterfaceElement()

  struct OpCalculateGap: public FlatPrismElementForcesAndSurcesCore::UserDataOperator {

    CommonData &commonData;
    OpCalculateGap(const string field_name,commonData &common_data):
      FlatPrismElementForcesAndSurcesCore::UserDataOperator(field_name),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

 	int nb_dofs = data.getFieldData().size();
	if(nb_dofs == 0) {
	  PetscFunctionReturn(0);
	}
	int nb_gauss_pts = data.getN().size1();

	if(type == MBVERTEX) {
	  commonData.gAp.resize(nb_gauss_pts,3);
	  commonData.gAp.clear();
	}

	for(int gg = 0;gg<nb_gauss_pts;gg++) {
	  for(int dd = 0;dd<3;dd++) {
	    commonData.gAp(gg,dd) += cblas_ddot(3,data.getFieldData()

	  }
	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "throw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };


};


