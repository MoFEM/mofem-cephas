/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __BCSRVEVOLUME_HPP
#define __BCSRVEVOLUME_HPP

namespace MoFEM {
  struct BCs_RVEVolume {
    
    /// \brief  definition of volume element
    struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {
      MyVolumeFE(FieldInterface &_mField): VolumeElementForcesAndSourcesCore(_mField) {}
    };
    
    MyVolumeFE feLhs; //< calculate left hand side for tetrahedral elements
    MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element
    
    FieldInterface &mField;
    BCs_RVEVolume(FieldInterface &m_field):
    feLhs(m_field),
    mField(m_field) {}
    
    /** \biref operator to calculate left hand side of heat conductivity terms
     * \infroup mofem_thermal_elem
     */
    struct OpVolumeCal: public VolumeElementForcesAndSourcesCore::UserDataOperator {
      Vec RVE_volume_Vec;
      NonlinearElasticElement::BlockData &dAta;
      FieldInterface &mField;
      
      OpVolumeCal(FieldInterface &m_field, string field_name, Vec _RVE_volume_Vec, NonlinearElasticElement::BlockData &data):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name, UserDataOperator::OPROWCOL),
      mField(m_field), RVE_volume_Vec(_RVE_volume_Vec), dAta(data) { }
      
      PetscErrorCode doWork(
                            int row_side,int col_side,
                            EntityType row_type,EntityType col_type,
                            DataForcesAndSurcesCore::EntData &row_data,
                            DataForcesAndSurcesCore::EntData &col_data) {
        PetscFunctionBegin;
        if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
          PetscFunctionReturn(0);
        }
        
        //      cout<<"Hi from OpVolumeCal "<<endl;
        if(row_type == MBVERTEX && col_type==MBVERTEX) {
//        cout<<"Hi from MBVERTEX "<<endl;
          ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
          int Indices[1];  Indices[0]=pcomm->rank();
          double Vol_elm[1];  Vol_elm[0]=0;
          for(unsigned int gg = 0;gg<row_data.getN().size1();gg++) {
            if(getHoGaussPtsDetJac().size()>0) {
//              cout<<"getHoGaussPtsDetJac()[gg] "<<getHoGaussPtsDetJac()[gg]<<endl;
//              cout<<"High order geometry "<<endl;
              Vol_elm[0]+=getVolume()*getGaussPts()(3,gg)*getHoGaussPtsDetJac()[gg];
            }else{
//              cout<<"Low order geometry "<<endl;
              Vol_elm[0]+=getVolume()*getGaussPts()(3,gg);
            }
          }
          PetscErrorCode ierr;
          ierr = VecSetValues(RVE_volume_Vec,1,Indices,Vol_elm,ADD_VALUES); CHKERRQ(ierr);
//          cout<<"Indices[0] "<<Indices[0] << endl;
//          cout<<"Vol_elm[0] "<<Vol_elm[0] << endl;
        }
        PetscFunctionReturn(0);
      }
    };
  
  
    PetscErrorCode setRVEVolumeOperators(FieldInterface &mField, string field_name, Vec _RVE_volume_Vec, map<int,NonlinearElasticElement::BlockData> &setOfBlocks) {
      PetscFunctionBegin;
      
      ////    cout<<"Hi from setRVEVolumeOperators "<<endl;
      map<int,NonlinearElasticElement::BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
//        cout<<"Hi from loop "<<endl;
//        cout<<"sit->second.tEts ===   "<<sit->second.tEts<<endl;
        feLhs.getOpPtrVector().push_back(new OpVolumeCal(mField, field_name, _RVE_volume_Vec, sit->second));
      }
    }
    
  };
  
}

#endif //__RVEVolume_HPP__
