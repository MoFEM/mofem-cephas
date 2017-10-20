/** \file Tools.cpp
 * \brief Auxilairy tools
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

namespace MoFEM {

  PetscErrorCode Tools::query_interface(const MOFEMuuid& uuid, UnknownInterface** iface) const {
    MoFEMFunctionBeginHot;
    *iface = NULL;
    if(uuid == IDD_MOFEMNodeMerger) {
      *iface = const_cast<Tools*>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    MoFEMFunctionReturnHot(0);
  }

  double Tools::volumeLengthQuality(const double *coords) {
    double lrms = 0;
    for(int dd = 0;dd!=3;dd++) {
      lrms +=
      pow(coords[0*3+dd]-coords[1*3+dd],2)+
      pow(coords[0*3+dd]-coords[2*3+dd],2)+
      pow(coords[0*3+dd]-coords[3*3+dd],2)+
      pow(coords[1*3+dd]-coords[2*3+dd],2)+
      pow(coords[1*3+dd]-coords[3*3+dd],2)+
      pow(coords[2*3+dd]-coords[3*3+dd],2);
    }
    lrms = sqrt((1./6.)*lrms);
    double diff_n[12];
    ShapeDiffMBTET(diff_n);
    FTensor::Tensor1<double*,3> t_diff_n(&diff_n[0],&diff_n[1],&diff_n[2],3);
    FTensor::Tensor1<const double*,3> t_coords(&coords[0],&coords[1],&coords[2],3);
    FTensor::Tensor2<double,3,3> jac;
    FTensor::Index<'i',3> i;
    FTensor::Index<'j',3> j;
    jac(i,j) = 0;
    for(int nn = 0;nn!=4;nn++) {
      jac(i,j) += t_coords(i)*t_diff_n(j);
      ++t_coords;
      ++t_diff_n;
    }
    double volume = dEterminant(jac)/6.;
    return 6.*sqrt(2.)*volume/pow(lrms,3);
  }

  PetscErrorCode Tools::minTetsQuality(const Range& tets,double &min_quality) {
    MoFEM::Interface& m_field = cOre;
    moab::Interface& moab(m_field.get_moab());
    MoFEMFunctionBeginHot;
    const EntityHandle* conn;
    int num_nodes;
    double coords[12];
    for(Range::iterator tit=tets.begin();tit!=tets.end();tit++) {
      rval = m_field.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      rval = moab.get_coords(conn,num_nodes,coords); CHKERRQ_MOAB(rval);
      double q = Tools::volumeLengthQuality(coords);
      min_quality = fmin(q,min_quality);
    }
    MoFEMFunctionReturnHot(0);
  }

}
