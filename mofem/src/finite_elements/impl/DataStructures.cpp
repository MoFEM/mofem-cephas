/** \file DataStructures.cpp

\brief Implementation for Data Structures in Forces and Sources

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

#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <ForcesAndSurcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

template<class T>
void cOnstructor(DataForcesAndSurcesCore *data,EntityType type,T) {

  data->dataOnEntities[MBENTITYSET].push_back(new T());

  switch (type) {
    case MBTET:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for(int ee = 0;ee<6;ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    for(int ff = 0;ff<4;ff++) {
      data->dataOnEntities[MBTRI].push_back(new T());
    }
    data->dataOnEntities[MBTET].push_back(new T());
    break;
    case MBTRI:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for(int ee = 0;ee<3;ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    data->dataOnEntities[MBTRI].push_back(new T());
    break;
    case MBEDGE:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    data->dataOnEntities[MBEDGE].push_back(new T());
    break;
    case MBVERTEX:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    break;
    case MBPRISM:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for(int ee = 0;ee<9;ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    for(int ff = 0;ff<5;ff++) {
      data->dataOnEntities[MBQUAD].push_back(new T());
    }
    for(int ff = 0;ff<5;ff++) {
      data->dataOnEntities[MBTRI].push_back(new T());
    }
    data->dataOnEntities[MBPRISM].push_back(new T());
    break;
    default:
    throw MoFEMException(MOFEM_NOT_IMPLEMENTED);
  }

}

DataForcesAndSurcesCore::DataForcesAndSurcesCore(EntityType type) {
  cOnstructor(this,type,EntData());
}


DerivedDataForcesAndSurcesCore::DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data): DataForcesAndSurcesCore() {

  boost::ptr_vector<EntData>::iterator iit;

  boost::ptr_vector<EntData>::iterator it;
  for(it = data.dataOnEntities[MBVERTEX].begin();it!=data.dataOnEntities[MBVERTEX].end();it++) {
    dataOnEntities[MBVERTEX].push_back(new DerivedEntData(*it));
  }
  for(it = data.dataOnEntities[MBEDGE].begin();it!=data.dataOnEntities[MBEDGE].end();it++) {
    dataOnEntities[MBEDGE].push_back(new DerivedEntData(*it));
  }
  for(it = data.dataOnEntities[MBTRI].begin();it!=data.dataOnEntities[MBTRI].end();it++) {
    dataOnEntities[MBTRI].push_back(new DerivedEntData(*it));
  }
  for(it = data.dataOnEntities[MBQUAD].begin();it!=data.dataOnEntities[MBQUAD].end();it++) {
    dataOnEntities[MBQUAD].push_back(new DerivedEntData(*it));
  }
  for(it = data.dataOnEntities[MBTET].begin();it!=data.dataOnEntities[MBTET].end();it++) {
    dataOnEntities[MBTET].push_back(new DerivedEntData(*it));
  }
  for(it = data.dataOnEntities[MBPRISM].begin();it!=data.dataOnEntities[MBPRISM].end();it++) {
    dataOnEntities[MBPRISM].push_back(new DerivedEntData(*it));
  }
}

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e) {
  os <<
    "sEnse: " << e.getSense() << endl <<
    "oRder: " << e.getOrder() << endl <<
    "iNdices: " << e.getIndices() << endl;
  os.precision(2);
  os <<
    "fieldData: " << std::fixed << e.getFieldData() << endl;
  os <<
    "N: " << std::fixed << e.getN() << endl <<
    "diffN: " << std::fixed << e.getDiffN();
  return os;
}

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e) {
  for(unsigned int nn = 0;nn < e.dataOnEntities[MBVERTEX].size(); nn++) {
    os << "dataOnEntities[MBVERTEX][" << nn << "]" << endl << e.dataOnEntities[MBVERTEX][nn] << endl;
  }
  for(unsigned int ee = 0;ee < e.dataOnEntities[MBEDGE].size(); ee++) {
    os << "dataOnEntities[MBEDGE][" << ee << "]" << endl << e.dataOnEntities[MBEDGE][ee] << endl;
  }
  for(unsigned int ff = 0;ff < e.dataOnEntities[MBTRI].size(); ff++) {
    os << "dataOnEntities[MBTRI][" << ff << "] " << endl << e.dataOnEntities[MBTRI][ff] << endl;
  }
  for(unsigned int vv = 0;vv < e.dataOnEntities[MBTET].size(); vv++) {
    os << "dataOnEntities[MBTET][" << vv << "]" << endl << e.dataOnEntities[MBTET][vv] << endl;
  }
  return os;
}

}
