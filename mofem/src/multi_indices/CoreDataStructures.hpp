/** \file CoreDataStructures.hpp
 * \brief Myltindex containers, data structures and other low-level functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#ifndef __DATASTRUCTURES_HPP__
#define __DATASTRUCTURES_HPP__

namespace MoFEM {

extern const char *FieldSpaceNames[];

const int prism_adj_edges[] = { 6,7,8, -1,-1,-1, 0,1,2 };
const int prism_edges_conn[6][2] = { {0,1},{1,2},{2,0}, {3,4}, {4,5}, {5,3} };

inline int fNBENTITYSET_NOFIELD(int P) { (void)P; return 1; }
//
inline int fNBVERTEX_L2(int P) { (void)P; return 0; }
inline int fNBVOLUMETET_L2(int P) { return NBVOLUMETET_L2(P); }
inline int fNBFACETRI_L2(int P) { return NBFACETRI_L2(P); }
inline int fNBEDGE_L2(int P) { return NBEDGE_L2(P); }

//
/// number of approx. functions for H1 space on vertex
inline int fNBVERTEX_H1(int P) { return (P==1) ? 1 : 0; }
/// number of approx. functions for H1 space on edge
inline int fNBEDGE_H1(int P) { return NBEDGE_H1(P); }
/// number of approx. functions for H1 space on face
inline int fNBFACETRI_H1(int P) { return NBFACETRI_H1(P); }
inline int fNBFACEQUAD_H1(int P) { return NBFACEQUAD_H1(P); }
/// number of approx. functions for H1 space on volume
inline int fNBVOLUMETET_H1(int P) { return NBVOLUMETET_H1(P); }
inline int fNBVOLUMEPRISM_H1(int P) { return NBVOLUMEPRISM_H1(P); }

//
/// number of approx. functions for HCURL space on vertex
inline int fNBVERTEX_HCURL(int P) { (void)P; return 0; }
inline int fNBEDGE_HCURL(int P) { return NBEDGE_HCURL(P); }
inline int fNBFACETRI_HCURL(int P) { return NBFACETRI_HCURL(P); }
inline int fNBVOLUMETET_HCURL(int P) { return NBVOLUMETET_HCURL(P); }
//
/// \brief number of approx. functions for HDIV space on vertex
///
/// zero number of digrees of freedom on vertex for that space
inline int fNBVERTEX_HDIV(int P) { (void)P; return 0; }
/// number of approx. functions for HDIV space on edge
inline int fNBEDGE_HDIV(int P) { assert(P==P); (void)P; return NBEDGE_HDIV(P); }
/// number of approx. functions for HDIV space on face
inline int fNBFACETRI_HDIV(int P) { return NBFACETRI_HDIV(P); }
/// number of approx. functions for HDIV space on voulem
inline int fNBVOLUMETET_HDIV(int P) { return NBVOLUMETET_HDIV(P); }

struct ltstr
{ inline bool operator()(const string &s1, const string& s2) const
  { return strcmp(s1.c_str(), s2.c_str()) < 0; } };


template<typename Tag>
void get_vector_by_multi_index_tag(vector<DofMoFEMEntity> &vec_dof,const DofMoFEMEntity_multiIndex &dofs,Tag* = 0) {
  const typename boost::multi_index::index<DofMoFEMEntity_multiIndex,Tag>::type& i = get<Tag>(dofs);
  vec_dof.insert(vec_dof.end(),i.begin(),i.end());
}

/** \brief if moab entity hanled
  */
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent);

}

#endif //__DATASTRUCTURES_HPP__
