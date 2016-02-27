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

const int prism_adj_edges[] = { 6,7,8, -1,-1,-1, 0,1,2 };
const int prism_edges_conn[6][2] = { {0,1},{1,2},{2,0}, {3,4}, {4,5}, {5,3} };

inline int fNBENTITY_GENERIC(int P) { (void)P; return 0; }
inline int fNBENTITYSET_NOFIELD(int P) { (void)P; return 1; }
//
inline int fNBVERTEX_L2(int P) { (void)P; return 0; }
inline int fNBVOLUMETET_L2_AINSWORTH_COLE(int P) { return NBVOLUMETET_L2_AINSWORTH_COLE(P); }
inline int fNBFACETRI_L2_AINSWORTH_COLE(int P) { return NBFACETRI_L2_AINSWORTH_COLE(P); }
inline int fNBEDGE_L2_AINSWORTH_COLE(int P) { return NBEDGE_L2_AINSWORTH_COLE(P); }

//
/// number of approx. functions for H1 space on vertex
inline int fNBVERTEX_H1_AINSWORTH_COLE(int P) { return (P==1) ? 1 : 0; }
/// number of approx. functions for H1 space on edge
inline int fNBEDGE_H1_AINSWORTH_COLE(int P) { return NBEDGE_H1_AINSWORTH_COLE(P); }
/// number of approx. functions for H1 space on face
inline int fNBFACETRI_H1_AINSWORTH_COLE(int P) { return NBFACETRI_H1_AINSWORTH_COLE(P); }
inline int fNBFACEQUAD_H1_AINSWORTH_COLE(int P) { return NBFACEQUAD_H1_AINSWORTH_COLE(P); }
/// number of approx. functions for H1 space on volume
inline int fNBVOLUMETET_H1_AINSWORTH_COLE(int P) { return NBVOLUMETET_H1_AINSWORTH_COLE(P); }
inline int fNBVOLUMEPRISM_H1_AINSWORTH_COLE(int P) { return NBVOLUMEPRISM_H1_AINSWORTH_COLE(P); }

//
/// number of approx. functions for HCURL space on vertex
inline int fNBVERTEX_HCURL_AINSWORTH_COLE(int P) { (void)P; return 0; }
inline int fNBEDGE_HCURL_AINSWORTH_COLE(int P) { return NBEDGE_HCURL_AINSWORTH_COLE(P); }
inline int fNBFACETRI_HCURL_AINSWORTH_COLE(int P) { return NBFACETRI_HCURL_AINSWORTH_COLE(P); }
inline int fNBVOLUMETET_HCURL_AINSWORTH_COLE(int P) { return NBVOLUMETET_HCURL_AINSWORTH_COLE(P); }

/// \brief number of approx. functions for HDIV space on vertex
///
/// zero number of digrees of freedom on vertex for that space
inline int fNBVERTEX_HDIV_AINSWORTH_COLE(int P) { (void)P; return 0; }
/// number of approx. functions for HDIV space on edge
inline int fNBEDGE_HDIV_AINSWORTH_COLE(int P) { assert(P==P); (void)P; return NBEDGE_HDIV_AINSWORTH_COLE(P); }
/// number of approx. functions for HDIV space on face
inline int fNBFACETRI_HDIV_AINSWORTH_COLE(int P) { return NBFACETRI_HDIV_AINSWORTH_COLE(P); }
/// number of approx. functions for HDIV space on voulem
inline int fNBVOLUMETET_HDIV_AINSWORTH_COLE(int P) { return NBVOLUMETET_HDIV_AINSWORTH_COLE(P); }

/** \brief Test MoAB entity handle if has structure as is assumed by MoFEM
  */
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent);

}

#endif //__DATASTRUCTURES_HPP__
