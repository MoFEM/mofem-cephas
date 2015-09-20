/** \file NetGenInterface.hpp
 * \brief NetGen interface 
 * 
 * MoFEM interface for NetGen
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __NETGENINTERFACE_HPP__
#define __NETGENINTERFACE_HPP__

#include "FieldUnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMNetGegInterface = MOFEMuuid( BitIntefaceId(NETGEN_INTERFACE) );

/** \brief Use NetGen to generate mesh
  * \ingroup mofem
  */
struct NetGenInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  NetGenInterface(MoFEM::Core& core): cOre(core) {};

  PetscErrorCode stlSetSurfaceTriangles(Ng_STL_Geometry *stl_geom,Range &ents,double *nv = NULL,int verb = 0);
  PetscErrorCode stlSetSurfaceEdges(Ng_STL_Geometry *stl_geom,Range &ents,int verb = 0);

  PetscErrorCode setPoints(Ng_Mesh *mesh,vector<EntityHandle> &pts,int verb = 0);
  PetscErrorCode setSurfaceElements(Ng_Mesh *mesh,vector<EntityHandle> &pts,vector<EntityHandle> &elms,Range *tets = NULL,int verb = 0);

  PetscErrorCode getPoints(Ng_Mesh *mesh,vector<EntityHandle> &pts);
  PetscErrorCode getSurfaceElements(Ng_Mesh *mesh,vector<EntityHandle> &pts,vector<EntityHandle> &elms);
  PetscErrorCode getVolumeElements(Ng_Mesh *mesh,vector<EntityHandle> &pts,vector<EntityHandle> &elms);

};

}

#endif //__NETGENINTERFACE_HPP__
