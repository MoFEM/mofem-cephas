/** \file Simple.hpp
 * \brief NodeMerger interface
 *
 * Make simplified interface, to speedup problem setup and analysts.
 * See discussion here
 * <a href=https://groups.google.com/d/msg/mofem-group/Vkc00aia4dU/o9RF3ZmPAAAJ>link to google groups</a>
 *
 * \ingroup mofem_node_merger
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __SIMPLE_HPP__
#define __SIMPLE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMSimple = MOFEMuuid( BitIntefaceId(SIMPLE_INTERFACE) );

  struct Simple: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    Simple(const MoFEM::Core& core);

    /**
    * \brief Destructor
    */
    ~Simple();

    /**
     * \brief Add field
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addDomainField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_EXCL,
      int verb = -1
    );

    /**
     * \brief Add field
     * @param  name              name of the filed
     * @param  space             space (L2,H1,Hdiv,Hcurl)
     * @param  base              approximation base, see FieldApproximationBase
     * @param  nb_of_cooficients number of field coefficients
     * @param  tag_type          type of the tag MB_TAG_DENSE or MB_TAG_SPARSE (DENSE is faster and uses less memory, SPARSE is more flexible if you define field on subdomains)
     * @param  bh                if MF_EXCL throws error if field exits, MF_ZERO no error if field exist
     * @param  verb              verbosity leve
     * @return                   error code
     */
    PetscErrorCode addBoundaryField(
      const std::string& name,
      const FieldSpace space,
      const FieldApproximationBase base,
      const FieldCoefficientsNumber nb_of_cooficients,
      const TagType tag_type = MB_TAG_SPARSE,
      const enum MoFEMTypes bh = MF_EXCL,
      int verb = -1
    );

    /**
     * \brief Setup problem
     * @param  meshset If only problem is set-up on part of the mesh
     * @return         Error code
     */
    PetscErrorCode setUp(
      std::string names_of_problem = "SimpleProblem",
      const EntityHandle meshset = 0,
      bool add_skeleton = false
    );

    /**
     * \brief add element to assemble domain
     * @param  fe element pointer
     * @return    errr code
     */
    PetscErrorCode domainAssemble(boost::shared_ptr<FEMethod> &fe);

    /**
     * \brief add element to assemble boundary
     * @param  fe element pointer
     * @return    errr code
     */
    PetscErrorCode boundaryAssemble(boost::shared_ptr<FEMethod> &fe);


    /**
     * \brief add element to assemble skeleton
     * @param  fe element pointer
     * @return    errr code
     */
    PetscErrorCode skeletonAssemble(boost::shared_ptr<FEMethod> &fe);

    /**
     * \brief Solve linear problem
     * @return Error code
     */
    PetscErrorCode kspSolve();

    /**
     * \brief ppstProcProblem
     * @return [description]
     */
    PetscErrorCode postProc();

  private:

    PetscLogEvent USER_EVENT_Simple   ;
    EntityHandle meshSet;                     ///< domain meshset
    std::string nameOfProblem;                ///< problem name
    std::vector<std::string> domainFields;    ///< domain fields
    std::vector<std::string> boundaryFields;  ///< boundary fields

    std::vector<boost::shared_ptr<FEMethod> > feDomain;
    std::vector<boost::shared_ptr<FEMethod> > feBoundary;
    std::vector<boost::shared_ptr<FEMethod> > feSkeleton;

  };
}

#endif // __SIMPLE_HPP__
