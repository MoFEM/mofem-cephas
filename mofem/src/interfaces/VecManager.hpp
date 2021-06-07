/** \file VecManager.hpp
 * \brief Interface managing vectors
 * \ingroup mofem_vectors
 *
 * Managing problems, build and partitioning.
 *
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

#ifndef __VECMANAGER_HPP__
#define __VECMANAGER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMVEC =
    MOFEMuuid(BitIntefaceId(VECMANAGER_INTERFACE));

/**
 * \brief Vector manager is used to create vectors
 * \mofem_vectors
 *
 * Managing Vectors, creation, scatter, etc.
 *
 */
struct VecManager : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  const MoFEM::Interface &cOre;
  bool dEbug;

  VecManager(const MoFEM::Core &core);

  /**
   * \brief Destructor
   */
  ~VecManager();

  /** \brief create local vector for problem
   * \ingroup mofem_vectors
   *
   * \param name problem name
   * \param RowColData specify what data is taken from Row, Col or Data
   * \param Vec the vector where data is stored
   */
  MoFEMErrorCode vecCreateSeq(const std::string name, RowColData rc,
                              Vec *V) const;

  /**
   * \brief create ghost vector for problem (collective)
   * \ingroup mofem_vectors

  collective - need to be run on all processors in communicator

   * \param name problem name
   * \param RowColData specify what data is taken from Row, Col or Data
   * \param Vec the vector where data is stored
   */
  MoFEMErrorCode vecCreateGhost(const std::string name, RowColData rc,
                                Vec *V) const;

  /** @copydoc MoFEM::VecManager::vecCreateGhost
   */
  MoFEMErrorCode vecCreateGhost(const std::string name, RowColData rc,
                                SmartPetscObj<Vec> &v_ptr) const;

  /**
    * \brief create scatter for vectors form one to another problem (collective)
    * \ingroup mofem_vectors
    *
    * User specify what name of field on one problem is scattered to another.
    *
    * \ingroup mofem_vectors
    *
    * \param xin vector
    * \param x_proble problem name
    * \param x_field name
    * \param yin vector
    * \param y_problem problem name
    * \param y_field_name
    * \retval newctx scatter

    */
  MoFEMErrorCode vecScatterCreate(Vec xin, const std::string x_problem,
                                  const std::string x_field_name,
                                  RowColData x_rc, Vec yin,
                                  const std::string y_problem,
                                  const std::string y_field_name,
                                  RowColData y_rc, VecScatter *newctx) const;

  /** @copydoc MoFEM::VecManager::vecScatterCreate
   */
  MoFEMErrorCode
  vecScatterCreate(Vec xin, const std::string x_problem,
                   const std::string x_field_name, RowColData x_rc, Vec yin,
                   const std::string y_problem, const std::string y_field_name,
                   RowColData y_rc,
                   SmartPetscObj<VecScatter> &smart_newctx) const;

  /**
    * \brief create scatter for vectors from one to another problem (collective)
    * \ingroup mofem_vectors
    *
    * \param xin vector
    * \param x_proble problem name
    * \param yin vector
    * \param y_problem problem name
    * \retval newctx scatter

    */
  MoFEMErrorCode vecScatterCreate(Vec xin, const std::string x_problem,
                                  RowColData x_rc, Vec yin,
                                  const std::string y_problem, RowColData y_rc,
                                  VecScatter *newctx) const;

  /** @copydoc MoFEM::VecManager::vecScatterCreate
   */
  MoFEMErrorCode
  vecScatterCreate(Vec xin, const std::string x_problem, RowColData x_rc,
                   Vec yin, const std::string y_problem, RowColData y_rc,
                   SmartPetscObj<VecScatter> &smart_newctx) const;

  /**
   * \brief set values of vector from/to meshdatabase
   * \ingroup mofem_vectors
   *
   * \param pointer to problem struture
   * \param RowColData for row or column:e (i.e. Row,Col)
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities from V vector.
   *
   * SCATTER_FORWARD set vector V from data field entities
   *
   */
  MoFEMErrorCode setLocalGhostVector(const Problem *problem_ptr, RowColData rc,
                                     Vec V, InsertMode mode,
                                     ScatterMode scatter_mode) const;

  /**
   * \brief set values of vector from/to meshdatabase
   * \ingroup mofem_vectors
   *
   * \param name of the problem
   * \param RowColData for row or column:e (i.e. Row,Col)
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities from V vector.
   *
   * SCATTER_FORWARD set vector V from data field entities
   *
   */
  MoFEMErrorCode setLocalGhostVector(const std::string name, RowColData rc,
                                     Vec V, InsertMode mode,
                                     ScatterMode scatter_mode) const;

  /**
    * \brief set values of vector from/to mesh database (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param pointer to porblem struture
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes
    are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  MoFEMErrorCode setGlobalGhostVector(const Problem *problem_ptr, RowColData rc,
                                      Vec V, InsertMode mode,
                                      ScatterMode scatter_mode) const;

  /**
    * \brief set values of vector from/to mesh database (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param name of the problem
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes
    are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  MoFEMErrorCode setGlobalGhostVector(const std::string name, RowColData rc,
                                      Vec V, InsertMode mode,
                                      ScatterMode scatter_mode) const;

  /**
   * \brief Copy vector to field which is not part of the problem
   * \ingroup mofem_vectors
   *
   * \param pointer to problem multi_index
   * \param field_name field name used for indexing petsc vectors used in the
   * problem \param cpy_field field name where data from vector are stored
   * \param RowColData for row or column
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities form V vector.
   *
   */
  MoFEMErrorCode setOtherLocalGhostVector(const Problem *problem_ptr,
                                          const std::string field_name,
                                          const std::string cpy_field_name,
                                          RowColData rc, Vec V, InsertMode mode,
                                          ScatterMode scatter_mode) const;

  /**
   * \brief Copy vector to field which is not part of the problem
   * \ingroup mofem_vectors
   *
   * \param name problem name
   * \param field_name field name used for indexing petsc vectors used in the
   * problem \param cpy_field field name where data from vector are stored
   * \param RowColData for row or column
   * \param V vector
   * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
   * \param scatter_mode see petsc manual for ScatterMode (The available modes
   * are: SCATTER_FORWARD or SCATTER_REVERSE)
   *
   * SCATTER_REVERSE set data to field entities form V vector.
   *
   */
  MoFEMErrorCode setOtherLocalGhostVector(const std::string name,
                                          const std::string field_name,
                                          const std::string cpy_field_name,
                                          RowColData rc, Vec V, InsertMode mode,
                                          ScatterMode scatter_mode) const;

  /** \brief Copy vector to field which is not part of the problem (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param problem_ptr pointer to problem
    * \param field_name field name used for indexing petsc vectors used in the
    problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes
    are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  MoFEMErrorCode setOtherGlobalGhostVector(const Problem *problem_ptr,
                                           const std::string field_name,
                                           const std::string cpy_field_name,
                                           RowColData rc, Vec V,
                                           InsertMode mode,
                                           ScatterMode scatter_mode) const;

  /** \deprecated VecManager
    * \brief Copy vector to field which is not part of the problem (collective)
    * \ingroup mofem_vectors

    collective - need tu be run on all processors in communicator

    * \param name problem name
    * \param field_name field name used for indexing petsc vectors used in the
    problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes
    are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  MoFEMErrorCode setOtherGlobalGhostVector(const std::string name,
                                           const std::string field_name,
                                           const std::string cpy_field_name,
                                           RowColData rc, Vec V,
                                           InsertMode mode,
                                           ScatterMode scatter_mode) const;
};

} // namespace MoFEM

#endif //__VECMANAGER_HPP__

/**
 * \defgroup mofem_vectors Vectors (Vec)
 * \brief Creating and scattering vectors on the mesh for given problem
 *
 * \ingroup mofem
 **/
