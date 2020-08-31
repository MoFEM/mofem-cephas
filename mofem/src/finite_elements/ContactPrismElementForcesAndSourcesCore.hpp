/** \file ContactPrismElementForcesAndSourcesCore.hpp
  \brief Implementation of the contact prism element.

  These elements are used to enforce contact constraints in the interface
  between two solids.

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

#ifndef __CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__
#define __CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__

namespace MoFEM {
template <int SWITCH>
struct VolumeElementForcesAndSourcesCoreOnContactPrismSideSwitch;

/** \brief ContactPrism finite element
 \ingroup mofem_forces_and_sources_prism_element

 User is implementing own operator at Gauss points level, by own class
 derived from ContactPrismElementForcesAndSourcesCoreL::UserDataOperator.
 Arbitrary number of operator added pushing instances to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct ContactPrismElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  ContactPrismElementForcesAndSourcesCore(Interface &m_field);

  /** \brief default operator for Contact Prism element
   * \ingroup mofem_forces_and_sources_prism_element
   */
  struct UserDataOperator : public ForcesAndSourcesCore::UserDataOperator {

    UserDataOperator(const FieldSpace space)
        : ForcesAndSourcesCore::UserDataOperator(space) {}

    UserDataOperator(const std::string &field_name, const char type)
        : ForcesAndSourcesCore::UserDataOperator(field_name, type) {}

    UserDataOperator(const std::string &row_field_name,
                     const std::string &col_field_name, const char type)
        : ForcesAndSourcesCore::UserDataOperator(row_field_name, col_field_name,
                                                 type) {}

    UserDataOperator(const std::string &row_field_name,
                     const std::string &col_field_name, const char type,
                     const char face_type)
        : ForcesAndSourcesCore::UserDataOperator(row_field_name, col_field_name,
                                                 type),
          faceType(face_type) {}

    UserDataOperator(const std::string &field_name, const char type,
                     const char face_type)
        : ForcesAndSourcesCore::UserDataOperator(field_name, type),
          faceType(face_type) {}

    /**
           * \brief Controls loop over faces and face combination on element
            *
            * FACEMASTER is used when column or row data needs to be accessed
      located at master face
            * FACESLAVE is used when column or row data needs to be accessed
      located at slave face
            * FACEMASTERMASTER is used for accessing simultaneously row and col
      data located at master face.
            * FACEMASTERSLAVE is used for accessing simultaneously row data that
      is located on master face and col data located at slave face.
            * FACESLAVEMASTER is used for accessing simultaneously row data that
      is located on slave face and col data located at master face.
      * FACESLAVESLAVE is used for accessing simultaneously row and col
      data located at slave face.
            *
      */
    enum FaceType {
      FACEMASTER = 1 << 0,
      FACESLAVE = 1 << 1,
      FACEMASTERMASTER = 1 << 2,
      FACEMASTERSLAVE = 1 << 3,
      FACESLAVEMASTER = 1 << 4,
      FACESLAVESLAVE = 1 << 5,
      FACELAST = 1 << 6
    };

    char faceType;

    /**
     * \brief Get operator types
     * @return Return operator type
     */
    inline int getFaceType() const;

    inline boost::shared_ptr<const NumeredEntFiniteElement>
    getNumeredEntFiniteElementPtr() const;

    /** \brief get face aRea Master
     */
    inline double getAreaMaster();

    /** \brief get face aRea Slave
     */
    inline double getAreaSlave();

    /** \brief get face normal vector to Master face
     */
    inline VectorAdaptor getNormalMaster();

    /** \brief get first face tangent vector to Master face
     */
    inline VectorAdaptor getTangentMasterOne();

    /** \brief get second face tangent vector to Master face
     */
    inline VectorAdaptor getTangentMasterTwo();

    /** \brief get face normal vector to Slave face
     */
    inline VectorAdaptor getNormalSlave();

    /** \brief get first face tangent vector to Slave face
     */
    inline VectorAdaptor getTangentSlaveOne();

    /** \brief get second face tangent vector to Slave face
     */
    inline VectorAdaptor getTangentSlaveTwo();

    /** \brief get Gauss point at Master face
     */
    inline MatrixDouble &getGaussPtsMaster();

    /** \brief get Gauss point at Slave face
     */
    inline MatrixDouble &getGaussPtsSlave();

    /**
     * @brief Get integration weights for slave side
     *
     * \code
     * auto t_w = getFTensor0IntegrationWeight();
     * for(int gg = 0; gg!=getGaussPts.size2(); ++gg) {
     *  // integrate something
     *  ++t_w;
     * }
     * \endcode
     *
     * @return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
     */
    inline auto getFTensor0IntegrationWeightSlave();

    /**
     * @brief Get integration weights for master side
     *
     * \code
     * auto t_w = getFTensor0IntegrationWeight();
     * for(int gg = 0; gg!=getGaussPts.size2(); ++gg) {
     *  // integrate something
     *  ++t_w;
     * }
     * \endcode
     *
     * @return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
     */
    inline auto getFTensor0IntegrationWeightMaster();

    /** \brief get triangle coordinates

      Vector has 9 elements, i.e. coordinates on Master face

     */
    inline VectorDouble getCoordsMaster();

    /** \brief get triangle coordinates

      Vector has 9 elements, i.e. coordinates on Slave face

     */
    inline VectorDouble getCoordsSlave();

    /** \brief get coordinates at Gauss pts on full prism.

      Matrix has size (nb integration points on master)x(3),
      i.e. coordinates on face Master

     */
    inline MatrixDouble &getCoordsAtGaussPtsMaster();

    /** \brief get coordinates at Gauss pts on full prism.

      Matrix has size (nb integration points on slave)x(3),
      i.e. coordinates on face Slave

     */
    inline MatrixDouble &getCoordsAtGaussPtsSlave();

    /** \brief return pointer to triangle finite element object
     */
    inline const ContactPrismElementForcesAndSourcesCore *
    getContactPrismElementForcesAndSourcesCore();

    /**
     *
     * User call this function to loop over elements on the side of face. This
     * function calls MoFEM::VolumeElementForcesAndSourcesCoreOnContactPrismSide with
     * is operator to do calculations.
     *
     * @param  fe_name Name of the element
     * @param  method  Finite element object
     * @param  side_type  states the side from which side element will work (0
     * for master 1 for slave)
     * @return         error code
     */
    template <int SWITCH>
    MoFEMErrorCode loopSideVolumes(
        const string &fe_name,
        VolumeElementForcesAndSourcesCoreOnContactPrismSideSwitch<SWITCH> &fe_method,
        const int side_type, const EntityHandle ent_for_side);

  protected:
    inline ForcesAndSourcesCore *getSidePtrFE() const;
  };

  MoFEMErrorCode operator()();

  inline const std::array<boost::shared_ptr<DataForcesAndSourcesCore>,
                          LASTSPACE>
  getDataOnMasterFromEleSide() {
    return dataOnMaster;
  }

  inline const std::array<boost::shared_ptr<DataForcesAndSourcesCore>,
                          LASTSPACE>
  getDataOnSlaveFromEleSide() {
    return dataOnSlave;
  }

  inline MatrixDouble &getGaussPtsMasterFromEleSide() { return gaussPtsMaster; }

  inline MatrixDouble &getGaussPtsSlaveFromEleSide() { return gaussPtsSlave; }

protected:
  std::array<double, 2> aRea; ///< Array storing master and slave faces areas

  VectorDouble
      normal; ///< vector storing vector normal to master or slave element
  VectorDouble coords;
  MatrixDouble coordsAtGaussPtsMaster; ///< matrix storing master Gauss points
                                       ///< global coordinates
  MatrixDouble coordsAtGaussPtsSlave;  ///< matrix storing slave Gauss points
                                       ///< global coordinates

  MatrixDouble gaussPtsMaster; ///< matrix storing master Gauss points local
                               ///< coordinates and weights
  MatrixDouble gaussPtsSlave;  ///< matrix storing slave Gauss points local
                               ///< coordinates and weights

  VectorDouble nOrmalSlave, tangentSlaveOne, tangentSlaveTwo;
  VectorDouble nOrmalMaster, tangentMasterOne, tangentMasterTwo;

  MatrixDouble normalsAtGaussPtsSlave;
  MatrixDouble tangentOneAtGaussPtsSlave;
  MatrixDouble tangentTwoAtGaussPtsSlave;
  
  OpSetContravariantPiolaTransformOnFace opContravariantTransform;
  // OpSetInvJacHdivAndHcurl opInvJacHdivAndHcurl;

  /**
   * @brief Entity data on element entity rows fields
   *
   *
   * FIXME: that should be moved to private class data and acessed only by
   * member function
   */
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      dataOnMaster;
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      dataOnSlave;

  /**
   * @brief Entity data on element entity columns fields
   *
   * FIXME: that should be moved to private class data and acessed only by
   * member function
   */
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      derivedDataOnMaster;
  const std::array<boost::shared_ptr<DataForcesAndSourcesCore>, LASTSPACE>
      derivedDataOnSlave;

  DataForcesAndSourcesCore &dataH1Master;
  DataForcesAndSourcesCore &dataH1Slave;

  DataForcesAndSourcesCore &dataNoFieldMaster;
  DataForcesAndSourcesCore &dataNoFieldSlave;
  DataForcesAndSourcesCore &dataHcurlMaster;
  DataForcesAndSourcesCore &dataHcurlSlave;
  DataForcesAndSourcesCore &dataHdivMaster;
  DataForcesAndSourcesCore &dataHdivSlave;
  DataForcesAndSourcesCore &dataL2Master;
  DataForcesAndSourcesCore &dataL2Slave;

  MoFEMErrorCode setDefaultGaussPts(const int rule);

  /**
   * @brief Iterate user data operators
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode loopOverOperators();

  /**
   * @brief Iterate user data operators
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts,
                                           FieldApproximationBase m_s_base,
                                           DataForcesAndSourcesCore &m_s_data);

  /** \brief function that gets entity field data.
   *
   * \param master_data data fot master face
   * \param slave_data data fot master face
   * \param field_name field name of interest
   * \param type_lo lowest dimension entity type to be searched
   * \param type_hi highest dimension entity type to be searched
   */
  MoFEMErrorCode
  getEntityFieldData(DataForcesAndSourcesCore &master_data,
                     DataForcesAndSourcesCore &slave_data,
                     const std::string &field_name,
                     const EntityType type_lo = MBVERTEX,
                     const EntityType type_hi = MBPOLYHEDRON) const;

  /** \brief function that gets entity indices.
   *
   * \param master_data data fot master face
   * \param slave_data data fot master face
   * \param field_name field name of interest
   * \param dofs MultiIndex container keeping FENumeredDofEntity.
   * \param type_lo lowest dimension entity type to be searched
   * \param type_hi highest dimension entity type to be searched
   */
  MoFEMErrorCode getEntityIndices(
      DataForcesAndSourcesCore &master_data,
      DataForcesAndSourcesCore &slave_data, const std::string &field_name,
      FENumeredDofEntity_multiIndex &dofs, const EntityType type_lo = MBVERTEX,
      const EntityType type_hi = MBPOLYHEDRON) const;

  /** \brief function that gets nodes indices.
   *
   * \param field_name field name of interest
   * \param dofs MultiIndex container keeping FENumeredDofEntity.
   * \param master_nodes_indices vector containing global master nodes indices
   * \param master_local_nodes_indices vector containing local master nodes
   * indices
   * \param slave_nodes_indices vector containing global master nodes indices
   * \param slave_local_nodes_indices vector containing local master nodes
   * indices
   */
  MoFEMErrorCode getNodesIndices(const boost::string_ref field_name,
                                 FENumeredDofEntity_multiIndex &dofs,
                                 VectorInt &master_nodes_indices,
                                 VectorInt &master_local_nodes_indices,
                                 VectorInt &slave_nodes_indices,
                                 VectorInt &slave_local_nodes_indices) const;

  /** \brief function that gets nodes field data.
   *
   * \param field_name field name of interest
   * \param dofs MultiIndex container keeping FENumeredDofEntity.
   * \param master_nodes_data vector containing master nodes data
   * \param slave_nodes_data vector containing master nodes data
   * \param master_nodes_dofs vector containing master nodes dofs
   * \param slave_nodes_dofs vector containing slave nodes dofs
   * \param master_space approximation energy space at master
   * \param slave_space approximation energy space at slave
   * \param master_base base for master face
   * \param slave_base base for slave face
   */
  MoFEMErrorCode getNodesFieldData(
      const boost::string_ref field_name, FEDofEntity_multiIndex &dofs,
      VectorDouble &master_nodes_data, VectorDouble &slave_nodes_data,
      VectorDofs &master_nodes_dofs, VectorDofs &slave_nodes_dofs,
      FieldSpace &master_space, FieldSpace &slave_space,
      FieldApproximationBase &master_base,
      FieldApproximationBase &slave_base) const;

private:
  int nbGaussPts;

};

boost::shared_ptr<const NumeredEntFiniteElement>
ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getNumeredEntFiniteElementPtr() const {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
      ->numeredEntFiniteElementPtr;
};

template <int SWITCH>
MoFEMErrorCode
ContactPrismElementForcesAndSourcesCore::UserDataOperator::loopSideVolumes(
    const string &fe_name,
    VolumeElementForcesAndSourcesCoreOnContactPrismSideSwitch<SWITCH> &fe_method,
    const int side_type, const EntityHandle ent_for_side) {
  return loopSide(fe_name, &fe_method, side_type, ent_for_side);
}

inline int
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getFaceType() const {
  return faceType;
}

inline double
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getAreaMaster() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->aRea[0];
}

inline double
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getAreaSlave() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->aRea[1];
}

inline VectorAdaptor
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getNormalMaster() {
  double *data = &(
      static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->normal[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangentMasterOne() {
  double *data = &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                       ->tangentMasterOne[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangentMasterTwo() {
  double *data = &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                       ->tangentMasterTwo[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getNormalSlave() {
  double *data = &(
      static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->normal[3]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangentSlaveOne() {
  double *data = &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                       ->tangentSlaveOne[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangentSlaveTwo() {
  double *data = &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                       ->tangentSlaveTwo[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline MatrixDouble &
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getGaussPtsMaster() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
      ->gaussPtsMaster;
}

inline MatrixDouble &
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getGaussPtsSlave() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
      ->gaussPtsSlave;
}

auto ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getFTensor0IntegrationWeightSlave() {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &(getGaussPtsSlave()(getGaussPtsSlave().size1() - 1, 0)));
}

auto ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getFTensor0IntegrationWeightMaster() {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &(getGaussPtsMaster()(getGaussPtsMaster().size1() - 1, 0)));
}

inline VectorDouble
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getCoordsMaster() {
  double *data = &(
      static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->coords[0]);
  return VectorAdaptor(9, ublas::shallow_array_adaptor<double>(9, data));
}

inline VectorDouble
ContactPrismElementForcesAndSourcesCore::UserDataOperator::getCoordsSlave() {
  double *data = &(
      static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->coords[9]);
  return VectorAdaptor(9, ublas::shallow_array_adaptor<double>(9, data));
}

inline MatrixDouble &ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getCoordsAtGaussPtsMaster() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
      ->coordsAtGaussPtsMaster;
}

inline MatrixDouble &ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getCoordsAtGaussPtsSlave() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
      ->coordsAtGaussPtsSlave;
}

inline const ContactPrismElementForcesAndSourcesCore *
ContactPrismElementForcesAndSourcesCore::UserDataOperator::
    getContactPrismElementForcesAndSourcesCore() {
  return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE);
}

} // namespace MoFEM

#endif //__CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__
