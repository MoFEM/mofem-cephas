/** \file ForcesAndSourcesCore.hpp
  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

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

using namespace boost::numeric;

#ifndef __CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__
#define __CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__

namespace MoFEM {

/** \brief ContactPrism finite element
 \ingroup mofem_forces_and_sources_prism_element

 User is implementing own operator at Gauss points level, by own object
 derived from ContactPrismElementForcesAndSourcesCoreL::UserDataOperator. Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct ContactPrismElementForcesAndSourcesCore : public ForcesAndSourcesCore {

  double aRea[2];
  VectorDouble normal;
  VectorDouble coords;
  MatrixDouble coordsAtGaussPtsMaster;
  MatrixDouble coordsAtGaussPtsSlave;

  MatrixDouble gaussPtsMaster;
  MatrixDouble gaussPtsSlave;

  // DataForcesAndSourcesCore &dataH1Master;
  // DataForcesAndSourcesCore &dataH1Slave;

  std::string meshPositionsFieldName;


   /**
     * @brief Entity data on element entity rows fields
     *
     *
     * FIXME: that should be moved to private class data and acessed only by
     * member function
     */
    const boost::shared_ptr<DataForcesAndSourcesCore> dataOnMaster[LASTSPACE];
    const boost::shared_ptr<DataForcesAndSourcesCore> dataOnSlave[LASTSPACE];
  
    /**
     * @brief Entity data on element entity columns fields
     *
     * FIXME: that should be moved to private class data and acessed only by
     * member function
     */
    const boost::shared_ptr<DataForcesAndSourcesCore>
        derivedDataOnMaster[LASTSPACE];
    const boost::shared_ptr<DataForcesAndSourcesCore>
        derivedDataOnSlave[LASTSPACE];

    DataForcesAndSourcesCore &dataH1Master;
    DataForcesAndSourcesCore &dataH1Slave;
    // boost::shared_ptr<DataForcesAndSourcesCore> dataH1Master;
    // boost::shared_ptr<DataForcesAndSourcesCore> dataH1Slave;

    DataForcesAndSourcesCore &dataNoFieldMaster;
    DataForcesAndSourcesCore &dataNoFieldSlave;
    DataForcesAndSourcesCore &dataHcurlMaster;
    DataForcesAndSourcesCore &dataHcurlSlave;
    DataForcesAndSourcesCore &dataHdivMaster;
    DataForcesAndSourcesCore &dataHdivSlave;
    DataForcesAndSourcesCore &dataL2Master;
    DataForcesAndSourcesCore &dataL2Slave;


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

    /** \brief get face aRea Master
    */
    
    inline double getAreaMaster() {
      return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
          ->aRea[0];
    }

    /** \brief get face aRea Slave
    */
    
    inline double getAreaSlave() {
      return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
          ->aRea[1];
    }


    inline VectorAdaptor getNormalMaster() {
      double *data =
          &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                ->normal[0]);
      return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
    }

    inline VectorAdaptor getNormalSlave() {
      double *data =
          &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
                ->normal[3]);
      return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
    }

    /** \brief get triangle coordinates

      Vector has 9 elements, i.e. coordinates on Master face

     */
    inline VectorDouble getCoordsMaster() {
     double *data =
          &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->coords[0]);
      return VectorAdaptor(9, ublas::shallow_array_adaptor<double>(9, data));
    }

    /** \brief get triangle coordinates

      Vector has 9 elements, i.e. coordinates on Slave face

     */

    inline VectorDouble getCoordsSlave() {
     double *data =
          &(static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)->coords[9]);
      return VectorAdaptor(9, ublas::shallow_array_adaptor<double>(9, data));
    }

    /** \brief get coordinates at Gauss pts on full prism.

      Matrix has size (nb integration points on master)x(3),
      i.e. coordinates on face Master

     */
    inline MatrixDouble &getCoordsAtGaussPtsMaster() {
      return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
          ->coordsAtGaussPtsMaster;
    }

    /** \brief get coordinates at Gauss pts on full prism.

      Matrix has size (nb integration points on master)x(3),
      i.e. coordinates on face Master

     */
    inline MatrixDouble &getCoordsAtGaussPtsSlave() {
      return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE)
          ->coordsAtGaussPtsSlave;
    }


     /** \brief return pointer to triangle finite element object
     */
    inline const ContactPrismElementForcesAndSourcesCore *
    getContactPrismElementForcesAndSourcesCore() {
      return static_cast<ContactPrismElementForcesAndSourcesCore *>(ptrFE);
    }
  };

  MoFEMErrorCode operator()();
};

} // namespace MoFEM

#endif //__CONTACTPRISMELEMENTFORCESANDSURCESCORE_HPP__

/***************************************************************************/ /**
* \defgroup mofem_forces_and_sources_prism_element Prism Element
* \ingroup mofem_forces_and_sources
******************************************************************************/
