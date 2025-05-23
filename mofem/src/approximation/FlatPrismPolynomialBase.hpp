/** \file FlatPrismPolynomialBase.hpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

*/



#ifndef __FLATPRISMPOLYNOMIALBASE_HPP__
#define __FLATPRISMPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Class used to pass element data to calculate base functions on flat
 * prism
 *
 * \ingroup mofem_base_functions
 * FIXME: Need moab and mofem finite element structure to work (that not
 * perfect)
 */
struct FlatPrismPolynomialBaseCtx : public EntPolynomialBaseCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  moab::Interface &mOab;
  const NumeredEntFiniteElement *fePtr;

  FlatPrismPolynomialBaseCtx(
      EntitiesFieldData &data, moab::Interface &moab,
      const NumeredEntFiniteElement *fe_ptr, const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE);

  ~FlatPrismPolynomialBaseCtx();
};

/**
 * \brief Calculate base functions on tetrahedral
 * \ingroup mofem_base_functions
 * FIXME: Need moab and mofem finite element structure to work (that not
 * perfect)
 */
struct FlatPrismPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  FlatPrismPolynomialBase();
  ~FlatPrismPolynomialBase();

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  FlatPrismPolynomialBaseCtx *cTx;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  int numNodes;
  const EntityHandle *connPrism;
  const EntityHandle *connFace3;
  const EntityHandle *connFace4;
  int faceNodes[2][3];
  MatrixDouble N;
  MatrixDouble diffN;
};

} // namespace MoFEM

#endif //__FLATPRISMPOLYNOMIALBASE_HPP__
