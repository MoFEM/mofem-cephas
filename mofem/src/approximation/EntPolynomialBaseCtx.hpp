/** \file EntPolynomialBaseCtx.hpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

*/



#ifndef __ENTPOLYNOMIALCTX_HPP__
#define __ENTPOLYNOMIALCTX_HPP__

namespace MoFEM {

struct EntitiesFieldData;
struct FEMethod;

/**
 * \brief Class used to pass element data to calculate base functions on
 * tet,triangle,edge
 *
 * \ingroup mofem_base_functions
 */
struct EntPolynomialBaseCtx : public BaseFunctionCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  PetscErrorCode (*basePolynomialsType0)(int p, double s, double *diff_s,
                                         double *L, double *diffL,
                                         const int dim);

  PetscErrorCode (*basePolynomialsType1)(int p, double alpha, double x,
                                         double t, double *diff_x,
                                         double *diff_t, double *L,
                                         double *diffL, const int dim);

  EntitiesFieldData &dAta;
  const FieldSpace sPace;
  const FieldApproximationBase bAse;
  const std::string fieldName;
  const FieldApproximationBase copyNodeBase;

  EntPolynomialBaseCtx(EntitiesFieldData &data, const FieldSpace space,
                       const FieldApproximationBase base,
                       const FieldApproximationBase copy_node_base = LASTBASE);

  EntPolynomialBaseCtx(EntitiesFieldData &data,
                       const std::string field_name, const FieldSpace space,
                       const FieldApproximationBase base,
                       const FieldApproximationBase copy_node_base = LASTBASE);

protected:
  MoFEMErrorCode setBase();
};

} // namespace MoFEM

#endif //__ENTPOLYNOMIALCTX_HPP__
