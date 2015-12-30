/** \file definitions.h
 * \brief useful compiler directives and definitions
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

#ifndef __DEFINITONS_H__
#define __DEFINITONS_H__

/** \brief Interfaces IDs
  *
  * To manage different complexities related to field, finite elements mesh
  * refinements, etc. a appropriate interfaces related to each complexities are
  * created. Interfaces by itself could vary by functionality or the same function
  * can me managed with two interfaces with waring level of abstraction.
  *
  */
enum MoFEMInterfaces {
  FIELD_UNKNOWNINTERFACE = 1<<0,
  //Field Interface
  FIELD_INTERFACE = 1<<0|1<<1,
  MESH_REFINE = 1<<1|1<<2,
  PRISM_INTEFACE = 1<<1|1<<3,
  SERIES_RECORDER = 1<<1|1<<4,
  //Loop Methods
  BASIC_METHOD = 1<<2|1<<3|1<<4,
  FE_METHOD = 1<<2|1<<3|1<<4|1<<5,
  ENT_METHOD = 1<<2|1<<3|1<<4|1<<6,
  //Independet Interfaces
  TETGEN_INTERFACE = 1<<3|1<<4,		///< used to generate mesh using TetGen
  NETGEN_INTERFACE = 1<<3|1<<5,		///< used to generate mesh using NetGen
  NODEMERGER_INTERFACE = 1<<3|1<<6,	///< used to merge nodes
  BITLEVELCOUPLER_INTERFACE = 1<<3|1<<7, ///< used to couple bit levels by finding parent children relation
  PRISMSFROMSURFACE_INTERFACE = 1<<3|1<<8 ///< create prisms from surface elements
};

/** \brief Error handling
  *
  * This is complementary to PETSC error codes. The numerical values for
  * these are defined in include/petscerror.h. The names are defined in err.c
  *
  * MoAB error messages are defined in naob/Types.hpp
  *
  */
enum MoFEMErrorCode {
  MOFEM_SUCESS = 0,
  MOFEM_DATA_INCONSISTENCY = 100,
  MOFEM_NOT_IMPLEMENTED = 101,
  MOFEM_NOT_FOUND = 102,
  MOFEM_OPERATION_UNSUCCESSFUL = 103,
  MOFEM_IMPOSIBLE_CASE = 104,
  MOFEM_MOFEMEXCEPTION_THROW = 105,
  MOFEM_STD_EXCEPTION_THROW = 106,
  MOFEM_INVALID_DATA = 107,
  MOFEM_ATOM_TEST_INVALID = 108,
  MOFEM_MOAB_ERROR = 109
};

/// \brief approximation spaces
enum FieldSpace {
  NOFIELD = 1, 	///< scalar or vector of scalars describe (no true field)
  H1, 		///< continuous field
  HDIV,		///< field with continuous normal traction
  HCURL,	///< field with continuous tangents
  L2,		///< field with C-1 continuity
  LASTSPACE 	///< FieldSpace in [ 0, LASTSPACE )
};

/// \brief Those types control how functions respond on arguments, f.e. error handling
enum MoFEMTypes {
  MF_ZERO = 0,
  MF_EXCL = 1<<0
};

/// \brief RowColData
enum RowColData {
  ROW,COL,DATA,LASTROWCOLDATA
};

enum ByWhat {
  BYROW = 1<<0, BYCOL = 1<<1, BYDATA = 1<<2,
  BYROWDATA = 1<<0|1<<2, BYCOLDATA = 1<<1|1<<2, BYROWCOL = 1<<0|1<<1,
  BYALL = 1<<0|1<<1|1<<2
};

/**
  * Types of sets and boundary conditions
  *
  */
enum CubitBC {
  UNKNOWNSET = 0,
  NODESET = 1<<0,
  SIDESET = 1<<1,
  BLOCKSET = 1<<2,
  MATERIALSET = 1<<3,
  DISPLACEMENTSET = 1<<4,
  FORCESET = 1<<5,
  PRESSURESET = 1<<6,
  VELOCITYSET = 1<<7,
  ACCELERATIONSET = 1<<8,
  TEMPERATURESET = 1<<9,
  HEATFLUXSET = 1<<10,
  INTERFACESET = 1<<11,
  UNKNOWNCUBITNAME = 1<< 12,
  MAT_ELASTICSET = 1<<13,	///< block name is "MAT_ELASTIC"
  MAT_INTERFSET = 1 <<14,
  MAT_THERMALSET = 1<<15,	///< block name is "MAT_THERMAL"
  BODYFORCESSET = 1<<16,	///< block name is "BODY_FORCES"
  MAT_MOISTURESET = 1<<17, 	///< block name is "MAT_MOISTURE"
  LASTCUBITSET
};

//taken from http://stackoverflow.com/questions/295120/c-mark-as-deprecated
#ifdef __GNUC__
  #define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
  #define DEPRECATED __declspec(deprecated)
#else
  #pragma message("WARNING: You need to implement DEPRECATED for this compiler")
  #define DEPRECATED
#endif

#define BITREFEDGES_SIZE 6 /*number of edges on tetrahedral*/
#define BITREFLEVEL_SIZE 128 /*max number of refinements*/
#define BITFIELDID_SIZE 32 /*max number of fields*/
#define BITFEID_SIZE 32 /*max number of finite elements*/
#define BITPROBLEMID_SIZE 32 /*max number of problems*/
#define BITINTERFACEUID_SIZE 32
#define BITCOORDSYS_SIZE 8

//// default comunicator number
#define MYPCOMM_INDEX 0

//This Is form MOAB
#define MB_TYPE_WIDTH 4
#define MB_ID_WIDTH (8*sizeof(EntityHandle)-MB_TYPE_WIDTH)
#define MB_TYPE_MASK ((EntityHandle)0xF << MB_ID_WIDTH)
//             2^MB_TYPE_WIDTH-1 ------^

#define MB_START_ID ((EntityID)1)        //!< All entity id's currently start at 1
#define MB_END_ID ((EntityID)MB_ID_MASK) //!< Last id is the complement of the MASK
#define MB_ID_MASK (~MB_TYPE_MASK)

#define NOT_USED(x) ( (void)(x) )

/** \brief set barrier start
 *
 * Run code in sequence, starting from process 0, and ends on last process.
 */
#define BARRIER_RANK_START(PCMB) \
  { for(unsigned int i = 0; \
  i<PCMB->proc_config().proc_rank(); i++) MPI_Barrier(PCMB->proc_config().proc_comm()); };
/// set barrier end
#define BARRIER_RANK_END(PCMB) \
  { for(unsigned int i = PCMB->proc_config().proc_rank(); \
  i<PCMB->proc_config().proc_size(); i++) MPI_Barrier(PCMB->proc_config().proc_comm()); };


//ERROR
/// check moab error
#define CHKERR(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    cerr << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    assert(1); \
  } \
} while (false)


/// check moab error and communicate it using petsc interface
#define CHKERR_PETSC(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::ostringstream ss; \
    ss << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    std::string str(ss.str()); \
    SETERRQ(PETSC_COMM_SELF,MOFEM_MOAB_ERROR,str.c_str()); \
  } \
} while (false)

#define CHKERR_THROW(a) do { \
  ErrorCode val = (a); \
  if (MB_SUCCESS != val) { \
    std::ostringstream ss; \
    ss << "Error code  " << val << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
    throw MoFEMException(MOFEM_MOAB_ERROR,ss.str().c_str() ); \
  } \
} while (false)

#define THROW_AT_LINE(a) { \
  std::ostringstream ss; \
  ss << a << " " << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
  throw MoFEMException( MOFEM_MOFEMEXCEPTION_THROW,ss.str().c_str() ); \
}

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
  ( std::ostringstream() << std::dec << x ) ).str()

#endif //__DEFINITONS_H__
