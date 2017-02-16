/** \file MedInterface.cpp
 * \brief Med file interface interface
 *
 * Interface loading mesh and data on mesh directly to mofem & moab
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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#ifdef WITH_MED

extern "C" {
#include <med.h>
}

#if (MED_MAJOR_NUM == 3)
  // To avoid too many ifdefs below we use defines for the bits of the
  // API that did not change too much between MED2 and MED3. If we
  // remove MED2 support at some point, please remove these defines and
  // replace the symbols accordingly.
  #define med_geometrie_element med_geometry_type
  #define med_entite_maillage med_entity_type
  #define med_type_champ med_field_type
  #define MED_LECTURE MED_ACC_RDONLY
  #define MED_LECTURE_AJOUT MED_ACC_RDEXT
  #define MED_NOEUD MED_NODE
  #define MED_MAILLE MED_CELL
  #define MED_NOEUD_MAILLE MED_NODE_ELEMENT
  #define MED_NOPFL MED_NO_PROFILE
  #define MEDouvrir MEDfileOpen
  #define MEDfermer MEDfileClose
  #define MEDnChamp MEDfieldnComponent
  #define MEDnValProfil MEDprofileSizeByName
#else
  #error "MED file is not MED_MAJOR_NUM == 3"
#endif

#include <MedInterface.hpp>

namespace MoFEM {

  PetscErrorCode MedInterface::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMMedInterface) {
      *iface = dynamic_cast<MedInterface*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  MedInterface::MedInterface(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  flgFile(PETSC_FALSE) {}

  PetscErrorCode MedInterface::getFileNameFromCommandLine(int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscErrorCode ierr;
    char mesh_file_name[255];
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(
      m_field.get_comm(),"","Shell prism configure","none"
    ); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-med_file",
      "med file name","", "mesh.med",mesh_file_name, 255, &flgFile
    ); CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    if(flgFile) {
      medFileName = std::string(mesh_file_name);
    } else {
      medFileName = std::string("mesh.bed");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::medGetFieldNames(const string &file,int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
    if(fid < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to open file '%s'",file.c_str());
    }
    med_int num_fields = MEDnField(fid);
    for(int index = 0; index < num_fields; index++){

      med_int num_comp = MEDnChamp(fid, index + 1);
      if(num_comp <= 0){
        SETERRQ(m_field.get_comm(),MOFEM_IMPOSIBLE_CASE,"Could not get number of components for MED field");
      }

      char name[MED_NAME_SIZE + 1], mesh_name[MED_NAME_SIZE + 1];
      char dt_unit[MED_SNAME_SIZE + 1];
      std::vector<char> comp_name(num_comp * MED_SNAME_SIZE + 1);
      std::vector<char> comp_unit(num_comp * MED_SNAME_SIZE + 1);
      med_int num_steps = 0;
      med_type_champ type;
      med_bool local_mesh;
      if(MEDfieldInfo(
        fid, index + 1, name, mesh_name, &local_mesh, &type,
        &comp_name[0], &comp_unit[0], dt_unit, &num_steps
      ) < 0) {
        SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not get MED field info");
      }

      fieldNames.push_back(name);
    }
    if(MEDfileClose(fid) < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to close file '%s'",file.c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::medGetFieldNames(int verb) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if(medFileName.empty()) {
      ierr = getFileNameFromCommandLine(verb); CHKERRQ(ierr);
    }
    ierr = medGetFieldNames(medFileName,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


  PetscErrorCode MedInterface::readMed(const string &file,int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;

    med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
    if(fid < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to open file '%s'",file.c_str());
    }

    med_int v[3], vf[3];
    MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
    MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);

    if(verb>0) {
      PetscPrintf(
        m_field.get_comm(),
        "Reading MED file V%d.%d.%d using MED library V%d.%d.%d\n",
        vf[0], vf[1], vf[2], v[0], v[1], v[2]
      );
    }
    if(vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)){
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Cannot read MED file older than V2.2");
    }

    for(int i = 0; i < MEDnMesh(fid); i++){
      char mesh_name[MED_NAME_SIZE + 1], mesh_desc[MED_COMMENT_SIZE + 1];
      med_int space_dim;
      med_mesh_type mesh_type;
      med_int mesh_dim, n_step;
      char dt_unit[ MED_SNAME_SIZE + 1];
      char axis_name[3 * MED_SNAME_SIZE + 1], axis_unit[3 * MED_SNAME_SIZE + 1];
      med_sorting_type sorting_type;
      med_axis_type axis_type;
      if(MEDmeshInfo(
        fid,i+1,mesh_name,&space_dim,&mesh_dim,&mesh_type,mesh_desc,
        dt_unit,&sorting_type,&n_step,&axis_type,axis_name,axis_unit
      ) < 0) {
        SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to read mesh information");
      }
      meshNames.push_back(mesh_name);
    }

    std::map<int,Range> family_elem_map;
    std::map<string,Range> group_elem_map;

    for(unsigned int ii = 0;ii!=meshNames.size();ii++) {
      ierr = readMesh(file,ii,family_elem_map,verb); CHKERRQ(ierr);
      ierr = readFamily(file,ii,family_elem_map,group_elem_map,verb); CHKERRQ(ierr);
      ierr = makeBlockSets(group_elem_map,verb); CHKERRQ(ierr);
    }

    if(MEDfileClose(fid) < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to close file '%s'",file.c_str());
    }

    PetscFunctionReturn(0);
  }

  static med_geometrie_element moab2med_element_type(const EntityType type) {
    switch(type) {
      case MBEDGE: return MED_SEG2;
      case MBTRI:  return MED_TRIA3;
      case MBQUAD: return MED_QUAD4;
      case MBTET: return MED_TETRA4;
      case MBHEX: return MED_HEXA8;
      case MBPRISM: return MED_PENTA6;
      case MBPYRAMID: return MED_PYRA5;
      case MBVERTEX: return MED_POINT1;
      default: return MED_NONE;
    }
  }

  PetscErrorCode MedInterface::readMesh(
    const string &file,
    const int index,
    std::map<int,Range> &family_elem_map,
    int verb
  ) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;

    med_idt fid = MEDfileOpen(file.c_str(), MED_ACC_RDONLY);
    if(fid < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to open file '%s'",file.c_str());
    }
    med_int v[3], vf[3];
    MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
    MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);
    // if(verb>1) {
    //   PetscPrintf(
    //     m_field.get_comm(),
    //     "Reading MED file V%d.%d.%d using MED library V%d.%d.%d\n",
    //     vf[0], vf[1], vf[2], v[0], v[1], v[2]
    //   );
    // }
    if(vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)){
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Cannot read MED file older than V2.2");
    }

    char mesh_name[MED_NAME_SIZE + 1], mesh_desc[MED_COMMENT_SIZE + 1];
    med_int space_dim;
    med_mesh_type mesh_type;
    med_int mesh_dim, n_step;
    char dt_unit[ MED_SNAME_SIZE + 1];
    char axis_name[3 * MED_SNAME_SIZE + 1], axis_unit[3 * MED_SNAME_SIZE + 1];
    med_sorting_type sorting_type;
    med_axis_type axis_type;
    if(MEDmeshInfo(
      fid,index+1,mesh_name,&space_dim,&mesh_dim,&mesh_type,mesh_desc,
      dt_unit,&sorting_type,&n_step,&axis_type,axis_name,axis_unit) < 0
    ) {
      SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to read mesh information");
    }
    if(verb>0) {
      PetscPrintf(m_field.get_comm(),"Reading mesh %s\n",mesh_name);
    }

    switch (axis_type) {
      case MED_CARTESIAN:
      break;
      case MED_SPHERICAL:
      case MED_CYLINDRICAL:
      default:
      SETERRQ(m_field.get_comm(),MOFEM_NOT_IMPLEMENTED,"Curvilinear coordinate system implemented");
    }
    if(space_dim!=3) {
      SETERRQ1(m_field.get_comm(),MOFEM_NOT_IMPLEMENTED,"Not implemented for space dim %d",space_dim);
    }

    EntityHandle mesh_meshset;
    {
      MeshsetsManager *meshsets_manager_ptr;
      ierr = m_field.query_interface(meshsets_manager_ptr); CHKERRQ(ierr);
      int max_id = 0;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {
        max_id = (max_id < cit->getMeshsetId()) ? cit->getMeshsetId() : max_id;
      }
      max_id++;
      ierr = meshsets_manager_ptr->addMeshset(BLOCKSET,max_id++,string(mesh_name)); CHKERRQ(ierr);
      CubitMeshsetByName::iterator cit =
      meshsets_manager_ptr->getMeshsetsMultindex().get<CubitMeshSets_name>().find(string(mesh_name));
      mesh_meshset = cit->getMeshset();
    }

    med_bool change_of_coord, geo_transform;
    med_int num_nodes = MEDmeshnEntity(
      fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_NODE,
      MED_NO_GEOTYPE, MED_COORDINATE, MED_NO_CMODE,
      &change_of_coord, &geo_transform
    );
    if(num_nodes < 0){
      SETERRQ(
        m_field.get_comm(),
        MOFEM_OPERATION_UNSUCCESSFUL,
        "Could not read number of MED nodes"
      );
    }
    if(num_nodes == 0) {
      SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"No nodes in MED mesh");
    }
    if(verb>1) {
      PetscPrintf(m_field.get_comm(),"Read number of node %d\n",num_nodes);
    }

    std::vector<med_float> coord_med(space_dim * num_nodes);
    if(MEDmeshNodeCoordinateRd(
      fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_NO_INTERLACE,&coord_med[0]
    ) < 0 ) {
      SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not read MED node coordinates");
    }

    // Add vertices to moab
    ReadUtilIface* iface;
    vector<double*> arrays_coord;
    EntityHandle startv;
    rval = m_field.get_moab().query_interface(iface); CHKERRQ_MOAB(rval);
    rval = iface->get_node_coords(3,num_nodes,0,startv,arrays_coord); CHKERRQ_MOAB(rval);
    Range verts(startv,startv+num_nodes-1);
    std::copy(&coord_med[0*num_nodes],&coord_med[1*num_nodes],arrays_coord[0]);
    std::copy(&coord_med[1*num_nodes],&coord_med[2*num_nodes],arrays_coord[1]);
    std::copy(&coord_med[2*num_nodes],&coord_med[3*num_nodes],arrays_coord[2]);
    ierr = m_field.get_moab().add_entities(mesh_meshset,verts); CHKERRQ(ierr);
    family_elem_map.clear();

    // get family for vertices
    {
      std::vector<med_int> fam(num_nodes,0);
      if(MEDmeshEntityFamilyNumberRd(
        fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_NODE,MED_NO_GEOTYPE,&fam[0]
      ) < 0) {
        SETERRQ(
          m_field.get_comm(),
          MOFEM_OPERATION_UNSUCCESSFUL,
          "No family number for elements: using 0 as default family number"
        );
      }
      for(int i = 0; i < num_nodes; i++) {
        // cerr << verts[i] << " " /*<< ele_tags[j] << " "*/ << fam[i] << endl;
        family_elem_map[fam[i]].insert(verts[i]);
      }
    }

    // read elements (loop over all possible MSH element types)
    for(EntityType ent_type = MBVERTEX; ent_type < MBMAXTYPE; ent_type++){

      med_geometrie_element type = moab2med_element_type(ent_type);
      if(type == MED_NONE) continue;

      // get number of cells
      med_bool change_of_coord;
      med_bool geo_transform;
      med_int num_ele = MEDmeshnEntity(
        fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_CELL,
        type, MED_CONNECTIVITY, MED_NODAL, &change_of_coord,
        &geo_transform
      );

      // get connectivity
      if(num_ele <= 0) continue;
      int num_nod_per_ele = type % 100;
      std::vector<med_int> conn_med(num_ele*num_nod_per_ele);
      if(MEDmeshElementConnectivityRd(
        fid, mesh_name, MED_NO_DT, MED_NO_IT, MED_CELL,
        type, MED_NODAL, MED_FULL_INTERLACE, &conn_med[0]
      ) < 0){
        SETERRQ(
          m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not read MED elements"
        );
      }

      // cerr << "type " << ent_type << " ";
      // cerr << "num_ele " << num_ele << " " << num_nod_per_ele << endl;;

      EntityHandle* conn_moab;
      EntityHandle starte;
      rval = iface->get_element_connect(
        num_ele,num_nod_per_ele,ent_type,0,starte,conn_moab
      ); CHKERRQ_MOAB(rval);
      switch (ent_type) {
        // FIXME: Some connectivity could not work, need to run and test
        default:
        {
          int ii = 0;
          for(int ee = 0;ee!=num_ele;ee++) {
            for(int nn = 0;nn!=num_nod_per_ele;nn++,ii++) {
              // cerr << conn_med[ii] << " ";
              conn_moab[ii] = verts[conn_med[ii]-1];
            }
            // cerr << endl;
          }
        }
      }
      rval = iface->update_adjacencies(
        starte,num_ele,num_nod_per_ele,conn_moab
      ); CHKERRQ_MOAB(rval);

      Range ents(starte,starte+num_ele-1);
      ierr = m_field.get_moab().add_entities(mesh_meshset,ents); CHKERRQ(ierr);

      // get family for cells
      {
        std::vector<med_int> fam(num_ele,0);
        if(MEDmeshEntityFamilyNumberRd(
          fid,mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,type,&fam[0]
        ) < 0) {
          SETERRQ(
            m_field.get_comm(),
            MOFEM_OPERATION_UNSUCCESSFUL,
            "No family number for elements: using 0 as default family number"
          );
        }
        // std::vector<med_int> ele_tags(num_ele);
        // if(MEDmeshEntityNumberRd(
        //   fid, mesh_name,MED_NO_DT,MED_NO_IT,MED_CELL,type,&ele_tags[0]) < 0
        // ) {
        //   ele_tags.clear();
        // }
        for(int j = 0; j < num_ele; j++) {
          // cerr << ents[j] << " " /*<< ele_tags[j] << " "*/ << fam[j] << endl;
          family_elem_map[fam[j]].insert(ents[j]);
        }
      }

    }

    if(MEDfileClose(fid) < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to close file '%s'",file.c_str());
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::readFamily(
    const string &file,
    const int index,
    const std::map<int,Range> &family_elem_map,
    std::map<string,Range> &group_elem_map,
    int verb
  ) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;

    med_idt fid = MEDfileOpen(file.c_str(),MED_ACC_RDONLY);
    if(fid < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to open file '%s'",file.c_str());
    }
    med_int v[3], vf[3];
    MEDlibraryNumVersion(&v[0], &v[1], &v[2]);
    MEDfileNumVersionRd(fid, &vf[0], &vf[1], &vf[2]);
    // if(verb>1) {
    //   PetscPrintf(
    //     m_field.get_comm(),
    //     "Reading MED file V%d.%d.%d using MED library V%d.%d.%d\n",
    //     vf[0], vf[1], vf[2], v[0], v[1], v[2]
    //   );
    // }
    if(vf[0] < 2 || (vf[0] == 2 && vf[1] < 2)){
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Cannot read MED file older than V2.2");
    }

    // clear groups
    group_elem_map.clear();

    med_int num_families = MEDnFamily(fid,meshNames[index].c_str());
    if(num_families < 0){
      SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not read MED families");
    }
    for(int i = 0;i<num_families;i++) {
      med_int num_attrib = (vf[0] == 2) ? MEDnFamily23Attribute(fid,meshNames[index].c_str(),i+1) : 0;
      med_int num_groups = MEDnFamilyGroup(fid,meshNames[index].c_str(),i+1);
      if(num_attrib < 0 || num_groups < 0){
        SETERRQ(
          m_field.get_comm(),
          MOFEM_OPERATION_UNSUCCESSFUL,
          "Could not read MED groups or attributes"
        );
      }

      std::vector<med_int> attribId(num_attrib + 1);
      std::vector<med_int> attrib_val(num_attrib + 1);
      std::vector<char> attrib_des(MED_COMMENT_SIZE * num_attrib + 1);
      std::vector<char> group_names(MED_LNAME_SIZE * num_groups + 1);
      char family_name[MED_NAME_SIZE + 1];
      med_int family_num;

      if(vf[0] == 2) { // MED2 file
        if(MEDfamily23Info(
          fid,meshNames[index].c_str(), i + 1, family_name, &attribId[0],
          &attrib_val[0], &attrib_des[0], &family_num,
          &group_names[0]) < 0
        ) {
          SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not read info for MED2 family %d",i+1);
        }
      }
      else{
        if(MEDfamilyInfo(
          fid,meshNames[index].c_str(),i+1,family_name,&family_num,
          &group_names[0]
        ) < 0) {
          SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Could not read info for MED3 family %d",i+1);
        }
      }

      // cerr << family_name << " " << family_num  << " " << num_groups << endl;
      for(int g = 0;g!=num_groups;g++) {
        group_elem_map[string(&group_names[MED_LNAME_SIZE*g])].merge(family_elem_map.at(family_num));
        // cerr << string(&group_names[MED_LNAME_SIZE*g]) << endl;
      }

    }

    if(MEDfileClose(fid) < 0) {
      SETERRQ1(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"Unable to close file '%s'",file.c_str());
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::makeBlockSets(
    const std::map<string,Range> &group_elem_map,
    int verb
  ) {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    MeshsetsManager *meshsets_manager_ptr;
    ierr = m_field.query_interface(meshsets_manager_ptr); CHKERRQ(ierr);

    int max_id = 0;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {
      max_id = (max_id < cit->getMeshsetId()) ? cit->getMeshsetId() : max_id;
    }
    max_id++;

    // cerr << group_elem_map.size() << endl;
    for(
      std::map<string,Range>::const_iterator git = group_elem_map.begin();
      git!=group_elem_map.end();git++
    ) {
      // cerr << "AAA\n";
      ierr = meshsets_manager_ptr->addMeshset(BLOCKSET,max_id++,git->first); CHKERRQ(ierr);
      CubitMeshsetByName::iterator cit =
      meshsets_manager_ptr->getMeshsetsMultindex().get<CubitMeshSets_name>().find(git->first);
      EntityHandle meshsets = cit->getMeshset();
      rval = m_field.get_moab().add_entities(meshsets,git->second); CHKERRQ_MOAB(rval);
      // cerr << git->second << endl;
    }

    // if(verb>0) {
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,cit)) {
        PetscPrintf(
          m_field.get_comm(),"%s\n",
          static_cast<std::ostringstream&>(
            std::ostringstream().seekp(0) << *cit
          ).str().c_str()
        );
      }
    // }

    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::readMed(int verb) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    if(medFileName.empty()) {
      ierr = getFileNameFromCommandLine(verb); CHKERRQ(ierr);
    }
    ierr = readMed(medFileName,verb); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MedInterface::writeMed(const string &file,int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    SETERRQ(m_field.get_comm(),MOFEM_NOT_IMPLEMENTED,"Not yet implemented");
    PetscFunctionReturn(0);
  }



}


#endif //WITH_MED
