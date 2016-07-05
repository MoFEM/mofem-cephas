/** \file ElasticMaterials.hpp
 * \ingroup nonlinear_elastic_elem
 * \brief Elastic materials
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

#ifndef __ELASTICMATERIALS_HPP__
#define __ELASTICMATERIALS_HPP__

#include <Hooke.hpp>
#include <NeoHookean.hpp>

#define MAT_KIRCHOFF "KIRCHOFF"
#define MAT_HOOKE "HOOKE"
#define MAT_NEOHOOKEAN "NEOHOOKEAN"

/** \brief Manage setting parameters and constitutive equations for nonlinear/linear elastic materials
  * \ingroup nonlinear_elastic_elem
  */
struct ElasticMaterials {

  FieldInterface &mField;
  string defMaterial;
  string configFile;

  bool iNitialized;

  ElasticMaterials(FieldInterface &m_field):
    mField(m_field),defMaterial(MAT_KIRCHOFF),
    configFile("elastic_material.in"),
    iNitialized(false) {}


  boost::ptr_map<std::string,NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble> > aDoubleMaterialModel;
  boost::ptr_map<std::string,NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double> > doubleMaterialModel;

  struct BlockOptionData {
    string mAterial;
    int oRder;
    double yOung;
    double pOisson;
    double dEnsity;
    double dashG;
    double dashPoisson;
    double aX,aY,aZ;
    BlockOptionData():
    mAterial(MAT_KIRCHOFF),
    oRder(-1),
    yOung(-1),
    pOisson(-2),
    dEnsity(-1),
    dashG(-1),
    dashPoisson(-1),
    aX(0),
    aY(0),
    aZ(0) {
    }
  };
  std::map<int,BlockOptionData> blockData;

  PetscBool isConfigFileSet;
  po::variables_map vM;


  virtual PetscErrorCode iNit() {
    PetscFunctionBegin;
    //add new material below
    string mat_name;
    mat_name = MAT_KIRCHOFF;
    aDoubleMaterialModel.insert(mat_name,new NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble>());
    doubleMaterialModel.insert(mat_name,new NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double>());
    mat_name = MAT_HOOKE;
    aDoubleMaterialModel.insert(mat_name,new Hooke<adouble>());
    doubleMaterialModel.insert(mat_name,new Hooke<double>());
    mat_name = MAT_NEOHOOKEAN;
    aDoubleMaterialModel.insert(mat_name,new NeoHookean<adouble>());
    doubleMaterialModel.insert(mat_name,new NeoHookean<double>());
    std::ostringstream avilable_materials;
    avilable_materials << "set elastic material < ";
    boost::ptr_map<std::string,NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double> >::iterator mit;
    mit = doubleMaterialModel.begin();
    for(;mit!=doubleMaterialModel.end();mit++) {
      avilable_materials << mit->first << " ";
    }
    avilable_materials << ">";
    PetscErrorCode ierr;
    ierr = PetscOptionsBegin(mField.get_comm(),"","Elastic Materials Configuration","none"); CHKERRQ(ierr);
    char default_material[255];
    PetscBool def_mat_set;
    ierr = PetscOptionsString("-default_material",avilable_materials.str().c_str(),"",MAT_KIRCHOFF,default_material,255,&def_mat_set); CHKERRQ(ierr);
    if(def_mat_set) {
      defMaterial = default_material;
      if(aDoubleMaterialModel.find(defMaterial)==aDoubleMaterialModel.end()) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"material <%s> not implemented",default_material);
      }
    }
    char config_file[255];
    ierr = PetscOptionsString("-elastic_material_configuration","elastic materials configure file name","",configFile.c_str(),config_file,255,&isConfigFileSet); CHKERRQ(ierr);
    if(isConfigFileSet) {
      configFile = config_file;

    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \brief read Elastic materials declaration for blocks and meshsets

    User has to include in file header:
    \code
    #include <boost/program_options.hpp>
    using namespace std;
    namespace po = boost::program_options;
    \endcode

    File parameters:
    \code
    [block_1]
    displacemet_order = 1/2 .. N
    material = KIRCHOFF/HOOKE/NEOHOOKEAN
    young_modulus = 1
    poisson_ratio = 0.25
    density = 1
    a_x = 0
    a_y = 0
    a_z = 10
    \endcode

    To read material configuration file you need to use option:
    \code
    -elastic_material_configuration name_of_config_file
    \endcode

    */
  virtual PetscErrorCode readConfigFile() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    try {

      po::options_description config_file_options;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {

        std::ostringstream str_order;
        str_order << "block_" << it->get_msId() << ".displacemet_order";
        config_file_options.add_options()
        (str_order.str().c_str(),po::value<int>(&blockData[it->get_msId()].oRder)->default_value(-1));

        std::ostringstream str_material;
        str_material << "block_" << it->get_msId() << ".material";
        config_file_options.add_options()
        (str_material.str().c_str(),po::value<std::string>(&blockData[it->get_msId()].mAterial)->default_value(defMaterial));

        std::ostringstream str_ym;
        str_ym << "block_" << it->get_msId() << ".young_modulus";
        config_file_options.add_options()
        (str_ym.str().c_str(),po::value<double>(&blockData[it->get_msId()].yOung)->default_value(-1));

        std::ostringstream str_pr;
        str_pr << "block_" << it->get_msId() << ".poisson_ratio";
        config_file_options.add_options()
        (str_pr.str().c_str(),po::value<double>(&blockData[it->get_msId()].pOisson)->default_value(-2));

        std::ostringstream str_density;
        str_density << "block_" << it->get_msId() << ".density";
        config_file_options.add_options()
        (str_density.str().c_str(),po::value<double>(&blockData[it->get_msId()].dEnsity)->default_value(-1));

        std::ostringstream str_dashG;
        str_dashG << "block_" << it->get_msId() << ".dashG";
        config_file_options.add_options()
        (str_dashG.str().c_str(),po::value<double>(&blockData[it->get_msId()].dashG)->default_value(-1));

        std::ostringstream str_dashPoisson;
        str_dashPoisson << "block_" << it->get_msId() << ".dashPoisson";
        config_file_options.add_options()
        (str_dashPoisson.str().c_str(),po::value<double>(&blockData[it->get_msId()].dashPoisson)->default_value(-2));

        std::ostringstream str_ax;
        str_ax << "block_" << it->get_msId() << ".a_x";
        config_file_options.add_options()
        (str_ax.str().c_str(),po::value<double>(&blockData[it->get_msId()].aX)->default_value(0));

        std::ostringstream str_ay;
        str_ay << "block_" << it->get_msId() << ".a_y";
        config_file_options.add_options()
        (str_ay.str().c_str(),po::value<double>(&blockData[it->get_msId()].aY)->default_value(0));

        std::ostringstream str_az;
        str_az << "block_" << it->get_msId() << ".a_z";
        config_file_options.add_options()
        (str_az.str().c_str(),po::value<double>(&blockData[it->get_msId()].aZ)->default_value(0));
      }
      ifstream file(configFile.c_str());
      if(isConfigFileSet) {
        if(!file.good()) {
          SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"file < %s > not found",configFile.c_str());
        }
      }
      po::parsed_options parsed = parse_config_file(file,config_file_options,true);
      store(parsed,vM);
      po::notify(vM);
      std::vector<std::string> additional_parameters;
      additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
      for(std::vector<std::string>::iterator vit = additional_parameters.begin();
      vit!=additional_parameters.end();vit++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARNING Unrecognized option %s\n",vit->c_str()); CHKERRQ(ierr);
      }
    } catch (std::exception& ex) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,"error parsing material elastic configuration file");
    }
    PetscFunctionReturn(0);

  }

  PetscErrorCode setBlocksOrder() {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
    //set app. order
    PetscBool flg = PETSC_TRUE;
    PetscInt disp_order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-order",&disp_order,&flg); CHKERRQ(ierr);
    if(flg!=PETSC_TRUE) {
      disp_order = 1;
    }
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(blockData[it->get_msId()].oRder == -1) continue;
      if(blockData[it->get_msId()].oRder == disp_order) continue;
      PetscPrintf(mField.get_comm(),"Set block %d oRder to %d\n",it->get_msId(),blockData[it->get_msId()].oRder);
      Range block_ents;
      rval = mField.get_moab().get_entities_by_handle(it->meshset,block_ents,true); CHKERR_MOAB(rval);
      Range ents_to_set_order;
      ierr = mField.get_moab().get_adjacencies(block_ents,3,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
      ierr = mField.get_moab().get_adjacencies(block_ents,2,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      ierr = mField.get_moab().get_adjacencies(block_ents,1,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      if(mField.check_field("DISPLACEMENT")) {
        ierr = mField.set_field_order(ents_to_set_order,"DISPLACEMENT",blockData[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      if(mField.check_field("SPATIAL_POSITION")) {
        ierr = mField.set_field_order(ents_to_set_order,"SPATIAL_POSITION",blockData[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      if(mField.check_field("DOT_SPATIAL_POSITION")) {
        ierr = mField.set_field_order(ents_to_set_order,"DOT_SPATIAL_POSITION",blockData[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      if(mField.check_field("SPATIAL_VELOCITY")) {
        ierr = mField.set_field_order(ents_to_set_order,"SPATIAL_VELOCITY",blockData[it->get_msId()].oRder); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  #ifdef __NONLINEAR_ELASTIC_HPP

  virtual PetscErrorCode setBlocks(std::map<int,NonlinearElasticElement::BlockData> &set_of_blocks) {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
    if(!iNitialized) {
      ierr = iNit(); CHKERRQ(ierr);
      ierr = readConfigFile(); CHKERRQ(ierr);
      ierr = setBlocksOrder(); CHKERRQ(ierr);
      iNitialized = true;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      int id = it->get_msId();
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      EntityHandle meshset = it->getMeshSet();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,set_of_blocks[id].tEts,true); CHKERRQ_MOAB(rval);
      set_of_blocks[id].iD = id;
      set_of_blocks[id].E = mydata.data.Young;
      if(blockData[id].yOung >= 0) set_of_blocks[id].E = blockData[id].yOung;
      if(blockData[id].pOisson >= -1) set_of_blocks[id].PoissonRatio = blockData[id].pOisson;
      PetscPrintf(mField.get_comm(),"Block Id %d Young Modulus %3.2g Poisson Ration %3.2f Material model %s\n",
      id,set_of_blocks[id].E,set_of_blocks[id].PoissonRatio,blockData[id].mAterial.c_str());
      if(blockData[id].mAterial.compare(MAT_KIRCHOFF)==0) {
        set_of_blocks[id].materialDoublePtr = &doubleMaterialModel.at(MAT_KIRCHOFF);
        set_of_blocks[id].materialAdoublePtr = &aDoubleMaterialModel.at(MAT_KIRCHOFF);
      } else
      if(blockData[id].mAterial.compare(MAT_HOOKE)==0) {
        set_of_blocks[id].materialDoublePtr = &doubleMaterialModel.at(MAT_HOOKE);
        set_of_blocks[id].materialAdoublePtr = &aDoubleMaterialModel.at(MAT_HOOKE);
      } else
      if(blockData[id].mAterial.compare(MAT_NEOHOOKEAN)==0) {
        set_of_blocks[id].materialDoublePtr = &doubleMaterialModel.at(MAT_NEOHOOKEAN);
        set_of_blocks[id].materialAdoublePtr = &aDoubleMaterialModel.at(MAT_NEOHOOKEAN);
      } else {
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }
    }
    PetscFunctionReturn(0);
  }

  #endif //__NONLINEAR_ELASTIC_HPP

  #ifdef __CONVECTIVE_MASS_ELEMENT_HPP

  PetscErrorCode setBlocks(std::map<int,ConvectiveMassElement::BlockData> &set_of_blocks) {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
    if(!iNitialized) {
      ierr = iNit(); CHKERRQ(ierr);
      ierr = readConfigFile(); CHKERRQ(ierr);
      ierr = setBlocksOrder(); CHKERRQ(ierr);
      iNitialized = true;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|BODYFORCESSET,it)) {
      int id = it->get_msId();
      EntityHandle meshset = it->getMeshSet();
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,set_of_blocks[id].tEts,true); CHKERRQ_MOAB(rval);
      Block_BodyForces mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      set_of_blocks[id].rho0 = mydata.data.density;
      set_of_blocks[id].a0.resize(3);
      set_of_blocks[id].a0[0] = mydata.data.acceleration_x;
      set_of_blocks[id].a0[1] = mydata.data.acceleration_y;
      set_of_blocks[id].a0[2] = mydata.data.acceleration_z;
      if(blockData[id].dEnsity>=0) {
        set_of_blocks[id].rho0 = blockData[id].dEnsity;
        std::ostringstream str_ax;
        str_ax << "block_" << it->get_msId() << ".a_x";
        std::ostringstream str_ay;
        str_ay << "block_" << it->get_msId() << ".a_y";
        std::ostringstream str_az;
        str_az << "block_" << it->get_msId() << ".a_z";
        if(vM.count(str_ax.str().c_str())) {
          set_of_blocks[id].a0[0] = blockData[id].aX;
        }
        if(vM.count(str_ay.str().c_str())) {
          set_of_blocks[id].a0[1] = blockData[id].aY;
        }
        if(vM.count(str_az.str().c_str())) {
          set_of_blocks[id].a0[2] = blockData[id].aZ;
        }
      }
      PetscPrintf(
        mField.get_comm(),"Block Id %d Density %3.2g a_x = %3.2g a_y = %3.2g a_z = %3.2g\n",
        id,set_of_blocks[id].rho0,set_of_blocks[id].a0[0],set_of_blocks[id].a0[1],set_of_blocks[id].a0[2]
      );
    }

    PetscFunctionReturn(0);
  }

  #endif //__CONVECTIVE_MASS_ELEMENT_HPP

  #ifdef __KELVIN_VOIGT_DAMPER_HPP__

  PetscErrorCode setBlocks(std::map<int,KelvinVoigtDamper::BlockMaterialData> &set_of_blocks) {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;

    if(!iNitialized) {
      ierr = iNit(); CHKERRQ(ierr);
      ierr = readConfigFile(); CHKERRQ(ierr);
      ierr = setBlocksOrder(); CHKERRQ(ierr);
      iNitialized = true;
    }

    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      bool set = false;
      int id = it->get_msId();
      EntityHandle meshset = it->getMeshSet();
      if(it->getName().compare(0,6,"DAMPER") == 0) {
        set = true;
        std::vector<double> data;
        ierr = it->get_attributes(data); CHKERRQ(ierr);
        if(data.size()<2) {
          SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
        }
        rval = mField.get_moab().get_entities_by_type(
          it->meshset,MBTET,set_of_blocks[it->get_msId()].tEts,true
        ); CHKERRQ_MOAB(rval);
        set_of_blocks[it->get_msId()].gBeta = data[0];
        set_of_blocks[it->get_msId()].vBeta = data[1];
      }
      if(blockData[id].dashG > 0) {
        set = true;
        Range tEts;
        rval = mField.get_moab().get_entities_by_type(meshset,MBTET,tEts,true); CHKERRQ_MOAB(rval);
        if(tEts.empty()) continue;
        set_of_blocks[it->get_msId()].tEts = tEts;
        set_of_blocks[it->get_msId()].gBeta = blockData[id].dashG;
        set_of_blocks[it->get_msId()].vBeta = blockData[id].dashPoisson;
      }
      if(set) {
        PetscPrintf(
          mField.get_comm(),
          "Block Id %d Damper Shear Modulus = %3.2g Poisson ratio = %3.2g\n",
          id,set_of_blocks[id].gBeta,set_of_blocks[id].vBeta
        );
      }
    }

    PetscFunctionReturn(0);
  }

  #endif //__KELVIN_VOIGT_DAMPER_HPP__

};

#endif //__ELASTICMATERIALS_HPP__
