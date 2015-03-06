/** \file ElasticMaterials.hpp 
 * \brief Elastic materials
 *
 * This file is part of MoFEM.
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

struct ElasticMaterials {

  FieldInterface &mField;
  boost::ptr_map<string,NonlinearElasticElement> elasticElements;
  boost::ptr_map<string,FunctionsToCalulatePiolaKirchhoffI<adouble> > adoubleMaterialModel;
  boost::ptr_map<string,PostPocOnRefinedMesh> postProcesElements;
  boost::ptr_map<string,PostPorcStress> postProcesStress;
  boost::ptr_map<string,FunctionsToCalulatePiolaKirchhoffI<double> > doubleMaterialModel;

  ElasticMaterials(FieldInterface &m_field):
    mField(m_field) {}

  struct BlockOptionData {
    string mAterial;
    int oRder;
    double yOung;
    double pOisson;
    double initTemp;
    BlockOptionData():
      mAterial(MAT_KIRCHOFF),
      oRder(-1),
      yOung(-1),
      pOisson(-1),
      initTemp(0) {}
  };
  map<int,BlockOptionData> blockData;

  /** \brief init Elastic materials declaration for blocks and meshsets

    User has to include in file header:
    \code 
    #include <boost/program_options.hpp>
    using namespace std;
    namespace po = boost::program_options;
    \endcode

    */
  PetscErrorCode readConfigFile() {
    char config_file[255];
    PetscBool is_config_set;
    ierr = PetscOptionsGetString(PETSC_NULL,"-elastic_material_config",config_file,255,&is_config_set); CHKERRQ(ierr);
    po::variables_map vm;
    po::options_description config_file_options;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      ostringstream str_material;
      str_material << "block_" << it->get_msId() << ".material";
      config_file_options.add_options()
	(str_cond.str().c_str(),po::value<string>(&blockData[it->get_msId()].material)->default_value(MAT_KIRCHOFF));
      ostringstream str_cond;
      str_cond << "block_" << it->get_msId() << ".young_modulus";
      config_file_options.add_options()
	(str_cond.str().c_str(),po::value<double>(&blockData[it->get_msId()].yOung)->default_value(-1));
      ostringstream str_capa;
      str_capa << "block_" << it->get_msId() << ".poisson_ratio";
      config_file_options.add_options()
	(str_capa.str().c_str(),po::value<double>(&blockData[it->get_msId()].pOisson)->default_value(-1));
      //ostringstream str_init_temp;
      //str_init_temp << "block_" << it->get_msId() << ".initail_temperature";
      //config_file_options.add_options()
	//(str_init_temp.str().c_str(),po::value<double>(&blockData[it->get_msId()].initTemp)->default_value(0));
    }
    po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
    store(parsed,vm);
    po::notify(vm); 
    vector<string> additional_parameters;
    additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
    for(vector<string>::iterator vit = additional_parameters.begin();
      vit!=additional_parameters.end();vit++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARRNING Unrecognised option %s\n",vit->c_str()); CHKERRQ(ierr);
    }
  }

  PetscErrorCode setBlockOrder() {
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      if(block_data[it->get_msId()].oRder == -1) continue;
      if(block_data[it->get_msId()].oRder == order) continue;
      PetscPrintf(mField.get_comm(),"Set block %d oRder to %d\n",it->get_msId(),block_data[it->get_msId()].oRder);
      Range block_ents;
      rval = mField.get_moab().get_entities_by_handle(it->meshset,block_ents,true); CHKERR(rval);
      Range ents_to_set_order;
      ierr = mField.get_moab().get_adjacencies(block_ents,3,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
      ierr = mField.get_moab().get_adjacencies(block_ents,2,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      ierr = mField.get_moab().get_adjacencies(block_ents,1,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
      if(mField.check_field("DISPLACEMENT")) {
	ierr = mField.set_field_order(ents_to_set_order,"DISPLACEMENT",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      if(mField.check_field("SPATIAL_POSITION")) {
	ierr = mField.set_field_order(ents_to_set_order,"DISPLACEMENT",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      if(mField.check_field("DOT_SPATIAL_POSITION")) {
	ierr = mField.set_field_order(ents_to_set_order,"DISPLACEMENT",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErroCode setElasticElements(int tag_start = 1000) {
    PetscFunctionBegin;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      int id = it->get_msId();
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      EntityHandle meshset = it->get_meshset();
      if(elasticElements.find(blockData.mAterial)==elasticElements.end()) {
	elasticElements.insert(blockData.mAterial,new NonlinearElasticElement(mField,tag_start+id));
	postProcesElasticElements.insert(blockData.mAterial,new 
      }
      rval = mField.get_moab().get_entities_by_type(meshset,MBTET,
	elasticElements[blockData[id].mAterial]->setOfBlocks[id].tEts,true); CHKERR_PETSC(rval);
      elasticElements[blockData[id].mAterial]->setOfBlocks[id].iD = id;
      elasticElements[blockData[id].mAterial]->setOfBlocks[id].E = blockData[id].yOung == -1 ? mydata.data.Young : blockData[id].yOung;
      elasticElements[blockData[id].mAterial]->setOfBlocks[id].E = blockData[id].PoissonRatio == -1 ? mydata.data.Poisson : blockData[id].pOisson;
    }
    boost::ptr_map<string,NonlinearElasticElement>::iterator mit = elasticElements.begin();
    for(;mit!=elasticElements.end();mit++) {
      ierr = mit->addElement(string("ELASTIC_")+mit->first,"SPATIAL_POSITION"); CHKERRQ(ierr);
      if(mit->first.compare(MAT_KIRCHOFF)==0) {
	adoubleMaterialModel.insert(MAT_KIRCHOFF,new NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<adouble>());
	doubleMaterialModel.insert(MAT_KIRCHOFF,new NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<double>());
	ierr = setOperators(adoubleMaterialModel[MAT_KIRCHOFF],"SPATIAL_POSITION"); CHKERRQ(ierr);
      } else
      if(mit->first.compare(MAT_HOOKE)==0) {
	adoubleMaterialModel.insert(MAT_HOOKE,new Hooke<adouble>());
	doubleMaterialModel.insert(MAT_HOOKE,new Hooke<double>());
	ierr = setOperators(adoubleMaterialModel[MAT_HOOKE],"SPATIAL_POSITION"); CHKERRQ(ierr);
      } else 
      if(mit->first.compare(MAT_NEOHOOKEAN)==0) {
	adoubleMaterialModel.insert(MAT_NEOHOOKEAN,new Hooke<adouble>());
	doubleMaterialModel.insert(MAT_NEOHOOKEAN,new Hooke<double>());
	ierr = setOperators(adoubleMaterialModel[MAT_NEOHOOKEAN],"SPATIAL_POSITION"); CHKERRQ(ierr);
      } else {
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }
    }
    PetscFunctionReturn(0);
  }

};

#endif //__ELASTICMATERIALS_HPP__


