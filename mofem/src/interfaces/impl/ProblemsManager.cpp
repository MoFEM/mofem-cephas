/** \file ProblemCore.cpp
 * \brief Managing complexities for problem
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
#include <Core.hpp>

#include <ProblemsManager.hpp>

#include <moab/MeshTopoUtil.hpp>

namespace MoFEM {

  PetscErrorCode ProblemsManager::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMProblemsManager) {
      *iface = dynamic_cast<ProblemsManager*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  ProblemsManager::ProblemsManager(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)) {
    PetscLogEventRegister("MoFEMProblemsManager",0,&USER_EVENT_ProblemsManager);
  }
  ProblemsManager::~ProblemsManager() {
  }

  PetscErrorCode ProblemsManager::partitionMesh(
    const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb
  ) {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;

    //get layout
    int rstart,rend,nb_elems;
    {
      PetscLayout layout;
      ierr = PetscLayoutCreate(m_field.get_comm(),&layout); CHKERRQ(ierr);
      ierr = PetscLayoutSetBlockSize(layout, 1); CHKERRQ(ierr);
      ierr = PetscLayoutSetSize(layout,ents.size()); CHKERRQ(ierr);
      ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
      ierr = PetscLayoutGetSize(layout,&nb_elems); CHKERRQ(ierr);
      ierr = PetscLayoutGetRange(layout,&rstart,&rend); CHKERRQ(ierr);
      ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
      if (verb > 0) {
        PetscSynchronizedPrintf(
          m_field.get_comm(),
          "Finite elements partition in problem: row lower %d row upper %d nb elems %d\n",
          rstart,rend,nb_elems
        );
        PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
      }
    }

    std::map<EntityHandle,int> problem_fe_ents;
    {
      Range::iterator eit = ents.begin();
      for(int ii = 0;eit!=ents.end();eit++,ii++) {
        problem_fe_ents[*eit] = ii;
      }
    }

    int *_i;
    int *_j;
    {
      MeshTopoUtil mtu(&m_field.get_moab());
      std::vector<int> i(rend-rstart+1,0),j;
      {
        int jj = 0;
        Range::iterator fe_it = ents.begin();
        for(int ii = 0;fe_it!=ents.end();fe_it++,ii++) {
          if(ii < rstart) continue;
          if(ii >= rend) break;
          if(m_field.get_moab().type_from_handle(*fe_it)==MBENTITYSET) {
            SETERRQ(
              m_field.get_comm(),
              MOFEM_NOT_IMPLEMENTED,
              "not yet implemented, don't know what to do for meshset element"
            );
          } else {
            Range adj_ents;
            rval = mtu.get_bridge_adjacencies(*fe_it,adj_dim,dim,adj_ents); CHKERRQ_MOAB(rval);
            adj_ents = intersect(adj_ents,ents);
            i[jj] = j.size();
            for(Range::iterator eit = adj_ents.begin();eit!=adj_ents.end();eit++) {
              if(*eit==*fe_it) continue; // no diagonal
              j.push_back(problem_fe_ents[*eit]);
            }
          }
          jj++;
        }
        i[jj] = j.size();
      }
      ierr = PetscMalloc(i.size()*sizeof(int),&_i); CHKERRQ(ierr);
      ierr = PetscMalloc(j.size()*sizeof(int),&_j); CHKERRQ(ierr);
      copy(i.begin(),i.end(),_i);
      copy(j.begin(),j.end(),_j);
    }

    {
      Mat Adj;
      // Adjacency matrix used to partition problems, f.e. METIS
      ierr = MatCreateMPIAdj(
        m_field.get_comm(),rend-rstart,nb_elems,_i,_j,PETSC_NULL,&Adj
      ); CHKERRQ(ierr);
      ierr = MatSetOption(Adj,MAT_STRUCTURALLY_SYMMETRIC,PETSC_TRUE); CHKERRQ(ierr);

      // if(1) {
      //   Mat A;
      //   MatConvert(Adj,MATMPIAIJ,MAT_INITIAL_MATRIX,&A);
      //   MatView(A,PETSC_VIEWER_DRAW_WORLD);
      //   std::string wait;
      //   std::cin >> wait;
      //   MatDestroy(&A);
      // }

      // run pets to do partitioning
      MatPartitioning part;
      IS is;
      ierr = MatPartitioningCreate(m_field.get_comm(),&part); CHKERRQ(ierr);
      ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
      ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
      ierr = MatPartitioningSetNParts(part,n_parts); CHKERRQ(ierr);
      ierr = MatPartitioningApply(part,&is); CHKERRQ(ierr);

      //gather
      IS is_gather,is_num,is_gather_num;
      ierr = ISAllGather(is,&is_gather); CHKERRQ(ierr);
      ierr = ISPartitioningToNumbering(is,&is_num); CHKERRQ(ierr);
      ierr = ISAllGather(is_num,&is_gather_num); CHKERRQ(ierr);

      const int *part_number,*gids;
      ierr = ISGetIndices(is_gather,&part_number);  CHKERRQ(ierr);
      ierr = ISGetIndices(is_gather_num,&gids);  CHKERRQ(ierr);

      // set partition tag and gid tag to entities
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
      Tag gid_tag;
      Tag part_tag = pcomm->part_tag();
      {
        const int zero =  0;
        rval = m_field.get_moab().tag_get_handle(
          GLOBAL_ID_TAG_NAME,1,MB_TYPE_INTEGER,gid_tag,MB_TAG_DENSE|MB_TAG_CREAT,&zero
        ); CHKERRQ_MOAB(rval);
        // get any sets already with this tag, and clear them
        rval = m_field.get_moab().tag_set_data(part_tag,ents,part_number); CHKERRQ_MOAB(rval);
        // rval = moab.tag_set_data(gid_tag,ents,&gids[0]); CHKERRQ_MOAB(rval);
        // std::vector<int> add_one(ents.size());
        // for(int ii = 0;ii<ents.size();ii++) {
        //   add_one[ii] = gids[ii]+1;
        // }
        // rval = moab.tag_set_data(gid_tag,ents,&add_one[0]); CHKERRQ_MOAB(rval);
      }

      std::map<int,Range> parts_ents;
      {
        // get entities on each part
        Range::iterator eit = ents.begin();
        for(int ii = 0;eit!=ents.end();eit++,ii++) {
          parts_ents[part_number[ii]].insert(*eit);
        }
        Range tagged_sets;
        rval = m_field.get_moab().get_entities_by_type_and_tag(
          0,MBENTITYSET,&part_tag,NULL,1,tagged_sets,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        if(!tagged_sets.empty()) {
          rval = m_field.get_moab().tag_delete_data(part_tag,tagged_sets); CHKERRQ_MOAB(rval);
          // rval = moab.delete_entities(tagged_sets); CHKERRQ_MOAB(rval);
        }
        if(n_parts > (int) tagged_sets.size()) {
          // too few partition sets - create missing ones
          int num_new = n_parts - tagged_sets.size();
          for(int i = 0;i < num_new;i++) {
            EntityHandle new_set;
            rval = m_field.get_moab().create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER, new_set); CHKERR_MOAB(rval);
            tagged_sets.insert(new_set);
          }
        } else if (n_parts < (int)tagged_sets.size()) {
          // too many partition sets - delete extras
          int num_del = tagged_sets.size() - n_parts;
          for(int i = 0; i < num_del; i++) {
            EntityHandle old_set = tagged_sets.pop_back();
            rval = m_field.get_moab().delete_entities(&old_set, 1); CHKERR_MOAB(rval);
          }
        }
        // write a tag to those sets denoting they're partition sets, with a value of the
        // proc number
        std::vector<int> dum_ids(n_parts);
        for(int i = 0;i<n_parts;i++) dum_ids[i] = i;
        rval = m_field.get_moab().tag_set_data(part_tag,tagged_sets,&*dum_ids.begin()); CHKERR_MOAB(rval);
        // for(int i = 0;i<n_parts;i++) {
        //   rval = m_field.get_moab().add_entities(tagged_sets[i],parts_ents[i]); CHKERR_MOAB(rval);
        // }

        // get lower dimension entities on each part
        for(int pp = 0;pp!=n_parts;pp++) {
          Range dim_ents = parts_ents[pp].subset_by_dimension(dim);
          for(int dd = dim-1;dd!=-1;dd--) {
            Range adj_ents;
            if(dim > 0 ) {
              rval = m_field.get_moab().get_adjacencies(
                dim_ents,dd,true,adj_ents,moab::Interface::UNION
              ); CHKERRQ_MOAB(rval);
            } else {
              rval = m_field.get_moab().get_connectivity(dim_ents,adj_ents,true); CHKERRQ_MOAB(rval);
            }
            parts_ents[pp].merge(adj_ents);
            // std::cerr << pp << " add " << parts_ents[pp].size() << std::endl;
          }
        }
        for(int pp = 1;pp!=n_parts;pp++) {
          for(int ppp = 0;ppp!=pp;ppp++) {
            // std::cerr << pp << "<-" << ppp << " " << parts_ents[pp].size() << " " << parts_ents[ppp].size();
            parts_ents[pp] = subtract(parts_ents[pp],parts_ents[ppp]);
            // std::cerr << " " << parts_ents[pp].size() << std::endl;
          }
        }
        for(int pp = 0;pp!=n_parts;pp++) {
          rval = m_field.get_moab().add_entities(tagged_sets[pp],parts_ents[pp]); CHKERR_MOAB(rval);
        }

        // set gid to lower dimension entities
        for(int dd = 0;dd<=dim;dd++) {
          int gid = 1; // moab indexing from 1
          for(int pp = 0;pp!=n_parts;pp++) {
            Range dim_ents = parts_ents[pp].subset_by_dimension(dd);
            // std::cerr << dim_ents.size() << " " << dd  << " " << pp << std::endl;
            for(Range::iterator eit = dim_ents.begin();eit!=dim_ents.end();eit++) {
              if(dd>0) {
                rval = m_field.get_moab().tag_set_data(part_tag,&*eit,1,&pp); CHKERRQ_MOAB(rval);
              }
              rval = m_field.get_moab().tag_set_data(gid_tag,&*eit,1,&gid); CHKERRQ_MOAB(rval);
              gid++;
            }
          }
        }

      }

      ierr = ISRestoreIndices(is_gather,&part_number);  CHKERRQ(ierr);
      ierr = ISRestoreIndices(is_gather_num,&gids);  CHKERRQ(ierr);
      ierr = ISDestroy(&is_num); CHKERRQ(ierr);
      ierr = ISDestroy(&is_gather_num); CHKERRQ(ierr);
      ierr = ISDestroy(&is_gather); CHKERRQ(ierr);
      ierr = ISDestroy(&is); CHKERRQ(ierr);
      ierr = MatPartitioningDestroy(&part); CHKERRQ(ierr);
      ierr = MatDestroy(&Adj); CHKERRQ(ierr);
    }

    cOre.getBuildMoFEM() |= Core::PARTITION_MESH;

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::buildProblem(const std::string &name,const bool square_matrix,int verb) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    if(!(cOre.getBuildMoFEM()&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
    if(!(cOre.getBuildMoFEM()&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
    if(!(cOre.getBuildMoFEM()&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
    const MoFEMProblem *problem_ptr;
    ierr = m_field.get_problem(name,&problem_ptr); CHKERRQ(ierr);
    ierr = buildProblem(const_cast<MoFEMProblem*>(problem_ptr),square_matrix,verb); CHKERRQ(ierr);
    cOre.getBuildMoFEM() |= 1<<3; // It is assumed that user who uses this function knows what he is doing
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::buildProblem(MoFEMProblem *problem_ptr,const bool square_matrix,int verb) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const EntFiniteElement_multiIndex *fe_ent_ptr;
    const DofEntity_multiIndex *dofs_field_ptr;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_ProblemsManager,0,0,0,0);

    // Note: Only allowed changes on problem_ptr structure which not influence multi-index.

    if(problem_ptr->getBitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->getName().c_str());
    }
    if(problem_ptr->getBitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->getName().c_str());
    }
    ierr = m_field.clear_problem(problem_ptr->getName()); CHKERRQ(ierr);
    ierr = m_field.get_ents_finite_elements(&fe_ent_ptr); CHKERRQ(ierr);
    ierr = m_field.get_dofs(&dofs_field_ptr); CHKERRQ(ierr);

    //zero finite elements
    problem_ptr->numeredFiniteElements.clear();

    DofEntity_multiIndex_active_view dofs_rows,dofs_cols;
    {
      EntFiniteElement_multiIndex::iterator miit = fe_ent_ptr->begin();
      EntFiniteElement_multiIndex::iterator hi_miit = fe_ent_ptr->end();
      //iterate all finite element entities in database
      for(;miit!=hi_miit;miit++) {
        //if element is in problem
        if(((*miit)->getId()&problem_ptr->getBitFEId()).any()) {
          //if finite element bit level has all refined bits sets
          if(
            ((*miit)->getBitRefLevel()&problem_ptr->getBitRefLevel()) ==
            problem_ptr->getBitRefLevel()
          ) {
            //get dof uids for rows and columns
            ierr = (*miit)->getRowDofView(*dofs_field_ptr,dofs_rows); CHKERRQ(ierr);
            if(!square_matrix) {
              ierr = (*miit)->getColDofView(*dofs_field_ptr,dofs_cols); CHKERRQ(ierr);
            }
          }
        }
      }
    }

    // Add row dofs to problem
    {
      // zero rows
      *problem_ptr->tag_nbdof_data_row = 0;
      *problem_ptr->tag_local_nbdof_data_row = 0;
      *problem_ptr->tag_ghost_nbdof_data_row = 0;

      //add dofs for rows
      DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,hi_miit;
      hi_miit = dofs_rows.get<0>().end();

      int count_dofs = 0;
      miit = dofs_rows.get<0>().begin();
      for(;miit!=hi_miit;miit++) {
        if(
          !(*miit)->getActive()||
          ((*miit)->getBitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=
          (*miit)->getBitRefLevel()
        ) {
          continue;
        }
        ++count_dofs;
      }

      boost::shared_ptr<std::vector<NumeredDofEntity> > dofs_array =
      boost::shared_ptr<std::vector<NumeredDofEntity> >(new std::vector<NumeredDofEntity>());
      problem_ptr->getRowDofsSeqence() = dofs_array;
      dofs_array->reserve(count_dofs);
      std::vector<boost::shared_ptr<NumeredDofEntity> > dofs_shared_array;
      dofs_shared_array.reserve(count_dofs);
      miit = dofs_rows.get<0>().begin();
      for(;miit!=hi_miit;miit++) {
        if(
          !(*miit)->getActive()||
          ((*miit)->getBitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=
          (*miit)->getBitRefLevel()
        ) {
          continue;
        }
        dofs_array->push_back(NumeredDofEntity(*miit));
        dofs_shared_array.push_back(
          boost::shared_ptr<NumeredDofEntity>(dofs_array,&dofs_array->back())
        );
        dofs_array->back().dofIdx = (*problem_ptr->tag_nbdof_data_row)++;
      }
      problem_ptr->numered_dofs_rows->insert(dofs_shared_array.begin(),dofs_shared_array.end());
    }

    // Add col dofs to problem
    if(!square_matrix) {
      //zero cols
      *problem_ptr->tag_nbdof_data_col = 0;
      *problem_ptr->tag_local_nbdof_data_col = 0;
      *problem_ptr->tag_ghost_nbdof_data_col = 0;

      //add dofs for cols
      DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,hi_miit;
      hi_miit = dofs_cols.get<0>().end();

      int count_dofs = 0;
      miit = dofs_cols.get<0>().begin();
      for(;miit!=hi_miit;miit++) {
        if(
          !(*miit)->getActive()||
          ((*miit)->getBitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=
          (*miit)->getBitRefLevel()
        ) {
          continue;
        }
        count_dofs++;
      }

      boost::shared_ptr<std::vector<NumeredDofEntity> > dofs_array =
      boost::shared_ptr<std::vector<NumeredDofEntity> >(new std::vector<NumeredDofEntity>());
      problem_ptr->getColDofsSeqence() = dofs_array;
      dofs_array->reserve(count_dofs);
      std::vector<boost::shared_ptr<NumeredDofEntity> > dofs_shared_array;
      dofs_shared_array.reserve(count_dofs);
      miit = dofs_cols.get<0>().begin();
      for(;miit!=hi_miit;miit++) {
        if(
          !(*miit)->getActive()||
          ((*miit)->getBitRefLevel()&problem_ptr->get_DofMask_BitRefLevel())!=
          (*miit)->getBitRefLevel()
        ) {
          continue;
        }
        dofs_array->push_back(NumeredDofEntity(*miit));
        dofs_shared_array.push_back(
          boost::shared_ptr<NumeredDofEntity>(dofs_array,&dofs_array->back())
        );
        dofs_array->back().dofIdx = (*problem_ptr->tag_nbdof_data_col)++;
      }
      problem_ptr->numered_dofs_cols->insert(
        dofs_shared_array.begin(),dofs_shared_array.end()
      );
    } else {
      problem_ptr->numered_dofs_cols = problem_ptr->numered_dofs_rows;
      *(problem_ptr->tag_local_nbdof_data_col) = *(problem_ptr->tag_local_nbdof_data_row);
      *(problem_ptr->tag_nbdof_data_col) = *(problem_ptr->tag_nbdof_data_row);
    }

    //job done, some debugging and postprocessing
    if(verb>0) {
      PetscSynchronizedPrintf(m_field.get_comm(),"Problem %s Nb. rows %u Nb. cols %u\n",
      problem_ptr->getName().c_str(),
      problem_ptr->numered_dofs_rows->size(),problem_ptr->numered_dofs_cols->size());
    }
    if(verb>1) {
      EntFiniteElement_multiIndex::iterator miit = fe_ent_ptr->begin();
      EntFiniteElement_multiIndex::iterator hi_miit = fe_ent_ptr->end();
      std::ostringstream ss;
      ss << "rank " << m_field.getCommRank() << " ";
      ss << "FEs data for problem " << *problem_ptr << std::endl;
      for(;miit!=hi_miit;miit++) {
        ss << "rank " << m_field.getCommRank() << " ";
        ss << **miit << std::endl;
      }
      ss << "rank " << m_field.getCommRank() << " ";
      ss << "FEs row dofs "<< *problem_ptr << " Nb. row dof " << problem_ptr->getNbDofsRow() << std::endl;
      NumeredDofEntity_multiIndex::iterator miit_dd_row = problem_ptr->numered_dofs_rows->begin();
      for(;miit_dd_row!=problem_ptr->numered_dofs_rows->end();miit_dd_row++) {
        ss << "rank " << m_field.getCommRank() << " ";
        ss<<**miit_dd_row<<std::endl;
      }
      ss << "rank " << m_field.getCommRank() << " ";
      ss << "FEs col dofs "<< *problem_ptr << " Nb. col dof " << problem_ptr->getNbDofsCol() << std::endl;
      NumeredDofEntity_multiIndex::iterator miit_dd_col = problem_ptr->numered_dofs_cols->begin();
      for(;miit_dd_col!=problem_ptr->numered_dofs_cols->end();miit_dd_col++) {
        ss << "rank " << m_field.getCommRank() << " ";
        ss<<**miit_dd_col<<std::endl;
      }
      PetscSynchronizedPrintf(m_field.get_comm(),ss.str().c_str());
    }

    if(verb>0) {
      PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
    }
    cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM; // It is assumed that user who uses this function knows what he is doing

    PetscLogEventEnd(USER_EVENT_ProblemsManager,0,0,0,0);

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::buildProblemOnDistributedMesh(
    const std::string &name,const bool square_matrix,int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;

    if(!((cOre.getBuildMoFEM())&Core::BUILD_FIELD)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
    if(!((cOre.getBuildMoFEM())&Core::BUILD_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
    if(!((cOre.getBuildMoFEM())&Core::BUILD_ADJ)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");

    const MoFEMProblem *problem_ptr;
    ierr = m_field.get_problem(name,&problem_ptr); CHKERRQ(ierr);
    ierr = buildProblemOnDistributedMesh(
      const_cast<MoFEMProblem*>(problem_ptr),square_matrix,verb
    ); CHKERRQ(ierr);

    cOre.getBuildMoFEM() |= Core::BUILD_PROBLEM;
    cOre.getBuildMoFEM() |= Core::PARTITION_PROBLEM;

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::buildProblemOnDistributedMesh(
    MoFEMProblem *problem_ptr,const bool square_matrix,int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const EntFiniteElement_multiIndex *fe_ent_ptr;
    const DofEntity_multiIndex *dofs_field_ptr;
    PetscFunctionBegin;
    PetscLogEventBegin(USER_EVENT_ProblemsManager,0,0,0,0);

    ierr = m_field.get_ents_finite_elements(&fe_ent_ptr); CHKERRQ(ierr);
    ierr = m_field.get_dofs(&dofs_field_ptr); CHKERRQ(ierr);

    if(problem_ptr->getBitRefLevel().none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem <%s> refinement level not set",problem_ptr->getName().c_str());
    }

    ierr = m_field.clear_problem(problem_ptr->getName()); CHKERRQ(ierr);

    int loop_size = 2;
    if(square_matrix) {
      loop_size = 1;
      problem_ptr->numered_dofs_cols = problem_ptr->numered_dofs_rows;
    } else if(problem_ptr->numered_dofs_cols == problem_ptr->numered_dofs_rows) {
      problem_ptr->numered_dofs_cols = boost::shared_ptr<NumeredDofEntity_multiIndex>(
        new NumeredDofEntity_multiIndex()
      );
    }

    // get rows and cols dofs view based on data on elements
    DofEntity_multiIndex_active_view dofs_rows,dofs_cols;
    {
      //fe_miit iterator for finite elements
      EntFiniteElement_multiIndex::iterator fe_miit = fe_ent_ptr->begin();
      EntFiniteElement_multiIndex::iterator hi_fe_miit = fe_ent_ptr->end();
      //iterate all finite elemen entities in database
      for(;fe_miit!=hi_fe_miit;fe_miit++) {
        //if element is in problem
        if(((*fe_miit)->getId()&problem_ptr->getBitFEId()).any()) {
          //if finite element bit level has all refined bits sets
          if(((*fe_miit)->getBitRefLevel()&problem_ptr->getBitRefLevel())==problem_ptr->getBitRefLevel()) {
            //get dof uids for rows and columns
            ierr = (*fe_miit)->getRowDofView(*dofs_field_ptr,dofs_rows); CHKERRQ(ierr);
            if(!square_matrix) {
              ierr = (*fe_miit)->getColDofView(*dofs_field_ptr,dofs_cols); CHKERRQ(ierr);
            }
          }
        }
      }
    }

    // get problem bit level
    const BitRefLevel &problem_bit_level = problem_ptr->get_DofMask_BitRefLevel();

    //add dofs for rows and cols and set ownership
    DofEntity_multiIndex_active_view* dofs_ptr[] = { &dofs_rows, &dofs_cols };
    boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_ptr[] = {
      problem_ptr->numered_dofs_rows, problem_ptr->numered_dofs_cols
    };
    int* nbdof_ptr[] = {
      problem_ptr->tag_nbdof_data_row, problem_ptr->tag_nbdof_data_col
    };
    int* local_nbdof_ptr[] = {
      problem_ptr->tag_local_nbdof_data_row, problem_ptr->tag_local_nbdof_data_col
    };
    int* ghost_nbdof_ptr[] = {
      problem_ptr->tag_ghost_nbdof_data_row, problem_ptr->tag_ghost_nbdof_data_col
    };
    for(int ss = 0;ss<2;ss++) {
      *(nbdof_ptr[ss]) = 0;
      *(local_nbdof_ptr[ss]) = 0;
      *(ghost_nbdof_ptr[ss]) = 0;
    }

    //Loop over dofs on rows and columns and add to multi-indices in dofs problem structure,
    //set partition for each dof
    int nb_local_dofs[] = { 0,0 };
    for(int ss = 0;ss<loop_size;ss++) {
      DofEntity_multiIndex_active_view::nth_index<1>::type::iterator miit,hi_miit;
      miit = dofs_ptr[ss]->get<1>().lower_bound(1);
      hi_miit = dofs_ptr[ss]->get<1>().upper_bound(1);
      for(;miit!=hi_miit;miit++) {
        const BitRefLevel &dof_bit_level = (*miit)->getBitRefLevel();
        if((dof_bit_level&problem_bit_level)!=dof_bit_level) {
          continue;
        }
        int owner_proc = (*miit)->getOwnerProc();
        if(owner_proc == m_field.getCommRank()) {
          nb_local_dofs[ss]++;
        }
      }
    }

    //get layout
    int start_ranges[2],end_ranges[2];
    for(int ss = 0;ss!=loop_size;ss++) {
      PetscLayout layout;
      ierr = PetscLayoutCreate(m_field.get_comm(),&layout); CHKERRQ(ierr);
      ierr = PetscLayoutSetBlockSize(layout,1); CHKERRQ(ierr);
      ierr = PetscLayoutSetLocalSize(layout,nb_local_dofs[ss]); CHKERRQ(ierr);
      ierr = PetscLayoutSetUp(layout); CHKERRQ(ierr);
      ierr = PetscLayoutGetSize(layout,&*nbdof_ptr[ss]); CHKERRQ(ierr); // get global size
      ierr = PetscLayoutGetRange(layout,&start_ranges[ss],&end_ranges[ss]); CHKERRQ(ierr); //get ranges
      ierr = PetscLayoutDestroy(&layout); CHKERRQ(ierr);
    }
    if(square_matrix) {
      nbdof_ptr[1] = nbdof_ptr[0];
      nb_local_dofs[1] = nb_local_dofs[0];
      start_ranges[1] = start_ranges[0];
      end_ranges[1] = end_ranges[0];
    }

    // if(sizeof(UId) != SIZEOFUID) {
    //   SETERRQ2(
    //     PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
    //     "check size of UId, size of UId is %u != %u",
    //     sizeof(UId),SIZEOFUID
    //   );
    // }

    // set local and global indices on own dofs
    const size_t idx_data_type_size = sizeof(IdxDataType);
    const size_t data_block_size = idx_data_type_size/sizeof(int);
    if(sizeof(IdxDataType) % sizeof(int)) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    std::vector<std::vector<IdxDataType> > ids_data_packed_rows(m_field.getCommSize()),ids_data_packed_cols(m_field.getCommSize());

    // used to keep shared_ptr
    std::vector<boost::shared_ptr<NumeredDofEntity> > dofs_shared_array;

    // ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    // Loop over dofs on this processor and prepare those dofs to send on another proc
    for(int ss = 0;ss<loop_size;ss++) {

      DofEntity_multiIndex_active_view::nth_index<0>::type::iterator miit,hi_miit;
      hi_miit = dofs_ptr[ss]->get<0>().end();

      boost::shared_ptr<std::vector<NumeredDofEntity> > dofs_array =
      boost::shared_ptr<std::vector<NumeredDofEntity> >(new std::vector<NumeredDofEntity>());
      int nb_dofs_to_add = 0;
      miit = dofs_ptr[ss]->get<0>().begin();
      for(;miit!=hi_miit;miit++) {
        // Only set global idx for dofs on this processor part
        if(!(miit->get()->getActive())) continue;
        const BitRefLevel &dof_bit_level = (*miit)->getBitRefLevel();
        if((dof_bit_level&problem_bit_level)!=dof_bit_level) {
          continue;
        }
        ++nb_dofs_to_add;
      }
      dofs_array->reserve(nb_dofs_to_add);
      dofs_shared_array.clear();
      dofs_shared_array.reserve(nb_dofs_to_add);
      if(ss == 0) {
        problem_ptr->getRowDofsSeqence() = dofs_array;
      } else {
        problem_ptr->getColDofsSeqence() = dofs_array;
      }

      int &local_idx = *local_nbdof_ptr[ss];
      miit = dofs_ptr[ss]->get<0>().begin();
      for(;miit!=hi_miit;miit++) {

        // Only set global idx for dofs on this processor part
        if(!(miit->get()->getActive())) continue;
        const BitRefLevel &dof_bit_level = (*miit)->getBitRefLevel();
        if((dof_bit_level&problem_bit_level)!=dof_bit_level) {
          continue;
        }

        dofs_array->push_back(NumeredDofEntity(*miit));
        dofs_shared_array.push_back(
          boost::shared_ptr<NumeredDofEntity>(dofs_array,&dofs_array->back())
        );

        int owner_proc = dofs_array->back().getOwnerProc();
        if(owner_proc<0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }

        if(owner_proc!=m_field.getCommRank()) {
          dofs_array->back().pArt = owner_proc;
          dofs_array->back().dofIdx = -1;
          dofs_array->back().petscGloablDofIdx = -1;
          dofs_array->back().petscLocalDofIdx = -1;
        } else {

          // Set part and indexes
          int glob_idx = start_ranges[ss]+local_idx;
          dofs_array->back().pArt = owner_proc;
          dofs_array->back().dofIdx = glob_idx;
          dofs_array->back().petscGloablDofIdx = glob_idx;
          dofs_array->back().petscLocalDofIdx = local_idx;
          local_idx++;

          unsigned char pstatus = dofs_array->back().getPStatus();
          // check id dof is shared, if that is a case global idx added to data structure and is
          // sended to other processors
          if(pstatus>0) {

            for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != dofs_array->back().getSharingProcsPtr()[proc]; proc++) {
              // make it different for rows and columns when needed
              if(ss == 0) {
                ids_data_packed_rows[dofs_array->back().getSharingProcsPtr()[proc]].push_back(
                  IdxDataType(dofs_array->back().getGlobalUniqueId(),glob_idx)
                );
              } else {
                ids_data_packed_cols[dofs_array->back().getSharingProcsPtr()[proc]].push_back(
                  IdxDataType(dofs_array->back().getGlobalUniqueId(),glob_idx)
                );
              }
              if(!(pstatus&PSTATUS_MULTISHARED)) {
                break;
              }
            }
          }

        }
      }

      numered_dofs_ptr[ss]->insert(dofs_shared_array.begin(),dofs_shared_array.end());

    }
    if(square_matrix) {
      local_nbdof_ptr[1] = local_nbdof_ptr[0];
    }

    // // If columns and rows have the same dofs indexing should be the same on both of them
    // {
    //   NumeredDofEntity_multiIndex::iterator dit_row = numered_dofs_ptr[0]->begin();
    //   NumeredDofEntity_multiIndex::iterator hi_dit_row = numered_dofs_ptr[0]->end();
    //   NumeredDofEntity_multiIndex::iterator dit_col = numered_dofs_ptr[1]->begin();
    //   NumeredDofEntity_multiIndex::iterator hi_dit_col = numered_dofs_ptr[1]->end();
    //   for(;dit_row!=hi_dit_row;dit_row++,dit_col++) {
    //     if(dit_row->get()->getPetscGlobalDofIdx()!=dit_col->get()->getPetscGlobalDofIdx()) {
    //       cerr << **dit_row << endl;
    //       cerr << **dit_col << endl;
    //     }
    //     if(dit_row->get()->getGlobalUniqueId()!=dit_col->get()->getGlobalUniqueId()) {
    //       cerr << "C " << **dit_row << endl;
    //       cerr << "R " << **dit_col << endl;
    //     }
    //
    //   }
    // }


    int nsends_rows = 0,nsends_cols = 0;
    // Non zero lengths[i] represent a message to i of length lengths[i].
    std::vector<int> lengths_rows(m_field.getCommSize()),lengths_cols(m_field.getCommSize());
    lengths_rows.clear();
    lengths_cols.clear();
    for(int proc = 0;proc<m_field.getCommSize();proc++) {
      lengths_rows[proc] = ids_data_packed_rows[proc].size()*data_block_size;
      lengths_cols[proc] = ids_data_packed_cols[proc].size()*data_block_size;
      if(!ids_data_packed_rows[proc].empty()) nsends_rows++;
      if(!ids_data_packed_cols[proc].empty()) nsends_cols++;
    }

    MPI_Status *status;
    ierr = PetscMalloc1(m_field.getCommSize(),&status);CHKERRQ(ierr);

    // Do rows
    int nrecvs_rows;	// number of messages received
    int *onodes_rows;	// list of node-ids from which messages are expected
    int *olengths_rows;	// corresponding message lengths
    int **rbuf_row;	// must bee freed by user

    {

      // make sure it is a PETSc comm
      ierr = PetscCommDuplicate(m_field.get_comm(),&m_field.get_comm(),NULL); CHKERRQ(ierr);

      //rows

      // Computes the number of messages a node expects to receive
      ierr = PetscGatherNumberOfMessages(m_field.get_comm(),NULL,&lengths_rows[0],&nrecvs_rows); CHKERRQ(ierr);
      //std::cerr << nrecvs_rows << std::endl;

      // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
      ierr = PetscGatherMessageLengths(
        m_field.get_comm(),
        nsends_rows,nrecvs_rows,
        &lengths_rows[0],&onodes_rows,&olengths_rows
      );  CHKERRQ(ierr);

      // Gets a unique new tag from a PETSc communicator. All processors that share
      // the communicator MUST call this routine EXACTLY the same number of times.
      // This tag should only be used with the current objects communicator; do NOT
      // use it with any other MPI communicator.
      int tag_row;
      ierr = PetscCommGetNewTag(m_field.get_comm(),&tag_row); CHKERRQ(ierr);

      // Allocate a buffer sufficient to hold messages of size specified in
      // olengths. And post Irecvs on these buffers using node info from onodes
      MPI_Request *r_waits_row; // must bee freed by user
      // rbuf has a pointers to messeges. It has size of of nrecvs (number of
      // messages) +1. In the first index a block is allocated,
      // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
      ierr = PetscPostIrecvInt(
        m_field.get_comm(),tag_row,nrecvs_rows,onodes_rows,olengths_rows,&rbuf_row,&r_waits_row
      ); CHKERRQ(ierr);
      ierr = PetscFree(onodes_rows); CHKERRQ(ierr);


      MPI_Request *s_waits_row; // status of sens messages
      ierr = PetscMalloc1(nsends_rows,&s_waits_row);CHKERRQ(ierr);

      // Send messeges
      for(int proc=0,kk=0; proc<m_field.getCommSize(); proc++) {
        if(!lengths_rows[proc]) continue; 	// no message to send to this proc
        ierr = MPI_Isend(
          &(ids_data_packed_rows[proc])[0], // buffer to send
          lengths_rows[proc], 		// message length
          MPIU_INT,proc, 			// to proc
          tag_row,m_field.get_comm(),s_waits_row+kk
        ); CHKERRQ(ierr);
        kk++;
      }

      if(nrecvs_rows) {
        ierr = MPI_Waitall(nrecvs_rows,r_waits_row,status);CHKERRQ(ierr);
      }
      if(nsends_rows) {
        ierr = MPI_Waitall(nsends_rows,s_waits_row,status);CHKERRQ(ierr);
      }

      ierr = PetscFree(r_waits_row); CHKERRQ(ierr);
      ierr = PetscFree(s_waits_row); CHKERRQ(ierr);


    }

    // cols
    int nrecvs_cols = nrecvs_rows;
    int *olengths_cols = olengths_rows;
    PetscInt **rbuf_col = rbuf_row;
    if(!square_matrix) {

      // Computes the number of messages a node expects to receive
      ierr = PetscGatherNumberOfMessages(m_field.get_comm(),NULL,&lengths_cols[0],&nrecvs_cols); CHKERRQ(ierr);

      // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
      int *onodes_cols;
      ierr = PetscGatherMessageLengths(m_field.get_comm(),
        nsends_cols,nrecvs_cols,
        &lengths_cols[0],&onodes_cols,&olengths_cols
      );  CHKERRQ(ierr);

      // Gets a unique new tag from a PETSc communicator.
      int tag_col;
      ierr = PetscCommGetNewTag(m_field.get_comm(),&tag_col); CHKERRQ(ierr);

      // Allocate a buffer sufficient to hold messages of size specified in
      // olengths. And post Irecvs on these buffers using node info from onodes
      MPI_Request *r_waits_col; // must bee freed by user
      ierr = PetscPostIrecvInt(m_field.get_comm(),tag_col,nrecvs_cols,onodes_cols,olengths_cols,&rbuf_col,&r_waits_col); CHKERRQ(ierr);
      ierr = PetscFree(onodes_cols); CHKERRQ(ierr);

      MPI_Request *s_waits_col; // status of sens messages
      ierr = PetscMalloc1(nsends_cols,&s_waits_col);CHKERRQ(ierr);

      // Send messeges
      for(int proc=0,kk=0; proc<m_field.getCommSize(); proc++) {
        if(!lengths_cols[proc]) continue; 	// no message to send to this proc
        ierr = MPI_Isend(
          &(ids_data_packed_cols[proc])[0],	// buffer to send
          lengths_cols[proc], 			// message length
          MPIU_INT,proc, 				// to proc
          tag_col,m_field.get_comm(),s_waits_col+kk); CHKERRQ(ierr);
          kk++;
        }

      if(nrecvs_cols) {
        ierr = MPI_Waitall(nrecvs_cols,r_waits_col,status);CHKERRQ(ierr);
      }
      if(nsends_cols) {
        ierr = MPI_Waitall(nsends_cols,s_waits_col,status);CHKERRQ(ierr);
      }

      ierr = PetscFree(r_waits_col); CHKERRQ(ierr);
      ierr = PetscFree(s_waits_col); CHKERRQ(ierr);

    }

    ierr = PetscFree(status); CHKERRQ(ierr);

    // set values received from other processors
    for(int ss = 0;ss<loop_size;ss++) {

      int nrecvs;
      int *olengths;
      int **data_procs;
      if(ss == 0) {
        nrecvs = nrecvs_rows;
        olengths = olengths_rows;
        data_procs = rbuf_row;
      } else {
        nrecvs = nrecvs_cols;
        olengths = olengths_cols;
        data_procs = rbuf_col;
      }

      IdxDataType *idx_data;
      GlobalUId uid;

      NumeredDofEntity_multiIndex::iterator dit;
      for(int kk=0; kk<nrecvs; kk++) {
        int len = olengths[kk];
        int *data_from_proc = data_procs[kk];
        for(int dd = 0;dd<len;dd+=data_block_size) {
          idx_data = (IdxDataType*)(&data_from_proc[dd]);
          bcopy(idx_data->uId,&uid,sizeof(UId));
          dit = numered_dofs_ptr[ss]->find(uid);
          if(dit == numered_dofs_ptr[ss]->end()) {
            // Dof is shared to this processor, however there is no element which
            // have this dof
            // continue;
            // DofEntity_multiIndex::iterator ddit = dofsField.find(uid);
            // if(ddit!=dofsField.end()) {
            //   std::cerr << **ddit << std::endl;
            // } else {
            //   std::ostringstream zz;
            //   zz << uid << std::endl;
            //   SETERRQ1(
            //     PETSC_COMM_SELF,
            //     MOFEM_OPERATION_UNSUCCESSFUL,
            //     "no such dof %s in mofem database",zz.str().c_str()
            //   );
            // }
            std::ostringstream zz;
            zz << uid << std::endl;
            SETERRQ1(
              PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"dof %s not found",zz.str().c_str()
            );
          }
          int global_idx = idx_data->globalDof;
          if(global_idx<0) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"received negative dof");
          }
          bool success;
          success = numered_dofs_ptr[ss]->modify(dit,NumeredDofEntity_mofem_index_change(global_idx));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          success = numered_dofs_ptr[ss]->modify(dit,NumeredDofEntity_part_change((*dit)->getPart(),global_idx));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
      }

    }

    if(square_matrix) {
      *(problem_ptr->tag_nbdof_data_col) = *(problem_ptr->tag_nbdof_data_row);
      *(problem_ptr->tag_local_nbdof_data_col) = *(problem_ptr->tag_local_nbdof_data_row);
    }

    ierr = PetscFree(olengths_rows); CHKERRQ(ierr);
    ierr = PetscFree(rbuf_row[0]); CHKERRQ(ierr);
    ierr = PetscFree(rbuf_row); CHKERRQ(ierr);
    if(!square_matrix) {
      ierr = PetscFree(olengths_cols); CHKERRQ(ierr);
      ierr = PetscFree(rbuf_col[0]); CHKERRQ(ierr);
      ierr = PetscFree(rbuf_col); CHKERRQ(ierr);
    }

    if(square_matrix) {
      if(numered_dofs_ptr[0]->size()!=numered_dofs_ptr[1]->size()) {
        SETERRQ2(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency for square_matrix %d!=%d",
          numered_dofs_ptr[0]->size(),numered_dofs_ptr[1]->size()
        );
      }
      if(problem_ptr->numered_dofs_rows!=problem_ptr->numered_dofs_cols) {
        SETERRQ(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency for square_matrix"
        );
      }
    }

    ierr = printPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);
    ierr = debugPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);

    PetscLogEventEnd(USER_EVENT_ProblemsManager,0,0,0,0);

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::buildSubProblem(
    const std::string &out_name,
    const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col,
    const std::string &main_problem,
    const bool square_matrix,
    int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;

    ierr = m_field.clear_problem(out_name); CHKERRQ(ierr);
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);

    // get reference to all problems
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblemByName;
    MoFEMProblemByName &problems_by_name =
    const_cast<MoFEMProblemByName&>(problems_ptr->get<Problem_mi_tag>());

    // get iterators to out problem, i.e. build problem
    MoFEMProblemByName::iterator out_problem_it = problems_by_name.find(out_name);
    if(out_problem_it==problems_by_name.end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "problem with name < %s > not defined (top tip check spelling)",
        out_name.c_str());
    }
    // get iterator to main problem, i.e. out problem is sub problem of main problem
    MoFEMProblemByName::iterator main_problem_it = problems_by_name.find(main_problem);
    if(main_problem_it==problems_by_name.end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "problem with name < %s > not defined (top tip check spelling)",
        main_problem.c_str());
    }

    // get dofs for row & columns for out problem,
    boost::shared_ptr<NumeredDofEntity_multiIndex> out_problem_dofs[] = {
      out_problem_it->numered_dofs_rows,
      out_problem_it->numered_dofs_cols
    };
    // get dofs for row & columns for main problem
    boost::shared_ptr<NumeredDofEntity_multiIndex> main_problem_dofs[] = {
      main_problem_it->numered_dofs_rows,
      main_problem_it->numered_dofs_cols
    };
    // get local indices counter
    int* nb_local_dofs[] = {
      out_problem_it->tag_local_nbdof_data_row,
      out_problem_it->tag_local_nbdof_data_col
    };
    // get global indices counter
    int* nb_dofs[] = {
      out_problem_it->tag_nbdof_data_row,
      out_problem_it->tag_nbdof_data_col
    };

    // set number of ghost nodes to zero
    {
      *out_problem_it->tag_ghost_nbdof_data_row = 0;
      *out_problem_it->tag_ghost_nbdof_data_col = 0;
    }

    // put rows & columns field names in array
    std::vector<std::string> fields[] = { fields_row,fields_col };

    // make data structure fos sub-problem data
    out_problem_it->subProblemData = boost::make_shared<MoFEMProblem::SubProblemData>();

    // use to keep shared_ptr
    std::vector<boost::shared_ptr<NumeredDofEntity> > dofs_shared_array;

    // Loop over rows and columns
    for(int ss = 0;ss!=(square_matrix ? 1 : 2);ss++) {

      // reset dofs and colimns counters
      (*nb_local_dofs[ss]) = 0;
      (*nb_dofs[ss]) = 0;
      // clear arrays
      out_problem_dofs[ss]->clear();

      // If DOFs are cleared clear finite elements too.
      out_problem_it->numeredFiniteElements.clear();

      // get dofs by field name and insert them in out problem multi-indices
      for(
        std::vector<std::string>::iterator fit = fields[ss].begin();
        fit!=fields[ss].end();fit++
      ) {
        NumeredDofEntityByFieldName::iterator dit = main_problem_dofs[ss]->get<FieldName_mi_tag>().lower_bound(*fit);
        NumeredDofEntityByFieldName::iterator hi_dit = main_problem_dofs[ss]->get<FieldName_mi_tag>().upper_bound(*fit);

        // Following reserve memory in sequences, only two allocations are here,
        // once for array of objects, next for array of shared pointers

        // reserve memory for field  dofs
        boost::shared_ptr<std::vector<NumeredDofEntity> > dofs_array =
        boost::shared_ptr<std::vector<NumeredDofEntity> >(new std::vector<NumeredDofEntity>());

        if(ss == 0) {
          out_problem_it->getRowDofsSeqence() = dofs_array;
        } else {
          out_problem_it->getColDofsSeqence() = dofs_array;
        }
        dofs_array->reserve(std::distance(dit,hi_dit));

        // create elements objects
        for(;dit!=hi_dit;dit++) {
          dofs_array->push_back(
            NumeredDofEntity(
              dit->get()->getDofEntityPtr(),
              dit->get()->getPetscGlobalDofIdx(),
              dit->get()->getPetscGlobalDofIdx(),
              dit->get()->getPetscLocalDofIdx(),
              dit->get()->getPart()
            )
          );
        }

        // reserve memory for shared pointers now
        dofs_shared_array.clear();
        dofs_shared_array.reserve(dofs_array->size());
        std::vector<NumeredDofEntity>::iterator vit = dofs_array->begin();
        for( ;vit!=dofs_array->end(); vit++ ) {
          dofs_shared_array.push_back(
            boost::shared_ptr<NumeredDofEntity>(dofs_array,&*vit)
          );
        }
        // fill multi-index
        out_problem_dofs[ss]->insert(
          dofs_shared_array.begin(),dofs_shared_array.end()
        );

      }
      // Set local indexes
      {
        NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type::iterator dit,hi_dit;
        dit = out_problem_dofs[ss]->get<Idx_mi_tag>().begin();
        hi_dit = out_problem_dofs[ss]->get<Idx_mi_tag>().end();
        for(;dit!=hi_dit;dit++) {
          int idx = -1; // if dof is not part of partition, set local index to -1
          if(dit->get()->getPart()==m_field.getCommRank()) {
            idx = (*nb_local_dofs[ss])++;
          }
          bool success = out_problem_dofs[ss]->modify(
            out_problem_dofs[ss]->project<0>(dit),NumeredDofEntity_local_idx_change(idx)
          );
        };
      }
      // Set global indexes, compress global indices
      {
        NumeredDofEntityByLocalIdx::iterator dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        NumeredDofEntityByLocalIdx::iterator hi_dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().end();
        const int nb = distance(dit,hi_dit);
        // get main problem global indices
        std::vector<int> main_indices(nb);
        for(std::vector<int>::iterator it = main_indices.begin();dit!=hi_dit;dit++,it++) {
          *it = dit->get()->getPetscGlobalDofIdx();
        }
        // crate is with global dofs
        IS is;
        ierr = ISCreateGeneral(
          m_field.get_comm(),nb,&*main_indices.begin(),PETSC_USE_POINTER,&is
        ); CHKERRQ(ierr);
        // create map form main problem global indices to out problem global indices
        AO ao;
        ierr = AOCreateMappingIS(is,PETSC_NULL,&ao); CHKERRQ(ierr);
        if(ss == 0) {
          ierr = ISDuplicate(is,&(out_problem_it->getSubData()->rowIs));
          // ierr = ISSort(out_problem_it->getSubData()->rowIs); CHKERRQ(ierr);
          out_problem_it->getSubData()->rowMap = ao;
          ierr = PetscObjectReference((PetscObject)ao); CHKERRQ(ierr);
        } else {
          ierr = ISDuplicate(is,&(out_problem_it->getSubData()->colIs));
          // ierr = ISSort(out_problem_it->getSubData()->colIs); CHKERRQ(ierr);
          out_problem_it->getSubData()->colMap = ao;
          ierr = PetscObjectReference((PetscObject)ao); CHKERRQ(ierr);
        }
        ierr = AOApplicationToPetscIS(ao,is); CHKERRQ(ierr);
        // set global number of DOFs
        ierr = ISGetSize(is,nb_dofs[ss]); CHKERRQ(ierr);
        ierr = ISDestroy(&is); CHKERRQ(ierr);
        // set out problem global indices after applying map
        dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(0);
        for(std::vector<int>::iterator it = main_indices.begin();dit!=hi_dit;dit++,it++) {
          bool success = out_problem_dofs[ss]->modify(
            out_problem_dofs[ss]->project<0>(dit),
            NumeredDofEntity_mofem_part_and_all_index_change(
              dit->get()->getPart(),*it,*it,dit->get()->getPetscLocalDofIdx()
            )
          );
        }
        // set global indices to nodes not on this part
        {
          NumeredDofEntityByLocalIdx::iterator dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(-1);
          NumeredDofEntityByLocalIdx::iterator hi_dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().upper_bound(-1);
          const int nb = distance(dit,hi_dit);
          std::vector<int> main_indices_non_local(nb);
          for(std::vector<int>::iterator it = main_indices_non_local.begin();dit!=hi_dit;dit++,it++) {
            *it = dit->get()->getPetscGlobalDofIdx();
          }
          IS is;
          ierr = ISCreateGeneral(
            m_field.get_comm(),nb,&*main_indices_non_local.begin(),PETSC_USE_POINTER,&is
          ); CHKERRQ(ierr);
          ierr = AOApplicationToPetscIS(ao,is); CHKERRQ(ierr);
          ierr = ISDestroy(&is); CHKERRQ(ierr);
          dit = out_problem_dofs[ss]->get<PetscLocalIdx_mi_tag>().lower_bound(-1);
          for(std::vector<int>::iterator it = main_indices_non_local.begin();dit!=hi_dit;dit++,it++) {
            bool success = out_problem_dofs[ss]->modify(
              out_problem_dofs[ss]->project<0>(dit),
              NumeredDofEntity_mofem_part_and_all_index_change(
                dit->get()->getPart(),dit->get()->getDofIdx(),*it,dit->get()->getPetscLocalDofIdx()
              )
            );
          }
        }
        ierr = AODestroy(&ao); CHKERRQ(ierr);
      }
    }

    if(square_matrix) {
      out_problem_it->numered_dofs_cols = out_problem_it->numered_dofs_rows;
      *(out_problem_it->tag_local_nbdof_data_col) = *(out_problem_it->tag_local_nbdof_data_row);
      *(out_problem_it->tag_nbdof_data_col) = *(out_problem_it->tag_nbdof_data_row);
      out_problem_it->getSubData()->colIs = out_problem_it->getSubData()->rowIs;
      out_problem_it->getSubData()->colMap = out_problem_it->getSubData()->rowMap;
      ierr = PetscObjectReference((PetscObject)out_problem_it->getSubData()->rowIs); CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)out_problem_it->getSubData()->rowMap); CHKERRQ(ierr);
    }

    ierr = printPartitionedProblem(&*out_problem_it,verb); CHKERRQ(ierr);
    ierr = debugPartitionedProblem(&*out_problem_it,verb); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

   PetscErrorCode ProblemsManager::buildCompsedProblem(
    const std::string &out_name,
    const std::vector<std::string> add_row_problems,
    const std::vector<std::string> add_col_problems,
    const bool square_matrix,
    int verb
  ) {

    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;

    ierr = m_field.clear_problem(out_name); CHKERRQ(ierr);
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    // get reference to all problems
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblemByName;
    MoFEMProblemByName &problems_by_name =
    const_cast<MoFEMProblemByName&>(problems_ptr->get<Problem_mi_tag>());

    // Get iterators to out problem, i.e. build problem
    MoFEMProblemByName::iterator out_problem_it = problems_by_name.find(out_name);
    if(out_problem_it==problems_by_name.end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "problem with name < %s > not defined (top tip check spelling)",
        out_name.c_str()
      );
    }
    // Make data structure for composed-problem data
    out_problem_it->composedProblemsData = boost::make_shared<ComposedProblemsData>();
    boost::shared_ptr<ComposedProblemsData> cmp_prb_data = out_problem_it->getComposedProblemsData();

    const std::vector<std::string>* add_prb[] = {
      &add_row_problems,&add_col_problems
    };
    std::vector<const MoFEMProblem*>* add_prb_ptr[] = {
      &cmp_prb_data->rowProblemsAdd,&cmp_prb_data->colProblemsAdd
    };
    std::vector<IS>* add_prb_is[] = {&cmp_prb_data->rowIs,&cmp_prb_data->colIs};

    // Get DOFs for row & columns for out problem,
    boost::shared_ptr<NumeredDofEntity_multiIndex> out_problem_dofs[] = {
      out_problem_it->numered_dofs_rows,
      out_problem_it->numered_dofs_cols
    };
    // Get local indices counter
    int* nb_local_dofs[] = {
      out_problem_it->tag_local_nbdof_data_row,
      out_problem_it->tag_local_nbdof_data_col
    };
    // Get global indices counter
    int* nb_dofs[] = {
      out_problem_it->tag_nbdof_data_row,
      out_problem_it->tag_nbdof_data_col
    };

    // Set number of ghost nodes to zero
    {
      *out_problem_it->tag_nbdof_data_row = 0;
      *out_problem_it->tag_nbdof_data_col = 0;
      *out_problem_it->tag_local_nbdof_data_row = 0;
      *out_problem_it->tag_local_nbdof_data_col = 0;
      *out_problem_it->tag_ghost_nbdof_data_row = 0;
      *out_problem_it->tag_ghost_nbdof_data_col = 0;
    }
    int nb_dofs_reserve[] = {0,0};

    // cerr << "AAA " << ((square_matrix)?1:2) << "\n";
    for(int ss = 0;ss!=((square_matrix)?1:2);ss++) {
      // cerr << "SS " << ss << endl;
      // cerr << add_prb[ss]->size() << endl;
      // cerr << add_prb_ptr[ss]->size() << endl;
      add_prb_ptr[ss]->reserve(add_prb[ss]->size());
      add_prb_is[ss]->reserve(add_prb[ss]->size());
      for(
        std::vector<std::string>::const_iterator
        vit=add_prb[ss]->begin();vit!=add_prb[ss]->end();vit++
      ) {
        MoFEMProblemByName::iterator prb_it = problems_by_name.find(*vit);
        if(prb_it==problems_by_name.end()) {
          SETERRQ1(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "problem with name < %s > not defined (top tip check spelling)",
            vit->c_str()
          );
        }
        add_prb_ptr[ss]->push_back(&*prb_it);
        if(ss==0) {
          *nb_dofs[ss] += add_prb_ptr[ss]->back()->getNbDofsRow();
          *nb_local_dofs[ss] += add_prb_ptr[ss]->back()->getNbLocalDofsRow();
          nb_dofs_reserve[ss] += add_prb_ptr[ss]->back()->numered_dofs_rows->size();
        } else {
          *nb_dofs[ss] += add_prb_ptr[ss]->back()->getNbDofsCol();
          *nb_local_dofs[ss] += add_prb_ptr[ss]->back()->getNbLocalDofsRow();
          nb_dofs_reserve[ss] += add_prb_ptr[ss]->back()->numered_dofs_cols->size();
        }
      }
    }
    if(square_matrix) {
      add_prb_ptr[1]->reserve(add_prb_ptr[0]->size());
      add_prb_is[1]->reserve(add_prb_ptr[0]->size());
      out_problem_it->numered_dofs_cols = out_problem_it->numered_dofs_rows;
      *nb_dofs[1] = *nb_dofs[0];
      *nb_local_dofs[1] = *nb_local_dofs[0];
    }

    boost::shared_ptr<std::vector<NumeredDofEntity> > dofs_array[2];

    // Reserve memory
    for(int ss = 0;ss!=((square_matrix)?1:2);ss++) {
      dofs_array[ss] = boost::make_shared<std::vector<NumeredDofEntity> >();
      dofs_array[ss]->reserve(nb_dofs_reserve[ss]);
      if(ss==0) {
        out_problem_it->getRowDofsSeqence() = dofs_array[ss];
      }
      if(ss==0) {
        out_problem_it->getColDofsSeqence() = dofs_array[ss];
      }
    }

    // Push back DOFs
     for(int ss = 0;ss!=((square_matrix)?1:2);ss++) {
      NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit,hi_dit;
      int shift_glob = 0;
      int shift_loc = 0;
      for(int pp = 0;pp!=add_prb_ptr[ss]->size();pp++) {
        PetscInt *dofs_out_idx_ptr;
        int nb_local_dofs = (*add_prb_ptr[ss])[pp]->getNbLocalDofsRow();
        ierr = PetscMalloc(nb_local_dofs*sizeof(int),&dofs_out_idx_ptr); CHKERRQ(ierr);
        if(ss == 0) {
          dit = (*add_prb_ptr[ss])[pp]->numered_dofs_rows->get<PetscGlobalIdx_mi_tag>().begin();
          hi_dit = (*add_prb_ptr[ss])[pp]->numered_dofs_rows->get<PetscGlobalIdx_mi_tag>().end();
        } else {
          dit = (*add_prb_ptr[ss])[pp]->numered_dofs_cols->get<PetscGlobalIdx_mi_tag>().begin();
          hi_dit = (*add_prb_ptr[ss])[pp]->numered_dofs_cols->get<PetscGlobalIdx_mi_tag>().end();
        }
        int is_nb = 0;
        for(;dit!=hi_dit;dit++) {
          // cerr << **dit << endl;d
          const int glob_idx = shift_glob+dit->get()->getPetscGlobalDofIdx();
          const int loc_idx =
          (dit->get()->getPart()==m_field.getCommRank())?(shift_loc+dit->get()->getPetscLocalDofIdx()):-1;
          // if(m_field.getCommRank()==1) {
          //   cerr << dit->get()->getPart() << " " << loc_idx << endl;
          // }
          dofs_array[ss]->push_back(
            NumeredDofEntity(
              dit->get()->getDofEntityPtr(),glob_idx,glob_idx,loc_idx,dit->get()->getPart()
            )
          );
          if(dit->get()->getPart()==m_field.getCommRank()) {
            dofs_out_idx_ptr[is_nb++] = glob_idx;
          }
        }
        if(is_nb!=nb_local_dofs) {
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
        }
        IS is;
        ierr = ISCreateGeneral(
          m_field.get_comm(),is_nb,dofs_out_idx_ptr,PETSC_OWN_POINTER,&is
        ); CHKERRQ(ierr);
        // cerr << "Push " << ss << " " << pp << endl;
        (*add_prb_is[ss]).push_back(is);
        if(ss==0) {
          shift_glob+=(*add_prb_ptr[ss])[pp]->getNbDofsRow();
          shift_loc+=(*add_prb_ptr[ss])[pp]->getNbLocalDofsRow();
        } else {
          shift_glob+=(*add_prb_ptr[ss])[pp]->getNbDofsCol();
          shift_loc+=(*add_prb_ptr[ss])[pp]->getNbLocalDofsCol();
        }
        if(square_matrix) {
          (*add_prb_ptr[1]).push_back((*add_prb_ptr[0])[pp]);
          // cerr << "Push square " << ss << " " << pp << endl;
          (*add_prb_is[1]).push_back(is);
          ierr = PetscObjectReference((PetscObject)is); CHKERRQ(ierr);
        }
      }
    }

    if((*add_prb_is[1]).size()!=(*add_prb_is[0]).size()) {
      SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
    }

    // Insert DOFs to problem multi-index
    std::vector<boost::shared_ptr<NumeredDofEntity> > dofs_shared_array;
    for(int ss = 0;ss!=((square_matrix)?1:2);ss++) {
      dofs_shared_array.clear();
      dofs_shared_array.reserve(dofs_array[ss]->size());
      for(
        std::vector<NumeredDofEntity>::iterator
        dit = dofs_array[ss]->begin();dit!=dofs_array[ss]->end();dit++
      ) {
        dofs_shared_array.push_back(boost::shared_ptr<NumeredDofEntity>(dofs_array[ss],&*dit));
      }
      if(ss==0) {
        out_problem_it->numered_dofs_rows->insert(dofs_shared_array.begin(),dofs_shared_array.end());
      } else {
        out_problem_it->numered_dofs_cols->insert(dofs_shared_array.begin(),dofs_shared_array.end());
      }
    }

    // Compress DOFs
    for(int ss = 0;ss!=((square_matrix)?1:2);ss++) {
      // if(ss == 0 && !enumerate_row) continue;
      // if(ss == 1 && !enumerate_col) continue;
      *nb_local_dofs[ss] = 0;
      boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_ptr;
      NumeredDofEntity_multiIndex::index<PetscGlobalIdx_mi_tag>::type::iterator dit,hi_dit;
      if(ss == 0) {
        dofs_ptr = out_problem_it->numered_dofs_rows;
      } else {
        dofs_ptr = out_problem_it->numered_dofs_cols;
      }
      dit = dofs_ptr->get<PetscGlobalIdx_mi_tag>().begin();
      hi_dit = dofs_ptr->get<PetscGlobalIdx_mi_tag>().end();
      std::vector<int> idx;
      idx.reserve(std::distance(dit,hi_dit));
      for(;dit!=hi_dit;dit++) {
        if(dit->get()->getPart()==m_field.getCommRank()) {
          bool success =
          dofs_ptr->get<PetscGlobalIdx_mi_tag>().modify(
            dit,NumeredDofEntity_local_idx_change((*nb_local_dofs[ss])++)
          );
          if(!success) {
            SETERRQ(
              m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
            );
          }
          idx.push_back(dit->get()->getPetscGlobalDofIdx());
        }
      }
      IS is;
      ierr = ISCreateGeneral(
        m_field.get_comm(),idx.size(),&*idx.begin(),PETSC_USE_POINTER,&is
      ); CHKERRQ(ierr);
      AO ao;
      ierr = AOCreateMappingIS(is,PETSC_NULL,&ao); CHKERRQ(ierr);
      for(int pp = 0;pp!=(*add_prb_is[ss]).size();pp++) {
        ierr = AOApplicationToPetscIS(ao,(*add_prb_is[ss])[pp]); CHKERRQ(ierr);
      }
      // Set DOFs numeration
      {
        std::vector<int> idx_new;
        idx_new.reserve(dofs_ptr->size());
        for(
          NumeredDofEntity_multiIndex::iterator
          dit = dofs_ptr->begin();dit!=dofs_ptr->end();dit++
        ) {
          idx_new.push_back(dit->get()->getPetscGlobalDofIdx());
        }
        IS is_new;
        ierr = ISCreateGeneral(
          m_field.get_comm(),idx_new.size(),&*idx_new.begin(),PETSC_USE_POINTER,&is_new
        ); CHKERRQ(ierr);
        ierr = AOApplicationToPetscIS(ao,is_new); CHKERRQ(ierr);
        std::vector<int>::iterator vit = idx_new.begin();
        for(
          NumeredDofEntity_multiIndex::iterator
          dit = dofs_ptr->begin();dit!=dofs_ptr->end();dit++
        ) {
          bool success =
          dofs_ptr->modify(
            dit,NumeredDofEntity_part_change(dit->get()->getPart(),*(vit++))
          );
          if(!success) {
            SETERRQ(
              m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
            );
          }
        }
        ierr = ISDestroy(&is_new); CHKERRQ(ierr);
      }
      ierr = ISDestroy(&is); CHKERRQ(ierr);
      ierr = AODestroy(&ao); CHKERRQ(ierr);

    }

    ierr = printPartitionedProblem(&*out_problem_it,verb); CHKERRQ(ierr);
    ierr = debugPartitionedProblem(&*out_problem_it,verb); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::partitionSimpleProblem(const std::string &name,int verb) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;
    if(!(cOre.getBuildMoFEM()&Core::BUILD_FIELD)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
    if(!(cOre.getBuildMoFEM()&Core::BUILD_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
    if(!(cOre.getBuildMoFEM()&Core::BUILD_ADJ)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
    if(!(cOre.getBuildMoFEM()&Core::BUILD_PROBLEM)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
    if(verb>0) {
      PetscPrintf(m_field.get_comm(),"Simple partition problem %s\n",name.c_str());
    }
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    // find p_miit
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblemByName;
    MoFEMProblemByName &problems_set = const_cast<MoFEMProblemByName&>(problems_ptr->get<Problem_mi_tag>());
    MoFEMProblemByName::iterator p_miit = problems_set.find(name);
    if(p_miit==problems_set.end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"problem < %s > is not found (top tip: check spelling)",name.c_str());
    }
    typedef boost::multi_index::index<NumeredDofEntity_multiIndex,Idx_mi_tag>::type NumeredDofEntitysByIdx;
    NumeredDofEntitysByIdx &dofs_row_by_idx = p_miit->numered_dofs_rows->get<Idx_mi_tag>();
    NumeredDofEntitysByIdx &dofs_col_by_idx = p_miit->numered_dofs_cols->get<Idx_mi_tag>();
    boost::multi_index::index<NumeredDofEntity_multiIndex,Idx_mi_tag>::type::iterator miit_row,hi_miit_row;
    boost::multi_index::index<NumeredDofEntity_multiIndex,Idx_mi_tag>::type::iterator miit_col,hi_miit_col;
    DofIdx &nb_row_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_row);
    DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
    nb_row_local_dofs = 0;
    nb_row_ghost_dofs = 0;
    DofIdx &nb_col_local_dofs = *((DofIdx*)p_miit->tag_local_nbdof_data_col);
    DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
    nb_col_local_dofs = 0;
    nb_col_ghost_dofs = 0;

    bool square_matrix = false;
    if(p_miit->numered_dofs_rows==p_miit->numered_dofs_cols) {
      square_matrix = true;
    }

    //get row range of local indices
    PetscLayout layout_row;
    const int *ranges_row;

    DofIdx nb_dofs_row = dofs_row_by_idx.size();
    ierr = PetscLayoutCreate(m_field.get_comm(),&layout_row); CHKERRQ(ierr);
    ierr = PetscLayoutSetBlockSize(layout_row,1); CHKERRQ(ierr);
    ierr = PetscLayoutSetSize(layout_row,nb_dofs_row); CHKERRQ(ierr);
    ierr = PetscLayoutSetUp(layout_row); CHKERRQ(ierr);
    ierr = PetscLayoutGetRanges(layout_row,&ranges_row); CHKERRQ(ierr);
    //get col range of local indices
    PetscLayout layout_col;
    const int *ranges_col;
    if(!square_matrix) {
      DofIdx nb_dofs_col = dofs_col_by_idx.size();
      ierr = PetscLayoutCreate(m_field.get_comm(),&layout_col); CHKERRQ(ierr);
      ierr = PetscLayoutSetBlockSize(layout_col,1); CHKERRQ(ierr);
      ierr = PetscLayoutSetSize(layout_col,nb_dofs_col); CHKERRQ(ierr);
      ierr = PetscLayoutSetUp(layout_col); CHKERRQ(ierr);
      ierr = PetscLayoutGetRanges(layout_col,&ranges_col); CHKERRQ(ierr);
    }
    for(unsigned int part = 0;part<(unsigned int)m_field.getCommSize();part++) {
      miit_row = dofs_row_by_idx.lower_bound(ranges_row[part]);
      hi_miit_row = dofs_row_by_idx.lower_bound(ranges_row[part+1]);
      if(distance(miit_row,hi_miit_row) != ranges_row[part+1]-ranges_row[part]) {
        SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
          "data inconsistency, distance(miit_row,hi_miit_row) != rend - rstart (%d != %d - %d = %d) ",
          distance(miit_row,hi_miit_row),ranges_row[part+1],ranges_row[part],ranges_row[part+1]-ranges_row[part]
        );
      }
      // loop rows
      for(;miit_row!=hi_miit_row;miit_row++) {
        bool success = dofs_row_by_idx.modify(miit_row,NumeredDofEntity_part_change(part,(*miit_row)->dofIdx));
        if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        if(part == (unsigned int)m_field.getCommRank()) {
          success = dofs_row_by_idx.modify(miit_row,NumeredDofEntity_local_idx_change(nb_row_local_dofs++));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
      }
      if(!square_matrix) {
        miit_col = dofs_col_by_idx.lower_bound(ranges_col[part]);
        hi_miit_col = dofs_col_by_idx.lower_bound(ranges_col[part+1]);
        if(distance(miit_col,hi_miit_col) != ranges_col[part+1]-ranges_col[part]) {
          SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
            "data inconsistency, distance(miit_col,hi_miit_col) != rend - rstart (%d != %d - %d = %d) ",
            distance(miit_col,hi_miit_col),ranges_col[part+1],ranges_col[part],ranges_col[part+1]-ranges_col[part]
          );
        }
        // loop cols
        for(;miit_col!=hi_miit_col;miit_col++) {
          bool success = dofs_col_by_idx.modify(miit_col,NumeredDofEntity_part_change(part,(*miit_col)->dofIdx));
          if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          if(part == (unsigned int)m_field.getCommRank()) {
            success = dofs_col_by_idx.modify(miit_col,NumeredDofEntity_local_idx_change(nb_col_local_dofs++));
            if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
        }
      }
    }
    ierr = PetscLayoutDestroy(&layout_row); CHKERRQ(ierr);
    if(!square_matrix) {
      ierr = PetscLayoutDestroy(&layout_col); CHKERRQ(ierr);
    }
    if(square_matrix) {
      nb_col_local_dofs = nb_row_local_dofs;
      nb_col_ghost_dofs = nb_row_ghost_dofs;
    }
    ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
    cOre.getBuildMoFEM() |= Core::PARTITION_PROBLEM;
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::partitionProblem(const std::string &name,int verb) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;

    if(!(cOre.getBuildMoFEM()&(1<<0))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
    if(!(cOre.getBuildMoFEM()&(1<<1))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
    if(!(cOre.getBuildMoFEM()&(1<<2))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"entFEAdjacencies not build");
    if(!(cOre.getBuildMoFEM()&(1<<3))) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"pRoblems not build");

    if(verb>0) {
      PetscPrintf(m_field.get_comm(),"Partition problem %s\n",name.c_str());
    }

    typedef NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofEntitysByIdx;
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;

    // Find problem pointer by name
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    ProblemsByName &problems_set = const_cast<ProblemsByName&>(problems_ptr->get<Problem_mi_tag>());
    ProblemsByName::iterator p_miit = problems_set.find(name);
    if(p_miit==problems_set.end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_INVALID_DATA,
        "problem with name %s not defined (top tip check spelling)",
        name.c_str());
    }
    int nb_dofs_row = p_miit->getNbDofsRow();

    Mat Adj;
    ierr = m_field.MatCreateMPIAdj_with_Idx_mi_tag(name,&Adj,verb); CHKERRQ(ierr);

    int m,n;
    ierr = MatGetSize(Adj,&m,&n); CHKERRQ(ierr);
    if(verb>2) {
      MatView(Adj,PETSC_VIEWER_STDOUT_WORLD);
    }

    //partitioning
    MatPartitioning part;
    IS is;
    ierr = MatPartitioningCreate(m_field.get_comm(),&part); CHKERRQ(ierr);
    //#ifdef __APPLE__
    // ierr = PetscBarrier((PetscObject)Adj); CHKERRQ(ierr);
    //#endif
    ierr = MatPartitioningSetAdjacency(part,Adj); CHKERRQ(ierr);
    ierr = MatPartitioningSetFromOptions(part); CHKERRQ(ierr);
    ierr = MatPartitioningSetNParts(part,m_field.getCommSize()); CHKERRQ(ierr);
    //#ifdef __APPLE__
    // ierr = PetscBarrier((PetscObject)part); CHKERRQ(ierr);
    //#endif
    ierr = MatPartitioningApply(part,&is); CHKERRQ(ierr);
    if(verb>2) {
      ISView(is,PETSC_VIEWER_STDOUT_WORLD);
    }
    // #ifdef __APPLE__
    // ierr = PetscBarrier((PetscObject)is); CHKERRQ(ierr);
    // #endif

    //gather
    IS is_gather,is_num,is_gather_num;
    ierr = ISAllGather(is,&is_gather); CHKERRQ(ierr);
    ierr = ISPartitioningToNumbering(is,&is_num); CHKERRQ(ierr);
    ierr = ISAllGather(is_num,&is_gather_num); CHKERRQ(ierr);
    const int *part_number,*petsc_idx;
    ierr = ISGetIndices(is_gather,&part_number);  CHKERRQ(ierr);
    ierr = ISGetIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
    int size_is_num,size_is_gather;
    ISGetSize(is_gather,&size_is_gather);
    if(size_is_gather != (int)nb_dofs_row) {
      SETERRQ2(
        PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,"data inconsistency %d != %d",
        size_is_gather,nb_dofs_row
      );
    }
    ISGetSize(is_num,&size_is_num);
    if(size_is_num != (int)nb_dofs_row) {
      SETERRQ2(
        PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ,
        "data inconsistency %d != %d",
        size_is_num,nb_dofs_row
      );
    }

    bool square_matrix = false;
    if(p_miit->numered_dofs_rows==p_miit->numered_dofs_cols) {
      square_matrix = true;
    }

    if(!square_matrix) {

      // FIXME: This is for back compatibility, if deprecate interface function
      // build interfaces is removed, this part of the code will be obsolete
      NumeredDofEntitysByIdx::iterator mit_row,hi_mit_row;
      mit_row = p_miit->numered_dofs_rows->get<Idx_mi_tag>().begin();
      hi_mit_row = p_miit->numered_dofs_rows->get<Idx_mi_tag>().end();
      NumeredDofEntitysByIdx::iterator mit_col,hi_mit_col;
      mit_col = p_miit->numered_dofs_cols->get<Idx_mi_tag>().begin();
      hi_mit_col = p_miit->numered_dofs_cols->get<Idx_mi_tag>().end();
      for(;mit_row!=hi_mit_row;mit_row++,mit_col++) {
        if(mit_col==hi_mit_col) {
          SETERRQ(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "check finite element definition, nb. of rows is not equal to number for columns"
          );
        }
        if(mit_row->get()->getGlobalUniqueId()!=mit_col->get()->getGlobalUniqueId()) {
          SETERRQ(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "check finite element definition, nb. of rows is not equal to number for columns"
          );
        }
      }

      // SETERRQ(
      //   PETSC_COMM_SELF,
      //   MOFEM_DATA_INCONSISTENCY,
      //   "check finite element definition, nb. of rows is not equal to number for columns"
      // );

    }

    if(verb>1) {
      PetscPrintf(m_field.get_comm(),"\tloop problem dofs");
    }

    try {

      // Set petsc global indices
      NumeredDofEntitysByIdx &dofs_row_by_idx_no_const =
      const_cast<NumeredDofEntitysByIdx&>(p_miit->numered_dofs_rows->get<Idx_mi_tag>());
      int &nb_row_local_dofs = *((int*)p_miit->tag_local_nbdof_data_row);
      int &nb_row_ghost_dofs = *((int*)p_miit->tag_ghost_nbdof_data_row);
      nb_row_local_dofs = 0;
      nb_row_ghost_dofs = 0;

      NumeredDofEntitysByIdx::iterator miit_dofs_row = dofs_row_by_idx_no_const.begin();
      for(;miit_dofs_row!=dofs_row_by_idx_no_const.end();miit_dofs_row++) {
        bool success = dofs_row_by_idx_no_const.modify(
          miit_dofs_row,NumeredDofEntity_part_change(
            part_number[(*miit_dofs_row)->dofIdx],
            petsc_idx[(*miit_dofs_row)->dofIdx]
          )
        );
        if(!success) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
        }
        if((*miit_dofs_row)->pArt == (unsigned int)m_field.getCommRank()) {
          success = dofs_row_by_idx_no_const.modify(
            miit_dofs_row,NumeredDofEntity_local_idx_change(nb_row_local_dofs++)
          );
          if(!success) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
        }
      }

      int &nb_col_local_dofs = *((int*)p_miit->tag_local_nbdof_data_col);
      int &nb_col_ghost_dofs = *((int*)p_miit->tag_ghost_nbdof_data_col);
      if(square_matrix) {
        nb_col_local_dofs = nb_row_local_dofs;
        nb_col_ghost_dofs = nb_row_ghost_dofs;
      } else {
        NumeredDofEntitysByIdx &dofs_col_by_idx_no_const =
        const_cast<NumeredDofEntitysByIdx&>(p_miit->numered_dofs_cols->get<Idx_mi_tag>());
        NumeredDofEntitysByIdx::iterator miit_dofs_col = dofs_col_by_idx_no_const.begin();
        nb_col_local_dofs = 0;
        nb_col_ghost_dofs = 0;
        for(;miit_dofs_col!=dofs_col_by_idx_no_const.end();miit_dofs_col++) {
          bool success = dofs_col_by_idx_no_const.modify(
            miit_dofs_col,NumeredDofEntity_part_change(
              part_number[(*miit_dofs_col)->dofIdx],
              petsc_idx[(*miit_dofs_col)->dofIdx]
            )
          );
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
          if((*miit_dofs_col)->pArt == (unsigned int)m_field.getCommRank()) {
            success = dofs_col_by_idx_no_const.modify(
              miit_dofs_col,NumeredDofEntity_local_idx_change(nb_col_local_dofs++)
            );
            if(!success) {
              SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
          }
        }
      }

    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    if(verb>1) {
      PetscPrintf(m_field.get_comm()," <- done\n");
    }

    ierr = ISRestoreIndices(is_gather,&part_number);  CHKERRQ(ierr);
    ierr = ISRestoreIndices(is_gather_num,&petsc_idx);  CHKERRQ(ierr);
    ierr = ISDestroy(&is_num); CHKERRQ(ierr);
    ierr = ISDestroy(&is_gather_num); CHKERRQ(ierr);
    ierr = ISDestroy(&is_gather); CHKERRQ(ierr);
    ierr = ISDestroy(&is); CHKERRQ(ierr);
    ierr = MatPartitioningDestroy(&part); CHKERRQ(ierr);
    ierr = MatDestroy(&Adj); CHKERRQ(ierr);
    ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);

    cOre.getBuildMoFEM() |= 1<<4;
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::inheretPartition(
    const std::string &name,
    const std::string &problem_for_rows,
    bool copy_rows,
    const std::string &problem_for_cols,
    bool copy_cols,
    int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;

    if(!(cOre.getBuildMoFEM()&Core::BUILD_PROBLEM))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"pRoblems not build");

    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type MoFEMProblemByName;

    //find p_miit
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    MoFEMProblemByName &problems_by_name = const_cast<MoFEMProblemByName&>(problems_ptr->get<Problem_mi_tag>());
    MoFEMProblemByName::iterator p_miit = problems_by_name.find(name);
    if(p_miit==problems_by_name.end()) {
      SETERRQ1(
        m_field.get_comm(),MOFEM_NOT_FOUND,
        "problem with name < %s > not defined (top tip check spelling)",
        name.c_str()
      );
    }
    if(verb>0) {
      PetscPrintf(
        m_field.get_comm(),
        "Compose problem %s from rows of %s and columns of %s\n",
        p_miit->getName().c_str(),
        problem_for_rows.c_str(),
        problem_for_cols.c_str()
      );
    }

    //find p_miit_row
    MoFEMProblemByName::iterator p_miit_row = problems_by_name.find(problem_for_rows);
    if(p_miit_row==problems_by_name.end()) {
      SETERRQ1(
        m_field.get_comm(),
        MOFEM_DATA_INCONSISTENCY,
        "problem with name < %s > not defined (top tip check spelling)",
        problem_for_rows.c_str()
      );
    }
    const boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_row = p_miit_row->numered_dofs_rows;

    //find p_mit_col
    MoFEMProblemByName::iterator p_miit_col = problems_by_name.find(problem_for_cols);
    if(p_miit_col==problems_by_name.end()) {
      SETERRQ1(
        m_field.get_comm(),
        MOFEM_DATA_INCONSISTENCY,
        "problem with name < %s > not defined (top tip check spelling)",
        problem_for_cols.c_str()
      );
    }
    boost::shared_ptr<NumeredDofEntity_multiIndex> dofs_col = p_miit_col->numered_dofs_cols;

    bool copy[] = { copy_rows, copy_cols };
    boost::shared_ptr<NumeredDofEntity_multiIndex> composed_dofs[] = {
      p_miit->numered_dofs_rows, p_miit->numered_dofs_cols
    };

    int* nb_local_dofs[] = { p_miit->tag_local_nbdof_data_row, p_miit->tag_local_nbdof_data_col };
    int* nb_dofs[] = { p_miit->tag_nbdof_data_row, p_miit->tag_nbdof_data_col };
    boost::shared_ptr<NumeredDofEntity_multiIndex> copied_dofs[] = { dofs_row, dofs_col };

    for(int ss = 0; ss<2;ss++) {

      // build indices
      *nb_local_dofs[ss] = 0;
      if(!copy[ss]) {

        // only copy indices which are belong to some elements if this problem
        std::vector<int> is_local,is_new;

        NumeredDofEntityByUId &dofs_by_uid = copied_dofs[ss]->get<Unique_mi_tag>();
        for(
          NumeredDofEntity_multiIndex::iterator
          dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++
        ) {
          NumeredDofEntityByUId::iterator diit = dofs_by_uid.find((*dit)->getGlobalUniqueId());
          if(diit==dofs_by_uid.end()) {
            SETERRQ(
              m_field.get_comm(),
              MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, could not find dof in composite problem"
            );
          }
          int part_number = (*diit)->getPart(); // get part number
          int petsc_global_dof = (*diit)->getPetscGlobalDofIdx();
          bool success;
          success = composed_dofs[ss]->modify(dit,NumeredDofEntity_part_change(part_number,petsc_global_dof));
          if(!success) {
            SETERRQ(
              m_field.get_comm(),
              MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful"
            );
          }
          if((*dit)->getPart() == (unsigned int)m_field.getCommRank()) {
            success =composed_dofs[ss]->modify(dit,NumeredDofEntity_local_idx_change((*nb_local_dofs[ss])++));
            if(!success) {
              SETERRQ(
                m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
              );
            }
            is_local.push_back(petsc_global_dof);
          }
        }

        AO ao;
        ierr = AOCreateMapping(m_field.get_comm(),is_local.size(),&is_local[0],NULL,&ao); CHKERRQ(ierr);

        // apply local to global mapping
        is_local.resize(0);
        for(
          NumeredDofEntity_multiIndex::iterator
          dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++
        ) {
          is_local.push_back((*dit)->getPetscGlobalDofIdx());
        }
        ierr = AOPetscToApplication(ao,is_local.size(),&is_local[0]); CHKERRQ(ierr);
        int idx2 = 0;
        for(
          NumeredDofEntity_multiIndex::iterator
          dit = composed_dofs[ss]->begin();dit!=composed_dofs[ss]->end();dit++
        ) {
          int part_number = (*dit)->getPart(); // get part number
          int petsc_global_dof = is_local[idx2++];
          bool success;
          success = composed_dofs[ss]->modify(dit,NumeredDofEntity_part_change(part_number,petsc_global_dof));
          if(!success) {
            SETERRQ(
              m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
            );
          }
        }

        ierr = AODestroy(&ao); CHKERRQ(ierr);

      } else {

        for(
          NumeredDofEntity_multiIndex::iterator
          dit = copied_dofs[ss]->begin();dit!=copied_dofs[ss]->end();dit++
        ) {
          std::pair<NumeredDofEntity_multiIndex::iterator,bool> p;
          p = composed_dofs[ss]->insert(
            boost::shared_ptr<NumeredDofEntity>(new NumeredDofEntity((*dit)->getDofEntityPtr()))
          );
          if(p.second) {
            (*nb_dofs[ss])++;
          }
          int dof_idx = (*dit)->getDofIdx();
          int part_number = (*dit)->getPart(); // get part number
          int petsc_global_dof = (*dit)->getPetscGlobalDofIdx();
          bool success;
          success = composed_dofs[ss]->modify(p.first,NumeredDofEntity_mofem_index_change(dof_idx));
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
          success = composed_dofs[ss]->modify(p.first,NumeredDofEntity_part_change(part_number,petsc_global_dof));
          if(!success) {
            SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
          if((*p.first)->getPart() == (unsigned int)m_field.getCommRank()) {
            success =composed_dofs[ss]->modify(p.first,NumeredDofEntity_local_idx_change((*nb_local_dofs[ss])++));
            if(!success) {
              SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
            }
          }
        }

      }
    }

    ierr = printPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);
    ierr = debugPartitionedProblem(&*p_miit,verb); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::printPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    if(verb>0) {
      std::ostringstream ss;
      ss << "partition_problem: rank = " << m_field.getCommRank() << " FEs row ghost dofs "<< *problem_ptr
      << " Nb. local dof " << problem_ptr->getNbLocalDofsRow() << " nb global row dofs " << problem_ptr->getNbDofsRow() << std::endl;
      ss << "partition_problem: rank = " << m_field.getCommRank() << " FEs col ghost dofs " << *problem_ptr
      << " Nb. local dof " << problem_ptr->getNbLocalDofsCol() << " nb global col dofs " << problem_ptr->getNbDofsCol() << std::endl;
      PetscSynchronizedPrintf(m_field.get_comm(),ss.str().c_str());
      PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
    }
    if(verb>1) {
      // FIXME mess if printed on more than one processors
      // std::ostringstream ss;
      std::cout << "rank = " << m_field.getCommRank() << " FEs row dofs "<< *problem_ptr << " Nb. row dof " << problem_ptr->getNbDofsRow()
      << " Nb. local dof " << problem_ptr->getNbLocalDofsRow() << std::endl;
      NumeredDofEntity_multiIndex::iterator miit_dd_row = problem_ptr->numered_dofs_rows->begin();
      for(;miit_dd_row!=problem_ptr->numered_dofs_rows->end();miit_dd_row++) {
        std::cout<<**miit_dd_row<<std::endl;
      }
      std::cout << "rank = " << m_field.getCommRank() << " FEs col dofs "<< *problem_ptr << " Nb. col dof " << problem_ptr->getNbDofsCol()
      << " Nb. local dof " << problem_ptr->getNbLocalDofsCol() << std::endl;
      NumeredDofEntity_multiIndex::iterator miit_dd_col = problem_ptr->numered_dofs_cols->begin();
      for(;miit_dd_col!=problem_ptr->numered_dofs_cols->end();miit_dd_col++) {
        std::cout<<**miit_dd_col<<std::endl;
      }
      // PetscSynchronizedPrintf(comm,ss.str().c_str());
      // PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::debugPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    if(debug>0) {

      typedef NumeredDofEntity_multiIndex::index<Idx_mi_tag>::type NumeredDofEntitysByIdx;
      NumeredDofEntitysByIdx::iterator dit,hi_dit;
      const NumeredDofEntitysByIdx* numered_dofs_ptr[] = {
        &(problem_ptr->numered_dofs_rows->get<Idx_mi_tag>()),
        &(problem_ptr->numered_dofs_cols->get<Idx_mi_tag>())
      };

      int* nbdof_ptr[] = {
        problem_ptr->tag_nbdof_data_row, problem_ptr->tag_nbdof_data_col
      };
      int* local_nbdof_ptr[] = {
        problem_ptr->tag_local_nbdof_data_row, problem_ptr->tag_local_nbdof_data_col
      };

      for(int ss = 0;ss<2;ss++) {

        dit = numered_dofs_ptr[ss]->begin();
        hi_dit = numered_dofs_ptr[ss]->end();
        for(;dit!=hi_dit;dit++) {
          if((*dit)->getPart()==(unsigned int)m_field.getCommRank()) {
            if((*dit)->getPetscLocalDofIdx()<0) {
              std::ostringstream zz;
              zz << "rank " << m_field.getCommRank() << " " << **dit;
              SETERRQ2(
                PETSC_COMM_SELF,
                MOFEM_IMPOSIBLE_CASE,
                "local dof index for %d (0-row, 1-col) not set, i.e. has negative value\n %s",
                ss,zz.str().c_str()
              );
            }
            if((*dit)->getPetscLocalDofIdx()>=*local_nbdof_ptr[ss]) {
              std::ostringstream zz;
              zz << "rank " << m_field.getCommRank() << " " << **dit;
              SETERRQ2(
                PETSC_COMM_SELF,
                MOFEM_IMPOSIBLE_CASE,
                "local dofs for %d (0-row, 1-col) out of range\n %s",
                ss,zz.str().c_str()
              );
            }
          } else {
            if((*dit)->getPetscGlobalDofIdx()<0) {
              std::ostringstream zz;
              zz << "rank " << m_field.getCommRank() << " " << **dit;
              SETERRQ2(
                PETSC_COMM_SELF,
                MOFEM_IMPOSIBLE_CASE,
                "global dof index for %d (0-row, 1-col) row not set, i.e. has negative value\n %s",
                ss,zz.str().c_str()
              );
            }
            if((*dit)->getPetscGlobalDofIdx()>=*nbdof_ptr[ss]) {
              std::ostringstream zz;
              zz << "rank " << m_field.getCommRank() << " " << **dit;
              SETERRQ2(
                PETSC_COMM_SELF,
                MOFEM_IMPOSIBLE_CASE,
                "global dofs for %d (0-row, 1-col) out of range\n %s",ss,zz.str().c_str()
              );
            }
          }
        }

      }

    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::partitionFiniteElements(
    const std::string &name,
    bool part_from_moab,
    int low_proc,
    int hi_proc,
    int verb
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    const EntFiniteElement_multiIndex *fe_ent_ptr;
    PetscFunctionBegin;

    if(!(cOre.getBuildMoFEM()&Core::BUILD_FIELD))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"fields not build");

    if(!(cOre.getBuildMoFEM()&Core::BUILD_FE))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"FEs not build");

    if(!(cOre.getBuildMoFEM()&Core::BUILD_ADJ))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
    if(!(cOre.getBuildMoFEM()&Core::BUILD_PROBLEM))

    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"problem not build");

    if(!(cOre.getBuildMoFEM()&Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"problem not partitioned");

    if(low_proc == -1) low_proc = m_field.getCommRank();
    if(hi_proc == -1) hi_proc = m_field.getCommRank();

    // Find pointer to problem of given name
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemByName;
    // Get pointer to problem data struture
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    ProblemByName &problems = const_cast<ProblemByName&>(problems_ptr->get<Problem_mi_tag>());
    ProblemByName::iterator p_miit = problems.find(name);
    if(p_miit == problems.end()) {
      SETERRQ1(
        m_field.get_comm(),
        MOFEM_NOT_FOUND,
        "problem < %s > not found (top tip: check spelling)",
        name.c_str()
      );
    }

    // Get reference on finite elements multi-index on the problem
    NumeredEntFiniteElement_multiIndex& problem_finite_elements = p_miit->numeredFiniteElements;

    // Clear all elements and data, build it again
    problem_finite_elements.clear();

    // Check if dofs and columns are the same, i.e. structurally symmetric problem
    bool do_cols_prob = true;
    if(p_miit->numered_dofs_rows == p_miit->numered_dofs_cols) {
      do_cols_prob = false;
    }

    // Loop over all elements in database and if right element is there add it
    // to problem finite element multi-index
    ierr = m_field.get_ents_finite_elements(&fe_ent_ptr); CHKERRQ(ierr);
    EntFiniteElement_multiIndex::iterator efit = fe_ent_ptr->begin();
    EntFiniteElement_multiIndex::iterator hi_efit = fe_ent_ptr->end();
    for(;efit!=hi_efit;efit++) {

      // if element is not part of problem
      if(((*efit)->getId()&p_miit->getBitFEId()).none()) continue;
      // if entity is not problem refinement level
      if(((*efit)->getBitRefLevel()&p_miit->getBitRefLevel())!=p_miit->getBitRefLevel()) continue;
      // create element
      boost::shared_ptr<NumeredEntFiniteElement> numered_fe(new NumeredEntFiniteElement(*efit));

      // check if rows and columns are the same on this element
      bool do_cols_fe = true;
      if(
        (numered_fe->sPtr->row_dof_view == numered_fe->sPtr->col_dof_view)
        && !do_cols_prob
      ) {
        do_cols_fe = false;
        numered_fe->cols_dofs = numered_fe->rows_dofs;
      } else {
        // different dofs on rows and columns
        numered_fe->cols_dofs = boost::shared_ptr<FENumeredDofEntity_multiIndex>(
          new FENumeredDofEntity_multiIndex()
        );
      }
      // get pointer to dofs multi-index on rows and columns
      boost::shared_ptr<FENumeredDofEntity_multiIndex> rows_dofs = numered_fe->rows_dofs;
      boost::shared_ptr<FENumeredDofEntity_multiIndex> cols_dofs = numered_fe->cols_dofs;
      // clear multi-indices
      rows_dofs->clear();
      if(do_cols_fe) {
        cols_dofs->clear();
      }
      NumeredDofEntity_multiIndex_uid_view_ordered rows_view;
      NumeredDofEntity_multiIndex_uid_view_ordered cols_view;

      // set partition to the element
      {
        if(part_from_moab) {
          // if partition is taken from moab partition
          int proc = (*efit)->getOwnerProc();
          NumeredEntFiniteElement_change_part(proc).operator()(numered_fe);
        } else {
          // count partition of the dofs in row, the larges dofs with given partition
          // is used to set partition of the element
          ierr = (*efit)->getRowDofView(
              *(p_miit->numered_dofs_rows),rows_view,moab::Interface::UNION
          ); CHKERRQ(ierr);
          std::vector<int> parts(m_field.getCommSize(),0);
          NumeredDofEntity_multiIndex_uid_view_ordered::iterator viit_rows;
          viit_rows = rows_view.begin();
          for(;viit_rows!=rows_view.end();viit_rows++) {
            parts[(*viit_rows)->pArt]++;
          }
          std::vector<int>::iterator pos = max_element(parts.begin(),parts.end());
          unsigned int max_part = distance(parts.begin(),pos);
          NumeredEntFiniteElement_change_part(max_part).operator()(numered_fe);
        }
      }

      // used to keep shared_ptr before inserting them to multi-index
      std::vector<boost::shared_ptr<FENumeredDofEntity> > dofs_shared_array;

      // set dofs on rows and columns (if are different)
      if(
        (numered_fe->getPart()>=(unsigned int)low_proc)&&
        (numered_fe->getPart()<=(unsigned int)hi_proc)
      ) {

        NumeredDofEntity_multiIndex_uid_view_ordered *dofs_view[] = {
          &rows_view, &cols_view
        };
        FENumeredDofEntity_multiIndex *fe_dofs[] = {
          rows_dofs.get(), cols_dofs.get()
        };

        for(int ss = 0;ss!=(do_cols_fe ? 2 : 1);ss++) {

          if(ss == 0) {
            if(part_from_moab) {
              // get row_view
              ierr = (*efit)->getRowDofView(
                *(p_miit->numered_dofs_rows),*dofs_view[ss],moab::Interface::UNION
              ); CHKERRQ(ierr);
            }
          } else {
            // get cols_views
            ierr = (*efit)->getColDofView(
              *(p_miit->numered_dofs_cols),*dofs_view[ss],moab::Interface::UNION
            ); CHKERRQ(ierr);
          }

          NumeredDofEntity_multiIndex_uid_view_ordered::iterator vit,hi_vit;
          vit = dofs_view[ss]->begin();
          hi_vit = dofs_view[ss]->end();

          // Following reserve memory in sequences, only two allocations are here,
          // once for array of objects, next for array of shared pointers

          // reserve memory for field  dofs
          boost::shared_ptr<std::vector<FENumeredDofEntity> > dofs_array =
          boost::make_shared<std::vector<FENumeredDofEntity> >();
          if(ss == 0) {
            numered_fe->getRowDofsSeqence() = dofs_array;
            if(!do_cols_fe) {
              numered_fe->getColDofsSeqence() = dofs_array;
            }
          } else {
            numered_fe->getColDofsSeqence() = dofs_array;
          }
          dofs_array->reserve(std::distance(vit,hi_vit));
          // reserve memory for shared pointers now
          dofs_shared_array.clear();
          dofs_shared_array.reserve(std::distance(vit,hi_vit));

          // create elements objects
          for(;vit!=hi_vit;vit++) {
            boost::shared_ptr<SideNumber> side_number_ptr;
            side_number_ptr = (*efit)->getSideNumberPtr((*vit)->getEnt());
            dofs_array->push_back(FENumeredDofEntity(side_number_ptr,*vit));
            dofs_shared_array.push_back(
              boost::shared_ptr<FENumeredDofEntity>(dofs_array,&dofs_array->back())
            );
          }

          // finally add DoFS to multi-indices
          fe_dofs[ss]->insert(dofs_shared_array.begin(),dofs_shared_array.end());

        }

      }
      std::pair<NumeredEntFiniteElement_multiIndex::iterator,bool> p;
      p = problem_finite_elements.insert(numered_fe);
      if(!p.second) {
        SETERRQ(m_field.get_comm(),MOFEM_NOT_FOUND,"element is there");
      }
      if(verb>1) {
        std::ostringstream ss;
        ss << *p_miit << std::endl;
        ss << *p.first << std::endl;
        typedef FENumeredDofEntityByUId FENumeredDofEntityByUId;
        FENumeredDofEntityByUId::iterator miit = (*p.first)->rows_dofs->get<Unique_mi_tag>().begin();
        for(;miit!= (*p.first)->rows_dofs->get<Unique_mi_tag>().end();miit++) ss << "rows: " << *(*miit) << std::endl;
        miit = (*p.first)->cols_dofs->get<Unique_mi_tag>().begin();
        for(;miit!=(*p.first)->cols_dofs->get<Unique_mi_tag>().end();miit++) ss << "cols: " << *(*miit) << std::endl;
        PetscSynchronizedPrintf(m_field.get_comm(),ss.str().c_str());
      }
    }
    if(verb>0) {
      typedef NumeredEntFiniteElement_multiIndex::index<Part_mi_tag>::type
      NumeredEntFiniteElementPart;
      NumeredEntFiniteElementPart::iterator miit,hi_miit;
      miit = problem_finite_elements.get<Part_mi_tag>().lower_bound(m_field.getCommRank());
      hi_miit = problem_finite_elements.get<Part_mi_tag>().upper_bound(m_field.getCommRank());
      int count = distance(miit,hi_miit);
      std::ostringstream ss;
      ss << *p_miit;
      ss << " Nb. elems " << count << " on proc " << m_field.getCommRank() << std::endl;
      PetscSynchronizedPrintf(m_field.get_comm(),ss.str().c_str());
      PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
    }

    cOre.getBuildMoFEM() |= Core::PARTITION_FE;
    PetscFunctionReturn(0);
  }

  PetscErrorCode ProblemsManager::partitionGhostDofs(const std::string &name,int verb) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    const MoFEMProblem_multiIndex *problems_ptr;
    PetscFunctionBegin;

    if(!(cOre.getBuildMoFEM()&Core::PARTITION_PROBLEM))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"partition of problem not build");
    if(!(cOre.getBuildMoFEM()&Core::PARTITION_FE))
    SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"partitions finite elements not build");

    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
    //find p_miit
    ierr = m_field.get_problems(&problems_ptr); CHKERRQ(ierr);
    ProblemsByName &problems_set = const_cast<ProblemsByName&>(problems_ptr->get<Problem_mi_tag>());
    ProblemsByName::iterator p_miit = problems_set.find(name);
    //
    DofIdx &nb_row_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_row);
    DofIdx &nb_col_ghost_dofs = *((DofIdx*)p_miit->tag_ghost_nbdof_data_col);
    //
    nb_row_ghost_dofs = 0;
    nb_col_ghost_dofs = 0;
    if(m_field.getCommSize()>1) {
      NumeredDofEntity_multiIndex_uid_view_ordered ghost_idx_col_view,ghost_idx_row_view;
      NumeredEntFiniteElement_multiIndex::index<Part_mi_tag>::type::iterator fe_it,hi_fe_it;
      fe_it = p_miit->numeredFiniteElements.get<Part_mi_tag>().lower_bound(m_field.getCommRank());
      hi_fe_it = p_miit->numeredFiniteElements.get<Part_mi_tag>().upper_bound(m_field.getCommRank());
      for(;fe_it!=hi_fe_it;fe_it++) {
        typedef FENumeredDofEntity_multiIndex::iterator dof_it;
        if((*fe_it)->rows_dofs->size()>0) {
          dof_it rowdofit,hi_rowdofit;
          rowdofit = (*fe_it)->rows_dofs->begin();
          hi_rowdofit = (*fe_it)->rows_dofs->end();
          for(;rowdofit!=hi_rowdofit;rowdofit++) {
            if((*rowdofit)->getPart()==(unsigned int)m_field.getCommRank()) continue;
            ghost_idx_row_view.insert((*rowdofit)->getNumeredDofEntityPtr());
          }
        }
        if((*fe_it)->cols_dofs->size()>0) {
          dof_it coldofit,hi_coldofit;
          coldofit = (*fe_it)->cols_dofs->begin();
          hi_coldofit = (*fe_it)->cols_dofs->end();
          for(;coldofit!=hi_coldofit;coldofit++) {
            if((*coldofit)->getPart()==(unsigned int)m_field.getCommRank()) continue;
            ghost_idx_col_view.insert((*coldofit)->getNumeredDofEntityPtr());
          }
        }
      }
      DofIdx *nb_ghost_dofs[2] = { &nb_col_ghost_dofs, &nb_row_ghost_dofs };
      DofIdx nb_local_dofs[2] = {
        *((DofIdx*)p_miit->tag_local_nbdof_data_col),
        *((DofIdx*)p_miit->tag_local_nbdof_data_row)
      };
      NumeredDofEntity_multiIndex_uid_view_ordered *ghost_idx_view[2] = { &ghost_idx_col_view, &ghost_idx_row_view };
      NumeredDofEntityByUId *dof_by_uid_no_const[2] = {
        &p_miit->numered_dofs_cols->get<Unique_mi_tag>(),
        &p_miit->numered_dofs_rows->get<Unique_mi_tag>()
      };
      int loop_size = 2;
      if(p_miit->numered_dofs_cols==p_miit->numered_dofs_rows) {
        loop_size = 1;
      }
      for(int ss = 0;ss<loop_size;ss++) {
        NumeredDofEntity_multiIndex_uid_view_ordered::iterator ghost_idx_miit = ghost_idx_view[ss]->begin();
        for(;ghost_idx_miit!=ghost_idx_view[ss]->end();ghost_idx_miit++) {
          NumeredDofEntityByUId::iterator diit = dof_by_uid_no_const[ss]->find((*ghost_idx_miit)->getGlobalUniqueId());
          if((*diit)->petscLocalDofIdx!=(DofIdx)-1) {
            SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"inconsistent data, ghost dof already set");
          }
          bool success = dof_by_uid_no_const[ss]->modify(diit,NumeredDofEntity_local_idx_change(nb_local_dofs[ss]++));
          if(!success) SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          (*nb_ghost_dofs[ss])++;
        }
      }
      if(loop_size==1) {
        (*nb_ghost_dofs[1]) = (*nb_ghost_dofs[0]);
      }
    }
    if(verb>0) {
      std::ostringstream ss;
      ss << "partition_ghost_col_dofs: rank = " << m_field.getCommRank()
      << " FEs col ghost dofs "<< *p_miit
      << " Nb. col ghost dof " << p_miit->getNbGhostDofsCol()
      << " Nb. local dof " << p_miit->getNbLocalDofsCol() << std::endl;
      ss << "partition_ghost_row_dofs: rank = " << m_field.getCommRank()
      << " FEs row ghost dofs "<< *p_miit
      << " Nb. row ghost dof " << p_miit->getNbGhostDofsRow()
      << " Nb. local dof " << p_miit->getNbLocalDofsRow() << std::endl;
      if(verb>1) {
        NumeredDofEntity_multiIndex::iterator miit_dd_col = p_miit->numered_dofs_cols->begin();
        for(;miit_dd_col!=p_miit->numered_dofs_cols->end();miit_dd_col++) {
          if((*miit_dd_col)->pArt==(unsigned int)m_field.getCommRank()) continue;
          if((*miit_dd_col)->petscLocalDofIdx==(DofIdx)-1) continue;
          ss<<*(*miit_dd_col)<<std::endl;
        }
      }
      PetscSynchronizedPrintf(m_field.get_comm(),ss.str().c_str());
      PetscSynchronizedFlush(m_field.get_comm(),PETSC_STDOUT);
    }

    cOre.getBuildMoFEM() |= Core::PARTITION_GHOST_DOFS;
    PetscFunctionReturn(0);
  }


}
