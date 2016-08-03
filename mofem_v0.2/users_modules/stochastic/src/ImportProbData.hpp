/* \brief ImportProbData.hpp
 *
 * Edited and modified by Xiaoyi Zhou.
 *
 * This is used to read probability information for basic random variables from
 * text file.
 *
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

struct ImportProbData {
  
  ublas::matrix<double> MargProb;    // Marginal probability distribution function
  ublas::matrix<double> CorrMat;     // Correlation matrix
  ublas::vector<double> MatStrength; // Material strength
  ublas::vector<double> PlyAngle;    // Ply angle of orientation
  vector<string> NameVars;           // Name of random variables
  vector<string> NameMatVars;       // Name of random variables for material properties
  string HomoMethod;                 // Homogenization method
  int NumVars;                       // Number of variables
  int NumMatVars =0;                 // Number of variables for material properties
  int NumLayers;                     // Number of layers
  int ExaminedLayer;                 // Examined layer
  int TimePoint;                     // Time point for selection of wt in degradation analysis
  int FailureCriterion;              // Type of failure criterion
  int AnalysisType;                  // Analysis type 10: FORM 20: SORM 30: MCS 31: MCIS
  double SearchStep;                 // Search step size
  
  
  //------------------------------------------------------------------------------
  // To count number of a specific character or a substring in a string;
  //
  
  virtual PetscErrorCode str_cnt(char *mstr,char substr,int &cnt) {
    PetscFunctionBegin;
    
    int j=0;
    for (int i=0;i<strlen(mstr);i++) {
      if (*(mstr+i)==substr) {
        j++;
      }
    }
    cnt=j;
    
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // To determine the position of a specific charater in a string
  //
  
  virtual PetscErrorCode str_pos(char *mstr,char substr,ublas::vector<int> &pos) {
    PetscFunctionBegin;
    int j=1;
    for (int i=0;i<strlen(mstr);i++) {
      if (*(mstr+i)==substr) {
        //cout<<*(mstr+i)<<'\t'<<i<<'\n';
        pos(j) = i;
        j++;
      }
    }
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // To read probability data from file
  //
  
  virtual PetscErrorCode ProbdataFileIn() {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    
    char prob_data_file_name[255];
    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetString(PETSC_NULL,"-my_prob_data_file",prob_data_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_prob_data_file (PROB_DATA FILE NEEDED)!!!");
    }
    
    ifstream ProbDataFile;
    //ProbDataFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Input_probdata.txt",ifstream::in);
    ProbDataFile.open(prob_data_file_name,ifstream::in);
    if (!ProbDataFile) {
      cout << "\n\nProbability data file does not exists!\n" << endl;
      exit(EXIT_FAILURE);
    }
    
    char   buffer[256];
    string stringbuf;
    string substringbuf;
    string datatype;
    //int    *pos;
    ublas::vector<int> pos;
    int    cnt;
    int    MAR_IX, COR_IX, STR_IX, ANG_IX;
    MAR_IX = 0; COR_IX = 0; STR_IX = 0; ANG_IX = 0;
    
    while (!ProbDataFile.eof()) {
      ProbDataFile.getline(buffer,200);
      if (strlen(buffer)>0) {
        if (isalpha(buffer[0])) {//(isdigit(buffer[0]) == 0) {
          stringbuf = (string)buffer;
          if (stringbuf.compare(0,3,"NUM") == 0) {
            // cout<<"Next line is data for number of variables"<<endl;
            datatype = "NUMBER";
          }
          else if (stringbuf.compare(0,4,"STEP") == 0) {
            // cout<<"Next line is data for size of search step"<<endl;
            datatype = "STEP";
          }
          else if (stringbuf.compare(0,4,"HOMO") == 0) {
            // cout<<"Next line is data for correlation matrix"<<endl;
            datatype = "HOMOGENIZATION";
          }
          else if (stringbuf.compare(0,3,"COR") == 0) {
            // cout<<"Next line is data for correlation matrix"<<endl;
            datatype = "CORRELATION";
          }
          else if (stringbuf.compare(0,3,"MAR") == 0) {
            // cout<<"Next line is data for marginal distribution"<<endl;
            datatype = "MARGINAL";
          }
          else if (stringbuf.compare(0,3,"STR") == 0) {
            // cout<<"Next lines are data for material strength"<<endl;
            datatype = "STRENGTH";
          }
          else if (stringbuf.compare(0,4,"NAME") == 0) {
            // cout<<"Next lines are data for names of variables"<<endl;
            datatype = "NAME";
          }
          else if (stringbuf.compare(0,6,"CHOICE") == 0) {
            // cout<<"Next lines are data for names of variables"<<endl;
            datatype = "FAILURE";
          }
          else if (stringbuf.compare(0,3,"LAY") == 0) {
           // cout<<"Next line is data for correlation matrix"<<endl;
           datatype = "LAYER";
           }
          else if (stringbuf.compare(0,3,"EXA") == 0) {
            // cout<<"Next line is data for correlation matrix"<<endl;
            datatype = "EXAMINED";
          }
          else if (stringbuf.compare(0,3,"ANG") == 0) {
           // cout<<"Next lines are data for angle of orietations"<<endl;
           datatype = "ANGLE";
           }
          else if (stringbuf.compare(0,4,"TIME") == 0) {
            // cout<<"Next lines are data for angle of orietations"<<endl;
            datatype = "TIME";
          }
          else if (stringbuf.compare(0,8,"ANALYSIS") == 0) {
            // cout<<"Next lines are data for angle of orietations"<<endl;
            datatype = "ANALYSIS";
          }
          stringbuf.clear();
        }
        else {
          ierr = str_cnt(buffer,',',cnt); CHKERRQ(ierr);
          //pos = new int[cnt];
          //cout<<"the number of specified character: \t"<<cnt<<endl;
          pos.resize(cnt+1);pos.clear();
          ierr = str_pos(buffer,',',pos); CHKERRQ(ierr);
          stringbuf = (string)buffer;//cout<<buffer<<endl;
          
          if (datatype.compare(0,3,"NUM") == 0) { // number of variables
            substringbuf = stringbuf.substr(0,pos(1));
            NumVars = atoi(substringbuf.c_str());
          }
          else if (datatype.compare(0,4,"STEP") == 0) { // number of variables
            substringbuf = stringbuf.substr(0,pos(1));
            SearchStep = atof(substringbuf.c_str());
          }
          else if (datatype.compare(0,3,"MAR") == 0) { // marginal distribution
            // Declaration
            if (MAR_IX ==0) {
              MargProb.resize(NumVars,10);
            }
            MAR_IX ++;
            
            // Insert data into probdata
            for (int i=1; i<=cnt;i++) {
              if (i==1) {
                substringbuf = stringbuf.substr(0,pos(i));
              }
              else {
                substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
              }
              MargProb(MAR_IX-1,i-1) = atof(substringbuf.c_str());
              substringbuf.clear();
            }
          }
          else if (datatype.compare(0,3,"COR") == 0) { // correlation matrix
            // Declaration
            if (COR_IX ==0) {
              CorrMat.resize(NumVars,NumVars);
            }
            COR_IX ++;
            
            // Insert data into probdata
            for (int i=1; i<=cnt;i++) {
              if (i == 1) {
                substringbuf = stringbuf.substr(0,pos(i));
              }
              else {
                substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
              }
              CorrMat(COR_IX-1, i-1) = atof(substringbuf.c_str());
              substringbuf.clear();
            }
          }
          else if (datatype.compare(0,3,"STR") == 0) { // material strength
            // Declaration
            if (STR_IX ==0) {
              MatStrength.resize(6);
            }
            STR_IX ++;
            
            // Insert data into probdata
            for (int i=1; i<=cnt;i++) {
              if (i == 1) {
                substringbuf = stringbuf.substr(0,pos(i));
              }
              else {
                substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
              }
              MatStrength(STR_IX-1) = atof(substringbuf.c_str());
              substringbuf.clear();
            }
          }
          else if (datatype.compare(0,4,"HOMO") == 0) {
            substringbuf = stringbuf.substr(pos(1)+1,pos(2)-pos(1)-1);
            cout<<"Homogenization method: "<<substringbuf<<endl;
            HomoMethod = substringbuf;
            substringbuf.clear();
          }
          else if (datatype.compare(0,4,"NAME") == 0) {
            // Insert data into probdata
            for (int i=1; i<=cnt;i++) {
              substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
              cout<<substringbuf<<endl;
              NameVars.push_back(substringbuf);
              // Identify random variables for material properties
              if (substringbuf.compare(0,2,"Em") == 0) {       // Young's modulus of matrix - isotropic
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,3,"NUm") == 0) { // Poisson's ratio of matrix - isotropic
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,3,"NUp") == 0) { // Poisson's ratio in p-direction of fibre
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,3,"NUz") == 0) { // Poisson's ratio in z-direction of fibre
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,2,"Ep") == 0) {  // Young's modulus in p-direction of fibre
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,2,"Ez") == 0) {  // Young's modulus in longitudinal direction
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              else if (substringbuf.compare(0,3,"Gzp") == 0) { // Shear modulus in z-direction of fibre
                NumMatVars ++;
                NameMatVars.push_back(substringbuf);
              }
              substringbuf.clear();
            }
          }
          else if (datatype.compare(0,7,"FAILURE") == 0) {
            substringbuf = stringbuf.substr(0,pos(1));
            FailureCriterion = atoi(substringbuf.c_str());
          }
          else if (datatype.compare(0,3,"LAY") == 0) { // number of layers
           substringbuf = stringbuf.substr(0,pos(1));
           NumLayers = atoi(substringbuf.c_str());
          }
          else if (datatype.compare(0,3,"EXA") == 0) { // number of layers
            substringbuf = stringbuf.substr(0,pos(1));
            ExaminedLayer = atoi(substringbuf.c_str());
          }
          else if (datatype.compare(0,8,"ANALYSIS") == 0) { // number of layers
            substringbuf = stringbuf.substr(0,pos(1));
            AnalysisType = atoi(substringbuf.c_str());
          }          else if (datatype.compare(0,4,"TIME") == 0) { // examined time point
            substringbuf = stringbuf.substr(0,pos(1));
            TimePoint = atoi(substringbuf.c_str());
          }
          else if (datatype.compare(0,5,"ANGLE") == 0) { // material strength
           // Declaration
           if (ANG_IX ==0) {
             PlyAngle.resize(NumLayers);
           }
           ANG_IX ++;
           
           // Insert data into probdata
           for (int i=1; i<=cnt;i++) {
             if (i == 1) {
               substringbuf = stringbuf.substr(0,pos(i));
             }
             else {
               substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
             }
             PlyAngle(ANG_IX-1) = atof(substringbuf.c_str());//cout<<"\n Layer angle \n"<<substringbuf<<endl;
             substringbuf.clear();
           }
          }
          // free dynamically allocated memory
          //delete [] pos;
          stringbuf.clear();
          //cout<<"size: "<<pos.size()<<datatype<<"\t"<<cnt<<"\t"<<buffer<<endl;
        }
      }
    }
    //delete [] pos;
    //cout<<"Marginal probability data:\t"<<MargProb<<endl;
    ProbDataFile.close();
    //cout<<"Number of names: "<<NameVars.size()<<endl;
    PetscFunctionReturn(0);
  }
  
};

