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
  ublas::matrix<double> MatStrength; // Material strength
  vector<string> NameVars;           // Name of random variables
  int NumVars;                       // Number of variables
  
  
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
    
    ErrorCode rval;
    PetscErrorCode ierr;
    
    char prob_data_file_name[255];
    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetString(PETSC_NULL,"-my_prob_data_file",prob_data_file_name,255,&flg); CHKERRQ(ierr);
    
    
    ifstream ProbDataFile;
    //ProbDataFile.open("//mnt//home//Dropbox//DURACOMP_Cal//009_MoFEM//04_ReliabilityAnalysis//Input_probdata.txt",ifstream::in);
    ProbDataFile.open(prob_data_file_name,ifstream::in);
    
    char   buffer[256];
    string stringbuf;
    string substringbuf;
    string datatype;
    //int    *pos;
    ublas::vector<int> pos;
    int    cnt;
    int    MAR_IX, COR_IX, STR_IX;
    MAR_IX = 0; COR_IX = 0; STR_IX = 0;
    
    while (!ProbDataFile.eof()) {
      ProbDataFile.getline(buffer,200);
      if (strlen(buffer)>0) {
        if (isdigit(buffer[0]) == 0) {
          stringbuf = (string)buffer;
          if (stringbuf.compare(0,3,"NUM") == 0) {
            // cout<<"Next line is data for number of variables"<<endl;
            datatype = "NUMBER";
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
          stringbuf.clear();
        }
        else {
          ierr = str_cnt(buffer,',',cnt); CHKERRQ(ierr);
          //pos = new int[cnt];
          cout<<"the number of specified character: \t"<<cnt<<endl;
          pos.resize(cnt+1);pos.clear();
          ierr = str_pos(buffer,',',pos); CHKERRQ(ierr);
          stringbuf = (string)buffer;//cout<<buffer<<endl;
          
          if (datatype.compare(0,3,"NUM") == 0) { // number of variables
            substringbuf = stringbuf.substr(0,pos(1));
            NumVars = atoi(substringbuf.c_str());
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
          else if (datatype.compare(0,3,"STR") == 0) { // material strenth
            // Declaration
            if (STR_IX ==0) {
              MatStrength.resize(6,cnt);
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
              MatStrength(STR_IX-1, i-1) = atof(substringbuf.c_str());
              substringbuf.clear();
            }
          }
          else if (datatype.compare(0,4,"NAME") == 0) {
            // Insert data into probdata
            for (int i=1; i<=cnt;i++) {
              substringbuf = stringbuf.substr(pos(i-1)+1,(pos(i)-pos(i-1))-1);
              //cout<<substringbuf<<endl;
              NameVars.push_back(substringbuf);
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
    cout<<"Marginal probability data:\t"<<MargProb<<endl;
    ProbDataFile.close();
    //cout<<"Number of names: "<<NameVars.size()<<endl;
    PetscFunctionReturn(0);
  }
  
};

