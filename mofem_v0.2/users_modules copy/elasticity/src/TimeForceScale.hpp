/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Edited and modified by Hassan.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
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

struct TimeForceScale: public MethodsForOp {
//Hassan: This function to read data file (once) and save it in a pair vector ts
   
    map<double,double> tSeries;
    int readFile,debug;

    TimeForceScale(): readFile(0),debug(1) {
      PetscErrorCode ierr;
      ierr = timeData(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    };

    ErrorCode rval;
    PetscErrorCode ierr;

    PetscErrorCode timeData() {
      PetscFunctionBegin;
      char time_file_name[255];
      PetscBool flg = PETSC_TRUE;
      ierr = PetscOptionsGetString(PETSC_NULL,"-my_time_data_file",time_file_name,255,&flg); CHKERRQ(ierr);
      if(flg != PETSC_TRUE) {
	SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_time_data_file (time_data FILE NEEDED)");
      }
      FILE *time_data = fopen(time_file_name,"r");
      if(time_data == NULL) {
	SETERRQ1(PETSC_COMM_SELF,1,"*** ERROR data file < %s > open unsucessfull",time_file_name);
      }
      double no1 = 0.0, no2 = 0.0;
      tSeries[no1] = no2;
      while(! feof (time_data)){
        int n = fscanf(time_data,"%lf %lf",&no1,&no2);
        if((n <= 0)||((no1==0)&&(no2==0))) {
          fgetc(time_data);
	  continue;
	}
	if(n != 2){
	  SETERRQ(PETSC_COMM_SELF,1,"*** ERROR read data file error (check input time data file)");
	}
	tSeries[no1] = no2;
      }
      int r = fclose(time_data);      
      if(debug) {
	map<double, double>::iterator tit = tSeries.begin();
	for(;tit!=tSeries.end();tit++) {
	  PetscPrintf(PETSC_COMM_WORLD,"** read time series %3.2e time %3.2e\n",tit->first,tit->second);
	}
      }
      if(r!=0) {
	SETERRQ(PETSC_COMM_SELF,1,"*** ERROR file cloase unsucessfull");
      }
      readFile=1;
      PetscFunctionReturn(0);
    }

    //Hassan: this fuction will loop over data in pair vector ts to find load scale based on ts_t
    PetscErrorCode scaleNf(const FEMethod *fe,ublas::vector<FieldData> &Nf) {
      PetscFunctionBegin;
      if(readFile==0) {
	SETERRQ(PETSC_COMM_SELF,1,"data file not readed");
      }
      double ts_t = fe->ts_t;
      double scale = 0;
      double t0 = 0,t1,s0 = 0,s1,dt;
      map<double, double>::iterator tit = tSeries.begin();
      for(;tit!=tSeries.end();tit++) {
	if(tit->first > ts_t) {
	  t1 = tit->first;
	  s1 = tit->second;
	  dt = ts_t - t0;
	  scale = s0 + ( (s1-s0)/(t1-t0) )*dt;
	  break;
	}
	t0 = tit->first;
	s0 = tit->second;
	scale = s0;
      }
      //Hassan : Here you can define time function rather than read from a file
      //Triangular loading over 10s (maximum at 5)
      /*double scale = 0;
      if(ts_t < 3.) scale = ts_t/5.;
      if(ts_t > 5.) scale = 1.+(5.-ts_t)/5.;
      if(ts_t > 10.) scale = 0;*/
      Nf *= scale;
      PetscFunctionReturn(0);
    }
};

