/* \file TimeSeries.hpp
  \ingroup mofem_helmholtz_elem

*/

/*
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/** \brief Read impulse and apply FFT
  */
struct TimeSeries {

  FieldInterface& mField;
  HelmholtzElement& helmholtzElement;
  AnalyticalDirihletBC::DirichletBC& analyticalDitihletBcReal;
  AnalyticalDirihletBC::DirichletBC& analyticalDitihletBcImag;

  bool dirichletBcSet;
  int readFile,debug;

  TimeSeries(FieldInterface &m_field,
    HelmholtzElement& helmholtz_element,
    AnalyticalDirihletBC::DirichletBC &analytical_ditihlet_bc_real,
    AnalyticalDirihletBC::DirichletBC &analytical_ditihlet_bc_imag,
    bool dirihlet_bc_set):
      mField(m_field),helmholtzElement(helmholtz_element),
      analyticalDitihletBcReal(analytical_ditihlet_bc_real),
      analyticalDitihletBcImag(analytical_ditihlet_bc_imag),
      dirichletBcSet(dirihlet_bc_set),
      readFile(0),debug(1) {}

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscErrorCode readData(const char* str,map<double,double> &series) {
    PetscFunctionBegin;

    char time_file_name[255];
    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetString(PETSC_NULL,str,time_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ1(PETSC_COMM_SELF,1,"*** ERROR %s (DATA FILE NEEDED)",str);
    }
    FILE *time_data = fopen(time_file_name,"r");
    if(time_data == NULL) {
      SETERRQ1(PETSC_COMM_SELF,1,"*** ERROR data file < %s > open unsucessfull",time_file_name);
    }
    double no1 = 0.0, no2 = 0.0;
    series[no1] = no2;
    while(! feof (time_data)){
      int n = fscanf(time_data,"%lf %lf",&no1,&no2);
      if((n <= 0)||((no1==0)&&(no2==0))) {
        fgetc(time_data);
        continue;
      }
      if(n != 2){
        SETERRQ1(PETSC_COMM_SELF,1,"*** ERROR read data file error (check input time data file) { n = %d }",n);
      }
      series[no1] = no2;
    }
    int r = fclose(time_data);

    if(debug) {
      map<double, double>::iterator tit = series.begin();
      for(;tit!=series.end();tit++) {
        PetscPrintf(PETSC_COMM_WORLD,"** read time series %3.2e time %3.2e\n",tit->first,tit->second);
      }
    }
    if(r!=0) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR file close unsuccessful");
    }
    readFile=1;

    PetscFunctionReturn(0);
  }

  map<double,double> tSeries;
  map<double,double> sSeries;

  PetscErrorCode readTimeData() {
    PetscFunctionBegin;
    ierr = readData("-time_data",tSeries); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode readSpaceData() {
    PetscFunctionBegin;
    ierr = readData("-space_data",sSeries); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode readData() {
    PetscFunctionBegin;

    ierr = readTimeData(); CHKERRQ(ierr);
    ierr = readSpaceData(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  boost::shared_array<kiss_fft_cpx> complexSpaceIn;  ///< Impulse amplitudes in space (imaginary)
  boost::shared_array<kiss_fft_cpx> complexSpaceOut; ///< Impulse amplitude in space (real)
  boost::shared_array<kiss_fft_cpx> complexTimeIn;  ///< Impulse amplitudes in time (imaginary)
  boost::shared_array<kiss_fft_cpx> complexTimeOut; ///< Impulse amplitude in time (real)

  kiss_fft_cfg forwardCfg;
  kiss_fft_cfg inverseCfg;

  PetscErrorCode forwardDft(
    map<double,double> series,
    boost::shared_array<kiss_fft_cpx> &complex_in,
    boost::shared_array<kiss_fft_cpx> &complex_out
  ) {
    PetscFunctionBegin;

    int n = tSeries.size();
    if(!(n%2)) {
      SETERRQ(PETSC_COMM_SELF,1,"odd number of number of points, should be even");
    }

    complex_in = boost::shared_array<kiss_fft_cpx>(new kiss_fft_cpx[n]);
    map<double,double>::iterator mit = tSeries.begin();
    for(int ii = 0;mit!=tSeries.end();mit++,ii++) {
      complex_in[ii].r = mit->second;
      complex_in[ii].i = 0;
    }

    forwardCfg = kiss_fft_alloc(n, 0, NULL, NULL);
    complex_out = boost::shared_array<kiss_fft_cpx>(new kiss_fft_cpx[n]);
    kiss_fft(forwardCfg,complex_in.get(),complex_out.get());

    PetscFunctionReturn(0);
  }

  PetscErrorCode forwardDft() {
    PetscFunctionBegin;


    ierr = forwardDft(tSeries,complexSpaceIn,complexSpaceOut); CHKERRQ(ierr);
    ierr = forwardDft(sSeries,complexTimeIn,complexTimeOut); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  ublas::vector<Vec> pSeriesReal,pSeriesImag;
  VecScatter scatterImag,scatterReal;

  PetscErrorCode createTimeVectorSeries(Vec T) {
    PetscFunctionBegin;

    int n = tSeries.size();

    pSeriesReal.resize(n);
    pSeriesImag.resize(n);

    ierr = mField.VecCreateGhost("PRESSURE_IN_TIME",ROW,&pSeriesReal[0]); CHKERRQ(ierr);
    for(int k = 1;k<n;k++) {
      ierr = VecDuplicate(pSeriesReal[0],&pSeriesReal[k]); CHKERRQ(ierr);
    }

    for(int k = 0;k<n;k++) {

      ierr = VecDuplicate(pSeriesReal[0],&pSeriesImag[k]); CHKERRQ(ierr);

    }

    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","rePRES",ROW,pSeriesReal[0],"PRESSURE_IN_TIME","P",ROW,&scatterReal
    ); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","imPRES",ROW,pSeriesImag[0],"PRESSURE_IN_TIME","P",ROW,&scatterImag
    ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode destroyTimeVectorSeries() {
    PetscFunctionBegin;

    int n = tSeries.size();

    for(int k = 0;n<k;k++) {

      ierr = VecDestroy(&pSeriesReal[k]); CHKERRQ(ierr);
      ierr = VecDestroy(&pSeriesImag[k]); CHKERRQ(ierr);

    }

    ierr = VecScatterDestroy(&scatterReal); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatterImag); CHKERRQ(ierr);

    kiss_fft_cleanup();

    PetscFunctionReturn(0);
  }

  /** \brief Solve problem for each wave number

  */
  PetscErrorCode solveForwardDFT(KSP solver,Mat A,Vec F,Vec T) {
    PetscFunctionBegin;

    const bool only_incident_wave = true;
    int nt = tSeries.size();
    int ns = sSeries.size();

    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    PetscBool is_partitioned = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-my_is_partitioned",&is_partitioned,PETSC_NULL); CHKERRQ(ierr);

    KSP approx_incindent_wave_solver;

    Mat A_approx_incident_wave;
    Vec D_approx_incident_wave;
    vector<Vec> vec_F;
    vec_F.resize(2);
    VecScatter scatter_incident_wave;

    if(add_incident_wave) {

      ierr = mField.MatCreateMPIAIJWithArrays("INCIDENT_WAVE",&A_approx_incident_wave); CHKERRQ(ierr);

      ierr = KSPCreate(PETSC_COMM_WORLD,&approx_incindent_wave_solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(approx_incindent_wave_solver,A_approx_incident_wave,A_approx_incident_wave); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(approx_incindent_wave_solver); CHKERRQ(ierr);
      ierr = KSPSetUp(approx_incindent_wave_solver); CHKERRQ(ierr);

      ierr = mField.VecCreateGhost("INCIDENT_WAVE",ROW,&vec_F[0]); CHKERRQ(ierr); /* real */
      ierr = mField.VecCreateGhost("INCIDENT_WAVE",ROW,&vec_F[1]); CHKERRQ(ierr); /* imag */
      ierr = mField.VecCreateGhost("INCIDENT_WAVE",COL,&D_approx_incident_wave); CHKERRQ(ierr);

      ierr = mField.VecScatterCreate(
        D_approx_incident_wave,"INCIDENT_WAVE","rePRES",ROW,pSeriesReal[0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave
      ); CHKERRQ(ierr);

    }

    for(int k = 0;k<nt;k++) {

      // Zero vectors
      ierr = VecZeroEntries(T); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);

      // Set wave number
      double reference_wave_speed = 1; // wave speed
      ierr = PetscOptionsGetScalar(NULL,"-reference_wave_speed",&reference_wave_speed,NULL); CHKERRQ(ierr);

      const double wave_number = 2*M_PI*k/reference_wave_speed;
      PetscPrintf(PETSC_COMM_WORLD,"\n\nCalculate frequency %d for wave number %3.4g out of %d\n",k,wave_number,nt);

      map<int,HelmholtzElement::VolumeData>::iterator vit = helmholtzElement.volumeData.begin();
      for(;vit != helmholtzElement.volumeData.end();vit++) {
        vit->second.waveNumber = wave_number;
      }
      map<int,HelmholtzElement::SurfaceData>::iterator sit =  helmholtzElement.sommerfeldBcData.begin();
      for(;sit != helmholtzElement.sommerfeldBcData.end(); sit++) {
        sit->second.aDmittance_imag = -wave_number;
      }
      sit =  helmholtzElement.baylissTurkelBcData.begin();
      for(;sit != helmholtzElement.baylissTurkelBcData.end(); sit++) {
        sit->second.aDmittance_imag = -wave_number;
      }

      helmholtzElement.globalParameters.waveNumber.first = wave_number;
      //helmholtzElement.globalParameters.powerOfIncidentWaveReal.first = complexOut[k].r;
      //helmholtzElement.globalParameters.powerOfIncidentWaveImag.first = complexOut[k].i;
      //PetscPrintf(PETSC_COMM_WORLD,"Complex amplitude %6.4e + i%6.4e\n",complexOut[k].r,complexOut[k].i);

      if(!only_incident_wave) {

        // Assemble problem
        if(dirichletBcSet) {
          ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analyticalDitihletBcReal); CHKERRQ(ierr);
          ierr = mField.problem_basic_method_preProcess("ACOUSTIC_PROBLEM",analyticalDitihletBcImag); CHKERRQ(ierr);
        }

        ierr = helmholtzElement.calculateA("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
        ierr = helmholtzElement.calculateF("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);

        if(dirichletBcSet) {
          ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analyticalDitihletBcReal); CHKERRQ(ierr);
          ierr = mField.problem_basic_method_postProcess("ACOUSTIC_PROBLEM",analyticalDitihletBcImag); CHKERRQ(ierr);
        }

        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = VecScale(F,-1); CHKERRQ(ierr);

        // Solve problem
        ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
        ierr = KSPSetUp(solver); CHKERRQ(ierr);

        ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);

      }

      //ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      //ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecScatterBegin(scatterReal,T,pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(scatterReal,T,pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterBegin(scatterImag,T,pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecScatterEnd(scatterImag,T,pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      if(add_incident_wave) {

        double incident_wave_speed = 1; // wave speed
        ierr = PetscOptionsGetScalar(NULL,"-incident_wave_speed",&incident_wave_speed,NULL); CHKERRQ(ierr);

        // loop over frequencies in space
        for (int s = 0; s < ns; s++) {

          double real_power = complexSpaceOut[s].r*complexTimeOut[k].r - complexSpaceOut[s].i*complexTimeOut[k].i;
          double imag_power = complexSpaceOut[s].r*complexTimeOut[k].i + complexSpaceOut[s].i*complexTimeOut[k].r;

          IncidentWave function_evaluator(
            s/incident_wave_speed,
            helmholtzElement.globalParameters.waveDirection.first,
            real_power,imag_power
          );

          if(s==0 && k == 0) {
            ierr = calculate_matrix_and_vector(
              mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRESS",A_approx_incident_wave,vec_F,function_evaluator
            ); CHKERRQ(ierr);
          } else {
            ierr = calculate_matrix_and_vector(
              mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRESS",PETSC_NULL,vec_F,function_evaluator
            ); CHKERRQ(ierr);
          }

        }

        for (size_t ss = 0; ss < 2; ss++) {
          // Solve incident wave approximation problem
          ierr = KSPSolve(approx_incindent_wave_solver,vec_F[ss],D_approx_incident_wave); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          if(ss) {
            ierr = VecScatterBegin(
              scatter_incident_wave,D_approx_incident_wave,pSeriesReal[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
            ierr = VecScatterEnd(
              scatter_incident_wave,D_approx_incident_wave,pSeriesReal[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
          } else {
            ierr = VecScatterBegin(
              scatter_incident_wave,D_approx_incident_wave,pSeriesImag[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
            ierr = VecScatterEnd(
              scatter_incident_wave,D_approx_incident_wave,pSeriesImag[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
          }
        }

      }

      ierr = VecGhostUpdateBegin(pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    }

    if(add_incident_wave) {

      ierr = KSPDestroy(&approx_incindent_wave_solver); CHKERRQ(ierr);

      ierr = MatDestroy(&A_approx_incident_wave); CHKERRQ(ierr);
      ierr = VecDestroy(&vec_F[0]); CHKERRQ(ierr);
      ierr = VecDestroy(&vec_F[1]); CHKERRQ(ierr);
      ierr = VecDestroy(&D_approx_incident_wave); CHKERRQ(ierr);

      ierr = VecScatterDestroy(&scatter_incident_wave); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode pressureInTimeDomainInverseDft() {
    PetscFunctionBegin;

    PostPocOnRefinedMesh post_proc(mField);
    ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("P"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);

    int n = tSeries.size();
    inverseCfg = kiss_fft_alloc(n, 1, NULL, NULL);

    int size;
    ierr = VecGetLocalSize(pSeriesReal[0],&size); CHKERRQ(ierr);

    for(int ii = 0;ii<size;ii++) {

      for(int k = 0;k<n;k++) {

        double *p_real,*p_imag;
        ierr = VecGetArray(pSeriesReal[k],&p_real); CHKERRQ(ierr);
        ierr = VecGetArray(pSeriesImag[k],&p_imag); CHKERRQ(ierr);

        //complexOut[k].r = p_real[ii];
        //complexOut[k].i = p_imag[ii];

        ierr = VecRestoreArray(pSeriesReal[k],&p_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(pSeriesImag[k],&p_imag); CHKERRQ(ierr);

      }

      //kiss_fft(inverseCfg,complexOut.get(),complexIn.get());

      for(int k = 0;k<n;k++) {

        double *a_p;
        ierr = VecGetArray(pSeriesReal[k],&a_p); CHKERRQ(ierr);
        //a_p[ ii ] = complexIn[k].r;
        ierr = VecRestoreArray(pSeriesReal[k],&a_p); CHKERRQ(ierr);

        ierr = VecGetArray(pSeriesImag[k],&a_p); CHKERRQ(ierr);
        //a_p[ ii ] = complexIn[k].i;
        ierr = VecRestoreArray(pSeriesImag[k],&a_p); CHKERRQ(ierr);

      }

    }

    for(int k = 0;k<n;k++) {

      ierr = VecGhostUpdateBegin(pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeriesReal[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeriesImag[k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = mField.set_local_ghost_vector(
        "PRESSURE_IN_TIME",ROW,pSeriesReal[k],INSERT_VALUES,SCATTER_REVERSE
      ); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("PRESSURE_IN_TIME","PRESSURE_FE",post_proc); CHKERRQ(ierr);

      {
        ostringstream ss;
        ss << "pressure_real_time_step_" << k << ".h5m";
        rval = post_proc.postProcMesh.write_file(ss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
        PetscPrintf(PETSC_COMM_WORLD,"Saved %s\n",ss.str().c_str());
      }

      ierr = mField.set_local_ghost_vector(
        "PRESSURE_IN_TIME",ROW,pSeriesImag[k],INSERT_VALUES,SCATTER_REVERSE
      ); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("PRESSURE_IN_TIME","PRESSURE_FE",post_proc); CHKERRQ(ierr);

      {
        ostringstream ss;
        ss << "pressure_imag_time_step_" << k << ".h5m";
        rval = post_proc.postProcMesh.write_file(ss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
        PetscPrintf(PETSC_COMM_WORLD,"Saved %s\n",ss.str().c_str());
      }

    }

    PetscFunctionReturn(0);
  }

};
