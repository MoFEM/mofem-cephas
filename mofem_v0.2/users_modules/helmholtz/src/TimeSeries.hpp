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

  map<double,double> sSeries;

  PetscErrorCode readSpaceData() {
    PetscFunctionBegin;
    ierr = readData("-space_data",sSeries); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode readData() {
    PetscFunctionBegin;
    ierr = readSpaceData(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  boost::shared_array<kiss_fft_cpx> complexIn;  ///< Impulse amplitudes in space (imaginary)
  boost::shared_array<kiss_fft_cpx> complexOut; ///< Impulse amplitude in space (real)

  kiss_fft_cfg forwardCfg;
  kiss_fft_cfg inverseCfg;

  PetscErrorCode forwardDft(
    map<double,double> series,
    boost::shared_array<kiss_fft_cpx> &complex_in,
    boost::shared_array<kiss_fft_cpx> &complex_out
  ) {
    PetscFunctionBegin;

    int n = sSeries.size();

    complex_in = boost::shared_array<kiss_fft_cpx>(new kiss_fft_cpx[n]);
    map<double,double>::iterator mit = sSeries.begin();
    for(int ii = 0;mit!=sSeries.end();mit++,ii++) {
      complex_in[ii].r = mit->second;
      complex_in[ii].i = 0;
    }

    complex_out = boost::shared_array<kiss_fft_cpx>(new kiss_fft_cpx[n]);
    kiss_fft(forwardCfg,complex_in.get(),complex_out.get());

    PetscFunctionReturn(0);
  }


  vector<ublas::matrix<Vec> > pTimeFrequencySpaceFrequency;

  /**

  */
  PetscErrorCode forwardDftIncidentWave() {
    PetscFunctionBegin;

    PetscBool is_partitioned = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-my_is_partitioned",&is_partitioned,PETSC_NULL); CHKERRQ(ierr);
    double signal_length = 1; // wave speed
    ierr = PetscOptionsGetScalar(NULL,"-signal_length",&signal_length,NULL); CHKERRQ(ierr);
    double signal_duration = 1; // wave speed
    ierr = PetscOptionsGetScalar(NULL,"-signal_duration",&signal_duration,NULL); CHKERRQ(ierr);

    KSP approx_incindent_wave_solver;

    Mat A_approx_incident_wave;
    Vec D_approx_incident_wave;
    vector<Vec> F(2);
    vector<VecScatter> scatter_incident_wave(2);

    ierr = mField.MatCreateMPIAIJWithArrays("INCIDENT_WAVE",&A_approx_incident_wave); CHKERRQ(ierr);
    ierr = KSPCreate(PETSC_COMM_WORLD,&approx_incindent_wave_solver); CHKERRQ(ierr);
    ierr = KSPSetOperators(approx_incindent_wave_solver,A_approx_incident_wave,A_approx_incident_wave); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(approx_incindent_wave_solver); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("INCIDENT_WAVE",ROW,&F[0]); CHKERRQ(ierr);
    ierr = VecDuplicate(F[0],&F[1]); CHKERRQ(ierr);
    ierr = mField.VecCreateGhost("INCIDENT_WAVE",COL,&D_approx_incident_wave); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      D_approx_incident_wave,"INCIDENT_WAVE","rePRES",ROW,pSeries[0][0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave[0]
    ); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      D_approx_incident_wave,"INCIDENT_WAVE","imPRES",ROW,pSeries[1][0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave[1]
    ); CHKERRQ(ierr);

    pTimeFrequencySpaceFrequency.resize(2);

    int n = sSeries.size();
    for(int ss = 0;ss<2;ss++) {
      pTimeFrequencySpaceFrequency[ss].resize(n,n);
      for(int f = 0;f<n;f++) {
        for(int t = 0;t<n;t++) {
          if(f == 0 && t == 0) {
            ierr = mField.VecCreateGhost("PRESSURE_IN_TIME",ROW,&pTimeFrequencySpaceFrequency[ss](f,t)); CHKERRQ(ierr);
          } else {
            ierr = VecDuplicate(pTimeFrequencySpaceFrequency[ss](0,0),&pTimeFrequencySpaceFrequency[ss](f,t)); CHKERRQ(ierr);
          }
        }
      }
    }

    IncidentWave function_evaluator(
      0,
      helmholtzElement.globalParameters.waveDirection.first,
      0,0,0
    );

    ierr = MatZeroEntries(A_approx_incident_wave); CHKERRQ(ierr);
    ierr = calculate_matrix_and_vector(
      mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",A_approx_incident_wave,F,function_evaluator
    ); CHKERRQ(ierr);
    ierr = KSPSetUp(approx_incindent_wave_solver); CHKERRQ(ierr);

    // For this it is real time and space is in frequencies. This loop is done to calculate values of
    // DOFs at dictate time steps. Note that wave is moving with given speed, knowing duration of signal,
    // phase shift in space as function of time can be calculated. It is assumed that wave is
    // periodic.
    for(int f = 0;f<n;f++) {

      double wave_number = 2*M_PI*f/signal_length;
      function_evaluator.wAvenumber = wave_number;
      function_evaluator.amplitudeReal = complexOut[f].r;
      function_evaluator.amplitudeImag = complexOut[f].i;

      /*int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      if(!rank) {
        const complex< double > i( 0.0, 1.0 );
        complex< double > v = 0.0;

        for(int k = 0;k<n;k++) {
          v = v + (complexOut[k].r+i*complexOut[k].i)*exp(i*2.0*M_PI*(double)k*(double)f/(double)n);
        }
        cerr << v*1./(double)n << endl;

      }*/

      for(int t = 0;t<n;t++) {

        double speed = signal_length/signal_duration;
        double delta_t = signal_duration/n;
        double distance = speed*delta_t*t;
        function_evaluator.pHase = 2*M_PI*f*(distance/signal_length);

        for(int ss = 0;ss<2;ss++) {
          ierr = VecZeroEntries(F[ss]); CHKERRQ(ierr);
        }

        // Calculate matrix and factor it only once
        ierr = calculate_matrix_and_vector(
          mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",PETSC_NULL,F,function_evaluator
        ); CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD,"\n\nWhave number %d on time step %d, out of %d\n",f,t,n);

        // solve problem for complex and real values
        for(int ss = 0;ss<2;ss++) {
          // Solve incident wave approximation problem
          ierr = VecZeroEntries(D_approx_incident_wave); CHKERRQ(ierr);
          ierr = KSPSolve(approx_incindent_wave_solver,F[ss],D_approx_incident_wave); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecScatterBegin(
            scatter_incident_wave[ss],D_approx_incident_wave,pTimeFrequencySpaceFrequency[ss](f,t),INSERT_VALUES,SCATTER_FORWARD
          ); CHKERRQ(ierr);
          ierr = VecScatterEnd(
            scatter_incident_wave[ss],D_approx_incident_wave,pTimeFrequencySpaceFrequency[ss](f,t),INSERT_VALUES,SCATTER_FORWARD
          ); CHKERRQ(ierr);
        }

      }

    }

    int size;
    ierr = VecGetLocalSize(pTimeFrequencySpaceFrequency[0](0,0),&size); CHKERRQ(ierr);

    for(int ss = 0;ss<2;ss++) {
      for(int t = 0;t<n;t++) {
        ierr = VecZeroEntries(pSeries[ss][t]); CHKERRQ(ierr);
      }
    }

    for(int ss = 0;ss<2;ss++) {
      for(int t = 0;t<n;t++) {
        for(int f = 0;f<n;f++) {
          ierr = VecAXPY(pSeries[ss][t],1./(double)n,pTimeFrequencySpaceFrequency[ss](f,t)); CHKERRQ(ierr);
        }
        ierr = VecGhostUpdateBegin(pSeries[ss][t],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(pSeries[ss][t],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
    }


    for(int i = 0;i<size;i++) {
      double *a_real,*a_imag;
      for(int t = 0;t<n;t++) {
        ierr = VecGetArray(pSeries[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecGetArray(pSeries[1](t),&a_imag); CHKERRQ(ierr);
        complexIn[t].r = a_real[i];
        complexIn[t].i = a_imag[i];
        ierr = VecRestoreArray(pSeries[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(pSeries[1](t),&a_imag); CHKERRQ(ierr);
      }
      // from now space is in frequencies and time is in frequencies
      kiss_fft(forwardCfg,complexIn.get(),complexOut.get());
      for(int t= 0;t<n;t++) {
        ierr = VecGetArray(pSeries[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecGetArray(pSeries[1](t),&a_imag); CHKERRQ(ierr);
        a_real[i] = complexOut[t].r;
        a_imag[i] = complexOut[t].i;
        ierr = VecRestoreArray(pSeries[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(pSeries[1](t),&a_imag); CHKERRQ(ierr);
      }
    }

    ierr = KSPDestroy(&approx_incindent_wave_solver); CHKERRQ(ierr);
    ierr = MatDestroy(&A_approx_incident_wave); CHKERRQ(ierr);
    ierr = VecDestroy(&F[0]); CHKERRQ(ierr);
    ierr = VecDestroy(&F[1]); CHKERRQ(ierr);
    ierr = VecDestroy(&D_approx_incident_wave); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_incident_wave[0]); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_incident_wave[1]); CHKERRQ(ierr);

    for(int ss = 0;ss<2;ss++) {
      for(int f = 0;f<n;f++) {
        for(int t = 0;t<n;t++) {
          ierr = VecDestroy(&pTimeFrequencySpaceFrequency[ss](f,t)); CHKERRQ(ierr);
        }
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode forwardDft() {
    PetscFunctionBegin;
    int n = sSeries.size();
    forwardCfg = kiss_fft_alloc(n, 0, NULL, NULL);
    ierr = forwardDft(sSeries,complexIn,complexOut); CHKERRQ(ierr);
    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = forwardDftIncidentWave(); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0)  ;
  }

  vector<ublas::vector<Vec> > pSeries;
  vector<VecScatter> scatterPressure;

  PetscErrorCode createTimeVectorSeries(Vec T) {
    PetscFunctionBegin;

    int n = sSeries.size();

    pSeries.resize(2);
    pSeries[0].resize(n);
    pSeries[1].resize(n);

    for(int ss = 0;ss<2;ss++) {
      ierr = mField.VecCreateGhost("PRESSURE_IN_TIME",ROW,&pSeries[ss][0]); CHKERRQ(ierr);
      for(int k = 1;k<n;k++) {
        ierr = VecDuplicate(pSeries[ss][0],&pSeries[ss][k]); CHKERRQ(ierr);
      }
    }

    scatterPressure.resize(2);
    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","rePRES",ROW,pSeries[0][0],"PRESSURE_IN_TIME","P",ROW,&scatterPressure[0]
    ); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","imPRES",ROW,pSeries[0][0],"PRESSURE_IN_TIME","P",ROW,&scatterPressure[1]
    ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode destroyTimeVectorSeries() {
    PetscFunctionBegin;

    int n = sSeries.size();

    for(int ss = 0;ss<2;ss++) {
      for(int k = 0;n<k;k++) {
        ierr = VecDestroy(&pSeries[ss][k]); CHKERRQ(ierr);
      }
      ierr = VecScatterDestroy(&scatterPressure[ss]); CHKERRQ(ierr);
    }

    kiss_fft_cleanup();

    PetscFunctionReturn(0);
  }

  /** \brief Solve problem for each wave number

  */
  PetscErrorCode solveForwardDFT(KSP solver,Mat A,Vec F,Vec T) {
    PetscFunctionBegin;

    /*const bool only_incident_wave = true;
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
    VecScatter scatter_incident_wave_real;

    if(add_incident_wave) {

      ierr = mField.MatCreateMPIAIJWithArrays("INCIDENT_WAVE",&A_approx_incident_wave); CHKERRQ(ierr);

      ierr = KSPCreate(PETSC_COMM_WORLD,&approx_incindent_wave_solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(approx_incindent_wave_solver,A_approx_incident_wave,A_approx_incident_wave); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(approx_incindent_wave_solver); CHKERRQ(ierr);

      ierr = mField.VecCreateGhost("INCIDENT_WAVE",ROW,&vec_F[0]); CHKERRQ(ierr);
      ierr = mField.VecCreateGhost("INCIDENT_WAVE",ROW,&vec_F[1]); CHKERRQ(ierr);
      ierr = mField.VecCreateGhost("INCIDENT_WAVE",COL,&D_approx_incident_wave); CHKERRQ(ierr);

      ierr = mField.VecScatterCreate(
        D_approx_incident_wave,"INCIDENT_WAVE","rePRES",ROW,pSeriesReal[0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave_real
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

        // loop over frequencies in space and sum them up, for given wave number in time frequency
        for (int s = 0; s < ns; s++) {

          double real_power = complexOut[s].r*complexTimeOut[k].r - complexOut[s].i*complexTimeOut[k].i;
          double imag_power = complexOut[s].r*complexTimeOut[k].i + complexOut[s].i*complexTimeOut[k].r;

          IncidentWave function_evaluator(
            s/incident_wave_speed,
            helmholtzElement.globalParameters.waveDirection.first,
            real_power,imag_power
          );

          if(s==0) {
            for(int ss = 0;ss<2;ss++) {
              ierr = VecZeroEntries(vec_F[ss]); CHKERRQ(ierr);
            }
          }

          if(s == 0 && k == 0) {

            ierr = MatZeroEntries(A_approx_incident_wave); CHKERRQ(ierr);
            ierr = calculate_matrix_and_vector(
              mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",A_approx_incident_wave,vec_F,function_evaluator
            ); CHKERRQ(ierr);

            //{
              //Matrix View
            //  MatView(A_approx_incident_wave,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
            //  std::string wait;
            //  std::cin >> wait;
            //}
            ierr = KSPSetUp(approx_incindent_wave_solver); CHKERRQ(ierr);

          } else {

            ierr = calculate_matrix_and_vector(
              mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",PETSC_NULL,vec_F,function_evaluator
            ); CHKERRQ(ierr);
          }

        }

        for (size_t ss = 0; ss < 2; ss++) {
          // Solve incident wave approximation problem
          ierr = VecZeroEntries(D_approx_incident_wave); CHKERRQ(ierr);
          ierr = KSPSolve(approx_incindent_wave_solver,vec_F[ss],D_approx_incident_wave); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          if(ss) {
            ierr = VecScatterBegin(
              scatter_incident_wave_real,D_approx_incident_wave,pSeriesReal[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
            ierr = VecScatterEnd(
              scatter_incident_wave_real,D_approx_incident_wave,pSeriesReal[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
          } else {
            ierr = VecScatterBegin(
              scatter_incident_wave_real,D_approx_incident_wave,pSeriesImag[k],ADD_VALUES,SCATTER_FORWARD
            ); CHKERRQ(ierr);
            ierr = VecScatterEnd(
              scatter_incident_wave_real,D_approx_incident_wave,pSeriesImag[k],ADD_VALUES,SCATTER_FORWARD
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

      ierr = VecScatterDestroy(&scatter_incident_wave_real); CHKERRQ(ierr);
    }*/

    PetscFunctionReturn(0);
  }

  PetscErrorCode pressureInTimeDomainInverseDft() {
    PetscFunctionBegin;

    PostPocOnRefinedMesh post_proc(mField);
    ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("P"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);

    int n = sSeries.size();
    inverseCfg = kiss_fft_alloc(n, 1, NULL, NULL);

    int size;
    ierr = VecGetLocalSize(pSeries[0][0],&size); CHKERRQ(ierr);

    for(int ii = 0;ii<size;ii++) {

      for(int k = 0;k<n;k++) {

        double *p_real,*p_imag;
        ierr = VecGetArray(pSeries[0][k],&p_real); CHKERRQ(ierr);
        ierr = VecGetArray(pSeries[1][k],&p_imag); CHKERRQ(ierr);

        complexIn[k].r = p_real[ii];
        complexIn[k].i = p_imag[ii];

        ierr = VecRestoreArray(pSeries[0][k],&p_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(pSeries[1][k],&p_imag); CHKERRQ(ierr);

      }
      kiss_fft(inverseCfg,complexIn.get(),complexOut.get());
      for(int k = 0;k<n;k++) {

        double *a_p;
        ierr = VecGetArray(pSeries[0][k],&a_p); CHKERRQ(ierr);
        a_p[ ii ] = complexOut[k].r;
        ierr = VecRestoreArray(pSeries[0][k],&a_p); CHKERRQ(ierr);

        ierr = VecGetArray(pSeries[1][k],&a_p); CHKERRQ(ierr);
        a_p[ ii ] = complexOut[k].i;
        ierr = VecRestoreArray(pSeries[1][k],&a_p); CHKERRQ(ierr);

      }

    }

    for(int k = 0;k<n;k++) {

      ierr = VecGhostUpdateBegin(pSeries[0][k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeries[0][k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecGhostUpdateBegin(pSeries[1][k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(pSeries[1][k],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = mField.set_local_ghost_vector(
        "PRESSURE_IN_TIME",ROW,pSeries[0][k],INSERT_VALUES,SCATTER_REVERSE
      ); CHKERRQ(ierr);
      ierr = mField.loop_finite_elements("PRESSURE_IN_TIME","PRESSURE_FE",post_proc); CHKERRQ(ierr);

      {
        ostringstream ss;
        ss << "pressure_real_time_step_" << k << ".h5m";
        rval = post_proc.postProcMesh.write_file(ss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
        PetscPrintf(PETSC_COMM_WORLD,"Saved %s\n",ss.str().c_str());
      }

      ierr = mField.set_local_ghost_vector(
        "PRESSURE_IN_TIME",ROW,pSeries[1][k],INSERT_VALUES,SCATTER_REVERSE
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
