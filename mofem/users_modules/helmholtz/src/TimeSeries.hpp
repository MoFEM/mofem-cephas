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
  \ingroup mofem_helmholtz_elem

 First signal is read from text file. Currently is assumed that this signal
gives profile for plane wave. Spherical or other wave types can be easily
added and implemented.

 Before signal is applied the Discrete Fourier Transform (in space) is applied. Next
for each wavelength in space a phase shift is applied according to wave
speed. With that at hand incident wave and boundary conditions are applied
for each discrete point in time.

 Note that wave in space is superposition of incident waves. This
superposition makes wave with approximated wave profile which is read from
text file. Moreover obtained wave fulfill wave propagation equation a priori.

 Next Dictate Fourier Transform in time domain is applied to get frequency
(time) space. Having boundary conditions in frequency domain, Helmholtz
problem is solved for each wave length.

 At this point plane wave (or other wave) and scattering wave in frequency
space is known. The last and final step is to apply inverse Discrete Fourier
Transform to get acoustic wave pressures in discrete times.

Note that analysis is with assumption that all signals are periodic, as
result of essential property of Dictate Fourier Transform.

In following procedure a KISS FTP library <http://sourceforge.net/projects/kissfft/>
is used to do Fast Fourier Transform.

  */
struct TimeSeries {

  FieldInterface& mField;
  HelmholtzElement& helmholtzElement;
  AnalyticalDirichletBC::DirichletBC& analyticalDitihletBcReal;
  AnalyticalDirichletBC::DirichletBC& analyticalDitihletBcImag;

  bool dirichletBcSet;
  int readFile,debug;

  PostProcVolumeOnRefinedMesh postProc;

  TimeSeries(FieldInterface &m_field,
    HelmholtzElement& helmholtz_element,
    AnalyticalDirichletBC::DirichletBC &analytical_ditihlet_bc_real,
    AnalyticalDirichletBC::DirichletBC &analytical_ditihlet_bc_imag,
    bool Dirichlet_bc_set):
      mField(m_field),helmholtzElement(helmholtz_element),
      analyticalDitihletBcReal(analytical_ditihlet_bc_real),
      analyticalDitihletBcImag(analytical_ditihlet_bc_imag),
      dirichletBcSet(Dirichlet_bc_set),
      readFile(0),
      debug(1),
      postProc(m_field)
      {}

  ErrorCode rval;
  PetscErrorCode ierr;

  /** \brief Read signal from text file
   */
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
    while(! feof (time_data)){ //check if end-of-file is reached
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

  /** \brief This is wrapper to reading signal profile in space
   */
  PetscErrorCode readSpaceData() {
    PetscFunctionBegin;
    ierr = readData("-space_data",sSeries); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /** \brief Read all text files

    At this point only signal profile in space is read.

  */
  PetscErrorCode readData() {
    PetscFunctionBegin;
    ierr = readSpaceData(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  boost::shared_array<kiss_fft_cpx> complexIn;  ///< Impulse amplitudes in space (imaginary)
  boost::shared_array<kiss_fft_cpx> complexOut; ///< Impulse amplitude in space (real)

  kiss_fft_cfg forwardCfg;
  kiss_fft_cfg inverseCfg;

  /** \brief Apply Forward FFT for signal.
  */
  PetscErrorCode forwardSignalDft(
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
    /* forwardCfg returned from kiss_fft_alloc */
    kiss_fft(forwardCfg,complex_in.get(),complex_out.get());

    /*const complex< double > i( 0.0, 1.0 );
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    for(int f = 0;f<n;f++) {
      double wave_number = 2*M_PI*f/(double)n;
      complex<double> p = 0;
      for(int k = 0;k<n;k++) {
        p += (complex_out[k].r+i*complex_out[k].i)*exp(i*(wave_number*k) - i*(2*M_PI*k*0.21) );
      }
      p /= n;
      if(!rank) {
        cerr << std::real(p) << " " << std::imag(p) << endl;
      }
    }*/

    PetscFunctionReturn(0);
  }

  /** \brief Apply Inverse FFT for incident wave
   */
  PetscErrorCode forwardDftIncidentWave() {
    PetscFunctionBegin;

    PetscBool is_partitioned = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-my_is_partitioned",&is_partitioned,PETSC_NULL); CHKERRQ(ierr);

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
      D_approx_incident_wave,"INCIDENT_WAVE","rePRES",ROW,pSeriesIncidentWave[0][0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave[0]
    ); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      D_approx_incident_wave,"INCIDENT_WAVE","imPRES",ROW,pSeriesIncidentWave[1][0],"PRESSURE_IN_TIME","P",ROW,&scatter_incident_wave[1]
    ); CHKERRQ(ierr);

    if(!helmholtzElement.globalParameters.signalLength.second) {
      SETERRQ(PETSC_COMM_SELF,1,"Signal length not set, -signal_length");
    }
    if(!helmholtzElement.globalParameters.signalDuration.second) {
      SETERRQ(PETSC_COMM_SELF,1,"Signal duration not set, -signal_duration");
    }

    int n = sSeries.size();
    IncidentWaveDFT function_evaluator(
      helmholtzElement.globalParameters.signalLength.first,
      helmholtzElement.globalParameters.signalDuration.first,
      helmholtzElement.globalParameters.waveDirection.first,
      complexOut,n,0
    );

    ierr = MatZeroEntries(A_approx_incident_wave); CHKERRQ(ierr);
    ierr = calculate_matrix_and_vector(
      mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",A_approx_incident_wave,F,function_evaluator
    ); CHKERRQ(ierr);
    ierr = KSPSetUp(approx_incindent_wave_solver); CHKERRQ(ierr);


    // For this it is real time and space is in frequencies.
    // It is assumed that wave is periodic.
    for(int t = 0;t<n;t++) {

      PetscPrintf(PETSC_COMM_WORLD,"\n\nTime step %d, out of %d\n",t,n);

      for(int ss = 0;ss<2;ss++) {
        ierr = VecZeroEntries(F[ss]); CHKERRQ(ierr);
      }

      function_evaluator.timeStep = t;
      // Calculate matrix and factor it only once
      ierr = calculate_matrix_and_vector(
        mField,"INCIDENT_WAVE","HELMHOLTZ_RERE_FE","rePRES",PETSC_NULL,F,function_evaluator
      ); CHKERRQ(ierr);


      // Solve incident wave approximation problem
      for(int ss = 0;ss<2;ss++) {
        ierr = VecZeroEntries(D_approx_incident_wave); CHKERRQ(ierr);
        ierr = KSPSolve(approx_incindent_wave_solver,F[ss],D_approx_incident_wave); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D_approx_incident_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecScatterBegin(
          scatter_incident_wave[ss],D_approx_incident_wave,pSeriesIncidentWave[ss](t),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
        ierr = VecScatterEnd(
          scatter_incident_wave[ss],D_approx_incident_wave,pSeriesIncidentWave[ss](t),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
      }

    }

    ierr = KSPDestroy(&approx_incindent_wave_solver); CHKERRQ(ierr);
    ierr = MatDestroy(&A_approx_incident_wave); CHKERRQ(ierr);
    ierr = VecDestroy(&F[0]); CHKERRQ(ierr);
    ierr = VecDestroy(&F[1]); CHKERRQ(ierr);
    ierr = VecDestroy(&D_approx_incident_wave); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_incident_wave[0]); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_incident_wave[1]); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  /** \brief Apply Froward FFT and compose wave

   */
  PetscErrorCode forwardSpaceDft() {
    PetscFunctionBegin;
    int n = sSeries.size();
    forwardCfg = kiss_fft_alloc(n, 0, NULL, NULL);

    // Transform impulse in physical space to frequency space domain
    ierr = forwardSignalDft(sSeries,complexIn,complexOut); CHKERRQ(ierr);

    // Add incident wave
    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = forwardDftIncidentWave(); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0)  ;
  }

  vector<ublas::vector<Vec> > pSeriesIncidentWave;
  vector<ublas::vector<Vec> > pSeriesScatterWave;
  vector<VecScatter> scatterPressure;


  PetscErrorCode createSeries(Vec T,vector<ublas::vector<Vec> > &p_series) {
    PetscFunctionBegin;

    int n = sSeries.size();

    p_series.resize(2);
    p_series[0].resize(n);
    p_series[1].resize(n);

    for(int ss = 0;ss<2;ss++) {
      ierr = mField.VecCreateGhost("PRESSURE_IN_TIME",ROW,&p_series[ss][0]); CHKERRQ(ierr);
      for(int k = 1;k<n;k++) {
        //Duplicate the ghost format to all the vectors inside real and imaginary containers.
        ierr = VecDuplicate(p_series[ss][0],&p_series[ss][k]); CHKERRQ(ierr);
      }
    }


    PetscFunctionReturn(0);
  }

  /** \brief Destroy pressure series.
  */
  PetscErrorCode destroySeries(vector<ublas::vector<Vec> > &p_series) {
    PetscFunctionBegin;

    int n = sSeries.size();

    for(int ss = 0;ss<2;ss++) {
      for(int k = 0;n<k;k++) {
        ierr = VecDestroy(&p_series[ss][k]); CHKERRQ(ierr);
      }
    }

    PetscFunctionReturn(0);
  }

  /** \brief Create vectors for each wave lengths

 There are two series, for real and imaginary values. Each series keeps
pressures at time steps or wave length depending on stage of algorithm.

  */
  PetscErrorCode createPressureSeries(Vec T) {
    PetscFunctionBegin;

    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = createSeries(T,pSeriesIncidentWave); CHKERRQ(ierr);
    }

    ierr = createSeries(T,pSeriesScatterWave); CHKERRQ(ierr);

    scatterPressure.resize(2);
    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","rePRES",ROW,pSeriesScatterWave[0][0],"PRESSURE_IN_TIME","P",ROW,&scatterPressure[0]
    ); CHKERRQ(ierr);
    ierr = mField.VecScatterCreate(
      T,"ACOUSTIC_PROBLEM","imPRES",ROW,pSeriesScatterWave[0][0],"PRESSURE_IN_TIME","P",ROW,&scatterPressure[1]
    ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  /** \brief Destroy pressure series.
  */
  PetscErrorCode destroyPressureSeries() {
    PetscFunctionBegin;


    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = destroySeries(pSeriesIncidentWave); CHKERRQ(ierr);
    }

    ierr = destroySeries(pSeriesScatterWave); CHKERRQ(ierr);

    for(int ss = 0;ss<2;ss++) {
      ierr = VecScatterDestroy(&scatterPressure[ss]); CHKERRQ(ierr);
    }

    kiss_fft_cleanup();

    PetscFunctionReturn(0);
  }

  /** \brief Solve problem for each wave number.

 Calculate right hand vector for boundary condition at each time step. Each
boundary condition is calculated as superposition of plane waves, having given
Fourier Transform of signal and phase shift in space according to wave speed in medium.

Next Forward FFT is applied to right hand vector and for each wave length
Helmholtz problem is solved. This is most time consuming stage of algorithm.

Currently only boundary conditions are implemented for hard and mix surface. Other types
of boundary conditions could be easily implemented.

  */
  PetscErrorCode solveForwardDFT(KSP solver,Mat A,Vec F,Vec T) {
    PetscFunctionBegin;

    int n = sSeries.size();
    helmholtzElement.globalParameters.complexOut = complexOut;
    helmholtzElement.globalParameters.complexOutSize = n;

    // Zero vectors
    ierr = VecZeroEntries(T); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    vector<ublas::vector<Vec> > rhs_series;
    ierr = createSeries(T,rhs_series); CHKERRQ(ierr);

    double signal_length = helmholtzElement.globalParameters.signalLength.first;
    double signal_duration = helmholtzElement.globalParameters.signalDuration.first;
    double speed = signal_length/signal_duration;

    // Calculate boundary conditions at time steps
    for(int t = 0;t<n;t++) {

      double wave_number = 2*M_PI*t/signal_length;

      PetscPrintf(PETSC_COMM_WORLD,"Calculate rhs for frequency %d for wave number (space) %3.4g out of %d\n",t,wave_number,n);

      helmholtzElement.globalParameters.waveNumber.first = wave_number;
      helmholtzElement.globalParameters.timeStep = t;

      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = helmholtzElement.calculateF("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

      for(int ss = 0;ss<2;ss++) {
        ierr = VecScatterBegin(
          scatterPressure[ss],F,rhs_series[ss](t),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
        ierr = VecScatterEnd(
          scatterPressure[ss],F,rhs_series[ss](t),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
      }

    }

    // From now transform right hand vector to frequency domain
    ierr = seriesForwardDft(rhs_series); CHKERRQ(ierr);

    // Solve problem in frequency domain
    for(int f = 0;f<n;f++) {

      double wave_number = 2*M_PI*f/speed;
      PetscPrintf(PETSC_COMM_WORLD,"\n\nSolve Helmholtz problem for frequency %d for wave number (time) %3.4g out of %d\n",f,wave_number,n);

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

      // Zero matrix
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      // Calculate matrix
      ierr = helmholtzElement.calculateA("ACOUSTIC_PROBLEM"); CHKERRQ(ierr);
      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

      // Solve problem
      ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);


      for(int ss = 0;ss<2;ss++) {
        ierr = VecScatterBegin(
          scatterPressure[ss],rhs_series[ss](f),F,INSERT_VALUES,SCATTER_REVERSE
        ); CHKERRQ(ierr);
        ierr = VecScatterEnd(
          scatterPressure[ss],rhs_series[ss](f),F,INSERT_VALUES,SCATTER_REVERSE
        ); CHKERRQ(ierr);
      }

      // Zero vectors
      ierr = VecZeroEntries(T); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);

      for(int ss = 0;ss<2;ss++) {
        ierr = VecScatterBegin(
          scatterPressure[ss],T,pSeriesScatterWave[ss](f),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
        ierr = VecScatterEnd(
          scatterPressure[ss],T,pSeriesScatterWave[ss](f),INSERT_VALUES,SCATTER_FORWARD
        ); CHKERRQ(ierr);
      }

    }

    ierr = destroySeries(rhs_series); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  /** \brief Apply Forward FFT to data series.

  Note that series is vector of vectors consisting DOFs.

  */
  PetscErrorCode seriesForwardDft(vector<ublas::vector<Vec> > &series) {
    PetscFunctionBegin;

    int n = sSeries.size();

    // Transform physical time domain to frequency time domain
    int size;
    ierr = VecGetLocalSize(series[0](0),&size); CHKERRQ(ierr);

    /*for(int ss = 0;ss<2;ss++) {
      for(int t = 0;t<n;t++) {
        ierr = VecScale(series[ss](t),1./(double)n); CHKERRQ(ierr);
      }
    }*/

    for(int i = 0;i<size;i++) {
      double *a_real,*a_imag;
      for(int t = 0;t<n;t++) {
        ierr = VecGetArray(series[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecGetArray(series[1](t),&a_imag); CHKERRQ(ierr);
        complexIn[t].r = a_real[i];
        complexIn[t].i = a_imag[i];
        ierr = VecRestoreArray(series[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(series[1](t),&a_imag); CHKERRQ(ierr);
      }
      kiss_fft(forwardCfg,complexIn.get(),complexOut.get());
      for(int t= 0;t<n;t++) {
        ierr = VecGetArray(series[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecGetArray(series[1](t),&a_imag); CHKERRQ(ierr);
        a_real[i] = complexOut[t].r;
        a_imag[i] = complexOut[t].i;
        ierr = VecRestoreArray(series[0](t),&a_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(series[1](t),&a_imag); CHKERRQ(ierr);
      }
    }

    PetscFunctionReturn(0);
  }


  /** \brief Apply Inverse FFT to data series.

  */
  PetscErrorCode seriesInverseDft(vector<ublas::vector<Vec> > &series) {
    PetscFunctionBegin;

    int n = sSeries.size();
    inverseCfg = kiss_fft_alloc(n, 1, NULL, NULL);

    int size;
    ierr = VecGetLocalSize(series[0][0],&size); CHKERRQ(ierr);

    for(int ii = 0;ii<size;ii++) {

      for(int k = 0;k<n;k++) {

        double *p_real,*p_imag;
        ierr = VecGetArray(series[0][k],&p_real); CHKERRQ(ierr);
        ierr = VecGetArray(series[1][k],&p_imag); CHKERRQ(ierr);

        complexIn[k].r = p_real[ii];
        complexIn[k].i = p_imag[ii];

        ierr = VecRestoreArray(series[0][k],&p_real); CHKERRQ(ierr);
        ierr = VecRestoreArray(series[1][k],&p_imag); CHKERRQ(ierr);

      }
      kiss_fft(inverseCfg,complexIn.get(),complexOut.get());
      for(int k = 0;k<n;k++) {

        double *a_p;
        ierr = VecGetArray(series[0][k],&a_p); CHKERRQ(ierr);
        a_p[ ii ] = complexOut[k].r;
        ierr = VecRestoreArray(series[0][k],&a_p); CHKERRQ(ierr);

        ierr = VecGetArray(series[1][k],&a_p); CHKERRQ(ierr);
        a_p[ ii ] = complexOut[k].i;
        ierr = VecRestoreArray(series[1][k],&a_p); CHKERRQ(ierr);

      }

    }

    PetscFunctionReturn(0);
  }

  /// Aplly Forward FFT to pressure degrees of freedom
  PetscErrorCode pressureForwardDft() {
    PetscFunctionBegin;
    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = seriesForwardDft(pSeriesIncidentWave); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  /// Aplly Inverse FFT to pressure degrees of freedom
  PetscErrorCode pressureInverseDft() {
    PetscFunctionBegin;
    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = seriesInverseDft(pSeriesIncidentWave); CHKERRQ(ierr);
    }
    ierr = seriesInverseDft(pSeriesScatterWave); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode generateReferenceElementMesh() {
    PetscFunctionBegin;
    ierr = postProc.generateReferenceElementMesh(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  /// Save data on mesh
  PetscErrorCode saveResults() {
    PetscFunctionBegin;

    ierr = postProc.clearOperators(); CHKERRQ(ierr);
    ierr = postProc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);

    Vec p_inicent_wave;
    PetscBool add_incident_wave = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,"-add_incident_wave",&add_incident_wave,NULL); CHKERRQ(ierr);
    if(add_incident_wave) {
      ierr = VecDuplicate(pSeriesIncidentWave[0](0),&p_inicent_wave); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("P","P_INCIDENT_WAVE",p_inicent_wave); CHKERRQ(ierr);
    }

    Vec p_scatter_wave_real;
    ierr = VecDuplicate(pSeriesScatterWave[0](0),&p_scatter_wave_real); CHKERRQ(ierr);
    ierr = postProc.addFieldValuesPostProc("P","P_SCATTER_WAVE_REAL",p_scatter_wave_real); CHKERRQ(ierr);

    Vec p_scatter_wave_imag;
    ierr = VecDuplicate(pSeriesScatterWave[1](0),&p_scatter_wave_imag); CHKERRQ(ierr);
    ierr = postProc.addFieldValuesPostProc("P","P_SCATTER_WAVE_IMAG",p_scatter_wave_imag); CHKERRQ(ierr);

    int n = sSeries.size();
    for(int k = 0;k<n;k++) {

      if(add_incident_wave) {
        ierr = VecCopy(pSeriesIncidentWave[0](k),p_inicent_wave); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(p_inicent_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(p_inicent_wave,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }

      ierr = VecCopy(pSeriesScatterWave[0](k),p_scatter_wave_real); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(p_scatter_wave_real,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(p_scatter_wave_real,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = VecCopy(pSeriesScatterWave[1](k),p_scatter_wave_imag); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(p_scatter_wave_imag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(p_scatter_wave_imag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      ierr = mField.loop_finite_elements("PRESSURE_IN_TIME","PRESSURE_FE",postProc); CHKERRQ(ierr);

      {
        ostringstream ss;
        ss << "pressure_real_time_step_" << k << ".h5m";
        rval = postProc.postProcMesh.write_file(ss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
        PetscPrintf(PETSC_COMM_WORLD,"Saved %s\n",ss.str().c_str());
      }

    }

    if(add_incident_wave) {
      ierr = VecDestroy(&p_inicent_wave); CHKERRQ(ierr);
    }
    ierr = VecDestroy(&p_scatter_wave_real); CHKERRQ(ierr);
    ierr = VecDestroy(&p_scatter_wave_imag); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

};
