/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
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

#ifndef __GROUND_TIME_DATA_HPP
#define __GROUND_TIME_DATA_HPP

//#include <boost/program_options.hpp>
//#include <iostream>
//#include <fstream>
//#include <iterator>

using namespace std;
namespace po = boost::program_options;

/**
  *
  * This is based on:
  * NEW MODELS FOR SIMULATING DAILY MINIMUM, DAILY MAXIMUM AND HOURLY OUTDOOR TEMPERATURES
  * H Bulut1, O. Buyukalaca2 and T. Yılmaz2
  *
  * However some modifications are introduced, that model can be set using data
  * easly found on internet.
  *
  */
struct CrudeClimateModel: public GenricClimateModel {

  CrudeClimateModel(const char file[]): GenricClimateModel() {
    PetscErrorCode ierr;
    ierr = setDefault(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = readFile(file); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  }

  double TdaySummer,TdayWinter;
  double TnigthSummer,TnigthWinter;
  int Low,High;
  
  double fractinalDayTime(double t) {
    const double nb_days = 365;
    double n = t/(60*60*24);
    return n - floor(n/nb_days)*nb_days;
  }
  double wetherTime(double n) {
    const double nb_days = 365;
    double a = -0.1e1 / (double) nb_days / (double) (High - Low - nb_days) * M_PI * (double) (2 * High - 2 * Low - nb_days) / (double) (High - Low);
    double b = M_PI * (double) (2 * High * High - 2 * Low * Low - 2 * nb_days * Low - nb_days * nb_days) / (double) nb_days / (double) (High - Low - nb_days) / (double) (High - Low);
    double c = -0.1e1 / (double) nb_days / (double) (High - Low - nb_days) * M_PI * (double) (2 * High - 2 * Low - nb_days) / (double) (High - Low);
    return a+b*n+c*n*n;
  }
  double funTday(double t) {
    double n = fractinalDayTime(t);
    return 0.5*(TdayWinter+TdaySummer+(TdayWinter-TdaySummer)*cos( wetherTime(n) ));
  }
  double funTnigth(double t) {
    double n = fractinalDayTime(t);
    return 0.5*(TnigthWinter+TnigthSummer+(TnigthWinter-TnigthSummer)*cos( wetherTime(n) ));
  }

  double Td; // dwe point
  struct tm startTime;
  time_t t0;

  
  PetscErrorCode readFile(const char file[]) {
    PetscFunctionBegin;
    try {
      ifstream ini_file(file);  
      po::variables_map vm;
      po::options_description config_file_options;
      config_file_options.add_options()
	("longitude",po::value<double>(&spaData.longitude)->default_value(0.127))
	("latitude",po::value<double>(&spaData.latitude)->default_value(51.5072))
	("year",po::value<int>(&spaData.year)->default_value(2014))
	("month",po::value<int>(&spaData.month)->default_value(6))
	("day",po::value<int>(&spaData.day)->default_value(21))
	("hour",po::value<int>(&spaData.hour)->default_value(0))
	("minute",po::value<int>(&spaData.minute)->default_value(0))
	("second",po::value<double>(&spaData.second)->default_value(0))
	("timezone",po::value<double>(&spaData.timezone)->default_value(0))
	("elevation",po::value<double>(&spaData.elevation)->default_value(10))
	("TdayAtSummer",po::value<double>(&TdaySummer)->default_value(10),"Day temerature at summer")
	("TnigthAtSummer",po::value<double>(&TnigthSummer)->default_value(8),"Nigth temerature at summer")
  	("TdayAtWinter",po::value<double>(&TdayWinter)->default_value(-8),"Day temerature at winter")
	("TnigthAtWinter",po::value<double>(&TnigthWinter)->default_value(-10),"Nigth temerature at winter")
	("DayOfLowTemperature",po::value<int>(&Low)->default_value(0))
	("DayOfHighTemperature",po::value<int>(&High)->default_value(182))
	("DewPoint",po::value<double>(&Td)->default_value(0),"Dew point in Celsius degrees")
	("u10",po::value<double>(&u10)->default_value(2.7))
	("CR",po::value<double>(&CR)->default_value(0))
	("Pressure",po::value<double>(&P)->default_value(101325))
	("Rs",po::value<double>(&Rs)->default_value(0));
      store(parse_config_file(ini_file,config_file_options), vm);
      po::notify(vm); 

      startTime.tm_sec = spaData.second;
      startTime.tm_min = spaData.minute;
      startTime.tm_hour = spaData.hour;
      startTime.tm_mday = spaData.day;
      startTime.tm_mon = spaData.month-1;
      startTime.tm_year = spaData.year-1900;
      t0 = mktime(&startTime);

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  double funT(double t,double t0 = 0) {
    double Tmax = funTday(t);
    double Tmin = funTnigth(t);

    //cerr << "Tmax " << Tmax << " Tmin " << Tmin << endl;

    //spaData.suntransit,spaData.sunrise,spaData.sunset
    double td = (spaData.sunset-spaData.sunrise);

    double tmin = 12-td/2;
    double tmax = 12+(tmin*(12-tmin))/13.5;

    double h = (t-floor(t/(3600*24))*3600*24)/3600.;

    double T;
    if((h<tmin)||(h>tmax)) {
      double h_ = h;
      if(h<tmin) h_+=24;
      T = (Tmax-(Tmax-Tmin)*pow(sin(0.5*M_PI*(h_-tmax)/(24+tmin-tmax)),1.2));
    } else {
      T = Tmin+(Tmax-Tmin)*pow(sin(0.5*M_PI*(h-tmin)/(tmax-tmin)),1.4);
    }

    return T;
  }

  PetscErrorCode testCode(double duration,double step) {
    PetscFunctionBegin;

    /*cerr << spaData.second << endl;
    cerr << spaData.minute << endl;
    cerr << spaData.hour << endl;
    cerr << spaData.day << endl;
    cerr << spaData.month-1 << endl;
    cerr << spaData.year << endl;*/
    //cerr << TdaySummer << endl;
    //cerr << TdayWinter << endl;
    //cerr << TnigthSummer << endl;
    //cerr << TnigthWinter << endl;

    #ifdef BOOST_IOSTREAMS_TEE_HPP_INCLUDED
      typedef tee_device<ostream, ofstream> TeeDevice;
      typedef stream<TeeDevice> TeeStream;
      ofstream ofs("time_variation_model.txt");
      TeeDevice my_tee(cout, ofs); 
      TeeStream my_split(my_tee);
    #endif //BOOST_IOSTREAMS_TEE_HPP_INCLUDED

    time_t t = t0;
    for(;t<t0+duration;t+=step) {

      struct tm current_time;
      current_time = *gmtime(&t);//localtime(&t);
      spaData.second = current_time.tm_sec;
      spaData.minute = current_time.tm_min;
      spaData.hour = current_time.tm_hour;
      spaData.day = current_time.tm_mday;
      spaData.month = current_time.tm_mon+1;
      spaData.year = current_time.tm_year+1970;
              
      spaData.function = SPA_ZA_RTS;
      int r;
      r = spa_calculate(&spaData);
      if(r) {
	SETERRQ1(PETSC_COMM_SELF,1,"wrong input data for solar position calulator error codde %d",r);
      }
      
      double T = funT(t,t0);

      #ifndef BOOST_IOSTREAMS_TEE_HPP_INCLUDED
	PetscPrintf(PETSC_COMM_WORLD,"%3.4f %3.4f\n",((double)t-(double)t0)/(60*60*24),T);
      #else //BOOST_IOSTREAMS_TEE_HPP_INCLUDED
	my_split.precision(3);
	my_split << ((double)t-(double)t0)/(60*60*24) << " " << T << endl;
      #endif //BOOST_IOSTREAMS_TEE_HPP_INCLUDED

    }

    PetscFunctionReturn(0);

  }

  PetscErrorCode setDefault() {
    PetscFunctionBegin;

    spaData.delta_ut1 = 0;    	// Fractional second difference between UTC and UT which is used
				// to adjust UTC for earth's irregular rotation rate and is derived
				// from observation only and is reported in this bulletin:
				// http://maia.usno.navy.mil/ser7/ser7.dat,
				// where delta_ut1 = DUT1
				// valid range: -1 to 1 second (exclusive), error code 17

    spaData.delta_t = 0;     	// Difference between earth rotation time and terrestrial time
				// It is derived from observation only and is reported in this
				// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
				// where delta_t = 32.184 + (TAI-UTC) - DUT1
				// valid range: -8000 to 8000 seconds, error code: 7

    spaData.pressure = 1013.25; // Annual average local pressure [millibars]
				// valid range:    0 to 5000 millibars,       error code: 12

    spaData.temperature = 20;  	// Annual average local temperature [degrees Celsius]
				// valid range: -273 to 6000 degrees Celsius, error code; 13


    spaData.slope = 0;        	// Surface slope (measured from the horizontal plane)
				// valid range: -360 to 360 degrees, error code: 14

    spaData.azm_rotation = 0; 	// Surface azimuth rotation (measured from south to projection of
				// surface normal on horizontal plane, negative east)
				// valid range: -360 to 360 degrees, error code: 15

    spaData.atmos_refract = 0.5667; // Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
				// valid range: -5   to   5 degrees, error code: 16

    PetscFunctionReturn(0);
  }

  spa_data spaData;

  /** \brief set time in days
    */
  PetscErrorCode set(double time = 0) {
    PetscFunctionBegin;

    time_t t = t0;
    t += time*60*60*24;

    struct tm current_time;
    current_time = *gmtime(&t);
    spaData.second = current_time.tm_sec;
    spaData.minute = current_time.tm_min;
    spaData.hour = current_time.tm_hour;
    spaData.day = current_time.tm_mday;
    spaData.month = current_time.tm_mon+1;
    spaData.year = current_time.tm_year+1970;

    spaData.function = SPA_ZA_RTS;
    int r;
    r = spa_calculate(&spaData);
    if(r) {
      SETERRQ1(PETSC_COMM_SELF,1,"wrong input data for solar position calulator error codde %d",r);
    }
    zenith = spaData.zenith;
    azimuth = spaData.azimuth;
    
    Ta = funT(t,t0); // sets air temperature

    PetscPrintf(PETSC_COMM_WORLD,
      "Suntransit %3.2f Sunrise %3.2f Sunset %3.2f Ta = %3.2f \n" ,
      spaData.suntransit,spaData.sunrise,spaData.sunset,
      Ta);
    PetscPrintf(PETSC_COMM_WORLD,"Time %s",asctime(&current_time));


    PetscFunctionReturn(0);
  }

};


#endif //__GROUND_TIME_DATA_HPP

