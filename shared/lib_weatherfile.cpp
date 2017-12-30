/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (�Alliance�) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as �System Advisor Model� or �SAM�. Except
*  to comply with the foregoing, the terms �System Advisor Model�, �SAM�, or any confusingly similar
*  designation may not be used to refer to any modified version of this software or any modified
*  version of the underlying software originally provided by Alliance without the prior written consent
*  of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/

#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <numeric>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

#if defined(__WINDOWS__)||defined(WIN32)||defined(_WIN32)
#define CASECMP(a,b) _stricmp(a,b)
#define CASENCMP(a,b,n) _strnicmp(a,b,n)
#else
#define CASECMP(a,b) strcasecmp(a,b) 
#define CASENCMP(a,b,n) strncasecmp(a,b,n)
#endif

#include "lib_util.h"
#include "lib_weatherfile.h"

using std::stof;
using std::stoi;

#ifdef _MSC_VER
#define my_isnan(x) ::_isnan( x )
#else
#define my_isnan(x) std::isnan( x )
#endif

static std::string trimboth(std::string &buf)
{
  	const auto strBegin = buf.find_first_not_of(" \t");
	if (strBegin == std::string::npos)
	  return std::string();

	const auto strEnd = buf.find_last_not_of(" \t\r\n");
	const auto strRange = strEnd - strBegin + 1;
	return buf.substr(strBegin, strRange);
}

static std::vector<std::string> split(const std::string &buf, char delim = ',')
{
	std::string token;
	std::vector<std::string> tokens;
	std::istringstream tokenStream(buf);
	while (std::getline(tokenStream, token, delim))
		tokens.push_back(token);
	return tokens;
}

static double conv_deg_min_sec(double degrees,
	double minutes,
	double seconds,
	char direction)
{
	double dd = degrees + minutes / 60.0 + seconds / 3600.0;
	if (tolower((int)direction) == 's' || tolower((int)direction) == 'w')
		dd = 0 - dd;
	return dd;
}

static int cmp_ext(const char *file, const char *ext)
{
	size_t len_ext, len_file;
	const char *extp;

	if (!file || !ext) return 0;
	len_ext = strlen(ext);
	len_file = strlen(file);
	extp = file + len_file - len_ext;

	if (extp < file) return 0;
	return CASENCMP(extp, ext, len_ext) == 0 ? 1 : 0;
}


std::string weatherfile::normalize_city(const std::string &in)
{
	std::string city = util::lower_case(in);
	util::replace(city, "?", " ");
	util::replace(city, "_", " ");
	util::replace(city, "\"", "");
	util::replace(city, "/", " ");
	util::replace(city, "\\", " ");

	for (size_t i = 0; i<city.length(); i++)
	{
		if (i == 0 || city[i - 1] == ' ')
			city[i] = (char)toupper((int)city[i]);
	}
	return city;
}


/***************************************************************************\

Function humidity()
This function calculates the relative humidity(%) based on the drybulb
temperature(C) and the dewpoint temperature.  It uses equations and
procedures presented in the 1993 ASHRAE Fundamentals Handbook, p6.7-10.
If humidity cannot be calculated an error value of 999 is returned.
1/4/00
List of Parameters Passed to Function:
db     = dry bulb temperature in degrees C
dpt    = dew point temperature in degrees C

Variable Returned
rh    = relative humidity in %, or error value of 999

\***************************************************************************/

int calc_humidity(float db, float dpt)  /* Function to find relative humidity */
{
	/* Local variables */
	double c1 = -5.6745359e+3, c2 = -0.51523058, c3 = -9.6778430e-3, c4 = 6.2215701e-7,
		c5 = 2.0747825e-9, c6 = -9.484024e-13, c7 = 4.1635019, c8 = -5.8002206e+3,
		c9 = -5.516256, c10 = -4.8640239e-2, c11 = 4.1764768e-5, c12 = -1.4452093e-8,
		c13 = 6.5459673;
	double arg, t, pres, pres_dew;
	int rh;

	if (db > 90.0 || dpt > 90.0 || dpt > db)    /* Check for valid input data */
		rh = 999;                   /* Missing data value */
	else
	{  /* Find water vapor saturation pressure (kPa) for drybulb temperature */
		t = db + 273.15;           /* Absolute temperature (deg K) */
		if (db < 0.0)             /* Use equation 3 if drybulb less than 0 */
		{
			arg = c1 / t + c2 + c3*t + c4*pow(t, 2.0) + c5*pow(t, 3.0) +
				c6*pow(t, 4.0) + c7*log(t);
			pres = exp(arg);
		}
		else
		{
			arg = c8 / t + c9 + c10*t + c11*pow(t, 2.0) + c12*pow(t, 3.0) +
				c13*log(t);
			pres = exp(arg);
		}
		/* Find water vapor saturation pressure (kPa) for dewpoint temperature */
		t = dpt + 273.15;           /* Absolute temperature (deg K) */
		if (dpt < 0.0)             /* Use equation 3 if dewpoint less than 0 */
		{
			arg = c1 / t + c2 + c3*t + c4*pow(t, 2.0) + c5*pow(t, 3.0) +
				c6*pow(t, 4.0) + c7*log(t);
			pres_dew = exp(arg);
		}
		else
		{
			arg = c8 / t + c9 + c10*t + c11*pow(t, 2.0) + c12*pow(t, 3.0) +
				c13*log(t);
			pres_dew = exp(arg);
		}

		rh = (int)(100.0*pres_dew / pres + 0.5); /* Relative humidity */
	}
	return(rh);
}

float calc_dewpt(float db, float rh)  /* Function to find dewpoint temperature */
{
	/* This function calculates the dewpoint temperature(C) based on the drybulb
	temperature(C) and the relative humidity(%).  It uses equations and
	procedures presented in the 1993 ASHRAE Fundamentals Handbook, p6.7-10.
	If dewpoint cannot be calculated an error value of 99.9 is returned.

	List of Parameters Passed to Function:
	db     = dry bulb temperature in degrees C
	rh     = relative humidity in %

	Variable Returned
	dpt    = dew point temperature in degrees C, or error value of 99.9   */

	/* Local variables */
	double c1 = -5.6745359e+3, c2 = -0.51523058, c3 = -9.6778430e-3, c4 = 6.2215701e-7,
		c5 = 2.0747825e-9, c6 = -9.484024e-13, c7 = 4.1635019, c8 = -5.8002206e+3,
		c9 = -5.516256, c10 = -4.8640239e-2, c11 = 4.1764768e-5, c12 = -1.4452093e-8,
		c13 = 6.5459673, c14 = 6.54, c15 = 14.526, c16 = 0.7389, c17 = 0.09486,
		c18 = 0.4569;
	double arg, t, pres, pres_dew, pta, ptb, ptc;
	float dpt =0;

	if (db > 90.0 || rh > 100.0 || rh < 1.0)    /* Check for valid input data */
		dpt = (float)99.9;                   /* Missing data value */
	else
	{  /* Find water vapor saturation pressure (kPa) for drybulb temperature */
		t = db + 273.15;           /* Absolute temperature (deg K) */
		if (db < 0.0)             /* Use equation 3 if drybulb less than 0 */
		{
			arg = c1 / t + c2 + c3*t + c4*pow(t, 2.0) + c5*pow(t, 3.0) +
				c6*pow(t, 4.0) + c7*log(t);
			pres = exp(arg);
		}
		else
		{
			arg = c8 / t + c9 + c10*t + c11*pow(t, 2.0) + c12*pow(t, 3.0) +
				c13*log(t);
			pres = exp(arg);
		}
		pres = pres*rh / 100.0;      /* Partial pressure (kPa) of water vapor */
		arg = log(pres);
		if (db >= 0.0)            /* Use equation 35 from ASHRAE */
			dpt = (float)(c14 + c15*arg + c16*pow(arg, 2.0) + c17*pow(arg, 3.0) +
			c18*pow(pres, 0.1984));
		if (db < 0.0 || dpt < 0)  /* Use eqn 36 if drybulb or dewpoint < 0 */
			dpt = (float)(6.09 + 12.608*arg + 0.4959*arg*arg);

		/* For dewpoint temperatures below -20C, check to see that the dewpoint
		temperature gives the correct vapor saturation pressure.  If not, iterate
		the correct dew point temperature */

		if (dpt < -20.0)
		{
			t = dpt + 273.15;       /* Absolute temperature (deg K) */
			arg = c1 / t + c2 + c3*t + c4*pow(t, 2.0) + c5*pow(t, 3.0) +
				c6*pow(t, 4.0) + c7*log(t);
			pres_dew = exp(arg);
			if (pres < pres_dew)   /* Set initial iteration points */
			{
				pta = t - 10.0;
				ptb = t;
				ptc = (pta + ptb) / 2.0;
			}
			else
			{
				pta = t;
				ptb = t + 10.0;
				ptc = (pta + ptb) / 2.0;
			}
			while (fabs(pres - pres_dew) > 0.00001 && fabs(pta - ptb) > 0.05)
			{
				dpt = (float)(ptc - 273.15);
				t = ptc;
				arg = c1 / t + c2 + c3*t + c4*pow(t, 2.0) + c5*pow(t, 3.0) +
					c6*pow(t, 4.0) + c7*log(t);
				pres_dew = exp(arg);
				if (pres < pres_dew)   /* Reset iteration points */
				{
					ptb = ptc;
					ptc = (pta + ptb) / 2.0;
				}
				else
				{
					pta = ptc;
					ptc = (pta + ptb) / 2.0;
				}
			}
		}
	}
	return(dpt);
}


double calc_twet(double T, double RH, double P)
{
	//	function [Twet] = calctwet(T, RH, P)
	//% calculate wet bulb temperature from T (dry bulb, 'C), RH (%), Pressure
	//% (mbar)
	//% see http://www.ejournal.unam.mx/atm/Vol07-3/ATM07304.pdf for eqns.

	/*
	Mike Wagner:
	There is a units error here! The original reference specifies that pressure should be provided in
	hPa (hectoPascals), which is equivalent with millibar. However, the units SHOULD BE in kPa, or mbar/10.
	Correct for the units issue here.

	IMPACT:
	This subroutine has been returning wet bulb temperatures much too high. This could adversely affect any
	model that calls the method and whose performance is sensitive to the wet bulb temperature.
	*/
	volatile double Pkpa = P / 10.;	//Correct for units problem

	//volatile double Twet = T*0.7;// initial guess
	volatile double Twet = T - 5.;	//Initial guess [mjw -- negative values of T were causing problems here]

	//[mjw] Use a bisection method to solve for Twet. The previous iteration method is unstable.
	bool
		hiflag = false,
		lowflag = false;
	double
		hival = 0, lowval = 0, err;
	const double tol = 0.05;

	int i = 0;
	while (i++ < 250)
	{
		err = exp((21.3 * Twet + 494.41) / (Twet + 273.15)) - RH / 100 * exp((21.3 * T + 494.41) / (T + 273.15)) - (6.53*10e-4) * Pkpa * (T - Twet);
		//double G = exp( (21.3 * Twet + 494.41) / (Twet+273.15) ) * ( (21.4 * (Twet+273.15) - (21.4 * Twet+494.41)) / pow(Twet+273.15, 2) ) + 6.53*10e-4 * Pkpa * Twet;
		if (err < 0.){
			lowval = Twet;
			lowflag = true;
		}
		else if (err > 0.){
			hival = Twet;
			hiflag = true;
		}

		if (fabs(err) < tol) break;

		//If the error is still too high, guess new values
		if (hiflag && lowflag){
			//Bisect
			Twet = (hival + lowval) / 2.;
		}
		else if (hiflag){
			//A lower bound hasn't yet been found. Try decreasing by 5 C
			Twet += -5;
		}
		else if (lowflag){
			//An upper bound hasn't yet been found. Bisect the current Twet and the Tdry
			Twet = (Twet + T) / 2.;
		}
		else{
			//Neither flags have been set. Guess a lower temp.
			Twet += -5.;
		}

	}

	if (Twet != Twet) // check for NaN
	{
		/*
		from biopower, Jennie Jorgenson:
		For estimating the dew point (first line of code), I used this very simple relation from wikipedia: http://en.wikipedia.org/wiki/Dew_point#Simple_approximation
		The second line is from a slightly sketchier source (http://www.theweatherprediction.com/habyhints/170/), meteorologist Jeff Haby. His procedure is for temperatures in F.
		*/

		double dp_est = T - ((1 - RH / 100) / 0.05);
		Twet = T - ((T - dp_est) / 3.0);
	}

	return Twet;
}

void weather_header::reset()
{
	location = city = state = country = source = description = url = "";
	hasunits = false;
	tz = lat = lon = elev = std::numeric_limits<double>::quiet_NaN();
}

void weather_record::reset()
{
	year = month = day = hour = 0;
	minute = std::numeric_limits<double>::quiet_NaN();
	gh = dn = df = poa = wspd = wdir = std::numeric_limits<double>::quiet_NaN();
	tdry = twet = tdew = rhum = pres = snow = alb =  aod = std::numeric_limits<double>::quiet_NaN();
}



#define NBUF 2048


weatherfile::weatherfile()
{
	reset();
}

weatherfile::weatherfile(const std::string &file, bool header_only)
{
	reset();
	m_ok = open(file, header_only);
}

weatherfile::~weatherfile()
{
	// virtual destructor, nothing to do..
}

void weatherfile::reset()
{
	m_startSec = m_stepSec = m_nRecords = 0;

	m_message.clear();
	m_ok = false;
	m_type = INVALID;
	m_startYear = 1900;
	m_time = 0;
	m_index = 0;

	m_type = INVALID;
	m_file.clear();
	m_startYear = 1900;

	m_hdr.reset();
	//m_rec.reset();
}


int weatherfile::type()
{
	return m_type;
}

std::string weatherfile::filename()
{
	return m_file;
}

bool weatherfile::open(const std::string &file, bool header_only)
{
	if (file.empty())
	{
		m_message = "no file name given to weather file reader";
		return false;
	}

	if (cmp_ext(file.c_str(), "tm2") || cmp_ext(file.c_str(), "tmy2"))
		m_type = TMY2;
	else if (cmp_ext(file.c_str(), "tm3") || cmp_ext(file.c_str(), "tmy3"))
		m_type = TMY3;
	else if (cmp_ext(file.c_str(), "csv"))
		m_type = WFCSV;
	else if (cmp_ext(file.c_str(), "epw"))
		m_type = EPW;
	else if (cmp_ext(file.c_str(), "smw"))
		m_type = SMW;
	else
	{
		m_message = "could not detect weather data file format from file extension (.csv,.tm2,.tm2,.epw)";
		return false;
	}

	std::string buf, buf1;
	std::ifstream ifs(file);

	if ( !ifs.is_open() )
	{
		m_message = "could not open file for reading: " + file;
		m_type = INVALID;
		return false;
	}

	if (m_type == WFCSV)
	{
		// if we opened a csv file, it could be SAM/WFCSV format or TMY3
		// try to autodetect a TMY3
		getline(ifs, buf);
		getline(ifs, buf1);
		int ncols = split(buf).size();
		int ncols1 = split(buf1).size();

		if (ncols == 7 && (ncols1 == 68 || ncols1 == 71))
			m_type = TMY3;

		ifs.seekg(0);
	}


	m_startYear = 1900;
	m_time = 1800;

	/* read header information */
	if (m_type == TMY2)
	{
		/*  93037 COLORADO_SPRINGS       CO  -7 N 38 49 W 104 43  1881 */
		char slat[10], slon[10];
		char pl[256], pc[256], ps[256];
		int dlat, mlat, dlon, mlon, ielv;

		getline(ifs, buf);
		sscanf(buf.c_str(),
		       "%s %s %s %lg %s %d %d %s %d %d %d",
			pl, pc, ps,
			&m_hdr.tz,
			slat, &dlat, &mlat,
			slon, &dlon, &mlon,
			&ielv);

		m_hdr.lat = conv_deg_min_sec(dlat, mlat, 0, slat[0]);
		m_hdr.lon = conv_deg_min_sec(dlon, mlon, 0, slon[0]);
		m_hdr.location = pl;
		m_hdr.city = pc;
		m_hdr.state = ps;
		m_hdr.elev = ielv;
		m_startSec = 1800;
		m_stepSec  = 3600;
		m_nRecords = 8760;
	}
	else if (m_type == TMY3)
	{
		/*  724699,"BROOMFIELD/JEFFCO [BOULDER - SURFRAD]",CO,-7.0,40.130,-105.240,1689 */
	  	getline(ifs, buf);
		auto cols = split(buf);
		if (cols.size() != 7)
		{
			m_message = "invalid TMY3 header: must contain 7 fields.  station,city,state,tz,lat,lon,elev";
			m_ok = false;
			return false;
		}

		m_hdr.location = cols[0];
		m_hdr.city = cols[1];
		m_hdr.state = cols[2];
		m_hdr.tz = stof(cols[3]);
		m_hdr.lat = stof(cols[4]);
		m_hdr.lon = stof(cols[5]);
		m_hdr.elev = stof(cols[6]);
		
		m_startSec = 1800;
		m_stepSec  = 3600;
		m_nRecords = 8760;

		getline(ifs, buf); // skip over labels line
	}
	else if (m_type == EPW)
	{
		/*  LOCATION,Cairo Intl Airport,Al Qahirah,EGY,ETMY,623660,30.13,31.40,2.0,74.0 */
		/*  LOCATION,Alice Springs Airport,NT,AUS,RMY,943260,-23.80,133.88,9.5,547.0 */

	  	getline(ifs, buf);
		auto cols = split(buf);

		if (cols.size() != 10)
		{
			m_message = "invalid EPW header: must contain 10 fields. LOCATION,city,state,country,source,station,lat,lon,tz,elev";
			m_ok = false;
			return false;
		}

		m_hdr.city = cols[1];
		m_hdr.state = cols[2];
		m_hdr.country = cols[3];
		m_hdr.source = cols[4];
		m_hdr.location = cols[5];
		m_hdr.lat = stof(cols[6]);
		m_hdr.lon = stof(cols[7]);
		m_hdr.tz = stof(cols[8]);
		m_hdr.elev = stof(cols[9]);

		/* skip over excess header lines */

		getline(ifs, buf);  // DESIGN CONDITIONS
		getline(ifs, buf);  // TYPICAL/EXTREME PERIODS
		getline(ifs, buf);  // GROUND TEMPERATURES
		getline(ifs, buf);  // HOLIDAY/DAYLIGHT SAVINGS
		getline(ifs, buf);  // COMMENTS 1
		getline(ifs, buf);  // COMMENTS 2
		getline(ifs, buf);  // DATA PERIODS
		
		m_startSec = 1800;
		m_stepSec  = 3600;
		m_nRecords = 8760;

	}
	else if (m_type == SMW)
	{
	  	getline(ifs, buf);
		auto cols = split(buf);

		if (cols.size() != 10)
		{
			m_message = "invalid SMW header format, 10 fields required";
			m_ok = false;
			return false;
		}

		m_hdr.location = cols[0];
		m_hdr.city = cols[1];
		m_hdr.state = cols[2];

		m_hdr.tz = stof(cols[3]);
		m_hdr.lat = stof(cols[4]);
		m_hdr.lon = stof(cols[5]);
		m_hdr.elev = stof(cols[6]);
		m_stepSec = (size_t)stof(cols[7]); // time step in seconds
		m_startYear = stoi(cols[8]);
		char *p = const_cast<char *>(cols[9].c_str());
		
		double start_hour = 0;
		double start_min = 30;
		double start_sec = 0;

		start_hour = atoi(p);

		p = strchr(p, ':');
		if (p && *p) p++;
		if (p && *p) start_min = atoi(p);

		p = strchr(p, ':');
		if (p && *p) p++;
		if (p && *p) start_sec = atoi(p);


		if (!header_only)
		{
			m_time = start_hour * 3600 + start_min * 60 + start_sec;
			m_startSec = (size_t)m_time;

			m_nRecords = 0;
			while (getline(ifs, buf))
				m_nRecords++;

			ifs.seekg(0);
			getline(ifs, buf);
			
			if ( m_nRecords%8784==0 )
			{
				// Check if the weather file contains a leap day
				// if so, exit out with an error 
				m_message = "could not determine timestep in CSV weather file. Does the file contain a leap day?";
				m_ok = false;
				return false;
			}
		}
	}
	else if (m_type == WFCSV)
	{
	  	getline(ifs, buf);
		auto cols = split(buf);
		int ncols = cols.size();
	  	getline(ifs, buf1);
		auto cols1 = split(buf1);
		int ncols1 = split(buf1).size();

		int hdr_step_sec = -1;

		if (ncols != ncols1)
		{
			m_message = "first two header lines must have same number of columns";
			return false;
		}

		for (size_t i = 0; (int)i<ncols; i++)
		{
		  	const std::string name = util::lower_case(trimboth(cols[i]));
			const std::string value = trimboth(cols1[i]);

			if (name == "lat" || name == "latitude")
			{
				m_hdr.lat = stof(value);
			}
			else if (name == "lon" || name == "long" || name == "longitude" || name == "lng")
			{
				m_hdr.lon = stof(value);
			}
			else if (name == "tz" || name == "timezone" || name == "time zone")
			{
				m_hdr.tz = stof(value);
			}
			else if (name == "el" || name == "elev" || name == "elevation" || name == "site elevation")
			{
				m_hdr.elev = stof(value);
			}
			else if (name == "year")
			{
			  	m_startYear = stoi(value);
			}
			else if (name == "id" || name == "location" || name == "location id" || name == "station" || name == "station id" || name == "wban" || name == "wban#")
			{
				m_hdr.location = value;
			}
			else if (name == "city")
			{
				m_hdr.city = value;
			}
			else if (name == "state" || name == "province" || name == "region")
			{
				m_hdr.state = value;
			}
			else if (name == "country")
			{
				m_hdr.country = value;
			}
			else if (name == "source" || name == "src")
			{
				m_hdr.source = value;
			}
			else if (name == "description" || name == "desc")
			{
				m_hdr.description = value;
			}
			else if (name == "url")
			{
				m_hdr.url = value;
			}
			else if (name == "hasunits" || name == "units")
			{
			  	m_hdr.hasunits = (util::lower_case(value) == "yes" || stoi(value) != 0);
			}
			else if (name == "step")
			{
			  	hdr_step_sec = stoi(value);
			}
		}

		if ( !std::isfinite(m_hdr.lat) || !std::isfinite(m_hdr.lon) )
		{
			m_message = "latitude and longitude required but not specified";
			return false;
		}

		if (!std::isfinite(m_hdr.tz))
		{
			m_message = "time zone required but not specified";
			return false;
		}


		// only scan to determine # of records
		// if we actually plan to read in the whole file
		if (!header_only)
		{
			m_startSec = 1800;
			m_stepSec = 3600;
			m_nRecords = 8760;

			getline(ifs, buf);  // col names
			if (m_hdr.hasunits)
			  	getline(ifs, buf);  // col units

			m_nRecords = 0; // figure out how many records there are
			while (getline(ifs, buf) && buf.length() > 0)
				m_nRecords++;

			// reposition to where we were
			ifs.seekg(0);
			getline(ifs, buf);  // header names
			getline(ifs, buf);  // header values

			// now determine timestep as best as possible
			int nmult = (int)m_nRecords / 8760;
			// divide by zero error 2/20/19
			if (nmult <= 0)
			{
				m_message = "could not determine number of records in CSV weather file";
				m_ok = false;
				return false;
			}



			if (hdr_step_sec > 0)
			{  // if explicitly specified in header?
				m_stepSec = hdr_step_sec;
				m_startSec = m_stepSec / 2;
			}
			else if (nmult * 8760 == (int)m_nRecords)
			{
				// multiple of 8760 records: assume 1 year of data
				m_stepSec = 3600 / nmult;
				m_startSec = m_stepSec / 2;
			}
			else if ( m_nRecords%8784==0 )
			{ 
				// Check if the weather file contains a leap day
				// if so, correct the number of nrecords 
				m_nRecords = m_nRecords/8784*8760;
				nmult = (int)m_nRecords/8760;
				m_stepSec = 3600 / nmult;
				m_startSec = m_stepSec / 2;
			}
			else
			{
				m_message = "could not determine timestep in CSV weather file";
				m_ok = false;
				return false;
			}
		}

	}
	else
	{
		m_message = "could not detect file format";
		return false;
	}

	if (header_only)
	{
		return true;
	}

	// preallocate memory for data
	for (size_t i = 0; i<_MAXCOL_; i++)
	{
		m_columns[i].index = -1;
		m_columns[i].data.resize(m_nRecords, std::numeric_limits<float>::quiet_NaN());
	}

	if (m_type == WFCSV)
	{
		// if it's a WFCSV format file, we need to determine which columns of data exist
	  	getline(ifs, buf);  // read column names
		if (ifs.eof())
		{
			m_message = "could not read column names";
			return false;
		}

		auto cols = split(buf);
		int ncols = cols.size();

		if (m_hdr.hasunits)
		{
		  	getline(ifs, buf);  // read column units
			if (ifs.eof())
			{
				m_message = "could not read column units";
				return false;
			}
			auto cols1 = split(buf1);
			int ncols1 = cols1.size();

			if (ncols != ncols1) {
				m_message = "column names and units must have the same number of fields";
				return false;
			}
		}

		// determine columns
		for (int i = 0; i<ncols; i++)
		{
		    	const std::string name = trimboth(cols[i]);
			if (name.length() > 0)
			{
			  	std::string lowname = util::lower_case(name);

				if (lowname == "yr" || lowname == "year") m_columns[YEAR].index = i;
				else if (lowname == "mo" || lowname == "month") m_columns[MONTH].index = i;
				else if (lowname == "day") m_columns[DAY].index = i;
				else if (lowname == "hour" || lowname == "hr") m_columns[HOUR].index = i;
				else if (lowname == "min" || lowname == "minute") m_columns[MINUTE].index = i;
				else if (lowname == "ghi" || lowname == "gh" || lowname == "global" || lowname == "global horizontal" || lowname == "global horizontal irradiance") m_columns[GHI].index = i;
				else if (lowname == "dni" || lowname == "dn" || lowname == "beam" || lowname == "direct normal" || lowname == "direct normal irradiance") m_columns[DNI].index = i;
				else if (lowname == "dhi" || lowname == "df" || lowname == "diffuse" || lowname == "diffuse horizontal" || lowname == "diffuse horizontal irradiance") m_columns[DHI].index = i;
				else if (lowname == "poa" || lowname == "pa" || lowname == "plane" || lowname == "plane of array" || lowname == "plane of array irradiance") m_columns[POA].index = i;
				else if (lowname == "tdry" || lowname == "dry bulb" || lowname == "dry bulb temp" || lowname == "temperature" || lowname == "ambient" || lowname == "ambient temp") m_columns[TDRY].index = i;
				else if (lowname == "twet" || lowname == "wet bulb" || lowname == "wet bulb temperature") m_columns[TWET].index = i;
				else if (lowname == "tdew" || lowname == "dew point" || lowname == "dew point temperature") m_columns[TDEW].index = i;
				else if (lowname == "wspd" || lowname == "wind speed") m_columns[WSPD].index = i;
				else if (lowname == "wdir" || lowname == "wind direction") m_columns[WDIR].index = i;
				else if (lowname == "rh" || lowname == "rhum" || lowname == "relative humidity" || lowname == "humidity") m_columns[RH].index = i;
				else if (lowname == "pres" || lowname == "pressure") m_columns[PRES].index = i;
				else if (lowname == "snow" || lowname == "snow cover" || lowname == "snow depth") m_columns[SNOW].index = i;
				else if (lowname == "alb" || lowname == "albedo") m_columns[ALB].index = i;
				else if (lowname == "aod" || lowname == "aerosol" || lowname == "aerosol optical depth") m_columns[AOD].index = i;
			}
		}
	}
	else if ( m_type == TMY2 )
	{
		// indicate which columns are available in TMY2 files
		m_columns[YEAR].index 
			= m_columns[MONTH].index
			= m_columns[DAY].index
			= m_columns[HOUR].index
			= m_columns[GHI].index
			= m_columns[DNI].index
			= m_columns[DHI].index
			= m_columns[TDRY].index
			= m_columns[TDEW].index
			= m_columns[WSPD].index
			= m_columns[WDIR].index
			= m_columns[RH].index
			= m_columns[PRES].index
			= m_columns[SNOW].index
			= 1;
	}
	else if ( m_type == TMY3 )
	{
		// indicate which columns are available in TMY3 files
		m_columns[YEAR].index 
			= m_columns[MONTH].index
			= m_columns[DAY].index
			= m_columns[HOUR].index
			= m_columns[GHI].index
			= m_columns[DNI].index
			= m_columns[DHI].index
			= m_columns[TDRY].index
			= m_columns[TDEW].index
			= m_columns[WSPD].index
			= m_columns[WDIR].index
			= m_columns[RH].index
			= m_columns[PRES].index
			= m_columns[ALB].index
			= 1;
	}
	else if ( m_type == EPW )
	{
		// indicate which columns are available in EPW files
		m_columns[YEAR].index 
			= m_columns[MONTH].index
			= m_columns[DAY].index
			= m_columns[HOUR].index
			= m_columns[GHI].index
			= m_columns[DNI].index
			= m_columns[DHI].index
			= m_columns[TDRY].index
			= m_columns[TWET].index
			= m_columns[WSPD].index
			= m_columns[WDIR].index
			= m_columns[RH].index
			= m_columns[PRES].index
			= m_columns[SNOW].index
			= 1;
	}
	else if ( m_type == SMW )
	{				
		// indicate which columns are available in SMW files
		m_columns[YEAR].index 
			= m_columns[MONTH].index
			= m_columns[DAY].index
			= m_columns[HOUR].index
			= m_columns[GHI].index
			= m_columns[DNI].index
			= m_columns[DHI].index
			= m_columns[TDRY].index
			= m_columns[TWET].index
			= m_columns[WSPD].index
			= m_columns[WDIR].index
			= m_columns[RH].index
			= m_columns[PRES].index
			= m_columns[SNOW].index
			= 1;
	}


	// by default, subtract 1 from hour of TMY3 files to switch
	// from 1-24 standard to 0-23
	int tmy3_hour_shift = 1;
	int n_leap_data_removed = 0;
	
	for (int i = 0; i<(int)m_nRecords; i++)
	{
		if (m_type == TMY2)
		{
			
			int yr, mn, dy, hr, ethor, etdn;
			int d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21;      /* which of these are used? d3, d10, d15 & d20 */
			int u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15, u16, u17, u18, u19, u20, u21;  /* are any of these ever used?? */
			int w1, w2, w3, w4, w5, w6, w7, w8, w9, w10;
			char f1[2], f2[2], f3[2], f4[2], f5[2], f6[2], f7[2], f8[2], f9[2], f10[2], f11[2], f12[2], f13[2], f14[2], f15[2], f16[2], f17[2], f18[2], f19[2], f20[2], f21[2];

			int nread = 0;

			for(;;)
			{
			  	getline(ifs, buf);
				nread = sscanf(buf.c_str(),
					"%2d%2d%2d%2d"
					"%4d%4d"
					"%4d%1s%1d%4d%1s%1d%4d%1s%1d%4d%1s%1d%4d%1s%1d%4d%1s%1d%4d%1s%1d"
					"%2d%1s%1d%2d%1s%1d%4d%1s%1d%4d%1s%1d%3d%1s%1d%4d%1s%1d%3d%1s%1d"
					"%3d%1s%1d%4d%1s%1d%5d%1s%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%3d%1s%1d%3d%1s%1d%3d%1s%1d%2d%1s%1d\n",
					&yr, &mn, &dy, &hr,
					&ethor, /* extraterrestrial horizontal radiation */
					&etdn, /* extraterrestrial direct normal radiation */
					&d1, f1, &u1, /* GH data value 0-1415 Wh/m2, Source, Uncertainty */
					&d2, f2, &u2, /* DN data value 0-1200 Wh/m2, Source, Uncertainty */
					&d3, f3, &u3, /* DF data value 0-700 Wh/m2, Source, Uncertainty */
					&d4, f4, &u4, /* GH illum data value, Source, Uncertainty */
					&d5, f5, &u5, /* DN illum data value, Source, Uncertainty */
					&d6, f6, &u6, /* DF illum data value, Source, Uncertainty */
					&d7, f7, &u7, /* Zenith illum data value, Source, Uncertainty */
					&d8, f8, &u8, /* Total sky cover */
					&d9, f9, &u9, /* opaque sky cover */
					&d10, f10, &u10, /* dry bulb temp -500 to 500 = -50.0 to 50.0 'C */
					&d11, f11, &u11, /* dew point temp -600 to 300 = -60.0 to 30.0 'C */
					&d12, f12, &u12, /* relative humidity 0-100 */
					&d13, f13, &u13, /* pressure millibars */
					&d14, f14, &u14, /* wind direction */
					&d15, f15, &u15, // wind speed 0 to 400 = 0.0 to 40.0 m/s
					&d16, f16, &u16, // visibility
					&d17, f17, &u17, // ceiling height
					&w1, &w2, &w3, &w4, &w5, &w6, &w7, &w8, &w9, &w10, // present weather
					&d18, f18, &u18, // precipitable water
					&d19, f19, &u19, // aerosol optical depth
					&d20, f20, &u20, // snow depth 0-150 cm
					&d21, f21, &u21); // days since last snowfall 0-88

				if ( mn == 2 && dy == 29 )
				{
					// skip data lines for february 29th if they exist in the file
					n_leap_data_removed++;
					continue;
				}

				m_columns[YEAR].data[i] = (float)yr + 1900;
				m_columns[MONTH].data[i] = (float)mn;
				m_columns[DAY].data[i] = (float)dy;
				m_columns[HOUR].data[i] = (float)hr - 1;  // hour goes 0-23, not 1-24
				m_columns[MINUTE].data[i] = 30;
				m_columns[GHI].data[i] = (float)(d1*1.0);
				m_columns[DNI].data[i] = (float)d2;           /* Direct radiation */
				m_columns[DHI].data[i] = (float)d3;           /* Diffuse radiation */
				m_columns[POA].data[i] = (float)(-999);       /* No POA in TMY2 */
				m_columns[TDRY].data[i] = (float)(d10 / 10.0);       /* Ambient dry bulb temperature(C) */
				m_columns[TDEW].data[i] = (float)(d11 / 10.0); /* dew point temp */
				m_columns[WSPD].data[i] = (float)(d15 / 10.0);       /* Wind speed(m/s) */
				m_columns[WDIR].data[i] = (float)d14; /* wind dir */
				m_columns[RH].data[i] = (float)d12;
				m_columns[PRES].data[i] = (float)d13;
				m_columns[SNOW].data[i] = (float)d20;
				m_columns[ALB].data[i] = -999; /* no albedo in TMY2 */
				m_columns[AOD].data[i] = -999; /* no AOD in TMY2 */
				m_columns[TWET].data[i] 
					= (float)calc_twet( 
					(double)m_columns[TDRY].data[i],
					(double)m_columns[RH].data[i],
					(double)m_columns[PRES].data[i] ); /* must calculate wet bulb */

				break;
			}


			if ( nread != 79 || ifs.eof() )
			{
				m_message = "TMY2: data line does not have at exactly 79 characters at record " + util::to_string(i);
				return false;
			}

		}
		else if (m_type == TMY3)
		{
			for(;;)
			{
				getline(ifs, buf);
				auto cols = split(buf);
				if (cols.size() < 68)
				{
					m_message = "TMY3: data line does not have at least 68 fields at record " + util::to_string(i);
					return false;
				}

				const char *p = cols[0].c_str();

				int month = stoi(p);
				p = strchr(p, '/');
				if (!p)
				{
					m_message = "TMY3: invalid date format at record " + util::to_string(i);
					return false;
				}
				p++;
				int day = stoi(p);
				p = strchr(p, '/');
				if (!p)
				{
					m_message = "TMY3: invalid date format at record " + util::to_string(i);
					return false;
				}
				p++;
				int year = stoi(p);

				int hour = stoi(cols[1]) - tmy3_hour_shift;  // hour goes 0-23, not 1-24
				if (i == 0 && hour < 0)
				{
					// this was a TMY3 file but with hours going 0-23 (against the tmy3 spec)
					// handle it anyway by NOT subtracting from the hour to convert from 1-24
					tmy3_hour_shift = 0;
					hour = 0;
				}

				if ( month == 2 && day == 29 )
				{
					n_leap_data_removed++;
					continue;
				}

				m_columns[YEAR].data[i] = (float)year;
				m_columns[MONTH].data[i] = (float)month;
				m_columns[DAY].data[i] = (float)day;
				m_columns[HOUR].data[i] = (float)hour;
				m_columns[MINUTE].data[i] = 30;

				m_columns[GHI].data[i] = (float)stof(cols[4]);
				m_columns[DNI].data[i] = (float)stof(cols[7]);
				m_columns[DHI].data[i] = (float)stof(cols[10]);
				m_columns[POA].data[i] = (float)(-999);       /* No POA in TMY3 */

				m_columns[TDRY].data[i] = (float)stof(cols[31]);
				m_columns[TDEW].data[i] = (float)stof(cols[34]);
				
				m_columns[WSPD].data[i] = (float)stof(cols[46]);
				m_columns[WDIR].data[i] = (float)stof(cols[43]);

				m_columns[RH].data[i] = (float)stof(cols[37]);
				m_columns[PRES].data[i] = (float)stof(cols[40]);
				m_columns[SNOW].data[i] = -999.0; // no snowfall in TMY3
				m_columns[ALB].data[i] = (float)stof(cols[61]);
				m_columns[AOD].data[i] = -999; /* no AOD in TMY3 */

				m_columns[TWET].data[i] 
					= (float)calc_twet( 
					(double)m_columns[TDRY].data[i],
					(double)m_columns[RH].data[i],
					(double)m_columns[PRES].data[i]); /* must calculate wet bulb */

				break;
			}

			if (ifs.eof())
			{
				m_message = "TMY3: data line formatting error at record " + util::to_string(i);
				return false;
			}
		}
		else if (m_type == EPW)
		{
			for(;;)
			{
			  	getline(ifs, buf);
				auto cols = split(buf);

				if (cols.size() < 32)
				{
					m_message = "EPW: data line does not have at least 32 fields at record " + util::to_string(i);
					return false;
				}

				int month = stoi(cols[1]);
				int day = stoi(cols[2] );

				if ( month == 2 && day == 29 )
				{
					n_leap_data_removed++;
					continue;
				}

				m_columns[YEAR].data[i] = (float)stoi(cols[0]);
				m_columns[MONTH].data[i] = (float)stoi(cols[1]);
				m_columns[DAY].data[i] = (float)stoi(cols[2]);
				m_columns[HOUR].data[i] = (float)stoi(cols[3]) - 1;  // hour goes 0-23, not 1-24;
				m_columns[MINUTE].data[i] = 30;

				m_columns[GHI].data[i] = (float)stof(cols[13]);
				m_columns[DNI].data[i] = (float)stof(cols[14]);
				m_columns[DHI].data[i] = (float)stof(cols[15]);
				m_columns[POA].data[i] = (float)(-999);       /* No POA in EPW */

				m_columns[WSPD].data[i] = (float)stof(cols[21]);
				m_columns[WDIR].data[i] = (float)stof(cols[20]);

				m_columns[TDRY].data[i] = (float)stof(cols[6]);
				m_columns[TWET].data[i] = (float)stof(cols[7]);

				m_columns[RH].data[i] = (float)stof(cols[8]);
				m_columns[PRES].data[i] = (float)(stof(cols[9]) * 0.01); /* convert Pa in to mbar */
				m_columns[SNOW].data[i] = (float)stof(cols[30]); // snowfall
				m_columns[ALB].data[i] = -999; /* no albedo in EPW file */
				m_columns[AOD].data[i] = -999; /* no AOD in EPW */

				m_columns[TDEW].data[i] = (float)wiki_dew_calc(m_columns[TDRY].data[i], m_columns[RH].data[i]);

				break;
			}

			if (ifs.eof())
			{
				m_message = "EPW: data line formatting error at record " + util::to_string(i);
				return false;
			}
		}
		else if (m_type == SMW)
		{
			getline(ifs, buf);
			auto cols = split(buf);

			if (cols.size() < 12)
			{
				m_message = "SMW: data line does not have at least 12 fields at record " + util::to_string(i);
				return false;
			}

			double T = m_time;

			m_columns[YEAR].data[i] = (float)m_startYear; // start year
			m_columns[MONTH].data[i] = (float)util::month_of(T / 3600.0); // 1-12
			m_columns[DAY].data[i] = (float)util::day_of_month((int)m_columns[MONTH].data[i], T / 3600.0); // 1-nday
			m_columns[HOUR].data[i] = (float)(((int)(T / 3600.0)) % 24);  // hour goes 0-23, not 1-24;
			m_columns[MINUTE].data[i] = (float)fmod(T / 60.0, 60.0);      // minute goes 0-59

			m_time += m_stepSec; // increment by step

			m_columns[GHI].data[i] = (float)stof(cols[7]);
			m_columns[DNI].data[i] = (float)stof(cols[8]);
			m_columns[DHI].data[i] = (float)stof(cols[9]);
			m_columns[POA].data[i] = (double)(-999);       /* No POA in SMW */

			m_columns[WSPD].data[i] = (float)stof(cols[4]);
			m_columns[WDIR].data[i] = (float)stof(cols[5]);

			m_columns[TDRY].data[i] = (float)stof(cols[0]);
			m_columns[TDEW].data[i] = (float)stof(cols[1]);
			m_columns[TWET].data[i] = (float)stof(cols[2]);

			m_columns[RH].data[i] = (float)stof(cols[3]);
			m_columns[PRES].data[i] = (float)stof(cols[6]);
			m_columns[SNOW].data[i] = (float)stof(cols[11]);
			m_columns[ALB].data[i] = (float)stof(cols[10]);
			m_columns[AOD].data[i] = -999; /* no AOD in SMW */

			if (ifs.eof())
			{
				m_message = "SMW: data line formatting error at record " + util::to_string(i);
				return false;
			}
		}
		else if (m_type == WFCSV)
		{	

			for(;;)
			{
			  	getline(ifs, buf);
				buf = trimboth(buf);
				if (buf.length() == 0)
				{
					m_message = "CSV: data line formatting error at record " + util::to_string(i);
					return false;
				}

				auto cols = split(buf);
				int ncols = cols.size();
				for (size_t k = 0; k < _MAXCOL_; k++)
				{
					if (m_columns[k].index >= 0
						&& m_columns[k].index < ncols)
					{
						m_columns[k].data[i] = (float)stof(trimboth(cols[m_columns[k].index]));
					} 
				}

				if ( m_columns[MONTH].data[i] == 2
					&& m_columns[DAY].data[i] == 29 )
				{
					n_leap_data_removed++;
					continue;
				}
				else
					break;
			}


		}

	}

	if( n_leap_data_removed > 0 )
		m_message = util::format("Skipped %d data lines for February 29th (leap day).", n_leap_data_removed );

	if (m_type == WFCSV)
	{
		// special handling for certain columns that we can calculate from others
		// if the data doesn't exist

		if (m_columns[TWET].index < 0
			&& m_columns[TDRY].index >= 0
			&& m_columns[PRES].index >= 0
			&& m_columns[RH].index >= 0)
		{
			for (size_t i = 0; i<m_nRecords; i++)
				m_columns[TWET].data[i] = (float)calc_twet(m_columns[TDRY].data[i], m_columns[RH].data[i], m_columns[PRES].data[i]);
		}

		if (m_columns[TDEW].index < 0
			&& m_columns[TDRY].index >= 0
			&& m_columns[RH].index >= 0)
		{
			for (size_t i = 0; i<m_nRecords; i++)
				m_columns[TDEW].data[i] = (float)wiki_dew_calc(m_columns[TDRY].data[i], m_columns[RH].data[i]);
		}

		if (m_columns[YEAR].index < 0)
		{
			for (size_t i = 0; i<m_nRecords; i++)
				m_columns[YEAR].data[i] = (float)m_startYear;
		}

		if (m_columns[MONTH].index < 0
			&& m_stepSec == 3600 && m_nRecords == 8760)
		{
			for (size_t i = 0; i<m_nRecords; i++)
				m_columns[MONTH].data[i] = (float)util::month_of((double)i);
		}

		if (m_columns[DAY].index < 0
			&& m_stepSec == 3600 && m_nRecords == 8760)
		{
			for (size_t i = 0; i<m_nRecords; i++)
			{
				int month = util::month_of((double)i);
				m_columns[DAY].data[i] = (float)util::day_of_month(month, (double)i);
			}
		}

		if (m_columns[HOUR].index < 0
			&& m_stepSec == 3600 && m_nRecords == 8760)
		{
			for (size_t i = 0; i<m_nRecords; i++)
			{
				size_t day = i / 24;
				size_t start_of_day = day * 24;
				m_columns[HOUR].data[i] = (float)(i - start_of_day);
			}
		}

		if (m_columns[MINUTE].index < 0 && (int)m_columns[HOUR].data[1] == m_columns[HOUR].data[1])
		{
			for (size_t i = 0; i<m_nRecords; i++)
				m_columns[MINUTE].data[i] = (float)((m_stepSec / 2) / 60);
		}
        else if( m_columns[MINUTE].index < 0 )  //implies fractional hours are provided
        {
            for (size_t i = 0; i<m_nRecords; i++)
            {
                float hr = m_columns[HOUR].data[i];
                m_columns[MINUTE].data[i] = (float)((hr - (int)hr)*60.);
                m_columns[HOUR].data[i] = (float)(int)hr;
            }
        }
	}

	return true;
}

bool weatherfile::read( weather_record *r )
{
	if ( r && m_index < m_nRecords)
	{
		r->year = (int)m_columns[YEAR].data[m_index];
		r->month = (int)m_columns[MONTH].data[m_index];
		r->day = (int)m_columns[DAY].data[m_index];
		r->hour = (int)m_columns[HOUR].data[m_index];
		r->minute = m_columns[MINUTE].data[m_index];
		r->gh = m_columns[GHI].data[m_index];
		r->dn = m_columns[DNI].data[m_index];
		r->df = m_columns[DHI].data[m_index];
		r->poa = m_columns[POA].data[m_index];
		r->wspd = m_columns[WSPD].data[m_index];
		r->wdir = m_columns[WDIR].data[m_index];
		r->tdry = m_columns[TDRY].data[m_index];
		r->twet = m_columns[TWET].data[m_index];
		r->tdew = m_columns[TDEW].data[m_index];
		r->rhum = m_columns[RH].data[m_index];
		r->pres = m_columns[PRES].data[m_index];
		r->snow = m_columns[SNOW].data[m_index];
		r->alb = m_columns[ALB].data[m_index];
		r->aod = m_columns[AOD].data[m_index];

		m_index++;
		return true;
	}
	else
		return false;
}

bool weatherfile::has_data_column( size_t id )
{
	return m_columns[id].index >= 0;
}

bool weatherfile::convert_to_wfcsv( const std::string &input, const std::string &output )
{
	weatherfile wf( input );
	if ( !wf.ok() ) return false;

	util::stdfile fp( output, "w" );
	if ( !fp.ok() ) return false;

	weather_header hdr; 
	wf.header( &hdr );
	weather_record rec;

	if ( wf.type() == weatherfile::TMY2 )
	{
		fprintf(fp, "Source,Location ID,City,State,Country,Latitude,Longitude,Time Zone,Elevation\n");
		fprintf(fp, "TMY2,%s,%s,%s,USA,%.6lf,%.6lf,%lg,%lg\n", hdr.location.c_str(),
			normalize_city(hdr.city).c_str(), hdr.state.c_str(), hdr.lat, hdr.lon, hdr.tz, hdr.elev );
		fprintf(fp, "Year,Month,Day,Hour,GHI,DNI,DHI,Tdry,Tdew,RH,Pres,Wspd,Wdir,Snow Depth\n" );
		for( size_t i=0;i<8760;i++ )
		{
			if (!wf.read( &rec ) ) return false;
			fprintf(fp, "%d,%d,%d,%d,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
				rec.year, rec.month, rec.day, rec.hour,
				rec.gh, rec.dn, rec.df, rec.tdry, rec.tdew, rec.rhum, rec.pres, rec.wspd, rec.wdir, rec.snow );
		}
	}
	else if ( wf.type() == weatherfile::TMY3 )
	{
		fprintf(fp, "Source,Location ID,City,State,Country,Latitude,Longitude,Time Zone,Elevation\n");
		fprintf(fp, "TMY3,%s,%s,%s,USA,%.6lf,%.6lf,%lg,%lg\n", hdr.location.c_str(),
			normalize_city(hdr.city).c_str(), hdr.state.c_str(), hdr.lat, hdr.lon, hdr.tz, hdr.elev );
		fprintf(fp, "Year,Month,Day,Hour,GHI,DNI,DHI,Tdry,Tdew,RH,Pres,Wspd,Wdir,Albedo\n" );
		for( size_t i=0;i<8760;i++ )
		{
			if (!wf.read( &rec )) return false;
			fprintf(fp, "%d,%d,%d,%d,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
				rec.year, rec.month, rec.day, rec.hour,
				rec.gh, rec.dn, rec.df, rec.tdry, rec.tdew, rec.rhum, rec.pres, rec.wspd, rec.wdir, rec.alb );
		}
	}
	else if ( wf.type() == weatherfile::EPW )
	{
		fprintf(fp, "Source,Location ID,City,State,Country,Latitude,Longitude,Time Zone,Elevation\n");
		fprintf(fp, "EPW,%s,%s,%s,%s,%.6lf,%.6lf,%lg,%lg\n", hdr.location.c_str(),
			normalize_city(hdr.city).c_str(), hdr.state.c_str(), hdr.country.c_str(), hdr.lat, hdr.lon, hdr.tz, hdr.elev );
		fprintf(fp, "Year,Month,Day,Hour,GHI,DNI,DHI,Tdry,Twet,RH,Pres,Wspd,Wdir,Albedo\n" );
		for( size_t i=0;i<8760;i++ )
		{
			if (!wf.read( &rec )) return false;
			fprintf(fp, "%d,%d,%d,%d,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
				rec.year, rec.month, rec.day, rec.hour,
				rec.gh, rec.dn, rec.df, rec.tdry, rec.twet, rec.rhum, rec.pres, rec.wspd, rec.wdir, rec.alb );
		}
	}
	else if ( wf.type() == weatherfile::SMW )
	{
		fprintf(fp, "Source,Location ID,City,State,Latitude,Longitude,Time Zone,Elevation\n");
		fprintf(fp, "SMW,%s,%s,%s,%s,%.6lf,%.6lf,%lg,%lg\n", hdr.location.c_str(),
			normalize_city(hdr.city).c_str(), hdr.state.c_str(), hdr.country.c_str(), hdr.lat, hdr.lon, hdr.tz, hdr.elev );
		fprintf(fp, "Month,Day,Hour,GHI,DNI,DHI,Tdry,Twet,Tdew,RH,Pres,Wspd,Wdir,Snow,Albedo\n" );
		for( size_t i=0;i<8760;i++ )
		{
			if (!wf.read( &rec )) return false;
			fprintf(fp, "%d,%d,%d,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n",
				rec.month, rec.day, rec.hour,
				rec.gh, rec.dn, rec.df, rec.tdry, rec.twet, rec.tdew, rec.rhum, rec.pres, rec.wspd, rec.wdir, rec.snow, rec.alb );
		}
	}
	else
		return false;


	return true;

}

