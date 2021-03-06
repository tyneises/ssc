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

#include "core.h"

static var_info _cm_vtab_timeseq[] = 
{	
/*   VARTYPE           DATATYPE         NAME                         LABEL                              UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT,        SSC_NUMBER,      "start_time",                 "Start time",                     "seconds", "0=jan1st 12am",         "Time Sequence", "*",                       "MIN=0,MAX=31536000",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "end_time",                   "End time",                       "seconds", "0=jan1st 12am",         "Time Sequence", "*",                       "MIN=0,MAX=31536000",                     "" },
	{ SSC_INPUT,        SSC_NUMBER,      "time_step",                  "Time step",                      "seconds", "",                     "Time Sequence", "*",                       "MIN=1,MAX=3600",                         "" },

	{ SSC_OUTPUT,       SSC_ARRAY,       "time",                       "Time",                           "secs",   "0=jan1st 12am",        "Time",          "*",                       "",                                       "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "timehr",                     "HourTime",                       "hours",  "0=jan1st 12am",        "Time",          "*",                       "",                                       "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "month",                      "Month",                          "",       "1-12",                 "Time",          "*",                       "",                                       "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "day",                        "Day",                            "",       "1-{28,30,31}",         "Time",          "*",                       "",                                       "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "hour",                       "Hour",                           "",       "0-23",                 "Time",          "*",                       "",                                       "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "minute",                     "Minute",                         "",       "0-59",                 "Time",          "*",                       "",                                       "" },

var_info_invalid };

class cm_timeseq : public compute_module
{
private:
public:
	cm_timeseq()
	{
		add_var_info( _cm_vtab_timeseq );
	}

	void exec( ) throw( general_error )
	{
		double t_start = as_double("start_time");
		double t_end = as_double("end_time");
		double t_step = as_double("time_step"); // seconds

		size_t num_steps = check_timestep_seconds( t_start, t_end, t_step );

		ssc_number_t *time = allocate("time", num_steps);
		ssc_number_t *timehr = allocate("timehr", num_steps);
		ssc_number_t *month = allocate("month", num_steps);
		ssc_number_t *day = allocate("day", num_steps);
		ssc_number_t *hour = allocate("hour", num_steps);
		ssc_number_t *minute = allocate("minute", num_steps);

		double T = t_start;
		size_t idx = 0;
		while (T < t_end && idx < num_steps)
		{
			double Thr = T / 3600.0;

			time[idx] = (float) T;
			timehr[idx] = (float) Thr;
						
			int m = util::month_of(Thr);
			month[idx] = (ssc_number_t) m ;              // month goes 1-12
			day[idx] = (ssc_number_t) util::day_of_month(m,Thr) ;   // day goes 1-nday_in_month
			hour[idx] = (ssc_number_t) ((int)(Thr)%24);		         // hour goes 0-23
			minute[idx] = (ssc_number_t) ((int)( (Thr-floor(Thr))*60  + t_step/3600.0*30));      // minute goes 0-59
	
			T += t_step;
			idx++;
		}

	}
};

DEFINE_MODULE_ENTRY( timeseq, "Time sequence generator", 1 )
