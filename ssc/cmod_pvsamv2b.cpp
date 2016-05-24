//PVSAMV2
//B!!!!!!!!!!!!


//THIS VERSION STEPS THROUGH THE SYSTEM RATHER THAN TIME. HOWEVER, THE SYSTEM'S YEAR IS STILL WRAPPED IN A LIFETIME LOOP.
//NEED TO REMOVE THIS OR PVSAMV2 FROM SSCAPI.CPP BEFORE WE MOVE FORWARD!!!++++++++++++++++++++++++++


/* A row of asterisks *************** indicates the beginning of a new section of code.
   A row of plus signs ++++++++++ indicates things that need to be modified or filled in.*/

#include "core.h"
#include "lib_soiling_b.h"
#include "lib_sunpos_b.h"

//VARIABLE TABLE COULD BE MOVED TO BE WITH SECTIONS OF CODE++++++++++++
static var_info _cm_vtab_pvsamv2b[] = {
	/* VARTYPE           DATATYPE         NAME                                   LABEL                                      UNITS     META      GROUP        REQUIRED_IF           CONSTRAINTS               UI_HINTS*/
	//{ SSC_INPUT,        SSC_NUMBER,      "pv_lifetime_simulation",            "PV lifetime simulation",                     "0/1",    "",       "pvsamv2",   "+=0",               "INTEGER,MIN=0,MAX=1",          "" },
	//INPUTS																																												    
	{ SSC_INPUT,        SSC_ARRAY,       "soiling",                            "Monthly soiling loss",                       "%",      "",      "pvsamv1",                         "*",             "LENGTH=12",     "" },         
	{ SSC_INPUT,        SSC_NUMBER,      "acwiring_loss",                      "AC wiring loss",                             "%",      "",      "pvsamv1",                         "*",             "MIN=0,MAX=100", "" },
	{ SSC_INPUT,        SSC_NUMBER,      "transformer_loss",                   "AC step-up transformer loss",                "%",      "",      "pvsamv1",                         "*",             "MIN=0,MAX=100", "" },

	//OUTPUTS														          														   										    
	{ SSC_OUTPUT,       SSC_ARRAY,       "year",                               "Year",                                       "",       "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "month",                              "Month",                                      "",       "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "day",                                "Day",                                        "",       "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "hour",                               "Hour",                                       "",       "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "minute",                             "Minute",                                     "",       "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "solazi",                             "Solar azimuth",                              "deg",    "",      "Time Series",                     "",              "",              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "solzen",                             "Solar zenith",                               "deg",    "",      "Time Series",                     "",              "",              "" },
																																							                     
	{ SSC_OUTPUT,       SSC_ARRAY,       "pretend_poa_eff",                    "POA total irradiance after soiling",         "W/m2",   "",      "Time Series",                     "",              "",              "" },
																																	   
	// AC LOSSES								
	{ SSC_OUTPUT,       SSC_NUMBER,      "ac_loss",                            "Interconnection AC loss",                    "%",      "",      "Annual",                          "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_clip_loss_percent",    "AC inverter power clipping loss",            "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_pso_loss_percent",     "AC inverter power consumption loss",         "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_pnt_loss_percent",     "AC inverter night tare loss",                "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_eff_loss_percent",     "AC inverter efficiency loss",                "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_wiring_loss_percent",      "AC wiring loss",                             "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_transformer_loss_percent", "AC step-up transformer loss",                "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_perf_adj_loss_percent",    "AC performance adjustment loss",             "%",      "",      "Loss",                            "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_wiring_loss",              "AC wiring loss",                            "kWh",     "",      "Annual",                          "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_transformer_loss",         "AC step-up transformer loss",               "kWh",     "",      "Annual",                          "*",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_dc_optimizer_loss",           "DC power optimizer loss",                   "kWh",     "",      "Annual",                          "*",             "",              "" },

	var_info_invalid };

class cm_pvsamv2b : public compute_module
{
public:
	cm_pvsamv2b()
	{
		add_var_info(_cm_vtab_pvsamv2b);
	}


	//MAIN FUNCTION***********************************************************************************************************************************************************************************
	//************************************************************************************************************************************************************************************************
	void exec() throw(general_error)
	{
		//READ INPUTS*********************************************************************************************************************************************************************************
		//Inputs are read in from the variable table here. Error checking is performed here as applicable. Like inputs are grouped together.
		//IF THIS BECOMES A META-COMPUTE MODULE, it's possible that inputs are read in by the functions individually+++++++++++++++++++++

		//Assume you're given these variables++++++++++++++++++++++++
		int nyears = 1; //this would probably be 1 for non-lifetime++++++++++++++++++++++++
		int nrec = 8760; //how do we know what nrec is up here? we don't actually know until we read the weather file++++++++++++++++++

		// location info would come from weatherfile
		double lat = 40;
		double lng = -107;
		double tz = -1;

		//Soiling
		size_t soiling_length = 0;
		ssc_number_t *soiling = as_array("soiling", &soiling_length); //don't need to error check length of soiling array because it is constrained to 12 in variable table

		// AC Losses
		double ac_derate = (1 - as_double("acwiring_loss") / 100) * (1 - as_double("transformer_loss") / 100);	
		double ac_loss_percent = (1 - ac_derate) * 100;
		assign("ac_loss", var_data((ssc_number_t)ac_loss_percent));

		//MODEL OVERRIDES*****************************************************************************************************************************************************************************
		bool user_soiling_model = false; //eventually read in from UI++++++++++++++++++++++++

		//ALLOCATE OUTPUTS****************************************************************************************************************************************************************************
		//Again, if this becomes a meta-compute module, outputs are allocated individually by functions++++ Somehow we have to know what has been created by what point, though.
		ssc_number_t *year = allocate_fill( "year", nrec, 1990 );
		ssc_number_t *month = allocate("month", nrec );
		ssc_number_t *day = allocate_fill( "day", nrec, 17 );
		ssc_number_t *hour = allocate_fill( "hour", nrec, 11 );
		ssc_number_t *minute = allocate_fill( "minute", nrec, 15 );

		std::vector<double> ibeam, idiff, ipoa;
		ibeam.reserve(nrec);
		idiff.reserve(nrec);
		ipoa.reserve(nrec);
		
		ssc_number_t *p_solazi = allocate( "solazi", nrec );
		ssc_number_t *p_solzen = allocate( "solzen", nrec );

		ssc_number_t *pretend_poa_eff = allocate("pretend_poa_eff", nrec);
		ssc_number_t *p_gen = allocate("p_gen", nrec);

		double acpwr_gross = 800.;
		for (int i = 0; i != nrec; i++)
			p_gen[i] = 100;

		//LIFETIME LOOP*******************************************************************************************************************************************************************************
		//This loop goes through the number of years in the simulation. 
		//The loop is broken up by sections of code. Inputs for each section of code are assigned prior to the lifetime loop, but outputs are assigned as they occur++++++++++++
		for (int current_year = 0; current_year < nyears; current_year++)
		{
			//YEAR 1 ONLY: WEATHER FILE, IRRADIANCE PROCESSING, SHADING, ETC... EVERYTHING PRIOR TO THE MODULE ITSELF.********************************************************************************
			if (current_year == 0)
			{
				//READ WEATHER FILE
				//faking this for now++++++++++++
				for (int i = 0; i < 8760; i++)
				{
					ibeam.push_back(500);
					idiff.push_back(200);
					//faking this for now++++++++++++++++++++
					month[i] = floor(i / 720) + 1;
					if (month[i] > 12) month[i] = 12;
				}

				//SUBARRAY LOOP WOULD GO HERE++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
				//CALCULATE SUN POSITION
				calculate_sun_positions( nrec, year, month, day, hour, minute, 
					lat, lng, tz,
					p_solazi, p_solzen );



				//CALCULATE TRACKER POSITION AND SURFACE ANGLES

				//CALCULATE POA
				//faking this for now++++++++++++
				for (int i = 0; i < 8760; i++)
					ipoa.push_back(600);

				//APPLY SHADING (3D, SHADING FACTORS, SHADING DATABASE, SELF-SHADING)

				//APPLY SNOW LOSSES

				//APPLY SOILING
				if (user_soiling_model)
					ibeam[0] = 400; //a real link to a user model would go here++++++++++++++++++++++++
				else
				{
					//soiling losses are applied to all components of irradiance for the entire year
					apply_soiling_loss_b(nrec, month, soiling, &ibeam);
					apply_soiling_loss_b(nrec, month, soiling, &idiff);
					apply_soiling_loss_b(nrec, month, soiling, &ipoa);
				}
				//assign pretend ipoa output for irradiance after soiling
				//there's probably a better way to store and assign these arrays...++++++++++++++++++++++++++++++++++++++++++++++++++++
				for (int i = 0; i < 8760; i++)
					pretend_poa_eff[i] = (ssc_number_t) ipoa[i];
			}

			//ALL YEARS (YEAR 1 AND FORWARD): MODULE TEMP, POWER, MPPT, INVERTER, BATTERY, LOSSES, ETC... EVERYTHING FROM THE MODULE AND DOWNSTREAM.**************************************************

			//MODULE TEMPERATURE +++++++++++potentially iterative with module power

			//MODULE POWER ++++++++++potentially iterative with MPPT

			//STRING POWER/MPPT/MISMATCH?++++++++++

			//DC LOSSES

			//INVERTER MODEL
			
			//BATTERY MODEL

			//AC LOSSES
			apply_ac_loss(current_year, p_gen, acpwr_gross, ac_loss_percent);

		}
	}
};

DEFINE_MODULE_ENTRY(pvsamv2b, "Photovoltaic Performance Model, Version 2b", 1)