//PVSAMV2
//B!!!!!!!!!!!!


//THIS VERSION STEPS THROUGH THE SYSTEM RATHER THAN TIME. HOWEVER, THE SYSTEM'S YEAR IS STILL WRAPPED IN A LIFETIME LOOP.
//NEED TO REMOVE THIS OR PVSAMV2 FROM SSCAPI.CPP BEFORE WE MOVE FORWARD!!!++++++++++++++++++++++++++


/* A row of asterisks *************** indicates the beginning of a new section of code.
   A row of plus signs ++++++++++ indicates things that need to be modified or filled in.*/

#include "core.h"
#include "lib_soiling_b.h"

//VARIABLE TABLE COULD BE MOVED TO BE WITH SECTIONS OF CODE++++++++++++
static var_info _cm_vtab_pvsamv2b[] = {
	/* VARTYPE           DATATYPE         NAME                                   LABEL                                                         UNITS     META                                                   GROUP                      REQUIRED_IF           CONSTRAINTS               UI_HINTS*/
	//{ SSC_INPUT,        SSC_NUMBER,      "pv_lifetime_simulation",              "PV lifetime simulation",                                      "0/1",    "",                                                   "pvsamv2",                  "+=0",               "INTEGER,MIN=0,MAX=1",          "" },
	//INPUTS
	{ SSC_INPUT,        SSC_ARRAY,       "soiling",                   "Monthly soiling loss",                      "%",       "",                              "pvsamv1",              "*",                        "LENGTH=12",                      "" },         
	
	//OUTPUTS
	{ SSC_OUTPUT,       SSC_ARRAY,       "month",                     "Month",                                      "",   "",                          "Time Series (Subarray 4)",       "",                    "",                              "" },
	{ SSC_OUTPUT,       SSC_ARRAY,       "pretend_poa_eff",           "POA total irradiance after soiling",         "W/m2",   "",                      "Time Series (Subarray 4)",       "",                    "",                              "" },

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

		//Soiling
		size_t soiling_length = 0;
		ssc_number_t *soiling = as_array("soiling", &soiling_length); //don't need to error check length of soiling array because it is constrained to 12 in variable table

		//MODEL OVERRIDES*****************************************************************************************************************************************************************************
		bool user_soiling_model = false; //eventually read in from UI++++++++++++++++++++++++

		//ALLOCATE OUTPUTS****************************************************************************************************************************************************************************
		//Again, if this becomes a meta-compute module, outputs are allocated individually by functions++++ Somehow we have to know what has been created by what point, though.
		ssc_number_t *month = allocate("month", nrec);
		std::vector<double> ibeam, idiff, ipoa;
		ibeam.reserve(nrec);
		idiff.reserve(nrec);
		ipoa.reserve(nrec);
		ssc_number_t *pretend_poa_eff = allocate("pretend_poa_eff", nrec);

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

		}
	}
};

DEFINE_MODULE_ENTRY(pvsamv2b, "Photovoltaic Performance Model, Version 2b", 1)