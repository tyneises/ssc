//PVSAMV2!!!!!!!!!!!!
/* A row of asterisks *************** indicates the beginning of a new section of code.
   A row of question marks ?????????? indicates things that need to be modified or filled in.*/

#include "core.h"
#include "lib_soiling.h"

static var_info _cm_vtab_pvsamv2[] = {
	/* VARTYPE           DATATYPE         NAME                                   LABEL                                                         UNITS     META                                                   GROUP                      REQUIRED_IF           CONSTRAINTS               UI_HINTS*/
	//{ SSC_INPUT,        SSC_NUMBER,      "pv_lifetime_simulation",              "PV lifetime simulation",                                      "0/1",    "",                                                   "pvsamv2",                  "?=0",               "INTEGER,MIN=0,MAX=1",          "" },
	{ SSC_INPUT,        SSC_ARRAY,       "soiling",                   "Monthly soiling loss",                      "%",       "",                              "pvsamv1",              "*",                        "LENGTH=12",                      "" },         

	var_info_invalid };

class cm_pvsamv2 : public compute_module
{
public:
	cm_pvsamv2()
	{
		add_var_info(_cm_vtab_pvsamv2);
	}


	//MAIN FUNCTION***********************************************************************************************************************************************************************************
	//************************************************************************************************************************************************************************************************
	void exec() throw(general_error)
	{
		//READ INPUTS*********************************************************************************************************************************************************************************
		//Inputs are read in from the variable table here. Error checking is performed here as applicable. Like inputs are grouped together.
		
		//Assume you're given these variables????????????????????????
		int nyears = 25;
		int ntimesteps = 1;
		int nsubarrays = 4;
		int radmode = 0;

		//Soiling
		size_t soiling_length = 0;
		ssc_number_t *soiling = as_array("soiling", &soiling_length); //don't need to error check length of soiling array because it is constrained to 12 in variable table


		//LIFETIME LOOP*******************************************************************************************************************************************************************************
		//This loop goes through the number of years in the simulation. If 
		int idx = 0; //keep track of where you are in the simulation
		for (int year_counter = 0; year_counter < nyears; year_counter++)
		{
			//Hourly loop*****************************************************************************************************************************************************************************
			for (int hour_counter = 0; hour_counter < 8760; hour_counter++)
			{
				//Sub-hourly loop*********************************************************************************************************************************************************************
				for (int step = 0; step < ntimesteps; step++)
				{
					//READ WEATHER FILE***************************************************************************************************************************************************************
					//time values would all be assigned from weather file values????????????????????????
					int year = 2003;
					int month = 3;
					int day = 14;
					int hour = 15;
					int min = 0;


					//Subarray loop*******************************************************************************************************************************************************************
					for (int subarray = 0; subarray < nsubarrays; subarray++)
					{


						//Soiling*********************************************************************************************************************************************************************
						//Soiling is the same for all four subarrays. It is passed in as an array of 12 values (monthly values).
						//assume you already have irrad components and know the irrad mode?????????????????????????
						double ibeam = 500; 
						double idiff = 200;
						double ipoa = 430;
						int radmode = 0;

						double soiling_loss = soiling[month];
						//if (radmode > 2 && soiling_loss)
							
						
						apply_soiling_loss(&ibeam, soiling_loss);
						apply_soiling_loss(&idiff, soiling_loss);
						apply_soiling_loss(&ipoa, soiling_loss);


					}
				}
			}
		}

		 



		
	}
};

DEFINE_MODULE_ENTRY(pvsamv2, "Photovoltaic Performance Model, Version 2", 1)