//PVSAMV2
//B!!!!!!!!!!!!


//THIS VERSION STEPS THROUGH THE SYSTEM RATHER THAN TIME. HOWEVER, THE SYSTEM'S YEAR IS STILL WRAPPED IN A LIFETIME LOOP.
//NEED TO REMOVE THIS OR PVSAMV2 FROM SSCAPI.CPP BEFORE WE MOVE FORWARD!!!++++++++++++++++++++++++++


/* A row of asterisks *************** indicates the beginning of a new section of code.
   A row of plus signs ++++++++++ indicates things that need to be modified or filled in.*/

#include "core.h"
#include "lib_util.h"
#include "lib_soiling_b.h"
#include "lib_sunpos_b.h"
#include "lib_sandia_b.h"
#include "lib_pvinv_b.h"

//VARIABLE TABLE COULD BE MOVED TO BE WITH SECTIONS OF CODE++++++++++++
static var_info _cm_vtab_pvsamv2b[] = {
	/* VARTYPE           DATATYPE         NAME                                   LABEL                                      UNITS     META      GROUP        REQUIRED_IF           CONSTRAINTS               UI_HINTS*/
	//{ SSC_INPUT,        SSC_NUMBER,      "pv_lifetime_simulation",            "PV lifetime simulation",                     "0/1",    "",       "pvsamv2",   "+=0",               "INTEGER,MIN=0,MAX=1",          "" },
	//INPUTS																																												    
	{ SSC_INPUT,        SSC_ARRAY,       "soiling",                            "Monthly soiling loss",                       "%",      "",      "pvsamv2b",                         "*",             "LENGTH=12",     "" },         
	{ SSC_INPUT,        SSC_NUMBER,      "acwiring_loss",                      "AC wiring loss",                             "%",      "",      "pvsamv2b",                         "*",             "MIN=0,MAX=100", "" },
	{ SSC_INPUT,        SSC_NUMBER,      "transformer_loss",                   "AC step-up transformer loss",                "%",      "",      "pvsamv2b",                         "*",             "MIN=0,MAX=100", "" },


	// Inverter inputs - general
	{ SSC_INPUT, SSC_NUMBER, "inverter_model", "Inverter model specifier", "", "0=cec,1=datasheet,2=partload", "pvsamv2b", "*", "INTEGER,MIN=0,MAX=2", "" },
	{ SSC_INPUT, SSC_NUMBER, "mppt_low_inverter", "Minimum inverter MPPT voltage window", "Vdc", "", "pvsamv2b", "", "?=0", "" },
	{ SSC_INPUT, SSC_NUMBER, "mppt_hi_inverter", "Maximum inverter MPPT voltage window", "Vdc", "", "pvsamv2b", "", "?=0", "" },
	// this will need to be an array for different inverters
	{ SSC_INPUT, SSC_NUMBER, "inverter_count", "Number of inverters", "", "", "pvsamv2b", "*", "INTEGER,POSITIVE", "" },

	// Invverter inputs - cec database (note - same as datasheet below - we should consolidate with the CEC parameter compute module in the works)
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_c0", "Curvature between ac-power and dc-power at ref", "1/W", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_c1", "Coefficient of Pdco variation with dc input voltage", "1/V", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_c2", "Coefficient of Pso variation with dc input voltage", "1/V", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_c3", "Coefficient of Co variation with dc input voltage", "1/V", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_paco", "AC maximum power rating", "Wac", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_pdco", "DC input power at which ac-power rating is achieved", "Wdc", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_pnt", "AC power consumed by inverter at night", "Wac", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_pso", "DC power required to enable the inversion process", "Wdc", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_vdco", "DC input voltage for the rated ac-power rating", "Vdc", "", "pvsamv2b", "inverter_model=0", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_snl_vdcmax", "Maximum dc input operating voltage", "Vdc", "", "pvsamv2b", "inverter_model=0", "", "" },

	// Inverter inputs - datasheet model specific
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_paco", "AC maximum power rating", "Wac", "", "pvsamv2b", "inverter_model=1", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_eff", "Weighted or Peak or Nominal Efficiency", "Wdc", "", "pvsamv2b", "inverter_model=1", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_pnt", "AC power consumed by inverter at night", "Wac", "", "pvsamv2b", "inverter_model=1", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_pso", "DC power required to enable the inversion process", "Wdc", "", "pvsamv2b", "inverter_model=1", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_vdco", "DC input voltage for the rated ac-power rating", "Vdc", "", "pvsamv2b", "inverter_model=1", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_ds_vdcmax", "Maximum dc input operating voltage", "Vdc", "", "pvsamv2b", "inverter_model=1", "", "" },

	// Inverter inputs - part load model specific
	{ SSC_INPUT, SSC_NUMBER, "inv_pd_paco", "AC maximum power rating", "Wac", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_pd_pdco", "DC input power at which ac-power rating is achieved", "Wdc", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_ARRAY, "inv_pd_partload", "Partload curve partload values", "%", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_ARRAY, "inv_pd_efficiency", "Partload curve efficiency values", "%", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_pd_pnt", "AC power consumed by inverter at night", "Wac", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_pd_vdco", "DC input voltage for the rated ac-power rating", "Vdc", "", "pvsamv2b", "inverter_model=2", "", "" },
	{ SSC_INPUT, SSC_NUMBER, "inv_pd_vdcmax", "Maximum dc input operating voltage", "Vdc", "", "pvsamv2b", "inverter_model=2", "", "" },



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
	{ SSC_OUTPUT,       SSC_NUMBER,      "ac_loss",                            "Interconnection AC loss",                    "%",      "",      "Annual",                          "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_clip_loss_percent",    "AC inverter power clipping loss",            "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_pso_loss_percent",     "AC inverter power consumption loss",         "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_pnt_loss_percent",     "AC inverter night tare loss",                "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_inv_eff_loss_percent",     "AC inverter efficiency loss",                "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_wiring_loss_percent",      "AC wiring loss",                             "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_transformer_loss_percent", "AC step-up transformer loss",                "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_perf_adj_loss_percent",    "AC performance adjustment loss",             "%",      "",      "Loss",                            "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_wiring_loss",              "AC wiring loss",                            "kWh",     "",      "Annual",                          "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_ac_transformer_loss",         "AC step-up transformer loss",               "kWh",     "",      "Annual",                          "",             "",              "" },
	{ SSC_OUTPUT,       SSC_NUMBER,      "annual_dc_optimizer_loss",           "DC power optimizer loss",                   "kWh",     "",      "Annual",                          "",             "",              "" },


	{ SSC_OUTPUT, SSC_ARRAY, "dc_power", "(AC) Power output from inverter", "kW", "", "Time Series", "*", "", "" },
	{ SSC_OUTPUT, SSC_ARRAY, "ac_power", "(AC) Power output from inverter", "kW", "", "Time Series", "*", "", "" },

	// Lifetime output
	{ SSC_OUTPUT, SSC_ARRAY, "gen", "Power generated by system", "kW", "", "Time Series", "*", "", "" },

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


		//Inverter initialization
		::sandia_inverter_t_b snlinv;
		//::partload_inverter_t_b plinv;

		int num_inverters = as_integer("inverter_count");
		int inv_type = as_integer("inverter_model");

		if (inv_type == 0) // cec database
		{
			snlinv.Paco = as_double("inv_snl_paco");
			snlinv.Pdco = as_double("inv_snl_pdco");
			snlinv.Vdco = as_double("inv_snl_vdco");
			snlinv.Pso = as_double("inv_snl_pso");
			snlinv.Pntare = as_double("inv_snl_pnt");
			snlinv.C0 = as_double("inv_snl_c0");
			snlinv.C1 = as_double("inv_snl_c1");
			snlinv.C2 = as_double("inv_snl_c2");
			snlinv.C3 = as_double("inv_snl_c3");

		}
		else if (inv_type == 1) // datasheet data
		{
			double eff_ds = as_double("inv_ds_eff") / 100.0;
			snlinv.Paco = as_double("inv_ds_paco");
			if (eff_ds != 0)
				snlinv.Pdco = snlinv.Paco / eff_ds;
			else
				snlinv.Pdco = 0;
			snlinv.Vdco = as_double("inv_ds_vdco");
			snlinv.Pso = as_double("inv_ds_pso");
			snlinv.Pntare = as_double("inv_ds_pnt");
			snlinv.C0 = 0;
			snlinv.C1 = 0;
			snlinv.C2 = 0;
			snlinv.C3 = 0;
		}
		/*
		else if (inv_type == 2) // partload curve
		{
			plinv.Paco = as_double("inv_pd_paco");
			plinv.Pdco = as_double("inv_pd_pdco");
			plinv.Pntare = as_double("inv_pd_pnt");

			std::vector<double> pl_pd = as_doublevec("inv_pd_partload");
			std::vector<double> eff_pd = as_doublevec("inv_pd_efficiency");

			plinv.Partload = pl_pd;
			plinv.Efficiency = eff_pd;
		}
		*/
		else
		{
			throw exec_error("pvsamv2b", "invalid inverter model type");
		}

		//MODEL OVERRIDES*****************************************************************************************************************************************************************************
		bool user_soiling_model = false; //eventually read in from UI++++++++++++++++++++++++

		// ALLOCATE intermediates (not used for outputs**************************************************************************************************************************************************************************** 
		// discuss what data structure to use for intermediates+++++++++++++++++++++++

		// inverter model input - will come from upstream
		util::matrix_t<ssc_number_t> dcpwr_net(nrec);
		util::matrix_t<ssc_number_t> inv_vdc(nrec);
		// inverter model outputs
		util::matrix_t<ssc_number_t> acpwr_gross(nrec);
		util::matrix_t<ssc_number_t> ac_parasitics(nrec);
		util::matrix_t<ssc_number_t> pl_ratio(nrec);
		util::matrix_t<ssc_number_t> aceff(nrec);
		util::matrix_t<ssc_number_t> cliploss(nrec);
		util::matrix_t<ssc_number_t> psoloss(nrec);
		util::matrix_t<ssc_number_t> pntloss(nrec);


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


		ssc_number_t *p_dcpwr_net = allocate("dc_power", nrec);
		ssc_number_t *p_acpwr_gross = allocate("ac_power", nrec);



		// lifetime outputs
		ssc_number_t *p_gen = allocate("gen", nrec * nyears);


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
			//faking this for now++++++++++++
			for (size_t i = 0; i < (nrec); i++)
			{
				dcpwr_net.at(i) = snlinv.Pdco - 0.5*snlinv.Pdco * ((float)rand() / (float)RAND_MAX);
				inv_vdc.at(i) = snlinv.Vdco - 0.1*snlinv.Vdco * ((float)rand() / (float)RAND_MAX);
			}

			if ((inv_type == 0) || (inv_type == 1))
			{
				snlinv.acpower(dcpwr_net, inv_vdc, num_inverters,
					acpwr_gross, ac_parasitics, pl_ratio, aceff, cliploss, psoloss, pntloss);
			}
			/*
			else if (inv_type == 2)
			{
				plinv.acpower(dcpwr_net / num_inverters, &acpwr_gross, &_par, &_plr, &aceff, &cliploss, &pntloss);
			}
			*/
			
			//BATTERY MODEL

			//AC LOSSES
			apply_ac_loss(nrec, current_year, p_gen, acpwr_gross, ac_loss_percent);

			if (current_year == 0)
			{
				for (size_t ii = 0; ii < nrec; ii++)
				{
					p_dcpwr_net[ii] = dcpwr_net.at(ii);
					p_acpwr_gross[ii] = acpwr_gross.at(ii);
				}
			}
		}
	}
};

DEFINE_MODULE_ENTRY(pvsamv2b, "Photovoltaic Performance Model, Version 2b", 1)