﻿/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  (“Alliance”) under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
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
*  the underlying software originally provided by Alliance as “System Advisor Model” or “SAM”. Except
*  to comply with the foregoing, the terms “System Advisor Model”, “SAM”, or any confusingly similar
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
#include "lib_weatherfile.h"
#include "RegenHX.h"
// for adjustment factors
#include "common.h"


static var_info _cm_vtab_regenhx[] = {
	//	  VARTYPE          DATATYPE          NAME                                        LABEL                                            UNITS            META     GROUP                REQUIRED_IF        CONSTRAINTS           UI_HINTS
		{ SSC_INPUT,       SSC_STRING,       "file_name",                                "Local weather file path",                        "",              "",      "dmitrii",          "",               "LOCAL_FILE",          "" },



		{ SSC_INPUT,       SSC_NUMBER,       "inletStates.T_H_in",					"",                    "K",					"",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "inletStates.P_H",						"",                    "kPa",				"",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "inletStates.T_C_in",					"",                    "K",					"",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "inletStates.P_C",						"",                    "kPa",				"",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.m_dot_H",					"",                    "kg/s",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.m_dot_C",					"",                    "kg/s",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.Q_dot_loss",				"",                    "kW",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.P_0",						"",                    "s",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.D_s",						"",                    "m",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "parameters.e_v",						"",                    "-",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.target_2_value",				"",                    "",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.target_1",					"",                    "-",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.target_2",					"",                    "-",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.target_1_value",				"",                    "",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.operationMode",				"",                    "-",              "",      "",          "*",               "",                    "" },
		{ SSC_INPUT,       SSC_NUMBER,       "targets.valveMode",					"",                    "-",              "",      "",          "*",               "",                    "" },



		//    OUTPUTS ----------------------------------------------------------------------------								      														   
		//	  VARTYPE           DATATYPE         NAME                          LABEL                                                   UNITS            META     GROUP                REQUIRED_IF        CONSTRAINTS           UI_HINTS

			{ SSC_OUTPUT,       SSC_NUMBER,      "epsilon_c",             "",                         "-",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "D_fr_c",				"",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "L_c",					"",                         "m",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "wallThickness_c",       "",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "UA_c",					"",                         "kW/K",        "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "totalCost_c",           "",                         "$",           "",      "",           "",               "",                    "" },

			{ SSC_OUTPUT,       SSC_NUMBER,      "epsilon_e",             "",                         "-",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "D_fr_e",				"",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "L_e",					"",                         "m",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "wallThickness_e",       "",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "UA_e",					"",                         "kW/K",        "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "totalCost_e",           "",                         "$",           "",      "",           "",               "",                    "" },

			{ SSC_OUTPUT,       SSC_NUMBER,      "epsilon_u",             "",                         "-",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "D_fr_u",				"",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "L_u",					"",                         "m",           "",      "",           "*",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "wallThickness_u",       "",                         "m",           "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "UA_u",					"",                         "kW/K",        "",      "",           "",               "",                    "" },
			{ SSC_OUTPUT,       SSC_NUMBER,      "totalCost_u",           "",                         "$",           "",      "",           "",               "",                    "" },


		var_info_invalid };

class cm_regenhx : public compute_module
{
private:
public:
	cm_regenhx()
	{
		add_var_info(_cm_vtab_regenhx);
		// performance adjustment factors
		//add_var_info(vtab_adjustment_factors);
		//add_var_info(vtab_technology_outputs);
	}

	void exec() throw(general_error)
	{
		double T_H_in = as_double("inletStates.T_H_in");
		double P_H = as_double("inletStates.P_H");
		double T_C_in = as_double("inletStates.T_C_in");
		double P_C = as_double("inletStates.P_C");

		double m_dot_H = as_double("parameters.m_dot_H");
		double m_dot_C = as_double("parameters.m_dot_C");
		double Q_dot_loss = as_double("parameters.Q_dot_loss");
		double P_0 = as_double("parameters.P_0");
		double D_s = as_double("parameters.D_s");
		double e_v = as_double("parameters.e_v");

		int target_1 = as_integer("targets.target_1");
		int target_2 = as_integer("targets.target_2");
		int operationMode = as_integer("targets.operationMode");
		int valveMode = as_integer("targets.valveMode");

		double target_2_value = as_double("targets.target_2_value");
		double target_1_value = as_double("targets.target_1_value");

		ofstream uaFile;
		string SSCDIR1(std::getenv("SSCDIR"));
		uaFile.open(SSCDIR1 + "/build_sdk/examples/Regen-UA.log");
		//uaFile << "Epsilon,\tUA,\t\tNTU,\tD_fr,\tL,\t\tL/D,\t\tT_H_out,\tComass,\t\tCost,\t\tdP_max,\tms" << endl;
		//uaFile << "Epsilon,UA,Q_dot,NTU,D_fr,L,L/D,T_H_out,Comass,Cost,dP_max,ms" << endl;
		char* format = "%.5f,\t%.0f,\t%.2f,\t%.2f,\t%.2f,\t%.2f,\t%.2f,\t\t%.2f,\t\t%.2f,\t\t%.0f,\t\t%.2f,\t\t\t\t%.0f,\t%.0f,\t%.0f,\t%.0f,\t%.0f";
		//char* format = "%.5f,%.0f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.0f,%.2f,%.0f";

		RegenHX* HT_regen = new RegenHX();

		double* results = new double[9];
		char* output = new char[150];
		clock_t begin, end;


		double q_dot, T_c_out, T_h_out, comass;
		//for (double cost = 650000; cost <= 750000; cost += 100) {
			comass = q_dot = T_c_out = T_h_out = std::numeric_limits<double>::quiet_NaN();

			begin = clock();
			try {
				HT_regen->set_params(target_1, target_2, operationMode, valveMode, target_2_value, P_0, D_s, e_v, Q_dot_loss);
				HT_regen->design_fix_TARGET_calc_outlet(target_1, target_1_value, 0.99, T_C_in, P_C, m_dot_C, P_C, T_H_in, P_H, m_dot_H, P_H - target_2_value, q_dot, T_c_out, T_h_out);
				end = clock();
				sprintf(output, format,
					HT_regen->ms_des_solved.m_eff_design,
					HT_regen->ms_des_solved.m_aUA_design_total,
					HT_regen->ms_des_solved.m_Q_dot_design,
					HT_regen->ms_des_solved.m_NTU_design,
					HT_regen->getD_fr(),
					HT_regen->getL(),
					HT_regen->getAspectRatio(),
					T_h_out,
					HT_regen->ms_des_solved.m_m_dot_carryover,
					HT_regen->getCost(),
					HT_regen->ms_des_solved.m_DP_hot_des,
					HT_regen->ms_des_solved.m_HTR_valve_HTHP_cv,
					HT_regen->ms_des_solved.m_HTR_valve_LTHP_cv,
					HT_regen->ms_des_solved.m_HTR_valve_HTLP_cv,
					HT_regen->ms_des_solved.m_HTR_valve_LTLP_cv,
					double(end - begin) / CLOCKS_PER_SEC * 1000.0);
				uaFile << output << endl;
				uaFile.flush();
			}
			catch(C_csp_exception &){
				end = clock();
				uaFile << "-,\t\t\t" << target_1_value << ",\t-,\t\t-,\t\t-,\t\t-,\t\t\t-,\t\t\t-,\t\t\t-,\t\t\t\t" << double(end - begin) / CLOCKS_PER_SEC * 1000.0 << endl;
				//continue;
			}
		//}

		uaFile.flush();
		uaFile.close();
	}

};

DEFINE_MODULE_ENTRY(regenhx, "Let's find out where this text goes", 2);

