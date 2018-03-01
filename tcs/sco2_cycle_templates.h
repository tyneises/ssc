#ifndef __SCO2_CYCLE_TEMPLATES_
#define __SCO2_CYCLE_TEMPLATES_

#include "sco2_cycle_components.h"
#include "heat_exchangers.h"

class C_sco2_cycle_core
{
public:

	enum E_cycle_state_points
	{
		// index values for c++ 0-based vectors for temperature, pressure, etc.
		MC_IN = 0,		// Main compressor inlet
		MC_OUT,			// Main compressor outlet
		LTR_HP_OUT,		// Low temp recuperator high pressure outlet
		MIXER_OUT,		// Mixer: LTR_HP_OUT + Recompressor outlet
		HTR_HP_OUT,		// High temp recuperator high pressure outlet
		TURB_IN,		// Turbine inlet
		TURB_OUT,		// Turbine outlet
		HTR_LP_OUT,		// High temp recuperator low pressure outlet
		LTR_LP_OUT,		// Low temp recuperator low pressure outlet
		RC_OUT,			// Recompresor outlet
		PC_IN,			// Precompressor inlet (partial cooling cycle)
		PC_OUT,			// Precompressor outlet (partial cooling cycle)

		END_SCO2_STATES
	};

	struct S_design_solved
	{
		std::vector<double> m_temp, m_pres, m_enth, m_entr, m_dens;		// thermodynamic states (K, kPa, kJ/kg, kJ/kg-K, kg/m3)
		double m_eta_thermal;	//[-]
		double m_W_dot_net;		//[kWe]
		double m_m_dot_mc;		//[kg/s]
		double m_m_dot_rc;		//[kg/s]
		double m_m_dot_pc;		//[kg/s]
		double m_m_dot_t;		//[kg/s]
		double m_recomp_frac;	//[-]
		double m_UA_LTR;			//[kW/K]
		double m_UA_HTR;			//[kW/K]

		bool m_is_rc;

		C_comp_multi_stage::S_des_solved ms_mc_ms_des_solved;
		C_comp_multi_stage::S_des_solved ms_rc_ms_des_solved;
		C_comp_multi_stage::S_des_solved ms_pc_ms_des_solved;
		C_turbine::S_design_solved ms_t_des_solved;
		C_HX_counterflow::S_des_solved ms_LTR_des_solved;
		C_HX_counterflow::S_des_solved ms_HTR_des_solved;

		S_design_solved()
		{
			m_eta_thermal = m_W_dot_net = m_m_dot_mc = m_m_dot_rc = m_m_dot_t = m_recomp_frac =
				m_UA_LTR = m_UA_HTR = std::numeric_limits<double>::quiet_NaN();

			m_is_rc = true;
		}
	};

	struct S_auto_opt_design_parameters
	{
		double m_W_dot_net;					//[kWe] Target net cycle power
		double m_T_mc_in;					//[K] Main compressor inlet temperature
		double m_T_pc_in;					//[K] Pre-compressor inlet temperature
		double m_T_t_in;					//[K] Turbine inlet temperature
		std::vector<double> m_DP_LTR;		//(cold, hot) positive values are absolute [kPa], negative values are relative (-)
		std::vector<double> m_DP_HTR;		//(cold, hot) positive values are absolute [kPa], negative values are relative (-)
		std::vector<double> m_DP_PC_pre;    //(cold, hot) positive values are absolute [kPa], negative values are relative (-)
		std::vector<double> m_DP_PC_main;   //(cold, hot) positive values are absolute [kPa], negative values are relative (-)
		std::vector<double> m_DP_PHX;		//(cold, hot) positive values are absolute [kPa], negative values are relative (-)
		double m_UA_rec_total;				//[kW/K] Total design-point recuperator UA
		double m_LTR_eff_max;				//[-] Maximum allowable effectiveness in LT recuperator
		double m_HTR_eff_max;				//[-] Maximum allowable effectiveness in HT recuperator
		double m_eta_mc;					//[-] design-point efficiency of the main compressor; isentropic if positive, polytropic if negative
		double m_eta_rc;					//[-] design-point efficiency of the recompressor; isentropic if positive, polytropic if negative
		double m_eta_pc;					//[-] design-point efficiency of the pre-compressor; 
		double m_eta_t;						//[-] design-point efficiency of the turbine; isentropic if positive, polytropic if negative
		int m_N_sub_hxrs;					//[-] Number of sub-heat exchangers to use when calculating UA value for a heat exchanger
		double m_P_high_limit;				//[kPa] maximum allowable pressure in cycle
		double m_tol;						//[-] Convergence tolerance
		double m_opt_tol;					//[-] Optimization tolerance
		double m_N_turbine;					//[rpm] Turbine shaft speed (negative values link turbine to compressor)
		
		int m_is_recomp_ok;					//[-] 1 = yes, 0 = no, other = invalid

		double m_PR_mc_guess;				//[-] Initial guess for ratio of P_mc_out to P_mc_in
		bool m_fixed_PR_mc;					//[-] if true, ratio of P_mc_out to P_mc_in is fixed at PR_mc_guess
		
		int m_des_objective_type;		//[2] = min phx deltat then max eta, [else] max eta
		double m_min_phx_deltaT;		//[C]

		// Callback function only log
		bool(*mf_callback_log)(std::string &log_msg, std::string &progress_msg, void *data, double progress, int out_type);
		void *mp_mf_active;


		S_auto_opt_design_parameters()
		{
			m_W_dot_net = m_T_mc_in = m_T_pc_in = m_T_t_in =
				m_UA_rec_total = m_LTR_eff_max = m_HTR_eff_max =
				m_eta_mc = m_eta_rc = m_eta_pc = m_eta_t = m_P_high_limit = m_tol = m_N_turbine = 
				m_fixed_PR_mc = std::numeric_limits<double>::quiet_NaN();
			m_N_sub_hxrs = -1;

			m_is_recomp_ok = 1;

			m_fixed_PR_mc = false;		//[-] If false, then should default to optimizing this parameter
			
			// Default to standard optimization to maximize cycle efficiency
			m_des_objective_type = 1;
			m_min_phx_deltaT = 0.0;		//[C]

			mf_callback_log = 0;
			mp_mf_active = 0;

			m_DP_LTR.resize(2);
			std::fill(m_DP_LTR.begin(), m_DP_LTR.end(), std::numeric_limits<double>::quiet_NaN());
			m_DP_HTR.resize(2);
			std::fill(m_DP_HTR.begin(), m_DP_HTR.end(), std::numeric_limits<double>::quiet_NaN());
			m_DP_PC_pre.resize(2);
			std::fill(m_DP_PC_pre.begin(), m_DP_PC_pre.end(), std::numeric_limits<double>::quiet_NaN());
			m_DP_PC_main.resize(2);
			std::fill(m_DP_PC_main.begin(), m_DP_PC_main.end(), std::numeric_limits<double>::quiet_NaN());
			m_DP_PHX.resize(2);
			std::fill(m_DP_PHX.begin(), m_DP_PHX.end(), std::numeric_limits<double>::quiet_NaN());
		}
	};
	
protected:

	S_design_solved ms_des_solved;

	S_auto_opt_design_parameters ms_auto_opt_des_par;

public:

	const S_design_solved * get_design_solved()
	{
		return &ms_des_solved;
	}

	virtual int auto_opt_design(S_auto_opt_design_parameters & auto_opt_des_par_in) = 0;
	
};


#endif // !__SCO2_CYCLE_TEMPLATES_
