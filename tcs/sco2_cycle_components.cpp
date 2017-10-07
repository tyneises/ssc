#include "sco2_cycle_components.h"
#include "CO2_properties.h"
#include <limits>
#include <algorithm>

#include "numeric_solvers.h"
#include "csp_solver_core.h"

const double C_turbine::m_nu_design = 0.7476;
const double C_comp_single_stage::m_snl_phi_design = 0.02971;		//[-] Design-point flow coef. for Sandia compressor (corresponds to max eta)
const double C_comp_single_stage::m_snl_phi_min = 0.02;				//[-] Approximate surge limit for SNL compressor
const double C_comp_single_stage::m_snl_phi_max = 0.05;				//[-] Approximate x-intercept for SNL compressor


void calculate_turbomachinery_outlet_1(double T_in /*K*/, double P_in /*kPa*/, double P_out /*kPa*/, double eta /*-*/, bool is_comp, int & error_code, double & spec_work /*kJ/kg*/)
{
	double enth_in, entr_in, dens_in, temp_out, enth_out, entr_out, dens_out;

	calculate_turbomachinery_outlet_1(T_in, P_in, P_out, eta, is_comp, error_code, enth_in, entr_in, dens_in, temp_out, enth_out, entr_out, dens_out, spec_work);
}

void calculate_turbomachinery_outlet_1(double T_in /*K*/, double P_in /*kPa*/, double P_out /*kPa*/, double eta /*-*/, bool is_comp, int & error_code, double & enth_in /*kJ/kg*/, double & entr_in /*kJ/kg-K*/,
	double & dens_in /*kg/m3*/, double & temp_out /*K*/, double & enth_out /*kJ/kg*/, double & entr_out /*kJ/kg-K*/, double & dens_out /*kg/m3*/, double & spec_work /*kJ/kg*/)
{
	/*Calculates the outlet state of a compressor or turbine using its isentropic efficiency.
	is_comp = .true.means the turbomachine is a compressor(w = w_s / eta)
	is_comp = .false.means the turbomachine is a turbine(w = w_s * eta) */

	CO2_state co2_props;

	error_code = 0;

	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);		// properties at the inlet conditions
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_in = co2_props.enth;
	double s_in = co2_props.entr;
	dens_in = co2_props.dens;

	prop_error_code = CO2_PS(P_out, s_in, &co2_props);			// outlet enthalpy if compression/expansion is isentropic
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_s_out = co2_props.enth;

	double w_s = h_in - h_s_out;			// specific work if process is isentropic (negative for compression, positive for expansion)

	double w = 0.0;
	if (is_comp)
		w = w_s / eta;						// actual specific work of compressor (negative)
	else
		w = w_s * eta;						// actual specific work of turbine (positive)

	double h_out = h_in - w;

	prop_error_code = CO2_PH(P_out, h_out, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}

	enth_in = h_in;
	entr_in = s_in;
	temp_out = co2_props.temp;
	enth_out = h_out;
	entr_out = co2_props.entr;
	dens_out = co2_props.dens;
	spec_work = w;

	return;
};

void isen_eta_from_poly_eta(double T_in /*K*/, double P_in /*kPa*/, double P_out /*kPa*/, double poly_eta /*-*/, bool is_comp, int & error_code, double & isen_eta)
{
	/* 9.3.14: code written by John Dyreby, translated to C++ by Ty Neises
	! Calculate the isentropic efficiency that corresponds to a given polytropic efficiency
	! for the expansion or compression from T_in and P_in to P_out.
	!
	! Inputs:
	!   T_in -- inlet temperature (K)
	!   P_in -- inlet pressure (kPa)
	!   P_out -- outlet pressure (kPa)
	!   poly_eta -- polytropic efficiency (-)
	!   is_comp -- if .true., model a compressor (w = w_s / eta); if .false., model a turbine (w = w_s * eta)
	!
	! Outputs:
	!   error_trace -- an ErrorTrace object
	!   isen_eta -- the equivalent isentropic efficiency (-)
	!
	! Notes:
	!   1) Integration of small DP is approximated numerically by using 200 stages.
	!   2) No error checking is performed on the inlet and outlet pressures; valid pressure ratios are assumed. */

	CO2_state co2_props;

	// Properties at the inlet conditions
	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_in = co2_props.enth;
	double s_in = co2_props.entr;

	// Outlet enthalpy if compression/expansion is isentropic
	prop_error_code = CO2_PS(P_out, s_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_s_out = co2_props.enth;

	double stage_P_in = P_in;		// Initialize first stage inlet pressure
	double stage_h_in = h_in;		// Initialize first stage inlet enthalpy
	double stage_s_in = s_in;		// Initialize first stage inlet entropy

	int N_stages = 200;

	double stage_DP = (P_out - P_in) / (double)N_stages;

	double stage_P_out = -999.9;
	double stage_h_out = -999.9;

	for (int i = 1; i <= N_stages; i++)
	{
		stage_P_out = stage_P_in + stage_DP;

		// Outlet enthalpy if compression/expansion is isentropic
		prop_error_code = CO2_PS(stage_P_out, stage_s_in, &co2_props);
		if (prop_error_code != 0)
		{
			error_code = prop_error_code;
			return;
		}
		double stage_h_s_out = co2_props.enth;

		double w_s = stage_h_in - stage_h_s_out;		// specific work if process is isentropic
		double w = std::numeric_limits<double>::quiet_NaN();
		if (is_comp)
			w = w_s / poly_eta;
		else
			w = w_s * poly_eta;
		stage_h_out = stage_h_in - w;

		// Reset next stage inlet values
		stage_P_in = stage_P_out;
		stage_h_in = stage_h_out;

		prop_error_code = CO2_PH(stage_P_in, stage_h_in, &co2_props);
		if (prop_error_code != 0)
		{
			error_code = prop_error_code;
			return;
		}
		stage_s_in = co2_props.entr;
	}

	// Note: last stage outlet enthalpy is equivalent to turbomachinery outlet enthalpy
	if (is_comp)
		isen_eta = (h_s_out - h_in) / (stage_h_out - h_in);
	else
		isen_eta = (stage_h_out - h_in) / (h_s_out - h_in);
}


void C_HeatExchanger::initialize(const S_design_parameters & des_par_in)
{
	ms_des_par = des_par_in;
	return;
}

void C_HeatExchanger::hxr_pressure_drops(const std::vector<double> & m_dots, std::vector<double> & hxr_deltaP)
{
	int N = (int)m_dots.size();
	hxr_deltaP.resize(N);
	for (int i = 0; i < N; i++)
		hxr_deltaP[i] = ms_des_par.m_DP_design[i] * pow((m_dots[i] / ms_des_par.m_m_dot_design[i]), 1.75);
}

void C_HeatExchanger::hxr_conductance(const std::vector<double> & m_dots, double & hxr_UA)
{
	double m_dot_ratio = 0.5*(m_dots[0] / ms_des_par.m_m_dot_design[0] + m_dots[1] / ms_des_par.m_m_dot_design[1]);
	hxr_UA = ms_des_par.m_UA_design*pow(m_dot_ratio, 0.8);
}

void C_turbine::turbine_sizing(const S_design_parameters & des_par_in, int & error_code)
{
	/* 9.4.14: code from John Dyreby, converted to C++ by Ty Neises
	! Determine the turbine rotor diameter, effective nozzle area, and design-point shaft
	! speed and store values in recomp_cycle%t.
	!
	! Arguments:
	!   recomp_cycle -- a RecompCycle object that defines the simple/recompression cycle at the design point
	!   error_trace -- an ErrorTrace object
	!
	! Notes:
	!   1) The value for recomp_cycle%t%N_design is required to be set.  If it is <= 0.0 then
	!      the value for recomp_cycle%mc%N_design is used (i.e., link the compressor and turbine
	!      shafts).  For this reason, turbine_sizing must be called after compressor_sizing if
	!      the shafts are to be linked. */

	CO2_state co2_props;

	ms_des_par = des_par_in;

	// Check that a design-point shaft speed is available
	if (ms_des_par.m_N_design <= 0.0)	// Link shafts
	{
		ms_des_solved.m_N_design = ms_des_par.m_N_comp_design_if_linked;
		if (ms_des_par.m_N_design <= 0.0)
		{
			error_code = 7;
			return;
		}
	}
	else
		ms_des_solved.m_N_design = ms_des_par.m_N_design;

	// Get speed of sound at inlet
	int prop_error_code = CO2_TD(ms_des_par.m_T_in, ms_des_par.m_D_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double ssnd_in = co2_props.ssnd;

	// Outlet specific enthalpy after isentropic expansion
	prop_error_code = CO2_PS(ms_des_par.m_P_out, ms_des_par.m_s_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_s_out = co2_props.enth;

	// Determine necessary turbine parameters
	ms_des_solved.m_nu_design = m_nu_design;
	double w_i = ms_des_par.m_h_in - h_s_out;			//[kJ/kg] Isentropic specific work of turbine
	double C_s = sqrt(2.0*w_i*1000.0);					//[m/s] Spouting velocity
	double U_tip = ms_des_solved.m_nu_design*C_s;		//[m/s] Tip speed
	ms_des_solved.m_D_rotor = U_tip / (0.5*ms_des_solved.m_N_design*0.104719755);	//[m]
	ms_des_solved.m_A_nozzle = ms_des_par.m_m_dot / (C_s*ms_des_par.m_D_in);		//[m^2]

	// Set other turbine variables
	ms_des_solved.m_w_tip_ratio = U_tip / ssnd_in;				//[-]
	ms_des_solved.m_eta = (ms_des_par.m_h_in - ms_des_par.m_h_out) / w_i;	//[-] Isentropic efficiency
}

void C_turbine::off_design_turbine(double T_in, double P_in, double P_out, double N, int & error_code, double & m_dot, double & T_out)
{
	/* 9.4.14: code from John Dyreby, converted to C++ by Ty Neises
	! Solve for the outlet state of 'turb' given its inlet conditions, outlet pressure, and shaft speed.
	!
	! Inputs:
	!   turb -- a Turbine object, with design-point values and sizing set
	!   T_in -- turbine inlet temperature (K)
	!   P_in -- turbine inlet pressure (kPa)
	!   P_out -- turbine outlet pressure (kPa)
	!   N -- shaft speed of turbine (rpm)
	!
	! Outputs:
	!   error_trace -- an ErrorTrace object
	!   m_dot -- allowable mass flow rate through the turbine (kg/s)
	!   T_out -- turbine outlet temperature (K)
	!
	! Notes:
	!   1) This subroutine also sets the following values in 'turb': nu, eta, m_dot, w, w_tip_ratio */

	CO2_state co2_props;

	// Get properties at turbine inlet
	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double D_in = co2_props.dens;
	double h_in = co2_props.enth;
	double s_in = co2_props.entr;
	double ssnd_in = co2_props.ssnd;

	prop_error_code = CO2_PS(P_out, s_in, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	double h_s_out = co2_props.enth;

	// Apply the radial turbine equations for efficiency
	double C_s = sqrt(2.0*(h_in - h_s_out)*1000.0);				//[m/s] spouting velocity
	double U_tip = ms_des_solved.m_D_rotor*0.5*N*0.104719755;	//[m/s] tip speed
	ms_od_solved.m_nu = U_tip / C_s;							//[-] ratio of tip speed to spouting velocity

	double eta_0 = (((1.0626*ms_od_solved.m_nu - 3.0874)*ms_od_solved.m_nu + 1.3668)*ms_od_solved.m_nu + 1.3567)*ms_od_solved.m_nu + 0.179921180;
	eta_0 = std::max(eta_0, 0.0);
	eta_0 = std::min(eta_0, 1.0);
	ms_od_solved.m_eta = eta_0*ms_des_solved.m_eta;		//[-] Actual turbine efficiency

	// Calculate the outlet state and allowable mass flow rate
	double h_out = h_in - ms_od_solved.m_eta*(h_in - h_s_out);		//[kJ/kg] Enthalpy at turbine outlet
	prop_error_code = CO2_PH(P_out, h_out, &co2_props);
	if (prop_error_code != 0)
	{
		error_code = prop_error_code;
		return;
	}
	T_out = co2_props.temp;

	m_dot = C_s*ms_des_solved.m_A_nozzle*D_in;			//[kg/s] Mass flow rate through turbine
	ms_od_solved.m_w_tip_ratio = U_tip / ssnd_in;		//[-] Ratio of the tip speed to the local (turbine inlet) speed of sound
	ms_od_solved.m_N = N;
	ms_od_solved.m_W_dot_out = m_dot*(h_in - h_out);		//[kW] Turbine power output
}

void C_turbine::od_turbine_at_N_des(double T_in, double P_in, double P_out, int & error_code, double & m_dot, double & T_out)
{
	double N = ms_des_solved.m_N_design;		//[rpm]

	off_design_turbine(T_in, P_in, P_out, N, error_code, m_dot, T_out);

	return;
}

int C_comp_single_stage::design_given_shaft_speed(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/,
	double N_rpm /*rpm*/, double eta_isen /*-*/, double & P_out /*kPa*/, double & T_out /*K*/, double & tip_ratio /*-*/)
{
	CO2_state co2_props;

	// Get inlet state
	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);
	if (prop_error_code != 0)
	{
		return prop_error_code;
	}
	double h_in = co2_props.enth;	//[kJ/kg]
	double s_in = co2_props.entr;	//[kJ/kg-K]
	double rho_in = co2_props.dens;	//[kg/m^3]

	// Convert shaft speed to rad/s
	double N_rad_s = N_rpm / 9.549296590;		//[rad/s]

	// Solve for the diameter that gives the design flow coefficient
	double D_rotor = std::pow(m_dot / (m_snl_phi_design * rho_in * 0.5 * N_rad_s), 1.0 / 3.0);		//[m]

	// Calculate psi at the design-point phi using Horner's method
	double psi_design = ((((-498626.0*m_snl_phi_design) + 53224.0) * m_snl_phi_design - 2505.0) * m_snl_phi_design + 54.6) *
		m_snl_phi_design + 0.04049;		// from dimensionless modified head curve(at design - point, psi and modified psi are equal)

	// Solve for idea head
	double U_tip = 0.5 * D_rotor * N_rad_s;		//[m/s]
	double w_i = psi_design * std::pow(U_tip, 2) * 0.001;		//[kJ/kg]

	// Solve for isentropic outlet enthalpy
	double h_out_isen = h_in + w_i;		//[kJ/kg]

	// Get isentropic outlet state
	prop_error_code = CO2_HS(h_out_isen, s_in, &co2_props);
	if (prop_error_code != 0)
	{
		return prop_error_code;
	}
	P_out = co2_props.pres;		//[kPa]

	// Get actual outlet state
	double h_out = h_in + w_i / eta_isen;	//[kJ/kg]
	prop_error_code = CO2_PH(P_out, h_out, &co2_props);
	if (prop_error_code != 0)
	{
		return prop_error_code;
	}
	T_out = co2_props.temp;		//[K]
	double ssnd_out = co2_props.ssnd;	//[m/s]

	// Solve for tip ratio
	tip_ratio = U_tip / ssnd_out;

	ms_des_solved.m_T_in = T_in;	//[K]
	ms_des_solved.m_P_in = P_in;	//[kPa]
	ms_des_solved.m_D_in = rho_in;	//[kg/m^3]
	ms_des_solved.m_h_in = h_in;	//[kJ/kg]
	ms_des_solved.m_s_in = s_in;	//[kJ/kg-K]

	ms_des_solved.m_T_out = T_out;	//[K]
	ms_des_solved.m_P_out = P_out;	//[kPa]
	ms_des_solved.m_h_out = h_out;	//[kJ/kg]
	ms_des_solved.m_D_out = co2_props.dens;	//[kg/m^3]

	ms_des_solved.m_m_dot = m_dot;	//[kg/s]

	ms_des_solved.m_D_rotor = D_rotor;	//[m]
	ms_des_solved.m_N_design = N_rpm;	//[rpm]
	ms_des_solved.m_tip_ratio = tip_ratio;	//[-]
	ms_des_solved.m_eta_design = eta_isen;		//[-]

	ms_des_solved.m_phi_des = m_snl_phi_design;
	ms_des_solved.m_phi_surge = m_snl_phi_min;
	ms_des_solved.m_phi_max = m_snl_phi_max;

	return 0;
}

int C_comp_single_stage::design_single_stage_comp(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/,
	double T_out /*K*/, double P_out /*K*/)
{
	CO2_state in_props;
	int prop_err_code = CO2_TP(T_in, P_in, &in_props);
	if (prop_err_code != 0)
	{
		return -1;
	}
	double s_in = in_props.entr;	//[kJ/kg-K]
	double h_in = in_props.enth;	//[kJ/kg]
	double rho_in = in_props.dens;	//[kg/m^3]

	CO2_state isen_out_props;
	prop_err_code = CO2_PS(P_out, s_in, &isen_out_props);
	if (prop_err_code != 0)
	{
		return -1;
	}
	double h_isen_out = isen_out_props.enth;	//[kJ/kg]

	CO2_state out_props;
	prop_err_code = CO2_TP(T_out, P_out, &out_props);
	if (prop_err_code != 0)
	{
		return -1;
	}
	double h_out = out_props.enth;

	// Calculate psi at the design-point phi using Horner's method
	double psi_design = ((((-498626.0*m_snl_phi_design) + 53224.0) * m_snl_phi_design - 2505.0) * m_snl_phi_design + 54.6) *
		m_snl_phi_design + 0.04049;		// from dimensionless modified head curve(at design - point, psi and modified psi are equal)

	// Determine required size and speed of compressor
	double w_i = h_isen_out - h_in;						//[kJ/kg] positive isentropic specific work of compressor
	double U_tip = sqrt(1000.0*w_i / psi_design);		//[m/s]
	double D_rotor = sqrt(m_dot / (m_snl_phi_design*rho_in*U_tip));
	double N_rad_s = U_tip*2.0 / D_rotor;				//[rad/s] shaft speed

	double ssnd_out = out_props.ssnd;	//[m/s]

	// Solve for tip ratio
	double tip_ratio = U_tip / ssnd_out;

	ms_des_solved.m_T_in = T_in;	//[K]
	ms_des_solved.m_P_in = P_in;	//[kPa]
	ms_des_solved.m_D_in = rho_in;	//[kg/m^3]
	ms_des_solved.m_h_in = h_in;	//[kJ/kg]
	ms_des_solved.m_s_in = s_in;	//[kJ/kg-K]

	ms_des_solved.m_T_out = T_out;	//[K]
	ms_des_solved.m_P_out = P_out;	//[kPa]
	ms_des_solved.m_h_out = h_out;	//[kJ/kg]
	ms_des_solved.m_D_out = out_props.dens;	//[kg/m^3]

	ms_des_solved.m_m_dot = m_dot;	//[kg/s]

	ms_des_solved.m_D_rotor = D_rotor;	//[m]
	ms_des_solved.m_N_design = N_rad_s * 9.549296590;			//[rpm] shaft speed
	ms_des_solved.m_tip_ratio = tip_ratio;	//[-]
	ms_des_solved.m_eta_design = (h_isen_out - h_in) / (h_out - h_in);		//[-]

	ms_des_solved.m_phi_des = m_snl_phi_design;
	ms_des_solved.m_phi_surge = m_snl_phi_min;
	ms_des_solved.m_phi_max = m_snl_phi_max;

	return 0;
}

int C_comp_single_stage::off_design_given_N(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/, double N_rpm /*rpm*/,
	double & T_out /*K*/, double & P_out /*kPa*/)
{
	CO2_state co2_props;

	ms_od_solved.m_N = N_rpm;		//[rpm]

	// Fully define the inlet state of the compressor
	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);
	if (prop_error_code != 0)
	{
		return prop_error_code;
	}
	double rho_in = co2_props.dens;	//[kg/m^3]
	double h_in = co2_props.enth;	//[kJ/kg]
	double s_in = co2_props.entr;	//[kJ/kg-K]

	// Calculate the modified flow and head coefficients and efficiency
	double U_tip = ms_des_solved.m_D_rotor*0.5*ms_od_solved.m_N*0.104719755;				//[m/s]
	double phi = m_dot / (rho_in*U_tip*pow(ms_des_solved.m_D_rotor, 2));	//[-]
	if (phi < m_snl_phi_min)
	{
		ms_od_solved.m_surge = true;
		phi = m_snl_phi_min;
	}
	else
		ms_od_solved.m_surge = false;

	double phi_star = phi*pow(ms_od_solved.m_N / ms_des_solved.m_N_design, 0.2);		//[-] modified flow coefficient
	double psi_star = ((((-498626.0*phi_star) + 53224.0)*phi_star - 2505.0)*phi_star + 54.6)*phi_star + 0.04049;	// from dimensionless modified head curve
	double eta_star = ((((-1.638e6*phi_star) + 182725.0)*phi_star - 8089.0)*phi_star + 168.6)*phi_star - 0.7069;	// from dimensionless modified efficiency curve
	double psi = psi_star / pow(ms_des_solved.m_N_design / ms_od_solved.m_N, pow(20.0*phi_star, 3.0));
	double eta_0 = eta_star*1.47528 / pow(ms_des_solved.m_N_design / ms_od_solved.m_N, pow(20.0*phi_star, 5.0));		// Efficiency is normalized so it equals 1.0 at snl_phi_design
	ms_od_solved.m_eta = std::max(eta_0*ms_des_solved.m_eta_design, 0.0);		//[-] Actual compressor efficiency, not allowed to go negative

	// Check that the specified mass flow rate is possible with the compressor's current shaft speed
	if (psi <= 0.0)
	{
		return 1;
	}

	// Calculate the compressor outlet state
	double dh_s = psi * pow(U_tip, 2.0) * 0.001;			//[kJ/kg] Ideal enthalpy rise in compressor, from definition of head coefficient
	double dh = dh_s / ms_od_solved.m_eta;					//[kJ/kg] Actual enthalpy rise in compressor
	double h_s_out = h_in + dh_s;							//[kJ/kg] Ideal enthalpy at compressor outlet
	double h_out = h_in + dh;								//[kJ/kg] Actual enthalpy at compressor outlet

	// Get the compressor outlet pressure
	prop_error_code = CO2_HS(h_s_out, s_in, &co2_props);
	if (prop_error_code != 0)
	{
		return 2;
	}
	P_out = co2_props.pres;

	// Determine compressor outlet temperature and speed of sound
	prop_error_code = CO2_PH(P_out, h_out, &co2_props);
	if (prop_error_code != 0)
	{
		return 2;
	}
	T_out = co2_props.temp;
	double ssnd_out = co2_props.ssnd;

	// Set a few compressor variables
	ms_od_solved.m_P_in = P_in;		//[kPa]
	ms_od_solved.m_h_in = h_in;		//[kJ/kg]
	ms_od_solved.m_T_in = T_in;		//[K]
	ms_od_solved.m_s_in = s_in;		//[kJ/kg-K]

	ms_od_solved.m_P_out = P_out;	//[kPa]
	ms_od_solved.m_h_out = h_out;	//[kJ/kg]
	ms_od_solved.m_T_out = T_out;	//[K]
	ms_od_solved.m_s_out = co2_props.entr;	//[kJ/kg-K]

	ms_od_solved.m_phi = phi;
	ms_od_solved.m_surge_safety = phi / m_snl_phi_min;	//[-] If > 1, then not in surge
	ms_od_solved.m_w_tip_ratio = U_tip / ssnd_out;
	ms_od_solved.m_W_dot_in = m_dot*(h_out - h_in);	//[kWe]

	return 0;
}

int C_comp_single_stage::calc_N_from_phi(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/, double phi_in /*-*/, double & N_rpm /*rpm*/)
{
	CO2_state co2_props;

	// Fully define the inlet state of the compressor
	int prop_error_code = CO2_TP(T_in, P_in, &co2_props);
	if (prop_error_code != 0)
	{
		return prop_error_code;
	}
	double rho_in = co2_props.dens;	//[kg/m^3]
	double U_tip = m_dot / (phi_in*rho_in*std::pow(ms_des_solved.m_D_rotor, 2));		//[m/s]
	N_rpm = (U_tip*2.0 / ms_des_solved.m_D_rotor)*9.549296590;		//[rpm]

	return 0;
}

int C_comp_multi_stage::C_MEQ_eta_isen__h_out::operator()(double eta_isen /*-*/, double *h_comp_out /*kJ/kg*/)
{
	C_MEQ_N_rpm__P_out c_stages(mpc_multi_stage, m_T_in, m_P_in, m_m_dot, eta_isen);
	C_monotonic_eq_solver c_solver(c_stages);

	// Set lowr bound
	double N_rpm_lower = 0.0;
	double N_rpm_upper = std::numeric_limits<double>::quiet_NaN();

	// Generate guess values
	double N_rpm_guess_1 = 3000.0;
	double N_rpm_guess_2 = 30000.0;

	c_solver.settings(1.E-4, 50, N_rpm_lower, N_rpm_upper, true);

	// Now solve for the shaft speed
	double N_rpm_solved = std::numeric_limits<double>::quiet_NaN();
	double tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

	int N_rpm_code = 0;
	try
	{
		N_rpm_code = c_solver.solve(N_rpm_guess_1, N_rpm_guess_2, m_P_out, N_rpm_solved, tol_solved, iter_solved);
	}
	catch (C_csp_exception)
	{
		throw(C_csp_exception("C_comp_multi_stage::C_MEQ_eta_isen__h_out threw an exception"));
	}

	if (N_rpm_code != C_monotonic_eq_solver::CONVERGED)
	{
		if (!(N_rpm_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.01))
		{
			throw(C_csp_exception("C_comp_multi_stage::C_MEQ_eta_isen__h_out failed to converge within a reasonable tolerance"));
		}
	}

	int n_stages = mpc_multi_stage->mv_stages.size();

	*h_comp_out = mpc_multi_stage->mv_stages[n_stages - 1].ms_des_solved.m_h_out;	//[kJ/kg]

	return 0;
}

int C_comp_multi_stage::C_MEQ_N_rpm__P_out::operator()(double N_rpm /*rpm*/, double *P_comp_out /*kPa*/)
{
	int n_stages = mpc_multi_stage->mv_stages.size();

	double T_in = m_T_in;	//[K]
	double P_in = m_P_in;	//[kPa]

	double P_out = std::numeric_limits<double>::quiet_NaN();
	double T_out = std::numeric_limits<double>::quiet_NaN();
	double tip_ratio = std::numeric_limits<double>::quiet_NaN();

	for (int i = 0; i < n_stages; i++)
	{
		if (i > 0)
		{
			T_in = T_out;	//[K]
			P_in = P_out;	//[kPa]
		}

		mpc_multi_stage->mv_stages[i].design_given_shaft_speed(T_in, P_in, m_m_dot, N_rpm, m_eta_isen, P_out, T_out, tip_ratio);
	}

	*P_comp_out = P_out;	//[kPa]

	return 0;
}

int C_comp_multi_stage::design_given_outlet_state(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/,
	double T_out /*K*/, double P_out /*K*/)
{
	mv_stages.resize(1);
	mv_stages[0].design_single_stage_comp(T_in, P_in, m_dot, T_out, P_out);

	double max_calc_tip_speed = mv_stages[0].ms_des_solved.m_tip_ratio;

	if (mv_stages[0].ms_des_solved.m_tip_ratio > 0.9)
	{
		CO2_state co2_props;

		double h_in = mv_stages[0].ms_des_solved.m_h_in;
		double s_in = mv_stages[0].ms_des_solved.m_s_in;

		int prop_err_code = CO2_PS(P_out, s_in, &co2_props);
		if (prop_err_code != 0)
		{
			return -1;
		}
		double h_out_isen = co2_props.enth;

		double h_out = mv_stages[0].ms_des_solved.m_h_out;

		double eta_isen_total = (h_out_isen - h_in) / (h_out - h_in);

		bool is_add_stages = true;
		int n_stages = 1;

		while (is_add_stages)
		{
			n_stages++;

			mv_stages.resize(n_stages);

			C_MEQ_eta_isen__h_out c_stages(this, T_in, P_in, P_out, m_dot);
			C_monotonic_eq_solver c_solver(c_stages);

			// Set bounds on isentropic efficiency
			double eta_isen_lower = 0.1;
			double eta_isen_upper = 1.0;

			// Generate guess values
			double eta_isen_guess_1 = eta_isen_total;
			double eta_isen_guess_2 = 0.95*eta_isen_total;

			c_solver.settings(1.E-4, 50, eta_isen_lower, eta_isen_upper, true);

			// Now solve for the isentropic efficiency for each stage that results in the total compressor design isentropic efficiency
			double eta_isen_solved = std::numeric_limits<double>::quiet_NaN();
			double tol_solved = std::numeric_limits<double>::quiet_NaN();
			int iter_solved = -1;

			int eta_isen_code = 0;

			try
			{
				eta_isen_code = c_solver.solve(eta_isen_guess_1, eta_isen_guess_2, h_out, eta_isen_solved, tol_solved, iter_solved);
			}
			catch (C_csp_exception)
			{
				throw(C_csp_exception("C_comp_multi_stage::design_given_outlet_state threw an exception"));
			}

			if (eta_isen_code != C_monotonic_eq_solver::CONVERGED)
			{
				if (!(eta_isen_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.01))
				{
					throw(C_csp_exception("C_comp_multi_stage::design_given_outlet_state failed to converge within a reasonable tolerance"));
				}
			}

			max_calc_tip_speed = 0.0;
			for (int i = 0; i < n_stages; i++)
			{
				max_calc_tip_speed = std::max(max_calc_tip_speed, mv_stages[i].ms_des_solved.m_tip_ratio);
			}

			if (max_calc_tip_speed < 0.9)
			{
				is_add_stages = false;
			}

			if (n_stages > 20)
			{
				return -1;
			}
		}
	}

	int n_stages = mv_stages.size();

	ms_des_solved.m_T_in = T_in;	//[K]
	ms_des_solved.m_P_in = P_in;	//[kPa]
	ms_des_solved.m_D_in = mv_stages[0].ms_des_solved.m_D_in;	//[kg/m^3]
	ms_des_solved.m_s_in = mv_stages[0].ms_des_solved.m_s_in;	//[kJ/kg-K]
	ms_des_solved.m_h_in = mv_stages[0].ms_des_solved.m_h_in;	//[kJ/kg]

	ms_des_solved.m_T_out = mv_stages[n_stages - 1].ms_des_solved.m_T_out;	//[K]
	ms_des_solved.m_P_out = mv_stages[n_stages - 1].ms_des_solved.m_P_out;	//[kPa]
	ms_des_solved.m_h_out = mv_stages[n_stages - 1].ms_des_solved.m_h_out;	//[kJ/kg]
	ms_des_solved.m_D_out = mv_stages[n_stages - 1].ms_des_solved.m_D_out;	//[kg/m^3]

	ms_des_solved.m_m_dot = m_dot;					//[kg/s]

	ms_des_solved.m_N_design = mv_stages[n_stages - 1].ms_des_solved.m_N_design;		//[rpm]
	ms_des_solved.m_phi_des = mv_stages[0].ms_des_solved.m_phi_des;		//[-]
	ms_des_solved.m_w_tip_ratio = max_calc_tip_speed;					//[-]
	ms_des_solved.m_n_stages = n_stages;								//[-]
	ms_des_solved.m_D_rotor = mv_stages[0].ms_des_solved.m_D_rotor;		//[m]
	ms_des_solved.m_phi_surge = mv_stages[0].m_snl_phi_min;				//[-]

	return 0;
}

void C_comp_multi_stage::off_design_at_N_des(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/,
	int & error_code, double & T_out /*K*/, double & P_out /*kPa*/)
{
	double N = ms_des_solved.m_N_design;	//[rpm]

	off_design_given_N(T_in, P_in, m_dot, N, error_code, T_out, P_out);
}

void C_comp_multi_stage::off_design_given_N(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/, double N_rpm /*rpm*/,
	int & error_code, double & T_out /*K*/, double & P_out /*kPa*/)
{
	int n_stages = mv_stages.size();

	double T_stage_in = T_in;	//[K]
	double P_stage_in = P_in;	//[kPa]

	double T_stage_out = std::numeric_limits<double>::quiet_NaN();
	double P_stage_out = std::numeric_limits<double>::quiet_NaN();

	double tip_ratio_max = 0.0;
	bool is_surge = false;
	double surge_safety_min = 10.0;
	double phi_min = 10.0;

	for (int i = 0; i < n_stages; i++)
	{
		if (i > 0)
		{
			T_stage_in = T_stage_out;
			P_stage_in = P_stage_out;
		}

		error_code = mv_stages[i].off_design_given_N(T_stage_in, P_stage_in, m_dot, N_rpm, T_stage_out, P_stage_out);
		if (error_code != 0)
		{
			return;
		}

		if (mv_stages[i].ms_od_solved.m_w_tip_ratio > tip_ratio_max)
		{
			tip_ratio_max = mv_stages[i].ms_od_solved.m_w_tip_ratio;
		}

		if (mv_stages[i].ms_od_solved.m_surge)
		{
			is_surge = true;
		}

		if (mv_stages[i].ms_od_solved.m_surge_safety < surge_safety_min)
		{
			surge_safety_min = mv_stages[i].ms_od_solved.m_surge;
		}

		phi_min = std::min(phi_min, mv_stages[i].ms_od_solved.m_phi);
	}

	P_out = mv_stages[n_stages - 1].ms_od_solved.m_P_out;		//[kPa]
	T_out = mv_stages[n_stages - 1].ms_od_solved.m_T_out;	//[K]

	double h_in = mv_stages[0].ms_od_solved.m_h_in;					//[kJ/kg]
	double s_in = mv_stages[0].ms_od_solved.m_s_in;					//[kJ/kg-K]

	CO2_state co2_props;
	int prop_err_code = CO2_PS(P_out, s_in, &co2_props);
	if (prop_err_code != 0)
	{
		error_code = prop_err_code;
		return;
	}
	double h_out_isen = co2_props.enth;

	double h_out = mv_stages[n_stages - 1].ms_od_solved.m_h_out;	//[kJ/kg]

	ms_od_solved.m_P_in = P_in;
	ms_od_solved.m_T_in = T_in;

	ms_od_solved.m_P_out = P_out;
	ms_od_solved.m_T_out = T_out;

	ms_od_solved.m_surge = is_surge;
	ms_od_solved.m_eta = (h_out_isen - h_in) / (h_out - h_in);
	ms_od_solved.m_phi = mv_stages[0].ms_od_solved.m_phi;
	ms_od_solved.m_w_tip_ratio = tip_ratio_max;

	ms_od_solved.m_N = N_rpm;

	ms_od_solved.m_W_dot_in = m_dot*(h_out - h_in);
	ms_od_solved.m_surge_safety = surge_safety_min;

}

int C_comp_multi_stage::C_MEQ_phi_od__P_out::operator()(double phi_od /*-*/, double *P_comp_out /*kPa*/)
{
	int error_code = 0;
	double N_rpm = std::numeric_limits<double>::quiet_NaN();
	error_code = mpc_multi_stage->mv_stages[0].calc_N_from_phi(m_T_in, m_P_in, m_m_dot, phi_od, N_rpm);
	if (error_code != 0)
	{
		*P_comp_out = std::numeric_limits<double>::quiet_NaN();
		return error_code;
	}

	double T_out = std::numeric_limits<double>::quiet_NaN();
	error_code = 0;
	mpc_multi_stage->off_design_given_N(m_T_in, m_P_in, m_m_dot, N_rpm, error_code, T_out, *P_comp_out);

	if (error_code != 0)
	{
		*P_comp_out = std::numeric_limits<double>::quiet_NaN();
		return error_code;
	}

	return 0;
}

void C_comp_multi_stage::off_design_given_P_out(double T_in /*K*/, double P_in /*kPa*/, double m_dot /*kg/s*/,
	double P_out /*kPa*/, int & error_code, double & T_out /*K*/)
{
	// Apply 1 var solver to find the phi that results in a converged recompressor
	C_MEQ_phi_od__P_out c_rc_od(this, T_in, P_in, m_dot);
	C_monotonic_eq_solver c_rd_od_solver(c_rc_od);

	// Set upper and lower bounds
	double phi_upper = mv_stages[0].ms_des_solved.m_phi_max;
	double phi_lower = mv_stages[0].ms_des_solved.m_phi_surge;

	// Generate first x-y pair
	double phi_guess_lower = ms_des_solved.m_phi_des;
	double P_solved_phi_guess_lower = std::numeric_limits<double>::quiet_NaN();
	int test_code = c_rd_od_solver.test_member_function(phi_guess_lower, &P_solved_phi_guess_lower);
	if (test_code != 0)
	{
		for (int i = 1; i < 9; i++)
		{
			phi_guess_lower = ms_des_solved.m_phi_des*(10 - i) / 10.0 + mv_stages[0].ms_des_solved.m_phi_max*i / 10.0;
			test_code = c_rd_od_solver.test_member_function(phi_guess_lower, &P_solved_phi_guess_lower);
			if (test_code == 0)
				break;
		}
	}
	if (test_code != 0)
	{
		// Can't find a RC phi guess value that returns an outlet pressure
		error_code = -20;
		return;
	}
	C_monotonic_eq_solver::S_xy_pair phi_pair_lower;
	phi_pair_lower.x = phi_guess_lower;
	phi_pair_lower.y = P_solved_phi_guess_lower;

	// Generate second x-y pair
	double phi_guess_upper = phi_guess_lower*0.5 + mv_stages[0].ms_des_solved.m_phi_max*0.5;
	double P_solved_phi_guess_upper = std::numeric_limits<double>::quiet_NaN();
	test_code = c_rd_od_solver.test_member_function(phi_guess_upper, &P_solved_phi_guess_upper);
	if (test_code != 0)
	{
		for (int i = 6; i < 10; i++)
		{
			phi_guess_upper = phi_guess_lower*i / 10.0 + mv_stages[0].ms_des_solved.m_phi_max*(10 - i) / 10.0;
			test_code = c_rd_od_solver.test_member_function(phi_guess_upper, &P_solved_phi_guess_upper);
			if (test_code == 0)
				break;
		}
		if (test_code != 0 && phi_guess_lower == ms_des_solved.m_phi_des)
		{
			for (int i = 6; i < 10; i++)
			{
				phi_guess_upper = phi_guess_lower*i / 10.0 + ms_des_solved.m_phi_surge*(10 - i) / 10.0;
				test_code = c_rd_od_solver.test_member_function(phi_guess_upper, &P_solved_phi_guess_upper);
				if (test_code == 0)
					break;
			}
		}
	}
	if (test_code != 0)
	{
		// Can't find a RC 2nd guess value (which, if we've found a first, means the solution space is really small?)
		error_code = -20;
		return;
	}
	C_monotonic_eq_solver::S_xy_pair phi_pair_upper;
	phi_pair_upper.x = phi_guess_upper;
	phi_pair_upper.y = P_solved_phi_guess_upper;

	// Set solver settings
	c_rd_od_solver.settings(0.001, 50, phi_lower, phi_upper, true);

	// Now, solve for the flow coefficient
	double phi_solved, tol_solved;
	phi_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

	int phi_code = 0;
	try
	{
		phi_code = c_rd_od_solver.solve(phi_pair_lower, phi_pair_upper, P_out, phi_solved, tol_solved, iter_solved);
	}
	catch (C_csp_exception)
	{
		error_code = -1;
		return;
	}

	if (phi_code != C_monotonic_eq_solver::CONVERGED)
	{
		int n_call_history = (int)c_rd_od_solver.get_solver_call_history()->size();

		if (n_call_history > 0)
			error_code = -(*(c_rd_od_solver.get_solver_call_history()))[n_call_history - 1].err_code;

		if (error_code == 0)
		{
			error_code = phi_code;
		}

		return;
	}

	T_out = ms_od_solved.m_T_out;		//[K]
}
