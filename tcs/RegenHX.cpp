#include "RegenHX.h"
#include <iostream>
#include <ctime>

//Root is at YOUR_SAM_DIR\ssc\examples

#ifdef _DEBUG

string const RegenHX::solver_exit_modes[] = {
	"REL_TOL_WITH_0_TARGET",
	"EQUAL_GUESS_VALUES",
	"NO_SOLUTION",

	"CONVERGED",

	"SLOPE_POS_NO_NEG_ERR",
	"SLOPE_NEG_NO_NEG_ERR",

	"SLOPE_POS_NO_POS_ERR",
	"SLOPE_NEG_NO_POS_ERR",

	"SLOPE_POS_BOTH_ERRS",
	"SLOPE_NEG_BOTH_ERRS",

	"MAX_ITER_SLOPE_POS_NO_NEG_ERR",
	"MAX_ITER_SLOPE_NEG_NO_NEG_ERR",

	"MAX_ITER_SLOPE_POS_NO_POS_ERR",
	"MAX_ITER_SLOPE_NEG_NO_POS_ERR",

	"MAX_ITER_SLOPE_POS_BOTH_ERRS",
	"MAX_ITER_SLOPE_NEG_BOTH_ERRS"
};

string const RegenHX::DEBUG_LOG_FILEPATH = "./";

#endif

string const RegenHX::LOG_FILEPATH = "./RegenHX_LOG.txt";
string const RegenHX::PROPERTY_FILES = "../tcs/PropertyFiles/";
string const RegenHX::SPHERES_RP_TABLE_PATH = PROPERTY_FILES + "Spheres_RP.csv";
string const RegenHX::BALANCED_REGENERATOR_TABLE_PATH = PROPERTY_FILES + "balanced-regenerator.csv";
string const RegenHX::FATIGUE_TABLE_PATH = PROPERTY_FILES + "fatigue.csv";

RegenHX::RegenHX()
{
	loadTables();
	initializeLOGs();
}

RegenHX::~RegenHX()
{
	LOG.flush();
	LOG.close();

#ifdef _DEBUG
	DEBUG_LOG.flush();
	DEBUG_LOG.close();
#endif

}

void RegenHX::loadTables() {
	bedMaterialTable = new LookupTable_1D(PROPERTY_FILES + bedMaterialName + ".csv");
	shellMaterialTable = new LookupTable_1D(PROPERTY_FILES + shellMaterialName + ".csv");
	insulationMaterialTable = new LookupTable_1D(PROPERTY_FILES + insulationMaterialName + ".csv");
	regeneratorTable = new LookupTable_2D(BALANCED_REGENERATOR_TABLE_PATH);
	spheresRPTable = new LookupTable_1D(SPHERES_RP_TABLE_PATH);
	fatigueTable = new LookupTable_1D(FATIGUE_TABLE_PATH);
}

/*
Re	[-]
f	[-]
j_H [-]
*/
void RegenHX::packedspheresNdFit(double Re, double * f, double * j_H)
{
	if (Re < 20 || Re >  50000) {
		(*f) = -1;
		(*j_H) = -1;
		throw invalid_argument("Re should be between 20 and 50,000!");
	}

	//Re is coloumn 0, f is coloumn 1
	(*f) = spheresRPTable->getValue("f", "Re", Re);
	(*j_H) = 0.23*pow(Re, -0.3);
}

/*
NTU [-]
C_1 [W/K]
C_2	[W/K]
However any matching units of C_1 and C_2 are okay.
*/
double RegenHX::hx(double NTU, double C_1, double C_2)
{
	if (C_1 <= 0 || C_2 <= 0) {
		throw invalid_argument("C_1 and C_2 should be greater than 0!");
	}

	double U = C_1 / C_2;	//[-]

	if (U > 1000) {
		throw invalid_argument("No value provided for U > 1000!");
	}

	return regeneratorTable->getValue(NTU, U);
}

/*
T	[K]
P	[kPa]
rho	[kg/m^3]
mu	[Pa-s]
k	[W/m-K]
Pr	[-]
Cp	[kJ/kg-K]
*/
void RegenHX::getpropsregenFitCO2(double T, double P, double * rho, double * mu, double * k, double * Pr, double * Cp)
{
	CO2_state CO2States;

	int error = CO2_TP(T, P, &CO2States);

	if (error != 0)
	{
		throw invalid_argument("Error fixing CO2 state with T and P provided!");
	}

	*Cp = CO2States.cp;
	*rho = CO2States.dens;
	*mu = CO2_visc(*rho, T) / 1000000; //Convert from uPa-s to Pa-s
	*k = CO2_cond(*rho, T);

	*Pr = (*Cp)*(*mu) / (*k)* 1000.0;	//Convert MPa to kPa
}

/*
m_dot		//[kg/s]
d			//[m]
A_fr		//[m^2]
L			//[m]
T			//[K]
P			//[kPa]
porosity	//[-]
f			//[-]
h			//[W/m^2-K]
NTU			//[-]
DP			//[kPa]
*/
void RegenHX::packedspheresFitCO2(double m_dot, double d, double A_fr, double L, double T, double P, double porosity, double * f, double * h, double * NTU, double * DP)
{
	double rho, mu, k, Pr, Cp;
	try {
		getpropsregenFitCO2(T, P, &rho, &mu, &k, &Pr, &Cp);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	double r_h = porosity*d / (6.0*(1 - porosity));	//[m]
	double G = m_dot / (porosity*A_fr);	//[kg/m^2-s]
	double Re = G*(4 * r_h) / mu;	//[-]

	double j_H;	//[-]
	try {
		packedspheresNdFit(Re, f, &j_H);
	}
	catch (const invalid_argument& e) {
		throw e;
	}


	double St = j_H / pow(Pr, (2.0 / 3.0));	//[-]
	*NTU = St*L / r_h;	//[-]
	*h = j_H*G*Cp / pow(Pr, (2.0 / 3.0)) * 1000.0;	//Convert kW to W. [W]

	
	double V = m_dot / rho / (A_fr);
	Re = d * rho*V / mu;
	double HELPER = (1 - porosity);
	double Re_m = Re / HELPER;
	double f_1L = 136 / pow(HELPER, 0.38);
	double f_1T = 29 / (pow(HELPER, 1.45) * pow(porosity, 2));
	double f_2 = 1.87 * pow(porosity, 0.75) / pow(HELPER, 0.26);
	double q = exp(-pow(porosity, 2) * HELPER / 12.6 * Re_m);
	double f_p = (q * f_1L / Re_m + (1 - q) * (f_2 + f_1T / Re_m)) * HELPER / pow(porosity, 3);
	*DP = f_p * rho * pow(V, 2) * L / d / 1000.0; //Pa to kPa


	//*DP = *f*L*pow(G, 2) / (2 * rho*r_h) / 1000.0;	//Convert Pa to kPa. [kPa]
}

//Explicitly calculates Q_dot_a and Q_dot_a_calc based on provided T_H_out
void RegenHX::CalculateThermoAndPhysicalModels()
{
	//							[SECTION 1]
	//Pressure at hot exit
	P_H_out = P_H - dP_H;	//[kPa]
	//Pressure at cold exit
	P_C_out = P_C - dP_C;	//[kPa]
	
	if (P_C_out <= N_co2_props::P_lower_limit || P_H_out <= N_co2_props::P_lower_limit || P_C_out >= N_co2_props::P_upper_limit || P_H_out >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Outlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) +"[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}
	if (T_H_in <= N_co2_props::T_lower_limit || T_C_in <= N_co2_props::T_lower_limit || T_H_in >= N_co2_props::T_upper_limit || T_C_in >= N_co2_props::T_upper_limit) {
		throw invalid_argument("Inlet temperatures are either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	//Enthalpy at hot inlet
	CO2_TP(T_H_in, P_H, &CO2State);
	h_H_in = CO2State.enth;	//[kJ/kg]

	//Enthalpy at cold inlet
	CO2_TP(T_C_in, P_C, &CO2State);
	h_C_in = CO2State.enth;	//[kJ/kg]

	CO2_TP(T_C_in, P_H_out, &CO2State);
	h_H_out_max = CO2State.enth;	//[kJ/kg]

	CO2_TP(T_H_in, P_C_out, &CO2State);
	h_C_out_max = CO2State.enth;	//[kJ/kg]

	//Maximum heat transfer if hot stream exiting at cold temp
	Q_dot_max_H = m_dot_H*(h_H_in - h_H_out_max);	//[kW]

	//Maximum heat transfer if cold stream exiting at hot temp
	Q_dot_max_C = m_dot_C*(h_C_out_max - h_C_in);	//[kW]

	//Actual maximum possible heat transfer
	Q_dot_max = min(Q_dot_max_H, Q_dot_max_C);	//[kW]

	if (T_H_out <= N_co2_props::T_lower_limit || T_H_out >= N_co2_props::T_upper_limit) {
		throw invalid_argument("T_H_out is either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}
	
	CO2_TP(T_H_out, P_H_out, &CO2State);
	h_H_out = CO2State.enth; //[kJ/kg]

	Q_dot_a = m_dot_H*(h_H_in - h_H_out);	//[kW]

	h_C_out = h_C_in + Q_dot_a / m_dot_C;	//[kJ/kg]
	CO2_PH(P_C_out, h_C_out, &CO2State);

	T_C_out = CO2State.temp;	//[K]

	Q_dot = Q_dot_a + Q_dot_loss;	//[kW]

	//							[SECTION 2]

	T_H_f = (T_H_in + T_H_out) / 2;	//[K]
	T_C_f = (T_C_in + T_C_out) / 2;	//[K]

	CO2_TP(T_C_in, P_H, &CO2State);
	h_H_out_max_p = CO2State.enth;	//[kJ/kg]

	CO2_TP(T_H_in, P_C, &CO2State);
	h_C_out_max_p = CO2State.enth;	//[kJ/kg]

	C_p_H = (h_H_in - h_H_out_max_p) / (T_H_in - T_C_in);	//[kJ/kg-K]
	C_p_C = (h_C_out_max_p - h_C_in) / (T_H_in - T_C_in);	//[kJ/kg-K]

	C_dot_H = C_p_H * m_dot_H;	//[kW/K]
	C_dot_C = C_p_C * m_dot_C;	//[kW/K]

	C_dot_min = min(C_dot_H, C_dot_C);	//[kW/K]
	C_dot_max = max(C_dot_H, C_dot_C);	//[kW/K]

	C_R = C_dot_min / C_dot_max;

	//							[SECTION 3]

	A_fr = PI*pow((D_fr) / 2.0, 2);	//[m^2]
	T_H_f = (T_H_in + T_H_out) / 2.0;	//[K]
	T_C_f = (T_C_in + T_C_out) / 2.0;	//[K]
	
	try {
		packedspheresFitCO2(m_dot_H, D_s, A_fr, L, T_H_f, P_H, e_v, &f_H, &h_H, &NTU_H, &dP_H_calc);
		packedspheresFitCO2(m_dot_C, D_s, A_fr, L, T_C_f, P_C, e_v, &f_C, &h_C, &NTU_C, &dP_C_calc);
	} catch (const invalid_argument& e) {
		throw e;
	}
	

	T_f = (T_H_in + T_C_in) / 2.0;	//[K]
	rho_s = bedMaterialTable->getValue("rho", "T", T_f);	//[kg/m^3]
	c_s = bedMaterialTable->getValue("c", "T", T_f);	//[kJ/kg-K]

	V_0 = A_fr * L;	//[m^3]

	m_s = V_0 * (1 - e_v) * rho_s;	//[kg]

	A_s = 6.0 * (1.0 - e_v)*V_0 / D_s;	//[m^2]

	NTU_R = 1.0 / (C_dot_min*1000.0)*(1.0 / ((1.0 / (h_H*A_s)) + (1.0 / (h_C*A_s))));

	UA = 1.0/((1.0 / (h_H*A_s)) + (1.0 / (h_C*A_s))) / 1000.0;	//[kW/K]

	f_0 = 1.0 / P_0;	//[1/s]

	C_m = 2.0 * m_s * c_s * f_0 / C_dot_min;

	NTU_R_e = 2.0 * C_R * NTU_R / (1.0 + C_R);

	C_m_e = 2.0 * C_R*C_m / (1.0 + C_R);

	try {
		epsilon_1 = hx(NTU_R_e * 2, C_dot_min, C_m_e*C_dot_min);
	}
	catch (exception e) {
		throw e;
	}

	X = (1.0 - pow(C_R, 2)) / (2.0 * C_R)*(epsilon_1 / (1.0 - epsilon_1));

	epsilon = (1.0 - exp(-X)) / (1.0 - C_R*exp(-X));

	if (epsilon < 0 || epsilon > 1) {
		throw invalid_argument("Epsilon is not between 0 and 1!");
	}

	Q_dot_calc = epsilon * Q_dot_max;	//[kW]
	Q_dot_a_calc = Q_dot_calc - Q_dot_loss;	//[kW]
	dP_max = max(dP_C, dP_H);
}

void RegenHX::setInletState(double T_H_in, double P_H, double T_C_in, double P_C)
{
	this->T_H_in = T_H_in;
	this->P_H = P_H;
	this->T_C_in = T_C_in;
	this->P_C = P_C;

	LOG << "Inlet states are set. T_H_in = " << this->T_H_in << ", P_H = " << this->P_H << ", T_C_in = " << this->T_C_in << ", P_C = " << this->P_C << endl;
}

void RegenHX::setParameters(string operationMode, double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->operationMode = operationMode;

	if(operationMode == "parallel") {
		this->m_dot_H = m_dot_H / numberOfSets;
		this->m_dot_C = m_dot_C / numberOfSets;
	}
	else if (operationMode == "redundant") {
		this->m_dot_H = m_dot_H;
		this->m_dot_C = m_dot_C;
	}
	this->Q_dot_loss = Q_dot_loss;
	this->P_0 = P_0;
	this->D_s = D_s;
	this->e_v = e_v;

	LOG << "Regenerator parameters are set. operationMode = " << this->operationMode << ", m_dot_H = " << this->m_dot_H << ", m_dot_C = " << this->m_dot_C;
	LOG << ", Q_dot_loss = " << this->Q_dot_loss << ", P_0 = " << this->P_0 << ", D_s = " << this->D_s << ", e_v = " << this->e_v;
	LOG << endl;
}

void RegenHX::setDesignTargets(string targetMode, double dP_max, double targetParameter)
{
	this->targetdP_max = dP_max;
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;

	LOG << "Design targets are set. targetdP_max = " << this->targetdP_max << ", targetMode = " << this->targetMode << ", targetParameter = " << this->targetParameter;
	LOG << endl;
}

void RegenHX::initializeLOGs() {

	LOG.open(LOG_FILEPATH, std::ofstream::app);

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, sizeof(buffer), "%d-%m-%Y %I_%M_%S", timeinfo);
	std::string timestamp(buffer);

	LOG << timestamp << " Run statred!\n";

#ifdef _DEBUG
	DEBUG_LOG.open(DEBUG_LOG_FILEPATH + "LOG_" + timestamp + "_" + std::to_string(clock()) + ".txt", std::ofstream::app);
#endif
}

void RegenHX::initialize(int N_sub_hx)
{
	
}

void RegenHX::solveSystem(double* results)
{
	
	BalanceQdotAs();
	figureOutL();
	calculateWallThickness();
	calculateCost();
	results[0] = epsilon;
	results[1] = costHXTotal;
	results[2] = UA;
	results[3] = T_H_out;
	results[4] = dP_H;
	results[5] = dP_C;
	results[6] = D_fr;
	results[7] = L;
	results[8] = wallThickness;
}

double RegenHX::od_delta_p_cold(double m_dot_c /*kg/s*/)
{
	return ms_des_solved.m_DP_cold_des*pow(m_dot_c / m_dot_C, 1.75);
}

double RegenHX::od_delta_p_hot(double m_dot_h /*kg/s*/)
{
	return ms_des_solved.m_DP_hot_des*pow(m_dot_h / m_dot_H, 1.75);
}

void RegenHX::design_fix_UA_calc_outlet(double Cost_target, double eff_target, double T_c_in, double P_c_in, double m_dot_c, double P_c_out, double T_h_in, double P_h_in, double m_dot_h, double P_h_out, double & q_dot, double & T_c_out, double & T_h_out)
{
	setInletState(T_h_in, P_h_in, T_c_in, P_c_in);
	double Q_dot_loss = 100;
	double P_0 = 45;
	double D_s = 0.003;
	double e_v = 0.37;
	double dP_H_guess = P_h_in - P_h_out;
	double dP_C_guess = P_c_in - P_c_out;
	double L_guess = 1;
	double D_fr_guess = 1;
	double dP_max = P_h_in - P_h_out;
	setParameters("parallel", m_dot_h, m_dot_c, Q_dot_loss, P_0, D_s, e_v);
	setDesignTargets("cost", dP_max, Cost_target);
	//Fix
//	solveSystem();
	q_dot = Q_dot_a;
	T_c_out = T_C_out;
	T_h_out = T_H_out;

	ms_des_solved.m_DP_cold_des = dP_C;
	ms_des_solved.m_DP_hot_des = dP_H;
	ms_des_solved.m_eff_design = epsilon;
	//???
	ms_des_solved.m_min_DT_design = min((T_H_in - T_H_out), (T_C_out - T_C_in));
	//???
	ms_des_solved.m_NTU_design = NTU_R_e;
	ms_des_solved.m_Q_dot_design = Q_dot_a;
	ms_des_solved.m_T_c_out = T_C_out;
	ms_des_solved.m_T_h_out = T_H_out;
	ms_des_solved.m_UA_design_total = UA;
}

void RegenHX::off_design_solution(double T_c_in, double P_c_in, double m_dot_c, double P_c_out, double T_h_in, double P_h_in, double m_dot_h, double P_h_out, double & q_dot, double & T_c_out, double & T_h_out)
{
	T_C_in = T_c_in;
	T_H_in = T_h_in;
	P_C = P_c_in;
	P_H = P_h_in;
	dP_H = P_h_in - P_h_out;
	dP_C = P_c_in - P_c_out;
	m_dot_C = m_dot_c;
	m_dot_H = m_dot_h;
	CalculateThermoAndPhysicalModels();
	q_dot = Q_dot_calc;
	T_c_out = T_C_out;
	T_h_out = T_H_out;

	ms_od_solved.m_eff = epsilon;
	ms_od_solved.m_min_DT = min((T_H_in - T_H_out), (T_C_out - T_C_in));
	ms_od_solved.m_NTU = NTU_R_e;
	ms_od_solved.m_P_c_out = P_C_out;
	ms_od_solved.m_P_h_out = P_H_out;
	ms_od_solved.m_q_dot = Q_dot_calc;
	ms_od_solved.m_T_c_out = T_C_out;
	ms_od_solved.m_T_h_out = T_H_out;
	ms_od_solved.m_UA_total = UA;
}

void RegenHX::BalanceQdotAs() {
	//For setting iteration parameters refer to header file
	BalanceQdotAsHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.T_H_out_lowerBound, helper.T_H_out_upperBound, false);

	double T_H_out_solved, tol_solved;
	T_H_out_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;
	double Q_dot_a_difference_target = 0;

	BalanceQdotAs_Status = eq_solver.solve(helper.T_H_out_guess1, helper.T_H_out_guess2, Q_dot_a_difference_target, T_H_out_solved, tol_solved, iter_solved);
}

void RegenHX::BalancedPHs()
{
	//For setting iteration parameters refer to header file
	BalancedPHsHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.dP_H_lowerBound, helper.dP_H_upperBound, false);

	double dP_H_solved, tol_solved;
	dP_H_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;
	double dP_H_difference_target = 0;

	BalancedPHs_Status = eq_solver.solve(helper.dP_H_guess1, helper.dP_H_guess2, dP_H_difference_target, dP_H_solved, tol_solved, iter_solved);
}

void RegenHX::BalancedPCs()
{
	//For setting iteration parameters refer to header file
	BalancedPCsHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.dP_C_lowerBound, helper.dP_C_upperBound, false);

	double dP_C_solved, tol_solved;
	dP_C_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;
	double dP_C_difference_target = 0;

	BalancedPCs_Status = eq_solver.solve(helper.dP_C_guess1, helper.dP_C_guess2, dP_C_difference_target, dP_C_solved, tol_solved, iter_solved);
}

void RegenHX::figureOutD_fr()
{
	//For setting iteration parameters refer to header file
	FigureOutD_frHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.D_fr_lowerBound, helper.D_fr_upperBound, true);

	double D_fr_solved, tol_solved;
	D_fr_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

#ifdef _DEBUG
	DEBUG_LOG << "\tD_fr solver entered: @L = " << std::to_string(L) <<
		"; D_fr_guess1: " + std::to_string(helper.D_fr_guess1) << "; D_fr_guess2: " + std::to_string(helper.D_fr_guess2) << endl;
#endif
		
	figureOutD_fr_Status = eq_solver.solve(helper.D_fr_guess1, helper.D_fr_guess2, targetParameter, D_fr_solved, tol_solved, iter_solved);

#ifdef _DEBUG
	DEBUG_LOG << "\tD_fr solver exited with: " << solver_exit_modes[figureOutD_fr_Status] << endl;
	DEBUG_LOG.flush();
#endif
}

void RegenHX::figureOutL()
{
	//For setting iteration parameters refer to header file
	FigureOutLHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.L_lowerBound, helper.L_upperBound, true);

	double L_solved, tol_solved;
	L_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

	double y1, y2;

	double L_guess1 = helper.L_guess1;
	double L_guess2 = helper.L_guess2;

	


	
#ifdef _DEBUG
	DEBUG_LOG << "targetParameter = " + std::to_string(targetParameter) + "\nL solver entered wiht: L_guess1: " + std::to_string(L_guess1);
	DEBUG_LOG << "; L_guess2: " + std::to_string(L_guess2) << endl;
#endif
	
		if (eq_solver.test_member_function(L_guess1, &y1) == 0 && eq_solver.test_member_function(L_guess2, &y2) == 0) {
			C_monotonic_eq_solver::S_xy_pair pair1;
			C_monotonic_eq_solver::S_xy_pair pair2;

			pair1.x = L_guess1;
			pair1.y = y1;
			pair2.x = L_guess2;
			pair2.y = y2;

			figureOutL_Status = eq_solver.solve(pair1, pair2, targetdP_max, L_solved, tol_solved, iter_solved);
		}
		else {
			L_guess1 = helper.L_guess3;
			L_guess2 = helper.L_guess4;

			figureOutL_Status = eq_solver.solve(L_guess1, L_guess2, targetdP_max, L_solved, tol_solved, iter_solved);
		}
		

#ifdef _DEBUG
		LOG << "L solver exited with: " << solver_exit_modes[figureOutL_Status] << endl;
		LOG.flush();
#endif


	if (figureOutL_Status != C_monotonic_eq_solver::CONVERGED)
	{
		if (figureOutL_Status != C_monotonic_eq_solver::NO_SOLUTION && figureOutL_Status != C_monotonic_eq_solver::EQUAL_GUESS_VALUES) {
			LOG << "L solver exited with worse tolerance of " << tol_solved << endl;
			LOG.flush();
		}
		else {
			throw invalid_argument("L solver did not converge!");
		}
	}
}

void RegenHX::calculateWallThickness()
{
	numberOfCycles = 1.2 * operationYears * 365 * operationHoursPerDay * 60 * 60 / P_0;
	D_shell = D_fr + insulationThickness * 2;

	//TODO
	stressAmplitude = fatigueTable->getValue("Sa", "N", numberOfCycles) * 6.89475729; //Convert ksi to MPA
	R_i = D_shell / 2.0;

	WallThicknessLHelper helper(this);
	C_monotonic_eq_solver eq_solver(helper);

	int iter_limit = 50;
	eq_solver.settings(helper.tolerance, iter_limit, helper.th_lowerBound, helper.th_upperBound, true);

	double th_solved, tol_solved;
	th_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

	WallThickness_Status = eq_solver.solve(helper.th_guess1, helper.th_guess2, stressAmplitude, th_solved, tol_solved, iter_solved);

	if (WallThickness_Status != C_monotonic_eq_solver::CONVERGED)
	{
		throw invalid_argument("Wall thickness solver did not converge!");
	}
}

void RegenHX::calculateCost()
{

	volumeShellPerModule = (PI * (pow(R_o, 2) - pow(R_i, 2)) * L +
					4 / 3.0 * PI * (pow(R_o, 3) - pow(R_i, 3))); // a cylinder plus two half spheres

	volumeInsulationPerModule = (PI * (pow(R_i, 2) - pow(R_i - insulationThickness, 2)) * L +
		4 / 3.0 * PI * pow(R_i, 3)*insulationParameter);

	volumeBedPerModule = V_0 * (1 - e_v);

	volumeShellTotal = volumeShellPerModule * numberOfModulesTotal;
	volumeInsulationTotal = volumeInsulationPerModule * numberOfModulesTotal;
	volumeBedTotal = volumeBedPerModule * numberOfModulesTotal;

	double shellMaterialDensity = shellMaterialTable->getValue("rho", "T", T_f);
	double insulationMaterialDensity = insulationMaterialTable->getValue("rho", "T", T_f);
	double bedMaterialDensity = bedMaterialTable->getValue("rho", "T", T_f);

	massShellPerModule = shellMaterialDensity * volumeShellPerModule;
	massInsulationPerModule = insulationMaterialDensity * volumeInsulationPerModule;
	massBedPerModule = bedMaterialDensity * volumeBedPerModule;

	massShellTotal = massShellPerModule * numberOfModulesTotal;
	massInsulationTotal = massInsulationPerModule * numberOfModulesTotal;
	massBedTotal = massBedPerModule * numberOfModulesTotal;

	costShellMaterialPerModule = specificCostShellMaterial * massShellPerModule;
	costInsulationMaterialPerModule = specificCostInsulationMaterial * massInsulationPerModule;
	costBedMaterialPerModule = specificCostBedMaterial * massBedPerModule;

	costShellMaterialTotal = specificCostShellMaterial * massShellTotal;
	costInsulationMaterialTotal = specificCostInsulationMaterial * massInsulationTotal;
	costBedMaterialTotal = specificCostBedMaterial * massBedTotal;

	costPerModuleMaterial = costShellMaterialPerModule + costInsulationMaterialPerModule + costBedMaterialPerModule;
	costPerModuleTotal = costPerModuleMaterial + priceCastingPerModule + priceCastingSteelPerModule + priceWeldingPerModule;

	costHXMaterial = costShellMaterialTotal + costInsulationMaterialTotal + costBedMaterialTotal;
	costHXTotal = costHXMaterial + (priceCastingPerModule + priceCastingSteelPerModule + priceWeldingPerModule) 
																		* numberOfModulesTotal;
}

int RegenHX::BalanceQdotAsHelper::operator()(double T_H_out, double * QdotAsDifference) {

	system->T_H_out = T_H_out;

	try {
		system->CalculateThermoAndPhysicalModels();
	}catch (const invalid_argument& e) {
#ifdef _DEBUG
		system->lastCoreError = e.what();
#endif
		return -1;
	}
	
	*QdotAsDifference = system->getQdotADifference();
	return 0;
}

int RegenHX::BalancedPHsHelper::operator()(double dP_H, double * dP_HsDifference)
{
	system->dP_H = dP_H;

	system->BalanceQdotAs();
	
	if (system->BalanceQdotAs_Status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_HsDifference = system->getdP_HsDifference();
	return 0;
}

int RegenHX::BalancedPCsHelper::operator()(double dP_C, double * dP_CsDifference)
{
	system->dP_C = dP_C;
	
	system->BalancedPHs();

	if (system->BalancedPHs_Status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_CsDifference = system->getdP_CsDifference();
	return 0;
}

int RegenHX::FigureOutD_frHelper::operator()(double D_fr, double * targetParameter)
{

#ifdef _DEBUG
	system->LOG << "\t\tL: " << system->L << ", D_fr: " << D_fr << " ";
#endif

	system->D_fr = D_fr;
	

	system->BalancedPCs();
#ifdef _DEBUG
	system->LOG << "epsilon -> " << system->epsilon;
	system->LOG << ", ua -> " << system->UA;
	system->LOG << ", dP_C -> " << system->dP_C << ", dP_C_calc -> " << system->dP_C_calc << ", dP_H -> " << system->dP_H << ", dP_H_calc -> " << system->dP_H_calc << endl;
	system->LOG << "\t\t\tdPC_Status -> " << system->solver_exit_modes[system->BalancedPCs_Status];
	system->LOG << ", dPH_Status -> " << system->solver_exit_modes[system->BalancedPHs_Status];
	system->LOG << ", Q_dot_Status -> " << system->solver_exit_modes[system->BalanceQdotAs_Status] << endl;
#endif


	if (system->BalancedPCs_Status != C_monotonic_eq_solver::CONVERGED) {
#ifdef _DEBUG
		system->LOG << "\t*Core Error " << system->lastCoreError << "*" << endl;
#endif
		return -1;
	}

	if (system->targetMode == "cost") {
		system->calculateWallThickness();
		system->calculateCost();
#ifdef _DEBUG
		system->LOG << "\t\tcost -> " << system->costHXTotal;
		system->LOG << ", thickness -> " << system->wallThickness;
		system->LOG << ", WT_Status -> " << system->solver_exit_modes[system->WallThickness_Status] << endl;
#endif
		*targetParameter = system->costHXTotal;
	}
	else if (system->targetMode == "epsilon") {
		*targetParameter = system->epsilon;
	}
	else if (system->targetMode == "ua") {
		*targetParameter = system->UA;
	}
	
	return 0;
}

int RegenHX::FigureOutLHelper::operator()(double L, double * dP_max)
{
	 system->L = L;

	 system->figureOutD_fr();

	 if (system->figureOutD_fr_Status != C_monotonic_eq_solver::CONVERGED) {
		 return -1;
	 }

	*dP_max = system->dP_max;

	return 0;
}

int RegenHX::WallThicknessLHelper::operator()(double th, double * stressAmplitude)
{
	system->wallThickness = th;

	system->R_o = system->R_i + system->wallThickness;

	double RiSquared = pow(system->R_i, 2);

	double RoSquared = pow(system->R_o, 2);

	double RSquareDifference = RoSquared - RiSquared;

	double RSquareProduct = RoSquared * RiSquared;

	double magicPieceRi = RiSquared * RSquareDifference;

	double magicPieceRo = RoSquared * RSquareDifference;

	double magicPieceHigh = RSquareProduct * (system->Patm - system->P_high);

	double magicPieceLow = RSquareProduct * (system->Patm - system->P_low);

	//sigma_a_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);		
	double sigma_a_h = (system->P_high * RiSquared - system->Patm * RoSquared) / RSquareDifference; //axial stress calculation

	//sigma_c_1_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at inner surface
	double sigma_c_1_h = sigma_a_h - magicPieceHigh / magicPieceRi;	//hoop stress at inner surface

	//sigma_c_2_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at outer surface
	double sigma_c_2_h = sigma_a_h - magicPieceHigh / magicPieceRo;	//hoop stress at outer surface
	
	//sigma_r_1_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at inner surface
	double sigma_r_1_h = sigma_a_h + magicPieceHigh / magicPieceRi;	//radial stress at inner surface
	
	//sigma_r_2_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_high) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at outer surface
	double sigma_r_2_h = sigma_a_h + magicPieceHigh / magicPieceRo;	//radial stress at outer surface
 
	double sigma_prime_1_h = sqrt((pow((sigma_r_1_h - sigma_c_1_h), 2) + pow((sigma_c_1_h - sigma_a_h), 2) + pow((sigma_a_h - sigma_r_1_h), 2)) / 2.0);	//equivalent stress at inner surface
	double sigma_prime_2_h = sqrt((pow((sigma_r_2_h - sigma_c_2_h), 2) + pow((sigma_c_2_h - sigma_a_h), 2) + pow((sigma_a_h - sigma_r_2_h), 2)) / 2.0);	//equivalent stress at outer surface
	double sigma_prime_h = max(sigma_prime_1_h, sigma_prime_2_h);		//largest equivalent stress


	// sigma_a_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);	//axial stress calculation
	double sigma_a_l = (system->P_low*RiSquared - system->Patm * RoSquared) / RSquareDifference;	//axial stress calculation

	//sigma_c_1_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at inner surface
	double sigma_c_1_l = sigma_a_l - magicPieceLow / magicPieceRi;	//hoop stress at inner surface
	
	//sigma_c_2_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) - r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//hoop stress at outer surface
	double sigma_c_2_l = sigma_a_l - magicPieceLow / magicPieceRo;	//hoop stress at outer surface
	
	//sigma_r_1_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_i ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at inner surface
	double sigma_r_1_l = sigma_a_l + magicPieceLow / magicPieceRi;	//radial stress at inner surface
	
	//sigma_r_2_l = (P_low*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2) + r_i ^ 2 * r_o ^ 2 * (P_o - P_low) / (r_o ^ 2 * (r_o ^ 2 - r_i ^ 2));	//radial stress at outer surface
	double sigma_r_2_l = sigma_a_l + magicPieceLow / magicPieceRo;	//radial stress at outer surface
 
	double sigma_prime_1_l = sqrt((pow((sigma_r_1_l - sigma_c_1_l), 2) + pow((sigma_c_1_l - sigma_a_l), 2) + pow((sigma_a_l - sigma_r_1_l), 2)) / 2.0);	//equivalent stress at inner surface
	double sigma_prime_2_l = sqrt((pow((sigma_r_2_l - sigma_c_2_l), 2) + pow((sigma_c_2_l - sigma_a_l), 2) + pow((sigma_a_l - sigma_r_2_l), 2)) / 2.0);	//equivalent stress at outer surface
	double sigma_prime_l = max(sigma_prime_1_l, sigma_prime_2_l);		//largest equivalent stress
	
	*stressAmplitude = (sigma_prime_h - sigma_prime_l) / 2.0;


	return 0;
}
