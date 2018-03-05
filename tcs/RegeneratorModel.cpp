#include "RegeneratorModel.h"

std::string SSCDIR1(std::getenv("SSCDIR"));

string const RegeneratorModel::LOG_FILEPATH = SSCDIR1 + "/build_sdk/examples/RegenHX_LOG.log";
string const RegeneratorModel::PROPERTY_FILES = SSCDIR1 + "/tcs/PropertyFiles/";
string const RegeneratorModel::SPHERES_RP_TABLE_PATH = PROPERTY_FILES + "Spheres_RP.csv";
string const RegeneratorModel::BALANCED_REGENERATOR_TABLE_PATH = PROPERTY_FILES + "balanced-regenerator.csv";
string const RegeneratorModel::FATIGUE_TABLE_PATH = PROPERTY_FILES + "fatigue.csv";


RegeneratorModel::RegeneratorModel()
{
	loadTables();

	if (spdlog::get("logger") == nullptr) {
		auto logger = spdlog::basic_logger_mt("logger", LOG_FILEPATH);
		spdlog::flush_on(spdlog::level::info);
		spdlog::drop("logger");
		spdlog::register_logger(logger);
		//spdlog::set_pattern("[%b %d %H:%M:%S] [%l] %v");
		spdlog::set_pattern("%v");
	}
}

RegeneratorModel::~RegeneratorModel()
{
	spdlog::get("logger")->flush();
	spdlog::drop("logger");
}

void RegeneratorModel::packedspheresNdFit(double Re, double * f, double * j_H)
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

double RegeneratorModel::hx(double NTU, double C_1, double C_2)
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

void RegeneratorModel::getpropsregenFitCO2(double T, double P, double * rho, double * mu, double * k, double * Pr, double * Cp)
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

void RegeneratorModel::packedspheresFitCO2(double m_dot, double d, double A_fr, double L, double T, double P, double porosity, 
																			double * f, double * h, double * NTU, double * DP)
{
	double rho, mu, k, Pr, Cp;
	try {
		getpropsregenFitCO2(T, P, &rho, &mu, &k, &Pr, &Cp);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	double r_h = porosity * d / (6.0*(1 - porosity));	//[m]
	double G = m_dot / (porosity*A_fr);	//[kg/m^2-s]
	double Re = G * (4 * r_h) / mu;	//[-]

	double j_H;	//[-]
	try {
		packedspheresNdFit(Re, f, &j_H);
	}
	catch (const invalid_argument& e) {
		throw e;
	}


	double St = j_H / pow(Pr, (2.0 / 3.0));	//[-]
	*NTU = St * L / r_h;	//[-]
	*h = j_H * G*Cp / pow(Pr, (2.0 / 3.0)) * 1000.0;	//Convert kW to W. [W]


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
}

void RegeneratorModel::massflowVariablesInit()
{
	C_dot_H = C_p_H * m_dot_H;
	C_dot_C = C_p_C * m_dot_C;

	C_dot_min = min(C_dot_H, C_dot_C);
	C_dot_max = max(C_dot_H, C_dot_C);

	C_R = C_dot_min / C_dot_max;
}

void RegeneratorModel::calculateModel()
{
	//							[SECTION 1]
	//Pressure at hot exit
	P_H_out = P_H - dP_H;	//[kPa]
							//Pressure at cold exit
	P_C_out = P_C - dP_C;	//[kPa]

	if (P_C_out <= N_co2_props::P_lower_limit || P_H_out <= N_co2_props::P_lower_limit || P_C_out >= N_co2_props::P_upper_limit || P_H_out >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Outlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}

	CO2_TP(T_C_in, P_H_out, &CO2State);
	h_H_out_max = CO2State.enth;

	CO2_TP(T_H_in, P_C_out, &CO2State);
	h_C_out_max = CO2State.enth;

	//Maximum heat transfer if hot stream exiting at cold temp
	Q_dot_max_H = m_dot_H * (h_H_in - h_H_out_max);

	//Maximum heat transfer if cold stream exiting at hot temp
	Q_dot_max_C = m_dot_C * (h_C_out_max - h_C_in);

	//Actual maximum possible heat transfer
	Q_dot_max = min(Q_dot_max_H, Q_dot_max_C);

	if (T_H_out <= N_co2_props::T_lower_limit || T_H_out >= N_co2_props::T_upper_limit) {
		throw invalid_argument("T_H_out is either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	CO2_TP(T_H_out, P_H_out, &CO2State);
	h_H_out = CO2State.enth;

	Q_dot_a = m_dot_H * (h_H_in - h_H_out);	

	h_C_out = h_C_in + Q_dot_a / m_dot_C;
	CO2_PH(P_C_out, h_C_out, &CO2State);

	T_C_out = CO2State.temp;

	Q_dot = Q_dot_a + Q_dot_loss;

	//							[SECTION 2]

	A_fr = PI * pow(D_fr, 2) / 4;
	T_H_f = (T_H_in + T_H_out) / 2;
	T_C_f = (T_C_in + T_C_out) / 2;

	try {
		packedspheresFitCO2(m_dot_H, D_s, A_fr, L, T_H_f, P_H, e_v, &f_H, &h_H, &NTU_H, &dP_H_calc);
		packedspheresFitCO2(m_dot_C, D_s, A_fr, L, T_C_f, P_C, e_v, &f_C, &h_C, &NTU_C, &dP_C_calc);
	}
	catch (const invalid_argument& e) {
		throw e;
	}

	V_0 = A_fr * L;

	m_s = V_0 * (1 - e_v) * rho_s;

	A_s = 6.0 * (1.0 - e_v)*V_0 / D_s;

	UA = 2 * A_s * h_H * h_C / (h_H + h_C) / 1000.0;

	NTU_R = UA / C_dot_min;

	C_m = 2.0 * m_s * c_s / P_0 / C_dot_min;

	NTU_R_e = 2.0 * C_R * NTU_R / (1.0 + C_R);

	C_m_e = 2.0 * C_R*C_m / (1.0 + C_R);

	try {
		if (C_m_e < 0.001) {
			throw invalid_argument("No value provided for U > 1000!");
		}

		epsilon_1 = regeneratorTable->getValue(NTU_R_e, 1 / C_m_e);
	}
	catch (exception e) {
		throw e;
	}

	X = (1.0 - pow(C_R, 2)) / (2.0 * C_R)*(epsilon_1 / (1.0 - epsilon_1));

	epsilon = (1.0 - exp(-X)) / (1.0 - C_R * exp(-X));

	if (epsilon < 0 || epsilon > 1) {
		throw invalid_argument("Epsilon is not between 0 and 1!");
	}

	Q_dot_calc = epsilon * Q_dot_max;
	Q_dot_a_calc = Q_dot_calc - Q_dot_loss;
	dP_max = max(dP_C, dP_H);
}

void RegeneratorModel::loadTables() {
	bedMaterialTable = new LookupTable_1D(PROPERTY_FILES + bedMaterialName + ".csv");
	shellMaterialTable = new LookupTable_1D(PROPERTY_FILES + shellMaterialName + ".csv");
	insulationMaterialTable = new LookupTable_1D(PROPERTY_FILES + insulationMaterialName + ".csv");
	regeneratorTable = new LookupTable_2D(BALANCED_REGENERATOR_TABLE_PATH);
	spheresRPTable = new LookupTable_1D(SPHERES_RP_TABLE_PATH);
	fatigueTable = new LookupTable_1D(FATIGUE_TABLE_PATH);
}

void RegeneratorModel::calculateCost()
{
	D_shell = D_fr + insulationThickness * 2;
	R_i = D_shell / 2.0;
	
	WallThickness->solve();

	volumeShell = 2 * (PI * (pow(R_o, 2) - pow(R_i, 2)) * L + 4 / 3.0 * PI * (pow(R_o, 3) - pow(R_i, 3))); // a cylinder plus two half spheres
	volumeInsulation = 2* (PI * (pow(R_i, 2) - pow(R_i - insulationThickness, 2)) * L + 4 / 3.0 * PI * pow(R_i, 3)*insulationParameter);
	volumeBed = 2* V_0 * (1 - e_v);

	double shellMaterialDensity = shellMaterialTable->getValue("rho", "T", T_f);
	double insulationMaterialDensity = insulationMaterialTable->getValue("rho", "T", T_f);
	double bedMaterialDensity = bedMaterialTable->getValue("rho", "T", T_f);

	massShell = shellMaterialDensity * volumeShell;
	massInsulation = insulationMaterialDensity * volumeInsulation;
	massBed = bedMaterialDensity * volumeBed;

	costShellMaterial = specificCostShellMaterial * massShell;
	costInsulationMaterial = specificCostInsulationMaterial * massInsulation;
	costBedMaterial = specificCostBedMaterial * massBed;

	costMaterial = costShellMaterial + costInsulationMaterial + costBedMaterial;
	costModule = costMaterial + 2 * (priceCasting + priceCastingSteel + priceWelding);
}

void RegeneratorModel::setInletStates(double T_H_in, double P_H, double m_dot_H, double T_C_in, double P_C, double m_dot_C)
{
	this->T_H_in = T_H_in;
	this->P_H = P_H;
	this->m_dot_H = m_dot_H;
	this->T_C_in = T_C_in;
	this->P_C = P_C;
	this->m_dot_C = m_dot_C;

	if (P_C <= N_co2_props::P_lower_limit || P_H <= N_co2_props::P_lower_limit || P_C >= N_co2_props::P_upper_limit || P_H >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Inlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
	}
	if (T_H_in <= N_co2_props::T_lower_limit || T_C_in <= N_co2_props::T_lower_limit || T_H_in >= N_co2_props::T_upper_limit || T_C_in >= N_co2_props::T_upper_limit) {
		throw invalid_argument("Inlet temperatures are either below " + std::to_string(N_co2_props::T_lower_limit) + "[K] or above " + std::to_string(N_co2_props::T_upper_limit) + "[K]!");
	}

	//Enthalpy at hot inlet
	CO2_TP(T_H_in, P_H, &CO2State);
	h_H_in = CO2State.enth;

	//Enthalpy at cold inlet
	CO2_TP(T_C_in, P_C, &CO2State);
	h_C_in = CO2State.enth;

	CO2_TP(T_C_in, P_H, &CO2State);
	h_H_out_max_p = CO2State.enth;

	CO2_TP(T_H_in, P_C, &CO2State);
	h_C_out_max_p = CO2State.enth;

	C_p_H = (h_H_in - h_H_out_max_p) / (T_H_in - T_C_in);
	C_p_C = (h_C_out_max_p - h_C_in) / (T_H_in - T_C_in);

	T_f = (T_H_in + T_C_in) / 2;
	rho_s = bedMaterialTable->getValue("rho", "T", T_f);
	c_s = bedMaterialTable->getValue("c", "T", T_f);

	massflowVariablesInit();
}

void RegeneratorModel::setParameters(double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->Q_dot_loss = Q_dot_loss;
	this->P_0 = P_0;
	this->D_s = D_s;
	this->e_v = e_v;

	numberOfCycles = 1.2 * operationYears * 365 * operationHoursPerDay * 60 * 60 / P_0;
	stressAmplitude = fatigueTable->getValue("Sa", "N", numberOfCycles) * 6.89475729; //Convert ksi to MPA
}

void RegeneratorModel::setDesignTargets(targetModes::targetModes targetMode, double targetParameter, double dP_max)
{
	this->targetdP_max = dP_max;
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;
}

double RegeneratorModel::densityIntegral(double T_low, double T_high, double P)
{
	double m = (T_high - T_low) / L;
	double dL = 0.001;
	double T_int = T_low - m * dL;
	double integral = 0;

	for (int i = 0; i <= ((int) (L / dL)); i++) {
		T_int += m * dL;
		CO2_TP(T_int, P, &CO2State);
		integral += CO2State.dens;
	}

	return integral * dL;
}

void RegeneratorModel::calcCarryoverMassFlow()
{
	CO2_TP(T_C_out, P_C, &CO2State);
	double rho_H_H_extra = CO2State.dens;

	CO2_TP(T_C_in, P_C, &CO2State);
	double rho_H_C_extra = CO2State.dens;

	CO2_TP(T_H_in, P_H, &CO2State);
	double rho_L_H_extra = CO2State.dens;

	CO2_TP(T_H_out, P_H, &CO2State);
	double rho_L_C_extra = CO2State.dens;

	double mass_H = vol_extra / 4 * (rho_H_H_extra - rho_L_H_extra);
	double mass_C = vol_extra / 4 * (rho_H_C_extra - rho_L_C_extra);

	double integral_H = densityIntegral(T_C_in, T_C_out, P_C);
	double integral_L = densityIntegral(T_H_out, T_H_in, P_H);

	double mass = 2 * (e_v * pow(D_fr, 2) / 4 * PI *(integral_H - integral_L) + (mass_C + mass_H)) / CO;
	m_dot_carryover = mass / P_0;
}

void RegeneratorModel::carryoverEnthDrop()
{
	h_H_out = (h_C_in * m_dot_carryover + h_H_out * (m_dot_H - m_dot_carryover)) / m_dot_H;
	CO2_PH(P_H_out, h_H_out, &CO2State);
	T_H_out = CO2State.temp;
}

int RegeneratorModel::HeatTransfer_ME(double T_H_out, double * QdotAsDifference) {

	this->T_H_out = T_H_out;

	if (T_H_out > T_H_in) {
		return -1;
	}

	try {
		calculateModel();
	}
	catch (const invalid_argument& e) {
		return -1;
	}

	*QdotAsDifference = Q_dot_a_calc - Q_dot_a;
	return 0;
}

int RegeneratorModel::HotPressureDrop_ME(double dP_H, double * dP_HsDifference)
{
	this->dP_H = dP_H;

	if (HeatTransfer->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_HsDifference = dP_H_calc - this->dP_H;
	return 0;
}

int RegeneratorModel::ColdPressureDrop_ME(double dP_C, double * dP_CsDifference)
{
	this->dP_C = dP_C;

	if (HotPressureDrop->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_CsDifference = dP_C_calc - this->dP_C;
	return 0;
}

int RegeneratorModel::Length_ME(double L, double * dP_max)
{
	this->L = L;

	int status = ColdPressureDrop->solve();
	

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_max = this->dP_max;

	return 0;
}

int RegeneratorModel::Diameter_ME(double D_fr, double * targetParameter)
{
	this->D_fr = D_fr;

	int status = Length->solve();

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	if (targetMode == targetModes::COST) {
		calculateCost();
		*targetParameter = costModule;
	}
	else if (targetMode == targetModes::EPSILON) {
		*targetParameter = epsilon;

	}
	else if (targetMode == targetModes::UA) {
		*targetParameter = UA;
	}

	return 0;
}

int RegeneratorModel::CarryoverMassFlow_ME(double m_dot_carryover, double *comass_difference)
{
	this->m_dot_C -= m_dot_carryover;
	this->m_dot_H -= m_dot_carryover;

	massflowVariablesInit();
	
	double status = Diameter->solve();

	this->m_dot_C += m_dot_carryover;
	this->m_dot_H += m_dot_carryover;
	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	calcCarryoverMassFlow();

	*comass_difference = m_dot_carryover - this->m_dot_carryover;
	return 0;
}

int RegeneratorModel::WallThickness_ME(double th, double * stressAmplitude)
{
	wallThickness = th;

	R_o = R_i + wallThickness;

	double RiSquared = pow(R_i, 2);

	double RoSquared = pow(R_o, 2);

	double RSquareDifference = RoSquared - RiSquared;

	double RSquareProduct = RoSquared * RiSquared;

	double magicPieceRi = RiSquared * RSquareDifference;

	double magicPieceRo = RoSquared * RSquareDifference;

	double magicPieceHigh = RSquareProduct * (Patm - P_C);

	double magicPieceLow = RSquareProduct * (Patm - P_H_out);

	//sigma_a_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);		
	double sigma_a_h = (P_C * RiSquared - Patm * RoSquared) / RSquareDifference; //axial stress calculation

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
	double sigma_a_l = (P_H_out*RiSquared - Patm * RoSquared) / RSquareDifference;	//axial stress calculation

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

int RegeneratorModel::solveSystem()
{
	HeatTransfer_SP.solverName = "Balance Heat Tansfer";
	HeatTransfer_SP.target = 0;
	HeatTransfer_SP.guessValue1 = T_C_in;						HeatTransfer_SP.guessValue2 = (T_C_in + T_H_in) / 2;
	HeatTransfer_SP.lowerBound = N_co2_props::T_lower_limit;	HeatTransfer_SP.upperBound = N_co2_props::T_upper_limit;
	HeatTransfer_SP.tolerance = 0.001;
	HeatTransfer_SP.iterationLimit = 50;
	HeatTransfer_SP.isErrorRel = false;
	HeatTransfer_SP.classInst = this;
	HeatTransfer_SP.monoEquation = &RegeneratorModel::HeatTransfer_ME;

	HotPressureDrop_SP.solverName = "Balance Hot Pressure Drop";
	HotPressureDrop_SP.target = 0;
	HotPressureDrop_SP.guessValue1 = targetdP_max - 10;		HotPressureDrop_SP.guessValue2 = targetdP_max;
	HotPressureDrop_SP.lowerBound = 0.1;					HotPressureDrop_SP.upperBound = P_H;
	HotPressureDrop_SP.tolerance = 0.001;
	HotPressureDrop_SP.iterationLimit = 50;
	HotPressureDrop_SP.isErrorRel = false;
	HotPressureDrop_SP.classInst = this;
	HotPressureDrop_SP.monoEquation = &RegeneratorModel::HotPressureDrop_ME;

	ColdPressureDrop_SP.solverName = "Balance Cold Pressure Drop";
	ColdPressureDrop_SP.target = 0;
	ColdPressureDrop_SP.guessValue1 = targetdP_max / 100 - 10;		ColdPressureDrop_SP.guessValue2 = targetdP_max / 100;
	ColdPressureDrop_SP.lowerBound = 0.1;							ColdPressureDrop_SP.upperBound = P_C;
	ColdPressureDrop_SP.tolerance = 0.001;
	ColdPressureDrop_SP.iterationLimit = 50;
	ColdPressureDrop_SP.isErrorRel = false;
	ColdPressureDrop_SP.classInst = this;
	ColdPressureDrop_SP.monoEquation = &RegeneratorModel::ColdPressureDrop_ME;
	
	Length_SP.solverName = "Length Solver";
	Length_SP.target = targetdP_max;
	Length_SP.guessValue1 = 0.8;		Length_SP.guessValue2 = 1;
	Length_SP.lowerBound = 0.1;			Length_SP.upperBound = 10;
	Length_SP.tolerance = 0.001;
	Length_SP.iterationLimit = 50;
	Length_SP.isErrorRel = true;
	Length_SP.classInst = this;
	Length_SP.monoEquation = &RegeneratorModel::Length_ME;

	Diameter_SP.solverName = "Diameter Solver";
	Diameter_SP.target = targetParameter;
	Diameter_SP.guessValue1 = 0.7;		Diameter_SP.guessValue2 = 0.9;
	Diameter_SP.lowerBound = 0.1;		Diameter_SP.upperBound = 10;
	Diameter_SP.tolerance = 0.001;
	Diameter_SP.iterationLimit = 50;
	Diameter_SP.isErrorRel = true;
	Diameter_SP.classInst = this;
	Diameter_SP.monoEquation = &RegeneratorModel::Diameter_ME;

	CarryoverMassFlow_SP.solverName = "Carryover Mass Solver";
	CarryoverMassFlow_SP.target = 0;
	CarryoverMassFlow_SP.lowerBound = 0;		CarryoverMassFlow_SP.upperBound = min(m_dot_C, m_dot_H);
	CarryoverMassFlow_SP.tolerance = 0.001;
	CarryoverMassFlow_SP.iterationLimit = 50;
	CarryoverMassFlow_SP.isErrorRel = false;
	CarryoverMassFlow_SP.classInst = this;
	CarryoverMassFlow_SP.monoEquation = &RegeneratorModel::CarryoverMassFlow_ME;

	WallThickness_SP.solverName = "Wall Thickness";
	WallThickness_SP.target = stressAmplitude;
	WallThickness_SP.guessValue1 = 0.04;		WallThickness_SP.guessValue2 = 0.05;
	WallThickness_SP.lowerBound = 0;			WallThickness_SP.upperBound = 1;
	WallThickness_SP.tolerance = 0.001;
	WallThickness_SP.iterationLimit = 50;
	WallThickness_SP.isErrorRel = true;
	WallThickness_SP.classInst = this;
	WallThickness_SP.monoEquation = &RegeneratorModel::WallThickness_ME;

	HeatTransfer = new MonoSolver<RegeneratorModel>(&HeatTransfer_SP);
	HotPressureDrop = new MonoSolver<RegeneratorModel>(&HotPressureDrop_SP);
	ColdPressureDrop = new MonoSolver<RegeneratorModel>(&ColdPressureDrop_SP);
	Length = new MonoSolver<RegeneratorModel>(&Length_SP);
	Diameter = new MonoSolver<RegeneratorModel>(&Diameter_SP);
	WallThickness = new MonoSolver<RegeneratorModel>(&WallThickness_SP);
	
	int statusSolver = Diameter->solve();
	if (statusSolver != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	calcCarryoverMassFlow();

	CarryoverMassFlow_SP.guessValue1 = m_dot_carryover;
	CarryoverMassFlow_SP.guessValue2 = m_dot_carryover + 1;
	CarryoverMassFlow = new MonoSolver<RegeneratorModel>(&CarryoverMassFlow_SP);
	HeatTransfer->updateGuesses(T_H_out - 1, T_H_out);
	HotPressureDrop->updateGuesses(dP_H - 10, dP_H);
	ColdPressureDrop->updateGuesses(dP_C - 10, dP_C);
	Diameter->updateGuesses(D_fr - 0.001, D_fr);
	Length->updateGuesses(L - 0.001, L);

	statusSolver = CarryoverMassFlow->solve();

	if (statusSolver == C_monotonic_eq_solver::CONVERGED) {
		carryoverEnthDrop();
		
		m_dot_H -= m_dot_carryover;
		m_dot_C -= m_dot_carryover;

		calculateCost();

		return 0;
	}

	return -1;
}

void RegeneratorModel::getSolution(RegeneratorSolution * solution)
{
	solution->dP_C = dP_C;
	solution->dP_H = dP_H;
	solution->epsilon = epsilon;
	solution->T_H_in = T_H_in;
	solution->T_C_in = T_C_in;
	solution->T_C_out = T_C_out;
	solution->T_H_out = T_H_out;
	solution->NTU_R_e = NTU_R_e;
	solution->Q_dot_a = Q_dot_a;
	solution->UA = UA;
	solution->m_dot_H = m_dot_H;
	solution->m_dot_C = m_dot_C;
	solution->m_dot_carryover = m_dot_carryover;
	solution->costModule = costModule;
	solution->L = L;
	solution->D_fr = D_fr;
	solution->wallThickness = wallThickness;
}
