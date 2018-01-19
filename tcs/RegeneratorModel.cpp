#include "RegeneratorModel.h"

string const RegeneratorModel::LOG_FILEPATH = "./RegenHX_LOG.log";
string const RegeneratorModel::PROPERTY_FILES = "../tcs/PropertyFiles/";
string const RegeneratorModel::SPHERES_RP_TABLE_PATH = PROPERTY_FILES + "Spheres_RP.csv";
string const RegeneratorModel::BALANCED_REGENERATOR_TABLE_PATH = PROPERTY_FILES + "balanced-regenerator.csv";
string const RegeneratorModel::FATIGUE_TABLE_PATH = PROPERTY_FILES + "fatigue.csv";

RegeneratorModel::RegeneratorModel()
{
	loadTables();

	auto logger = spdlog::basic_logger_mt("logger", LOG_FILEPATH);

#ifdef _DEBUG
	spdlog::set_level(spdlog::level::debug);
#endif
	
	spdlog::drop("logger");
	spdlog::register_logger(logger);
	spdlog::set_pattern("[%b %d %H:%M:%S] [%l] %v");
}

RegeneratorModel::~RegeneratorModel()
{
	spdlog::get("logger")->flush();
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

void RegeneratorModel::CalculateThermoAndPhysicalModels()
{
	//							[SECTION 1]
	//Pressure at hot exit
	P_H_out = P_H - dP_H;	//[kPa]
							//Pressure at cold exit
	P_C_out = P_C - dP_C;	//[kPa]

	if (P_C_out <= N_co2_props::P_lower_limit || P_H_out <= N_co2_props::P_lower_limit || P_C_out >= N_co2_props::P_upper_limit || P_H_out >= N_co2_props::P_upper_limit) {
		throw invalid_argument("Outlet pressures are either below " + std::to_string(N_co2_props::P_lower_limit) + "[kPa] or above " + std::to_string(N_co2_props::P_upper_limit) + "[kPa]!");
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

	T_H_f = (T_H_in + T_H_out) / 2;
	T_C_f = (T_C_in + T_C_out) / 2;

	CO2_TP(T_C_in, P_H, &CO2State);
	h_H_out_max_p = CO2State.enth;

	CO2_TP(T_H_in, P_C, &CO2State);
	h_C_out_max_p = CO2State.enth;

	C_p_H = (h_H_in - h_H_out_max_p) / (T_H_in - T_C_in);
	C_p_C = (h_C_out_max_p - h_C_in) / (T_H_in - T_C_in);

	C_dot_H = C_p_H * m_dot_H;
	C_dot_C = C_p_C * m_dot_C;

	C_dot_min = min(C_dot_H, C_dot_C);
	C_dot_max = max(C_dot_H, C_dot_C);

	C_R = C_dot_min / C_dot_max;

	//							[SECTION 3]

	A_fr = PI * pow((D_fr) / 2.0, 2);
	T_H_f = (T_H_in + T_H_out) / 2.0;
	T_C_f = (T_C_in + T_C_out) / 2.0;

	try {
		packedspheresFitCO2(m_dot_H, D_s, A_fr, L, T_H_f, P_H, e_v, &f_H, &h_H, &NTU_H, &dP_H_calc);
		packedspheresFitCO2(m_dot_C, D_s, A_fr, L, T_C_f, P_C, e_v, &f_C, &h_C, &NTU_C, &dP_C_calc);
	}
	catch (const invalid_argument& e) {
		throw e;
	}


	T_f = (T_H_in + T_C_in) / 2.0;
	rho_s = bedMaterialTable->getValue("rho", "T", T_f);
	c_s = bedMaterialTable->getValue("c", "T", T_f);

	V_0 = A_fr * L;

	m_s = V_0 * (1 - e_v) * rho_s;

	A_s = 6.0 * (1.0 - e_v)*V_0 / D_s;

	NTU_R = 1.0 / (C_dot_min*1000.0)*(1.0 / ((1.0 / (h_H*A_s)) + (1.0 / (h_C*A_s))));

	UA = 1.0 / ((1.0 / (h_H*A_s)) + (1.0 / (h_C*A_s))) / 1000.0;

	f_0 = 1.0 / P_0;

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
	
	solveForWallThickness->solve();

	volumeShell = (PI * (pow(R_o, 2) - pow(R_i, 2)) * L + 4 / 3.0 * PI * (pow(R_o, 3) - pow(R_i, 3))); // a cylinder plus two half spheres
	volumeInsulation = (PI * (pow(R_i, 2) - pow(R_i - insulationThickness, 2)) * L + 4 / 3.0 * PI * pow(R_i, 3)*insulationParameter);
	volumeBed = V_0 * (1 - e_v);

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
	costModule = costMaterial + priceCasting + priceCastingSteel + priceWelding;
}

void RegeneratorModel::setInletStates(double T_H_in, double P_H, double T_C_in, double P_C)
{
	this->T_H_in = T_H_in;
	this->P_H = P_H;
	this->T_C_in = T_C_in;
	this->P_C = P_C;

	//LOG << "Inlet states are set. T_H_in = " << this->T_H_in << ", P_H = " << this->P_H << ", T_C_in = " << this->T_C_in << ", P_C = " << this->P_C << endl;
}

void RegeneratorModel::setParameters(double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->m_dot_H = m_dot_H;
	this->m_dot_C = m_dot_C;
	this->Q_dot_loss = Q_dot_loss;
	this->P_0 = P_0;
	this->D_s = D_s;
	this->e_v = e_v;

	numberOfCycles = 1.2 * operationYears * 365 * operationHoursPerDay * 60 * 60 / P_0;
	stressAmplitude = fatigueTable->getValue("Sa", "N", numberOfCycles) * 6.89475729; //Convert ksi to MPA

	//LOG << "Regenerator parameters are set. operationMode = " << this->operationMode << ", m_dot_H = " << this->m_dot_H << ", m_dot_C = " << this->m_dot_C;
	//LOG << ", Q_dot_loss = " << this->Q_dot_loss << ", P_0 = " << this->P_0 << ", D_s = " << this->D_s << ", e_v = " << this->e_v;
	//LOG << endl;
}

void RegeneratorModel::setDesignTargets(targetModes::targetModes targetMode, double targetParameter, double dP_max)
{
	this->targetdP_max = dP_max;
	this->targetMode = targetMode;
	this->targetParameter = targetParameter;

	//LOG << "Design targets are set. targetdP_max = " << this->targetdP_max << ", targetMode = " << this->targetMode << ", targetParameter = " << this->targetParameter;
	//LOG << endl;
}

int RegeneratorModel::balanceHeatTransfer_Equation(double T_H_out, double * QdotAsDifference) {

	this->T_H_out = T_H_out;

	try {
		CalculateThermoAndPhysicalModels();
	}
	catch (const invalid_argument& e) {
		return -1;
	}

	*QdotAsDifference = Q_dot_a_calc - Q_dot_a;
	return 0;
}

int RegeneratorModel::balanceHotPressureDrop_Equation(double dP_H, double * dP_HsDifference)
{
	this->dP_H = dP_H;

	if (balanceHeatTransfer->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_HsDifference = dP_H_calc - this->dP_H;
	return 0;
}

int RegeneratorModel::balanceColdPressureDrop_Equation(double dP_C, double * dP_CsDifference)
{
	this->dP_C = dP_C;

	if (balanceHotPressureDrop->solve() != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_CsDifference = dP_C_calc - this->dP_C;
	return 0;
}

int RegeneratorModel::solveForL_Equation(double L, double * dP_max)
{
	this->L = L;

	int status = balanceColdPressureDrop->solve();

	spdlog::get("logger")->debug("\tL = " + std::to_string(L) + " epsilon -> " + std::to_string(epsilon) + " ua -> " + std::to_string(UA)
		+ " dP_max -> " + std::to_string(this->dP_max) + " status  -> " + std::to_string(status));
	

	if (status != C_monotonic_eq_solver::CONVERGED) {
		return -1;
	}

	*dP_max = this->dP_max;

	return 0;
}

int RegeneratorModel::solveForDfr_Equation(double D_fr, double * targetParameter)
{
	this->D_fr = D_fr;

	spdlog::get("logger")->debug(" --> Entered with D_fr = " + std::to_string(D_fr));

	int status = solveForL->solve();

	spdlog::get("logger")->debug("\tepsilon -> " + std::to_string(epsilon) + +" ua -> " + std::to_string(UA) +
		" dP_max -> " + std::to_string(this->dP_max) + " status  -> " + std::to_string(status));


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

int RegeneratorModel::solveForWallThickness_Equation(double th, double * stressAmplitude)
{
	wallThickness = th;

	R_o = R_i + wallThickness;

	double RiSquared = pow(R_i, 2);

	double RoSquared = pow(R_o, 2);

	double RSquareDifference = RoSquared - RiSquared;

	double RSquareProduct = RoSquared * RiSquared;

	double magicPieceRi = RiSquared * RSquareDifference;

	double magicPieceRo = RoSquared * RSquareDifference;

	double magicPieceHigh = RSquareProduct * (Patm - P_high);

	double magicPieceLow = RSquareProduct * (Patm - P_low);

	//sigma_a_h = (P_high*r_i ^ 2 - P_o*r_o ^ 2) / (r_o ^ 2 - r_i ^ 2);		
	double sigma_a_h = (P_high * RiSquared - Patm * RoSquared) / RSquareDifference; //axial stress calculation

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
	double sigma_a_l = (P_low*RiSquared - Patm * RoSquared) / RSquareDifference;	//axial stress calculation

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
	SolverParameters<RegeneratorModel> HT;
	HT.solverName = "Balance Heat Tansfer";
	HT.target = 0;
	HT.guessValue1 = (T_C_in + T_H_in) * 0.5;		HT.guessValue2 = (T_C_in + T_H_in) * 0.7;
	HT.lowerBound = N_co2_props::T_lower_limit;		HT.upperBound = N_co2_props::T_upper_limit;
	HT.tolerance = 0.001;
	HT.iterationLimit = 50;
	HT.isErrorRel = false;
	HT.classInst = this;
	HT.monoEquation = &RegeneratorModel::balanceHeatTransfer_Equation;

	SolverParameters<RegeneratorModel> HPD;
	HPD.solverName = "Balance Hot Pressure Drop";
	HPD.target = 0;
	HPD.guessValue1 = (P_H - 100) / 100.0;		HPD.guessValue2 = P_H / 100.0;
	HPD.lowerBound = 0.1;						HPD.upperBound = P_H;
	HPD.tolerance = 0.01;
	HPD.iterationLimit = 50;
	HPD.isErrorRel = false;
	HPD.classInst = this;
	HPD.monoEquation = &RegeneratorModel::balanceHotPressureDrop_Equation;

	SolverParameters<RegeneratorModel> CPD;
	CPD.solverName = "Balance Cold Pressure Drop";
	CPD.target = 0;
	CPD.guessValue1 = (P_C - 100) / 1000.0;		CPD.guessValue2 = P_C / 1000.0;
	CPD.lowerBound = 0.1;						CPD.upperBound = P_C;
	CPD.tolerance = 0.01;
	CPD.iterationLimit = 50;
	CPD.isErrorRel = false;
	CPD.classInst = this;
	CPD.monoEquation = &RegeneratorModel::balanceColdPressureDrop_Equation;
	
	SolverParameters<RegeneratorModel> LS;
	LS.solverName = "Length Solver";
	LS.target = targetdP_max;
	LS.guessValue1 = 0.8;	LS.guessValue2 = 1;
	LS.lowerBound = 0.1;	LS.upperBound = 10;
	LS.tolerance = 0.01;
	LS.iterationLimit = 50;
	LS.isErrorRel = true;
	LS.classInst = this;
	LS.monoEquation = &RegeneratorModel::solveForL_Equation;

	SolverParameters<RegeneratorModel> DS;
	DS.solverName = "Diameter Solver";
	DS.target = targetParameter;
	DS.guessValue1 = 0.7;	DS.guessValue2 = 0.9;
	DS.lowerBound = 0.1;	DS.upperBound = 10;
	DS.tolerance = 0.001;
	DS.iterationLimit = 50;
	DS.isErrorRel = true;
	DS.classInst = this;
	DS.monoEquation = &RegeneratorModel::solveForDfr_Equation;

	SolverParameters<RegeneratorModel> WT;
	WT.solverName = "Wall Thickness";
	WT.target = stressAmplitude;
	WT.guessValue1 = 0.04;	WT.guessValue2 = 0.05;
	WT.lowerBound = 0;	WT.upperBound = 1;
	WT.tolerance = 0.001;
	WT.iterationLimit = 50;
	WT.isErrorRel = true;
	WT.classInst = this;
	WT.monoEquation = &RegeneratorModel::solveForWallThickness_Equation;

	balanceHeatTransfer = new MonoSolver<RegeneratorModel>(&HT);
	balanceHotPressureDrop = new MonoSolver<RegeneratorModel>(&HPD);
	balanceColdPressureDrop = new MonoSolver<RegeneratorModel>(&CPD);
	solveForL = new MonoSolver<RegeneratorModel>(&LS);
	solveForDfr = new MonoSolver<RegeneratorModel>(&DS);
	solveForWallThickness = new MonoSolver<RegeneratorModel>(&WT);
	
	
	double toleranceD_fr;
	int statusSolverDfr = solveForDfr->solve(&toleranceD_fr);

	if (statusSolverDfr == C_monotonic_eq_solver::CONVERGED) {
		calculateCost();

		return 0;
	}
	else if(toleranceD_fr <= 0.01){

		spdlog::get("logger")->warn("Worse tolerance = " + std::to_string(toleranceD_fr) + " | TargetMode = " + std::to_string(targetMode) +
			", Target = " + std::to_string(targetParameter));

		calculateCost();
		return 1;
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
	solution->costModule = costModule;
	solution->L = L;
	solution->D_fr = D_fr;
	solution->wallThickness = wallThickness;
}
