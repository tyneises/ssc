#include "RegenHX.h"

//Root is at YOUR_SAM_DIR\ssc\examples

RegenHX::RegenHX()
{

}

RegenHX::~RegenHX()
{
}

void RegenHX::initialize(int N_sub_hx)
{
	
}

void RegenHX::solveSystem()
{
	figureOutL();
	calculateWallThickness();
	calculateCost();
}

void RegenHX::solveSystem(double* results)
{
	solveSystem();
	results[0] = epsilon;
	results[1] = costHXTotal;

	if (this->operationMode == PARALLEL_OM) {
		results[2] = UA * numberOfSets;
	}
	else {
		results[2] = UA;
	}

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

void RegenHX::design_fix_UA_calc_outlet(double UA_target /*kW/K*/, double eff_target /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
	double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
	double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/)
{
	setInletState(T_h_in, P_h_in, T_c_in, P_c_in);
	double Q_dot_loss = 100;
	double P_0 = 45;
	double D_s = 0.003;
	double e_v = 0.37;
	double dP_max = P_h_in - P_h_out;
	setParameters(PARALLEL_OM, m_dot_h, m_dot_c, Q_dot_loss, P_0, D_s, e_v);
	setDesignTargets(UA_TM, UA_target, dP_max);
	
	solveSystem();
	q_dot = Q_dot_a * numberOfSets;
	T_c_out = T_C_out;
	T_h_out = T_H_out;

	ms_des_solved.m_DP_cold_des = dP_C;
	ms_des_solved.m_DP_hot_des = dP_H;
	ms_des_solved.m_eff_design = epsilon;
	//???
	ms_des_solved.m_min_DT_design = min((T_H_in - T_H_out), (T_C_out - T_C_in));
	//???
	ms_des_solved.m_NTU_design = NTU_R_e * numberOfSets;
	ms_des_solved.m_Q_dot_design = Q_dot_a * numberOfSets;
	ms_des_solved.m_T_c_out = T_C_out;
	ms_des_solved.m_T_h_out = T_H_out;
	ms_des_solved.m_UA_design_total = UA * numberOfSets;
}

void RegenHX::off_design_solution(double T_c_in, double P_c_in, double m_dot_c, double P_c_out, double T_h_in, double P_h_in, double m_dot_h, double P_h_out, double & q_dot, double & T_c_out, double & T_h_out)
{
	/*T_C_in = T_c_in;
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
	ms_od_solved.m_UA_total = UA;*/
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
