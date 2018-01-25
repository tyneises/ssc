#include "RegenHX.h"

//Root is at YOUR_SAM_DIR\ssc\examples

RegenHX::RegenHX()
{
	regenModel = new RegeneratorModel();
}

RegenHX::~RegenHX()
{
}

void RegenHX::initialize(int N_sub_hx)
{
	
}

int RegenHX::solveSystem()
{
	RegeneratorSolution solution;
	int status = regenModel->solveSystem();
	//spdlog::get("logger")->info("Regenerator exited with status = " + std::to_string(status));
	if (status < 0) {
		return status;
	}

	regenModel->getSolution(&solution);

	costHX = solution.costModule * numberOfModules;
	T_H_out = solution.T_H_out;
	T_C_out = solution.T_C_out;
	dP_H = solution.dP_H;
	dP_C = solution.dP_C;
	epsilon = solution.epsilon;
	NTU_R_e = solution.NTU_R_e;
	Q_dot_a = solution.Q_dot_a;
	UA = solution.UA;
	L = solution.L;
	D_fr = solution.D_fr;
	wallThickness = solution.wallThickness;

	if (operationMode == operationModes::PARALLEL) {
		UA *= numberOfSets;
		Q_dot_a *= numberOfSets;
	}

	return status;
}

int RegenHX::solveSystem(double* results)
{
	int status = solveSystem();

	if (status < 0) {
		results[0] = -1;
		results[1] = -1;
		results[2] = -1;
		results[3] = -1;
		results[4] = -1;
		results[5] = -1;
		results[6] = -1;
		results[7] = -1;
		results[8] = -1;

		return status;
	}

	results[0] = epsilon;
	results[1] = costHX;
	results[2] = UA;
	results[3] = T_H_out;
	results[4] = dP_H;
	results[5] = dP_C;
	results[6] = D_fr;
	results[7] = L;
	results[8] = wallThickness;

	return status;
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
	double Q_dot_loss = 100;
	double P_0 = 45;
	double D_s = 0.003;
	double e_v = 0.37;
	double dP_max = P_h_in - P_h_out;
	if (dP_max < 220) {
		dP_max = 220;
	}

	setParameters(operationModes::PARALLEL, Q_dot_loss, P_0, D_s, e_v);
	setInletStates(T_h_in, P_h_in, m_dot_h, T_c_in, P_c_in, m_dot_c);
	setDesignTargets(targetModes::UA, UA_target, dP_max);
	
	int status = solveSystem();

	if (status < 0) {
		resetDesignStructure();
		return;
	}

	q_dot = Q_dot_a;
	T_c_out = T_C_out;
	T_h_out = T_H_out;

	ms_des_solved.m_DP_cold_des = 0;
	ms_des_solved.m_DP_hot_des = 0;
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

void RegenHX::resetDesignStructure()
{
	ms_des_solved.m_DP_cold_des = ms_des_solved.m_DP_hot_des = ms_des_solved.m_eff_design = ms_des_solved.m_min_DT_design = ms_des_solved.m_NTU_design =
		ms_des_solved.m_Q_dot_design = ms_des_solved.m_T_c_out = ms_des_solved.m_T_h_out = ms_des_solved.m_UA_design_total = std::numeric_limits<double>::quiet_NaN();
}

void RegenHX::setInletStates(double T_H_in, double P_H, double m_dot_H, double T_C_in, double P_C, double m_dot_C)
{
	this->T_C_in = T_C_in;
	this->T_H_in = T_H_in;
	this->m_dot_H = m_dot_H;
	this->m_dot_C = m_dot_C;

	if (this->operationMode == operationModes::PARALLEL) {
		regenModel->setInletStates(T_H_in, P_H, m_dot_H / numberOfSets, T_C_in, P_C, m_dot_C / numberOfSets);
	}
	else {
		regenModel->setInletStates(T_H_in, P_H, m_dot_H, T_C_in, P_C,  m_dot_C);
	}
}

void RegenHX::setParameters(operationModes::operationModes operationMode, double Q_dot_loss, double P_0, double D_s, double e_v)
{
	this->operationMode = operationMode;
	regenModel->setParameters(Q_dot_loss, P_0, D_s, e_v);
}

void RegenHX::setDesignTargets(targetModes::targetModes targetMode, double targetParameter, double dP_max)
{
	if (targetMode == targetModes::UA && this->operationMode == operationModes::PARALLEL) {
		regenModel->setDesignTargets(targetMode, targetParameter / numberOfSets, dP_max);
	}
	else if(targetMode == targetModes::COST && this->operationMode == operationModes::PARALLEL){
		regenModel->setDesignTargets(targetMode, targetParameter / numberOfModules, dP_max);
	}
	else {
		regenModel->setDesignTargets(targetMode, targetParameter, dP_max);
	}
}
