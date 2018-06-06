#pragma once
#include <string>
#include "heat_exchangers.h"
#include "definitions.h"
#include "RegeneratorModel.h"

using namespace std;

/*!
  \brief     Regenerative heat exchanger class.
  \details   How to use: 1) Create an instance of the class. 2) Set heat exchanger parameters with setParameters(). 3) Set inlet parameters using setInletState() 
				4) Set design targets with setDesignTargets() method. Now it is ready to go.
  \author    Dmitrii Turygin 	https://github.com/tvdmitrii
  \version   2.0
  \date      1/21/2018
*/
class RegenHX
{
private:

	RegeneratorModel* regenModel;
	valve* valves;

	double dP_C;
	double dP_H;
	double epsilon;
	double T_H_in;
	double T_C_in;
	double T_C_out;
	double T_H_out;
	double NTU_R_e;
	double Q_dot_a;
	double UA;
	double m_dot_H;
	double m_dot_C;
	double L;
	double D_fr;
	double AspectRatio;
	double wallThickness;
	double m_dot_carryover;
	double m_HTR_LP_dP, m_HTR_HP_dP;

	targetModes::targetModes target_1;
	targetModes::target2Modes target_2;
	double target_2_value;
	double P_0;
	double D_s;
	double Q_dot_loss;
	double e_v;

	/*! \brief Number of (hot module + cold module) sets in heat exchanger. [-]
		
		Heat exchanger needs two modules. One for the hot flow and one for cold flow.
		Hot flow + cold flow modules form a set.
		Usually number of sets is 2.
		\sa operationModes
	*/
	int modulesInParallel = 2;

	/*! \brief Sets regenerator operation mode. See operationModes. [-]		
		\sa operationModes
	*/
	operationModes::operationModes operationMode;

	/*! \brief Total cost of the regenerator. [$]
	
		RegeneratorMModel::costModule * numberOfModulesTotal bed modules.
	*/
	double costHX;

	double costValves;

	void resetDesignStructure();

	/*!	\brief Sets fluid states for hot and cold inlets
		
		Mass flows adjusted, depending on operationMode.
		setParameters() prior to calling this method!
		\param T_H_in Temperature of fluid at hot inlet in [K]
		\param P_H_in Pressure of fluid at hot inlet in [kPa]
		\param m_dot_H Mass flow rate of hot stream in [kg/s]
		\param T_C_in Temperature of fluid at cold inlet in [K]
		\param P_C Pressure of fluid at cold inlet in [kPa]
		\param m_dot_C Mass flow rate of cold stream in [kg/s]
		\sa setParameters(), setDesignTargets()
	*/
	void setInletStates(double T_H_in, double P_H_in, double m_dot_H, double T_C_in, double P_C, double m_dot_C);

	/*!	\brief Sets flow and regenerator parameters
	
		\param operationMode Please see operationMode.
		\param Q_dot_loss Heat loss rate of heat exchanger in [kW]
		\param P_0 Switching period of heat exchanger in [s]. Deremines how often flows through beds is switched
		\param D_s Diameter of spheres (particles) packed inside of the heat exchanger in [m]
		\param e_v Porosity or ratio of empty space inside of the heat exchanger to its total volume
		\sa setInletState(), setDesignTargets()
	*/
	void setParameters(operationModes::operationModes operationMod, double Q_dot_loss, double P_0, double D_s, double e_v);

	/*!	\brief Sets design parameters.
	

		\param targetMode Allows to choose between targetModes
		\param targetParameter Can be either effectiveness, cost or UA, depending on targetMode parameter
		\param dP_max Target maximum pressure drop in the system in [kPa]
		\sa setParameters(), setInletState()
	*/
	void setDesignTargets(targetModes::targetModes targetMode, targetModes::target2Modes secondTargetMode, double targetParameter, double secondTargetParameter);

	void setValves(valve* valves);

public:
	/*!	\brief Solves the model so that it meets design parameters.
	*
	*	Method runs initialize() method and then figureOutL() which is a higher iteration loop, which causes
	*	figureOutD_fr() to be called, which calls BalancedPCs(), which calls BalancedPHs(), which calls BalanceQdotAs(),
	*	which calls CalculateThermoAndPhysicalModels() that explicitly solves system of equations describing the system.
	*	These are 5 levels of nested iteration cycles with figureOutL() being the outter cycle and BalanceQdotAs() being the deepest inner cycle.
	*
	*	Call setParameters(), setGuesses() and setDesignTargets() prior to calling this function!
	*
	*	\sa setDesignTargets(), figureOutL(), figureOutD_fr(), BalancedPCs(), BalancedPHs(), BalanceQdotAs(), CalculateThermoAndPhysicalModels()
	*/
	int getDesignSolution();

	void set_params(int target_1, int target_2, int operation_mode, double target_2_value, double P_0, double D_s, double e_v, double Q_dot_loss);
	
	int getDesignSolution(double* results);

	C_HX_counterflow::S_des_solved ms_des_solved;
	C_HX_counterflow::S_od_solved ms_od_solved;
	void initialize(double N_sub_hx);

	double getD_fr() { return D_fr; }
	double getL() { return L; }
	double getAspectRatio() { return AspectRatio; }
	double getCost() { return costHX; }

	//! Constructor that calls RegenHX::loadTables()
	RegenHX();

	~RegenHX();

	double od_delta_p_cold(double m_dot_c /*kg/s*/);
	double od_delta_p_hot(double m_dot_h /*kg/s*/);

	void design_fix_UA_calc_outlet(double UA_target /*kW/K*/, double eff_limit /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);

	void design_fix_TARGET_calc_outlet(int targetType /*-*/, double targetValue /*kW/K or $*/, double eff_limit /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);

	void off_design_solution(double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);
};

