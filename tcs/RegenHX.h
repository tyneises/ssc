#pragma once
#include <string>
#include "heat_exchangers.h"
#include "definitions.h"
#include "RegeneratorModel.h"

using namespace std;

/*!
*  \brief     Regenerative heat exchanger class.
*  \details   How to use: 1) Create an instance of the class. 2) Set inlet parameters using setInletState() method. 3) Set heat exchanger parameters with setParameters().
*				4) Set design targets with setDesignTargets() method. Now it is ready to go.
*  \author    Dmitrii Turygin 	https://github.com/tvdmitrii
*  \version   1.0
*  \date      1/1/2018
*/
class RegenHX
{
private:

	RegeneratorModel* regenModel;

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
	double wallThickness;

	/*! \brief Number of hot+cold sets in heat exchenger. [-]
	*	
	*	Heat exchanger needs two bed modules. One for hot flow and one for cold flow.
	*	Hot flow + cold flow modules form a set.
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell
	*/
	int numberOfSets = 2;

	/*! \brief Sets bed operation mode. [-]
	*
	*	Two available options are either 'redundant' or 'parallel'.
	*	'redundant' - each hot and cold modules can handle full mass-flow and only one set out of
	*	numberOfSets operates at a time.
	*
	*	'parallel' - massflow is split equally between all numberOfSets
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	operationModes::operationModes operationMode;

	/*! \brief Number of bed modules that total (twice the number of numberOfSets). [-]
	*
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	int numberOfModules = 2 * numberOfSets;

	/*! \brief Total cost of the regenerative heat exchanger including static costs of welding and casting. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costHXMaterial, costPerModuleTotal, costPerModuleMaterial
	*/
	double costHX;

public:
	/*!	\brief Sets fluid state at hot and cold inlets
	*
	*	Without these values CalculateThermoAndPhysicalModels() cannot run.
	*	\param T_H_in Temperature of fluid at hot inlet in [K]
	*	\param P_H Pressure of fluid at hot inlet in [kPa]
	*	\param T_C_in Temperature of fluid at cold inlet in [K]
	*	\param P_C Pressure of fluid at cold inlet in [kPa]
	*	\sa T_H_in, P_H, T_C_in, P_C, CalculateThermoAndPhysicalModels(), setParameters(), setParameters(), initialize(), setDesignTargets()
	*/
	void setInletStates(double T_H_in, double P_H, double T_C_in, double P_C);

	/*!	\brief Sets flow and heat exchanger parameters
	*
	*	Without these values CalculateThermoAndPhysicalModels() cannot run.
	*	\param operationMode Please see operationMode.
	*	\param m_dot_H Mass flow rate of hot stream in [kg/s]
	*	\param m_dot_C Mass flow rate of cold stream in [kg/s]
	*	\param Q_dot_loss Heat loss rate of heat exchanger in [kW]
	*	\param P_0 Switching period of heat exchanger in [s]. Deremines how often flows through beds is switched
	*	\param D_s Diameter of spheres (particles) packed inside of the heat exchanger in [m]
	*	\param e_v Porosity or ratio of empty space inside of the heat exchanger to its total volume
	*	\sa m_dot_H, m_dot_C, Q_dot_loss, P_0, D_s, e_v, CalculateThermoAndPhysicalModels(), setguesses(), setInletState(), initialize(), setDesignTargets()
	*/
	void setParameters(operationModes::operationModes operationMode, double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v);

	/*!	\brief Sets design parameters such as dP_max and epsilon.
	*
	*	\param dP_max Target maximum pressure drop in the system in [kPa]
	*	\param epsilon Target effectiveness of the system
	*	\sa dP_max, epsilon, CalculateThermoAndPhysicalModels(), setParameters(), setInletState(), initialize(), setGuesses()
	*/
	void setDesignTargets(targetModes::targetModes targetMode, double targetParameter, double dP_max);

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
	int solveSystem();
	
	int solveSystem(double* results);

	C_HX_counterflow::S_des_solved ms_des_solved;
	C_HX_counterflow::S_od_solved ms_od_solved;

	void initialize(int N_sub_hx);

	//! Constructor that calls RegenHX::loadTables()
	RegenHX();

	~RegenHX();

	double od_delta_p_cold(double m_dot_c /*kg/s*/);
	double od_delta_p_hot(double m_dot_h /*kg/s*/);

	void design_fix_UA_calc_outlet(double UA_target /*kW/K*/, double eff_target /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);

	void off_design_solution(double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);
};

