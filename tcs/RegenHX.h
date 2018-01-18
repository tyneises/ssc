#pragma once
#include "LookupTable_1D.h"
#include "LookupTable_2D.h"
#include "CO2_properties.h"
#include <string>
#include "numeric_solvers.h"
#include "heat_exchangers.h"
#include "definitions.h"
#include <iostream>
#include <ctime>

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
public:
	
	enum operationModes {
		PARALLEL_OM,
		REDUNDANT_OM
	};

private:

	/*! \brief Number of hot+cold sets in heat exchenger. [-]
	*	
	*	Heat exchanger needs two bed modules. One for hot flow and one for cold flow.
	*	Hot flow + cold flow modules form a set.
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell
	*/
	double numberOfSets = 2;

	/*! \brief Sets bed operation mode. [-]
	*
	*	Two available options are either 'redundant' or 'parallel'.
	*	'redundant' - each hot and cold modules can handle full mass-flow and only one set out of
	*	numberOfSets operates at a time.
	*
	*	'parallel' - massflow is split equally between all numberOfSets
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	operationModes operationMode;

	/*! \brief Number of bed modules that total (twice the number of numberOfSets). [-]
	*
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	double numberOfModulesTotal = 2 * numberOfSets;

	/*! \brief Volume of the regenerative heat exchanger's shell. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeShellPerModule
	*/
	double volumeShellTotal;

	/*! \brief Volume of the regenerative heat exchanger's insulation. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeInsulationPerModule
	*/
	double volumeInsulationTotal;

	/*! \brief Volume of the regenerative heat exchanger's bed consisting of spheres. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeBedPerModule
	*/
	double volumeBedTotal;

	/*! \brief Mass of the regenerative heat exchanger's shell. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massShellPerModule
	*/
	double massShellTotal;

	/*! \brief Mass of the regenerative heat exchanger's insulation. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massInsulationPerModule
	*/
	double massInsulationTotal;

	/*! \brief Mass of the regenerative heat exchanger's bed consisting of spheres. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massBedPerModule
	*/
	double massBedTotal;

	/*! \brief Cost of the regenerative heat exchanger's shell. [$]
	*	\sa calculateCost()
	*/
	double costShellMaterialTotal;

	/*! \brief Cost of the regenerative heat exchanger's insulation. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costInsulationMaterialTotal;

	/*! \brief Cost of the regenerative heat exchanger's bed consisting of spheres. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costBedMaterialTotal;

	/*! \brief Cost of the regenerative heat exchanger material. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costHXTotal, costPerModuleMaterial, costPerModuleTotal
	*/
	double costHXMaterial;

	/*! \brief Total cost of the regenerative heat exchanger including static costs of welding and casting. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costHXMaterial, costPerModuleTotal, costPerModuleMaterial
	*/
	double costHXTotal;

	/*! \brief A helper class inherited from C_monotonic_equation class. This class is encloses the 1st layer for iteratively solved equations.
	*
	*	Next level is BalancedPHsHelper.
	*
	*	BalanceQdotAsHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is a wrapper for RegenHX::CalculateThermoAndPhysicalModels()
	*	\sa RegenHX::CalculateThermoAndPhysicalModels(), BalanceQdotAsHelper::operator()()
	*/
	class BalanceQdotAsHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
		RegenHX* system;

	public:
		/*! \brief Lower bound value for temperature of hot outlet in [K]
		*
		*	This limitation is set by CO2_state class.
		*	Used by RegenHX::BalanceQdotAs()
		*/
		double const T_H_out_lowerBound = N_co2_props::T_lower_limit;	//[K]

		/*! \brief Upper bound value for temperature of hot outlet in [K]
		*
		*	This limitation is set by CO2_state class.
		*	Used by RegenHX::BalanceQdotAs()
		*/
		double const T_H_out_upperBound = N_co2_props::T_upper_limit;	//[K]

		/*! \brief Smaller Temperature guess value in [K]
		*
		*	Used in RegenHX::BalanceQdotAs(). Set in BalanceQdotAsHelper().
		*	\sa T_H_out_guess2
		*/
		double T_H_out_guess1;	//[K]

		/*! \brief Larger Temperature guess value in [K]
		*
		*	Used in RegenHX::BalanceQdotAs(). Set in BalanceQdotAsHelper().
		*	\sa T_H_out_guess1
		*/
		double T_H_out_guess2;	//[K]

		/*! \brief Solver tolerance or allowed error
		*
		*	Used in RegenHX::BalanceQdotAs().
		*/
		double const tolerance = 0.001;

		/*! \brief Constructor sets guess values based on data provided from RegenHX class.
		*
		*	T_H_out_guess1 is set to (RegenHX::T_C_in + RegenHX::T_H_in) * 0.5 value.
		*	T_H_out_guess2 is set to (RegenHX::T_C_in + RegenHX::T_H_in) * 0.7 value.
		*
		*	\param *system Pointer to RegenHX class.
		*/
		BalanceQdotAsHelper(RegenHX* system)
		{
			this->system = system;
			T_H_out_guess1 = (this->system->T_C_in + this->system->T_H_in) * 0.5;
			T_H_out_guess2 = (this->system->T_C_in + this->system->T_H_in) * 0.7;
		}

		/*! \brief Monotonic equation that is a wrapper for RegenHX::CalculateThermoAndPhysicalModels()
		*
		*	Sets T_H_out, calls RegenHX::CalculateThermoAndPhysicalModels() and calculates difference Q_dot_a_calc - Q_dot_a.
		*	For more information please refer to RegenHX::CalculateThermoAndPhysicalModels()
		*	\param T_H_out Temperature of hot outlet in [K]
		*	\param *QdotAsDifference Difference between RegenHX::Q_dot_a_calc and RegenHX::Q_dot_a in [kW].
		*	This is a returned value
		*	\return -1 if RegenHX::CalculateThermoAndPhysicalModels() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double T_H_out, double * QdotAsDifference);
	};

	/*! \brief A helper class inherited from C_monotonic_equation class. This class is encloses the 2nd layer for iteratively solved equations.
	*
	*	Previous level is BalanceQdotAsHelper, next level is BalancedPCsHelper.
	*	Almost identical to BalancedPCsHelper class.
	*
	*	BalancedPHsHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is a wrapper for RegenHX::BalanceQdotAs()
	*	\sa RegenHX::BalanceQdotAs(), BalancedPHsHelper::operator()()
	*/
	class BalancedPHsHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
		RegenHX* system;

	public:

		/*! \brief Lower bound value for pressure drop in [kPa]
		*
		*	Used in RegenHX::BalancedPHs(). 
		*/
		double const dP_H_lowerBound = 0.01;	//[kPa]

		/*! \brief Upper bound value for pressure drop in [kPa]
		*
		*	Used in RegenHX::BalancedPHs(). Set in BalancedPHsHelper()
		*/
		double dP_H_upperBound;	//[kPa]

		/*! \brief Smaller Pressure drop guess value in [kPa]
		*
		*	Used in RegenHX::BalancedPHs(). Set in BalancedPHsHelper()
		*	\sa dP_H_guess2
		*/
		double dP_H_guess1;	  //[kPa]

		/*! \brief Larger Pressure drop guess value in [kPa]
		*
		*	Used in RegenHX::BalancedPHs(). Set in BalancedPHsHelper()
		*	\sa dP_H_guess1
		*/
		double dP_H_guess2;	//[kPa]

		/*! \brief Solver tolerance or allowed error
		*
		*	Used in RegenHX::BalanceQdotAs().
		*/
		double const tolerance = 0.01;

		/*! \brief Constructor sets guess values based on data provided from RegenHX class.
		*
		*	dP_H_guess1 is set to (RegenHX::dP_H_upperBound - 100) / 100.0 value.
		*	dP_H_guess2 is set to RegenHX::dP_H_upperBound / 100.0 value.
		*	dP_H_upperBound is set to RegenHX::P_H value. 
		*
		*	\param *system Pointer to RegenHX class.
		*/
		BalancedPHsHelper(RegenHX* system)
		{
			this->system = system;
			this->dP_H_upperBound = this->system->P_H;
			this->dP_H_guess1 = (this->dP_H_upperBound - 100) / 100.0;
			this->dP_H_guess2 = this->dP_H_upperBound / 100.0;

			
		}

		/*! \brief Monotonic equation that is a wrapper for RegenHX::BalanceQdotAs()
		*
		*	For more information please refer to RegenHX::BalanceQdotAs()
		*	\param dP_H Pressure drop in hot stream in [kPa]
		*	\param *dP_HsDifference Difference between RegenHX::dP_H_calc and RegenHX::dP_H in [kPa].
		*	This is a returned value
		*	\return -1 if RegenHX::BalanceQdotAs() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double dP_H, double * dP_HsDifference);
	};

	/*! \brief A helper class inherited from C_monotonic_equation class. This class is encloses the 3rd layer for iteratively solved equations.
	*
	*	Previous level is BalancedPHsHelper, next level is FigureOutD_frHelper.
	*	Almost identical to BalancedPHsHelper class.
	*
	*	BalancedPCsHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is a wrapper for RegenHX::BalancedPHs()
	*	\sa RegenHX::BalancedPHs(), BalancedPCsHelper::operator()()
	*/
	class BalancedPCsHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
		RegenHX* system;

	public:
		/*! \brief Lower bound value for pressure drop in [kPa]
		*
		*	Used in RegenHX::BalancedPCs().
		*/
		double const dP_C_lowerBound = 0.01;	//[kPa]

		/*! \brief Upper bound value for pressure drop in [kPa]
		*
		*	Used in RegenHX::BalancedCHs(). Set in BalancedPCsHelper()
		*/
		double dP_C_upperBound;	//[kPa]

		/*! \brief Smaller Pressure drop guess value in [kPa]
		*
		*	Used in RegenHX::BalancedPCs(). Set in BalancedPCsHelper()
		*	\sa dP_C_guess2
		*/
		double dP_C_guess1;	//[kPa]

		/*! \brief Larger Pressure drop guess value in [kPa]
		*
		*	Used in RegenHX::BalancedPCs(). Set in BalancedPCsHelper()
		*	\sa dP_C_guess1
		*/
		double dP_C_guess2;	//[kPa]

		/*! \brief Solver tolerance or allowed error
		*
		*	Used in RegenHX::BalancedPCs().
		*/
		double const tolerance = 0.01;

		/*! \brief Constructor sets guess values based on data provided from RegenHX class.
		*
		*	dP_C_guess1 is set to (RegenHX::dP_C_upperBound - 100) / 1000.0 value.
		*	dP_C_guess2 is set to RegenHX::dP_C_upperBound / 1000.0 value.
		*	dP_C_upperBound is set to RegenHX::P_C value.
		*
		*	\param *system Pointer to RegenHX class.
		*/
		BalancedPCsHelper(RegenHX* system)
		{
			this->system = system;
			this->dP_C_upperBound = this->system->P_C;

			this->dP_C_guess1 = (this->dP_C_upperBound - 100) / 1000.0;
			this->dP_C_guess2 = this->dP_C_upperBound / 1000.0;
			
			
		}

		/*! \brief Monotonic equation that is a wrapper for RegenHX::BalancedPHs()
		*
		*	For more information please refer to RegenHX::BalancedPHs()
		*	\param dP_C Pressure drop in cold stream in [kPa]
		*	\param *dP_CsDifference Difference between RegenHX::dP_C_calc and RegenHX::dP_C in [kPa].
		*	This is a returned value
		*	\return -1 if RegenHX::BalancedPHs() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double dP_C, double * dP_CsDifference);
	};

	/*! \brief A helper class inherited from C_monotonic_equation class. This class is encloses the 4th layer for iteratively solved equations.
	*
	*	Previous level is BalancedPCsHelper, next level is FigureOutLHelper.
	*
	*	FigureOutD_frHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is a wrapper for RegenHX::BalancedPCs()
	*	\sa RegenHX::BalancedPCs(), FigureOutD_frHelper::operator()()
	*/
	class FigureOutD_frHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
		RegenHX* system;
	
	public:
		/*! \brief Lower bound value for frontal diameter of heat exchanger in [m]
		*
		*	Used in RegenHX::figureOutD_fr().
		*/
		double const D_fr_lowerBound = 0.1;//[m]

		/*! \brief Upper bound value for frontal diameter of heat exchanger in [m]
		*
		*	Used in RegenHX::figureOutD_fr().
		*/
		double const D_fr_upperBound = 10;	//[m]

		/*! \brief Smaller Frontal diameter guess value in [m]
		*
		*	Used in RegenHX::figureOutD_fr().
		*	\sa D_fr_guess2
		*/
		const double D_fr_guess1 = 0.7;	//[m]

		/*! \brief Larger Frontal diameter guess value in [m]
		*
		*	Used in RegenHX::figureOutD_fr().
		*	\sa D_fr_guess1
		*/
		const double D_fr_guess2 = 0.9;	//[m]

		/*! \brief Solver tolerance or allowed error
		*
		*	Used in RegenHX::figureOutD_fr().
		*/
		double const tolerance = 0.001;

		/*! \brief Constructor sets FigureOutD_frHelper::system.
		*
		*	\param *system Pointer to RegenHX class.
		*/
		FigureOutD_frHelper(RegenHX* system)
		{
			this->system = system;
		}

		/*! \brief Monotonic equation that is a wrapper for RegenHX::BalancedPCs()
		*
		*	For more information please refer to RegenHX::BalancedPCs()
		*	\param D_fr Frontal diameter of the heat exchanger in [m]
		*	\param *epsilon Effectiveness of the heat exchanger
		*	This is a returned value
		*	\return -1 if RegenHX::BalancedPCs() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double D_fr, double * targetParameter);
	};

	/*! \brief A helper class inherited from C_monotonic_equation class. This class is encloses the 5th layer for iteratively solved equations.
	*
	*	Previous level is FigureOutD_frHelper.
	*
	*	FigureOutD_frHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is a wrapper for RegenHX::figureOutD_fr()
	*	\sa RegenHX::figureOutD_fr(), FigureOutLHelper::operator()()
	*/
	class FigureOutLHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
		RegenHX* system;
		
	public:
		/*! \brief Lower bound value for length of heat exchanger in [m]
		*
		*	Used in RegenHX::figureOutL().
		*/
		double const L_lowerBound = 0.1;//[m]

		/*! \brief Upper bound value for length of heat exchanger in [m]
		*
		*	Used in RegenHX::figureOutL().
		*/
		double const L_upperBound = 10;	//[m]


		/*! \brief Smaller Length guess value in "Large Guess Value Set" in [m]
		*
		*	Used in RegenHX::figureOutL().
		*	\sa L_guess2
		*/
		const double L_guess1 = 0.8;	//[m]

		/*! \brief Larger Length guess value in "Large Guess Value Set" in [m]
		*
		*	Used in RegenHX::figureOutL().
		*	\sa L_guess1
		*/
		const double L_guess2 = 1;	//[m]

		/*! \brief Smaller Length guess value in "Small Guess Value Set" in [m]
		*
		*	Used in RegenHX::figureOutL().
		*	\sa L_guess4
		*/
		const double L_guess3 = 0.4;	//[m]

		/*! \brief Larger Length guess value in "Small Guess Value Set" in [m]
		*
		*	Used in RegenHX::figureOutL().
		*	\sa L_guess3
		*/
		const double L_guess4 = 0.6;	//[m]

		/*! \brief Solver tolerance or allowed error
		*
		*	Used in RegenHX::figureOutL().
		*/
		double const tolerance = 0.01;

		/*! \brief Constructor sets FigureOutLHelper::system.
		*
		*	\param *system Pointer to RegenHX class.
		*/
		FigureOutLHelper(RegenHX* system)
		{
			this->system = system;
		}

		/*! \brief Monotonic equation that is a wrapper for RegenHX::figureOutD_fr()
		*
		*	For more information please refer to RegenHX::figureOutD_fr()
		*	\param L Length of the heat exchanger in [m]
		*	\param *dP_max Maximum pressure drop in heat exchanger. Maximum of RegenHX::dP_H and RegenHX::dP_C
		*	This is a returned value
		*	\return -1 if RegenHX::figureOutD_fr() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double L, double * dP_max);
	};

	/*! \brief This is a helper class. It is used to find wall thickness of the regenerator.
	*
	*	Previous level is FigureOutLHelper.
	*
	*	WallThicknessLHelper::operator()() is overloaded to be used with C_monotonic_eq_solver and is used in RegenHX::calculateWallThickness()
	*	\sa RegenHX::calculateWallThickness(), WallThicknessLHelper::operator()()
	*/
	class WallThicknessLHelper : public C_monotonic_equation
	{
	private:
		//! Pointer to a RegenHX class, containing data about the cycle and fluid states.
	RegenHX* system;

	public:
		
		//! Lower bound for wall thickness. Set to 0 [m]
		double const th_lowerBound = 0;//[m]

		//! Upper bound for wall thickness. Set to RegenHX::D_shell / 2.0			
		double th_upperBound;	//[m]

		//! Smaller guess value for wall thickness [m]			
		const double th_guess1 = 0.04;	//[m]

		//! Larger guess value for wall thickness [m]					
		const double th_guess2 = 0.05;	//[m]

		//! Solver tolerance or allowed error					
		double const tolerance = 0.001;

		
		WallThicknessLHelper(RegenHX* system)
		{
			this->system = system;
			th_upperBound = system->D_shell / 2.0;
		}

		/*! \brief Calculates stressAmplitude based on thickness provided
		*
		*	For more information please refer to RegenHX::figureOutD_fr()
		*	\param L Length of the heat exchanger in [m]
		*	\param *dP_max Maximum pressure drop in heat exchanger. Maximum of RegenHX::dP_H and RegenHX::dP_C
		*	This is a returned value
		*	\return -1 if RegenHX::figureOutD_fr() throwed an exception. 0 otherwise
		*/
		virtual int operator()(double th, double * stressAmplitude);
	};

	/*!	\brief Adjusts T_H_out until Q_dot_a_calc and Q_dot_a match.
	*
	*	Solves monotonic equation BalanceQdotAsHelper::operator()(). Solver parameters are set in BalanceQdotAsHelper class.
	*
	*	\exception invalid_argument exception is thrown in case solver did not converge.
	*	\sa T_H_out, Q_dot_a_calc, Q_dot_a
	*/
	void BalanceQdotAs();

	/*!	\brief Adjusts dP_H until it matches dP_H_calc.
	*
	*	Solves monotonic equation BalancedPHsHelper::operator()(). Solver parameters are set in BalancedPHsHelper class.
	*
	*	\exception invalid_argument exception is thrown in case solver did not converge.
	*	\sa dP_H, dP_H_calc
	*/
	void BalancedPHs();

	/*!	\brief Adjusts dP_C until it matches dP_C_calc.
	*
	*	Solves monotonic equation BalancedPCsHelper::operator()(). Solver parameters are set in BalancedPCsHelper class.
	*
	*	\exception invalid_argument exception is thrown in case solver did not converge.
	*	\sa dP_C, dP_C_calc
	*/
	void BalancedPCs();

	/*!	\brief Adjusts D_fr until epsilon matches targetEffectiveness.
	*
	*	Solves monotonic equation figureOutD_frHelper::operator()(). Solver parameters are set in figureOutD_frHelper class.
	*
	*	\exception invalid_argument exception is thrown in case solver did not converge.
	*	\sa D_fr, epsilon, targetEffectiveness
	*/
	void figureOutD_fr();

	/*!	\brief Adjusts L until dP_max matches targetdP_max.
	*
	*	Solves monotonic equation figureOutLHelper::operator()(). Solver parameters are set in figureOutLHelper class.
	*
	*	\exception invalid_argument exception is thrown in case solver did not converge.
	*	\sa L, dP_max, targetdP_max
	*/
	void figureOutL();

	//! Calculates thickness of walls of the regenerator, required to operate for numberOfCycles
	void calculateWallThickness();

	/*! \breif Calculates total cost of the regenerator.
	*	calculateWallThickness() must be called prior!
	*/
	void calculateCost();

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
	void setInletState(double T_H_in, double P_H, double T_C_in, double P_C);

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
	void setParameters(operationModes operationMode, double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v);

	/*!	\brief Sets design parameters such as dP_max and epsilon.
	*
	*	\param dP_max Target maximum pressure drop in the system in [kPa]
	*	\param epsilon Target effectiveness of the system
	*	\sa dP_max, epsilon, CalculateThermoAndPhysicalModels(), setParameters(), setInletState(), initialize(), setGuesses()
	*/
	void setDesignTargets(targetModes targetMode, double targetParameter, double dP_max);

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
	void solveSystem();
	
	void solveSystem(double* results);

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

