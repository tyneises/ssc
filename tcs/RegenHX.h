#pragma once
#include "LookupTable_1D.h"
#include "LookupTable_2D.h"
#include "CO2_properties.h"
#include <string>
#include "numeric_solvers.h"
#include "heat_exchangers.h"
#include "definitions.h"

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

#ifdef _DEBUG
	/*! \breif Compiled in _DEBUG only.
	*	Allows to convert C_monotonic_eq_solver::solver_exit_modes values back to strings
	*/
	static const string solver_exit_modes[];

	/*! \breif Compiled in _DEBUG only.
	*	Handle to debug log file
	*/
	ofstream DEBUG_LOG;

	/*! \breif Compiled in _DEBUG only.
	*	Debug log filepath
	*/
	static const string DEBUG_LOG_FILEPATH;

	/*! \breif Compiled in _DEBUG only.
	*	Last core error that caused BalanceQdotAs solver to return -1
	*/
	string lastCoreError;
#endif

	/*! \breif Shows status of BalanceQdotAs solver
	*	/sa BalanceQdotAs()
	*/
	int BalanceQdotAs_Status;

	/*! \breif Shows status of BalancedPHs solver
	*	/sa BalancedPHs()
	*/
	int BalancedPHs_Status;

	/*! \breif Shows status of BalancedPCs solver
	*	/sa BalancedPCs()
	*/
	int BalancedPCs_Status;

	/*! \breif Shows status of figureOutD_fr solver
	*	/sa figureOutD_fr()
	*/
	int figureOutD_fr_Status;

	/*! \breif Shows status of figureOutL solver
	*	/sa figureOutL()
	*/
	int figureOutL_Status;

	/*! \breif Shows status of WallThickness solver
	*	/sa calculateWallThickness()
	*/
	int WallThickness_Status;

	//! Handle to a log file
	ofstream LOG;
	
	//! Log filepath
	static const string LOG_FILEPATH;

	//! String containg path to a folder where property files are.
	static const string PROPERTY_FILES;

	//! String containg path to a lookup CSV file.
	static const string SPHERES_RP_TABLE_PATH;

	//! String containg path to a lookup CSV file.
	static const string BALANCED_REGENERATOR_TABLE_PATH;
	
	/*! \brief String containg path to a lookup CSV file.
	*
	*	First coloumn contains number of cycles. Second coloumn contains SigmaA stress in [ksi].
	*/
	static const string FATIGUE_TABLE_PATH;

	/*! \brief Pointer to a lookup table containing heat exchanger bed material properties.
	*	\sa bedMaterialName
	*/
	LookupTable_1D* bedMaterialTable;

	/*! \brief Pointer to a lookup table containing heat exchanger shell material properties.
	*	\sa shellMaterialName
	*/
	LookupTable_1D* shellMaterialTable;

	/*! \brief Pointer to a lookup table containing heat exchanger insulation material properties.
	*	\sa insulationMaterialName
	*/
	LookupTable_1D* insulationMaterialTable;

	//! Pointer to a lookup table.
	LookupTable_2D* regeneratorTable;

	//! Pointer to a lookup table.
	LookupTable_1D* spheresRPTable;

	//! Pointer to a lookup table. In [ksi]
	LookupTable_1D* fatigueTable;

	//! Contains property functions that are used to set CO2 states.
	CO2_state CO2State;

	//! Used to store error codes produced by CO2_state class methods. 
	int error;

	//! Temperature at hot inlet in [K]
	double T_H_in;			//[K]

	//! Pressure at hot inlet in [kPa]
	double P_H;				//[kPa]

	//! Enthalpy at hot inlet in [kJ/kg]
	double h_H_in;			//[kJ/kg]

	//! Mass flow rate of fluid at hot side in [kg/s]
	double m_dot_H;			//[kg/s]

	//! Temperature at cold inlet in [K]
	double T_C_in;			//[K]

	//! Pressure at cold inlet in [kPa] 
	double P_C;				//[kPa]

	//! Enthalpy at cold inlet in [kJ/kg]
	double h_C_in;			//[kJ/kg]

	//! Mass flow rate of fluid at cold side in [kg/s]
	double m_dot_C;			//[kg/s]


	//	[Outlet conditions]
	//! Temperature at hot outlet in [K]
	double T_H_out;			//[K]

	//! Pressure at hot outlet in [kPa]
	double P_H_out;			//[kPa]

	//! Enthalpy at hot outlet in [kJ/kg]
	double h_H_out;			//[kJ/kg]
	
	//! Temperature at cold outlet in [K]
	double T_C_out;			//[K]

	//! Pressure at cold outlet in [kPa]
	double P_C_out;			//[kPa]

	//! Enthalpy at cold outlet in [kJ/kg]
	double h_C_out;			//[kJ/kg]

	//	[Calculation values]
	//! Enthalpy at state fixed by T_C_in and P_H_out in [kJ/kg]
	double h_H_out_max;		//[kJ/kg]

	//! Enthalpy at state fixed by T_C_in and P_H in [kJ/kg]
	double h_H_out_max_p;	//[kJ/kg]

	/*! \brief Calculated pressure drop at hot side in [kPa]
	*
	*	Calculated by packedspheresFitCO2().
	*	\sa dP_H
	*/
	double dP_H_calc;		//[kPa]

	/*! \brief Friction factor for the hot side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double f_H;				//[-]

	/*! \brief Heat transfer coefficient for the hot side in [W/m^2-K]
	*
	*	Heat transfer coefficientbased on the surface area of the matrix material.
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double h_H;				//[-]

	/*! \brief Number of transfer units for the hot side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double NTU_H;			//[-]

	//! Maximum heat transfer if hot stream exiting at cold temperature in [kW]
	double Q_dot_max_H;		//[kW]

	/*! \brief Hot film temperature. [K]
	*
	*	Input to packedspheresFitCO2().
	*/
	double T_H_f;			//[K]

	/*! \brief Cp over possible temperature span in regenerator for the hot stream in [kJ/kg-K]
	*	\sa h_H_out_max_p
	*/
	double C_p_H;			//[kJ/kg-K]

	//! Capacitance rate of hot stream in [kW/K]
	double C_dot_H;			//[kW/K]

	//! Enthalpy at state fixed by T_H_in and P_C_out in [kJ/kg]
	double h_C_out_max;		//[kJ/kg]

	//! Enthalpy at state fixed by T_H_in and P_C in [kJ/kg]
	double h_C_out_max_p;	//[kJ/kg]

	/*! \brief Calculated pressure drop at cold side in [kPa]
	*
	*	Calculated by packedspheresFitCO2().
	*	\sa dP_C
	*/
	double dP_C_calc;		//[kPa]

	/*! \brief Friction factor at cold side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double f_C;				//[-]

	/*! \brief Heat transfer coefficient for the cold side in [W/m^2-K]
	*
	*	Heat transfer coefficientbased on the surface area of the matrix material.
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double h_C;				//[W/m^2-K]

	/*! \brief Number of transfer units for the cold side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double NTU_C;			//[-]

	//! Maximum heat transfer if cold stream exiting at hot temperature in [kW]
	double Q_dot_max_C;		//[kW]

	/*! \brief Cold film temperature. [K]
	*
	*	Input to packedspheresFitCO2().
	*/
	double T_C_f;			//[K]

	/*! \brief Cp over possible temperature span in regenerator for the cold stream in [kJ/kg-K]
	*	\sa h_H_out_max_p
	*/
	double C_p_C;			//[kJ/kg-K]

	//! Capacitance rate of cold stream in [kW/K]
	double C_dot_C;			//[kW/K]

	/*! \brief Maximum of C_dot_C and C_dot_H capacitance rates in [kW/K]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_dot_max;		//[kW/K]

	/*! \brief Minimum of C_dot_C and C_dot_H capacitance rates in [kW/K]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_dot_min;		//[kW/K]

	/*! \brief Capacitance ratio. [-]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_R;				//[-]

	/*! \brief Calculated heat transferred between cold and hot sides with losses substracted in [kW]
	*	\sa Q_dot, Q_dot_a_calc, Q_dot_loss, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_calc;		//[kW]

	/*! \brief Calculated heat transferred between cold and hot sides in [kW]
	*	\sa Q_dot_a, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_a_calc;	//[kW]
	
	// [Heat Exchanger Design Parameters]

	/*! \brief Target maximum pressure drop in [kPa]
	*	\sa dP_max
	*/
	double targetdP_max;	//[kPa]

	/*! \brief Flow switch period in [s]
	*	\sa f_0
	*/
	double P_0;				//[s]

	//! Volume of regenerator in [m^3]
	double V_0;				//[m^3]

	/*! \brief Flow switch frequency in [Hz]
	*	\sa P_0
	*/
	double f_0;				//[Hz]

	//! Diameter of bed spheres in [m]
	double D_s;				//[m]

	//! Ratio of empty volume inside of the regenerator to total regenerator volume. [-]
	double e_v;				//[-]

	/*! \brief Heat lost during heat exchange between hot stream and cold stream in [kW]
	*	\sa Q_dot, Q_dot_a, Q_dot_calc, Q_dot_a_calc
	*/
	double Q_dot_loss;		//[kW]

	//! Length of the heat exchanger in [m]
	double L;				//[m]

	/*! \brief This diameter is used in themodynamic model calculations [m].
	*
	*	Diameter of one of the numberOfModulesTotal.
	*	\sa numberOfModulesPerRX, D_shell
	*/
	double D_fr;			//[m]

	/*! \brief Diameter of the steel shell that accounts for extra insulationThickness added [m].
	*
	*	This diameter is used to calculate stresses and material used. D_shell = D_fr + insulationThickness * 2
	*	\sa numberOfModulesPerRX, D_fr, D_fr, calculateWallThickness()
	*/
	double D_shell;			//[m]

	/*! \brief Inner radius of the heat exchanger shell in [m]
	*	\sa D_shell
	*/
	double R_i;			//[m]

	/*! \brief Outer radius of the heat exchanger shell in [m]
	*	\sa D_shell
	*/
	double R_o;			//[m]

	/*! \brief Frontal area of the heat exchanger in [m^2]
	*	\sa D_fr
	*/
	double A_fr;			//[m^2]

	//	[Heat Exchanger Physical Material Properties]

	/*! \brief Film temperature in [K]
	*	\sa T_f_h, T_f_c
	*/
	double T_f;				//[K]

	//! Heat exchanger material density in [kg/m^3]
	double rho_s;			//[kg/m^3]

	//! Heat capacity of the heat exchanger material in [kJ/kg-K]
	double c_s;				//[kJ/kg-K]

	//! Mass of the heat exchanger in [kg]
	double m_s;				//[kg]

	//! Surface area of spheres packed in the bed in [m^2]
	double A_s;				//[m^2]

	/*! \brief Matrix capacity ratio. [-]
	*	\sa C_m_e
	*/
	double C_m;				//[-]

	/*! \brief Effective matrix capacity ratio. [-]
	*	\sa C_m
	*/
	double C_m_e;			//[-]

	// [Heat Exchanger Metrics]

	/*! \brief Actual thermodynamic efficiency. [-]
	*	\sa targetEffectivness
	*/
	double epsilon;			//[-]

	/*! \brief Actual maximum pressure drop in [kPa]
	*	It is a maximum pressure drop out of dP_C and dP_H.
	*	\sa targetdP_max
	*/
	double dP_max;			//[kPa]

	/*! \brief Heat transferred between cold and hot sides with losses substracted in [kW]
	*	\sa Q_dot_calc, Q_dot_a, Q_dot_loss, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot;			//[kW]

	/*! \brief Maximum possible heat transfer in [kW]
	*	\sa Q_dot_max_H, Q_dot_max_C, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_max;		//[kW]

	/*! \brief Heat transferred between cold and hot sides in [kW]
	*	\sa Q_dot_a_calc, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_a;			//[kW]

	/*! \brief Actual pressure drop at hot side in [kPa]
	*	\sa dP_H_calc
	*/
	double dP_H;			//[kPa]

	/*! \brief Actual pressure drop at cold side in [kPa].
	*	\sa dP_C_calc
	*/
	double dP_C;			//[kPa]

	/*! \brief NTU of the cycle. [-]
	*	\sa NTU_R_e, CalculateThermoAndPhysicalModels()
	*/
	double NTU_R;			//[-]

	/*! \brief Effective NTU of regenerator for unbalanced flow. [-]
	*	\sa CalculateThermoAndPhysicalModels()
	*/
	double NTU_R_e;			//[-]

	/*! \brief UA of the heat exchanger in [kW/K]
	*	\sa NTU_R, CalculateThermoAndPhysicalModels()
	*/
	double UA;				//[kW/K]

	/*! \brief Effective efficiency. [-]
	*	\sa epsilon, CalculateThermoAndPhysicalModels()
	*/
	double epsilon_1;		//[-]

	/*! \brief Correction parameter. [-]
	*	\sa epsilon, epsilon_1, CalculateThermoAndPhysicalModels()
	*/
	double X;				//[-]

	/*! \brief Number of cycles completed during operationYears of operationHoursPerDay operation. [-]
	*	\sa stressAmplitude, wallThickness, calculateWallThickness()
	*/
	double numberOfCycles;

	/*! \brief Number of years that regenerative heat exchanger should operate. Used to calculate numberOfCycles. [years]
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), operationHoursPerDay
	*/
	double operationYears = 30;

	/*! \brief Number of hours per day that regenerative heat exchanger should operate. Used to calculate numberOfCycles [hrs/day]
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), operationYears
	*/
	double operationHoursPerDay = 10;

	/*! \brief Number of hot+cold sets in heat exchenger. [-]
	*	
	*	Heat exchanger needs two bed modules. One for hot flow and one for cold flow.
	*	Hot flow + cold flow modules form a set.
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell
	*/
	double numberOfSets = 2;

	/*! \brief Sets bed operation mode. [-]
	*
	*	To available options are either 'redundant' or 'parallel'.
	*	'redundant' - each hot and cold modules can handle full mass-flow and only one set out of
	*	numberOfSets operates at a time.
	*
	*	'parallel' - massflow is split equally between all numberOfSets
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	string operationMode;

	/*! \brief Allows to choose which second design parameter to set.
	*
	*	Three possible modes are: "epsilon", "cost" or "ua"
	*	\sa targetParameter
	*/
	string targetMode;

	/*! \brief Number of bed modules that total (twice the number of numberOfSets). [-]
	*
	*	\sa stressAmplitude, wallThickness, calculateWallThickness(), D_fr, D_shell, numberOfSets
	*/
	double numberOfModulesTotal = 2 * numberOfSets;

	/*! \brief Thickness of the insulation on the inside of the D_shell pipe [m]
	*	\sa numberOfModulesTotal, D_fr, D_shell
	*/
	double insulationThickness = 0.05;

	/*! \brief Sets how much of insulation goes on the end half-spheres of the module. [-]
	*
	*	From 0 (no insulation) to 1 (all of the volume).
	*
	*	\sa numberOfModulesTotal, D_fr, D_shell
	*/
	double insulationParameter = 1;

	/*! \brief stressAmplitude is amlitude of stess. [MPa]
	*	\sa FATIGUE_TABLE_PATH, numberOfCycles, wallThickness, calculateWallThickness()
	*/
	double stressAmplitude;

	/*! \brief Wall thickness of regenerator. [m]
	*	\sa stressAmplitude, numberOfCycles, calculateWallThickness()
	*/
	double wallThickness;

	/*! \brief Atmospheric pressure. [MPa]
	*	\sa calculateWallThickness()
	*/
	double Patm = 0.101325;

	//TODO [MPa]
	double P_high = 24.81;
	
	//TODO [MPa]
	double P_low = 8.292;

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

	/*! \brief Volume of the steel shell per module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeShellTotal
	*/
	double volumeShellPerModule;

	/*! \brief Volume of the bed insulation per bed module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeInsulationTotal
	*/
	double volumeInsulationPerModule;

	/*! \brief Volume of the bed consisting of spheres per bed module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeBedTotal
	*/
	double volumeBedPerModule;

	/*! \brief Mass of the steel shell per moudle. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massShellTotal
	*/
	double massShellPerModule;

	/*! \brief Mass of the bed insulation per bed module. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massInsulationTotal
	*/
	double massInsulationPerModule;

	/*! \brief Mass of the bed consisting of spheres per bed module. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massBedTotal
	*/
	double massBedPerModule;

	/*! \brief Material of the regenerative heat exchanger's shell.
	*	\sa calculateCost()
	*/
	string shellMaterialName = "Carbon_steel_AISI1010";

	/*! \brief Material of the regenerative heat exchanger's insulation.
	*	\sa calculateCost()
	*/
	string insulationMaterialName = "Al_oxide_polycryst";

	/*! \brief Material of the regenerative heat exchanger's bed consisting of spheres.
	*	\sa calculateCost()
	*/
	string bedMaterialName = "Stainless_AISI304";

	/*! \brief Cost of the regenerative heat exchanger's shell material per kilogram. [$/kg]
	*	\sa calculateCost()
	*/
	double specificCostShellMaterial = 1.431;

	/*! \brief Cost of the regenerative heat exchanger's insulation material per kilogram. [$/kg]
	*	\sa calculateCost()
	*/
	double specificCostInsulationMaterial = 3.08;

	/*! \brief Cost of the regenerative heat exchanger's bed material per kilogram. [$/kg]
	*	\sa calculateCost()
	*/
	double specificCostBedMaterial = 4.171;

	/*! \brief Cost of the steel shell per bed module. [$]
	*	\sa calculateCost()
	*/
	double costShellMaterialPerModule;

	/*! \brief Cost of the bed insulation per bed module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costInsulationMaterialPerModule;

	/*! \brief Cost of the bed consisting of spheres per bed module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costBedMaterialPerModule;

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

	/*! \brief Static price of welding componets of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double priceWeldingPerModule = 1254.77 * 2;

	/*! \brief Static price of casting componets of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double priceCastingPerModule = 4984.65;

	/*! \brief Static price of casting steel components of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costHXTotal, costPerModuleTotal
	*/
	double priceCastingSteelPerModule = 11610.05 + 3285.91 * 2;

	/*! \brief Cost of one bed module material. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double costPerModuleMaterial;

	/*! \brief Total cost of the bed module including static costs of welding and casting. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleMaterial
	*/
	double costPerModuleTotal;

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

	/*!	\brief The value of the second target design parameter.
	*	\sa targetMode
	*/
	double targetParameter;

	/*! \brief Loads lookup tables into memory for quick access.
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa LookupTable, bedMaterialTable, spheresRPTable, regeneratorTable
	*/
	void loadTables();

	/*! \brief Calculates the friction factor and Colburn j factor.
	*
	*	PackedSpheres_ND returns the friction factor and Colburn j factor (St Pr^(2/3)) for a packed bed consisting of randomly packed spheres.  
	*	The data used in this procedure are from Kays and London, Compact Heat Exchanger, 2nd ed. , McGraw-Hill, 1964, Figure 7-10
	*
	*	\exception invalid_argument exception is thrown if Re value is outside of range.
	*
	*	\param Re Reynolds number. Re should be inside [20, 50 000] range. [-]
	*	\param *f Friction factor. Interpolated from SpheresRP lookup table. This is a returned value. [-]
	*	\param *j_H Colburn j factor. This is a returned value. [-]
	*	\sa spheresRPTable, packedspheresFitCO2()
	*/
	void packedspheresNdFit(double Re, double* f, double* j_H);

	/*! \brief Returns effectiveness of the heat exchanger.
	*
	*	Balanced flow is assumed. Effectivness is interpolated from Balanced Regenerator lookup table using NTU and utilization U.
	*
	*	C_1 is the product of the mass flow rate and specific heat of the fluid during the hot or cold flow period.
	*	The flow is assumed to be balanced so the capacitance rates and flow duration during the hot and cold flow periods are equal.
	*	C_2 is the thermal capacity of the regenerator divided by the duration of the hot or cold flow process.
	*	C_1 and C_2 are in units of [W/K]. However only ratio C_1/C_2 is utilized, so any units can be used as long as they are the same.
	*
	*	\exception invalid_argument exception is thrown if C_1 or C_2 is negative.
	*	\exception invalid_argument exception is thrown if U > 1000.
	*
	*	\param NTU Number of transfer units. [-]
	*	\param C_1 Capacitance rate of the stream in units of [W/K].
	*	\param C_2 Capacitance rate of the stream in units of [W/K].
	*	\return Effectiveness of the heat exchanger from 0 to 1. [-]
	*	\sa regeneratorTable, CalculateThermoAndPhysicalModels()
	*/
	double hx(double NTU, double C_1, double C_2);

	/*! \brief Calculates CO2 properties fixed by T and P.
	*
	*	\exception invalid_argument exception is thrown if CO2States class was not able to fix the state with T and P provided.
	*
	*	\param T Temperature of CO2. [K]
	*	\param P Pressure of CO2. [kPa]
	*	\param *rho Density in [kg/m^3]. This is a returned value.
	*	\param *mu Viscosity in [Pa-s]. This is a returned value.
	*	\param *k Conductivity in [W/m-K]. This is a returned value.
	*	\param *Pr Prandtl number [-]. This is a returned value.
	*	\param *Cp Specific heat at constant pressure in [kJ/kg-K]. This is a returned value.
	*	\sa packedspheresFitCO2()
	*/
	void getpropsregenFitCO2(double T, double P, double* rho, double* mu, double* k, double* Pr, double* Cp);

	/*! 
	*	\brief Calculates the pressure drop and heat transfer coefficient, for a packed bed consisting of randomly packed spheres.
	*
	*	This procedure relies on getpropsregenFitCO2() and packedspheresNdFit() methods.
	*
	*	\exception invalid_argument exception is thrown if getpropsregenFitCO2() or packedspheresNdFit() throws an exception
	*
	*	\param m_dot Fluid flow rate in [kg/s]
	*	\param d Particle (sphere) diameter in [m]
	*	\param A_fr Frontal area in [m^2]
	*	\param L Length of the regenerator in [m]
	*	\param T Fluid temperature in [K]
	*	\param P Fluid pressure in [kPa]
	*	\param porosity Fraction of the volume of voids over the total volume, between 0 and 1 in [-]
	*	\param *f Friction factor. This is a returned value. [-]
	*	\param *h Heat transfer coefficient based on the surface area of the matrix material in [W/m^2-K]. This is a returned value.
	*	\param *NTU Number of transfer units. This is a returned value. [-]
	*	\param *DP Pressure drop in [kPa]. This is a returned value.
	*	\sa getpropsregenFitCO2(), packedspheresNdFit()
	*/
	void packedspheresFitCO2(double m_dot, double d, double A_fr, double L, double T, double P, double porosity, double* f, double* h, double* NTU, double* DP);

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


	/*! \brief Adjusts T_H_out until Q_dot_a and Q_dot_a_calc are within tolerance.
	*
	*	This method is thermodynamical core of the program, containing a full set of equations that describe this regenerative heat exchanger.
	*	This method needs dP_H, dP_C, D_fr, L and T_H_out set. T_H_out is set by a BalanceQdotAsHelper::operator()() that guesses T_H_out 
	*	and represents a monotonic equation where x = T_H_out and y = (Q_dot_a_calc - Q_dot_a). This monotonic equation is solved by
	*	BalanceQdotAs().
	*
	*	This method calculates Q_dot_a using two different ways and also calculates all parameters of the entire model, including
	*	physical heat exchanger properties like mass.
	*
	*	[First way] Knowing pressure drop dP_H and inlet pressure P_H_in, P_H_out can be calculated.
	*	State of the fluid at hot outlet can be fixed by P_H_out and T_H_out. Energy balance gives us Q_dot_a and T_C_out.
	*
	*	[Second way] Knowing D_fr, L, T_H_out and T_C_out packedspheresFitCO2() can calculate properties that allow to find epsilon.
	*	Using epsilon and Q_dot_max, Q_dot_calc and Q_dot_a_calc can be calculated. packedspheresFitCO2() also provides dP_H_calc and
	*	dP_C_calc that are used by BalancedPHs() and BalancedPCs() to make dP_H = dP_H_calc and dP_C = dP_C_calc.
	*
	*	BalanceQdotAs() iteratively ajusts input T_H_out so (Q_dot_a_calc - Q_dot_a) = 0.
	*
	*
	*	\exception invalid_argument exception in case temperature, pressure or epsilon values do not make sense.
	*	The same exception is thrown in case packedspheresFitCO2() throwed an exception.
	*
	*	\sa T_H_out, T_C_out, dP_H, dP_C, D_fr, L, Q_dot_a, Q_dot_a_calc, dP_C_calc, dP_H_calc, BalanceQdotAs(), packedspheresFitCO2()
	*	
	*/
	void CalculateThermoAndPhysicalModels();

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

	//! Returns (RegenHX::Q_dot_a_calc - RegenHX::Q_dot_a) difference
	double getQdotADifference() { return Q_dot_a_calc - Q_dot_a; }

	//! Returns (RegenHX::dP_H_calc - RegenHX::dP_H) difference
	double getdP_HsDifference() { return dP_H_calc - dP_H; }

	//! Returns (RegenHX::dP_C_calc - RegenHX::dP_C) difference
	double getdP_CsDifference() { return dP_C_calc - dP_C; }

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
	void setParameters(string operationMode, double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v);

	/*!	\brief Sets design parameters such as dP_max and epsilon.
	*
	*	\param dP_max Target maximum pressure drop in the system in [kPa]
	*	\param epsilon Target effectiveness of the system
	*	\sa dP_max, epsilon, CalculateThermoAndPhysicalModels(), setParameters(), setInletState(), initialize(), setGuesses()
	*/
	void setDesignTargets(string targetMode, double dP_max, double targetParameter);

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
	void solveSystem(double* results);

	
	//public:

	C_HX_counterflow::S_des_solved ms_des_solved;
	C_HX_counterflow::S_od_solved ms_od_solved;

	void initializeLOGs();
	void initialize(int N_sub_hx);

	//! Constructor that calls RegenHX::loadTables()
	RegenHX();

	~RegenHX();

	double od_delta_p_cold(double m_dot_c /*kg/s*/);
	double od_delta_p_hot(double m_dot_h /*kg/s*/);

	void design_fix_UA_calc_outlet(double Cost_target /*kW/K*/, double eff_target /*-*/, double T_c_in /*K*/, double P_c_in /*kPa*/, 
		double m_dot_c /*kg/s*/, double P_c_out /*kPa*/, double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);

	void off_design_solution(double T_c_in /*K*/, double P_c_in /*kPa*/, double m_dot_c /*kg/s*/, double P_c_out /*kPa*/,
		double T_h_in /*K*/, double P_h_in /*kPa*/, double m_dot_h /*kg/s*/, double P_h_out /*kPa*/,
		double & q_dot /*kWt*/, double & T_c_out /*K*/, double & T_h_out /*K*/);
};

