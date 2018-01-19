#include "MonoSolver.h"
#include "LookupTable_1D.h"
#include "LookupTable_2D.h"
#include "CO2_properties.h"
#include <string>
#include "numeric_solvers.h"
#include "definitions.h"

namespace operationModes {
	enum operationModes {
		PARALLEL,
		REDUNDANT
	};
}
namespace targetModes {
	enum targetModes {
		EPSILON,
		COST,
		UA
	};
}

struct RegeneratorSolution {
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
	double costModule;
	double L;
	double D_fr;
	double wallThickness;
};

class RegeneratorModel
{
private:

	//! Path to the log file
	static string const RegeneratorModel::LOG_FILEPATH;

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
	double T_H_in;

	//! Pressure at hot inlet in [kPa]
	double P_H;

	//! Enthalpy at hot inlet in [kJ/kg]
	double h_H_in;

	//! Mass flow rate of fluid at hot side in [kg/s]
	double m_dot_H;

	//! Temperature at cold inlet in [K]
	double T_C_in;

	//! Pressure at cold inlet in [kPa] 
	double P_C;

	//! Enthalpy at cold inlet in [kJ/kg]
	double h_C_in;

	//! Mass flow rate of fluid at cold side in [kg/s]
	double m_dot_C;


	//	[Outlet conditions]
	//! Temperature at hot outlet in [K]
	double T_H_out;

	//! Pressure at hot outlet in [kPa]
	double P_H_out;

	//! Enthalpy at hot outlet in [kJ/kg]
	double h_H_out;

	//! Temperature at cold outlet in [K]
	double T_C_out;

	//! Pressure at cold outlet in [kPa]
	double P_C_out;

	//! Enthalpy at cold outlet in [kJ/kg]
	double h_C_out;

	//	[Calculation values]
	//! Enthalpy at state fixed by T_C_in and P_H_out in [kJ/kg]
	double h_H_out_max;

	//! Enthalpy at state fixed by T_C_in and P_H in [kJ/kg]
	double h_H_out_max_p;

	/*! \brief Calculated pressure drop at hot side in [kPa]
	*
	*	Calculated by packedspheresFitCO2().
	*	\sa dP_H
	*/
	double dP_H_calc;

	/*! \brief Friction factor for the hot side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double f_H;

	/*! \brief Heat transfer coefficient for the hot side in [W/m^2-K]
	*
	*	Heat transfer coefficientbased on the surface area of the matrix material.
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double h_H;

	/*! \brief Number of transfer units for the hot side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double NTU_H;

	//! Maximum heat transfer if hot stream exiting at cold temperature in [kW]
	double Q_dot_max_H;

	/*! \brief Hot film temperature. [K]
	*
	*	Input to packedspheresFitCO2().
	*/
	double T_H_f;

	/*! \brief Cp over possible temperature span in regenerator for the hot stream in [kJ/kg-K]
	*	\sa h_H_out_max_p
	*/
	double C_p_H;

	//! Capacitance rate of hot stream in [kW/K]
	double C_dot_H;

	//! Enthalpy at state fixed by T_H_in and P_C_out in [kJ/kg]
	double h_C_out_max;

	//! Enthalpy at state fixed by T_H_in and P_C in [kJ/kg]
	double h_C_out_max_p;

	/*! \brief Calculated pressure drop at cold side in [kPa]
	*
	*	Calculated by packedspheresFitCO2().
	*	\sa dP_C
	*/
	double dP_C_calc;

	/*! \brief Friction factor at cold side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double f_C;

	/*! \brief Heat transfer coefficient for the cold side in [W/m^2-K]
	*
	*	Heat transfer coefficientbased on the surface area of the matrix material.
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double h_C;

	/*! \brief Number of transfer units for the cold side. [-]
	*
	*	Calculated by packedspheresFitCO2().
	*/
	double NTU_C;

	//! Maximum heat transfer if cold stream exiting at hot temperature in [kW]
	double Q_dot_max_C;

	/*! \brief Cold film temperature. [K]
	*
	*	Input to packedspheresFitCO2().
	*/
	double T_C_f;

	/*! \brief Cp over possible temperature span in regenerator for the cold stream in [kJ/kg-K]
	*	\sa h_H_out_max_p
	*/
	double C_p_C;

	//! Capacitance rate of cold stream in [kW/K]
	double C_dot_C;

	/*! \brief Maximum of C_dot_C and C_dot_H capacitance rates in [kW/K]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_dot_max;

	/*! \brief Minimum of C_dot_C and C_dot_H capacitance rates in [kW/K]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_dot_min;

	/*! \brief Capacitance ratio. [-]
	*	\sa C_dot_C, C_dot_H
	*/
	double C_R;

	/*! \brief Calculated heat transferred between cold and hot sides with losses substracted in [kW]
	*	\sa Q_dot, Q_dot_a_calc, Q_dot_loss, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_calc;

	/*! \brief Calculated heat transferred between cold and hot sides in [kW]
	*	\sa Q_dot_a, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_a_calc;

	// [Heat Exchanger Design Parameters]

	/*! \brief Target maximum pressure drop in [kPa]
	*	\sa dP_max
	*/
	double targetdP_max;

	/*! \brief Flow switch period in [s]
	*	\sa f_0
	*/
	double P_0;

	//! Volume of regenerator in [m^3]
	double V_0;

	/*! \brief Flow switch frequency in [Hz]
	*	\sa P_0
	*/
	double f_0;

	//! Diameter of bed spheres in [m]
	double D_s;

	//! Ratio of empty volume inside of the regenerator to total regenerator volume. [-]
	double e_v;

	/*! \brief Heat lost during heat exchange between hot stream and cold stream in [kW]
	*	\sa Q_dot, Q_dot_a, Q_dot_calc, Q_dot_a_calc
	*/
	double Q_dot_loss;

	//! Length of the heat exchanger in [m]
	double L;

	/*! \brief This diameter is used in themodynamic model calculations [m].
	*
	*	Diameter of one of the numberOfModulesTotal.
	*	\sa numberOfModulesPerRX, D_shell
	*/
	double D_fr;

	/*! \brief Diameter of the steel shell that accounts for extra insulationThickness added [m].
	*
	*	This diameter is used to calculate stresses and material used. D_shell = D_fr + insulationThickness * 2
	*	\sa numberOfModulesPerRX, D_fr, D_fr, calculateWallThickness()
	*/
	double D_shell;

	/*! \brief Inner radius of the heat exchanger shell in [m]
	*	\sa D_shell
	*/
	double R_i;

	/*! \brief Outer radius of the heat exchanger shell in [m]
	*	\sa D_shell
	*/
	double R_o;

	/*! \brief Frontal area of the heat exchanger in [m^2]
	*	\sa D_fr
	*/
	double A_fr;

	//	[Heat Exchanger Physical Material Properties]

	/*! \brief Film temperature in [K]
	*	\sa T_f_h, T_f_c
	*/
	double T_f;

	//! Heat exchanger material density in [kg/m^3]
	double rho_s;

	//! Heat capacity of the heat exchanger material in [kJ/kg-K]
	double c_s;

	//! Mass of the heat exchanger in [kg]
	double m_s;

	//! Surface area of spheres packed in the bed in [m^2]
	double A_s;

	/*! \brief Matrix capacity ratio. [-]
	*	\sa C_m_e
	*/
	double C_m;

	/*! \brief Effective matrix capacity ratio. [-]
	*	\sa C_m
	*/
	double C_m_e;

	// [Heat Exchanger Metrics]

	/*! \brief Actual thermodynamic efficiency. [-]
	*	\sa targetEffectivness
	*/
	double epsilon;

	/*! \brief Actual maximum pressure drop in [kPa]
	*	It is a maximum pressure drop out of dP_C and dP_H.
	*	\sa targetdP_max
	*/
	double dP_max;

	/*! \brief Heat transferred between cold and hot sides with losses substracted in [kW]
	*	\sa Q_dot_calc, Q_dot_a, Q_dot_loss, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot;

	/*! \brief Maximum possible heat transfer in [kW]
	*	\sa Q_dot_max_H, Q_dot_max_C, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_max;

	/*! \brief Heat transferred between cold and hot sides in [kW]
	*	\sa Q_dot_a_calc, CalculateThermoAndPhysicalModels()
	*/
	double Q_dot_a;

	/*! \brief Actual pressure drop at hot side in [kPa]
	*	\sa dP_H_calc
	*/
	double dP_H;

	/*! \brief Actual pressure drop at cold side in [kPa].
	*	\sa dP_C_calc
	*/
	double dP_C;

	/*! \brief NTU of the cycle. [-]
	*	\sa NTU_R_e, CalculateThermoAndPhysicalModels()
	*/
	double NTU_R;

	/*! \brief Effective NTU of regenerator for unbalanced flow. [-]
	*	\sa CalculateThermoAndPhysicalModels()
	*/
	double NTU_R_e;

	/*! \brief UA of the heat exchanger in [kW/K]
	*	\sa NTU_R, CalculateThermoAndPhysicalModels()
	*/
	double UA;

	/*! \brief Effective efficiency. [-]
	*	\sa epsilon, CalculateThermoAndPhysicalModels()
	*/
	double epsilon_1;

	/*! \brief Correction parameter. [-]
	*	\sa epsilon, epsilon_1, CalculateThermoAndPhysicalModels()
	*/
	double X;

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

	/*! \brief Allows to choose which second design parameter to set.
	*
	*	Three possible modes are: "epsilon", "cost" or "ua"
	*	\sa targetParameter
	*/
	targetModes::targetModes targetMode;

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
	double P_high = 24.94;

	//TODO [MPa]
	double P_low = 8.218;

	/*! \brief Volume of the steel shell per module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeShellTotal
	*/
	double volumeShell;

	/*! \brief Volume of the bed insulation per bed module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeInsulationTotal
	*/
	double volumeInsulation;

	/*! \brief Volume of the bed consisting of spheres per bed module. [m^3]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), volumeBedTotal
	*/
	double volumeBed;

	/*! \brief Mass of the steel shell per moudle. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massShellTotal
	*/
	double massShell;

	/*! \brief Mass of the bed insulation per bed module. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massInsulationTotal
	*/
	double massInsulation;

	/*! \brief Mass of the bed consisting of spheres per bed module. [kg]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), massBedTotal
	*/
	double massBed;

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
	double costShellMaterial;

	/*! \brief Cost of the bed insulation per bed module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costInsulationMaterial;

	/*! \brief Cost of the bed consisting of spheres per bed module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost()
	*/
	double costBedMaterial;

	/*! \brief Static price of welding componets of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double priceWelding = 1254.77 * 2;

	/*! \brief Static price of casting componets of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double priceCasting = 4984.65;

	/*! \brief Static price of casting steel components of the module. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costHXTotal, costPerModuleTotal
	*/
	double priceCastingSteel = 11610.05 + 3285.91 * 2;

	/*! \brief Cost of one bed module material. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleTotal
	*/
	double costMaterial;

	/*! \brief Total cost of the bed module including static costs of welding and casting. [$]
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa calculateCost(), costPerModuleMaterial
	*/
	double costModule;

	/*!	\brief The value of the second target design parameter.
	*	\sa targetMode
	*/
	double targetParameter;

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

	/*! \brief Loads lookup tables into memory for quick access.
	*
	*	Regenerative heat exchanger consists of numberOfModulesTotal bed modules.
	*	\sa LookupTable, bedMaterialTable, spheresRPTable, regeneratorTable
	*/
	void loadTables();

	/*! \breif Calculates total cost of the regenerator.
	*	calculateWallThickness() must be called prior!
	*/
	void calculateCost();

	MonoSolver<RegeneratorModel>* balanceHeatTransfer;
	MonoSolver<RegeneratorModel>* balanceHotPressureDrop;
	MonoSolver<RegeneratorModel>* balanceColdPressureDrop;
	MonoSolver<RegeneratorModel>* solveForL;
	MonoSolver<RegeneratorModel>* solveForDfr;
	MonoSolver<RegeneratorModel>* solveForWallThickness;

public:
	RegeneratorModel();
	~RegeneratorModel();

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
	*	\param m_dot_H Mass flow rate of hot stream in [kg/s]
	*	\param m_dot_C Mass flow rate of cold stream in [kg/s]
	*	\param Q_dot_loss Heat loss rate of heat exchanger in [kW]
	*	\param P_0 Switching period of heat exchanger in [s]. Deremines how often flows through beds is switched
	*	\param D_s Diameter of spheres (particles) packed inside of the heat exchanger in [m]
	*	\param e_v Porosity or ratio of empty space inside of the heat exchanger to its total volume
	*	\sa m_dot_H, m_dot_C, Q_dot_loss, P_0, D_s, e_v, CalculateThermoAndPhysicalModels(), setguesses(), setInletState(), initialize(), setDesignTargets()
	*/
	void setParameters(double m_dot_H, double m_dot_C, double Q_dot_loss, double P_0, double D_s, double e_v);

	/*!	\brief Sets design parameters such as dP_max and epsilon.
	*
	*	\param dP_max Target maximum pressure drop in the system in [kPa]
	*	\param epsilon Target effectiveness of the system
	*	\sa dP_max, epsilon, CalculateThermoAndPhysicalModels(), setParameters(), setInletState(), initialize(), setGuesses()
	*/
	void setDesignTargets(targetModes::targetModes targetMode, double targetParameter, double dP_max);

	int balanceHeatTransfer_Equation(double T_H_out, double * QdotAsDifference);
	int balanceHotPressureDrop_Equation(double dP_H, double * dP_HsDifference);
	int balanceColdPressureDrop_Equation(double dP_C, double * dP_CsDifference);
	int solveForL_Equation(double L, double * dP_max);
	int solveForDfr_Equation(double D_fr, double * targetParameter);
	int solveForWallThickness_Equation(double th, double * stressAmplitude);

	int solveSystem();
	void getSolution(RegeneratorSolution* solution);
};

