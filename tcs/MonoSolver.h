#pragma once
#include "numeric_solvers.h"
#include "spdlog/spdlog.h"
#include <stdlib.h>

using namespace std;

template <class T>
class MonoSolver
{
private:
	
	string solverName;

	double target;

	double* guessValues;

	double lowerBound;

	double upperBound;

	double tolerance;

	int iterationLimit;

	bool isErrorRel;

	C_monotonic_equation* equation;

	C_monotonic_eq_solver* eqSolver;

	void setSettings()
	{
		eqSolver->reset();
		eqSolver->settings(tolerance, iterationLimit, lowerBound, upperBound, isErrorRel);
	}
	
public:
	MonoSolver(string solverName, double target, double* guessValues, double lowerBound, double upperBound, double tolerance,
		int iterationLimit, bool isErrorRel, T* chain, int (T::*monoEquation)(double x, double* y))
	{
		this->solverName = solverName;
		this->target = target;
		this->guessValues = guessValues;
		this->lowerBound = lowerBound;
		this->upperBound = upperBound;
		this->tolerance = tolerance;
		this->iterationLimit = iterationLimit;
		this->isErrorRel = isErrorRel;

		equation = new C_member_mono_eq<T>(chain, monoEquation);
		eqSolver = new C_monotonic_eq_solver(*equation);
	}
	
	int solve(double* tolSolved = nullptr) {
		if (sizeof(*guessValues) < 2) {
			throw invalid_argument("MonoSolver requires at least two guess values!");
		}

		setSettings();

		double solution, toleranceSolved;
		solution = toleranceSolved = numeric_limits<double>::quiet_NaN();
		int iter_solved = -1;

		double y1, y2;

		double guess1 = guessValues[0];
		double guess2 = guessValues[1];
		int guess1_Result, guess2_Result;
		int i = 0;

		/*while (i + 1 < sizeof(*guessValues)) {
			guess1_Result = eqSolver->test_member_function(guess1, &y1);
			guess2_Result = eqSolver->test_member_function(guess2, &y2);

			if (guess1_Result != 0 || guess2_Result != 0) {
				i += 2;
				continue;
			}

			C_monotonic_eq_solver::S_xy_pair pair1;
			C_monotonic_eq_solver::S_xy_pair pair2;

			pair1.x = guess1;
			pair1.y = y1;
			pair2.x = guess2;
			pair2.y = y2;

			int status = eqSolver->solve(pair1, pair2, target, solution, toleranceSolved, iter_solved);*/
		int status = eqSolver->solve(guess1, guess2, target, solution, toleranceSolved, iter_solved);

			if (tolSolved != nullptr && toleranceSolved != numeric_limits<double>::quiet_NaN()) {
				*tolSolved = toleranceSolved;
			}

			//spdlog::get("logger")->info(solverName + ". Result -> " + std::to_string(status) + ". Target = " + std::to_string(target) + ". Solution = " + std::to_string(solution)
			 //+ ". Tol = " + std::to_string(toleranceSolved));

			return status;
		/*}

		return C_monotonic_eq_solver::NO_SOLUTION;*/
	}

	~MonoSolver() {};
};
