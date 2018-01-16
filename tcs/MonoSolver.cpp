#include "MonoSolver.h"



void MonoSolver::setSettings()
{
	eqSolver->settings(tolerance, iterationLimit, lowerBound, upperBound, isErrorRel);
}

MonoSolver::MonoSolver(double target, double* guessValues, double lowerBound, double upperBound, double tolerance,
								int iterationLimit, bool isErrorRel, MonoChain* chain, int (MonoSolver::*nextSolver)(), int (MonoChain::*model)(double x, double* y))
{
	this->target = target;
	this->guessValues = guessValues;
	this->lowerBound = lowerBound;
	this->upperBound = upperBound;
	this->tolerance = tolerance;
	this->iterationLimit = iterationLimit;
	this->isErrorRel = isErrorRel;
	this->nextSolver = nextSolver;

	equation = new C_member_mono_eq<MonoChain>(chain, model);
	eqSolver = new C_monotonic_eq_solver(*equation);
	setSettings();
}

int MonoSolver::solve()
{
	double L_solved, tol_solved;
	L_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
	int iter_solved = -1;

	double y1, y2;

	double guess1 = guessValues[0];
	double guess2 = guessValues[1];
	int guess1_Result, guess2_Result;
	int i = 0;

	while (i + 1 < sizeof(*guessValues)){
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

		return eqSolver->solve(pair1, pair2, target, L_solved, tol_solved, iter_solved);
	}

}


MonoSolver::~MonoSolver()
{
}