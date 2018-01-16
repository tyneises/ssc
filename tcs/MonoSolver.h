#pragma once
#include "numeric_solvers.h"
#include "MonoChain.h"

class MonoSolver
{
private:
	
	double target;

	double* guessValues;

	double lowerBound;

	double upperBound;

	double tolerance;

	int iterationLimit;

	bool isErrorRel;

	C_monotonic_equation* equation;

	C_monotonic_eq_solver* eqSolver;

	int (MonoSolver::*nextSolver)();

	void setSettings();
	
public:
	MonoSolver(double target, double* guessValues, double lowerBound, double upperBound, double tolerance,
		int iterationLimit, bool isErrorRel, MonoChain* chain, int (MonoSolver::*nextSolver)(), int (MonoChain::*model)(double x, double* y));

	int solve();

	~MonoSolver();
};

