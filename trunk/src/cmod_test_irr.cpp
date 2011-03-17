#include "lib_financial.h"
#include "core.h"
#include <sstream>

static var_info _cm_vtab_test_irr[] = {

/*test and validation of irr function used in ssc financial models */
/*   VARTYPE           DATATYPE         NAME                         LABEL                              UNITS     META                      GROUP          REQUIRED_IF                 CONSTRAINTS                      UI_HINTS*/
	{ SSC_INPUT,        SSC_NUMBER,      "analysis_years",			"years in test",	"years",   "",                      "DHF",             "*",						   "INTEGER",                              "" },
	{ SSC_OUTPUT,        SSC_NUMBER,      "cf_length",			"length of cash flow items",	"",   "",                      "DHF",             "*",						   "INTEGER",                              "" },

	{ SSC_INPUT,        SSC_NUMBER,      "max_iterations",			"maximum number of iterations",	"",   "",                      "DHF",             "*",						   "INTEGER",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "tolerance",			"convergence tolerance",	"",   "",                      "DHF",             "*",						   "",                              "" },
	{ SSC_INPUT,        SSC_NUMBER,      "initial_guess",			"initial guess",	"",   "",                      "DHF",             "*",						   "",                              "" },
	{ SSC_INPUT,        SSC_ARRAY,      "cf_test",		"cash flow input values over which to test irr",	"$",   "",  "DHF",             "*",						   "",                              "" },


	//{ SSC_OUTPUT,        SSC_ARRAY,       "cf_test_scaled",            "scaled cash flow input",                     "",      "",                      "DHF",             "*",                      "LENGTH_EQUAL=cf_length",                             "" },
	{ SSC_OUTPUT,        SSC_ARRAY,       "cf_irr",            "calculated irr",                     "%",      "",                      "DHF",             "*",                      "LENGTH_EQUAL=cf_length",                             "" },


var_info_invalid };

extern var_info
	vtab_standard_financial[],
	vtab_oandm[],
	vtab_tax_credits[],
	vtab_payment_incentives[];

enum {
	CF_test,
	CF_irr,

	CF_max };



class cm_test_irr : public compute_module
{
private:
	util::matrix_t<double> cf;

public:
	cm_test_irr()
	{
		add_var_info( _cm_vtab_test_irr );
	}

	void exec( ) throw( general_error )
	{

		// cash flow initialization
		int nyears = as_integer("analysis_years");
		int cf_length = nyears+1;
		int max_iterations = as_integer("max_iterations");
		double tolerance = as_number("tolerance");
		double initial_guess = as_number("initial_guess");
		cf.resize_fill( CF_max, nyears+1, 0.0 );
		int i;
		// initialize energy and revenue
		size_t count = 0;
		ssc_number_t *arrp = 0;
		
		arrp = as_array("cf_test", &count);
		i=0;
		while ( i <= nyears && i < (int)count )
		{
			cf.at(CF_test, i) = arrp[i];
			i++;
		}

		for (i=1;i<=nyears;i++) 
		{
			cf.at(CF_irr,i) = irr(CF_test,i,initial_guess,tolerance,max_iterations)*100;
		}


	    assign("cf_length", var_data((ssc_number_t) cf_length ));
		save_cf( CF_irr, nyears, "cf_irr" );
	}


	void save_cf(int cf_line, int nyears, const std::string &name)
	{
		ssc_number_t *arrp = allocate( name, nyears+1 );
		for (int i=0;i<=nyears;i++)
			arrp[i] = (ssc_number_t)cf.at(cf_line, i);
	}



/* ported from http://code.google.com/p/irr-newtonraphson-calculator/ */
	bool is_valid_iter_bound(double estimated_return_rate)
	{
		return estimated_return_rate != -1 && (estimated_return_rate < INT_MAX) && (estimated_return_rate > INT_MIN);
	}

	double irr_poly_sum(double estimated_return_rate, int cf_line, int count)
	{
		double sum_of_polynomial = 0;
		if (is_valid_iter_bound(estimated_return_rate))
		{
			for (int j = 0; j <= count ; j++)
			{
				double val = (pow((1 + estimated_return_rate), j));
				if (val != 0.0)
					sum_of_polynomial += cf.at(cf_line,j)/val;
				else
					break;
			}
		}
		return sum_of_polynomial;
	}

	double irr_derivative_sum(double estimated_return_rate,int cf_line, int count)
	{
		double sum_of_derivative = 0;
		if (is_valid_iter_bound(estimated_return_rate))
			for (int i = 1; i <= count ; i++)
			{
				sum_of_derivative += cf.at(cf_line,i)*(i)/pow((1 + estimated_return_rate), i+1);
			}
		return sum_of_derivative*-1;
	}

	double irr_scale_factor( int cf_unscaled, int count)
	{
		// scale to max value for better irr convergence
		if (count<1) return 1.0;
		int i=0;
		double max=fabs(cf.at(cf_unscaled,0));
		for (i=0;i<=count;i++) 
			if (fabs(cf.at(cf_unscaled,i))> max) max =fabs(cf.at(cf_unscaled,i));
		return (max>0 ? max:1);
	}

	double irr( int cf_line, int count, double initial_guess=-2, double tolerance=1e-5, int max_iterations=200 )
	{
		int number_of_iterations=0;
		double calculated_irr=0;


		if (count < 1)
			return calculated_irr;

		// only possible for first value negative
		if ( (cf.at(cf_line,0) <= 0))
		{
			// initial guess from http://zainco.blogspot.com/2008/08/internal-rate-of-return-using-newton.html
			if ((initial_guess < -1) && (count > 1))// second order
			{
				if (cf.at(cf_line,0) !=0) 
				{
					double b = 2.0+ cf.at(cf_line,1)/cf.at(cf_line,0);
					double c = 1.0+cf.at(cf_line,1)/cf.at(cf_line,0)+cf.at(cf_line,2)/cf.at(cf_line,0);
					initial_guess = -0.5*b - 0.5*sqrt(b*b-4.0*c);
					if ((initial_guess <= 0) || (initial_guess >= 1)) initial_guess = -0.5*b + 0.5*sqrt(b*b-4.0*c);
				}
			}
			else if (initial_guess < 0) // first order
			{
				if (cf.at(cf_line,0) !=0) initial_guess = -(1.0 + cf.at(cf_line,1)/cf.at(cf_line,0));
			}
			if ((initial_guess <= 0) || (initial_guess >= 1)) initial_guess = 0.1;

			double deriv_sum = irr_derivative_sum(initial_guess,cf_line,count);
			if (deriv_sum != 0.0)
				calculated_irr = initial_guess - irr_poly_sum(initial_guess,cf_line,count)/deriv_sum;
			else
				return initial_guess;

			number_of_iterations++;

			double scale_factor = irr_scale_factor(cf_line,count);

			double residual = irr_poly_sum(calculated_irr,cf_line,count) / scale_factor;

			while (!(fabs(residual) <= tolerance) && (number_of_iterations < max_iterations))
			{
				deriv_sum = irr_derivative_sum(initial_guess,cf_line,count);
				if (deriv_sum != 0.0)
					calculated_irr = calculated_irr - irr_poly_sum(calculated_irr,cf_line,count)/deriv_sum;
				else
					break;

			//	std::stringstream outm;
			//	outm << "iteration=" << number_of_iterations << ", residual="  << residual << ", deriv_sum=" << deriv_sum  << ", calculated_irr=" << calculated_irr ;
			//	log( outm.str() );

				number_of_iterations++;
				residual = irr_poly_sum(calculated_irr,cf_line,count) / scale_factor;
			}
			//std::stringstream outm;
			//outm << "initial_guess=" << initial_guess << " iterations=" << number_of_iterations << " irr=" << calculated_irr;
			//log( outm.str() );

		}
		return calculated_irr;
	}


	double min( double a, double b )
	{
		return (a < b) ? a : b;
	}

	double max( double a, double b )
	{
		return (a > b) ? a : b;
	}

};




DEFINE_MODULE_ENTRY( test_irr, "Test of IRR function in ssc_", 1 );


