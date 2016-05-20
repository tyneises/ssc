
//turns out you can't include files from ssc in shared
//#include "core.h"
//in which case, where do we put the variable tables??

//typedef std::vector<double> doublevec; //fancy hints from Nick

//applies soiling to the entire array of irradiance based on an array of soiling losses, and months.
void apply_soiling_loss_b(int nrec, ssc_number_t *month, ssc_number_t *soiling_loss, std::vector<double> *irradiance)
{
	for (int i = 0; i < nrec; i++)
	{
		int m = month[i];
		(*irradiance)[i] *= (1 - soiling_loss[m-1] / 100);
	}

}

void apply_ac_loss(int current_year, ssc_number_t *p_gen, double acpwr_gross, double ac_loss_percent)
{
	// would also need to know the timestep
	int begin = current_year * 8760; 
	int end = (current_year + 1) * 8760;
	for (int i = begin; i != end; i++)
		p_gen[i] = (ssc_number_t)((acpwr_gross - (fabs(acpwr_gross) * ac_loss_percent / 100)) * 0.001);
}