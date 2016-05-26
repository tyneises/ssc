
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

void apply_ac_loss(int nrec, int current_year, ssc_number_t *p_gen, util::matrix_t<ssc_number_t> &acpwr_gross, double ac_loss_percent)
{
	// p_gen is lifetime (nrec*nyears) and p_acpwr_gross is for current year (nrec)
	// assumes current year is 0 based
	size_t begin = current_year * nrec; 
	for (size_t i = 0; i != nrec; i++)
		p_gen[i + begin] = (ssc_number_t)((acpwr_gross.at(i) - (fabs(acpwr_gross.at(i)) * ac_loss_percent / 100)) * 0.001);
}