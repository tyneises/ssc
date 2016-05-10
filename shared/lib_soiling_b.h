
//turns out you can't include files from ssc in shared
//#include "core.h"
//in which case, where do we put the variable tables??

//typedef std::vector<double> doublevec; //fancy hints from Nick

//applies soiling to the entire array of irradiance based on an array of soiling losses, and months.
void apply_soiling_loss_b(double nrec, std::vector<double> month, ssc_number_t *soiling_loss, std::vector<double> *irradiance)
{
	for (int i = 0; i < nrec; i++)
	{
		int m = month[i];
		(*irradiance)[i] *= (1 - soiling_loss[m] / 100);
	}

}


