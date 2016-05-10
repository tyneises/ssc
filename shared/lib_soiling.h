
//turns out you can't include files from ssc in shared
//#include "core.h"
//in which case, where do we put the variable tables??

//calculates soiling for a single timestep
void apply_soiling_loss(double *irradiance, double soiling_loss)
{
	*irradiance *= (1 - soiling_loss / 100); //convert soiling loss from a percent
}


