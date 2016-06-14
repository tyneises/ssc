#include "lib_irradproc.h"

void calculate_sun_positions( size_t Nrec, 
	ssc_number_t *year, ssc_number_t *month, ssc_number_t *day, ssc_number_t *hour, ssc_number_t *minute,
	ssc_number_t lat, ssc_number_t lng, ssc_number_t tz,
	ssc_number_t *solazi, ssc_number_t *solzen )
{
	double sundata[9];
	for( size_t i=0;i<Nrec;i++ )
	{	
		// comments from irradproc:
		//   sunn[0] = azm = sun azimuth in radians, measured east from north, 0 to 2*pi
		//   sunn[1] = 0.5*pi - elv = sun zenith in radians, 0 to pi
		solarpos( (int)year[i], (int)month[i], (int)day[i], (int)hour[i], (int)minute[i], 
			lat, lng, tz, sundata );
		solazi[i] = sundata[0]*180/3.1415926;
		solzen[i] = sundata[1]*180/3.1415926;
	}
}