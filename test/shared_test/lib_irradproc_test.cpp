#include <stdlib.h>
#include <gtest/gtest.h>

#include "lib_irradproc.h"

using namespace std;

/**
 * \class IrradTest
 *
 * Month: 1-12, Hour: 0-23, Minute: 0-59.
 *
 */

class IrradTest : public ::testing::Test{
protected:
	double lat, lon, tz, alb, tilt, azim, rotlim, gcr;
	int year, month, day, skymodel, tracking;
	bool backtrack_on;
	double calc_sunrise, calc_sunset;
	double e;

	void SetUp(){
		// parameters
		lat = 31.6340;
		lon = 74.8723;
		tz = 5.5;
		year = 2017;
		month = 7;
		day = 19;
		skymodel = 2;
		alb = 0.2;
		tracking = 0;
		tilt = 10;
		azim = 180;
		rotlim = 0;
		backtrack_on = false;
		gcr = 0;
		e = 0.0001;

		// correct sunrise and sunset times
		calc_sunrise = 5.70924; // 5:43 am
		calc_sunset = 19.5179;  // 7:31 pm
	}
};

class NightCaseIrradProc : public IrradTest{
protected:
	// Test time: 1:30 am
	irrad irr_hourly_night;
	// Test time: 1:15 am
	irrad irr_15m_night;

	void SetUp(){
		IrradTest::SetUp();
		int night_hr(1);
		irr_hourly_night.set_time(year, month, day, night_hr, 30, 1);
		irr_hourly_night.set_location(lat, lon, tz);
		irr_hourly_night.set_sky_model(skymodel, alb);
		irr_hourly_night.set_beam_diffuse(0, 0);
		irr_hourly_night.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
		irr_15m_night.set_time(year, month, day, night_hr, 15, -1);
		irr_15m_night.set_location(lat, lon, tz);
		irr_15m_night.set_sky_model(skymodel, alb);
		irr_15m_night.set_beam_diffuse(0, 0);
		irr_15m_night.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
	}
};

class SunriseCaseIrradProc : public IrradTest{
protected:
	// Test time: 5:30 am
	irrad irr_hourly_sunrise;
	// Test time: 5:30 am
	irrad irr_15m_sunrise;

	void SetUp(){
		IrradTest::SetUp();
		int sr_hr(5);
		irr_hourly_sunrise.set_time(year, month, day, sr_hr, 30, 1);
		irr_hourly_sunrise.set_location(lat, lon, tz);
		irr_hourly_sunrise.set_sky_model(skymodel, alb);
		irr_hourly_sunrise.set_beam_diffuse(0, 1);
		irr_hourly_sunrise.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
		irr_15m_sunrise.set_time(year, month, day, sr_hr, 30, 1);
		irr_15m_sunrise.set_location(lat, lon, tz);
		irr_15m_sunrise.set_sky_model(skymodel, alb);
		irr_15m_sunrise.set_beam_diffuse(0, 1);
		irr_15m_sunrise.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
	}
};

class DayCaseIrradProc : public IrradTest{
protected:
	// Test time: 12:30 pm
	irrad irr_hourly_day;
	// Test time: 12:45 pm
	irrad irr_15m_day;
	
	void SetUp(){
		IrradTest::SetUp();
		int day_hr(12);
		irr_hourly_day.set_time(year, month, day, day_hr, 30, 1);
		irr_hourly_day.set_location(lat, lon, tz);
		irr_hourly_day.set_sky_model(skymodel, alb);
		irr_hourly_day.set_beam_diffuse(2, 2);
		irr_hourly_day.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
		irr_15m_day.set_time(year, month, day, day_hr, 45, 1);
		irr_15m_day.set_location(lat, lon, tz);
		irr_15m_day.set_sky_model(skymodel, alb);
		irr_15m_day.set_beam_diffuse(2, 2);
		irr_15m_day.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
	}
};

class SunsetCaseIrradProc : public IrradTest{
protected:
	// Test time: 7:30 pm
	irrad irr_hourly_sunset;
	// Test time: 7:30 pm
	irrad irr_15m_sunset;
	
	virtual void SetUp(){
		IrradTest::SetUp();
		int ss_hr(19);
		irr_hourly_sunset.set_time(year, month, day, ss_hr, 30, 1);
		irr_hourly_sunset.set_location(lat, lon, tz);
		irr_hourly_sunset.set_sky_model(skymodel, alb);
		irr_hourly_sunset.set_beam_diffuse(0, 1);
		irr_hourly_sunset.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
		irr_15m_sunset.set_time(year, month, day, ss_hr, 30, 1);
		irr_15m_sunset.set_location(lat, lon, tz);
		irr_15m_sunset.set_sky_model(skymodel, alb);
		irr_15m_sunset.set_beam_diffuse(0, 1);
		irr_15m_sunset.set_surface(tracking, tilt, azim, rotlim, backtrack_on, gcr);
	}
};

/**
 * Solar Position Function Tests
 * Output: sun[] = azimuth (rad), zenith(rad), elevation(rad), declination(rad), sunrise time, sunset time, 
 * eccentricity correction factor, true solar time, extraterrestrial solar irradiance on horizontal (W/m2) 
 */

TEST_F(NightCaseIrradProc, solarposTest_lib_irradproc){
	
	double sun[9];
	vector<double> sunrise_times;
	vector<double> sunset_times;

	/* Just before sunrise test case */
	solarpos(year, month, day, 4, 30, lat, lon, tz, sun);
	vector<double> solution = { 0.95662, 1.79457, -0.223771, 0.363938, 5.70882, 19.5183, 0.968276, 3.88646, 0 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "hourly before-sunrise case, parameter " << i << " fail\n";
	}
	solarpos(year, month, day, 5, 15, lat, lon, tz, sun);
	solution = { 1.0744, 1.65255, -0.0817513, 0.363839, 5.7091, 19.518, 0.96828, 4.63642, 0 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "15m before-sunrise case, parameter " << i << " fail\n";
	}

	/* Just after sunset test case */
	solarpos(year, month, day, 20, 30, lat, lon, tz, sun);
	solution = { 5.28748, 1.75391, -0.183117, 0.361807, 5.71544, 19.5131, 0.968361, 19.8857, 0 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "hourly after-sunset case, parameter " << i << " fail\n";
	}
	solarpos(year, month, day, 19, 45, lat, lon, tz, sun);
	solution = { 5.17431, 1.60864, -0.0378397, 0.361908, 5.71513, 19.5133, 0.968357, 19.1358, 0 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "15m after-sunrise case, parameter " << i << " fail\n";
	}
}

TEST_F(SunriseCaseIrradProc, solarposTest_lib_irradproc){
	double sun[9];
	vector<double> sunrise_times;
	vector<double> sunset_times;

	solarpos(year, month, day, 5, 30, lat, lon, tz, sun);
	vector<double> solution = { 1.11047, 1.6031, -0.0323028, 0.363806, 5.70924, 19.5179, 0.968281, 4.88641, 0 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "sunrise case, parameter " << i << " fail\n";
	}
}

TEST_F(DayCaseIrradProc, solarposTest_lib_irradproc){
	double sun[9];
	vector<double> sunrise_times;
	vector<double> sunset_times;

	/* Just before sunset test case */
	solarpos(year, month, day, 18, 30, lat, lon, tz, sun);
	vector<double>solution = { 5.01022, 1.3584, 0.212397, 0.362076, 5.71461, 19.5137, 0.96835, 17.8858, 279.08756 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "hourly before-sunset case, parameter " << i << " fail\n";
	}
	solarpos(year, month, day, 19, 15, lat, lon, tz, sun);
	solution = { 5.10579, 1.51295, 0.0578472, 0.361975, 5.71492, 19.5135, 0.968354, 18.6358, 76.5423 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "15m before-sunset case, parameter " << i << " fail\n";
	}

	/* Sunset time test case */
	solarpos(year, month, day, 19, 30, lat, lon, tz, sun);
	solution = { 5.13947, 1.55886, 0.0119379, 0.361941, 5.71503, 19.5134, 0.968356, 18.8858, 15.8044 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "sunset case, parameter " << i << " fail\n";
	}
}

TEST_F(SunsetCaseIrradProc, solarposTest_lib_irradproc){
	double sun[9];
	vector<double> sunrise_times;
	vector<double> sunset_times;

	/* Sunset time test case */
	solarpos(year, month, day, 19, 30, lat, lon, tz, sun);
	vector<double>solution = { 5.13947, 1.55886, 0.0119379, 0.361941, 5.71503, 19.5134, 0.968356, 18.8858, 15.8044 };
	sunrise_times.push_back(solution[4]);
	sunset_times.push_back(solution[5]);
	for (int i = 0; i < 9; i++){
		EXPECT_NEAR((double)sun[i], solution[i], e) << "sunset case, parameter " << i << " fail\n";
	}
}

/**
* Solar Incidence Function Test
* Mode = 0 for fixed tilt.
* Output: angle[] = incident angle (rad), tilt angle (rad), surface azimuth (rad), tracking axis rotation angle for single axis tracker (rad),
* backtracking angle difference: rot - ideal_rot (rad)
*/

TEST_F(NightCaseIrradProc, incidenceTest_lib_irradproc){
	int mode = 0;
	double angle[5] = { 0 };
	double sun_zen, sun_azm;
	vector<double> solutions;

	/* Just before sunrise test case */
	sun_azm = 0.95662;
	sun_zen = 1.79457;
	incidence(mode, tilt, azim, rotlim, sun_zen, sun_azm, backtrack_on, gcr, angle);
	solutions = { 1.89243, 0.174533, 3.14159, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle[i], solutions[i], e) << "before-sunrise case";
	}
}

TEST_F(SunriseCaseIrradProc, incidenceTest_lib_irradproc){
	int mode = 0;
	double angle[5] = { 0 };
	double sun_zen, sun_azm;
	double solution;
	vector<double> solutions;

	sun_azm = 1.11047;
	sun_zen = 1.6031;
	incidence(mode, tilt, azim, rotlim, sun_zen, sun_azm, backtrack_on, gcr, angle);
	solution = 1.67992;
	EXPECT_NEAR(angle[0], solution, e) << "sunrise case";
}

TEST_F(DayCaseIrradProc, incidenceTest_lib_irradproc){
	int mode = 0;
	double angle[5] = { 0 };
	double sun_zen, sun_azm;
	double solution;
	vector<double> solutions;

	sun_azm = 0;
	sun_zen = 0;
	incidence(mode, tilt, azim, rotlim, sun_zen, sun_azm, backtrack_on, gcr, angle);
	solution = 0.174533;
	EXPECT_NEAR(angle[0], solution, e) << "noon case";
}

TEST_F(SunsetCaseIrradProc, incidenceTest_lib_irradproc){
	int mode = 0;
	double angle[5] = { 0 };
	double sun_zen, sun_azm;
	double solution;
	vector<double> solutions;

	sun_azm = 5.13947;
	sun_zen = 1.55886;
	incidence(mode, tilt, azim, rotlim, sun_zen, sun_azm, backtrack_on, gcr, angle);
	solution = 1.631;
	EXPECT_NEAR(angle[0], solution, e) << "sunset case";
}

/**
* Calc Function Tests
* Output:
* sun[] =	azimuth (rad), zenith(rad), elevation(rad), declination(rad), sunrise time, sunset time,
*			eccentricity correction factor, true solar time, extraterrestrial solar irradiance on horizontal (W/m2);
* angle_p[] = incident angle (rad), tilt angle (rad), surface azimuth (rad), tracking axis rotation angle for single axis tracker (rad),
*			backtracking angle difference: rot - ideal_rot (rad);
* poa_p[] = incident beam, incident sky diffuse, incident ground diffuse, diffuse isotropic, diffuse circumsolar, horizon brightening (W/m2);
* irrad parameters: ghi, dni, dhi
*/

TEST_F(NightCaseIrradProc, CalcTestRadMode0_lib_irradproc){
	vector<double> sun_p;
	sun_p.resize(10);
	int sunup = false;
	vector<double> angle_p;
	angle_p.resize(5);
	vector<double> poa_p;
	poa_p.resize(6);
	vector<double> rad_p = { 1, 1, 1 };

	irr_hourly_night.set_beam_diffuse(rad_p[1], rad_p[2]);
	irr_15m_night.set_beam_diffuse(rad_p[1], rad_p[2]);

	/* Hourly during the night */
	irr_hourly_night.calc();
	irr_hourly_night.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_hourly_night.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_hourly_night.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_hourly_night.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	vector<double> sun_solution = { -999, -999, -999, 20.795182, 5.711921, 19.515852, 0, 0.968315, 11.386113, 1286.786711 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "hourly_night, sun parameter " << i << " fail\n";
	}
	vector<double> angle_solution = { 0, 0, 0, 0, 0 };	// azim & tilt returned as 0 when sun is down
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "hourly_night, angle parameter " << i << " fail\n";
	}
	vector<double> poa_solution = { 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "hourly_night, poa parameter " << i << " fail\n";
	}
	vector<double>  rad_solution = { 0, 0, 0 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "hourly_night, irradiance parameter " << i << " fail\n";
	}

	/* 15m during the night */
	irr_15m_night.calc();
	irr_15m_night.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_15m_night.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_15m_night.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_15m_night.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	sun_solution = { -999, -999, -999, 20.795182, 5.711921, 19.515852, 0, 0.968315, 11.386113, 1286.786711 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "15m_night, sun parameter " << i << " fail\n";
	}
	angle_solution = { 0, 0, 0, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "15m_night, angle parameter " << i << " fail\n";
	}
	poa_solution = { 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "15m_night, poa parameter " << i << " fail\n";
	}
	rad_solution = { 0, 0, 0 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "15m_night, irradiance parameter " << i << " fail\n";
	}
}

TEST_F(SunriseCaseIrradProc, CalcTestRadMode0_lib_irradproc){
	vector<double> sun_p;
	sun_p.resize(10);
	int sunup = true;
	vector<double> angle_p;
	angle_p.resize(5);
	vector<double> poa_p;
	poa_p.resize(6);
	vector<double> rad_p = { 1, 1, 1 };

	irr_hourly_sunrise.set_beam_diffuse(rad_p[1], rad_p[2]);
	irr_15m_sunrise.set_beam_diffuse(rad_p[1], rad_p[2]);

	/* hourly during sunrise */
	irr_hourly_sunrise.calc();
	irr_hourly_sunrise.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_hourly_sunrise.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_hourly_sunrise.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_hourly_sunrise.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	vector<double> sun_solution = { 66.441256, 87.969757, 2.030243, 20.841844, 5.709384, 19.517827, 2.0, 0.968283, 5.242355, 46.902536 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "hourly_sunrise, sun parameter " << i << " fail\n";
	}
	vector<double> angle_solution = { 91.975545, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "hourly_sunrise, angle parameter " << i << " fail\n";
	}
	vector<double> poa_solution = { 0, 0.992404, 0, 0.992404, 0, 0 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "hourly_sunrise, poa parameter " << i << " fail\n";
	}
	vector<double> rad_solution = { -999, 1, 1 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "hourly_sunrise, irradiance parameter " << i << " fail\n";
	}

	/* 15m during sunrise */
	irr_15m_sunrise.calc();
	irr_15m_sunrise.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_15m_sunrise.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_15m_sunrise.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_15m_sunrise.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	sun_solution = { 66.441256, 87.969757, 2.030243, 20.841844, 5.709384, 19.517827, 2.0, 0.968283, 5.242355, 46.902536 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "15m_sunrise, sun parameter " << i << " fail\n";
	}
	angle_solution = { 91.975545, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "15m_sunrise, angle parameter " << i << " fail\n";
	}
	poa_solution = { 0, 0.992404, 0, 0.992404, 0, 0 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "15m_sunrise, poa parameter " << i << " fail\n";
	}
	rad_solution = { -999, 1, 1 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "15m_sunrise, irradiance parameter " << i << " fail\n";
	}
}

TEST_F(DayCaseIrradProc, CalcTestRadMode0_lib_irradproc){
	vector<double> sun_p;
	sun_p.resize(10);
	int sunup = true;
	vector<double> angle_p;
	angle_p.resize(5);
	vector<double> poa_p;
	poa_p.resize(6);
	vector<double> rad_p = { 1, 1, 1 };

	irr_hourly_day.set_beam_diffuse(rad_p[1], rad_p[2]);
	irr_15m_day.set_beam_diffuse(rad_p[1], rad_p[2]);

	/* Hourly during the day */
	irr_hourly_day.calc();
	irr_hourly_day.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_hourly_day.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_hourly_day.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_hourly_day.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	vector<double> sun_solution = { 171.563258, 10.938523, 79.061477, 20.791368, 5.712128, 19.515691, 1.0, 0.968318, 11.886091, 1299.866650 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "hourly_day, sun parameter " << i << " fail\n";
	}
	vector<double> angle_solution = { 1.795054, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "hourly_day, angle parameter " << i << " fail\n";
	}
	vector<double> poa_solution = { 0.999509, 1.052130, 0.003011, 0.194810, 0.818170, 0.039150 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "hourly_day, poa parameter " << i << " fail\n";
	}
	vector<double> rad_solution = { -999, 1.0, 1.0 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "hourly_day, irradiance parameter " << i << " fail\n";
	}

	/* 15m during the day */
	irr_15m_day.calc();
	irr_15m_day.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_15m_day.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_15m_day.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_15m_day.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	sun_solution = { 190.054756, 10.986005, 79.013994, 20.789459, 5.712231, 19.515691, 1.000000, 0.968318, 12.136079, 1299.658007 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "15m_day, sun parameter " << i << " fail\n";
	}
	angle_solution = { 2.075962, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "15m_day, angle parameter " << i << " fail\n";
	}
	poa_solution = { 0.999343, 1.052130, 0.003011, 0.195107, 0.817860, 0.039150 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "15m_day, poa parameter " << i << " fail\n";
	}
	rad_solution = { -999, 1.0, 1.0 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "15m_day, irradiance parameter " << i << " fail\n";
	}
}

TEST_F(SunsetCaseIrradProc, CalcTestRadMode0_lib_irradproc){
	vector<double> sun_p;
	sun_p.resize(10);
	int sunup = false;
	vector<double> angle_p;
	angle_p.resize(5);
	vector<double> poa_p;
	poa_p.resize(6); 
	vector<double> rad_p = { 1, 1, 1 };
	irr_hourly_sunset.set_beam_diffuse(rad_p[1], rad_p[2]);
	irr_15m_sunset.set_beam_diffuse(rad_p[1], rad_p[2]);
	
	/* hourly during sunset */
	irr_hourly_sunset.calc();
	irr_hourly_sunset.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_hourly_sunset.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_hourly_sunset.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_hourly_sunset.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	vector<double>sun_solution = { 292.600441, 86.774381, 3.225619, 20.739568, 5.714928, 19.513485, 3.0, 0.968355, 18.643720, 74.494280 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "hourly_sunset, sun parameter " << i << " fail\n";
	}
	vector<double>angle_solution = { 90.642562, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "hourly_sunset, angle parameter " << i << " fail\n";
	}
	vector<double>poa_solution = { 0, 0.981644, 0.001605, 0.992404, 0, -0.010760 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "hourly_sunset, poa parameter " << i << " fail\n";
	}
	vector<double>rad_solution = { -999, 1, 1 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "hourly_sunset, irradiance parameter " << i << " fail\n";
	}

	/* 15m during sunset */
	irr_15m_sunset.calc();
	irr_15m_sunset.get_sun(&sun_p[0], &sun_p[1], &sun_p[2], &sun_p[3], &sun_p[4], &sun_p[5], &sunup, &sun_p[7], &sun_p[8], &sun_p[9]);
	irr_15m_sunset.get_angles(&angle_p[0], &angle_p[1], &angle_p[2], &angle_p[3], &angle_p[4]);
	irr_15m_sunset.get_poa(&poa_p[0], &poa_p[1], &poa_p[2], &poa_p[3], &poa_p[4], &poa_p[5]);
	irr_15m_sunset.get_irrad(&rad_p[0], &rad_p[1], &rad_p[2]);

	sun_p[6] = (double)sunup;
	sun_solution = { 292.600441, 86.774381, 3.225619, 20.739568, 5.714928, 19.513485, 3.0, 0.968355, 18.643720, 74.494280 };
	for (int i = 0; i < 10; i++){
		EXPECT_NEAR(sun_p[i], sun_solution[i], e) << "15m_sunset, sun parameter " << i << " fail\n";
	}
	angle_solution = { 90.642562, tilt, azim, 0, 0 };
	for (int i = 0; i < 5; i++){
		EXPECT_NEAR(angle_p[i], angle_solution[i], e) << "15m_sunset, angle parameter " << i << " fail\n";
	}
	poa_solution = { 0, 0.981644, 0.001605, 0.992404, 0, -0.010760 };
	for (int i = 0; i < 6; i++){
		EXPECT_NEAR(poa_p[i], poa_solution[i], e) << "15m_sunset, poa parameter " << i << " fail\n";
	}
	rad_solution = { -999, 1, 1 };
	for (int i = 0; i < 3; i++){
		EXPECT_NEAR(rad_p[i], rad_solution[i], e) << "15m_sunset, irradiance parameter " << i << " fail\n";
	}

	/*
	printf("sun:%f, %f, %f, %f, %f, %f, %f, %f, %f, %f", sun_p[0], sun_p[1], sun_p[2], sun_p[3], sun_p[4], sun_p[5], (double)sunup, sun_p[6], sun_p[7], sun_p[8]);
	printf("angles: %f, %f, %f, %f, %f \n", angle_p[0], angle_p[1], angle_p[2], angle_p[3], angle_p[4]);
	printf("poa: %f, %f, %f, %f, %f, %f \n", poa_p[0], poa_p[1], poa_p[2], poa_p[3], poa_p[4], poa_p[5]);
	printf("irrad: %f, %f, %f \n", &rad_p[0], &rad_p[1], &rad_p[2]);
	*/
}