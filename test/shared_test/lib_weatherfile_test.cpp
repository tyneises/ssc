#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "lib_weatherfile.h"
#include "common.h"
#include "vartab.h"

using namespace std;

/**
* \class weatherfileTest
*
* For testing input as files. Possible files are SAMCSV, EPY, TMY2, TMY3...
*
*/

class weatherfileTest : public ::testing::Test{
protected:
	weatherfile wf;
	string file;
	double e;		//epsilon for double comparison

	virtual void SetUp(){}
};

class CSVCase_WeatherfileTest : public weatherfileTest{
protected:
	void SetUp(){
		e = 0.001;
#ifdef _MSC_VER	
		file = "../../../test/input_docs/weather-noRHum.csv";
#else	
		file = "../test/input_docs/weather-noRHum.csv";
#endif	
		ASSERT_TRUE(wf.open(file));
	}
};

/// Test some init actions
TEST_F(CSVCase_WeatherfileTest, initTest){

	EXPECT_EQ(wf.header().location, "875760") << "CSV Case: Init test\n";
	EXPECT_EQ(wf.header().city, "Buenos_Aires") << "CSV Case: Init test\n";
	EXPECT_EQ("", wf.message()) << "CSV Case: Init test\n";
	EXPECT_EQ(wf.type(), 5) << "CSV Case: Init test\n";
	EXPECT_FALSE(wf.ok()) << "CSV Case: Init test\n";
	wf.rewind();
	EXPECT_EQ(wf.get_counter_value(), 0) << "CSV Case: Init test\n";
	EXPECT_EQ(wf.start_sec(), 1800) << "CSV Case: Init test\n";
	EXPECT_EQ(wf.step_sec(), 3600) << "CSV Case: Init test\n";
	EXPECT_EQ(wf.nrecords(), 8760) << "CSV Case: Init test\n";
	EXPECT_TRUE(wf.has_data_column(0));
	EXPECT_FALSE(wf.has_data_column(4));
	EXPECT_FALSE(wf.has_calculated_data(10)) << "Twet not able to be calculated\n";
}

TEST_F(CSVCase_WeatherfileTest, normalizeCityTest_lib_weatherfile){
	EXPECT_EQ("Buenos Aires", wf.normalize_city("buenos aires"));
}

/// Test reading first, second and custom row
TEST_F(CSVCase_WeatherfileTest, readTest){
	weather_record r;
	/* read first row */
	wf.read(&r);
	EXPECT_EQ(r.year, 1988) << "CSV Case: 1st row\n";
	EXPECT_EQ(r.month, 1) << "CSV Case: 1st row\n";
	EXPECT_EQ(r.day, 1) << "CSV Case: 1st row\n";
	EXPECT_EQ(r.hour, 0) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.minute, 30, e) << "CSV Case: 1st row\n";
	EXPECT_TRUE(isnan(r.gh)) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.dn, 0, e) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.df, 0, e) << "CSV Case: 1st row\n";
	EXPECT_TRUE(isnan(r.poa)) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.wspd, 2.1, e) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.wdir, 20, e) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.tdry, 20.9, e) << "CSV Case: 1st row\n";
	EXPECT_TRUE(isnan(r.twet)) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.tdew, 19.3, e) << "CSV Case: 1st row\n";
	EXPECT_TRUE(isnan(r.rhum)) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.pres, 1010, e) << "CSV Case: 1st row\n";
	EXPECT_TRUE(isnan(r.snow)) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.alb, 0.17, e) << "CSV Case: 1st row\n";
	EXPECT_NEAR(r.aod, 0.291, e) << "CSV Case: 1st row\n";
	EXPECT_EQ(wf.get_counter_value(), 1);

	/* read second row */
	wf.read(&r);
	EXPECT_EQ(r.year, 1988) << "CSV Case: 2nd row\n";
	EXPECT_EQ(r.month, 1) << "CSV Case: 2nd row\n";
	EXPECT_EQ(r.day, 1) << "CSV Case: 2nd row\n";
	EXPECT_EQ(r.hour, 1) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.minute, 30, e) << "CSV Case: 2nd row\n";
	EXPECT_TRUE(isnan(r.gh)) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.dn, 0, e) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.df, 0, e) << "CSV Case: 2nd row\n";
	EXPECT_TRUE(isnan(r.poa)) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.wspd, 1.5, e) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.wdir, 360, e) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.tdry, 20.9, e) << "CSV Case: 2nd row\n";
	EXPECT_TRUE(isnan(r.twet)) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.tdew, 19.4, e) << "CSV Case: 2nd row\n";
	EXPECT_TRUE(isnan(r.rhum)) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.pres, 1007, e) << "CSV Case: 2nd row\n";
	EXPECT_TRUE(isnan(r.snow)) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.alb, 0.17, e) << "CSV Case: 2nd row\n";
	EXPECT_NEAR(r.aod, 0.291, e) << "CSV Case: 2nd row\n";
	EXPECT_EQ(wf.get_counter_value(), 2);

	/* setting counter to another step */
	wf.set_counter_to(0);
	wf.read(&r);
	EXPECT_EQ(r.year, 1988) << "CSV Case: Reset to 1st row\n";
	EXPECT_EQ(r.month, 1) << "CSV Case: Reset to 1st row\n";
	EXPECT_EQ(r.day, 1) << "CSV Case: Reset to 1st row\n";
	EXPECT_EQ(r.hour, 0) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.minute, 30, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_TRUE(isnan(r.gh)) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.dn, 0, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.df, 0, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_TRUE(isnan(r.poa)) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.wspd, 2.1, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.wdir, 20, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.tdry, 20.9, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_TRUE(isnan(r.twet)) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.tdew, 19.3, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_TRUE(isnan(r.rhum)) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.pres, 1010, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_TRUE(isnan(r.snow)) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.alb, 0.17, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_NEAR(r.aod, 0.291, e) << "CSV Case: Reset to 1st row\n";
	EXPECT_EQ(wf.get_counter_value(), 1);
}

/**
* \class weatherdataTest
*
* For testing input as var_data type. Base class contains input var_data which contains a var_table.
* After child SetUp() modifies var_data that go into the var_table, base SetUp() called to assign as input.
*
*/

class weatherdataTest : public ::testing::Test{
protected: 
	var_table* vt;
	var_data vd, time;
	var_data* input;
	double e;		//epsilon for double comparison

	void SetUp(){
		e = 0.001;
		vt->assign("lat", 1);
		vt->assign("lon", 2);
		vt->assign("tz", 3);
		vt->assign("elev", 4);
		vt->assign("year", 5);
		vt->assign("month", time);
		vt->assign("day", time);
		vt->assign("hour", time);
		vt->assign("minute", vd);
		vt->assign("dn", vd);
		vt->assign("df", vd);
		vt->assign("wspd", vd);
		vt->assign("wdir", vd);
		vt->assign("tdry", vd);
		vt->assign("tdew", vd);
		vt->assign("alb", vd);
		vt->assign("aod", vd);

		input = new var_data;
		input->type = SSC_TABLE;
		input->table = *vt;
	}
};

/// Hourly data Case
class Data8760CaseWeatherData : public weatherdataTest{
protected:
	void SetUp(){
		// set-up input var_table: vd contains only zeros and time will be 1, 2, 3 
		vt = new var_table;
		float empty[8760] = { 0 };
		vd = var_data(empty, 8760);
		float order[3] = { 1, 2, 3 };
		time = var_data(order, 3);
		weatherdataTest::SetUp();
		
	}
};


TEST_F(Data8760CaseWeatherData, initTest_lib_weatherfile){
	weatherdata wd(input);
	EXPECT_FALSE(wd.has_message()) << "Error message was found:" << wd.message(); 
	EXPECT_EQ(wd.header().lat, 1) << "Latitude?";
	EXPECT_EQ(wd.header().lon, 2) << "Longitude?";
	EXPECT_EQ(wd.get_counter_value(), 0) << "Counter at beginning";
	EXPECT_EQ(wd.nrecords(), 8760) << "Number of records?";
	EXPECT_FALSE(wd.has_data_column(0)) << "Should not have Year column";
	EXPECT_TRUE(wd.has_data_column(4)) << "Should have Minute column";
	EXPECT_TRUE(wd.has_data_column(6)) << "Should have DN column";
	EXPECT_FALSE(wd.has_data_column(8)) << "Should not have POA column";
	EXPECT_FALSE(wd.has_calculated_data(10)) << "Twet not able to be calculated\n";
}

TEST_F(Data8760CaseWeatherData, readTest_lib_weatherfile){
	weatherdata wd(input);
	weather_record r;
	EXPECT_TRUE(wd.read(&r));
	EXPECT_EQ(r.year, 2000) << "Data8760 Case: 1st row\n";
	EXPECT_EQ(r.month, 1) << "Data8760 Case: 1st row\n";
	EXPECT_EQ(r.day, 1) << "Data8760 Case: 1st row\n";
	EXPECT_EQ(r.hour, 1) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.minute, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.gh)) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.dn, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.df, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.poa)) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.wspd, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.wdir, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.tdry, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.twet)) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.tdew, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.rhum)) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.pres)) << "Data8760 Case: 1st row\n";
	EXPECT_TRUE(isnan(r.snow)) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.alb, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_NEAR(r.aod, 0, e) << "Data8760 Case: 1st row\n";
	EXPECT_EQ(wd.get_counter_value(), 1);

	wd.set_counter_to(2);
	EXPECT_TRUE(wd.read(&r));
	EXPECT_EQ(r.year, 2000) << "Data8760 Case: 3rd row\n";
	EXPECT_EQ(r.month, 3) << "Data8760 Case: 3rd row\n";
	EXPECT_EQ(r.day, 3) << "Data8760 Case: 3rd row\n";
	EXPECT_EQ(r.hour, 3) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.minute, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.gh)) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.dn, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.df, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.poa)) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.wspd, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.wdir, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.tdry, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.twet)) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.tdew, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.rhum)) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.pres)) << "Data8760 Case: 3rd row\n";
	EXPECT_TRUE(isnan(r.snow)) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.alb, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_NEAR(r.aod, 0, e) << "Data8760 Case: 3rd row\n";
	EXPECT_EQ(wd.get_counter_value(), 3);
	// poa is not an assigned var_data input, but r.poa returns 0 whereas pres, rhum & twet
	// are not assigned but are NULL
}

/// Error Case
class Data9999CaseWeatherData : public weatherdataTest{
protected:
	void SetUp(){
		// set-up input var_table
		vt = new var_table;
		float empty[9999] = { 0 };
		vd = var_data(empty, 9999);
		float order[3] = { 1, 2, 3 };
		time = var_data(order, 3);
		weatherdataTest::SetUp();
	}
};

TEST_F(Data9999CaseWeatherData, initTest_lib_weatherfile){
	weatherdata wd(input);
	EXPECT_EQ(wd.nrecords(), 9999);
	string error = "could not determine timestep in weatherdata";
	EXPECT_EQ(error, wd.message()) << "Timestep should be invalid";
}

TEST_F(Data9999CaseWeatherData, initTest2_lib_weatherfile){
	input->table.unassign("gh");
	input->table.unassign("dn");
	input->table.unassign("df");
	weatherdata wd(input);
	string error = "missing irradiance: could not find gh, dn, df, or poa";
	EXPECT_EQ(error, wd.message()) << "No irradiance provided error";
}

TEST_F(Data9999CaseWeatherData, readTest2_lib_weatherfile){
	float wrong_length[1000] = { 0 };
	var_data vd_err = var_data(wrong_length, 1000);
	input->table.unassign("dn");
	input->table.assign("dn", vd_err);
	weatherdata wd(input);
	string error = "aod number of entries doesn't match with other fields";
	EXPECT_EQ(error, wd.message()) << "Irradiance entry length mismatch";
}