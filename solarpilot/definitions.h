//File automatically generated by variable_map_makefile.py 2018/2/14


/*******************************************************************************************************
*  Copyright 2017 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "System Advisor Model" or "SAM". Except
*  to comply with the foregoing, the terms "System Advisor Model", "SAM", or any confusingly similar
*  designation may not be used to refer to any modified version of this software or any modified
*  version of the underlying software originally provided by Alliance without the prior written consent
*  of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#ifndef _VARDEFS_
#define _VARDEFS_ 1

#include "mod_base.h"
#include "Toolbox.h"

#ifdef _MSC_VER
#pragma warning(disable : 4267)
#endif

//Enumeration of data columns in the variable map file


//Custom module settings
#define _CUSTOM_REC 0		//If using custom geometry functions in the CustomReceiverWindow, define 1
//Sandbox mode
#define _SANDBOX 0
//demo only
#define _DEMO 0	//Is this a demo version?
const int _demo_date[] = {2014,8,1};
//Include Coretrace (relevant to fieldcore only! Disabling this option will cause SolarPILOT compilation to fail.).
#ifdef SP_STANDALONE
	#define SP_USE_SOLTRACE
	//Compile without threading functionality? Comment out to remove.
	#define SP_USE_THREADS
	//crete local make-dir functions
	#ifdef _WIN32 
	    #define SP_USE_MKDIR
	#endif
#endif

#ifndef PI
    #define PI 3.14159265358979311600
    #define R2D 57.29577951308232286465
    #define D2R 0.01745329251994329547
#endif



struct var_ambient
{
	spvar< matrix_t<double> > atm_coefs; 		//[none] Atmospheric attenuation coefficients for user-defined analysis
	spvar< std::string > atm_model; 		//[none] Atmospheric attenuation model {0=25km Barstow, 1 = 5km Barstow, 2 = user defined}
	struct ATM_MODEL{ enum EN{DELSOL3_CLEAR_DAY=0,DELSOL3_HAZY_DAY=1,USERDEFINED=2}; };
	spvar< std::string > class_name; 		//[none] Class name
	spvar< double > del_h2o; 		//[mm H2O] Atmospheric precipitable water depth for use in the Allen insolation model
	spvar< double > dni_layout; 		//[W/m2] DNI to use during all layout calculations. CONSTANT model only.
	spvar< double > dpres; 		//[atm] Local ambient pressure relative to sea-level pressure
	spvar< double > elevation; 		//[m] Plant mean elevation
	spvar< std::string > insol_type; 		//[none] Model used to determine insolation as a function of time
	struct INSOL_TYPE{ enum EN{WEATHER_FILE_DATA=-1,MEINEL_MODEL=0,HOTTEL_MODEL=1,CONSTANT_VALUE=2,ALLEN_MODEL=3,MOON_MODEL=4}; };
	spvar< double > latitude; 		//[deg] Plant latitude
	spvar< std::string > loc_city; 		//[none] City or place name for weather station (informational only)
	spvar< std::string > loc_state; 		//[none] State name for weather station (informational only)
	spvar< double > longitude; 		//[deg] Plant longitude
	spvar< double > sun_csr; 		//[none] Ratio of solar flux contained in the circumsolar ring over the solar disc flux
	spvar< matrix_t<double> > sun_pos_map; 		//[[deg,deg]] Map of sun positions to use for calculations
	spvar< double > sun_rad_limit; 		//[mrad] Half-angle of sunshape size (4.65mrad for Pillbox, 2.73mrad for Gaussian)
	spvar< std::string > sun_type; 		//[none] Sunshape model - {0=point sun, 1=limb darkened sun, 2=square wave sun, 3=user sun}
	struct SUN_TYPE{ enum EN{PILLBOX_SUN=2,GAUSSIAN_SUN=4,LIMBDARKENED_SUN=1,POINT_SUN=0,BUIE_CSR=5,USER_SUN=3}; };
	spvar< double > time_zone; 		//[hr] Time zone
	spvar< matrix_t<double> > user_sun; 		//[[deg,deg]] Array of intensity at various angles from the centroid of the sun
	spvar< std::string > weather_file; 		//[none] Weather file to use for analysis
	spvar< WeatherData > wf_data; 		//[none] Data entries in the weather file
	spout< double > atm_atten_est; 		//[%] Average solar field attenuation due to atmospheric scattering
	spout< double > sim_time_step; 		//[sec] Simulation weather data time step

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_financial
{
	spvar< std::string > class_name; 		//[none] Class name
	spvar< double > contingency_rate; 		//[%] Fraction of the direct capital costs added to account for contingency
	spvar< double > fixed_cost; 		//[$] Cost that does not scale with any plant parameter
	spvar< double > heliostat_spec_cost; 		//[$/m2] Cost per square meter of heliostat aperture area of the heliostat field
	spvar< bool > is_pmt_factors; 		//[] Enable or disable the use of weighting factors in determining field layout
	spvar< double > land_spec_cost; 		//[$/acre] Cost of land per acre including the footprint of the land occupied by the entire plant.
	spvar< double > plant_spec_cost; 		//[$/kWe] Cost of the power block and balance of plant equipment per kilowatt (electric) gross design power
	spvar< std::vector< double > > pmt_factors; 		//[none] Relative value of electricity produced during this period compared to the average
	spvar< double > rec_cost_exp; 		//[none] Exponent in the equation (total cost) = (ref. cost) * ( (area) / (ref. area) ) ^ X
	spvar< double > rec_ref_area; 		//[m2] Receiver surface area corresponding to the receiver reference cost
	spvar< double > rec_ref_cost; 		//[$] Cost of the receiver at the sizing indicated by the reference receiver area
	spvar< double > sales_tax_frac; 		//[%] Fraction of the direct capital costs for which sales tax applies
	spvar< double > sales_tax_rate; 		//[%] Sales tax rate applid to the total direct capital cost
	spvar< double > site_spec_cost; 		//[$/m2] Cost per square meter of heliostat aperture area of site improvements
	spvar< double > tes_spec_cost; 		//[$/kWht] Cost of thermal storage per kilowatt hour (thermal) capacity
	spvar< double > tower_exp; 		//[none] Exponent in the equation (total cost) = (fixed cost) * exp( X * (tower height) )
	spvar< double > tower_fixed_cost; 		//[$] Fixed tower cost - used as the basis for scaling tower cost as a function of height
	spvar< std::string > weekday_sched; 		//[] Weekday dispatch period schedule
	spvar< std::string > weekend_sched; 		//[] Weekend dispatch period schedule
	spvar< double > wiring_user_spec; 		//[$/m2] Cost of wiring per square meter of heliostat aperture area
	spout< double > contingency_cost; 		//[$] Contingency cost
	spout< double > cost_per_capacity; 		//[$/kWe] Estimated capital cost per capacity (net)
	spout< double > heliostat_cost; 		//[$] Heliostat field cost
	spout< double > land_cost; 		//[$] Land cost
	spout< double > plant_cost; 		//[$] Cost of the power block and balance of plant equipment
	spout< std::vector< double > > pricing_array; 		//[none] Yearly time series schedule of price multipliers to incentivize electricity sales at particular times
	spout< double > rec_cost; 		//[$] Receiver cost
	spout< double > sales_tax_cost; 		//[$] Sales tax cost
	spout< std::vector< int > > schedule_array; 		//[none] Yearly time series schedule of TOU periods
	spout< double > site_cost; 		//[$] Site improvements cost
	spout< double > tes_cost; 		//[$] Thermal storage cost
	spout< double > total_direct_cost; 		//[$] Sum of all direct costs
	spout< double > total_indirect_cost; 		//[$] Sum of all indirect costs
	spout< double > total_installed_cost; 		//[$] Sum of direct and indirect costs
	spout< double > tower_cost; 		//[$] Tower cost
	spout< double > wiring_cost; 		//[$] Wiring cost

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_fluxsim
{
	spvar< std::string > aim_method; 		//[] Method for determining the aim point for each heliostat
	struct AIM_METHOD{ enum EN{SIMPLE_AIM_POINTS=0,SIGMA_AIMING=1,PROBABILITY_SHIFT=2,IMAGE_SIZE_PRIORITY=3,KEEP_EXISTING=4,FREEZE_TRACKING=5}; };
	spvar< std::string > class_name; 		//[none] Class name
	spvar< double > cloud_depth; 		//[m] Depth of the cloud shape
	spvar< double > cloud_loc_x; 		//[m] Base location of the cloud(s) relative to the tower position - X dimension
	spvar< double > cloud_loc_y; 		//[m] Base location of the cloud(s) relative to the tower position - Y dimension
	spvar< double > cloud_opacity; 		//[none] Fraction of DNI obfuscated by a cloud shadow
	spvar< double > cloud_sep_depth; 		//[] Cloud pattern depth spacing
	spvar< double > cloud_sep_width; 		//[] Cloud pattern width spacing
	spvar< std::string > cloud_shape; 		//[] Shape used to model the cloud shadow
	struct CLOUD_SHAPE{ enum EN{ELLIPTICAL=0,RECTANGULAR=1,FRONT=2}; };
	spvar< double > cloud_skew; 		//[deg] Angle between North and the depth direction (-180 to +180 with clockwise positive)
	spvar< double > cloud_width; 		//[m] Width of the cloud shape
	spvar< std::string > flux_data; 		//[] 2D matrix of flux data
	spvar< int > flux_day; 		//[] Day of the month for the flux simulation
	spvar< std::string > flux_dist; 		//[] Sampling basis for random positioning. Non-uniform distributions are weighted away from the center.
	struct FLUX_DIST{ enum EN{TRIANGULAR=0,NORMAL=1,UNIFORM=2}; };
	spvar< double > flux_dni; 		//[W/m2] Direct Normal Irradiation at the specified simulation point
	spvar< double > flux_hour; 		//[hr] Hour of the day for the flux simulation
	spvar< std::string > flux_model; 		//[none] Desired flux simulation tool. Not all geometries can be simulated using the Hermite approximation.
	struct FLUX_MODEL{ enum EN{HERMITE_ANALYTICAL=0,SOLTRACE=1}; };
	spvar< int > flux_month; 		//[] Month of the year for the flux simulation
	spvar< double > flux_solar_az_in; 		//[] Solar azimuth angle to use for the flux simulation
	spvar< double > flux_solar_el_in; 		//[] Solar elevation angle to use for the flux simulation
	spvar< std::string > flux_time_type; 		//[none] Method for specifying the desired flux simulation time.
	struct FLUX_TIME_TYPE{ enum EN{SUN_POSITION=0,HOURDAY=1}; };
	spvar< bool > is_autoscale; 		//[none] Autoscale the Z-axis of the contour plot
	spvar< bool > is_cloud_pattern; 		//[] Create a pattern based on the specified cloud
	spvar< bool > is_cloud_symd; 		//[] Mirror the cloud pattern below the width axis
	spvar< bool > is_cloud_symw; 		//[] Mirror the cloud pattern to the left of the depth axis
	spvar< bool > is_cloudy; 		//[] Enable simulation for a cloud transient
	spvar< bool > is_load_raydata; 		//[none] Load heliostat field raytrace data from an existing file
	spvar< bool > is_optical_err; 		//[] Include the reflector optical error sources in the SolTrace simulation
	spvar< bool > is_save_raydata; 		//[none] Save heliostat field raytrace data to a file for future re-use
	spvar< bool > is_sunshape_err; 		//[] Include the sun shape error in the SolTrace simulation
	spvar< int > max_rays; 		//[none] The maximum number of generated rays allowed before terminating the simulation. Overrides the desired rays setting.
	spvar< int > min_rays; 		//[none] The minimum number of ray hits on the receiver before terminating the simulation.
	spvar< double > norm_dist_sigma; 		//[] Size of the standard distribution relative to half of the height of the receiver.
	spvar< double > plot_zmax; 		//[none] Z-axis maximum value
	spvar< double > plot_zmin; 		//[none] Z-axis minimum value
	spvar< std::string > raydata_file; 		//[none] Location and file of the ray data
	spvar< bool > save_data; 		//[] Save the results for each ray
	spvar< std::string > save_data_loc; 		//[] Choose a location to save the ray data
	spvar< int > seed; 		//[none] The seed for the random number generator
	spvar< double > sigma_limit_x; 		//[] Minimum distance (std. dev.) between optical center of heliostat image and the receiver edge in the receiver-X direction
	spvar< double > sigma_limit_y; 		//[none] Minimum distance (std. dev.) between optical center of heliostat image and the receiver edge in the receiver-Y direction
	spvar< int > x_res; 		//[none] Number of flux test points per panel (maximum) in the vertical direction for the flux simulation
	spvar< int > y_res; 		//[none] Number of flux test points per panel (maximum) in the horizontal direction for the flux simulation
	spout< double > flux_solar_az; 		//[deg] Solar azimuth angle to use for the flux simulation
	spout< double > flux_solar_el; 		//[deg] Solar elevation angle to use for the flux simulation

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_heliostat
{
	spvar< int > cant_day; 		//[day] Day of the year used for canting the heliostats (1-365)
	spvar< double > cant_hour; 		//[hr] Hours past noon at which the mirror panels are canted (-12 to 12)
	spvar< std::string > cant_method; 		//[none] Integer to specify the canting method {0=none, -1=Cant on-axis equal to slant range, 1=user-defined on-axis, 3=user-defined off-axis at hour + day}
	struct CANT_METHOD{ enum EN{NO_CANTING=0,ONAXIS_AT_SLANT=-1,ONAXIS_USERDEFINED=1,OFFAXIS_DAY_AND_HOUR=3,USERDEFINED_VECTOR=4}; };
	spvar< double > cant_rad_scaled; 		//[none] Canting radius value (absolute value if radius is not scaled, multiplied by tower height if scaled)
	spvar< double > cant_vect_i; 		//[none] Canting vector - x-component
	spvar< double > cant_vect_j; 		//[none] Canting vector y-component
	spvar< double > cant_vect_k; 		//[none] Canting vector z-component
	spvar< double > cant_vect_scale; 		//[m] Value to scale the canting unit vector to determine actual canting magnitude
	spvar< void* > cbdata; 		//[none] Data pointer for UI page
	spvar< std::string > class_name; 		//[none] Class name
	spvar< double > diameter; 		//[m] Diameter of the heliostat structure (round heliostats only)
	spvar< double > err_azimuth; 		//[rad] Standard deviation of the normal error dist. of the azimuth angle
	spvar< double > err_elevation; 		//[rad] Standard deviation of the normal error dist. of the elevation angle
	spvar< double > err_reflect_x; 		//[rad] error in reflected vector (horiz.) caused by atmospheric refraction, tower sway, etc.
	spvar< double > err_reflect_y; 		//[rad] error in reflected vector (vert.) caused by atmospheric refraction, tower sway, etc.
	spvar< double > err_surface_x; 		//[rad] Std.dev. of the normal error dist. of the reflective surface normal in the X (horizontal)
	spvar< double > err_surface_y; 		//[rad] Same as above, but in the vertical direction
	spvar< std::string > focus_method; 		//[none] The focusing method {0=Flat, 1=Each at slant, 2=Average of group, 3=User defined}
	struct FOCUS_METHOD{ enum EN{FLAT=0,AT_SLANT=1,GROUP_AVERAGE=2,USERDEFINED=3}; };
	spvar< double > height; 		//[m] Height of the heliostat structure
	spvar< std::string > helio_name; 		//[] Heliostat template name
	spvar< int > id; 		//[none] Unique ID number for the heliostat template
	spvar< bool > is_cant_rad_scaled; 		//[] The cant radius scales with tower height
	spvar< bool > is_cant_vect_slant; 		//[none] Multiply the canting vector by the slant range
	spvar< bool > is_enabled; 		//[] Is template enabled?
	spvar< bool > is_faceted; 		//[none] The number of reflective panels per heliostat is greater than 1
	spvar< bool > is_focal_equal; 		//[none] Both the X and Y focal lengths will use a single value as indicated by the X focal length
	spvar< std::string > is_round; 		//[none] Is the heliostat round (true) or rectangular (false)
	struct IS_ROUND{ enum EN{RECTANGULAR=0,ROUND=1}; };
	spvar< bool > is_xfocus; 		//[none] Reflector is focused in with respect to the heliostat X axis
	spvar< bool > is_yfocus; 		//[none] Reflector is focused in with respect to the heliostat Y axis
	spvar< int > n_cant_x; 		//[] Number of cant panels in the X direction
	spvar< int > n_cant_y; 		//[] Number of cant panels in the Y direction
	spvar< double > reflect_ratio; 		//[none] Ratio of mirror area to total area of the heliostat defined by wm x hm
	spvar< double > reflectivity; 		//[none] Average reflectivity (clean) of the mirrored surface
	spvar< double > rvel_max_x; 		//[rad/s] maximum rotational velocity about the x axis
	spvar< double > rvel_max_y; 		//[rad/s] maximum rotational velocity about the z axis
	spvar< double > soiling; 		//[none] Average soiling factor
	spvar< std::string > st_err_type; 		//[none] Error distribution of the reflected rays from the heliostat optical surface
	struct ST_ERR_TYPE{ enum EN{GAUSSIAN=0,PILLBOX=1}; };
	spvar< double > temp_az_max; 		//[deg] Angular boundary for heliostat geometry - on the clockwise side of the region
	spvar< double > temp_az_min; 		//[deg] Angular boundary for heliostat geometry - on the counter-clockwise side of the region
	spvar< double > temp_rad_max; 		//[none] Maximum radius at which this heliostat geometry can be used
	spvar< double > temp_rad_min; 		//[none] Minimum radius at which this heliostat geometry can be used
	spvar< int > template_order; 		//[none] template_order
	spvar< std::string > track_method; 		//[none] Specify how often heliostats update their tracking position 
	struct TRACK_METHOD{ enum EN{CONTINUOUS=0,PERIODIC=1}; };
	spvar< double > track_period; 		//[sec] The amount of time between tracking updates for each heliostat
	spvar< int > type; 		//[] Integer used to group heliostats into geometries within a field, (e.g. 5 different focal length designs)
	spvar< double > width; 		//[m] Width of the heliostat structure
	spvar< double > x_focal_length; 		//[m] Reflector focal length with respect to the heliostat X (horizontal) axis
	spvar< double > x_gap; 		//[m] Separation between panels in the horizontal direction
	spvar< double > y_focal_length; 		//[m] Reflector focal length with respect to the heliostat Y (vertical) axis
	spvar< double > y_gap; 		//[m] Separation between panels in the vertical direction
	spout< double > area; 		//[m2] Aperture area including geometry penalties and gaps in the structure
	spout< double > cant_mag_i; 		//[none] Total canting vector - i
	spout< double > cant_mag_j; 		//[none] Total canting vector - j
	spout< double > cant_mag_k; 		//[none] Total canting vector - k
	spout< double > cant_norm_i; 		//[none] Normalized canting vector - i
	spout< double > cant_norm_j; 		//[none] Normalized canting vector - j
	spout< double > cant_norm_k; 		//[none] Normalized canting vector - k
	spout< double > cant_radius; 		//[m] Radius for canting focal point assuming on-axis canting
	spout< double > cant_sun_az; 		//[deg] Sun azimuth angle at the moment the cant panels are focused on the receiver
	spout< double > cant_sun_el; 		//[deg] Sun elevation angle at the moment the cant panels are focused on the receiver
	spout< double > err_total; 		//[rad] Total convolved optical error in the reflected beam from the above sources
	spout< double > r_collision; 		//[m] Distance between heliostat center and maximum radial extent of structure
	spout< double > ref_total; 		//[none] Effective reflectance - product of the mirror reflectance and soiling

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_land
{
	spvar< std::string > class_name; 		//[none] Class name
	spvar< std::vector<std::vector<sp_point> > > exclusions; 		//[] Vector of arrays that specify the regions of land to exclude in the heliostat layout
	spvar< double > import_tower_lat; 		//[deg] Imported land boundary tower latitude
	spvar< double > import_tower_lon; 		//[deg] Imported land boundary tower longitude
	spvar< bool > import_tower_set; 		//[none] Has the tower location been set for imported land geometries?
	spvar< std::vector<std::vector<sp_point> > > inclusions; 		//[] Vector of arrays that specify the regions of land to include in the heliostat layout
	spvar< bool > is_bounds_array; 		//[none] Land boundary is specified by points array
	spvar< bool > is_bounds_fixed; 		//[none] Land boundary has fixed limits (not more than | not less than)
	spvar< bool > is_bounds_scaled; 		//[none] Land boundary scales with tower hight value
	spvar< bool > is_exclusions_relative; 		//[none] Shift the exclusion regions along with the tower offset values
	spvar< double > land_const; 		//[acre] Fixed land area that is added to the area occupied by heliostats
	spvar< double > land_mult; 		//[none] Factor multiplying the land area occupied by heliostats
	spvar< double > max_fixed_rad; 		//[m] Outer land boundary for circular land plot
	spvar< double > max_scaled_rad; 		//[none] Maximum radius (in units of tower height) for positioning of the heliostats
	spvar< double > min_fixed_rad; 		//[m] Inner land boundary for circular land plot
	spvar< double > min_scaled_rad; 		//[none] Minimum radius (in units of tower height) for positioning of the heliostats
	spvar< double > tower_offset_x; 		//[m] Displacement of the tower in X relative to the X-positions specified in the land table
	spvar< double > tower_offset_y; 		//[m] Displacement of the tower in Y relative to the Y-positions specified in the land table
	spout< double > bound_area; 		//[acre] Land area occupied by heliostats. This value is the area of a convex hull surrounding the heliostat positions.
	spout< double > land_area; 		//[acre] Land area, including solar field and multiplying factor
	spout< double > radmax_m; 		//[m] Calculated maximum distance between tower and last row of heliostats
	spout< double > radmin_m; 		//[m] Calculated minimum distance between tower and first row of heliostats

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_optimize
{
	spvar< std::string > algorithm; 		//[] Optimization algorithm
	struct ALGORITHM{ enum EN{BOBYQA=0,COBYLA=1,NEWOUA=2,NELDERMEAD=3,SUBPLEX=4,RSGS=5}; };
	spvar< std::string > class_name; 		//[] 
	spvar< double > converge_tol; 		//[none] Relative change in the objective function below which convergence is achieved
	spvar< double > flux_penalty; 		//[none] Relative weight in the objective function given to flux intensity over the allowable limit
	spvar< int > max_desc_iter; 		//[none] Maximum number of steps along the direction of steepest descent before recalculating the response surface
	spvar< int > max_gs_iter; 		//[none] Maximum number of golden section iterations to refine the position of a local minimum
	spvar< int > max_iter; 		//[none] Maximum number of times the optimization can iterate
	spvar< double > max_step; 		//[none] Maximum total relative step size during optimization
	spvar< double > power_penalty; 		//[none] Relative weight in the objective function given to power to the receiver below the required minimum
	spout< double > aspect_display; 		//[none] Current receiver aspect ratio (H/W)
	spout< double > gs_refine_ratio; 		//[none] The relative step size of the refined area during refinement simulations. More iterations will allow greater refinement

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_parametric
{
	spvar< std::string > class_name; 		//[none] Class name
	spvar< std::string > eff_file_name; 		//[] Name of the output file containing the efficiency matrix
	spvar< std::string > flux_file_name; 		//[] Name of the output file containing the fluxmap data
	spvar< std::string > fluxmap_format; 		//[] Dimensions of the fluxmap data (rows,cols)
	struct FLUXMAP_FORMAT{ enum EN{SAM_FORMAT=0,N12X10_ARRAY=1,SPECIFIED_DIMENSIONS=2}; };
	spvar< bool > is_fluxmap_norm; 		//[] Flux data is reported as normalized
	spvar< bool > par_save_field_img; 		//[none] Save field efficiency image
	spvar< bool > par_save_flux_dat; 		//[none] Save receiver flux data
	spvar< bool > par_save_flux_img; 		//[none] Save receiver flux image
	spvar< bool > par_save_helio; 		//[none] Save detailed heliostat performance data for each run
	spvar< bool > par_save_summary; 		//[none] Save detailed system performance data to a file for each run
	spvar< std::string > sam_grid_format; 		//[none] SAM data grid format
	struct SAM_GRID_FORMAT{ enum EN{AUTO_SPACING=0,EVEN_GRID=1}; };
	spvar< std::string > sam_out_dir; 		//[] Output directory
	spvar< bool > upar_save_field_img; 		//[none] Save field efficiency image
	spvar< bool > upar_save_flux_dat; 		//[none] Save receiver flux data
	spvar< bool > upar_save_flux_img; 		//[none] Save receiver flux image
	spvar< bool > upar_save_helio; 		//[none] Save detailed heliostat performance data for each run
	spvar< bool > upar_save_summary; 		//[none] Save detailed system performance data to a file for each run
	spvar< std::string > user_par_values; 		//[none] User parametric values

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_plant
{
	spvar< std::string > class_name; 		//[none] Class name
	spvar< double > eta_cycle; 		//[none] Thermodynamic efficiency of the power cycle, including feedwater pumps and cooling equipment parasitics
	spvar< double > hours_tes; 		//[hr] Capacity of Hours of thermal storage operation at full cycle output
	spvar< double > par_factor; 		//[none] Estimated ratio of net power output to gross power output at design
	spvar< double > solar_mult; 		//[none] Ratio of thermal power output from the solar field to power cycle thermal input at design
	spout< double > power_gross; 		//[MWe] Rated nameplate design gross turbine electric output, not accounting for parasitic losses
	spout< double > power_net; 		//[MWe] Estimated net electric power at design, accounting for all parasitic losses

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_receiver
{
	spvar< double > absorptance; 		//[none] Energy absorbed by the receiver surface before accounting for radiation/convection losses
	spvar< std::string > accept_ang_type; 		//[none] Receiver angular acceptance window defines angles about the aperture normal, can be rectangular or elliptical shape
	struct ACCEPT_ANG_TYPE{ enum EN{RECTANGULAR=0,ELLIPTICAL=1}; };
	spvar< double > accept_ang_x; 		//[deg] Acceptance angle of the receiver in the horizontal direction (in aperture coordinates)
	spvar< double > accept_ang_y; 		//[deg] Acceptance angle of the receiver in the vertical direction (in aperture coordinates)
	spvar< std::string > aperture_type; 		//[] The shape of the receiver aperture
	struct APERTURE_TYPE{ enum EN{RECTANGULAR=0}; };
	spvar< void* > cbdata; 		//[none] Data pointer for UI page
	spvar< std::string > class_name; 		//[none] Class name
	spvar< int > id; 		//[] Template ID
	spvar< bool > is_aspect_opt; 		//[] Optimize receiver aspect ratio (height / width)
	spvar< bool > is_enabled; 		//[] Is template enabled?
	spvar< bool > is_open_geom; 		//[] If true, the receiver is represented by an arc rather than a closed circle/polygon
	spvar< bool > is_polygon; 		//[] Receiver geometry is represented as discrete polygon of N panels rather than continuous arc
	spvar< int > n_panels; 		//[none] Number of receiver panels (polygon facets) for a polygonal receiver geometry
	spvar< double > panel_rotation; 		//[deg] Azimuth angle between the normal vector to the primary 'north' panel and North
	spvar< double > peak_flux; 		//[kW/m2] Maximum allowable flux intensity on any portion of the receiver surface
	spvar< double > piping_loss_coef; 		//[kW/m] Loss per meter of tower height
	spvar< double > piping_loss_const; 		//[kW] Constant thermal loss due to piping - doesn't scale with tower height
	spvar< double > rec_azimuth; 		//[deg] Receiver azimuth orientation: 0 deg is north, positive clockwise
	spvar< double > rec_cav_cdepth; 		//[m] Offset of centroid of cavity absorber surface from the aperture plane. (Positive->Increased depth)
	spvar< double > rec_cav_rad; 		//[m] Radius of the receiver cavity absorbing surface
	spvar< double > rec_diameter; 		//[m] Receiver diameter for cylindrical receivers
	spvar< double > rec_elevation; 		//[deg] Receiver elevation orientation: 0 deg to the horizon, negative rotating downward
	spvar< double > rec_height; 		//[m] Height of the absorbing component
	spvar< std::string > rec_name; 		//[] Receiver template name
	spvar< double > rec_offset_x; 		//[m] Offset of receiver center in the East(+)/West(-) direction from the tower
	spvar< double > rec_offset_y; 		//[m] Offset of receiver center in the North(+)/South(-) direction from the tower
	spvar< double > rec_offset_z; 		//[m] Offset of the receiver center in the vertical direction, positive upwards
	spvar< std::string > rec_type; 		//[none] Receiver geometrical configuration
	struct REC_TYPE{ enum EN{EXTERNAL_CYLINDRICAL=0,FLAT_PLATE=2}; };
	spvar< double > rec_width; 		//[m] Receiver width for cavity or flat receivers
	spvar< double > span_max; 		//[deg] Maximum (CW) bound of the arc defining the receiver surface
	spvar< double > span_min; 		//[deg] Minimum (CCW) bound of the arc defining the receiver surface
	spvar< double > therm_loss_base; 		//[kW/m2] Thermal loss from the receiver at design-point conditions
	spvar< matrix_t<double> > therm_loss_load; 		//[none] Temperature-dependant thermal loss
	spvar< matrix_t<double> > therm_loss_wind; 		//[none] Wind speed-dependant thermal loss
	spout< double > absorber_area; 		//[m2] Effective area of the receiver absorber panels
	spout< double > optical_height; 		//[m] Calculated height of the centerline of the receiver above the plane of the heliostats
	spout< double > piping_loss; 		//[MW] Thermal loss from non-absorber receiver piping
	spout< double > rec_aspect; 		//[none] Ratio of receiver height to width
	spout< double > therm_eff; 		//[none] Receiver calculated thermal efficiency
	spout< double > therm_loss; 		//[MW] Receiver thermal loss at design

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_solarfield
{
	spvar< double > accept_max; 		//[deg] Upper bound of the angular range containing the heliostat field
	spvar< double > accept_min; 		//[deg] Lower bound of the angular range containing the heliostat field
	spvar< double > az_spacing; 		//[none] Azimuthal spacing factor for the first row of heliostats after a reset. Heliostats separated by heliostat width times this factor.
	spvar< std::string > class_name; 		//[none] Class name
	spvar< std::string > des_sim_detail; 		//[none] Simulation detail for placing heliostats (see definitions in options spreadsheet)
	struct DES_SIM_DETAIL{ enum EN{SUBSET_OF_DAYSHOURS=2,SINGLE_SIMULATION_POINT=1,DO_NOT_FILTER_HELIOSTATS=0,ANNUAL_SIMULATION=3,LIMITED_ANNUAL_SIMULATION=4,REPRESENTATIVE_PROFILES=5,EFFICIENCY_MAP__ANNUAL=6}; };
	spvar< int > des_sim_ndays; 		//[none] For limited annual simulation, the number of evenly spaced days to simulate
	spvar< int > des_sim_nhours; 		//[none] Simulation will run with the specified hourly frequency (1=every hour, 2=every other hour...)
	spvar< double > dni_des; 		//[W/m2] DNI value at which the design-point receiver thermal power is achieved
	spvar< std::string > hsort_method; 		//[none] Select the criteria by which heliostats will be included in the solar field layout.
	struct HSORT_METHOD{ enum EN{POWER_TO_RECEIVER=0,TOTAL_EFFICIENCY=1,COSINE_EFFICIENCY=2,ATTENUATION_EFFICIENCY=3,INTERCEPT_EFFICIENCY=4,BLOCKING_EFFICIENCY=5,SHADOWING_EFFICIENCY=6,TOUWEIGHTED_POWER=7}; };
	spvar< bool > is_opt_zoning; 		//[none] Enables grouping of heliostats into zones for intercept factor calculation during layout only
	spvar< bool > is_prox_filter; 		//[none] Post-process the layout to select heliostats that are closer to the tower.
	spvar< bool > is_sliprow_skipped; 		//[none] Radial gap before first row after slip plane is sufficient to eliminate blocking
	spvar< bool > is_tht_opt; 		//[none] Vary the tower height during optimization to identify optimal level?
	spvar< std::string > layout_data; 		//[] Layout data in string form
	spvar< std::string > layout_method; 		//[none] Field layout method
	struct LAYOUT_METHOD{ enum EN{RADIAL_STAGGER=1,CORNFIELD=2,USERDEFINED=3}; };
	spvar< double > max_zone_size_az; 		//[tower-ht] Maximum zone size (azimuthal direction) for grouping optical intercept factor calculations
	spvar< double > max_zone_size_rad; 		//[tower-ht] Maximum zone size (radial direction) for grouping optical intercept factor calculations
	spvar< double > min_zone_size_az; 		//[tower-ht] Minimum zone size (azimuthal direction) for grouping optical intercept factor calculations
	spvar< double > min_zone_size_rad; 		//[tower-ht] Minimum zone size (radial direction) for grouping optical intercept factor calculations
	spvar< double > prox_filter_frac; 		//[none] Fraction of heliostats to subject to proximity filter.
	spvar< double > q_des; 		//[MWt] Design thermal power delivered from the solar field
	spvar< std::string > rad_spacing_method; 		//[none] Method for determining radial spacing during field layout for radial-stagger
	struct RAD_SPACING_METHOD{ enum EN{NO_BLOCKINGDENSE=3,ELIMINATE_BLOCKING=2,DELSOL_EMPIRICAL_FIT=1}; };
	spvar< double > row_spacing_x; 		//[none] Separation between adjacent heliostats in the X-direction, multiplies heliostat radius
	spvar< double > row_spacing_y; 		//[none] Separation between adjacent heliostats in the Y-direction, multiplies heliostat radius
	spvar< double > shadow_height; 		//[m] Effective tower height for shadowing calculations
	spvar< double > shadow_width; 		//[m] Effective tower diameter for shadowing calculations
	spvar< double > slip_plane_blocking; 		//[none] Allowable blocking in slip plane
	spvar< double > spacing_reset; 		//[none] For heliostat layout - ratio of maximum to initial azimuthal spacing before starting new compressed row
	spvar< double > sun_az_des_user; 		//[deg] Solar azimuth angle at the design point
	spvar< double > sun_el_des_user; 		//[deg] Solar elevation angle at the design point
	spvar< std::string > sun_loc_des; 		//[none] Sun location when thermal power rating is achieved
	struct SUN_LOC_DES{ enum EN{SUMMER_SOLSTICE=0,EQUINOX=1,WINTER_SOLSTICE=2,ZENITH=3,OTHER=4}; };
	spvar< std::string > temp_which; 		//[none] Select the heliostat geometry template that will be used in the layout
	struct TEMP_WHICH{ enum EN{}; };
	spvar< std::string > template_rule; 		//[] Method for distributing heliostat geometry templates in the field
	struct TEMPLATE_RULE{ enum EN{USE_SINGLE_TEMPLATE=0,SPECIFIED_RANGE=1,EVEN_RADIAL_DISTRIBUTION=2}; };
	spvar< double > tht; 		//[m] Average height of the tower receiver centerline above the base heliostat pivot point elevation
	spvar< double > trans_limit_fact; 		//[none] Determines the point at which close-packing switches to standard layout. =1 at no-blocking transition limit.
	spvar< std::string > xy_field_shape; 		//[] Enforced shape of the heliostat field
	struct XY_FIELD_SHAPE{ enum EN{HEXAGON=0,RECTANGLE=1,UNDEFINED=2}; };
	spvar< double > xy_rect_aspect; 		//[none] Aspect ratio of the rectangular field layout (height in Y / width in X)
	spvar< double > zone_div_tol; 		//[none] Allowable variation in optical intercept factor within a layout zone
	spout< double > rec_area; 		//[m2] Surface area from all receivers included in the solar field
	spout< double > sf_area; 		//[m2] The sum of all heliostat reflector area in the current layout
	spout< WeatherData > sim_step_data; 		//[none] Data used for design simulations
	spout< double > sun_az_des; 		//[deg] Calculated design-point solar azimuth
	spout< double > sun_el_des; 		//[deg] Calculated design-point solar elevation

    void addptrs(unordered_map<std::string, spbase*> &pmap);
};


struct var_map
{
    var_ambient amb;
    var_financial fin;
    var_fluxsim flux;
    var_land land;
    var_optimize opt;
    var_parametric par;
    var_plant plt;
    var_solarfield sf;
    
    std::vector< var_heliostat > hels;
    std::vector< var_receiver > recs;
   
    var_map();
    var_map( var_map &vc );   
    void copy( var_map &vc );
    void reset();

    unordered_map<std::string, spbase*> _varptrs;
    void add_receiver(int id);
    void add_heliostat(int id);
    void drop_receiver(int id);
    void drop_heliostat(int id);
};

#endif
