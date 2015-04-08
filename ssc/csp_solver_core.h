#ifndef __csp_solver_core_
#define __csp_solver_core_

#include "tcstype.h"
#include "core.h"

#include "lib_weatherfile.h"

class C_csp_component : public compute_module
{
private:
	const tcsvarinfo *params;
	std::vector<tcsvalue> param_values;

	void set_csp_component_value(const char *name, double val);
	void set_csp_component_value(const char *name, double *p_array, size_t l_array);
	void set_csp_component_value(const char *name, double *p_array, size_t nr, size_t nc);
	void set_csp_component_value(const char *name, const char *s);

	void set_csp_component_value(int idx, double val);
	void set_csp_component_value(int idx, double *p_array, size_t l_array);
	void set_csp_component_value(int idx, double *p_array, size_t nr, size_t nc);
	void set_csp_component_value(int idx, const char *s);

	int find_var(const char *name);

	void tcsvalue_free(int idx);
	void tcsvalue_set_number(int idx, double val);
	void tcsvalue_set_array(int idx, double *p_array, size_t l_array);
	void tcsvalue_set_matrix(int idx, double *p_array, size_t nr, size_t nc);
	void tcsvalue_set_string(int idx, const char *s);	

protected:
	void set_params_and_size_vector(const tcsvarinfo *params_in);
	double value(size_t idx);
	double *value(size_t idx, int *len);
	double *value(size_t idx, int *nr, int *nc);
	void value(size_t idx, double val);
	std::string value_str(size_t idx);
	double *allocate(size_t idx, int nr, int nc);

	void tcsmatrix_index(size_t idx, int row, int col, double val);
	double tcsmatrix_index(size_t idx, int row, int col);

	// Have to have this to access 'compute_module' public methods??
	virtual void exec() throw(general_error) = 0;
	

public:
	virtual ~C_csp_component(){};

	C_csp_component(){};

	void set_csp_component_value_ssc_double(const char *name, double val);
	void set_csp_component_value_ssc_array(const char *name, ssc_number_t *p_array, size_t *l_array_in);
	void set_csp_component_value_ssc_matrix(const char *name, ssc_number_t *p_array, size_t *nr_in, size_t *nc_in);
	void set_csp_component_value_ssc_string(const char *name, const char *s);

	virtual void init() {return;}

	virtual void converged() {return;}
	
};


class C_csp_weatherreader : public C_csp_component
{
private:
	weatherfile m_wf;
	bool m_first;		// flag to indicate whether this is the first call

	virtual void exec() throw(general_error)
	{

	}

protected:

public:

	C_csp_weatherreader();

	~C_csp_weatherreader(){};

	virtual void init();

	virtual void converged() {return;}

	void timestep_call(double time, double step);
};



class C_csp_collector_receiver : public C_csp_component
{
private:

protected:

	virtual void exec() throw(general_error) = 0;

public:
	C_csp_collector_receiver(){};

	~C_csp_collector_receiver(){};

	virtual void init() = 0;	// pure virtual function

	virtual void converged() = 0;

	// solve_field

	// optical_efficiency

	// internal_energy || time_to_startup

};

class C_csp_thermal_storage : public C_csp_component
{
private:

protected:

	virtual void exec() throw(general_error) = 0;

public:
	C_csp_thermal_storage(){};

	~C_csp_thermal_storage(){};

	virtual void init() = 0;		// pure virtual function

	virtual void converged() = 0;

	virtual double estimate_charging_availability() = 0;

	virtual double estimate_discharging_availability() = 0;
};

class C_csp_power_cycle : public C_csp_component
{
private:

protected:

	virtual void exec() throw(general_error) = 0;

public:

	C_csp_power_cycle(){};

	~C_csp_power_cycle(){};

	virtual void init() = 0;

	virtual void converged() = 0;

};




class C_csp_solver : public compute_module
{
private:
	C_csp_weatherreader *mp_weatherreader;
	C_csp_collector_receiver *mp_collector_receiver;
	C_csp_thermal_storage *mp_thermal_storage;
	C_csp_power_cycle *mp_power_cycle;

	virtual void exec() throw(general_error)
	{

	}


public:
	C_csp_solver();

	~C_csp_solver(){};

	void setup_technology_model(C_csp_weatherreader *p_weatherreader, C_csp_collector_receiver *p_collector_receiver,
		C_csp_thermal_storage *p_thermal_storage, C_csp_power_cycle *p_power_cycle);

	void timeseries_simulation();
};



#endif // !__csp_solver_core_