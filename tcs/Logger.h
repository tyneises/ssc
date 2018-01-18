#include <iostream>
#include "spdlog\spdlog.h"

using namespace std;
namespace spd = spdlog;

class Logger
{
private:
	Logger();

	static const string LOG_PATH;

	static shared_ptr<spd::logger> logger;

public:

	static shared_ptr<spd::logger> get();

	~Logger();
};
