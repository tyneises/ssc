#include "Logger.h"

const string Logger::LOG_PATH = "./RegenHX_LOG.log";
shared_ptr<spd::logger> Logger::logger = nullptr;

Logger::Logger()
{
}


shared_ptr<spd::logger> Logger::get()
{
	if (logger == nullptr) {
		logger = spd::basic_logger_mt("logger", LOG_PATH);
		spd::drop("logger");
		spd::register_logger(logger);
		spd::set_pattern("[%b %d %H:%M:%S] [%l] %v");
	}

	return logger;
}

Logger::~Logger()
{
}
