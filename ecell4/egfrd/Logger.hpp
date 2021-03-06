#ifndef ECELL4_EGFRD_LOGGER_HPP
#define ECELL4_EGFRD_LOGGER_HPP

#include <cstdarg>
#include <set>
#include <vector>
#include <string>
#include <memory>

namespace ecell4
{
namespace egfrd
{

class LogAppender;
class LoggerManager;
class LoggerManagerRegistry;

class Logger
{
public:
    enum level
    {
        L_OFF = 0,
        L_DEBUG = 1,
        L_INFO = 2,
        L_WARNING = 3,
        L_ERROR = 4,
        L_FATAL = 5
    };

public:
    ~Logger();

    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;


    LoggerManager const& logging_manager() const;

    void level(enum level level);

    enum level level() const;

    char const* name() const
    {
        return name_.c_str();
    }

    std::shared_ptr<LoggerManager> manager() const;

    void debug(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_DEBUG, format, ap);
        va_end(ap);
    }

    void info(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_INFO, format, ap);
        va_end(ap);
    }

    void warn(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_WARNING, format, ap);
        va_end(ap);
    }

    void error(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_ERROR, format, ap);
        va_end(ap);
    }

    void fatal(char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(L_FATAL, format, ap);
        va_end(ap);
    }

    void log(enum level lv, char const* format, ...)
    {
        va_list ap;
        va_start(ap, format);
        logv(lv, format, ap);
        va_end(ap);
    }

    void logv(enum level lv, char const* format, va_list ap);

    void flush();

    Logger(LoggerManagerRegistry const& registry, char const* name);

    static Logger& get_logger(char const* name);

    static char const* stringize_error_level(enum level lv);

private:
    void ensure_initialized();

protected:
    LoggerManagerRegistry const& registry_; 
    std::string const name_;
    std::shared_ptr<LoggerManager> manager_;
    enum level level_;
    std::vector<std::shared_ptr<LogAppender> > appenders_;
};

class LoggerManager
{
    friend class Logger;

public:

    LoggerManager(const LoggerManager&) = delete;
    LoggerManager& operator=(const LoggerManager&) = delete;

    void level(enum Logger::level level);

    enum Logger::level level() const;

    char const* name() const;

    std::vector<std::shared_ptr<LogAppender> > const& appenders() const;

    void add_appender(std::shared_ptr<LogAppender> const& appender);

    LoggerManager(char const* name, enum Logger::level level = Logger::L_WARNING);
    // LoggerManager(char const* name, enum Logger::level level = Logger::L_INFO);

    static void register_logger_manager(char const* logger_name_pattern,
                                        std::shared_ptr<LoggerManager> const& manager);

    static std::shared_ptr<LoggerManager> get_logger_manager(char const* logger_name_patern);

protected:
    void manage(Logger* logger);

protected:
    std::string const name_;
    enum Logger::level level_;
    std::set<Logger*> managed_loggers_;
    std::vector<std::shared_ptr<LogAppender> > appenders_;
};

class LogAppender
{
public:
    virtual ~LogAppender();

    virtual void flush() = 0;

    virtual void operator()(enum Logger::level lv,
                            char const* name, char const** chunks) = 0;
};

#define LOG_DEBUG(args) if (log_.level() == Logger::L_DEBUG) log_.debug args

#define LOG_INFO(args) if (enum Logger::level const level = log_.level()) if (level <= Logger::L_INFO) log_.info args

#define LOG_WARNING(args) if (enum Logger::level const level = log_.level()) if (level <= Logger::L_WARNING) log_.warn args

#define LOG_ERROR(args) if (enum Logger::level const level = log_.level()) if (level <= Logger::L_ERROR) log_.error args

} //egfrd
} //ecell4
#endif /* LOGGER_HPP */
