
#ifndef SCREAM_LOGGING_HPP
#define SCREAM_LOGGING_HPP

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <stdexcept>
#include "scream_config.h"

#if LOG_STACK_TRACE
#include <execinfo.h>
#endif

namespace scream {
namespace util {
struct Error : public std::runtime_error {
  explicit Error(const std::string &s) : std::runtime_error(s) {}
};
}  // namespace util
}  // namespace scream

#include <assert.h>
#include <iostream>
#include <sstream>
#include <ctime>

namespace scream {
namespace util   {

inline void InitLogging(const char*) {
  // DO NOTHING
}

class LogCheckError {
 public:
  LogCheckError() : str(nullptr) {}
  explicit LogCheckError(const std::string& str_) : str(new std::string(str_)) {}
  ~LogCheckError() { if (str != nullptr) delete str; }
  operator bool() {return str != nullptr; }
  std::string* str;
};


#define DEFINE_CHECK_FUNC(name, op)                               \
  template <typename X, typename Y>                               \
  inline LogCheckError LogCheck##name(const X& x, const Y& y) {   \
    if (x op y) return LogCheckError();                           \
    std::ostringstream os;                                        \
    os << " (" << x << " vs. " << y << ") ";                      \
    return LogCheckError(os.str());                               \
  }                                                               \
  inline LogCheckError LogCheck##name(int x, int y) {             \
    return LogCheck##name<int, int>(x, y);                        \
  }

#define CHECK_BINARY_OP(name, op, x, y)                                              \
  if (scream::util::LogCheckError _check_err = scream::util::LogCheck##name(x, y))   \
      scream::util::LogMessageFatal(__FILE__, __LINE__).stream()                     \
      << "Check failed: " << #x " " #op " " #y << *(_check_err.str)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
DEFINE_CHECK_FUNC(_LT, <)
DEFINE_CHECK_FUNC(_GT, >)
DEFINE_CHECK_FUNC(_LE, <=)
DEFINE_CHECK_FUNC(_GE, >=)
DEFINE_CHECK_FUNC(_EQ, ==)
DEFINE_CHECK_FUNC(_NE, !=)
#pragma GCC diagnostic pop

// Always-on checking
#define SCREAM_CHECK(x)                                           \
  if (!(x))                                                       \
    stream::util::LogMessageFatal(__FILE__, __LINE__).stream()    \
      << "Check failed: " #x << ' '
#define SCREAM_CHECK_LT(x, y) CHECK_BINARY_OP(_LT, <, x, y)
#define SCREAM_CHECK_GT(x, y) CHECK_BINARY_OP(_GT, >, x, y)
#define SCREAM_CHECK_LE(x, y) CHECK_BINARY_OP(_LE, <=, x, y)
#define SCREAM_CHECK_GE(x, y) CHECK_BINARY_OP(_GE, >=, x, y)
#define SCREAM_CHECK_EQ(x, y) CHECK_BINARY_OP(_EQ, ==, x, y)
#define SCREAM_CHECK_NE(x, y) CHECK_BINARY_OP(_NE, !=, x, y)
#define SCREAM_CHECK_NOTNULL(x) \
  ((x) == NULL ? scream::util::LogMessageFatal(__FILE__, __LINE__).stream() << "Check  notnull: "  #x << ' ', (x) : (x)) // NOLINT(*)

// Debug-only checking.
#ifdef NDEBUG
#define SCREAM_DCHECK(x) \
  while (false) SCREAM_CHECK(x)
#define SCREAM_DCHECK_LT(x, y) \
  while (false) SCREAM_CHECK((x) < (y))
#define SCREAM_DCHECK_GT(x, y) \
  while (false) SCREAM_CHECK((x) > (y))
#define SCREAM_DCHECK_LE(x, y) \
  while (false) SCREAM_CHECK((x) <= (y))
#define SCREAM_DCHECK_GE(x, y) \
  while (false) SCREAM_CHECK((x) >= (y))
#define SCREAM_DCHECK_EQ(x, y) \
  while (false) SCREAM_CHECK((x) == (y))
#define SCREAM_DCHECK_NE(x, y) \
  while (false) SCREAM_CHECK((x) != (y))
#else
#define SCREAM_DCHECK(x) SCREAM_CHECK(x)
#define SCREAM_DCHECK_LT(x, y) SCREAM_CHECK((x) < (y))
#define SCREAM_DCHECK_GT(x, y) SCREAM_CHECK((x) > (y))
#define SCREAM_DCHECK_LE(x, y) SCREAM_CHECK((x) <= (y))
#define SCREAM_DCHECK_GE(x, y) SCREAM_CHECK((x) >= (y))
#define SCREAM_DCHECK_EQ(x, y) SCREAM_CHECK((x) == (y))
#define SCREAM_DCHECK_NE(x, y) SCREAM_CHECK((x) != (y))
#endif  // NDEBUG

#if LOG_CUSTOMIZE
#define LOG_INFO scream::util::CustomLogMessage(__FILE__, __LINE__)
#else
#define LOG_INFO scream::util::LogMessage(__FILE__, __LINE__)
#endif

#define LOG_ERROR LOG_INFO
#define LOG_WARNING LOG_INFO
#define LOG_FATAL scream::util::LogMessageFatal(__FILE__, __LINE__)
#define LOG_QFATAL LOG_FATAL

#define LOG(severity) LOG_##severity.stream()
#define LG LOG_INFO.stream()
#define LOG_IF(severity, condition) \
  !(condition) ? (void)0 : scream::util::LogMessageVoidify() & LOG(severity)

#ifdef NDEBUG
#define LOG_DFATAL LOG_ERROR
#define DFATAL ERROR
#define DLOG(severity) true ? (void)0 : scream::util::LogMessageVoidify() & LOG(severity)
#define DLOG_IF(severity, condition) \
  (true || !(condition)) ? (void)0 : scream::util::LogMessageVoidify() & LOG(severity)
#else
#define LOG_DFATAL LOG_FATAL
#define DFATAL FATAL
#define DLOG(severity) LOG(severity)
#define DLOG_IF(severity, condition) LOG_IF(severity, condition)
#endif

class DateLogger {
 public:
  DateLogger() { }
  const char* HumanDate() {
    time_t time_value = time(NULL);
    struct tm *pnow;
    pnow = localtime(&time_value);  // NOLINT(*)
    snprintf(buffer_, sizeof(buffer_), "%02d:%02d:%02d",
             pnow->tm_hour, pnow->tm_min, pnow->tm_sec);
    return buffer_;
  }

 private:
  char buffer_[9];
};


class LogMessage {
 public:
  LogMessage(const char* file, int line)
      : log_stream_(std::cerr)
  {
    log_stream_ << "[" << pretty_date_.HumanDate() << "] " << file << ":"
                << line << ": ";
  }
  ~LogMessage() { log_stream_ << '\n'; }
  std::ostream& stream() { return log_stream_; }

 protected:
  std::ostream& log_stream_;

 private:
  DateLogger pretty_date_;
  LogMessage(const LogMessage&);
  void operator=(const LogMessage&);
};


class CustomLogMessage {
 public:
  CustomLogMessage(const char* file, int line) {
    log_stream_ << "[" << DateLogger().HumanDate() << "] " << file << ":"
                << line << ": ";
  }
  ~CustomLogMessage() {
    Log(log_stream_.str());
  }
  std::ostream& stream() { return log_stream_; }

  static void Log(const std::string& msg);

 private:
  std::ostringstream log_stream_;
};

#if LOG_FATAL_THROW == 0
class LogMessageFatal : public LogMessage {
 public:
  LogMessageFatal(const char* file, int line) : LogMessage(file, line) {}
  ~LogMessageFatal() {
#if LOG_STACK_TRACE
    const int MAX_STACK_SIZE = 10;
    void *stack[MAX_STACK_SIZE];

    int nframes = backtrace(stack, MAX_STACK_SIZE);
    log_stream_ << "\n\n" << "Stack trace returned " << nframes << " entries:\n";
    char **msgs = backtrace_symbols(stack, nframes);
    if (msgs != nullptr) {
      for (int i = 0; i < nframes; ++i) {
        log_stream_ << "[bt] (" << i << ") " << msgs[i] << "\n";
      }
    }
#endif

    log_stream_ << "\n";
    abort();
  }

 private:
  LogMessageFatal(const LogMessageFatal&);
  void operator=(const LogMessageFatal&);
};
#else
class LogMessageFatal {
 public:
  LogMessageFatal(const char* file, int line) {
    log_stream_ << "[" << pretty_date_.HumanDate() << "] " << file << ":"
                << line << ": ";
  }
  std::ostringstream &stream() { return log_stream_; }
  ~LogMessageFatal() THROW_EXCEPTION {
#if LOG_STACK_TRACE
    const int MAX_STACK_SIZE = 10;
    void *stack[MAX_STACK_SIZE];

    int nframes = backtrace(stack, MAX_STACK_SIZE);
    log_stream_ << "\n\n" << "Stack trace returned " << nframes << " entries:\n";
    char **msgs = backtrace_symbols(stack, nframes);
    if (msgs != nullptr) {
      for (int i = 0; i < nframes; ++i) {
        log_stream_ << "[bt] (" << i << ") " << msgs[i] << "\n";
      }
    }
#endif

#if LOG_BEFORE_THROW
    LOG(ERROR) << log_stream_.str();
#endif
    throw Error(log_stream_.str());
  }

 private:
  std::ostringstream log_stream_;
  DateLogger pretty_date_;
  LogMessageFatal(const LogMessageFatal&);
  void operator=(const LogMessageFatal&);
};
#endif

class LogMessageVoidify {
 public:
  LogMessageVoidify() {}

  void operator&(std::ostream&) {}
};

}  // namespace util
}  // namespace scream

#endif  // SCREAM_LOGGING_HPP
