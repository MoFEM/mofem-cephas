
#define BOOST_LOG_DYN_LINK 1
//#define BOOST_SYSTEM_NO_DEPRECATED

#include <MoFEM.hpp>
#include <LogManager.hpp>

#undef likely

#include <cstddef>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/core/null_deleter.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;

using namespace MoFEM;

static char help[] = "...\n\n";

// enum severity_level { normal, notification, warning, error, critical };

// namespace boost {
// namespace log {
// namespace expressions {
// BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
// BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
// BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
// BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope", attrs::named_scope::value_type)
// BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline", attrs::timer::value_type)
// } // namespace expressions
// } // namespace log
// } // namespace boost

// void logging_function() {
//   BOOST_LOG_NAMED_SCOPE("named_scope_logging");
//   src::severity_logger<severity_level> slg;

//   BOOST_LOG_SEV(slg, normal) << "A regular message";
//   BOOST_LOG_SEV(slg, warning)
//       << "Something bad is going on but I can handle it";
//   BOOST_LOG_SEV(slg, critical) << "Everything crumbles, shoot me now!";
// }

void named_scope_logging() {
  BOOST_LOG_NAMED_SCOPE("named_scope_logging");

  src::severity_logger<LogManager::SeverityLevel> slg;
  slg.add_attribute("Tag", attrs::constant<std::string>("My tag value"));

  BOOST_LOG_SEV(slg, LogManager::SeverityLevel::normal)
      << "Hello from the function named_scope_logging!";
}

// void tagged_logging() {
//   src::severity_logger<severity_level> slg;
//   slg.add_attribute("Tag", attrs::constant<std::string>("My tag value"));

//   // // BOOST_LOG_FUNC();
//   // BOOST_LOG_FUNCTION();
//   // clang-format off
//   BOOST_LOG_NAMED_SCOPE("aaa"); BOOST_LOG_SEV(slg, normal) << "Here goes the tagged record";
//   // clang-format on
// }

// void timed_logging() {
//   BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());

//   src::severity_logger<severity_level> slg;
//   BOOST_LOG_SEV(slg, normal) << "Starting to time nested functions";

//   logging_function();

//   BOOST_LOG_SEV(slg, normal) << "Stopping to time nested functions";
// }

// // The operator puts a human-friendly representation of the severity level to
// // the stream
// std::ostream &operator<<(std::ostream &strm, severity_level level) {
//   static const char *strings[] = {"normal", "notification", "warning", "error",
//                                   "critical"};

//   if (static_cast<std::size_t>(level) < sizeof(strings) / sizeof(*strings))
//     strm << strings[level];
//   else
//     strm << static_cast<int>(level);

//   return strm;
// }

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    // // LogManager::WorldStreamBuf sync_buf(m_field.get_comm());
    // // auto stream_ptr = boost::make_shared<std::ostream>(&sync_buf);
    // auto stream_ptr = m_field.getInterface<LogManager>()->getStrmWorld();

    // auto core_log = logging::core::get();
    // auto backend = boost::make_shared<sinks::text_ostream_backend>();
    // backend->add_stream(stream_ptr);
    // // backend->add_stream(
    //     // boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter()));
    // backend->auto_flush(true);

    // typedef sinks::synchronous_sink<sinks::text_ostream_backend> sink_t;
    // boost::shared_ptr<sink_t> sink(new sink_t(backend));

    // sink->set_formatter(

    //     expr::stream
    //     << std::hex << std::setw(8) << std::setfill('0')
    //     << boost::log::expressions::line_id << std::dec << std::setfill(' ')
    //     << ": <" << boost::log::expressions::severity << ">\t"
    //     << boost::log::expressions::format_named_scope(
    //            "Scope", keywords::format = "[%f:%l]")
    //     << "(" << boost::log::expressions::scope << ") "
    //     << expr::if_(expr::has_attr(boost::log::expressions::tag_attr))
    //            [expr::stream << "[" << boost::log::expressions::tag_attr
    //                          << "] "]
    //     << expr::if_(expr::has_attr(boost::log::expressions::timeline))
    //            [expr::stream << "[" << boost::log::expressions::timeline
    //                          << "] "]
    //     << expr::smessage

    // );

    // core_log->add_sink(sink);

    // logging::add_common_attributes();
    // core_log->add_global_attribute("LineID", attrs::counter<unsigned int>(1));
    // core_log->add_global_attribute("TimeStamp", attrs::local_clock());
    // core_log->add_global_attribute("Scope", attrs::named_scope());

    BOOST_LOG_FUNCTION();
    named_scope_logging();
    // logging_function();
    // tagged_logging();
    // timed_logging();

    // backend->flush();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}