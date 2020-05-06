
#define BOOST_LOG_DYN_LINK 1
//#define BOOST_SYSTEM_NO_DEPRECATED

#include <MoFEM.hpp>
#include <LogManager.hpp>

#include <boost/log/utility/setup/console.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
// #include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;

namespace boost {
namespace log {
namespace keywords {

// // We define our own severity levels
// enum severity_level
// {
//     normal,
//     notification,
//     warning,
//     error,
//     critical
// };


BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
// BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", severity_level)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope", attrs::named_scope::value_type)
BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline", attrs::timer::value_type)

} 
} 
}

#if 0

//[ example_tutorial_file_simple
void init()
{
    logging::add_file_log("sample.log");

    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::info
    );
}
//]

// We need this due to this bug: https://svn.boost.org/trac/boost/ticket/4416
//[ example_tutorial_file_advanced_no_callouts
void init()
{
    logging::add_file_log
    (
        keywords::file_name = "sample_%N.log",
        keywords::rotation_size = 10 * 1024 * 1024,
        keywords::time_based_rotation = sinks::file::rotation_at_time_point(0, 0, 0),
        keywords::format = "[%TimeStamp%]: %Message%"
    );

    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::info
    );
}
//]

#else

//[ example_tutorial_file_advanced
void init()
{
  // logging::add_file_log(keywords::format =
  //                           "[%TimeStamp%]: %Message%" /*< log record format >*/
  // );



  logging::core::get()->set_filter(logging::trivial::severity >=
                                   logging::trivial::info);
}
//]

#endif



using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;    

    // LogManager::SynchronizedStream synchronized_stream(m_field.get_comm());

    // std::ostringstream ss;
    // ss << "Hello world from proc " << m_field.get_comm_rank() << endl;
    // synchronized_stream << ss;
    // synchronized_stream.flush();

    init();
    logging::add_common_attributes();

    using namespace logging::trivial;
    using namespace boost::log::keywords;

    src::severity_logger< severity_level > lg;

    // Construct the sink
    // auto sink = logging::add_console_log();

    // sink->set_formatter(
    //     expr::stream << std::hex << std::setw(8) << std::setfill('0')
    //                  << line_id
    //                  //  << std::dec << std::setfill(' ') << ": <" << severity
    //                  //  << ">\t"
    //                  << "(" << scope << ") "
    //                  << expr::if_(expr::has_attr(
    //                         tag_attr))[expr::stream << "[" << tag_attr << "] "]
    //                 //  << expr::if_(expr::has_attr(
    //                 //         timeline))[expr::stream << "[" << timeline << "] "]
    //                  << expr::smessage);

    // Add a stream to write log to

    LogManager::WorldStreamBuf sync_buf(m_field.get_comm());
    auto stream_ptr  = boost::make_shared<std::ostream>(&sync_buf);
    // (*stream_ptr) << "AAAAA" << std::endl;

    // sink->locked_backend()->add_stream(
        // boost::make_shared<std::ostream>(&sync_buf));

    // sink->locked_backend()->add_stream(
    //     boost::make_shared<std::ofstream>("sample.log"));

    // Register the sink in the logging core
    // logging::core::get()->add_sink(sink);

    // We have to provide an empty deleter to avoid destroying the global stream object
    // boost::shared_ptr<std::ostream> stream(&stream, boost::null_deleter());
    // sink->locked_backend()->add_stream(stream_ptr);


    auto core_log = logging::core::get();
    auto backend = boost::make_shared<sinks::text_ostream_backend>();
    // backend->add_stream(
    //     boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter()));
    backend->add_stream(stream_ptr);
    backend->auto_flush(true);
    
    typedef sinks::synchronous_sink< sinks::text_ostream_backend > sink_t;
    boost::shared_ptr< sink_t > sink(new sink_t(backend));
    core_log->add_sink(sink);

    BOOST_LOG_SEV(lg, trace) << "A trace severity message";
    BOOST_LOG_SEV(lg, debug) << "A debug severity message";
    BOOST_LOG_SEV(lg, info) << "An informational severity message";
    BOOST_LOG_SEV(lg, warning) << "A warning severity message";
    BOOST_LOG_SEV(lg, error) << "An error severity message";
    BOOST_LOG_SEV(lg, fatal) << "A fatal severity message";

  
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}