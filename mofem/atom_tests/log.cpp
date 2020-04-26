


#include <MoFEM.hpp>
#include <LogManager.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;    

    LogManager::SynchronizedStream synchronized_stream(m_field.get_comm());

    std::ostringstream ss;
    ss << "Hello world from proc " << m_field.get_comm_rank() << endl;
    synchronized_stream << ss;
    synchronized_stream.flush();

  
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}