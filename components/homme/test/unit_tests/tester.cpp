#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include <Kokkos_Core.hpp>

#include <Hommexx_Session.hpp>
#include <ErrorDefs.hpp>

#include <Context.hpp>
#include <mpi/Comm.hpp>

#include "test_session.hpp"

#include <mpi.h>

// Make command-line arguments available to tests. Anything following "hommexx"
// in the command-line argument list is passed to the tests.
int hommexx_catch2_argc;
char** hommexx_catch2_argv;

int main(int argc, char **argv) {

  // Initialize mpi
  MPI_Init(&argc,&argv);

  Homme::Context::singleton().create<Homme::Comm>().reset_mpi_comm(MPI_COMM_WORLD);
  Homme::initialize_hommexx_session();

  // Filter arguments so catch2 doesn't try to interpret hommexx-specific ones.
  hommexx_catch2_argc = argc;
  hommexx_catch2_argv = argv;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "hommexx") {
      argc = i;
      hommexx_catch2_argc -= i + 1;
      hommexx_catch2_argv = argv + i + 1;
      break;
    }
  }

  // Read possible ekat-specific arguments
  auto const readCommaSeparatedParams = [] (const std::string& cmd_line_arg) {
    if (cmd_line_arg=="") {
      return;
    }
    auto& ts = Homme::TestSession::instance();

    std::stringstream input(cmd_line_arg);
    std::string option;
    while (getline(input,option,',')) {
      // Split option according to key=val
      auto pos = option.find('=');
      Homme::Errors::runtime_check(pos!=std::string::npos, "Error! Badly formatted command line options.\n");
      std::string key = option.substr(0,pos);
      std::string val = option.substr(pos+1);
      Homme::Errors::runtime_check(key!="", "Error! Empty key found in the list of command line options.\n");
      Homme::Errors::runtime_check(val!="", "Error! Empty value for key '" + key + "'.\n");
      ts.params[key] = val;
    }
  };
  auto const readCommaSeparatedOptions = [] (const std::string& cmd_line_arg) {
    if (cmd_line_arg=="") {
      return;
    }
    auto& ts = Homme::TestSession::instance();

    std::stringstream input(cmd_line_arg);
    std::string option;
    while (getline(input,option,',')) {
      // Split option according to key=val
      ts.flags[option] = true;
    }
  };

  Catch::Session catch_session;
  auto cli = catch_session.cli();
  cli |= Catch::clara::Opt(readCommaSeparatedParams, "key1=val1[,key2=val2[,...]]")
             ["--params"]
             ("list of parameters to forward to the test");
  cli |= Catch::clara::Opt(readCommaSeparatedOptions, "option1[,option2[,...]]")
             ["--flags"]
             ("list of flags to forward to the test");

  catch_session.cli(cli);

  Homme::Errors::runtime_check(catch_session.applyCommandLine(argc,argv)==0,
                     "Error! Something went wrong while parsing command line.\n");

  int result = catch_session.run();

  Homme::finalize_hommexx_session();

  // Finalize mpi
  MPI_Finalize();

  return result;
}
