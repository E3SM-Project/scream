#include <catch2/catch.hpp>

#include "share/io/scream_scorpio_interface.hpp"
#include <ekat/util/ekat_string_utils.hpp>

namespace scream {

using namespace scorpio;

TEST_CASE ("io_subsystem") {
  ekat::Comm comm (MPI_COMM_WORLD);

  init_pio_subsystem (comm);
  REQUIRE_THROWS (init_pio_subsystem(comm));

  finalize_pio_subsystem ();
  REQUIRE_THROWS (finalize_pio_subsystem());
}

TEST_CASE ("write_and_read") {
  ekat::Comm comm (MPI_COMM_WORLD);

  EKAT_REQUIRE_MSG (comm.size()<=4,
      "Error! This test is tailored for an MPI_Comm of size 1, 2, 3, or 4.\n"
      " - MPI_Comm size: " + std::to_string(comm.size()) + "\n");

  init_pio_subsystem (comm);

  std::string filename = "scorpio_interface_write_test_np" + std::to_string(comm.size()) + ".nc";

  const int dim1 = 2;
  const int dim2 = 4;
  const int dim3 = 12;
  const int ldim3 = dim3 / comm.size();

  // Offsets for dim3 decomp owned by this rank
  std::vector<offset_t> my_offsets;
  for (int i=0; i<ldim3; ++i) {
    my_offsets.push_back(ldim3*comm.rank() + i);
  }

  // Write phase
  {
    register_file (filename,Write);
    REQUIRE_THROWS (register_file (filename,Read)); // ERROR: cannot open in both read and write modes
    REQUIRE (is_file_open(filename));
    REQUIRE (is_file_open(filename,Write));
    REQUIRE (not is_file_open(filename,Read));

    define_dim (filename,"dim1",dim1);
    define_dim (filename,"dim2",dim2);
    define_dim (filename,"dim2",dim2); // OK, same specs
    REQUIRE_THROWS (define_dim (filename,"dim2",dim2+3)); // ERROR: changing dim length
    define_dim (filename,"dim3",dim3);

    REQUIRE_THROWS (set_dim_decomp (filename,"dim4",my_offsets)); // ERROR: dimension not found

    set_dim_decomp (filename,"dim3",my_offsets);

    REQUIRE_THROWS (define_var (filename,"var1",{"dim1"},"double",true)); // ERROR: no time dimension (yet)
    REQUIRE_THROWS (define_var (filename,"var1",{"dim0"},"double",false)); // ERROR: dim0 not found
    REQUIRE_THROWS (define_var (filename,"var1",{"dim1"},"complex",false)); // ERROR: unsupported dtype
    define_var (filename,"var1",{"dim1"},"double",false);
    define_var (filename,"var1",{"dim1"},"double",false); // OK, same specs
    REQUIRE_THROWS (define_var (filename,"var1",{"dim2"},"double",false)); // ERROR: changing var dimensions
    REQUIRE_THROWS (define_var (filename,"var1",{"dim1"},"int",false)); // ERROR: changing dtype

    define_time(filename,"some_units","the_time");
    REQUIRE_THROWS (define_time(filename,"","another_name"));

    define_var (filename,"var2",{"dim1","dim2"},"float",true);
    define_var (filename,"var3",{},"int",true);
    REQUIRE_THROWS (define_var (filename,"var3",{},"int",false)); // ERROR: changing time_dep flag
    define_var (filename,"var4",{"dim3","dim1"},"double",false);
    define_var (filename,"var5",{"dim3","dim1"},"double",true);
    set_vars_decomp (filename,{"var4","var5"},"dim3");

    enddef (filename);

    std::vector<double> var1 (dim1);
    std::vector<float> var2 (dim1*dim2);
    std::vector<int> var3 (1);
    std::vector<double> var45 (ldim3*dim1);

    // Write first time slice
    update_time (filename,0.0);
    std::iota (var1.begin(),var1.end(),100);
    std::iota (var2.begin(),var2.end(),100);
    std::iota (var3.begin(),var3.end(),100);
    std::iota (var45.begin(),var45.end(),100+comm.rank()*ldim3*dim1);
    write_var (filename,"var1",var1.data());
    write_var (filename,"var2",var2.data());
    write_var (filename,"var3",var3.data());
    write_var (filename,"var4",var45.data());
    write_var (filename,"var5",var45.data());

    REQUIRE_THROWS (write_var (filename,"var3",static_cast<int*>(nullptr))); // ERROR: invalid pointer

    // Write second time slice
    update_time (filename,0.5);
    std::iota (var2.begin(),var2.end(),  200);
    std::iota (var3.begin(),var3.end(),  200);
    std::iota (var45.begin(),var45.end(),200+comm.rank()*ldim3*dim1);
    write_var (filename,"var2",var2.data());
    write_var (filename,"var3",var3.data());
    double minus_one = -1;
    write_var (filename,"var5",var45.data(),&minus_one);

    // Cleanup
    release_file (filename);

    REQUIRE_THROWS (release_file (filename)); // ERROR: file not open
    REQUIRE (not is_file_open(filename));
  }
  return;
  // Read phase
  {
    register_file (filename,Read);
    REQUIRE (is_file_open(filename));

    // Check dims
    REQUIRE (has_dimension(filename,"the_time",2));
    REQUIRE (has_dimension(filename,"dim1",dim1));
    REQUIRE (has_dimension(filename,"dim2",dim2));
    REQUIRE (has_dimension(filename,"dim3",dim3));

    // Check vars
    REQUIRE (has_variable(filename,"the_time",{},"","yes"));
    REQUIRE (has_variable(filename,"var1",{"dim1"},"","no"));
    REQUIRE (has_variable(filename,"var2",{"dim1","dim2"},"","yes"));
    REQUIRE (has_variable(filename,"var3",{},"","yes"));
    REQUIRE (has_variable(filename,"var4",{"dim3","dim1"},"","no"));
    REQUIRE (has_variable(filename,"var5",{"dim3","dim1"},"","yes"));

    REQUIRE (not has_variable(filename,"var0"));                  // Var not in file
    REQUIRE (not has_variable(filename,"var3",{"dim1"}));         // Var has different dims
    REQUIRE (not has_variable(filename,"var3",{},"some_units"));  // Var has different units

    std::vector<double> var1 (dim1);
    std::vector<float> var2 (dim1*dim2);
    std::vector<int> var3 (1);
    std::vector<double> var45 (ldim3*dim2);

    std::vector<double> tgt_var1 (dim1);
    std::vector<float> tgt_var2 (dim1*dim2);
    std::vector<int> tgt_var3 (1);
    std::vector<double> tgt_var45 (ldim3*dim1);
    std::iota (tgt_var1.begin(),tgt_var1.end(),  100);
    std::iota (tgt_var2.begin(),tgt_var2.end(),  100);
    std::iota (tgt_var3.begin(),tgt_var3.end(),  100);
    std::iota (tgt_var45.begin(),tgt_var45.end(),100+comm.rank()*ldim3*dim1);

    read_var (filename,"var1",var1.data());
    REQUIRE (tgt_var1==var1);

    // Read first time slice
    REQUIRE_THROWS (read_var (filename,"var3",var3.data(),3)); // ERROR: time_idx out of bounds
    REQUIRE_THROWS (read_var (filename,"var3",static_cast<int*>(nullptr))); // ERROR: invalid pointer

    read_var (filename,"var2",var2.data(),0);
    REQUIRE (tgt_var2==var2);

    read_var (filename,"var3",var3.data(),0);
    REQUIRE (tgt_var3==var3);

    read_var (filename,"var4",var45.data(),0);
    REQUIRE (tgt_var45==var45);

    read_var (filename,"var5",var45.data(),0);
    REQUIRE (tgt_var45==var45);

    // Read second time slice
    std::iota (tgt_var2.begin(),tgt_var2.end(),  200);
    std::iota (tgt_var3.begin(),tgt_var3.end(),  200);
    std::iota (tgt_var45.begin(),tgt_var45.end(),200+comm.rank()*ldim3*dim2);

    read_var (filename,"var2",var2.data(),1);
    REQUIRE (tgt_var2==var2);

    read_var (filename,"var3",var3.data(),1);
    REQUIRE (tgt_var3==var3);

    read_var (filename,"var5",var45.data(),1);
    REQUIRE (tgt_var45==var45);

    // Cleanup
    release_file (filename);
  }

  finalize_pio_subsystem ();
}

} // namespace scream
