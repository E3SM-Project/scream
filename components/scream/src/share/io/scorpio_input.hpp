#ifndef SCREAM_SCORPIO_INPUT_HPP
#define SCREAM_SCORPIO_INPUT_HPP

#include "scream_config.h"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "ekat/ekat_pack_utils.hpp"

#include "share/io/scream_scorpio_interface.hpp"

#include "share/field/field_manager.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_identifier.hpp"

#include "share/grid/grids_manager.hpp"

/*  The AtmosphereInput class handles all input streams to SCREAM.
 *  It is important to note that there does not exist an InputManager,
 *  like in the case of output.  So all input streams have to be managed
 *  directly by the process that requires it.
 *
 *  The typical input call will be to the outward facing routine 'pull_input'.
 *
 *  Currently, input can only handle single timesnap input files.  In other words
 *  files that will be opened, read and closed within the same timestep and only
 *  store one timesnap of data.
 *
 *  Note: the init, and finalize are separate routines that are outward facing in
 *  this class to facilitate cases where reading input over some number of simulation
 *  timesteps is possible.
 *
 *  At construction time All input instances require
 *  1. an EKAT comm group and
 *  2. a EKAT parameter list
 *  3. a shared pointer to the field manager
 *  4. a shared pointer to the grids manager
 *  The parameter list contains all of the essential data regarding
 *  the input file name and the set of variables to be read.  The parameter list can be
 *  created localling in the process needing input, see src/share/io/tests/ for examples of
 *  setting up the input parameter list.
 *
 *  A typical input parameter list looks like:
 *  -----
 *  Input Parameters
 *    FILENAME: STRING
 *    GRID: STRING
 *    FIELDS
 *      Number of Fields: INT
 *      field_1: STRING
 *      ...
 *      field_N: STRING
 *  -----
 *  where,
 *  FILENAME: is a string value of the name of the input file to be read.
 *  GRID: is a string of the grid the input file is written on, currently only "Physics" is supported.
 *  FIELDS: designation of a sublist, so empty here
 *    Number of Fields: is an integer value>0 telling the number of fields
 *    field_x: is the xth field variable name.  Should match the name in the file and the name in the field manager.  TODO: add a rename option if variable names differ in file and field manager.
 *
 * Usage:
 * 1. Construct an instance of the AtmosphereInput class:
 *    AtmosphereInput in_object(comm,params);
 * 2. Use pull input to gather the desired data:
 *    in_object.pull_input(*grid_manager)
 *
 * The AtmosphereInput class will replace all fields in the field_manager that are part of the input with
 * data read from the input file.
 *    
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */

namespace scream
{

class AtmosphereInput 
{
public:
  using dofs_list_type = AbstractGrid::dofs_list_type;
  using view_type_host = typename KokkosTypes<DefaultDevice>::view_1d<Real>::HostMirror;

  // --- Constructor(s) & Destructor --- //
  AtmosphereInput (const ekat::Comm& comm, const ekat::ParameterList& params,
                   const std::shared_ptr<const FieldManager<Real>>& field_mgr,
                   const std::shared_ptr<const GridsManager>& grid_mgr)
    : m_params    (params)
    , m_comm      (comm)
    , m_field_mgr (field_mgr)
    , m_grid_mgr  (grid_mgr)
  {
    // Nothing to do here
  }

  AtmosphereInput (const ekat::Comm& comm, const std:: string grid_name,
                   const std::shared_ptr<const GridsManager>& grid_mgr)
    : m_comm      (comm)
    , m_grid_mgr  (grid_mgr)
    , m_grid_name (grid_name)
  {
    // Nothing to do here
  }

  virtual ~AtmosphereInput () = default;

  // --- Methods --- //

  void pull_input ();

  // Used by scorpio_output when handling restart history files.
  view_type_host pull_input (const std::string& name);

  // var_dims is a list of tags for each of the physical dimensions of this variable.
  // dim_lens is a vector of the physical dimension lengths without any padding, in other words
  // dim_lens is a vector of integer for the physical length of each of the tags in var_dims.
  // one last way to think of dim_lens would be the array dimensions lengths if ValueType=Real (no padding)
  void pull_input (const std::string& filename, const std::string& var_name,
                   const std::vector<std::string>& var_dims,
                   const bool has_columns, const std::vector<int>& dim_lens,
                   const int padding, Real* data);

  // Determine padding from the type of the variable.
  template<typename ValueType>
  void pull_input (const std::string& filename, const std::string& var_name,
                   const std::vector<std::string>& var_dims,
                   const bool has_columns, const std::vector<int>& dim_lens,
                   ValueType* data)
  {
    // Determine the padding for this data
    constexpr int pack_size = sizeof(ValueType) / sizeof(Real);
    const int padding = ekat::PackInfo<pack_size>::padding(dim_lens.back());
    // Make sure to pass the data as a Real pointer.
    auto data_real = reinterpret_cast<Real*>(data);
    pull_input(filename, var_name, var_dims, has_columns, dim_lens, padding, data_real);
  }

  void init();
  void finalize();

protected:
  // Internal functions
  void register_variables();
  void set_degrees_of_freedom();

  std::vector<std::string> get_vec_of_dims (const FieldLayout& layout);
  std::string get_io_decomp (const std::vector<std::string>& vec_of_dims);
  std::vector<Int> get_var_dof_offsets (const int dof_len, const bool has_cols);

  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;

  std::shared_ptr<const FieldManager<Real>>   m_field_mgr;
  std::shared_ptr<const GridsManager>         m_grid_mgr;
  
  std::string m_filename;
  std::string m_avg_type;
  std::string m_grid_name;

  std::vector<std::string>               m_fields_names;
  std::map<std::string,Int>              m_dofs_sizes;
  typename dofs_list_type::HostMirror    m_gids_host;

  bool m_is_rhist = false;

}; // Class AtmosphereInput

} //namespace scream

#endif // SCREAM_SCORPIO_INPUT_HPP
