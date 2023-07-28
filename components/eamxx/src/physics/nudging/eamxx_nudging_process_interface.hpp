#ifndef SCREAM_NUDGING_HPP
#define SCREAM_NUDGING_HPP

#include "share/util/eamxx_time_interpolation.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_lin_interp.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/util/scream_vertical_interpolation.hpp"
#include "share/util/scream_time_stamp.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the nudging of variables
*/

// enum to track how the source pressure levels are defined
enum SourcePresType {
  DYNAMIC = 0,  // DEFAULT - source data should include time/spatially varying p_mid
  STATIC  = 1,  // source data includes p_lev which is a static set of levels in both space and time.
};

class Nudging : public AtmosphereProcess
{
public:
  using mPack = ekat::Pack<Real,1>;
  using mMask = ekat::Mask<1>;
  using KT = KokkosTypes<DefaultDevice>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;

  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_Nd_host = typename KT::template view_ND<S,N>::HostMirror;

  template <typename S>
  using view_1d_host = view_Nd_host<S,1>;

  template <typename S>
  using view_2d_host = view_Nd_host<S,2>;

  // Constructors
  Nudging (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Nudging"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void run_impl        (const double dt);

protected:

  // The two other main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  // Internal function to apply nudging at specific timescale
  void apply_tendency(Field& base, const Field& next, const int dt);

  std::shared_ptr<const AbstractGrid>   m_grid;
  // Keep track of field dimensions and the iteration count
  int m_num_cols;
  int m_num_levs;
  int m_num_src_levs;
  int m_timescale;
  std::vector<std::string> m_datafiles;
  SourcePresType m_src_pres_type;
  

  std::vector<std::string> m_fields_nudge;

  util::TimeInterpolation m_time_interp;
}; // class Nudging

} // namespace scream

#endif // SCREAM_NUDGING_HPP
