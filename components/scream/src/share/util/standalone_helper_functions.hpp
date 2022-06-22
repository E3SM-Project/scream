#ifndef SCREAM_STANDALONE_TEST_UTILS_HPP
#define SCREAM_STANDALONE_TEST_UTILS_HPP

#include "share/field/field_manager.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

namespace scream {

// Helper function to check the total water mass
Real calculate_water_mass(const std::shared_ptr<const GridsManager>& grids_mgr,const std::shared_ptr<FieldManager>& field_mgr, const bool use_precip) {
  
  using PC = scream::physics::Constants<Real>;
  constexpr Real gravit = PC::gravit;

  const auto& grid = grids_mgr->get_grid("Point Grid");

  int ncol  = grid->get_num_local_dofs();
  int nlay  = grid->get_num_vertical_levels();
  int dof   = ncol*nlay;
  REQUIRE(field_mgr->has_field("pseudo_density"));
  auto d_pdel       = field_mgr->get_field("pseudo_density").get_view<Real**>();

  // Start adding water species to mass calculation
  Real total_mass = 0.0;
  if (field_mgr->has_field("qv")) {
    auto d_tmp         = field_mgr->get_field("qv").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qr")) {
    auto d_tmp         = field_mgr->get_field("qr").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qc")) {
    auto d_tmp         = field_mgr->get_field("qc").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qi")) {
    auto d_tmp         = field_mgr->get_field("qi").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }
  if (field_mgr->has_field("qs")) {
    auto d_tmp         = field_mgr->get_field("qs").get_view<Real**>();
    Real result;
    Kokkos::parallel_reduce("",dof,KOKKOS_LAMBDA(const int& i,Real& lsum) {
      int icol = i / nlay;
      int klev = i % nlay;
      lsum += d_tmp(icol,klev) * d_pdel(icol,klev) / gravit;
    },result);
    total_mass += result;
  }

  if (use_precip) {
    auto d_tmp_liq = field_mgr->get_field("precip_liq_surf").get_view<Real*>();
    auto d_tmp_ice = field_mgr->get_field("precip_ice_surf").get_view<Real*>();
    Real result;
    Kokkos::parallel_reduce("",ncol, KOKKOS_LAMBDA(const int& icol,Real& lsum) {
      lsum += (d_tmp_liq(icol) + d_tmp_ice(icol) ) * 1000.0;
    },result);
    total_mass -= result;
  }

  return total_mass;
}

Real calculate_water_mass(const std::shared_ptr<const GridsManager>& grids_mgr,const std::shared_ptr<FieldManager>& field_mgr) {
  return calculate_water_mass(grids_mgr,field_mgr,false);
}

} //namespace scream
#endif //SCREAM_STANDALONE_TEST_UTILS_HPP
