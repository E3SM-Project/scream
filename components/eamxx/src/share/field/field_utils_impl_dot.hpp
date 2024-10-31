#ifndef SCREAM_FIELD_UTILS_IMPL_DOT_HPP
#define SCREAM_FIELD_UTILS_IMPL_DOT_HPP

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "share/field/field.hpp"

namespace scream {

// dot product of rank 1 field with a rank N field, return a rank N-1 field
template <typename ST>
Field do_dot_along_rank1_dim(const int &pd, const Field &f1, const Field &f2,
                             const ekat::Comm *co) {
  using KT          = ekat::KokkosTypes<DefaultDevice>;
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  using TeamPolicy  = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember  = typename TeamPolicy::member_type;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const auto &l1 = f1.get_header().get_identifier().get_layout();

  const auto &n2 = f2.get_header().get_identifier().name();
  const auto &l2 = f2.get_header().get_identifier().get_layout();
  const auto &u2 = f2.get_header().get_identifier().get_units();
  const auto &g2 = f2.get_header().get_identifier().get_grid_name();

  EKAT_REQUIRE_MSG(pd <= 5 && pd >= 0,
                   "Error! First argument pd must be between 0 and 5.\n"
                   "The input pd is "
                       << pd << ", which is not accepted.\n");
  EKAT_REQUIRE_MSG(l1.rank() == 1,
                   "Error! Second argument f1 must be rank-1.\n"
                   "The input f1 rank is "
                       << l1.rank() << ", which is not accepted.\n");
  EKAT_REQUIRE_MSG(l2.rank() <= 6,
                   "Error! Third argument f2 must be at most rank-6.\n"
                   "The input f2 rank is "
                       << l2.rank() << ", which is not accepted.\n");
  EKAT_REQUIRE_MSG(
      l1.dim(0) == l2.dim(pd),
      "Error! The two input fields must have the same dimension along "
      "which we are taking the dot product.\n"
      "The first field f1 has dimension "
          << l1.dim(0)
          << " while "
             "the second field f2 has dimension "
          << l2.dim(pd)
          << " \n"
             "along the provided dimension pd "
          << pd << " .\n");

  auto v1 = f1.get_view<const ST *>();

  FieldIdentifier fo_id(n2, l2.clone().strip_dim(pd), u2, g2);
  Field fo(fo_id);
  fo.allocate_view();
  fo.deep_copy(0);

  const int d0 = l2.dim(0);

  switch(l2.rank()) {
    case 1: {
      auto v2 = f2.get_view<ST *>();
      auto vo = fo.get_view<ST>();
      Kokkos::parallel_reduce(
          fo.name(), Kokkos::RangePolicy<>(0, d0),
          KOKKOS_LAMBDA(const int i, Real &ls) { ls += v1(i) * v2(i); }, vo);
    } break;
    case 2: {
      auto v2      = f2.get_view<const ST **>();
      auto vo      = fo.get_view<ST *>();
      const int d1 = l2.dim(1);
      if(pd == 0) {
        auto p = ESU::get_default_team_policy(d1, d0);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int j = tm.league_rank();
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d0),
                  [&](int i, ST &ac) { ac += v1(i) * v2(i, j); }, vo(j));
            });
      } else if(pd == 1) {
        auto p = ESU::get_default_team_policy(d0, d1);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int i = tm.league_rank();
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d1),
                  [&](int j, ST &ac) { ac += v1(j) * v2(i, j); }, vo(i));
            });
      }
    } break;
    case 3: {
      auto v2      = f2.get_view<const ST ***>();
      auto vo      = fo.get_view<ST **>();
      const int d1 = l2.dim(1);
      const int d2 = l2.dim(2);
      if(pd == 0) {
        auto p = ESU::get_default_team_policy(d1 * d2, d0);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int j   = idx / d2;
              const int k   = idx % d2;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d0),
                  [&](int i, ST &ac) { ac += v1(i) * v2(i, j, k); }, vo(j, k));
            });
      } else if(pd == 1) {
        auto p = ESU::get_default_team_policy(d0 * d2, d1);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = idx / d2;
              const int k   = idx % d2;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d1),
                  [&](int j, ST &ac) { ac += v1(j) * v2(i, j, k); }, vo(i, k));
            });
      } else if(pd == 2) {
        auto p = ESU::get_default_team_policy(d0 * d1, d2);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = idx / d1;
              const int j   = idx % d1;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d2),
                  [&](int k, ST &ac) { ac += v1(k) * v2(i, j, k); }, vo(i, j));
            });
      }
    } break;
    case 4: {
      auto v2      = f2.get_view<const ST ****>();
      auto vo      = fo.get_view<ST ***>();
      const int d1 = l2.dim(1);
      const int d2 = l2.dim(2);
      const int d3 = l2.dim(3);
      if(pd == 0) {
        auto p = ESU::get_default_team_policy(d1 * d2 * d3, d0);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int j   = (idx / d3) / d2;
              const int k   = (idx / d3) % d2;
              const int l   = idx % d3;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d0),
                  [&](int i, ST &ac) { ac += v1(i) * v2(i, j, k, l); },
                  vo(j, k, l));
            });
      } else if(pd == 1) {
        auto p = ESU::get_default_team_policy(d0 * d2 * d3, d1);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (idx / d3) / d2;
              const int k   = (idx / d3) % d2;
              const int l   = idx % d3;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d1),
                  [&](int j, ST &ac) { ac += v1(j) * v2(i, j, k, l); },
                  vo(i, k, l));
            });
      } else if(pd == 2) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d3, d2);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (idx / d3) / d1;
              const int j   = (idx / d3) % d1;
              const int l   = idx % d3;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d2),
                  [&](int k, ST &ac) { ac += v1(k) * v2(i, j, k, l); },
                  vo(i, j, l));
            });
      } else if(pd == 3) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2, d3);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (idx / d2) / d1;
              const int j   = (idx / d2) % d1;
              const int k   = idx % d2;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d3),
                  [&](int l, ST &ac) { ac += v1(l) * v2(i, j, k, l); },
                  vo(i, j, k));
            });
      }
    } break;
    case 5: {
      auto v2      = f2.get_view<const ST *****>();
      auto vo      = fo.get_view<ST ****>();
      const int d1 = l2.dim(1);
      const int d2 = l2.dim(2);
      const int d3 = l2.dim(3);
      const int d4 = l2.dim(4);
      if(pd == 0) {
        auto p = ESU::get_default_team_policy(d1 * d2 * d3 * d4, d0);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int j   = ((idx / d4) / d3) / d2;
              const int k   = ((idx / d4) / d3) % d2;
              const int l   = (idx / d4) % d3;
              const int n   = idx % d4;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d0),
                  [&](int i, ST &ac) { ac += v1(i) * v2(i, j, k, l, n); },
                  vo(j, k, l, n));
            });
      } else if(pd == 1) {
        auto p = ESU::get_default_team_policy(d0 * d2 * d3 * d4, d1);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = ((idx / d4) / d3) / d2;
              const int k   = ((idx / d4) / d3) % d2;
              const int l   = (idx / d4) % d3;
              const int n   = idx % d4;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d1),
                  [&](int j, ST &ac) { ac += v1(j) * v2(i, j, k, l, n); },
                  vo(i, k, l, n));
            });
      } else if(pd == 2) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d3 * d4, d2);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = ((idx / d4) / d3) / d1;
              const int j   = ((idx / d4) / d3) % d1;
              const int l   = (idx / d4) % d3;
              const int n   = idx % d4;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d2),
                  [&](int k, ST &ac) { ac += v1(k) * v2(i, j, k, l, n); },
                  vo(i, j, l, n));
            });
      } else if(pd == 3) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2 * d4, d3);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = ((idx / d4) / d2) / d1;
              const int j   = ((idx / d4) / d2) % d1;
              const int k   = (idx / d4) % d2;
              const int n   = idx % d4;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d3),
                  [&](int l, ST &ac) { ac += v1(l) * v2(i, j, k, l, n); },
                  vo(i, j, k, n));
            });
      } else if(pd == 4) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2 * d3, d4);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = ((idx / d3) / d2) / d1;
              const int j   = ((idx / d3) / d2) % d1;
              const int k   = (idx / d3) % d2;
              const int l   = idx % d3;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d4),
                  [&](int n, ST &ac) { ac += v1(n) * v2(i, j, k, l, n); },
                  vo(i, j, k, l));
            });
      }
    } break;
    case 6: {
      auto v2      = f2.get_view<const ST ******>();
      auto vo      = fo.get_view<ST *****>();
      const int d1 = l2.dim(1);
      const int d2 = l2.dim(2);
      const int d3 = l2.dim(3);
      const int d4 = l2.dim(4);
      const int d5 = l2.dim(5);
      if(pd == 0) {
        auto p = ESU::get_default_team_policy(d1 * d2 * d3 * d4 * d5, d0);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int j   = (((idx / d5) / d4) / d3) / d2;
              const int k   = (((idx / d5) / d4) / d3) % d2;
              const int l   = ((idx / d5) / d4) % d3;
              const int n   = (idx / d5) % d4;
              const int m   = idx % d5;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d0),
                  [&](int i, ST &ac) { ac += v1(i) * v2(i, j, k, l, n, m); },
                  vo(j, k, l, n, m));
            });
      } else if(pd == 1) {
        auto p = ESU::get_default_team_policy(d0 * d2 * d3 * d4 * d5, d1);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (((idx / d5) / d4) / d3) / d2;
              const int k   = (((idx / d5) / d4) / d3) % d2;
              const int l   = ((idx / d5) / d4) % d3;
              const int n   = (idx / d5) % d4;
              const int m   = idx % d5;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d1),
                  [&](int j, ST &ac) { ac += v1(j) * v2(i, j, k, l, n, m); },
                  vo(i, k, l, n, m));
            });
      } else if(pd == 2) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d3 * d4 * d5, d2);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (((idx / d5) / d4) / d3) / d1;
              const int j   = (((idx / d5) / d4) / d3) % d1;
              const int l   = ((idx / d5) / d4) % d3;
              const int n   = (idx / d5) % d4;
              const int m   = idx % d5;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d2),
                  [&](int k, ST &ac) { ac += v1(k) * v2(i, j, k, l, n, m); },
                  vo(i, j, l, n, m));
            });
      } else if(pd == 3) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2 * d4 * d5, d3);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (((idx / d5) / d4) / d2) / d1;
              const int j   = (((idx / d5) / d4) / d2) % d1;
              const int k   = ((idx / d5) / d4) % d2;
              const int n   = (idx / d5) % d4;
              const int m   = idx % d5;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d3),
                  [&](int l, ST &ac) { ac += v1(l) * v2(i, j, k, l, n, m); },
                  vo(i, j, k, n, m));
            });
      } else if(pd == 4) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2 * d3 * d5, d4);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (((idx / d5) / d3) / d2) / d1;
              const int j   = (((idx / d5) / d3) / d2) % d1;
              const int k   = ((idx / d5) / d3) % d2;
              const int l   = (idx / d5) % d3;
              const int m   = idx % d5;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d4),
                  [&](int n, ST &ac) { ac += v1(n) * v2(i, j, k, l, n, m); },
                  vo(i, j, k, l, m));
            });
      } else if(pd == 5) {
        auto p = ESU::get_default_team_policy(d0 * d1 * d2 * d3 * d4, d5);
        Kokkos::parallel_for(
            fo.name(), p, KOKKOS_LAMBDA(const TeamMember &tm) {
              const int idx = tm.league_rank();
              const int i   = (((idx / d4) / d3) / d2) / d1;
              const int j   = (((idx / d4) / d3) / d2) % d1;
              const int k   = ((idx / d4) / d3) % d2;
              const int l   = (idx / d4) % d3;
              const int n   = idx % d4;
              Kokkos::parallel_reduce(
                  Kokkos::TeamVectorRange(tm, d5),
                  [&](int m, ST &ac) { ac += v1(m) * v2(i, j, k, l, n, m); },
                  vo(i, j, k, l, n));
            });
      }
    } break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank.\n");
  }
  Kokkos::fence();
  if(co) {
    fo.sync_to_host();
    co->all_reduce(fo.template get_internal_view_data<ST, Host>(),
                   l2.size() / l2.dim(pd), MPI_SUM);
    fo.sync_to_dev();
  }
  return fo;
}

}  // namespace scream

#endif  // SCREAM_FIELD_UTILS_IMPL_HPP
