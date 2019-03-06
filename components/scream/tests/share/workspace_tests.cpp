#include "catch2/catch.hpp"

#include "share/scream_workspace.hpp"
#include "share/util/scream_kokkos_utils.hpp"

namespace unit_test {

using namespace scream;

struct UnitWrap {

template <typename DeviceType>
struct UnitTest {

using Device     = DeviceType;
using MemberType = typename KokkosTypes<Device>::MemberType;
using TeamPolicy = typename KokkosTypes<Device>::TeamPolicy;
using ExeSpace   = typename KokkosTypes<Device>::ExeSpace;

template <typename S>
using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
template <typename S>
using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

static void unittest_workspace()
{
  using namespace scream;

  int nerr = 0;
  const int ints_per_ws = 37;
  static constexpr const int num_ws = 4;
  const int ni = 128;
  const int nk = 128;

  TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk));

  {
    WorkspaceManager<double, Device> wsmd(17, num_ws, policy);
    REQUIRE(wsmd.m_reserve == 1);
    REQUIRE(wsmd.m_size == 17);
  }
  {
    WorkspaceManager<char, Device> wsmc(16, num_ws, policy);
    REQUIRE(wsmc.m_reserve == 8);
    REQUIRE(wsmc.m_size == 16);
    Kokkos::parallel_for(
      "unittest_workspace char", policy,
      KOKKOS_LAMBDA(const MemberType& team) {
        auto ws = wsmc.get_workspace(team);
        const auto t1 = ws.take("t1");
        const auto t2 = ws.take("t1");
        ws.release(t1);
        ws.release(t2);
      });
  }
  {
    WorkspaceManager<short, Device> wsms(16, num_ws, policy);
    REQUIRE(wsms.m_reserve == 4);
    REQUIRE(wsms.m_size == 16);
  }

  // Test host-explicit WorkspaceMgr
  {
    using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
    typename KokkosTypes<HostDevice>::TeamPolicy policy_host(util::ExeSpaceUtils<typename KokkosTypes<HostDevice>::ExeSpace>::get_default_team_policy(ni, nk));
    WorkspaceManager<short, HostDevice> wsmh(16, num_ws, policy_host);
    wsmh.m_data(0, 0) = 0; // check on cuda machine
  }

  WorkspaceManager<int, Device> wsm(ints_per_ws, num_ws, policy);

  REQUIRE(wsm.m_reserve == 2);
  REQUIRE(wsm.m_size == ints_per_ws);

  Kokkos::parallel_reduce("unittest_workspace", policy, KOKKOS_LAMBDA(const MemberType& team, int& total_errs) {

    int nerrs_local = 0;
    auto ws = wsm.get_workspace(team);

    // Test getting workspaces of different type
    {
      const auto ws_int = ws.take("ints");
      // These nerrs_local increments are write race conditions among threads in
      // a team, but that's OK: nerrs_local doesn't have to be accurate. A 0
      // result will be a true 0 result.
      if (ws_int.extent(0) != ints_per_ws) ++nerrs_local;
      ws.release(ws_int);

      const auto ws_dlb = ws.template take<double>("doubles");
      if (ws_dlb.extent(0) != 18) ++nerrs_local;
      ws.release(ws_dlb);
    }
    team.team_barrier();

    Kokkos::Array<ko::Unmanaged<view_1d<int> >, num_ws> wssub;

    // Main test. Test different means of taking and release spaces.
    for (int r = 0; r < 8; ++r) {
      if (r % 4 == 0) {
        for (int w = 0; w < num_ws; ++w) {
          char buf[8] = "ws";
          buf[2] = 48 + w; // 48 is offset to integers in ascii
          wssub[w] = ws.take(buf);
        }
      }
      else {
        ko::Unmanaged<view_1d<int> > ws1, ws2, ws3, ws4;
        Kokkos::Array<ko::Unmanaged<view_1d<int> >*, num_ws> ptrs = { {&ws1, &ws2, &ws3, &ws4} };
        Kokkos::Array<const char*, num_ws> names = { {"ws0", "ws1", "ws2", "ws3"} };
        if (r % 4 == 1) {
          ws.take_many(names, ptrs);
        }
        else if (r % 4 == 2) {
          ws.take_many_contiguous_unsafe(names, ptrs);
        }
        else { // % 4 == 3
          ws.take_many_and_reset(names, ptrs);
        }

        for (int w = 0; w < num_ws; ++w) {
          wssub[w] = *ptrs[w];
        }
      }

      for (int w = 0; w < num_ws; ++w) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
          wssub[w](i) = i * w;
        });
      }

      team.team_barrier();

      for (int w = 0; w < num_ws; ++w) {
        // These spaces aren't free, but their metadata should be the same as it
        // was when they were initialized
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
          if (wsm.get_index(wssub[w]) != w) ++nerrs_local;
          if (wsm.get_next(wssub[w]) != w+1) ++nerrs_local;
          char buf[8] = "ws";
          buf[2] = 48 + w; // 48 is offset to integers in ascii
#ifndef NDEBUG
          if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
          if (ws.get_num_used() != 4) ++nerrs_local;
#endif
          for (int i = 0; i < ints_per_ws; ++i) {
            if (wssub[w](i) != i*w) ++nerrs_local;
          }
        });
      }

      team.team_barrier();

      if (r % 4 == 2) {
        // let take_and_reset do the reset
      }
      else if (r % 2 == 0) {
        ws.reset();
      }
      else {
        for (int w = num_ws - 1; w >= 0; --w) {
          ws.release(wssub[w]);
        }
      }

      team.team_barrier();
    }

#ifndef KOKKOS_ENABLE_CUDA
    // Test weird take/release permutations.
    for (int r = 0; r < 3; ++r) {
      int take_order[]    = {0, 1, 2, 3};
      int release_order[] = {-3, -2, -1, 0};

      do {
        for (int w = 0; w < num_ws; ++w) {
          char buf[8] = "ws";
          buf[2] = 48 + take_order[w]; // 48 is offset to integers in ascii
          wssub[take_order[w]] = ws.take(buf);
        }
        team.team_barrier();

        for (int w = 0; w < num_ws; ++w) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
            wssub[w](i) = i * w;
          });
        }

        team.team_barrier();

        // verify stuff
        for (int w = 0; w < num_ws; ++w) {
          Kokkos::single(Kokkos::PerTeam(team), [&] () {
            char buf[8] = "ws";
            buf[2] = 48 + w; // 48 is offset to integers in ascii
#ifndef NDEBUG
            if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
            if (ws.get_num_used() != 4) ++nerrs_local;
#endif
            for (int i = 0; i < ints_per_ws; ++i) {
              if (wssub[w](i) != i*w) ++nerrs_local;
            }
          });
        }

        team.team_barrier();

        for (int w = 0; w < num_ws; ++w) {
          ws.release(wssub[release_order[w] * -1]);
        }

        team.team_barrier();

        std::next_permutation(release_order, release_order+4);

      } while (std::next_permutation(take_order, take_order+4));
    }
    ws.reset();

    // Test weird take/release permutations.
    {
      int actions[] = {-3, -2, -1, 1, 2, 3};
      bool exp_active[] = {false, false, false, false};

      do {
        for (int a = 0; a < 6; ++a) {
          int action = actions[a];
          if (action < 0) {
            action *= -1;
            if (exp_active[action]) {
              ws.release(wssub[action]);
              exp_active[action] = false;
            }
          }
          else {
            if (!exp_active[action]) {
              char buf[8] = "ws";
              buf[2] = 48 + action; // 48 is offset to integers in ascii
              wssub[action] = ws.take(buf);
              exp_active[action] = true;
            }
          }
        }

        for (int w = 0; w < num_ws; ++w) {
          if (exp_active[w]) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ints_per_ws), [&] (Int i) {
              wssub[w](i) = i * w;
            });
          }
        }

        team.team_barrier();

        // verify stuff
        Kokkos::single(Kokkos::PerTeam(team), [&] () {
#ifndef NDEBUG
          int exp_num_active = 0;
#endif
          for (int w = 0; w < num_ws; ++w) {
            char buf[8] = "ws";
            buf[2] = 48 + w; // 48 is offset to integers in ascii
            if (exp_active[w]) {
#ifndef NDEBUG
              if (util::strcmp(ws.get_name(wssub[w]), buf) != 0) ++nerrs_local;
              ++exp_num_active;
              if (!ws.template is_active<int>(wssub[w])) ++nerrs_local;
#endif
              for (int i = 0; i < ints_per_ws; ++i) {
                if (wssub[w](i) != i*w) ++nerrs_local;
              }
            }
          }
#ifndef NDEBUG
          if (ws.get_num_used() != exp_num_active) ++nerrs_local;
#endif
        });

        team.team_barrier();

      } while (std::next_permutation(actions, actions + 6));
    }
    ws.reset();
#endif

    total_errs += nerrs_local;
    team.team_barrier();
  }, nerr);

#if 0
  wsm.report();
#endif

  REQUIRE(nerr == 0);
}

}; // struct UnitTest
}; // struct UnitWrap

} // namespace unit_test

namespace {

TEST_CASE("workspace_manager", "[utils]") {
  unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::unittest_workspace();
}

#ifdef KOKKOS_ENABLE_CUDA
// Force host testing on CUDA
TEST_CASE("workspace_manager_host", "[utils]") {
  unit_test::UnitWrap::UnitTest<scream::HostDevice>::unittest_workspace();
}
#endif

} // anonymous namespace
