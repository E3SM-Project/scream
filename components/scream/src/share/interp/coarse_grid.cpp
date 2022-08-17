#include <share/interp/coarse_grid.hpp>
#include <share/io/scream_scorpio_interface.hpp>

namespace scream {
namespace interpolators {

CoarseGrid::CoarseGrid(const std::string& filename) {
}

namespace {

// This function computes the Jacobian of the coordinate transformation
// (a, b) -> (lon, lat, 1) at the reference point (a, b).
// s = (lon, lat, 1)
// s_ab <-> jacobian (derivative of s w.r.t. (a, b)
void ref2sphere_deriv(Real corners[4][3], Real a, Real b,
                      Real s_ab[3][2], Real s[3]) {

  // compute shape functions
  Real q[4] = {
    0.25 * (1.0-a)*(1.0-b),
    0.25 * (1.0+a)*(1.0-b),
    0.25 * (1.0+a)*(1.0+b),
    0.25 * (1.0-a)*(1.0+b)
  };

  // compute spherical coordinates
  s[0] = s[1] = s[2] = 0.0;
  for (int v = 0; v < 4; ++v) {
    for (int i = 0; i < 3; ++i) {
      s[i] += q[v] * corners[v][i];
    }
  }

  // compute jacobian
  Real r2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
  Real r_inv = 1.0/std::sqrt(r2);
  Real q_ab[2][4] = {
    {-0.25*(1.0-b),  0.25*(1.0-b), 0.25*(1.0+b), -0.25*(1.0+b)},
    {-0.25*(1.0-a), -0.25*(1.0+a), 0.25*(1.0+a),  0.25*(1.0-a)}
  };
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 3; ++i) {
      s_ab[i][j] = 0.0;
      for (int v = 0; v < 4; ++v) {
        s_ab[i][j] += corners[v][i] * q_ab[j][4];
      }
    }
    for (int i = 0; i < 3; ++i) {
      s_ab[i][j] = r_inv * s_ab[i][j] -
        r_inv*r_inv*r_inv *
        ((s[0]*s_ab[0][j]+s[1]*s_ab[1][j]+s[2]*s_ab[2][j]))*s[i];
    }
  }

  // normalize coordinates
  s[0] *= r_inv, s[1] *= r_inv, s[2] *= r_inv;
}

}

std::pair<Real, Real> CoarseGrid::
ll_to_ref(int e, Real lon, Real lat, Real tol, int max_iter) const {

  // extract element corners (on a unit sphere)
  const auto& elem = elements[e];
  Real c[4][3];
  for (int v = 0; v < 4; ++v) {
    c[v][0] = longitudes[elem.vertices[v]];
    c[v][1] = latitudes[elem.vertices[v]];
    c[v][2] = 1.0;
  }
  Real s0[3] = {lon, lat, 1.0}; // spherical coordinates

  Real tol2 = tol*tol;

  // Find (a, b) using Newton iteration.
  Real a = 0.0, b = 0.0;
  Real s[3], s_ab[3][2], r[3], fac[3], x[2];
  for (int iter = 0; iter < max_iter; ++iter) {
    // compute the jacobian s_ab
    ref2sphere_deriv(c, a, b, s_ab, s);

    // evaluate the residual r
    Real r[3] = {
      s[0] - s0[0],
      s[1] - s0[1],
      s[2] - s0[2]
    };
    if ((r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) <= tol2) break;

    // solve s_ab x = r using QR factorization

    // Q
    fac[0] = std::sqrt(s_ab[0][0]*s_ab[0][0] +
                       s_ab[1][0]*s_ab[1][0] +
                       s_ab[2][0]*s_ab[2][0]);
    for (int i = 0; i < 3; ++i)
      s_ab[i][0] /= fac[0];

    fac[1] = s_ab[0][1]*s_ab[0][0] +
             s_ab[1][1]*s_ab[1][0] +
             s_ab[2][1]*s_ab[2][0];
    for (int i = 0; i < 3; ++i)
      s_ab[i][1] -= fac[1]*s_ab[i][0];

    fac[2] = std::sqrt(s_ab[0][1]*s_ab[0][1] +
                       s_ab[1][1]*s_ab[1][1] +
                       s_ab[2][1]*s_ab[2][1]);
    for (int i = 0; i < 3; ++i)
      s_ab[i][1] /= fac[2];

    // x = Q'r
    x[0] = s_ab[0][0]*r[0] + s_ab[1][0]*r[1] + s_ab[2][0]*r[2];
    x[1] = s_ab[0][1]*r[0] + s_ab[1][1]*r[1] + s_ab[2][1]*r[2];

    // x = R \ x
    x[1] /= fac[2];
    x[0] = (x[0] - fac[1]*x[1]) / fac[0];

    // newton update
    a -= x[0];
    b -= x[1];
  }
  return std::make_pair(a, b);
}

} // namespace interpolators
} // namespace scream


