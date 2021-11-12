#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include <string>
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "share/field/field.hpp"
#include "share/grid//abstract_grid.hpp"

namespace scream {

template<typename RT>
RT& access_2d_dof (const Field<RT>& f, const std::shared_ptr<const AbstractGrid>& grid, const int idof, const int icomp = -1) {
  using namespace ShortFieldTagsNames;

  auto i2s = [](const int i) -> std::string { return std::to_string(i); };

  const auto& fh  = f.get_header();
  const auto& fid = fh.get_identifier();
  EKAT_REQUIRE_MSG (grid->name()==fid.get_grid_name(),
      "Error! Input grid is not the same as the one associated with the input field.\n"
      "   - input grid name: " + grid->name() + "\n"
      "   - field grid name: " + fid.get_grid_name() + "\n");
  const auto& lt = fid.get_layout();
  const auto lt_type = get_layout_type(lt.tags());
  EKAT_REQUIRE_MSG (lt_type==LayoutType::Scalar2D ||
                    lt_type==LayoutType::Vector2D,
      "Error! Function 'access_2d_dof' called on a non-2d field.\n");
  const bool vector = lt_type==LayoutType::Vector2D;
  EKAT_REQUIRE_MSG ( (vector && icomp>=0) || (not vector && icomp==-1),
      "Error! Vector component must be specified if and only if the field is a vector field.\n");
  EKAT_REQUIRE_MSG ( not vector || icomp<lt.dim(CMP),
      "Error! Vector component (" + i2s(icomp) + ") out of bounds [0," + i2s(lt.dim(CMP)) + ")\n");

  const auto& lid2idx = grid->get_lid_to_idx_map();
  EKAT_REQUIRE_MSG (idof>=0 && idof<=lid2idx.extent_int(0),
      "Error! Dof index (" + i2s(idof) + ") out of bounds [0," + i2s(lid2idx.extent(0)) + ").\n");

  const auto& idx = ekat::subview(lid2idx,idof);
  switch (grid->type()) {
    case GridType::Point: 
      {
        if (vector) {
          const auto v = f.template get_view<RT**,Host>();
          return v(idx(0),icomp);
        } else {
          const auto v = f.template get_view<RT*,Host>();
          return v(idx(0));
        }
      }
    case GridType::SE:
      {
        if (vector) {
          const auto v = f.template get_view<RT****,Host>();
          return v(idx(0),icomp,idx(1),idx(2));
        } else {
          const auto v = f.template get_view<RT***,Host>();
          return v(idx(0),idx(1),idx(2));
        }
      }
    default:
      EKAT_ERROR_MSG("Error! Unexpected grid type.\n");
  }
}

template<typename RT>
auto access_3d_column (const Field<RT>& f, const std::shared_ptr<const AbstractGrid>& grid, const int idof, const int icomp = -1) 
-> decltype(f.template get_view<RT*,Host>())
{
  using namespace ShortFieldTagsNames;

  auto i2s = [](const int i) -> std::string { return std::to_string(i); };

  const auto& fh  = f.get_header();
  const auto& fid = fh.get_identifier();
  EKAT_REQUIRE_MSG (grid->name()==fid.get_grid_name(),
      "Error! Input grid is not the same as the one associated with the input field.\n"
      "   - input grid name: " + grid->name() + "\n"
      "   - field grid name: " + fid.get_grid_name() + "\n");
  const auto& lt = fid.get_layout();
  const auto lt_type = get_layout_type(lt.tags());
  EKAT_REQUIRE_MSG (lt_type==LayoutType::Scalar3D ||
                    lt_type==LayoutType::Vector3D,
      "Error! Function 'access_2d_dof' called on a non-2d field.\n");
  const bool vector = lt_type==LayoutType::Vector3D;
  EKAT_REQUIRE_MSG ( (vector && icomp>=0) || (not vector && icomp==-1),
      "Error! Vector component must be specified if and only if the field is a vector field.\n");
  EKAT_REQUIRE_MSG ( not vector || icomp<lt.dim(CMP),
      "Error! Vector component (" + i2s(icomp) + ") out of bounds [0," + i2s(lt.dim(CMP)) + ")\n");

  const auto& lid2idx = grid->get_lid_to_idx_map();
  EKAT_REQUIRE_MSG (idof>=0 && idof<=lid2idx.extent_int(0),
      "Error! Dof index (" + i2s(idof) + ") out of bounds [0," + i2s(lid2idx.extent(0)) + ").\n");

  const auto& idx = ekat::subview(lid2idx,idof);
  switch (grid->type()) {
    case GridType::Point: 
      {
        if (vector) {
          const auto v = f.template get_view<RT**,Host>();
          return ekat::subview(v,idx(0),icomp);
        } else {
          const auto v = f.template get_view<RT*,Host>();
          return ekat::subview(v,idx(0));
        }
      }
    case GridType::SE:
      {
        if (vector) {
          const auto v = f.template get_view<RT****,Host>();
          return ekat::subview(v,idx(0),icomp,idx(1),idx(2));
        } else {
          const auto v = f.template get_view<RT***,Host>();
          return ekat::subview(v,idx(0),idx(1),idx(2));
        }
      }
    default:
      EKAT_ERROR_MSG("Error! Unexpected grid type.\n");
  }
}

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
template<typename RT1, typename RT2>
bool views_are_equal(const Field<RT1>& f1, const Field<RT2>& f2) {
  static_assert(
      std::is_same<typename std::remove_cv<RT1>::type,
                   typename std::remove_cv<RT2>::type>::value,
      "Error! Real types must be the same (except possibly for cv qualifiers).\n");

  // Get physical layout (shoudl be the same for both fields)
  const auto& l1 = f1.get_header().get_identifier().get_layout();
  const auto& l2 = f2.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (l1==l2,
      "Error! Input fields have different layouts.\n");

  // For simplicity, we perform the check on Host only. This is not a big
  // limitation, since this code is likely used only in testing.
  f1.sync_to_host();
  f2.sync_to_host();

  // Reshape based on the rank, then loop over all entries.
  const auto& dims = l1.dims();
  switch (l1.rank()) {
    case 1:
      {
        auto v1 = f1.template get_view<RT1*,Host>();
        auto v2 = f2.template get_view<RT2*,Host>();
        for (int i=0; i<dims[0]; ++i) {
          if (v1(i) != v2(i)) {
            return false;
          }
        }
      }
      break;
    case 2:
      {
        auto v1 = f1.template get_view<RT1**,Host>();
        auto v2 = f2.template get_view<RT2**,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            if (v1(i,j) != v2(i,j)) {
              return false;
            }
        }}
      }
      break;
    case 3:
      {
        auto v1 = f1.template get_view<RT1***,Host>();
        auto v2 = f2.template get_view<RT2***,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              if (v1(i,j,k) != v2(i,j,k)) {
                return false;
              }
        }}}
      }
      break;
    case 4:
      {
        auto v1 = f1.template get_view<RT1****,Host>();
        auto v2 = f2.template get_view<RT2****,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                if (v1(i,j,k,l) != v2(i,j,k,l)) {
                  return false;
                }
        }}}}
      }
      break;
    case 5:
      {
        auto v1 = f1.template get_view<RT1*****,Host>();
        auto v2 = f2.template get_view<RT2*****,Host>();
        for (int i=0; i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                for (int m=0; m<dims[4]; ++m) {
                  if (v1(i,j,k,l,m) != v2(i,j,k,l,m)) {
                    return false;
                  }
        }}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  // If we get here, then all entries matched.
  return true;
}

template<typename RT, typename Engine, typename PDF>
void randomize (const Field<RT>& f, Engine& engine, PDF&& pdf)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  const auto& fl = f.get_header().get_identifier().get_layout();
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_view<RT*,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          v(i) = pdf(engine);
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<RT**,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            v(i,j) = pdf(engine);
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<RT***,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              v(i,j,k) = pdf(engine);
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<RT****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                v(i,j,k,l) = pdf(engine);
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_view<RT*****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  v(i,j,k,l,m) = pdf(engine);
        }}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  // Sync the dev view with the host view.
  f.sync_to_dev();
}

inline Real frobenius_norm(const Field<const Real>& f)
{
  const auto& fl = f.get_header().get_identifier().get_layout();
  f.sync_to_host();

  // Note: use Kahan algorithm to increase accuracy
  Real norm = 0;
  Real c = 0;
  Real temp,y;
  switch (fl.rank()) {
    case 1:
      {
        auto v = f.template get_view<Real*,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          y = std::pow(v(i),2) - c;
          temp = norm + y;
          c = (temp - norm) - y;
          norm = temp;
        }
      }
      break;
    case 2:
      {
        auto v = f.template get_view<Real**,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            y = std::pow(v(i,j),2) - c;
            temp = norm + y;
            c = (temp - norm) - y;
            norm = temp;
        }}
      }
      break;
    case 3:
      {
        auto v = f.template get_view<Real***,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              y = std::pow(v(i,j,k),2) - c;
              temp = norm + y;
              c = (temp - norm) - y;
              norm = temp;
        }}}
      }
      break;
    case 4:
      {
        auto v = f.template get_view<Real****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                y = std::pow(v(i,j,k,l),2) - c;
                temp = norm + y;
                c = (temp - norm) - y;
                norm = temp;
        }}}}
      }
      break;
    case 5:
      {
        auto v = f.template get_view<Real*****,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  y = std::pow(v(i,j,k,l,m),2) - c;
                  temp = norm + y;
                  c = (temp - norm) - y;
                  norm = temp;
        }}}}}
      }
      break;
    case 6:
      {
        auto v = f.template get_view<Real******,Host>();
        for (int i=0; i<v.extent_int(0); ++i) {
          for (int j=0; j<v.extent_int(1); ++j) {
            for (int k=0; k<v.extent_int(2); ++k) {
              for (int l=0; l<v.extent_int(3); ++l) {
                for (int m=0; m<v.extent_int(4); ++m) {
                  for (int n=0; n<v.extent_int(4); ++n) {
                    y = std::pow(v(i,j,k,l,m,n),2) - c;
                    temp = norm + y;
                    c = (temp - norm) - y;
                    norm = temp;
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  return std::sqrt(norm);
}


} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
