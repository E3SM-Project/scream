#include "simple_op.hpp"

#include "share/util/scream_array_utils.hpp"

namespace scream {

SimpleOp::
SimpleOp  (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic (comm,params)
{
  m_name = m_params.get<std::string>("diag_name",m_params.name());

  auto op = m_params.get<std::string>("op");
  if (op=="add") {
    m_op = ADD;
  } else if (op=="sub") {
    m_op = SUB;
  } else if (op=="prod") {
    m_op = PROD;
  } else if (op=="div") {
    m_op = DIV;
  } else if (op=="neg") {
    m_op = NEG;
  } else {
    EKAT_ERROR_MSG ("Error! Invalid value for 'op' in SimpleOp constructor.\n"
        "  - input value: '" + op + "'\n"
        "  - valid values: 'add', 'sub', 'prod', 'div', 'neg'\n");
  }

  m_f1_scale = m_params.get<Real>("f1_scale",1);

  // Note: you may ask why we can't just store f and be done, rather than go through the
  // add_field infrastructure. Unfortunately, there are some checks in the base classes
  // that require the diagnostic to have 1+ required fields. So we have to do this.
  // TODO: one day we may make atm diags *not* inherit from atm process...
  EKAT_REQUIRE_MSG (m_params.isParameter("f1"),
      "Error! BinaryOp constructor requires 'f1' parameter.\n");
  auto f1 = m_params.get<Field>("f1");
  add_field<Required>(f1.get_header().get_identifier());

  if (m_op!=NEG) {
    auto f2_type = m_params.get<std::string>("f2_type","field");
    if (f2_type=="field") {
      m_f2_type = FIELD;
      m_f2_scale = m_params.get<Real>("f2_scale",1);
      EKAT_REQUIRE_MSG (m_params.isParameter("f2"),
          "Error! BinaryOp constructor requires 'f2' parameter if 'f2_type' is 'field'.\n");
      auto f2 = m_params.get<Field>("f2");
      add_field<Required>(f2.get_header().get_identifier());
    } else if (f2_type=="value") {
      m_f2_type = VALUE;
      m_f2_value = m_params.get<Real>("f2_value");
    } else {
      EKAT_ERROR_MSG (
          "Error! Invalid choice for f2_type in BinaryOp constructor.\n"
          "  - input value: " + f2_type + "\n"
          "  - valid values: 'field', 'value'\n");
    }
  }
}

void SimpleOp::initialize_impl (const RunType /*run_type*/)
{
  // Run some checks
  const auto& fs = get_fields_in();
  EKAT_REQUIRE_MSG (fs.size()==2 or
      (fs.size()==1 and (m_op==NEG or m_f2_type==VALUE)),
      "Error! Wrong number of input fields for SimpleOp diagnostic.\n"
      "  - diag name: " + m_name + "\n"
      "  - simple op: " + op_str(m_op) + "\n"
      "  - f2 type: " + f2type_str(m_f2_type) + "\n"
      "  - num input fields: " + std::to_string(fs.size()) + "\n");

  auto f1 = fs.front();
  EKAT_REQUIRE_MSG (f1.rank()<=3,
      "Error! SimpleOp only supports fields up to rank 3.\n"
      "  - input fields rank: " + std::to_string(f1.rank()) + "\n");
  EKAT_REQUIRE_MSG (f1.data_type()==DataType::RealType,
      "Error! SimpleOp only implemented for RealType data type.\n"
      "  - input field: " + f1.name() + "\n"
      "  - data type  : " + e2str(f1.data_type()) + "\n");

  if (m_f2_type==FIELD) {
    auto f2 = fs.back();
    EKAT_REQUIRE_MSG (f1.get_header().get_identifier().get_layout()==f2.get_header().get_identifier().get_layout(),
        "Error! SimpleOp requires compatible layout input fields.\n"
        "  - f1 name: " + f1.name() + "\n"
        "  - f2 name: " + f1.name() + "\n"
        "  - f1 layout: " + to_string(f1.get_header().get_identifier().get_layout()) + "\n"
        "  - f2 layout: " + to_string(f1.get_header().get_identifier().get_layout()) + "\n");

    EKAT_REQUIRE_MSG (f1.data_type()==f2.data_type(),
        "Error! SimpleOp requires input fields with same data type.\n"
        "  - f1 name: " + f1.name() + "\n"
        "  - f2 name: " + f1.name() + "\n"
        "  - f1 data type: " + e2str(f1.data_type()) + "\n"
        "  - f2 data type: " + e2str(f1.data_type()) + "\n");
  }

  m_diagnostic_output = f1.clone(m_name);
}

void SimpleOp::compute_diagnostic_impl ()
{
  const auto f1 = get_fields_in().front();

  if (m_op==NEG) {
    m_diagnostic_output.scale(Real(-1));
    return;
  } else if (m_f2_type==VALUE) {
    if (m_op==PROD or m_op==DIV) {
      Real alpha = m_op==PROD ? m_f1_scale*m_f2_value : m_f1_scale/m_f2_value;
      m_diagnostic_output.update(f1,alpha,Real(0));
      return;
    } else if (m_op==ADD or m_op==SUB) {
      Real val = m_op==ADD ? 1 : -1;
      m_diagnostic_output.deep_copy(val*m_f2_value);
      m_diagnostic_output.update(f1,m_f1_scale,Real(1));
      return;
    } else {
      // Should never reach this, but just in case we forget to add an impl
      EKAT_ERROR_MSG ("Error! Unexpected value for m_op in BinaryOp::compute_diagnostic_impl\n");
    }
  }

  using KT = KokkosTypes<DefaultDevice>;

  int n = m_diagnostic_output.get_header().get_identifier().get_layout().size();
  typename KT::RangePolicy policy (0,n);

  const auto f2 = get_fields_in().back();
  const auto rank = f1.rank();

  using view_0d = decltype(f1.get_view<const Real>());
  using view_1d = decltype(f1.get_view<const Real*>());
  using view_2d = decltype(f1.get_view<const Real**>());
  using view_3d = decltype(f1.get_view<const Real***>());

  view_0d v1_0d,v2_0d;
  view_1d v1_1d,v2_1d;
  view_2d v1_2d,v2_2d;
  view_3d v1_3d,v2_3d;
  view_0d::non_const_type d_0d;
  view_1d::non_const_type d_1d;
  view_2d::non_const_type d_2d;
  view_3d::non_const_type d_3d;
  switch (rank) {
    case 0:
      v1_0d = f1.get_view<const Real>();
      v2_0d = f2.get_view<const Real>();
      d_0d  = m_diagnostic_output.get_view<Real>();
      break;
    case 1:
      v1_1d = f1.get_view<const Real*>();
      v2_1d = f2.get_view<const Real*>();
      d_1d  = m_diagnostic_output.get_view<Real*>();
      break;
    case 2:
      v1_2d = f1.get_view<const Real**>();
      v2_2d = f2.get_view<const Real**>();
      d_2d  = m_diagnostic_output.get_view<Real**>();
      break;
    case 3:
      v1_3d = f1.get_view<const Real***>();
      v2_3d = f2.get_view<const Real***>();
      d_3d  = m_diagnostic_output.get_view<Real***>();
      break;
  }

  auto op = m_op;
  auto extents = f1.get_header().get_identifier().get_layout().extents();
  auto scale1 = m_f1_scale;
  auto scale2 = m_f2_scale;
  Kokkos::parallel_for (m_name,policy,
                        KOKKOS_LAMBDA(const int idx) {
    int i,j,k;
    if (rank==2) {
      unflatten_idx(idx,extents,i,j);
    } else if (rank==3) {
      unflatten_idx(idx,extents,i,j,k);
    }
    const auto v1 = rank==0 ? v1_0d()
                            :rank==1 ? v1_1d(idx)
                                     : rank==2 ? v1_2d(i,j)
                                               : v1_3d(i,j,k);
    const auto v2 = rank==0 ? v2_0d()
                            :rank==1 ? v2_1d(idx)
                                     : rank==2 ? v2_2d(i,j)
                                               : v2_3d(i,j,k);
          auto& d = rank==0 ? d_0d()
                            :rank==1 ? d_1d(idx)
                                     : rank==2 ? d_2d(i,j)
                                               : d_3d(i,j,k);
    switch (op) {
      case ADD:
        d = scale1*v1 + scale2*v2; break;
      case SUB:
        d = scale1*v1 - scale2*v2; break;
      case PROD:
        d = scale1*v1 * scale2*v2; break;
      case DIV:
        d = scale1*v1 / scale2*v2; break;
      default:
        EKAT_KERNEL_ERROR_MSG (
            "Error! Unexpected op type in BinaryOp::compute_diagnostic_impl kernel\n");
    }
  });
}

} // namespace scream
