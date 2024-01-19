#include "eamxx_fvphys_model_init.hpp"

#include "eamxx_homme_fv_phys_helper.hpp"
#include "physics_dynamics_remapper.hpp"
#include "share/io/scorpio_input.hpp"

// Homme includes
#include "Context.hpp"
#include "Tracers.hpp"
#include "ElementsState.hpp"

namespace scream
{

FvPhysModelInit::
FvPhysModelInit (const strmap_t<Field>& eamxx_inputs,
                 const std::shared_ptr<const GridsManager>& gm,
                 const ekat::ParameterList& ic_pl,
                 const util::TimeStamp& run_t0)
 : ModelInit (eamxx_inputs,gm,ic_pl,run_t0)
{
  // FvPhys requires certain fields to be present
  for (const std::string& n : {"T_mid","horiz_winds","ps","phis","pseudo_density","tracers"}) {
    EKAT_REQUIRE_MSG (m_eamxx_inputs.count(n)==1,
        "Error! Using FvPhysModelInit but '" + n + "' is not a model input.\n");
  }
}

void FvPhysModelInit::
set_initial_conditions (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  set_constant_fields (logger);
  if (m_params.isParameter("Filename") and m_ic_fields.size()>0) {
    read_ic_file (logger);
  }
  if (m_params.isParameter("topography_filename") and m_topo_fields.size()>0) {
    read_topo_file (logger);
  }

  remap_fields ();
}

void FvPhysModelInit::
read_ic_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  const auto& filename = m_params.get<std::string>("Filename");

  logger->info("    [EAMxx] Reading initial conditions file ...");
  logger->info("        filename: " + filename);

  using namespace ShortFieldTagsNames;

  std::vector<Field> nc_fields;
  const auto cgll_grid = HommeFvPhysHelper::instance().m_cgll_grid;
  std::cout << "fvphys, reading with grid " << cgll_grid->name() << "\n";
  for (const auto& f : m_ic_fields) {

    const auto& fid = f.get_header().get_identifier();
    const auto& layout = fid.get_layout();
    const auto lt = get_layout_type(layout.tags());
          auto gll_layout = FieldLayout::invalid();
    const bool midpoints = layout.tags().back()==LEV;
    const int vec_dim = layout.is_vector_layout() ? layout.get_vector_dim() : -1;
    const int vec_extent = vec_dim>=0 ? layout.dims()[vec_dim] : 0;
    switch (lt) {
      case LayoutType::Scalar2D:
        gll_layout = cgll_grid->get_2d_scalar_layout();
        break;
      case LayoutType::Vector2D:
        gll_layout = cgll_grid->get_2d_vector_layout(CMP,vec_extent);
        break;
      case LayoutType::Scalar3D:
        gll_layout = cgll_grid->get_3d_scalar_layout(midpoints);
        break;
      case LayoutType::Vector3D:
        gll_layout = cgll_grid->get_3d_vector_layout(midpoints,CMP,vec_extent);
        break;
      default:
        EKAT_ERROR_MSG ("Unexpected/unsupported layout for EAMxx input.\n"
                        " - field name: " + fid.name() + "\n"
                        " - layout    : " + to_string(layout) + "\n");
    }

    FieldIdentifier gll_fid(fid.name(),gll_layout,fid.get_units(),cgll_grid->name());
    Field gll_f (gll_fid);
    const int pack_size = f.get_header().get_alloc_properties().get_largest_pack_size();
    gll_f.get_header().get_alloc_properties().request_allocation(pack_size);
    gll_f.allocate_view();
    m_gll_fields[f.name()] = gll_f;

    if (f.name()=="tracers") {
      // We do need to pass the tracers to the fvphys remapper, but we won't
      // find them in the IC file. So don't add it to nc_fields.
      continue;
    }

    // Also check whether this field was already init-ed to a constant.
    // If yes, simply copy the constant in the gll field, otherwise add
    // the gll field to the vector of fields to read in
    if (m_const_fields_values.count(f.name())==1) {
      const auto& fl = f.get_header().get_identifier().get_layout();
      const auto& val_str = m_const_fields_values.at(f.name());
      const bool vector_values = val_str.front()=='(' and val_str.back()==')';
      if (vector_values) {
        const auto vals_str = ekat::split(val_str.substr(0,val_str.size()-1),",");
        const int ncmp = vals_str.size();
        const int vec_dim = fl.get_vector_dim();
        for (int icmp=0; icmp<ncmp; ++icmp) {
          auto comp = f.subfield(vec_dim,icmp);
          comp.deep_copy(std::stod(vals_str[icmp]));
        }
      } else {
        gll_f.deep_copy(std::stod(val_str));
      }
    } else {
      nc_fields.push_back(gll_f);
    }
  }

  AtmosphereInput ic_reader (filename,cgll_grid,nc_fields);
  ic_reader.set_logger(logger);
  ic_reader.read_variables();
  ic_reader.finalize();

  logger->info("    [EAMxx] Reading initial conditions file ... done!");
}

void FvPhysModelInit::
read_topo_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  using namespace ShortFieldTagsNames;

  const auto& filename = m_params.get<std::string>("topography_filename");
  logger->info("    [EAMxx] Reading topography file ...");
  logger->info("        filename: " + filename);

  using namespace ShortFieldTagsNames;

  const auto cgll_grid = HommeFvPhysHelper::instance().m_cgll_grid;
  const auto phys_grid = HommeFvPhysHelper::instance().m_phys_grid;

  for (const auto& f : m_topo_fields) {
    if (f.name()=="phis") {
      // We only read phis on GLL grid. But the ncol dim for gll is 'ncol_d'
      // in the topo file, so we need to create a soft copy of cgll grid,
      // with a new name for the COL dim.
      auto io_grid = cgll_grid->clone(cgll_grid->name(),true);
      io_grid->reset_field_tag_name(COL,"ncol_d");
      const auto& fid = f.get_header().get_identifier();
      const auto layout = io_grid->get_2d_scalar_layout();
      FieldIdentifier phis_fid("PHIS_d",layout,fid.get_units(),io_grid->name());
      Field phis (phis_fid);
      phis.allocate_view();
      m_gll_fields["phis"] = phis;

      // Also check whether phis was already init-ed to a constant.
      // If yes, simply copy the constant in the gll field, otherwise
      // proceed to read from file
      if (m_const_fields_values.count(f.name())==1) {
        const auto& val_str = m_const_fields_values.at(f.name());
        phis.deep_copy(std::stod(val_str));
      } else {
        std::vector<Field> nc_fields = {phis};
        AtmosphereInput topo_reader (filename,io_grid,nc_fields);
        topo_reader.set_logger(logger);
        topo_reader.read_variables();
        topo_reader.finalize();
      }
    } else if (f.name()=="sgh30") {
      // Only read sgh30 on FV grid
      EKAT_REQUIRE_MSG (f.get_header().get_identifier().get_grid_name()==phys_grid->name(),
          "Error! Somehow, the field 'sgh30' is required on non-PG2 grid.\n");

      // If someone inits sgh30 to a constant, skip this part
      if (m_const_fields_names.count(f.name())==0) {
        std::vector<Field> nc_fields = {f.alias("SGH30") };

        AtmosphereInput topo_reader (filename,phys_grid,nc_fields);
        topo_reader.set_logger(logger);
        topo_reader.read_variables();
        topo_reader.finalize();
      }
    } else {
      EKAT_ERROR_MSG ("Unexpected topology-related field.\n"
          " - field name : " + f.name() + "\n"
          " - valid names: phis, sgh30\n");
    }
  }

  logger->info("    [EAMxx] Reading topography file ... done!");
}

void FvPhysModelInit::
remap_fields ()
{
  using kt = KokkosTypes<DefaultDevice>;
  using view_1d = typename kt::view_ND<Real,1>;
  using view_2d = typename kt::view_ND<Real,2>;
  using view_3d = typename kt::view_ND<Real,3>;
  using view_4d = typename kt::view_ND<Real,4>;
  using view_5d = typename kt::view_ND<Real,5>;
  using view_6d = typename kt::view_ND<Real,6>;

  auto as_real = [] (Homme::Scalar* p) {
    return reinterpret_cast<Real*>(p);
  };

  auto& fv_phys = HommeFvPhysHelper::instance();
  auto& c = Homme::Context::singleton();
  auto& state = c.get<Homme::ElementsState>();
  auto& tracers = c.get<Homme::Tracers>();

  const int np = HOMMEXX_NGP;
  const int nelem = m_dyn_grid->get_num_local_dofs() / (np*np);
  const int homme_num_lev  = HOMMEXX_PACK_SIZE * HOMMEXX_NUM_LEV;
  const int homme_num_ilev = HOMMEXX_PACK_SIZE * HOMMEXX_NUM_LEV_P;

  // Create dyn grid fields from Homme's views
  auto vTheta = view_5d(as_real(state.m_vtheta_dp.data()),

  // First, remap the fields from cgll to dyn grid
  PhysicsDynamicsRemapper gll2dyn(fv_phys.m_cgll_grid,fv_phys.m_dyn_grid);
  gll2dyn.registration_begins();
  for (const std::string& n : {"T_mid","horiz_winds","ps","phis","pseudo_density","tracers"}) {
    gll2dyn.register_field(m_gll_fields.at(n),fv_phys.m_dyn_fields.at(n));
  }
  gll2dyn.registration_ends();
  gll2dyn.remap(/*forward = */ true);

  // Second, use FvPhysHelper to remap dyn->phys
}

} // namespace scream
