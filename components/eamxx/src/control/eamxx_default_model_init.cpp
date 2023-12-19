#include "eamxx_default_model_init.hpp"
#include "intensive_observation_period.hpp"

#include "share/io/scorpio_input.hpp"

namespace scream
{

DefaultModelInit::
DefaultModelInit (const strmap_t<Field>& eamxx_inputs,
                  const std::shared_ptr<const GridsManager>& gm,
                  const ekat::ParameterList& ic_pl,
                  const util::TimeStamp& run_t0)
 : ModelInit (eamxx_inputs,gm,ic_pl,run_t0)
{
  // Nothing to do here
}

void DefaultModelInit::
set_initial_conditions (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  set_constant_fields (logger);
  if (m_params.isParameter("Filename") and m_ic_fields.size()>0) {
    read_ic_file (logger);
  }
  if (m_params.isParameter("topography_filename") and m_topo_fields.size()>0) {
    read_topo_file (logger);
  }
}

void DefaultModelInit::
read_ic_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  const auto& filename = m_params.get<std::string>("Filename");

  logger->info("    [EAMxx] Reading initial conditions file ...");
  logger->info("        filename: " + filename);

  // Use convert std::set to a std::vector
  std::vector<Field> ic_fields;
  for (const auto& f : m_ic_fields) {
    ic_fields.push_back(f);
  }
  auto grid = m_gm->get_grid("Physics");

  auto iop = m_params.get<std::shared_ptr<control::IntensiveObservationPeriod>>("iop");

  if (iop==nullptr) {
    AtmosphereInput ic_reader (filename,grid,ic_fields);
    ic_reader.set_logger(logger);
    ic_reader.read_variables();
    ic_reader.finalize();
  } else {
    // Setup io grid
    iop->setup_io_info(filename,grid);

    // Copy data to the closest lat/lon col in the io data to every other col
    iop->read_fields_from_file_for_iop (filename,
                                          ic_fields,
                                          m_run_t0);

    // Overwrite IC data with data from iop file
    iop->read_iop_file_data(m_run_t0);
    iop->set_fields_from_iop_data(m_eamxx_inputs);
  }

  // Update time stamp of ic fields
  for (auto& f : ic_fields) {
    f.get_header().get_tracking().update_time_stamp(m_run_t0);
  }

  logger->info("    [EAMxx] Reading initial conditions file ... done!");
}

void DefaultModelInit::
read_topo_file (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  using namespace ShortFieldTagsNames;
  const auto& filename = m_params.get<std::string>("topography_filename");
  logger->info("    [EAMxx] Reading topography file ...");
  logger->info("        filename: " + filename);

  // Store topo fields as std::vector. Alias name, to match what's in the file
  std::vector<Field> fields;
  for (const auto& f : m_topo_fields) {
    if (f.name()=="phis") {
      fields.push_back(f.alias("PHIS_d"));
    } else if (f.name()=="sgh30") {
      EKAT_ERROR_MSG(
        "[DefaultModelInit] Error! Requesting sgh30 field from topo file.\n"
        "  The DefaultModelInit class does not handle Physics PG2 grid,\n"
        "  but sgh30 should only be needed for Physics PG2 grid.\n");
    } else {
      EKAT_ERROR_MSG(
        "[DefaultModelInit] Error! Unexpected/unrecognized topo field name.\n"
        " - field name: " + f.name() + "\n");
    }
  }

  auto phys_grid = m_gm->get_grid("Physics");
  auto iop = m_params.get<std::shared_ptr<control::IntensiveObservationPeriod>>("iop");
  if (iop==nullptr) {
    // The topo file stores colum dim as ncol for PG2 and ncol_d for GLL.
    // So create a shallow clone of the gll grid, with different name for COL tag
    auto io_grid   = phys_grid->clone("Physics",true);
    io_grid->reset_field_tag_name(COL,"ncol_d");
    for (auto n : phys_grid->aliases()) {
      io_grid->add_alias(n);
    }

    // Load fields
    AtmosphereInput ic_reader (filename,io_grid,fields);
    ic_reader.set_logger(logger);
    ic_reader.read_variables();
    ic_reader.finalize();
  } else {
    iop->setup_io_info(filename,phys_grid);
    iop->read_fields_from_file_for_iop (filename,fields,m_run_t0);
  }

  // Update time stamp of ic fields
  for (auto& f : fields) {
    f.get_header().get_tracking().update_time_stamp(m_run_t0);
  }
  logger->info("    [EAMxx] Reading topography file ... done!");
}

} // namespace scream
