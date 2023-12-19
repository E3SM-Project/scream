#include "eamxx_default_model_init.hpp"

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

  AtmosphereInput ic_reader (filename,grid,ic_fields);
  ic_reader.set_logger(logger);
  ic_reader.read_variables();
  ic_reader.finalize();

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

  auto phys_grid = m_gm->get_grid("Physics");
  auto ic_grid   = phys_grid->clone("Physics",true);
  ic_grid->reset_field_tag_name(COL,"ncol_d");
  for (auto n : phys_grid->aliases()) {
    ic_grid->add_alias(n);
  }

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

  // Load fields
  const auto& filename = m_params.get<std::string>("topography_filename");
  AtmosphereInput ic_reader (filename,ic_grid,fields);
  ic_reader.set_logger(logger);
  ic_reader.read_variables();
  ic_reader.finalize();

  // Update time stamp of ic fields
  for (auto& f : fields) {
    f.get_header().get_tracking().update_time_stamp(m_run_t0);
  }
  logger->info("    [EAMxx] Reading topography file ... done!");
}

} // namespace scream
