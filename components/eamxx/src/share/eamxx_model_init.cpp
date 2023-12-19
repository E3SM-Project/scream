#include "eamxx_model_init.hpp"

namespace scream
{

ModelInit::
ModelInit (const strmap_t<Field>& eamxx_inputs,
           const std::shared_ptr<const GridsManager>& gm,
           const ekat::ParameterList& ic_pl,
           const util::TimeStamp& run_t0)
 : m_params       (ic_pl)
 , m_gm           (gm)
 , m_run_t0       (run_t0)
 , m_eamxx_inputs (eamxx_inputs)
{
  separate_inputs ();
}

void ModelInit::
set_constant_fields (const std::shared_ptr<ekat::logger::LoggerBase>& logger)
{
  logger->info("    [EAMxx] Initializing constant fields ...");

  auto str2double = [&] (const std::string& fname, const std::string& val_str) {
    double val;
    try {
      val = std::stod (val_str);
    } catch (std::exception&) {
      EKAT_ERROR_MSG (
        "Error! Bad value specification for constant field initialization.\n"
        " Field name: " + fname + "\n"
        "Could not convert the string '" + val_str + "' to double.\n");
    }
    return val;
  };

  const auto& const_fields = m_params.get<strvec_t>("constant_fields",{});
  for (const auto& s : const_fields) {
    logger->info("      " + s);

    // Split input string FNAME = VALUE(S)
    auto tokens = ekat::split(s,"=");
    EKAT_REQUIRE_MSG (tokens.size()==2,
        "Error! Bad syntax specifying constant input fields.\n"
        "  Valid syntax: 'field_name = field_value(s)'\n"
        "  Input string: '" + s + "'\n");
    const auto fname   = ekat::trim(tokens[0]);
    const auto val_str = ekat::trim(tokens[1]);
    EKAT_REQUIRE_MSG (fname.size()>0,
        "Error! Bad syntax specifying constant input fields (empty field name).\n"
        " Input string: '" + s + "'\n");
    EKAT_REQUIRE_MSG (val_str.size()>0,
        "Error! Bad syntax specifying constant input fields (empty value).\n"
        " Input string: '" + s + "'\n");

    auto get_fname = [&] (const strmap_t<Field>::value_type& it) {
      return it.first;
    };
    EKAT_REQUIRE_MSG (m_eamxx_inputs.count(fname)==1,
        "Error! Attempt to prescribe constant value for a field that is not an atm input.\n"
        " - field name: " + fname + "\n"
        " - atm inputs: " + ekat::join(m_eamxx_inputs,get_fname,",") + "\n");

    // Get the field, and ensure it's in the ic fields list
    auto& f = m_eamxx_inputs.at(fname);

    // Check if we are attempting to init a vec field with f=(v1,v2,...,vN) syntax
    bool vector_values = val_str.front()=='(' and val_str.back()==')';
    if (vector_values) {
      const auto& fl = f.get_header().get_identifier().get_layout();
      EKAT_REQUIRE_MSG (fl.is_vector_layout(),
          "Error! Attempt to assing a vector of values to a non-vector field.\n"
          "  - field name: " + fname + "\n"
          "  - field layout: " + to_string(fl) + "\n"
          "  - input string: " + s + "\n");
      const auto vals = ekat::split(val_str.substr(1,val_str.size()-2),",");
      const int vec_dim = fl.get_vector_dim();
      EKAT_REQUIRE_MSG (static_cast<int>(vals.size())==fl.dim(vec_dim),
          "Error! Vector of values does not match vector field extent.\n"
          "  - field name: " + fname + "\n"
          "  - field layout: " + to_string(fl) + "\n"
          "  - input string: " + s + "\n");
      for (int icmp=0; icmp<fl.dim(vec_dim); ++icmp) {
        auto f_i = f.subfield(vec_dim,icmp);
        auto val = str2double(fname,vals[icmp]);
        f_i.deep_copy(val);
      }
    } else {
      // Not a vector assignment, so just init the whole field
      // NOTE: if user had a typo, like "f=(v1,v2", notice missing ')',
      //       this fcn will throw, cause "(v1,v2" cannot be converted to double
      auto val = str2double(fname,val_str);
      f.deep_copy(val);
    }
    
    f.get_header().get_tracking().update_time_stamp(m_run_t0);

    // Move this field out of m_ic_fields/m_topo_fields;
    if (m_ic_fields.count(f)==1) {
      m_ic_fields.erase(f);
      m_ic_fields_names.erase(ekat::find(m_ic_fields_names,f.name()));
    }
    if (m_topo_fields.count(f)==1) {
      m_topo_fields.erase(f);
      m_topo_fields_names.erase(ekat::find(m_topo_fields_names,f.name()));
    }

    m_const_fields.insert(f);
  }
  logger->info("    [EAMxx] Initializing constant fields ... done!");
}

void ModelInit::
separate_inputs ()
{
  // Separate inputs between topo-related fields and non-topo-related fields
  for (const auto& it : m_eamxx_inputs) {
    // const auto& fid = f.get_header().get_identifier();
    // const auto& gname = fid.get_grid_name();
    const auto& name = it.first;
    const auto& f    = it.second;

    if (name=="phis" or name=="sgh30") {
      // Topography fields are set in a separate vector
      // Also, keep track of the field name in the topo file vs eamxx
      m_topo_fields.insert(f);
      m_topo_fields_names.push_back(name);
      // if (fid.name()=="phis") {
      //   m_topo_fields_names_eamxx[gname].push_back("phis");
      //   if (gname=="Physics PG2") {
      //     m_topo_fields_names_file[gname].push_back("PHIS");
      //   } else {
      //     m_topo_fields_names_file[gname].push_back("PHIS_d");
      //   }
      // } else {
      //   EKAT_ASSERT_MSG(gname == "Physics PG2",
      //       "Error! Requesting sgh30 field on " + gname +
      //       ", but topo file only has sgh30 for Physics PG2.\n");
      //   m_topo_fields_names_file[gname].push_back("SGH30");
      //   m_topo_fields_names_eamxx[gname].push_back("sgh30");
      // }
    // } else if (m_ic_grid->is_valid_alias(gname)) {
    } else {
      // We consider the field as an input only if it's on the IC grid
      m_ic_fields.insert(f);
      m_ic_fields_names.push_back(name);
    // } else {
    //   EKAT_ERROR_MSG (
    //       "Error! Unexpected grid for non-topography-related eamxx input field.\n"
    //       " - field name: " + fid.name() + "\n"
    //       " - field grid: " + gname + "\n"
    //       " - ic grid aliases: " + ekat::join(m_ic_grid->aliases(),",") + "\n");
    }
  };
}

} // namespace scream
