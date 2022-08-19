#include "share/io/scream_io_utils.hpp"
#include "share/util/scream_utils.hpp"

#include <fstream>

namespace scream {

std::string find_filename_in_rpointer (
    const std::string& casename,
    const bool model_restart,
    const ekat::Comm& comm,
    const util::TimeStamp& run_t0)
{
  std::string filename;
  bool found = false;
  std::string content;
  std::string suffix = model_restart ? ".r." : ".rhist.";
  if (comm.am_i_root()) {
    std::ifstream rpointer_file;
    std::string line;
    rpointer_file.open("rpointer.atm");

    // If the timestamp is in the filename, then the filename ends with "S.nc",
    // with S being the string representation of the timestamp
    auto extract_ts = [&] (const std::string& line) -> util::TimeStamp {
      auto ts_len = run_t0.to_string().size();
      auto min_size = ts_len+3;
      if (line.size()>=min_size) {
        auto ts_str = line.substr(line.size()-min_size);
        auto ts = util::str_to_time_stamp(ts_str);
        return ts;
      } else {
        return util::TimeStamp();
      }
    };

    // Note: keep swallowing line, even after the first match, since we never wipe
    //       rpointer.atm, so it might contain multiple matches, and we want to pick
    //       the last one (which is the last restart file that was written).
    while (rpointer_file >> line) {
      content += line + "\n";

      if (line.find(casename) != std::string::npos && line.find(suffix) != std::string::npos) {
        // Extra check: make sure the date in the filename (if present) precedes this->t0.
        // If there's no time stamp in one of the filenames, we assume this is some sort of
        // unit test, with no time stamp in the filename, and we accept the filename.
        auto ts = extract_ts(line);
        if (not ts.is_valid() || !(ts<=run_t0) ) {
          found = true;
          filename = line;
        }
      }
    }
  }
  int ifound = int(found);
  comm.broadcast(&ifound,1,0);
  found = bool(ifound);
  broadcast_string(content,comm,comm.root_rank());

  // If the history restart file is not found, it must be because the last
  // model restart step coincided with a model output step, in which case
  // a restart history file is not written.
  // If that's the case, *disable* output restart, by setting
  //   'Restart'->'Perform Restart' = false
  // in the input parameter list
  EKAT_REQUIRE_MSG (found,
      "Error! Restart requested, but no restart file found in 'rpointer.atm'.\n"
      "   restart case name: " + casename + "\n"
      "   restart file type: " + std::string(model_restart ? "model restart" : "history restart") + "\n"
      "   rpointer content:\n" + content);

  // Have the root rank communicate the nc filename
  broadcast_string(filename,comm,comm.root_rank());

  return filename;
}

std::string get_nc_tag_name (const FieldTag& t, const int extent) {
  using namespace ShortFieldTagsNames;

  std::string name = "";
  switch(t) {
    case EL:
      name = "elem";
      break;
    case LEV:
      name = "lev";
      break;
    case ILEV:
      name = "ilev";
      break;
    case TL:
      name = "tl";
      break;
    case COL:
      name = "ncol";
      break;
    case GP:
      name = "gp";
      break;
    case CMP:
      name = "dim" + std::to_string(extent);
      break;
    // Added for rrtmgp - TODO revisit this paradigm, see comment in field_tag.hpp
    case NGAS:
      name = "ngas";
      break;
    case SWBND:
      name = "swband";
      break;
    case LWBND:
      name = "lwband";
      break;
    default:
      EKAT_ERROR_MSG("Error! Field tag not supported in netcdf files.");
  }

  return name;
}

std::string get_io_tag_name (const std::shared_ptr<const AbstractGrid>& io_grid,
                             const FieldTag tag, const int extent) {
  // Check if this is a partitioned dim AND io_grid has a special name
  // for it, otherwise use general utility to get the name
  const auto is_partitioned = io_grid->get_partitioned_dim_tag()==tag;
  return is_partitioned && io_grid->has_partitioned_dim_name()
        ? io_grid->get_partitioned_dim_name()
        : get_nc_tag_name(tag,extent);
}

} // namespace scream
