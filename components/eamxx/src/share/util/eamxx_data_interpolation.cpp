#include "share/util/eamxx_data_interpolation.hpp"

#include <filesystem>
#include <regex>

namespace scream{

DataInterpolation::
DataInterpolation (const std::shared_ptr<const AbstractGrid>& grid,
                   const ekat::ParameterList& params)
 : m_params (params)
{
  namespace fs = std::filesystem;

  auto path = m_params.get<std::string>("input_files_path","./");
  EKAT_REQUIRE_MSG (fs::exists(path),
      "Error! Input value for 'input_files_path' does not exist.\n"
      " - path: " + path + "\n");
  EKAT_REQUIRE_MSG (fs::is_directory(path),
      "Error! Input value for 'input_files_path' is not a directory.\n"
      " - path: " + path + "\n");

  auto dir_perm_rx = [](const std::string& path) {
    for (const auto& entry : fs::directory_iterator(path) {
      return true; // If you can iterate, path has rx permissions
    }
    return false;
  };
  EKAT_REQUIRE_MSG (dir_perm_rx(path),
    "Error! Input files path does not have correct permissions.\n"
      " - path: " + path + "\n");

  auto file_readable = [] (const fs::path& filePath) {
    std::ifstream file(filePath.string());
    return file.good(); // Check if the file can be opened
  }

  // In order to form full path filenames
  if (path.back()!='/') path += '/';

  if (m_params.isParameter("input_files_names")) {
    EKAT_REQUIRE_MSG (not m_params.isParameter("input_files_pattern"),
        "Error! Cannot provide both 'input_file_pattern' and 'input_file_name');

    m_input_files = m_params.get<strvec_t>("input_files_names");
    for (auto& f : m_input_files) {
      f = path + f;
    }
  } else if (m_params.isParameter("input_files_pattern")) {
    EKAT_REQUIRE_MSG (not m_params.isParameter("input_file_name"),
        "Error! Cannot provide both 'input_file_pattern' and 'input_file_name');

    const auto& pattern = m_params.get<std::string>("input_files_pattern");
    std::regex pattern (p);
    for (const auto& entry : fs::directory_iterator(path)) {
      if (not entry.is_regular_file())
        continue;

      std::string filename = entry.path().filename().string();

      if (std::regex_match(filename,pattern)) {
        m_input_files.push_back(path+filename);
      }
    }
  } else {
    EKAT_ERROR_MSG ("Error! Missing 'input_files_names' or 'input_files_patterns'.\n");
  }

  for (const auto& f : m_input_files) {
    EKAT_REQUIRE_MSG (file_readable(entry.path()),
        "Error! One of the input files is not readable.\n"
        " - file   : " + f + "\n");
  }
}

void DataInterpolation::
complete_setup (const util::TimeStamp& t0)
{
  m_time_state.end = m_time_state.beg = t0;
  update_end_fields ();
}

void DataInterpolation::run (const util::TimeStamp& ts)
{
  if (m_time_state.end<ts) {
    std::swap(m_time_state.beg,m_time_state.end);
    std::swap(m_fields_beg,m_fields_end);
    update_end_fields ();
  }
}

void DataInterpolation::update_end_fields ()
{

}

void DataInterpolation::set_field (const Field& f)
{
  m_tgt_fields.push_back(f);
}

std::vector<Field>&
DataInterpolation::get_fields (Phase phase, bool beg)
{
  if (beg) {
    return m_fields_beg[phase];
  } else {
    return m_fields_end[phase];
  }
}

} // namespace scream
