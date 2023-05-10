#include "share/util/scream_time_interpolation.hpp"

namespace scream {
namespace util {

// Constructors:
TimeInterpolation::TimeInterpolation(const grid_ptr_type& grid)
{
  m_grid = grid;

  // Construct field manager to use with ScorpioInput interface
  m_fm_t0 = std::make_shared<FieldManager>(m_grid);
  m_fm_t1 = std::make_shared<FieldManager>(m_grid);
  m_fm_t0->registration_begins();
  m_fm_t0->registration_ends();
  m_fm_t1->registration_begins();
  m_fm_t1->registration_ends();

}

/* ----------------------------------------------------------------
 * At initialization we organize all of the internal structures that
 * will be used for time interpolation.
 * Inputs:
 *   1. list_of_variables -
 *      A vector of strings that lists each variable that we intend
 *      to interpolate.  This will be used to set up the interpolation
 *      data structures.
 *   2. list_of_files -
 *      A vector of strings that lists each file that contains
 *      interpolation data.  This will be used to load the individual
 *      time stamps of data.
 * We only need to store two timestamps of data at any given step.  So we don't
 * load all of the data at once.  Instead, we store an internal vector of timestamps
 * that will be used to figure out which file contains the needed data.  During initialization
 * we query all files and populate a map of timestamp to filename.
 * ---------------------------------------------------------------- */

void TimeInterpolation::add_field(const Field& field)
{
  EKAT_REQUIRE_MSG(std::find(m_fields_names.begin(),m_fields_names.end(),field.name())==m_fields_names.end(),"Error!! TimeInterpolation:add_field, field + " + field.name() + " has already been added");
  auto f0 = field.clone();
  auto f1 = field.clone();
  m_fm_t0->add_field(f0);
  m_fm_t1->add_field(f1);
  m_fields_t0.push_back(f0);
  m_fields_t1.push_back(f1);
  m_fields_names.push_back(field.name());
}

void TimeInterpolation::init(
  const std::vector<std::string>& files)
{
  // Construct map of timestamps and files.  Cycle through all files and check the
  // retrieve the start date and the 'time' variable to back our timestamps.
  set_list_of_files(files);

  // We have a sorted list of files and timestamps.  We can now populate the first two sets of interpolation
  // data using the first 2 indices.
  m_tstamp_ind = -1;
  advance_index_and_update_data();
  advance_index_and_update_data();
  // After that need a function to check if we need to load new data at all.  Check timestamp of the current max data
  // against the simulation timestamp and shift_data, load_new_data if needed.

}

void TimeInterpolation::advance_index_and_update_data()
{
  // Since we are loading new data, we need to shift the data stored for fm_t1 to fm_t0.
  shift_data();

  // Now we update the index to the next timesnap and load new data into t1.
  m_tstamp_ind++;
  const auto data_triplet = m_list_of_timestamps[m_tstamp_ind];
  if (data_triplet.m_filename != m_data_input.get_filename()) {
    // Then we need to close this input stream and open a new file.
    m_data_input.finalize();
    ekat::ParameterList input_params;
    input_params.set("Field Names",m_fields_names);
    input_params.set("Filename",data_triplet.m_filename);
    m_data_input = AtmosphereInput(input_params,m_fm_t1);
  }
  m_data_input.read_variables(data_triplet.m_snap);
}

void TimeInterpolation::shift_data()
{
  for (auto ff = m_fm_t1->begin(); ff != m_fm_t1->end(); ff++) {
    auto field1 = ff->second;
    auto name = field1->name();
    auto field0 = m_fm_t0->get_field(name);
    field0.deep_copy(*field1);
  }
}

void TimeInterpolation::set_list_of_files(const std::vector<std::string>& files)
{
  int it = -1;
  std::map<Real,int>  map_of_times; // A vector listing the timestamps associated with each file
  std::vector<std::string> list_of_files;
  std::vector<TimeStamp>   list_of_timestamps;
  std::vector<int>         list_of_snaps;
  for (int ii=0; ii<files.size(); ii++) {
    const auto filename = files[ii];
    // TODO:  The units of 'time' could also be parsed to get the reference date/time.
    const int date_start = scorpio::get_attribute<int>(filename,"start_date"); // Start date is in YYYYMMDD format
    const int time_start = scorpio::get_attribute<int>(filename,"start_time");         // Start time is in hhmmss format
    // Need to parse the start time and date into a timestamp
    const int YY = date_start/10000;
    const int MM = (date_start - YY*10000)/100;
    const int DD = (date_start - YY*10000 - MM*100);
    const int hh = time_start/10000;
    const int mm = (time_start - hh*10000)/100;
    const int ss = (time_start - hh*10000 - mm*100);
    TimeStamp ts_start(YY,MM,DD,hh,mm,ss);
    // Now we retrieve the 'time' variable and use it to update the map of timestamp to file
    // TODO: This might be good as a scream IO util function since getting time might happen in more than one place.
    scorpio::register_file(filename,scorpio::Read);
    const int ntime = scorpio::get_dimlen(filename,"time");
    for (int tt=0; tt<ntime; tt++) {
      it++;
      auto time_snap = scorpio::read_time_at_index_c2f(filename.c_str(),tt+1);
      TimeStamp ts_snap = ts_start;
      if (time_snap>0) {
        ts_snap += (time_snap*86400); // note, time is assumed to be in days.
      }
      map_of_times.emplace(time_snap,it);
      list_of_files.push_back(filename);
      list_of_snaps.push_back(tt);
      list_of_timestamps.push_back(ts_snap);
    } 
  }
  // Now that we have a sorted map of times from the file we populate an order list of files, timestamps and the
  // index for each time in the associated file.
  for (auto a = map_of_times.begin();a!=map_of_times.end();a++) {
    auto ind = a->second;
    TimesnapTriplet trip;
    trip.m_ts       = list_of_timestamps[ind];
    trip.m_filename = list_of_files[ind];
    trip.m_snap     = list_of_snaps[ind];
    m_list_of_timestamps.push_back(trip);
  }
}
//// Time interpolation routine
//// templated on P as a pack or scalar
//template<typename P>
//void perform_time_interpolation(
//  const std::string& field_name,
//  const TimeStamp&   time,
//        view_Nd<P>&  view_out
//)
//{
//  // Do nothing
//}


} // namespace util
} // namespace scream
