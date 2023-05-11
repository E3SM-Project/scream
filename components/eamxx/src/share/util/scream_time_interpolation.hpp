#ifndef SCREAM_TIME_INTERPOLATION_HPP
#define SCREAM_TIME_INTERPOLATION_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/io/scorpio_input.hpp"

namespace scream {
namespace util {

class TimeInterpolation {
public:
  using grid_ptr_type    = std::shared_ptr<const AbstractGrid>;

  TimeInterpolation(const grid_ptr_type& grid);

  void init(
    const TimeStamp& timestamp,
    const std::vector<std::string>& list_of_files);

  std::map<std::string,Field> perform_time_interpolation( const TimeStamp& time);

  void set_list_of_files(const std::vector<std::string>& files);
  void add_field(const Field& field);

protected:

  void read_data(const int idx);
  void shift_data();
  void advance_index_and_update_data(const TimeStamp& timestamp);

  // A simple structure to store all the time snap information to be used for interpolation
  struct TimesnapTriplet {
  public:
    TimeStamp   m_ts;
    std::string m_filename;
    int         m_snap;

    void print() {
      printf("TimesnapTriplet: timestamp = %s, file = %s, tstep = %d\n",m_ts.to_string().c_str(),m_filename.c_str(),m_snap);
    }
  };

  grid_ptr_type                 m_grid;
  std::vector<TimesnapTriplet>  m_list_of_timestamps; // A vector that stores the full list of timestamps available for interpolation
  std::shared_ptr<FieldManager> m_fm_t0;
  std::shared_ptr<FieldManager> m_fm_t1;
  TimeStamp                     m_t0;
  TimeStamp                     m_t1;
  std::vector<std::string>      m_fields_names;
  int                           m_tstamp_ind=-1;       // A tally of which timestamps have been loaded from files.
  AtmosphereInput               m_data_input;


}; // class TimeInterpolation



} // namespace util
} // namespace scream

#endif // SCREAM_TIME_INTERPOLATION_HPP
