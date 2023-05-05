#ifndef SCREAM_TIME_INTERPOLATION_HPP
#define SCREAM_TIME_INTERPOLATION_HPP

#include "share/util/scream_time_stamp.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace scream {
namespace util {

class TimeInterpolation {
public:

  TimeInterpolation();

  void init(
    const std::vector<std::string>& list_of_variables,
    const std::vector<std::string>& list_of_files);

  void set_list_of_files(const std::vector<std::string>& files);

protected:

  struct DataTriplet {
  public:
    TimeStamp   m_ts;
    std::string m_filename;
    int         m_snap;

    void print() {
      printf("DataTriplet: timestamp = %s, file = %s, tstep = %d\n",m_ts.to_string().c_str(),m_filename.c_str(),m_snap);
    }
  };

  std::vector<std::string> m_list_of_vars;       // A vector listing the variables that will be stored for interpolation
  std::vector<DataTriplet> m_list_of_timestamps;

//  template<typename P>
//  void perform_time_interpolation(
//    const std::string& field_name,
//    const TimeStamp&   time,
//          view_Nd<P>&  view_out);

}; // class TimeInterpolation



} // namespace util
} // namespace scream

#endif // SCREAM_TIME_INTERPOLATION_HPP
