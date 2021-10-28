#include <catch2/catch.hpp>

#include "share/util/scream_utils.hpp"
#include "share/util/scream_time_stamp.hpp"

TEST_CASE("field_layout") {
  using namespace scream;

  std::string A = "A";
  std::string B = "B";
  std::string C = "C";
  std::string D = "D";
  std::string E = "E";
  std::string F = "F";
  std::string G = "G";

  using LOLS_type = std::list<std::list<std::string>>;

  // These three lists do not allow a superset from which they can all be
  // contiguously subviewed.
  LOLS_type lol1 = { {A,B}, {B,C}, {A,C} };
  REQUIRE(contiguous_superset(lol1).size()==0);

  // Input inner lists are not sorted
  REQUIRE_THROWS(contiguous_superset(LOLS_type{ {B,A} }));

  // The following should both allow the superset (A,B,C,D,E,F,G)
  // Note: lol3 is simply a shuffled version of lol2
  LOLS_type lol2 = { {A,B,C}, {B,C,D,E}, {C,D}, {C,D,E,F}, {D,E,F,G} };
  LOLS_type lol3 = { {D,E,F,G}, {C,D,E,F}, {A,B,C}, {C,D}, {B,C,D,E} };

  // Flipping a list is still a valid solution, so consider both tgt and its reverse.
  std::list<std::string> tgt = {A,B,C,D,E,F,G};
  std::list<std::string> tgt_rev = tgt;
  tgt_rev.reverse();

  auto superset2 = contiguous_superset(lol2);
  auto superset3 = contiguous_superset(lol3);
  REQUIRE ( (superset2==tgt || superset2==tgt_rev) );
  REQUIRE ( (superset3==tgt || superset3==tgt_rev) );
}

TEST_CASE ("time_stamp") {
  using namespace scream;
  using TS = util::TimeStamp;

  TS ts1 (2021,10,12,17,8,30);
  REQUIRE (ts1.get_year()==2021);
  REQUIRE (ts1.get_month()==10);
  REQUIRE (ts1.get_day()==12);
  REQUIRE (ts1.get_hours()==17);
  REQUIRE (ts1.get_minutes()==8);
  REQUIRE (ts1.get_seconds()==30);

  // Julian day = frac_of_year_in_days.fraction_of_day, with frac_of_year_in_days=0 at Jan 1st.
  REQUIRE (ts1.frac_of_year_in_days()==(284 + (17*3600+8*60+30)/86400.0));

  REQUIRE (ts1.get_date_string()=="2021-10-12");
  REQUIRE (ts1.get_time_string()=="17:08:30");
  REQUIRE (ts1.to_string()=="2021-10-12.170830");

  // Comparisons
  TS ts2 = ts1;
  ts2 += 10;
  REQUIRE (ts1<ts2);
  REQUIRE (ts2<=ts2);
  REQUIRE (ts2==ts2);

  // Cannot rewind time
  REQUIRE_THROWS (ts2+=-10);

  // Update: check carries
  auto ts3 = ts1 + 1;
  REQUIRE (ts3.get_seconds()==(ts1.get_seconds()+1));
  REQUIRE (ts3.get_minutes()==ts1.get_minutes());
  REQUIRE (ts3.get_hours()==ts1.get_hours());
  REQUIRE (ts3.get_day()==ts1.get_day());
  REQUIRE (ts3.get_month()==ts1.get_month());
  REQUIRE (ts3.get_year()==ts1.get_year());

  ts3 += 60;
  REQUIRE (ts3.get_seconds()==(ts1.get_seconds()+1));
  REQUIRE (ts3.get_minutes()==(ts1.get_minutes()+1));
  REQUIRE (ts3.get_hours()==ts1.get_hours());
  REQUIRE (ts3.get_day()==ts1.get_day());
  REQUIRE (ts3.get_month()==ts1.get_month());
  REQUIRE (ts3.get_year()==ts1.get_year());

  ts3 += 3600;
  REQUIRE (ts3.get_seconds()==(ts1.get_seconds()+1));
  REQUIRE (ts3.get_minutes()==(ts1.get_minutes()+1));
  REQUIRE (ts3.get_hours()==ts1.get_hours()+1);
  REQUIRE (ts3.get_day()==ts1.get_day());
  REQUIRE (ts3.get_month()==ts1.get_month());
  REQUIRE (ts3.get_year()==ts1.get_year());

  ts3 += 86400;
  REQUIRE (ts3.get_seconds()==(ts1.get_seconds()+1));
  REQUIRE (ts3.get_minutes()==(ts1.get_minutes()+1));
  REQUIRE (ts3.get_hours()==(ts1.get_hours()+1));
  REQUIRE (ts3.get_day()==(ts1.get_day()+1));
  REQUIRE (ts3.get_month()==ts1.get_month());
  REQUIRE (ts3.get_year()==ts1.get_year());

  ts3 += 86400*20;
  REQUIRE (ts3.get_seconds()==(ts1.get_seconds()+1));
  REQUIRE (ts3.get_minutes()==(ts1.get_minutes()+1));
  REQUIRE (ts3.get_hours()==(ts1.get_hours()+1));
  REQUIRE (ts3.get_day()==(ts1.get_day()+1+20-31)); // Add 20 days, subtract Oct 31 days (carry)
  REQUIRE (ts3.get_month()==(ts1.get_month()+1));
  REQUIRE (ts3.get_year()==ts1.get_year());

  ts3 += 86400*365;
  REQUIRE (ts3.get_seconds()==ts1.get_seconds()+1);
  REQUIRE (ts3.get_minutes()==(ts1.get_minutes()+1));
  REQUIRE (ts3.get_hours()==(ts1.get_hours()+1));
  REQUIRE (ts3.get_day()==(ts1.get_day()+1+20-31)); // Add 20 days, subtract Oct 31 days (carry)
  REQUIRE (ts3.get_month()==(ts1.get_month()+1));
  REQUIRE (ts3.get_year()==(ts1.get_year()+1));

  // Check update across leap date
#ifdef EKAT_HAS_LEAP_YEAR
  TS ts4(2024,2,28,0,0,0);
  ts4 += 86400;
  REQUIRE (ts4.get_month()==2);
  REQUIRE (ts4.get_day()==29);
  ts4 += 86400;
  REQUIRE (ts4.get_month()==3);
  REQUIRE (ts4.get_day()==1);
#endif

  // Comparisons
  REQUIRE ( TS({2021,12,31},{23,59,59}) < TS({2022,1,1},{0,0,0}));
  REQUIRE ( TS({2022,1,1},{0,0,0}) <= TS({2022,1,1},{0,0,0}));
  REQUIRE ( (TS({2021,12,31},{23,59,59})+1) == TS({2022,1,1},{0,0,0}));
}
