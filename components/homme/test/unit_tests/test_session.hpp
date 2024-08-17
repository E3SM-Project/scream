#ifndef HOMMEXX_TEST_SESSION_HPP
#define HOMMEXX_TEST_SESSION_HPP

#include <map>
#include <string>

namespace Homme {

struct TestSession {

  static TestSession& instance() {
    static TestSession s;
    return s;
  }

  std::map<std::string,std::string> params;
  std::map<std::string,bool>        flags;
private:
  TestSession() = default;
};

} // namespace Homme

#endif // HOMMEXX_TEST_SESSION_HPP
