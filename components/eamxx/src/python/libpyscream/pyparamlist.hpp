#ifndef PYPARAMLIST_HPP
#define PYPARAMLIST_HPP

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/list.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <ekat/ekat_parameter_list.hpp>
#include <functional>

namespace scream {

struct PyParamList {
  ekat::ParameterList pl;
  std::reference_wrapper<ekat::ParameterList> pl_ref;

  PyParamList(ekat::ParameterList &src) : pl_ref(src) {}

  PyParamList(const nanobind::dict &d) : PyParamList(d, "") {}

  PyParamList(const nanobind::dict &d, const std::string &name)
      : pl(name), pl_ref(pl) {
    parse_dict(d, pl);
  }

  PyParamList sublist(const std::string &name) {
    PyParamList spl(pl.sublist(name));
    return spl;
  }

  bool get_bool(const std::string &name) const {
    return pl_ref.get().get<bool>(name);
  }
  int get_int(const std::string &name) const {
    return pl_ref.get().get<int>(name);
  }
  double get_dbl(const std::string &name) const {
    return pl_ref.get().get<double>(name);
  }
  std::string get_str(const std::string &name) const {
    return pl_ref.get().get<std::string>(name);
  }

  std::vector<int> get_int_vec(const std::string &name) const {
    return pl_ref.get().get<std::vector<int>>(name);
  }
  std::vector<double> get_dbl_vec(const std::string &name) const {
    return pl_ref.get().get<std::vector<double>>(name);
  }
  std::vector<std::string> get_str_vec(const std::string &name) const {
    return pl_ref.get().get<std::vector<std::string>>(name);
  }

  template <typename T>
  void set(const std::string &name, T val) {
    pl_ref.get().set(name, val);
  }

  void print() { pl_ref.get().print(); }

 private:
  void parse_dict(const nanobind::dict &d, ekat::ParameterList &p) {
    for(auto item : d) {
      auto key = nanobind::cast<std::string>(nanobind::str(item.first));
      if(nanobind::isinstance<nanobind::str>(item.second)) {
        auto pystr = nanobind::str(item.second);
        p.set<std::string>(key, nanobind::cast<std::string>(pystr));
      } else if(nanobind::isinstance<nanobind::bool_>(item.second)) {
        auto pyint = nanobind::cast<nanobind::bool_>(item.second);
        p.set(key, nanobind::cast<bool>(pyint));
        // p.set(key,pyint.cast<bool>());
      } else if(nanobind::isinstance<nanobind::int_>(item.second)) {
        auto pyint = nanobind::cast<nanobind::int_>(item.second);
        p.set(key, nanobind::cast<int>(pyint));
        // p.set(key,pyint.cast<int>());
      } else if(nanobind::isinstance<nanobind::float_>(item.second)) {
        auto pydouble = nanobind::cast<nanobind::float_>(item.second);
        p.set(key, nanobind::cast<double>(pydouble));
        // p.set(key,pydouble.cast<double>());
      } else if(nanobind::isinstance<nanobind::list>(item.second)) {
        auto pylist = nanobind::cast<nanobind::list>(item.second);
        parse_list(pylist, p, key);
      } else if(nanobind::isinstance<nanobind::dict>(item.second)) {
        auto pydict = nanobind::cast<nanobind::dict>(item.second);
        parse_dict(pydict, p.sublist(key));
      } else {
        EKAT_ERROR_MSG("Unsupported/unrecognized dict entry type.\n");
      }
    }
  }

  void parse_list(const nanobind::list &l, ekat::ParameterList &p,
                  const std::string &key) {
    EKAT_REQUIRE_MSG(
        nanobind::len(l) > 0,
        "Error! Cannot deduce type for dictionary list entry '" + key + "'\n");
    auto first       = l[0];
    bool are_ints    = nanobind::isinstance<nanobind::int_>(first);
    bool are_floats  = nanobind::isinstance<nanobind::float_>(first);
    bool are_strings = nanobind::isinstance<nanobind::str>(first);
    if(are_ints) {
      parse_list_impl<int, nanobind::int_>(l, p, key);
    } else if(are_floats) {
      parse_list_impl<double, nanobind::float_>(l, p, key);
    } else if(are_strings) {
      parse_list_impl<std::string, nanobind::str>(l, p, key);
    } else {
      EKAT_ERROR_MSG("Unrecognized/unsupported list entry type.\n");
    }
  }

  template <typename Txx, typename Tpy>
  void parse_list_impl(const nanobind::list &l, ekat::ParameterList &p,
                       const std::string &key) {
    std::vector<Txx> vals;
    for(auto item : l) {
      EKAT_REQUIRE_MSG(nanobind::isinstance<Tpy>(item),
                       "Error! Inconsistent types in list entries.\n");
      auto item_py = nanobind::cast<Tpy>(item);
      vals.push_back(nanobind::cast<Txx>(item_py));
    }
    p.set(key, vals);
  }
};

inline void nanobind_pyparamlist(nanobind::module_ &m) {
  // Param list
  nanobind::class_<PyParamList>(m, "ParameterList")
      .def(nanobind::init<const nanobind::dict &>())
      .def(nanobind::init<const nanobind::dict &, const std::string &>())
      .def("sublist", &PyParamList::sublist)
      .def("print", &PyParamList::print)
      .def("set", &PyParamList::set<bool>)
      .def("set", &PyParamList::set<int>)
      .def("set", &PyParamList::set<double>)
      .def("set", &PyParamList::set<std::string>)
      .def("set", &PyParamList::set<std::vector<int>>)
      .def("set", &PyParamList::set<std::vector<double>>)
      .def("set", &PyParamList::set<std::vector<std::string>>)
      .def("get_bool", &PyParamList::get_bool)
      .def("get_int", &PyParamList::get_int)
      .def("get_dbl", &PyParamList::get_dbl)
      .def("get_str", &PyParamList::get_str)
      .def("get_int_vec", &PyParamList::get_int_vec)
      .def("get_dbl_vec", &PyParamList::get_dbl_vec)
      .def("get_str_vec", &PyParamList::get_str_vec);
}

}  // namespace scream

#endif  // PYPARAMLIST_HPP
