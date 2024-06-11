#ifndef SCREAM_DATA_TYPE_HPP
#define SCREAM_DATA_TYPE_HPP

#include "share/util/scream_utils.hpp"

#include <ekat/ekat_assert.hpp>

#include <string>

namespace scream
{

// An enum for specifying fields data type
enum class DataType {
  Invalid,
  IntType,
  FloatType,
  DoubleType,
  EnsembleDouble,
  EnsembleFloat,
#ifdef SCREAM_DOUBLE_PRECISION
  RealType     = DoubleType,
  EnsembleReal = EnsembleDouble
#else
  RealType     = FloatType,
  EnsembleReal = EnsembleFloat
#endif
};

template<typename ST>
DataType get_data_type () {
  // if statements are compiled out
  if (std::is_same<ST,int>::value) {
    return DataType::IntType;
  } else if (std::is_same<ST,float>::value) {
    return DataType::FloatType;
  } else if (std::is_same<ST,double>::value) {
    return DataType::DoubleType;
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported data type.\n"
        " - typeid(ST): " + std::string(typeid(ST).name()) + "\n");
  }
  return DataType::Invalid;
}

inline bool is_ensemble (const DataType data_type) {
  switch (data_type) {
    case DataType::IntType:         [[fallthrough]];
    case DataType::FloatType:       [[fallthrough]];
    case DataType::DoubleType:      [[fallthrough]];
    case DataType::Invalid:         return false;
    case DataType::EnsembleDouble:  [[fallthrough]];
    case DataType::EnsembleFloat:   return true;
    default:
      EKAT_ERROR_MSG("Error! Unsupported DataType value.\n");
  }
}


inline bool is_narrowing_conversion (const DataType from, const DataType to) {
  return (from==DataType::FloatType || from==DataType::DoubleType) && to==DataType::IntType;
}

inline std::string e2str (const DataType data_type) {
  switch (data_type) {
    case DataType::IntType:         return "int";
    case DataType::FloatType:       return "float";
    case DataType::DoubleType:      return "double";
    case DataType::Invalid:         return "invalid";
    case DataType::EnsembleDouble:  return "ensemble(double)";
    case DataType::EnsembleFloat:   return "ensemble(float)";
    default:
      EKAT_ERROR_MSG("Error! Unsupported DataType value.\n");
  }
}

inline int get_type_size (const DataType data_type) {
  switch (data_type) {
    case DataType::IntType:         return sizeof(int);
    case DataType::FloatType:       return sizeof(float);
    case DataType::DoubleType:      return sizeof(double);
    case DataType::EnsembleDouble:  return sizeof(double)*SCREAM_ENSEMBLE_SIZE;
    case DataType::EnsembleFloat:   return sizeof(float)*SCREAM_ENSEMBLE_SIZE;
    default:
      EKAT_ERROR_MSG("Error! Unsupported DataType value.\n");
  }
}

} // namespace scream

#endif // SCREAM_DATA_TYPE_HPP
