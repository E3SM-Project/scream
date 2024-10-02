#include "scream_output_stream.hpp"

namespace scream
{

std::shared_ptr<OutputStreamBase>
create_output_stream (const ekat::ParameterList& params,
                      const ekat::Comm& comm)
{
  const auto& avg_type_str = params.get<std::string>("Averaging Type");

  if (avg_type_str

}

} // namespace scream
