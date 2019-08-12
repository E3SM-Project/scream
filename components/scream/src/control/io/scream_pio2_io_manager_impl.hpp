#include "control/io/scream_pio2_io_manager.hpp"

namespace scream {

// =================== PIO2_Data ================ //

template<typename ScalarType, typename DeviceType>
struct PIO2_IOManager<ScalarType,DeviceType>::PIO2_Data {

  PIO2_Data (const Comm& comm, ParameterList& params) {
    // Num of I/O tasks
    const int comm_size = comm.size();
    this->num_io_tasks = params.get<int>("Number of I/O tasks");
    scream_require_msg(this->num_io_tasks<=comm_size,
                       "Error! Requested more I/O tasks (" + std::to_string(this->num_io_tasks) + ")"
                      " than there are processes in the communicator (" + comm_size + ").\n");
    this->base = 0; // We always assume first I/O task is proc 0. TODO: should we relax this?

    // I/O tasks stride
    this->stride = params.get<int>("I/O tasks stride", comm_size / this->num_io_tasks);
    scream_require_msg(this->base+(this->num_io_tasks-1)*stride < comm_size,
                       "Error! Requested stride (" + std::to_string(this->stride) + ") is too large.\n");

    // I/O tasks rearrangement
    const util::CaseInsensitiveString rearr_type = params.get<std::string>("Rearranger","Subset");
    scream_require_msg(rearr_type=="Subset" || rearr_type=="Box",
                       "Error! Invalid choice (" + rearr_type + ") for parameter 'Rearranger'.\n");
    this->rearr = (rearr_type=="Box" ? PIO_REARR_BOX : PIO_REARR_SUBSET);

    // Setup input/output streams
    this->num_inputs  = params.get<int>("Number of Input Files",0);
    this->num_outputs = params.get<int>("Number of Output Files",0);
    for (int ifile=0; ifile<this->num_inputs; ++ifile) {
      auto& plist = params.sublist(util::strint("Input File",ifile));
      const std::string& fname = plist.get<std::string>("File Name");
    }
    for (int ofile=0; ofile<this->num_outputs; ++ofile) {
      auto& plist = params.sublist(util::strint("Input File",ofile));
      const std::string& fname = plist.get<std::string>("File Name");
    }
  }

  int num_io_tasks;
  int stride;
  int base;
  int rearr;
  int iosysidp;

  int num_inputs;
  int num_outputs;
};

// =================== PIO2_IOManager ================ //

template<typename ScalarType, typename DeviceType>
PIO2_IOManager<ScalarType,DeviceType>::
PIO2_IOManager (const Comm& comm, const ParameterList& params)
 : m_io_params(params)
 , m_atm_comm (comm)
{
  m_pio2_struct.reset(new PIO2_Data(m_atm_comm,m_io_params));

  PIOc_Init_Intracomm(m_atm_comm.mpi_comm(),
                      m_pio2_struct->num_io_tasks,
                      m_pio2_struct->stride,
                      m_pio2_struct->base,
                      m_pio2_struct->rearr,
                      &m_pio2_struct->iosysidp);
}

template<typename ScalarType, typename DeviceType>
util::TimeStamp PIO2_IOManager<ScalarType,DeviceType>::
get_next_restart_step (const util::TimeStamp& current_time) const {
}

template<typename ScalarType, typename DeviceType>
util::TimeStamp PIO2_IOManager<ScalarType,DeviceType>::
get_next_output_step (const util::TimeStamp& current_time) const {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
read_restart_fields (const util::TimeStamp& time) {
}

template<typename ScalarType, typename DeviceType>
bool PIO2_IOManager<ScalarType,DeviceType>::
read_input_fields (const util::TimeStamp& time, const bool must_match) {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
write_restart_fields (const util::TimeStamp& time) {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
write_output_fields  (const util::TimeStamp& time) {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
do_registration_begins () {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
do_registration_ends () {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
register_input_field (const field_type& f) {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
register_output_field (const field_type& f) {
}

template<typename ScalarType, typename DeviceType>
void PIO2_IOManager<ScalarType,DeviceType>::
register_io_field (const field_type& f) {
}

} // namespace scream
