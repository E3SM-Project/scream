#include <share/util/scream_interpolators.hpp>

#include <pio.hpp>

SourceDataFile::SourceDataFile(const ekat::Comm& comm,
                               const std::string& filename) {
  int err, num_io_tasks = 1, stride = 1, base = 1, rearr = 0;
  err = PIOc_Init_Intracom(comm.mpi_comm(), num_io_tasks, stride, base,
                           rearr, pio_id_);
  err = PIOc_openfile2(pio_id_, &file_id_, &io_type, filename, PIO_NOWRITE);
}

SourceDataFile::~SourceDataFile() {
  int err = PIOc_closefile(nc_id_);
  err = PIOc_finalize(pio_id_);
}

void SourceDataFile::get_times(std::vector<Real>& times) const {
  int err, time_id;
  err = PIOc_inq_dimid(file_id_, "time", &time_id);
  PIO_Offset n;
  err = PIOc_inq_dimlen(file_id_, time_id, &n);
  times.resize(n);
  int time_var_id;
  err = PIOc_inq_varid(file_id_, "time", &time_var_id);
  err = PIOc_read_darray(file_id_, time_var_id, pio_id_, n, &times[0]);
}

void SourceDataFile::get_array(const std::string& name,
                               std::vector<Real>& values) const {
  int var_id;
  err = PIOc_inq_varid(file_id_, name.c_str(), &var_id);
  int var_dim_id;
  err = PIOc_inq_vardimid(file_id_, var_id, &time_var_id);
  PIO_Offset n;
  err = PIOc_inq_dimlen(file_id_, var_id, &n);
  values.resize(n);
  err = PIOc_read_darray(file_id_, var_id, pio_id_, n, &values[0]);
}

