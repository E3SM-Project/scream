#include "scream_scorpio_interface.hpp"

#include "scream_config.h"

#include <ekat/ekat_assert.hpp>
#include <ekat/util/ekat_string_utils.hpp>

#include <pio.h>

#include <numeric>

namespace scream {
namespace scorpio {

// This class is an implementation detail, and therefore it is hidden inside
// a cpp file. All customers of IO capabilities must use the common interfaces
// exposed in the header file of this source file.
// This class simply serves as a container for persistent IO data.
struct ScorpioSession
{
public:
  static ScorpioSession& instance () {
    static ScorpioSession s;
    return s;
  }

  template<typename T>
  using strmap_t = std::map<std::string,T>;

  strmap_t<PIOFile>                     files;
  strmap_t<std::shared_ptr<PIODecomp>>  decomps;

  // In the above map, we would like to label decomps as dtype_dim1$N1_dim2$N2...
  // where N$i is the global length of dim$i. However, it *may* happen that we use
  // two different decompositions for the same global layout, which would clash
  // the names. For this reason, we append at the end an increasing counter,
  // which disambiguate between globally-equivalent partitions.
  // When adding a new decomp, we check this map to see if another decomp
  // already exists with the same global layout. If so, we first check to see
  // if that decomp is equivalent to the new one *on all ranks*. If yes, we
  // recycle it, otherwise we create a new PIO decomp.
  // strmap_t<std::list<std::string>>    decomp_global_layout_to_decomp_name;
  // strmap_t<int>                       decomp_global_layout_to_counter;

  int         pio_sysid      = -1;
  int         pio_rearranger = -1;
  int         pio_type       = -1;
  int         pio_format     = -1;

  ekat::Comm  comm;

private:

  ScorpioSession () = default;
};

// --------------------------------------------------------------------------------------------- //

// Utility for common IO operation failure
void check_scorpio_noerr (const int err,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: PIOc_" + pioc_func_name + "\n");
}

void check_scorpio_noerr (const int err, const std::string& filename,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - filename: " + filename + "\n" 
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: " + pioc_func_name + "\n");
}

void check_scorpio_noerr (const int err,
                          const std::string& filename,
                          const std::string& entity_type,
                          const std::string& entity_name,
                          const std::string& func_name,
                          const std::string& pioc_func_name)
{
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "Error! Something went wrong while performing a pio operation.\n"
      " - pio error code: " + std::to_string(err) + "\n"
      " - filename: " + filename + "\n" 
      " - " + entity_type + ": " + entity_name + "\n" 
      " - from interface function: scorpio::" + func_name + "\n"
      " - calling PIOc function: " + pioc_func_name + "\n");
}

// Return name of a shared ptr to PIO entity (to use inside ekat::join)
std::string get_entity_name (const std::shared_ptr<const PIOEntity>& e)
{
  return e->name;
}

template<typename T>
std::string print_map_keys (const std::map<std::string,T>& map) {
  std::string s;
  for (const auto& it : map) {
    s += it.first + ",";
  }
  s.pop_back();
  return s;
}

// Retrieve the int codes PIO uses to specify data types
int nctype (const std::string& type) {
  if (type=="int") {
    return PIO_INT;
  } else if (type=="int64") {
    return PIO_INT64;
  } else if (type=="float" || type=="single") {
    return PIO_FLOAT;
  } else if (type=="double") {
    return PIO_DOUBLE;
  } else if (type=="real") {
#if defined(SCREAM_DOUBLE_PRECISION)
    return PIO_DOUBLE;
#else
    return PIO_FLOAT;
#endif
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported data type '" + type + "'.\n");
  }
}

// ====================== Local utilities ========================== //
// Note: these utilities are used in this file to retrieve PIO entities,
//       so that we implement all checks once (rather than in every function)

PIOFile& get_file (const std::string& filename,
                   const std::string& context)
{
  auto& s = ScorpioSession::instance();

  EKAT_REQUIRE_MSG (s.files.count(filename)==1,
      "Error! Could not retrieve the file. File not open.\n"
      " - filename: " + filename + "\n"
      "Context:\n"
      " " + context + "\n");

  return s.files.at(filename);
}

PIODim& get_dim (const std::string& filename,
                 const std::string& dimname,
                 const std::string& context)
{
  const auto& f = get_file(filename,context);
  EKAT_REQUIRE_MSG (f.dims.count(dimname)==1,
      "Error! Could not retrieve dimension. Dimension not found.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n"
      " - dims on file: " + print_map_keys(f.dims) + "\n"
      "Context:\n"
      " " + context + "\n");

  return *f.dims.at(dimname);
}

PIOVar& get_var (const std::string& filename,
                 const std::string& varname,
                 const std::string& context)
{
  const auto& f = get_file(filename,context);
  EKAT_REQUIRE_MSG (f.vars.count(varname)==1,
      "Error! Could not retrieve variable. Variable not found.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n"
      " - vars on file : " + print_map_keys(f.vars) + "\n"
      "Context:\n"
      " " + context + "\n");

  return *f.vars.at(varname);
}

// ====================== Global IO operations ======================= // 

void init_pio_subsystem(const ekat::Comm& comm, const int atm_id)
{
  auto& s = ScorpioSession::instance();
  s.comm = comm;

  EKAT_REQUIRE_MSG (s.pio_sysid==-1,
      "Error! Attmept to re-initialize pio subsystem.\n");

#ifdef SCREAM_CIME_BUILD
  s.pio_sysid      = shr_get_iosysid_c2d(atm_id);
  s.pio_type       = shr_get_iotype_c2f(atm_id);
  s.pio_rearranger = shr_get_rearranger_c2f(atm_id);
  s.pio_format     = shr_get_ioformat_c2f(atm_id);
#else
  // Use some reasonable defaults for standalone EAMxx tests
  int stride = 1;
  int base = 0;

  s.pio_rearranger = PIO_REARR_SUBSET;
  s.pio_type       = PIO_IOTYPE_PNETCDF;
  s.pio_format     = PIO_64BIT_DATA;

  auto err = PIOc_Init_Intracomm(comm.mpi_comm(), comm.size(), stride, base, s.pio_rearranger, &s.pio_sysid);
  check_scorpio_noerr (err,"init_pio_subsystem", "Init_Intracomm");

  // Unused in standalone mode
  (void) atm_id;
#endif

  static_assert (sizeof(offset_t)==sizeof(PIO_Offset),
      "Error! PIO was configured with PIO_OFFSET not a 64-bit int.\n");
}

bool is_pio_subsystem_inited () {
  return ScorpioSession::instance().pio_sysid!=-1;
}

void finalize_pio_subsystem ()
{
  auto& s = ScorpioSession::instance();

  // TODO: should we simply return instead? I think trying to finalize twice
  //       *may* be a sign of possible bugs, though with Catch2 testing
  //       I *think* there may be some issue with how the code is run.
  EKAT_REQUIRE_MSG (s.pio_sysid!=-1,
      "Error! PIO subsystem was already finalized.\n");

  for (auto& it : s.files) {
    EKAT_REQUIRE_MSG (it.second.num_customers==0,
      "Error! ScorpioSession::finalize called, but a file is still in use elsewhere.\n"
      " - filename: " + it.first + "\n");
  }
  s.files.clear();
  for (auto& it : s.decomps) {
    EKAT_REQUIRE_MSG (it.second.use_count()==1,
      "Error! ScorpioSession::finalize called, but a decomp is still stored elsewhere.\n"
      " - decomp name: " + it.first + "\n");
  }
  s.decomps.clear();

  PIOc_finalize (s.pio_sysid);

  s.pio_type       = -1;
  s.pio_sysid      = -1;
  s.pio_format     = -1;
  s.pio_rearranger = -1;
}

// ========================= File operations ===================== //

void register_file (const std::string& filename, const FileMode mode)
{
  auto& s = ScorpioSession::instance();
  auto& f = s.files[filename];
  EKAT_REQUIRE_MSG (f.mode==Unset || f.mode==mode,
      "Error! File was already opened with a different mode.\n"
      " - filename: " + filename + "\n"
      " - old mode: " + e2str(f.mode) + "\n"
      " - new mode: " + e2str(f.mode) + "\n");

  if (f.mode == Unset) {
    // First time we ask for this file. Call PIO open routine(s)
    int err;
    if (mode & Read) {
      auto write = mode & Write ? PIO_WRITE : PIO_NOWRITE;
      err = PIOc_openfile(s.pio_sysid,&f.ncid,&s.pio_type,filename.c_str(),write);
      f.enddef = true;
    } else {
      err = PIOc_createfile(s.pio_sysid,&f.ncid,&s.pio_type,filename.c_str(),s.pio_format);
      f.enddef = false;
    }

    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "Error! Something went wrong while opening a file.\n"
        " - filename : " + filename + "\n"
        " - file mode: " + e2str(mode) + "\n"
        " - pio error: " + std::to_string(err) + "\n");

    f.mode = mode;
    f.name = filename;

    if (mode & Read) {
      // Read all dims/vars from file
      PIO_Offset len;
      int ndims, nvars, ngatts, unlimdimid;
      err = PIOc_inq(f.ncid, &ndims, &nvars, &ngatts, &unlimdimid);
      check_scorpio_noerr(err,f.name,"register_file","inq");

      char name[PIO_MAX_NAME];
      for (int idim=0; idim<ndims; ++idim) {
        err = PIOc_inq_dim(f.ncid,idim,name,&len);
        check_scorpio_noerr (err,f.name,"register_file","inq_dim");
        auto dim = f.dims[name] = std::make_shared<PIODim>();
        dim->name = name;
        dim->length = dim->length = len;
        dim->ncid = idim;
        dim->unlimited = idim==unlimdimid;

        if (dim->unlimited) {
          f.time_dim = dim;
        }
      }

      auto find_dim = [&](int dimid) -> std::shared_ptr<const PIODim> {
        std::shared_ptr<const PIODim> d;
        for (auto it : f.dims) {
          if (it.second->ncid==dimid) {
            d = it.second;
          }
        }
        EKAT_REQUIRE_MSG (d!=nullptr,
            "Error! Could not locat dimension id in the file.\n"
            " - filename: " + f.name + "\n"
            " - dim id  : " + std::to_string(dimid) + "\n");
        return d;
      };
      int dtype,natts;
      int dimids[PIO_MAX_DIMS];
      for (int ivar=0; ivar<nvars; ++ivar) {
        err = PIOc_inq_var(f.ncid,ivar,name,PIO_MAX_NAME,&dtype,&ndims,dimids,&natts);
        check_scorpio_noerr(err,f.name,"register_file","inq_var");

        auto var = f.vars[name] = std::make_shared<PIOVar>();
        var->name = name;
        var->ncid = ivar;
        for (int idim=0; idim<ndims; ++idim) {
          auto dim = find_dim(dimids[idim]);
          if (dim->unlimited) {
            var->time_dep = true;
            var->num_records = dim->length;
          } else {
            var->dims.push_back(find_dim(dimids[idim]));
          }
        }

        for (const auto& dt : {"int","int64","float","double"}) {
          if (nctype(dt)==dtype) {
            var->dtype = var->nc_dtype = dt;
            break;
          }
        }
        EKAT_REQUIRE_MSG (var->dtype!="",
            "Error! Variable data type not supported.\n"
            " - filename: " + filename + "\n"
            " - varname : " + var->name + "\n"
            " - nc dtype: " + std::to_string(dtype) + "\n");

        for (int iatt=0; iatt<natts; ++iatt) {
          err = PIOc_inq_attname(f.ncid,ivar,iatt,name);
          check_scorpio_noerr(err,f.name,"register_file","inq_attname");
          if (std::string(name)=="units") {
            err = PIOc_get_att_text(f.ncid,ivar,"units",name);
            var->units = name;
            break;
          }
        }
      }
    }
  }
  ++f.num_customers;
}

void release_file  (const std::string& filename)
{
  auto& f = get_file(filename,"scorpio::release_file");

  --f.num_customers;
  if (f.num_customers>0) {
    return;
  }

  int err;
  if (f.mode & Write) {
    err = PIOc_sync(f.ncid);
    check_scorpio_noerr (err,f.name,"release_file","sync");
  }

  err = PIOc_closefile(f.ncid);
  check_scorpio_noerr (err,f.name,"release_file","closefile");

  auto& s = ScorpioSession::instance();
  s.files.erase(filename);
}

void sync_file (const std::string &filename)
{
  auto& f = get_file(filename,"scorpio::sync_file");
  
  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Cannot call sync_file. File is read-only.\n"
      " - filename: " + filename + "\n");

  int err = PIOc_sync(f.ncid);
  check_scorpio_noerr (err,f.name,"sync_file","sync");
}

void redef(const std::string &filename)
{
  auto& f = get_file(filename,"scorpio::redef");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not call redef on the input file. File is read-only.\n"
      " - filename: " + filename + "\n");

  if (f.enddef) {
    int err = PIOc_redef(f.ncid);
    check_scorpio_noerr (err,f.name,"redef","redef");
    f.enddef = false;
  }
}

void enddef(const std::string &filename)
{
  auto& f = get_file(filename,"scorpio::enddef");

  if (not f.enddef) {
    int err = PIOc_enddef(f.ncid);
    check_scorpio_noerr (err,f.name,"enddef","enddef");
    f.enddef = true;
  }
}

bool is_file_open (const std::string& filename, const FileMode mode)
{
  auto& s = ScorpioSession::instance();
  auto it = s.files.find(filename);
  if (it==s.files.end()) return false;

  return mode==Unset || (mode & it->second.mode);
}

// =================== Dimensions operations ======================= //

void define_dim (const std::string& filename, const std::string& dimname, const int length)
{
  auto& f = get_file(filename,"scorpio::define_dim");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not define dimension. File is read-only.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  auto& dim = f.dims[dimname];

  bool unlimited = length==0;

  if (dim==nullptr) {
    // Create new dimension
    dim = std::make_shared<PIODim>();
    dim->name = dimname;
    dim->length = length;
    dim->unlimited = unlimited;

    // Define the dimension in PIO
    int err = PIOc_def_dim(f.ncid,dimname.c_str(),dim->length,&dim->ncid);
    check_scorpio_noerr (err,f.name,"dimension",dimname,"define_dim","def_dim");
  } else {
    // Already defined. Check that the dim specs are the same.
    EKAT_REQUIRE_MSG (unlimited==dim->unlimited,
        "Error! Redefining dimension with different unlimited flag.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - old unlimited:" + (dim->unlimited ? "yes" : "no" )+ "\n"
        " - new unlimited:" + (unlimited ? "yes" : "no") + "\n");

    EKAT_REQUIRE_MSG (unlimited || length==dim->length,
        "Error! Redefining dimension with a different (local) length.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - old length:" + std::to_string(dim->length)+ "\n"
        " - new length:" + std::to_string(length) + "\n");
  }
}

bool has_dimension (const std::string& filename, const std::string& dimname, const int length)
{
  const auto& f = get_file(filename,"scorpio::has_dimension");

  auto it = f.dims.find(dimname);
  if (it==f.dims.end()) {
    return false;
  } else if (length==0) {
    return it->second->unlimited;
  } else if (length>0) {
    return it->second->length==length;
  }
  return true;
}

int get_dimlen (const std::string& filename, const std::string& dimname)
{
  EKAT_REQUIRE_MSG (has_dimension(filename,dimname),
      "Error! Could not inquire dimension length. The dimension is not in the file.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");

  auto& s = ScorpioSession::instance();
  auto& f = s.files[filename];

  return f.dims.at(dimname)->length;
}

// =================== Decompositions operations ==================== //

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const std::vector<offset_t>& my_offsets,
                     const std::string& decomp_name)
{
  auto& d = get_dim(filename,dimname,"scorpio::set_dim_decomp");

  EKAT_REQUIRE_MSG (not d.unlimited,
      "Error! Cannot partition an unlimited dimension.\n"
      " - filename: " + filename + "\n"
      " - dimname : " + dimname + "\n");
  
  if (d.decomps.count(decomp_name)==1) {
    // Should we allow to redefine the same decomposition for the same dim?
    // I'm not sure, but for now I'll allow it.

    // Was already partitioned with this tag. Check that we're decomposing in the same way
    int same = d.decomps[decomp_name]==my_offsets;
    const auto& comm = ScorpioSession::instance().comm;
    comm.all_reduce(&same,1,MPI_MIN);
    EKAT_REQUIRE_MSG(same==1,
        "Error! Attempt to redefint a decomposition with a different dofs distribution.\n"
        " - filename: " + filename + "\n"
        " - dimname : " + dimname + "\n"
        " - decomp name: " + decomp_name + "\n");
  } else {
    // Check that offsets are less than the global dimension length
    for (auto o : my_offsets) {
      EKAT_REQUIRE_MSG (o>=0 && o<d.length,
          "Error! Offset for dimension decomposition is out of bounds.\n"
          " - filename: " + filename + "\n"
          " - dimname : " + dimname + "\n"
          " - dim glen: " + std::to_string(d.length) + "\n"
          " - offset  : " + std::to_string(o) + "\n");
    }

    d.decomps[decomp_name] = my_offsets;
  }
}

void set_dim_decomp (const std::string& filename,
                     const std::string& dimname,
                     const offset_t start, const offset_t count,
                     const std::string& decomp_name)
{
  std::vector<offset_t> offsets(count);
  std::iota(offsets.begin(),offsets.end(),start);
  set_dim_decomp(filename,dimname,offsets,decomp_name);
}

void set_var_decomp (const std::string& filename,
                     const std::string& varname,
                     const std::string& dimname,
                     const std::string& dim_decomp_name,
                     const bool throw_if_decomp_already_set,
                     const bool throw_if_var_does_not_have_decomp_dim)
{
  auto& f = get_file(filename,"scorpio::set_var_decomp");

  auto& var = get_var(filename,varname,"scorpio::set_var_decomp");

  EKAT_REQUIRE_MSG (not throw_if_decomp_already_set or var.decomp==nullptr,
      "Error! Variable decomposition was already set, and you requested to not allow a reset.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n"
      " - old decomp name: " + var.decomp->name + "\n"
      " - new dim decomp name: " + dim_decomp_name + "\n");

  // Create decomp name: dtype-dim1$N1_dim2$N2_..._dimk$Nk-$dim_decomp_name
  std::shared_ptr<const PIODim> decomp_dim;
  std::string decomp_tag = var.dtype + "-";
  int idecompdim = -1;
  int ndims = var.dims.size();
  for (int idim=0; idim<ndims; ++idim) {
    auto d = var.dims[idim];
    if (d->name==dimname) {
      decomp_dim = d;
      idecompdim = idim;
    }
    decomp_tag += d->name + std::to_string(d->length) + "_";
  }

  if (decomp_tag.back()=='_') {
    decomp_tag.pop_back();
  }

  if (decomp_dim==nullptr and not throw_if_var_does_not_have_decomp_dim) {
    // We are allowing to call this fucntion on not-decomposed vars.
    return;
  }

  EKAT_REQUIRE_MSG (decomp_dim!=nullptr,
      "Error! Could not locate the decomp dim among the variable dimensions.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n"
      " - decomp dim: " + dimname + "\n"
      " - var dims  : " + ekat::join(var.dims,get_entity_name,",") + "\n");

  if (dim_decomp_name=="DEFAULT") {
    EKAT_REQUIRE_MSG (decomp_dim->decomps.size()==1,
        "Error! 'DEFAULT' dim_decomp_name only allowed if there is one decomposition for the dimension.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n"
        " - dimname : " + dimname + "\n"
        " - dim decmops: " + print_map_keys(decomp_dim->decomps) + "\n");
  }
  
  // The decomp_tag is dtype_dim1$N1_dim2$N2_.._dimK$NK-$decomp_name, where decomp_name is
  // the tag passed to decomp_dimension.
  decomp_tag += "-" + dim_decomp_name;

  // Check if a decomp with this name already exists
  auto& s = ScorpioSession::instance();
  auto& decomp = s.decomps[decomp_tag];
#ifndef NDEBUG
  // Extra check: all ranks must agree on whether they have the decomposition!
  int my_found = decomp==nullptr ? 0 : 1;
  int found;
  const auto& comm = ScorpioSession::instance().comm;
  comm.all_reduce(&my_found,&found,1,MPI_MIN);
  EKAT_REQUIRE_MSG(found==my_found,
      "Error! Decomposition already present on some ranks but not all.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n"
      " - var dims: " + ekat::join(var.dims,get_entity_name,",") + "\n"
      " - decopm tag: " + decomp_tag + "\n");
#endif

  if (decomp==nullptr) {
    // We haven't create this decomp yet. Go ahead and create one
    decomp = std::make_shared<PIODecomp>();
    decomp->name = decomp_tag;

    EKAT_REQUIRE_MSG (idecompdim==0,
        "Error! We currently only allow decomposition on slowest-striding dimension.\n"
        "       Generalizing is not complicated, but it was not a priority.\n"
        " - filename  : " + filename + "\n"
        " - varname   : " + varname + "\n"
        " - var dims  : " + ekat::join(var.dims,get_entity_name,",") + "\n"
        " - decomp dim: " + decomp_dim->name + "\n");

    std::vector<int> gdimlen;
    int non_decomp_dim_prod = 1;

    for (int idim=0; idim<ndims; ++idim) {
      auto d = var.dims[idim];
      gdimlen.push_back(d->length);
      if (idim!=idecompdim) {
        non_decomp_dim_prod *= d->length;
      }
    }

    const auto& decomp_dim_offsets = decomp_dim->decomps.at(dim_decomp_name);

    int decom_dim_loc_len = decomp_dim_offsets.size();
    decomp->offsets.resize (non_decomp_dim_prod*decom_dim_loc_len);
    for (int idof=0; idof<decom_dim_loc_len; ++idof) {
      auto dof_offset = decomp_dim_offsets[idof];
      auto beg = decomp->offsets.begin()+ idof*non_decomp_dim_prod;
      auto end = beg + non_decomp_dim_prod;
      std::iota (beg,end,non_decomp_dim_prod*dof_offset);
    }

    // Create PIO decomp
    int maplen = decomp->offsets.size();
    PIO_Offset* compmap = reinterpret_cast<PIO_Offset*>(decomp->offsets.data());
    int err = PIOc_init_decomp(s.pio_sysid,nctype(var.dtype),ndims,gdimlen.data(),
                               maplen,compmap, &decomp->ncid,s.pio_rearranger,
                               nullptr,nullptr);
    check_scorpio_noerr(err,f.name,"decomp",decomp_tag,"set_var_decomp","InitDecomp");
  }
  var.decomp = decomp;
  var.decomp_dim = dimname;
  var.dim_decomp_name = dim_decomp_name;
}

void set_vars_decomp (const std::string& filename,
                      const std::vector<std::string>& varnames,
                      const std::string& dimname,
                      const std::string& dim_decomp_name,
                      const bool throw_if_decomp_already_set,
                      const bool throw_if_var_does_not_have_decomp_dim)
{
  for (const auto& vname : varnames) {
    set_var_decomp(filename,vname,dimname,dim_decomp_name,
                   throw_if_decomp_already_set,
                   throw_if_var_does_not_have_decomp_dim);
  }
}

// Clean up any currently unused decompositions (meaning the decomp is associated to *no* variables)
void free_unused_decomps ()
{
  auto& s = ScorpioSession::instance();
  for (auto it = s.decomps.begin(); it!=s.decomps.end(); )
  {
    if (it->second.use_count()==1) {
      it = s.decomps.erase(it);
    } else {
      ++it;
    }
  }
}

// ================== Variable operations ================== //

// Define var on output file (cannot call on Read/Append files)
void define_var (const std::string& filename, const std::string& varname,
                 const std::string& units, const std::vector<std::string>& dimensions,
                 const std::string& dtype, const std::string& nc_dtype,
                 const bool time_dep)
{
  auto& f = get_file(filename,"scorpio::define_var");

  EKAT_REQUIRE_MSG (f.mode & Write,
      "Error! Could not define variable. File is read-only.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  EKAT_REQUIRE_MSG (not time_dep || f.time_dim!=nullptr,
      "Error! Cannot define time-dependent variable: no time dimension defined.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  if (f.vars.count(varname)==0) {
    // Create new variable
    auto var = std::make_shared<PIOVar>();
    var->name = varname;
    var->units = units;
    var->dtype = dtype;
    var->nc_dtype = nc_dtype;
    var->time_dep = time_dep;
    int ndims = dimensions.size() + (time_dep ? 1 : 0);
    std::vector<int> dimids;
    if (time_dep) {
      dimids.push_back(f.time_dim->ncid);
    }
    for (const auto& dname : dimensions) {
      EKAT_REQUIRE_MSG (has_dimension(filename,dname),
          "Error! Cannot create variable. Dimension not found.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - var dims : " + ekat::join(dimensions,",") + "\n"
          " - file dims: " + print_map_keys(f.dims) + "\n");
      auto dim = f.dims.at(dname);
      var->dims.push_back(dim);
      dimids.push_back(dim->ncid);
    }

    // Define the variable in PIO
    int err = PIOc_def_var(f.ncid,varname.c_str(),nctype(nc_dtype),ndims,dimids.data(),&var->ncid);
    check_scorpio_noerr(err,f.name,"variable",varname,"define_var","def_var");

    f.vars[varname] = var;

    if (units!="") {
      // Add units attribute
      set_attribute(filename,varname,"units",units);
    }
  } else {
    const auto& var = f.vars.at(varname);
    // The variable was already defined. Check that important metadata is the same
    EKAT_REQUIRE_MSG (var->units==units,
        "Error! Attempt to redefine variable with different units.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old units: " + var->units + "\n"
          " - new units: " +      units + "\n");
    EKAT_REQUIRE_MSG (var->dtype==dtype,
        "Error! Attempt to redefine variable with different data type.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old dtype: " + var->dtype + "\n"
          " - new dtype: " +      dtype + "\n");
    EKAT_REQUIRE_MSG (var->nc_dtype==nc_dtype,
        "Error! Attempt to redefine variable with different PIO data type.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old pio dtype: " + var->nc_dtype + "\n"
          " - new pio dtype: " +      nc_dtype + "\n");
    EKAT_REQUIRE_MSG (var->time_dep==time_dep,
        "Error! Attempt to redefine variable with different time dep flag.\n"
          " - filename : " + filename + "\n"
          " - varname  : " + varname + "\n"
          " - old time_dep: " + (var->time_dep ? "yes" : "no") + "\n"
          " - new time_dep: " + (     time_dep ? "yes" : "no") + "\n");
    const auto var_dims = ekat::join(var->dims,get_entity_name,",");
    EKAT_REQUIRE_MSG (var_dims==ekat::join(dimensions,","),
        "Error! Attempt to redefine variable with different dimensions.\n"
          " - filename: " + filename + "\n"
          " - varname : " + varname + "\n"
          " - old dims: " + var_dims + "\n"
          " - new dims: " + ekat::join(dimensions,",") + "\n");
  }
}

void define_var (const std::string& filename, const std::string& varname,
                 const std::vector<std::string>& dimensions,
                 const std::string& dtype,
                 const bool time_dependent)
{
  define_var(filename,varname,"",dimensions,dtype,dtype,time_dependent);
}

void change_var_dtype (const std::string& filename, const std::string& varname,
                       const std::string& dtype)
{
  auto& var = get_var(filename,varname,"scorpio::change_var_dtype");

  if (dtype==var.dtype) {
    return;
  }

  var.dtype = dtype;
  if (var.decomp) {
    // Re-decompose the variable, with new data type
    set_var_decomp (filename,varname,var.decomp_dim,var.dim_decomp_name,false,true);
  }
}

bool has_variable (const std::string& filename, const std::string& varname,
                   const std::vector<std::string>& dims,
                   const std::string& units,
                   const std::string& time_dep)
{
  const auto& f = get_file(filename,"scorpio::has_variable");

  if (f.vars.count(varname)==0) {
    return false;
  }

  const auto& var = *f.vars.at(varname);
  if (units!="" && var.units!=units) {
    return false;
  }

  if (dims.size()!=0) {
    size_t ndims     = dims.size();
    if (var.dims.size()!=ndims) {
      return false;
    }

    for (size_t idx=0; idx!=ndims; ++idx) {
      if (dims[idx]!=var.dims[idx]->name) {
        return false;
      }
    }
  }

  if (ekat::upper_case(time_dep)=="YES" or
      ekat::upper_case(time_dep)=="TRUE") {
    return var.time_dep;
  } else if (ekat::upper_case(time_dep)=="NO" or
             ekat::upper_case(time_dep)=="FALSE") {
    return not var.time_dep;
  } else {
    return true;
  }
}

std::vector<std::string> get_vardims (const std::string& filename,
                                      const std::string& varname)
{
  const auto& v = get_var(filename,varname,"scorpio::get_vardims");
  std::vector<std::string> dims;
  for (auto d : v.dims) {
    dims.push_back(d->name);
  }
  return dims;
}

void define_time (const std::string& filename, const std::string& units, const std::string& time_name)
{
  auto& f = get_file(filename,"scorpio::define_time");
  EKAT_REQUIRE_MSG (f.time_dim==nullptr,
      "Error! Attempt to redeclare unlimited dimension.\n"
      " - filename: " + filename + "\n");

  define_dim(filename,time_name,0);
  f.time_dim = f.dims.at(time_name);

  define_var(filename,time_name,units,{},"double","double",true);
}

// Update value of time variable, increasing time dim length
void update_time(const std::string &filename, const double time) {
  const auto& f = get_file(filename,"scorpio::update_time");
        auto& time_dim = *f.time_dim;
  const auto& var = get_var(filename,time_dim.name,"scorpio::update_time");

  PIO_Offset index = time_dim.length;
  int err = PIOc_put_var1(f.ncid,var.ncid,&index,&time);
  check_scorpio_noerr (err,f.name,"update time","put_var1");
  ++time_dim.length;
}

double get_time (const std::string& filename, const int time_index)
{
  const auto& f = get_file(filename,"scorpio::get_time");
  const auto& time_name = f.time_dim->name;
  double t;
  read_var(filename,time_name,&t,time_index);
  return t;
}

std::vector<double> get_all_times (const std::string& filename)
{
  const auto& f = get_file(filename,"scorpio::get_all_times");
  const auto& dim = *f.time_dim;

  std::vector<double> times (dim.length);
  for (int i=0; i<dim.length; ++i) {
    read_var (filename, dim.name, &times[i], i);
  }
  return times;
}

// Read variable into user provided buffer.
// If time dim is present, read given time slice (time_index=-1 means "read last record).
// If time dim is not present, time_index must be -1 (error out otherwise)
template<typename T>
void read_var (const std::string &filename, const std::string &varname, T* buf, const int time_index)
{
  EKAT_REQUIRE_MSG (buf!=nullptr,
      "Error! Cannot read from provided pointer. Invalid buffer pointer.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  const auto& f = get_file(filename,"scorpio::read_var");
  const auto& var = get_var(filename,varname,"scorpio::read_var");

  int err, frame;
  if (var.time_dep) {
    frame = time_index>=0 ? time_index : f.time_dim->length-1;
    EKAT_REQUIRE_MSG (frame<f.time_dim->length,
        "Error! Time index out of bounds.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n"
        " - time idx: " + std::to_string(time_index) + "\n"
        " - time len: " + std::to_string(f.time_dim->length));
    err = PIOc_setframe(f.ncid,var.ncid,frame);
    check_scorpio_noerr (err,f.name,"variable",varname,"read_var","setframe");
  }

  std::string pioc_func;
  if (var.decomp) {
    // A decomposed variable, requires read_darray
    err = PIOc_read_darray(f.ncid,var.ncid,var.decomp->ncid,var.decomp->offsets.size(),buf);
    pioc_func = "read_darray";
  } else {
    // A non-decomposed variable, use PIOc_get_var(a)
    if (var.time_dep) {
      // We need to get the start/count for each dimension
      int ndims = var.dims.size();
      std::vector<PIO_Offset> start (ndims+1,0), count(ndims+1); // +1 for time
      start[0] = frame;
      count[0] = 1;
      for (int idim=0; idim<ndims; ++idim) {
        count[idim+1] = var.dims[idim]->length;
      }
      err = PIOc_get_vara(f.ncid,var.ncid,start.data(),count.data(),buf);
      pioc_func = "get_vara";
    } else {
      // Easy: just pass the buffer, and read all entries
      err = PIOc_get_var(f.ncid,var.ncid,buf);
      pioc_func = "get_var";
    }
  }
  check_scorpio_noerr (err,f.name,"variable",varname,"read_var",pioc_func);
}

// Write data from user provided buffer into the requested variable
template<typename T>
void write_var (const std::string &filename, const std::string &varname, const T* buf, const T* fillValue)
{
  EKAT_REQUIRE_MSG (buf!=nullptr,
      "Error! Cannot write in provided pointer. Invalid buffer pointer.\n"
      " - filename: " + filename + "\n"
      " - varname : " + varname + "\n");

  const auto& f = get_file(filename,"scorpio::write_var");
  auto& var = get_var(filename,varname,"scorpio::write_var");
  int err;

  if (var.time_dep) {
    ++var.num_records;
    EKAT_REQUIRE_MSG (var.num_records==f.time_dim->length,
        "Error! Number of records for variable does not match time length.\n"
        " - filename: " + filename + "\n"
        " - varname : " + varname + "\n"
        " - time len: " + std::to_string(f.time_dim->length) + "\n"
        " - nrecords: " + std::to_string(var.num_records) + "\n");
    err = PIOc_setframe (f.ncid,var.ncid,var.num_records-1);
    check_scorpio_noerr (err,f.name,"variable",varname,"write_var","setframe");
  }

  std::string pioc_func;
  if (var.decomp) {
    // A decomposed variable, requires write_darray
    err = PIOc_write_darray(f.ncid,var.ncid,var.decomp->ncid,var.decomp->offsets.size(),buf,fillValue);
    pioc_func = "write_darray";
  } else {
    // A non-decomposed variable, use PIOc_put_var(a)
    if (var.time_dep) {
      // We need to get the start/count for each dimension
      int ndims = var.dims.size();
      std::vector<PIO_Offset> start (ndims+1,0), count(ndims+1);
      start[0] = f.time_dim->length-1;
      count[0] = 1;
      for (int idim=0; idim<ndims; ++idim) {
        count[idim+1] = var.dims[idim]->length;
      }
      err = PIOc_put_vara(f.ncid,var.ncid,start.data(),count.data(),buf);
      pioc_func = "put_vara";
    } else {
      // Easy: just pass the buffer, and write all entries
      err = PIOc_put_var(f.ncid,var.ncid,buf);
      pioc_func = "put_var";
    }
  }
  check_scorpio_noerr (err,f.name,"variable",varname,"write_var",pioc_func);
}

// ========================== READ/WRITE ETI ========================== //

template void read_var<int>       (const std::string&, const std::string&, int*,       const int);
template void read_var<long long> (const std::string&, const std::string&, long long*, const int);
template void read_var<float>     (const std::string&, const std::string&, float*,     const int);
template void read_var<double>    (const std::string&, const std::string&, double*,    const int);

template void write_var<int>       (const std::string&, const std::string&, const int*,       const int*);
template void write_var<long long> (const std::string&, const std::string&, const long long*, const long long*);
template void write_var<float>     (const std::string&, const std::string&, const float*,     const float*);
template void write_var<double>    (const std::string&, const std::string&, const double*,    const double*);

// =============== Attributes operations ================== //

// To specify GLOBAL attributes, pass "GLOBAL" as varname
ekat::any get_any_attribute (const std::string& filename,
                             const std::string& varname,
                             const std::string& attname)
{
  register_file(filename,Read);
  const int ncid = get_file(filename,"scorpio::get_any_attribute").ncid;

  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    const auto& var = get_var(filename,varname,"scorpio::get_any_attribute");
    varid = var.ncid;
  }
  int err;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    err = PIOc_inq_varid(ncid, varname.c_str(), &varid);
    EKAT_REQUIRE_MSG (err==PIO_NOERR,
        "[get_any_attribute] Error! Something went wrong while inquiring variable id.\n"
          " - filename : " + filename + "\n"
          " - variable : " + varname + "\n"
          " - attribute: " + attname + "\n"
          " - pio error: " << err << "\n");
  }

  nc_type type;
  PIO_Offset len;
  err = PIOc_inq_att(ncid,varid,attname.c_str(),&type,&len);
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring attribute.\n"
        " - filename : " + filename + "\n"
        " - variable : " + varname + "\n"
        " - attribute: " + attname + "\n"
        " - pio error: " << err << "\n");

  EKAT_REQUIRE_MSG (len==1 || type==PIO_CHAR,
      "[get_any_attribute] Error! Only single value attributes allowed.\n"
        " - filename : " + filename + "\n"
        " - variable : " + varname + "\n"
        " - attribute: " + attname + "\n"
        " - nc type  : " << type << "\n"
        " - att len  : " << len << "\n");

  ekat::any att;
  if (type==PIO_INT) {
    int val;
    err = PIOc_get_att(ncid,varid,attname.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_DOUBLE) {
    double val;
    err = PIOc_get_att(ncid,varid,attname.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_FLOAT) {
    float val;
    err = PIOc_get_att(ncid,varid,attname.c_str(),&val);
    att.reset(val);
  } else if (type==PIO_CHAR) {
    std::string val(len,'\0');
    err = PIOc_get_att(ncid,varid,attname.c_str(),&val[0]);
    att.reset(val);
  } else {
    EKAT_ERROR_MSG ("[get_any_attribute] Error! Unsupported/unrecognized nc type.\n"
        " - filename : " + filename + "\n"
        " - variable : " + varname + "\n"
        " - attribute: " + attname + "\n"
        " - nc type  : " << type << "\n");
  }
  EKAT_REQUIRE_MSG (err==PIO_NOERR,
      "[get_any_attribute] Error! Something went wrong while inquiring attribute.\n"
        " - filename : " + filename + "\n"
        " - variable : " + varname + "\n"
        " - attribute: " + attname + "\n"
        " - pio error: " << err << "\n");

  release_file(filename);
  return att;
}

void set_any_attribute (const std::string& filename,
                        const std::string& varname,
                        const std::string& attname,
                        const ekat::any& att)
{
  const auto& f = get_file (filename,"scorpio::set_any_attribute");
  const int ncid = f.ncid;

  int varid;
  if (varname=="GLOBAL") {
    varid = PIO_GLOBAL;
  } else {
    const auto& var = get_var(filename,varname,"scorpio::set_any_attribute");
    varid = var.ncid;
  }
  int err;

  const bool needs_redef = f.enddef;
  if (needs_redef) {
    redef(filename);
  }

  if (att.isType<int>()) {
    const int& data = ekat::any_cast<int>(att);
    err = PIOc_put_att(ncid,varid,attname.c_str(),PIO_INT,1,&data);
  } else if (att.isType<double>()) {
    const double& data = ekat::any_cast<double>(att);
    err = PIOc_put_att(ncid,varid,attname.c_str(),PIO_DOUBLE,1,&data);
  } else if (att.isType<float>()) {
    const float& data = ekat::any_cast<float>(att);
    err = PIOc_put_att(ncid,varid,attname.c_str(),PIO_FLOAT,1,&data);
  } else if (att.isType<std::string>()) {
    const std::string& data = ekat::any_cast<std::string>(att);
    err = PIOc_put_att(ncid,varid,attname.c_str(),PIO_CHAR,data.size(),data.data());
  } else {
    EKAT_ERROR_MSG ("[set_any_attribute] Error! Unsupported/unrecognized att type.\n"
        " - filename : " + filename + "\n"
        " - att name : " + attname + "\n"
        " - att value: " << att << "\n"
        " - type info: " << att.content().type().name() << "\n");
  }

  check_scorpio_noerr(err,filename,"attribute",attname,"set_any_attribute","put_att");

  if (needs_redef) {
    enddef(filename);
  }
}

} // namespace scorpio
} // namespace scream
