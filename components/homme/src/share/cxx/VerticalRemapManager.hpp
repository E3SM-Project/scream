/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VERTICAL_REMAP_MANAGER_HPP
#define HOMMEXX_VERTICAL_REMAP_MANAGER_HPP

#include <memory>

namespace Homme {

class FunctorsBuffersManager;

struct VerticalRemapManager {
  VerticalRemapManager();

  VerticalRemapManager(const int num_elems);

  void run_remap(int np1, int np1_qdp, double dt) const;

  int requested_buffer_size () const;
  void init_buffers(const FunctorsBuffersManager& fbm);

  bool setup_needed () { return !is_setup; }

  void setup ();

private:
  struct Impl;
  std::shared_ptr<Impl> p_;

  int m_num_elems;
  bool is_setup;
};

}

#endif
