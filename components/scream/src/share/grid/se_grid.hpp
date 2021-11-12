#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"

namespace scream
{

/*
 * An enum for the type of SE grid.
 * For all SEGrid object, the dofs are the gauss points (GP) of a Spectral Element mesh.
 * For the DG case corresponding edge dofs on bordering elements are considered independent
 * (i.e., may have different values), while for the CG case they are assumed to be equal.
 * An important consequence, is that for CG grids, we are free to *read* the value of a dof
 * from any of the elements that share the dof.
 */
enum class SEType {
  DG,
  CG
};

class SEGrid : public AbstractGrid
{
public:

  // Creates a SE grid (either CG or DG)
  SEGrid (const std::string& grid_name,
          const int num_my_elements,
          const int num_gauss_pts,
          const int num_vertical_levels,
          const SEType se_type,
          const ekat::Comm& comm);

  virtual ~SEGrid () = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  FieldLayout get_2d_scalar_layout () const override;
  FieldLayout get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const override;
  FieldLayout get_3d_scalar_layout (const bool midpoints) const override;
  FieldLayout get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const override;


  SEType get_se_type () const { return m_se_type; }

  // Set/retrieve the CG grid
  void set_cg_grid (const std::shared_ptr<const SEGrid>& cg_grid);
  std::shared_ptr<const AbstractGrid> get_cg_grid () const;

protected:
  void check_dofs_list () const override;
  void check_lid_to_idx_map () const override;
  void check_geo_data (const std::string& name, const geo_view_type& data) const override;

  // SE dims
  int       m_num_local_elem;
  int       m_num_gp;

  // SE type
  SEType m_se_type;

  // A CG version of this grid
  std::shared_ptr<const AbstractGrid> m_cg_grid;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
