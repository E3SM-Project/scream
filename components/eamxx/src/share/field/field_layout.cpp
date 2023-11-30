#include "field_layout.hpp"

#include <ekat/util/ekat_string_utils.hpp>

namespace scream
{

FieldLayout::FieldLayout (const std::initializer_list<FieldTag>& tags)
 : FieldLayout(std::vector<FieldTag>(tags))
{
  // Nothing to do here
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags)
 : FieldLayout(tags,std::vector<int>(tags.size(),-1))
{
  // Nothing to do here
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : FieldLayout (tags,evec2str(tags),dims)
{
  // Nothing to do here
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<std::string>& names,
                          const std::vector<int>& dims)
 : m_rank (tags.size())
 , m_tags (tags)
 , m_names(names)
 , m_dims (dims)
 , m_extents ("",tags.size())
{
  EKAT_REQUIRE_MSG (dims.size()==tags.size(),
      "Error! Tags and dims vectors dimensions mismatch.\n"
      "  tags size: " + std::to_string(tags.size()) + "\n"
      "  dims size: " + std::to_string(dims.size()) + "\n");

  EKAT_REQUIRE_MSG (names.size()==tags.size(),
      "Error! Tags and names vectors dimensions mismatch.\n"
      "  tags size : " + std::to_string(tags.size()) + "\n"
      "  names size: " + std::to_string(names.size()) + "\n");

  set_extents ();
  compute_type ();
}

bool FieldLayout::is_vector_layout () const {
  const auto lt = type();
  return lt==LayoutType::Vector2D || lt==LayoutType::Vector3D;
}

int FieldLayout::get_vector_dim () const {
  EKAT_REQUIRE_MSG (is_vector_layout(),
      "Error! 'get_vector_dim' available only for vector layouts.\n"
      "       Current layout: " + e2str(type()) + "\n");

  using namespace ShortFieldTagsNames;
  int idim = -1;
  if (has_tag(CMP)) {
    idim = std::distance(m_tags.begin(),ekat::find(m_tags,CMP));
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized layout for a '" + e2str(type()) + "' quantity.\n");
  }
  return idim;
}

FieldLayout FieldLayout::strip_dim (const FieldTag tag) const {
  auto it = ekat::find(m_tags,tag);

  // Check if found
  EKAT_REQUIRE_MSG(it!=m_tags.end(), "Error! Tag '" + e2str(tag) + "' not found.\n");

  // Check only one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_tags,tag)==1,
                     "Error! Tag '" + e2str(tag) + "' appears multiple times.\n"
                     "       You must inspect tags() and dims() manually.\n");

  return strip_dim (std::distance(m_tags.begin(),it));
}

FieldLayout FieldLayout::strip_dim (const int idim) const {
  EKAT_REQUIRE_MSG (idim>=0 and idim<m_rank,
      "Error! Cannot strip dimension, because it is out of bounds.\n"
      "  - input dim index: " + std::to_string(idim) + "\n"
      "  - layout rank    : " + std::to_string(m_rank) + "\n");
  std::vector<FieldTag>    t = tags();
  std::vector<int>         d = dims();
  std::vector<std::string> n = names();
  t.erase(t.begin()+idim);
  d.erase(d.begin()+idim);
  n.erase(n.begin()+idim);
  return FieldLayout (t,n,d);
}

void FieldLayout::set_dimension (const int idim, const int dimension) {
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  EKAT_REQUIRE_MSG(dimension>=0, "Error! Dimensions must be non-negative.");
  m_dims[idim] = dimension;

  // Recompute device extents
  set_extents ();
}

void FieldLayout::rename_dim (const int idim, const std::string& n) {
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");

  m_names[idim] = n;
}
void FieldLayout::rename_dim (const FieldTag tag, const std::string& n) {
  rename_dim(dim(tag),n);
}

void FieldLayout::set_extents () {
  auto extents_h = Kokkos::create_mirror_view(m_extents);
  std::copy_n(m_dims.begin(),m_rank,extents_h.data());
  Kokkos::deep_copy(m_extents,extents_h);
}

void FieldLayout::compute_type () {
  using namespace ShortFieldTagsNames;

  using ekat::erase;
  using ekat::count;

  auto tags = this->tags();

  const int n_element = count(tags,EL);
  const int n_column  = count(tags,COL);
  const int ngp       = count(tags,GP);
  const int nvlevs    = count(tags,LEV) + count(tags,ILEV);
  const int ncomps    = count(tags,CMP);

  if (n_element==1 && ngp==2 && n_column==0) {
    // A Dynamics layout

    // Remove the Element and the two GaussPoint tags
    erase(tags,EL);
    erase(tags,GP);
    erase(tags,GP);
  } else if (n_element==0 && ngp==0 && n_column==1) {
    // A Physics layout

    // Remove the column tag
    erase(tags,COL);
  } else if (tags.size()==0) {
    m_type = LayoutType::Scalar0D; return;
  } else if (tags.size()==1 and tags[0]==CMP) {
    m_type = LayoutType::Vector0D; return;
  } else if (tags.size()==1 and nvlevs==1) {
    m_type = LayoutType::Scalar1D; return;
  } else if (tags.size()==2 and ncomps==1 and nvlevs==1) {
    m_type = LayoutType::Vector1D; return;
  } else {
    // Not a supported layout.
    m_type = LayoutType::Invalid; return;
  }

  // Get the size of what's left
  const auto size = tags.size();
  switch (size) {
    case 0:
      m_type = LayoutType::Scalar2D;
      break;
    case 1:
      // The only tag left should be 'CMP', 'TL', or 'LEV'/'ILEV'
      if (tags[0]==CMP || tags[0]==TL) {
        m_type = LayoutType::Vector2D;
      } else if (tags[0]==LEV || tags[0]==ILEV) {
        m_type = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible supported scenarios:
      //  1) <CMP|TL,LEV|ILEV>
      //  2) <TL,CMP>
      if ( (tags[1]==LEV || tags[1]==ILEV) && (tags[0]==CMP || tags[0]==TL)) {
        m_type = LayoutType::Vector3D;
      } else if (tags[0]==TL && tags[1]==CMP ) {
        m_type = LayoutType::Tensor2D;
      }
      break;
    case 3:
      // The only supported scenario is:
      //  1) <TL,  CMP, LEV|ILEV>
      if ( tags[0]==TL && tags[1]==CMP &&
          (tags[2]==LEV || tags[2]==ILEV)) {
        m_type = LayoutType::Tensor3D;
      }
    default:
      // If nothing worked, this type is not recognized
      m_type = LayoutType::Invalid;
  }
}

std::string to_string (const FieldLayout& layout)
{
  if (layout.rank()==0) {
    return "<>()";
  }

  std::string s;
  s += "<" + e2str(layout.tags()[0]);
  for (int dim=1; dim<layout.rank(); ++dim) {
    s += "," + e2str(layout.tags()[dim]);
  }
  s += ">(" + ekat::join(layout.dims(),",") + ")";

  return s;
}

} // namespace scream
