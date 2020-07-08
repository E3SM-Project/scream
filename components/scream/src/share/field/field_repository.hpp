#ifndef SCREAM_FIELD_REPOSITORY_HPP
#define SCREAM_FIELD_REPOSITORY_HPP

#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/string_utils.hpp"
#include "ekat/util/scream_std_utils.hpp"
#include "share/field/field.hpp"

#include <map>
#include <set>

namespace scream
{

 /*
  *  A database for all the persistent fields needed in an atm time step
  *  We template a field repository over the field value type and over
  *  the device type. These two are enough to fully deduce the type of
  *  the stored views.
  *
  *  The fields are internally organized by name. Within each name,
  *  there can be multiple fields, which differ by their layout.
  *  For instance, we could have two version of 'temperature',
  *  one on the "physics" grid (tags: Column, Level), and one on
  *  the "dynamics" grid (tags: Element, GaussPoint, GaussPoint, Level).
  *
  *  When you query the repo for its 'size', you get the number
  *  of different *names*. So if the repo is storing the two
  *  versions of 'temperature' above and nothing else, its size
  *  will be 1. To get the number of different Field objects
  *  stored, you need to use the 'internal_size' method.
  *
  */

template<typename ScalarType, typename Device>
class FieldRepository {
public:

  // Public types
  using device_type     = Device;
  using scalar_type     = ScalarType;
  using field_type      = Field<scalar_type,device_type>;
  using header_type     = typename field_type::header_type;
  using identifier_type = typename header_type::identifier_type;
  using ci_string       = typename identifier_type::ci_string;
  using alias_map_type  = std::map<identifier_type,field_type>;
  using repo_type       = std::map<ci_string,alias_map_type>;
  using groups_map_type = std::map<ci_string,std::set<ci_string>>;

  // Constructor(s)
  FieldRepository ();

  // No copies, cause the internal database is not a shared_ptr.
  // NOTE: you can change this if you find that copies are needed/useful.
  FieldRepository (const FieldRepository&) = delete;
  FieldRepository& operator= (const FieldRepository&) = delete;

  // Change the state of the database
  void registration_begins ();
  void registration_ends ();
  void clean_up ();

  // Deduce the pack size from the scalar type (which must be of type Pack<ScalarType,N>, for some int N>0, or ScalarType)
  template<typename RequestedValueType = scalar_type>
  void register_field (const identifier_type& identifier, const std::set<std::string>& groups_names);

  template<typename RequestedValueType = scalar_type>
  void register_field (const identifier_type& identifier, const std::string& field_group);

  template<typename RequestedValueType = scalar_type>
  void register_field (const identifier_type& identifier);

  // Get information about the state of the repo
  int size () const { return m_fields.size(); }
  int internal_size () const;
  RepoState repository_state () const { return m_repo_state; }

  // Query for a particular field or group of fields
  bool has_field (const identifier_type& identifier) const;
  const field_type& get_field (const identifier_type& identifier) const;
  const groups_map_type& get_field_groups () const { return m_field_groups; }

  // Iterators, to allow range for loops over the repo.
  typename repo_type::const_iterator begin() const { return m_fields.begin(); }
  typename repo_type::const_iterator end()   const { return m_fields.end(); }

  typename repo_type::iterator begin() { return m_fields.begin(); }
  typename repo_type::iterator end()   { return m_fields.end(); }

protected:

  // The state of the repository
  RepoState           m_repo_state;

  // The actual repo.
  repo_type           m_fields;

  // The groups
  groups_map_type     m_field_groups;

  // The reserved groups. Users are not allowed to manually add fields to these groups
  std::set<ci_string> m_reserved_groups;
};

// ============================== IMPLEMENTATION ============================= //

template<typename ScalarType, typename Device>
FieldRepository<ScalarType,Device>::FieldRepository ()
 : m_repo_state (RepoState::Clean)
{
  m_reserved_groups.insert("state");
  m_reserved_groups.insert("old state");

  // TODO: we should add names to the 'state' and 'old state' group.
  //       This means that the repo should know the name of the state variables,
  //       as well as the names of the old state variables (probably "blah old").
  //       This may require passing a parameter list or a list of strings to the constructor.
}

template<typename ScalarType, typename Device>
template<typename RequestedValueType>
void FieldRepository<ScalarType,Device>::register_field (const identifier_type& id) {
  std::set<std::string> empty_set;
  register_field<RequestedValueType>(id,empty_set);
}

template<typename ScalarType, typename Device>
template<typename RequestedValueType>
void FieldRepository<ScalarType,Device>::
register_field (const identifier_type& id, const std::string& group_name) {
  std::set<std::string> group_name_set;
  group_name_set.insert(group_name);
  register_field<RequestedValueType>(id,group_name_set);
}

template<typename ScalarType, typename Device>
template<typename RequestedValueType>
void FieldRepository<ScalarType,Device>::
register_field (const identifier_type& id, const std::set<std::string>& groups_names) {

  // Check that the requested value type is either ScalarType or Pack<ScalarType,N>, for some N>0.
  static_assert(std::is_same<ScalarType,RequestedValueType>::value ||
                std::is_same<ScalarType,typename util::ScalarProperties<RequestedValueType>::scalar_type>::value,
                "Error! The template argument 'RequestedValueType' of this function must either match "
                "the template argument 'ScalarType' of this class or be a Pack type based on ScalarType.\n");

  // Sanity checks
  scream_require_msg(m_repo_state==RepoState::Open,"Error! Registration of new fields not started or no longer allowed.\n");

  // Get the map of all fields with this name
  auto& map = m_fields[id.name()];

  if (map.size()>0) {
    using units::to_string;
    // Make sure the new field id stores the same units as all the other ones.
    // TODO: this is the easiest way to ensure everyone uses the same units.
    //       However, in the future, we *may* allow different units, providing
    //       the users with conversion routines perhaps.
    scream_require_msg(id.get_units()==map.begin()->first.get_units(),
                       "Error! Request to register field '" + id.name() + "' with units '" +
                       to_string(id.get_units()) + "',\n"
                       "       but there is already a request for the same field with units '" +
                       to_string(map.begin()->first.get_units()) + "'\n"
                       "       Please, check and make sure all atmosphere processes use the same units.\n");
  }

  // Try to create the field. Allow case where it is already existing.
  auto it_bool = map.emplace(id,field_type(id));

  // Make sure the field can accommodate the requested value type
  it_bool.first->second.get_header().get_alloc_properties().template request_value_type_allocation<RequestedValueType>();

  // Finally, add the field to the given groups
  for (const auto& group_name : groups_names) {
    // First, make sure it's not a reserved group
    scream_require_msg(!util::contains(m_reserved_groups,group_name),"");

    // Add the group name to the field tracking of all the fields with that name
    // Remember: fields with the same name can differ only because of tags/extents (i.e., different grids).
    //           Morally, they are different layouts of the same field.
    for (auto& f_it : map) {
      f_it.second.get_header().get_tracking().add_to_group(group_name);
    }

    // Add the field name to the set of fields belonging to this group
    m_field_groups[group_name].insert(id.name());
  }
}

template<typename ScalarType, typename Device>
int FieldRepository<ScalarType,Device>::
internal_size () const {
  int s = 0;
  for (auto x : m_fields) {
    s+= x.second.size();
  }
  return s;
}

template<typename ScalarType, typename Device>
bool FieldRepository<ScalarType,Device>::
has_field (const identifier_type& identifier) const {
  auto it = m_fields.find(identifier.name());
  return it!=m_fields.end() && it->second.find(identifier)!=it->second.end();
}

template<typename ScalarType, typename Device>
const typename FieldRepository<ScalarType,Device>::field_type&
FieldRepository<ScalarType,Device>::get_field (const identifier_type& id) const {
  scream_require_msg(m_repo_state==RepoState::Closed,"Error! You are not allowed to grab fields from the repo until after the registration phase is completed.\n");

  const auto& map = m_fields.find(id.name());
  scream_require_msg(map!=m_fields.end(), "Error! Field not found.\n");
  auto it = map->second.find(id);
  scream_require_msg(it!=map->second.end(), "Error! Field not found.\n");
  return it->second;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::registration_begins () {
  // Update the state of the repo
  m_repo_state = RepoState::Open;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::registration_ends () {
  // Proceed to allocate fields
  for (auto& map : m_fields) {
    for (auto& it : map.second) {
      it.second.allocate_view();
    }
  }

  // Prohibit further registration of fields
  m_repo_state = RepoState::Closed;
}

template<typename ScalarType, typename Device>
void FieldRepository<ScalarType,Device>::clean_up() {
  m_fields.clear();
  m_repo_state = RepoState::Clean;
}

} // namespace scream

#endif // SCREAM_FIELD_REPOSITORY_HPP
