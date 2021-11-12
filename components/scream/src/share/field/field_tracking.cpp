#include "share/field/field_tracking.hpp"

namespace scream {

void FieldTracking::add_provider (const std::weak_ptr<AtmosphereProcess>& provider) {
  // In case of 'owners', we do not allow any other atm proc to be a provider
  EKAT_REQUIRE_MSG (m_owner.expired(),
      "Error! Cannot add 'providers' if the field has a 'owner' already.\n");
  auto p = this->get_parent().lock();
  EKAT_REQUIRE_MSG(not p or p->get_owner().expired(),
      "Error! Cannot set the provider of a field that is the child of a field with an owner.\n");
  m_providers.insert(provider);
}

void FieldTracking::add_customer (const std::weak_ptr<AtmosphereProcess>& customer) {
  m_customers.insert(customer);
}

void FieldTracking::set_owner (const std::weak_ptr<AtmosphereProcess>& owner) {
  // In case of 'owners', we do not allow any other atm proc to be a provider
  EKAT_REQUIRE_MSG (m_providers.size()==0,
      "Error! Cannot set a field 'owner' if there are other 'providers'.\n");
  // Also, there can only be *one* owner
  EKAT_REQUIRE_MSG (m_owner.expired(),
      "Error! Cannot re-set a field owner once it's been set.\n");
  // Also, this cannot be the child of another field or the parent of a field with providers
  EKAT_REQUIRE_MSG (this->get_parent().expired(),
      "Error! Cannot set owner of a field that is the child of another field.\n");
  for (auto it : this->get_children()) {
    auto c = it.lock();
    EKAT_REQUIRE_MSG(c, "Error! A weak pointer of a child field expired.\n");
    EKAT_REQUIRE_MSG(c->get_providers().size()==0,
        "Error! Cannot set the owner of a field that is the parent of fields with providers.\n");
  }

  m_owner = owner;
}

void FieldTracking::
add_to_group (const std::shared_ptr<const FieldGroupInfo>& group) {
  m_groups.insert(group);
}

void FieldTracking::update_time_stamp (const TimeStamp& ts) {
  // We check that the given time stamp is not in the past.
  // This is to prevent users from tampering with time stamps (e.g., rewinding time).
  EKAT_REQUIRE_MSG(!m_time_stamp.is_valid() || !(ts<m_time_stamp),
      "Error! Input time stamp is in the past.\n");

  m_time_stamp = ts;

  // If you update a field, all its subviews will automatically be updated
  for (auto it : this->get_children()) {
    auto c = it.lock();
    EKAT_REQUIRE_MSG(c, "Error! A weak pointer of a child field expired.\n");
    c->update_time_stamp(ts);
  }
}

} // namespace scream
