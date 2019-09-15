/*! Definition of trigger results.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#include "TauTriggerTools/Common/interface/TriggerDescriptor.h"
#include "TauTriggerTools/Common/interface/TextIO.h"

namespace tau_trigger {
TriggerDescriptorCollection::TriggerDescriptorCollection(const edm::VParameterSet& trig_vpset)
{
    if(trig_vpset.size() > MaxNumberOfTriggers)
        throw analysis::exception("The max number of triggers is exceeded");
    for(const auto& pset : trig_vpset) {
        TriggerDescriptor desc;
        desc.path = pset.getParameter<std::string>("path");
        desc.is_tag = pset.getParameter<bool>("is_tag");
        const std::vector<int> leg_types = pset.getParameter<std::vector<int>>("leg_types");
        const std::vector<std::string> filters = pset.getParameter<std::vector<std::string>>("filters");
        if(desc_indices.count(desc.path))
            throw analysis::exception("Duplicated trigger path = '%1%'.") % desc.path;
        if(leg_types.size() != filters.size())
            throw analysis::exception("Inconsitent leg_types and filters for trigger path = '%1%'.") % desc.path;

        for(size_t n = 0; n < leg_types.size(); ++n) {
            TriggerLeg leg;
            leg.type = leg_types.at(n);
            leg.filters = analysis::SplitValueList(filters.at(n), false);
            desc.legs.push_back(leg);
        }
        descs.push_back(desc);
        desc_indices[desc.path] = descs.size() - 1;
    }
}

void TriggerDescriptorCollection::updateGlobalIndices(const std::vector<std::string>& triggerNames)
{
    // TODO
}

} // namespace tau_trigger
