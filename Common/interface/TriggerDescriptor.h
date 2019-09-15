/*! Definition of trigger results.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include <bitset>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace tau_trigger {

struct TriggerLeg {
    int type{0};
    std::vector<std::string> filters;
};

struct TriggerDescriptor {
    std::string path;
    int global_index{-1};
    std::vector<TriggerLeg> legs;
    bool is_tag{false};
};

using TriggerBitsContainer = unsigned long long;
constexpr size_t MaxNumberOfTriggers = std::numeric_limits<TriggerBitsContainer>::digits;
using TriggerResults = std::bitset<MaxNumberOfTriggers>;

class TriggerDescriptorCollection {
public:
    TriggerDescriptorCollection(const edm::VParameterSet& trig_pset);
    const std::vector<TriggerDescriptor>& getDescriptors() const { return descs; }
    const TriggerDescriptor& at(size_t n) const { return descs.at(n); }
    size_t getIndex(const std::string& path) const { return desc_indices.at(path); }
    size_t size() const { return descs.size(); }

    void updateGlobalIndices(const std::vector<std::string>& triggerNames);

private:
    std::vector<TriggerDescriptor> descs;
    std::map<std::string, size_t> desc_indices;
};

}
