/*! Definition of a tuple with summary information about production.
This file is part of https://github.com/cms-tau-pog/TauTriggerTools. */

#pragma once

#include "EventTuple.h"

#define SUMMARY_DATA() \
    /* Run statistics */ \
    VAR(UInt_t, exeTime) \
    VAR(ULong64_t, numberOfProcessedEvents) \
    VAR(Double_t, totalGenEventWeight) \
    /* Trigger information */ \
    VAR(std::vector<UInt_t>, trigger_index) \
    VAR(std::vector<std::string>, trigger_pattern) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(tau_trigger, ProdSummary, SummaryTuple, SUMMARY_DATA, "summary")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(tau_trigger, SummaryTuple, SUMMARY_DATA)
#undef VAR
#undef SUMMARY_DATA


#define EVENT_EXPRESS_DATA() \
    VAR(UInt_t, run) /* run */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, npv) /* number of primary vertices */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(tau_trigger, ExpressEvent, ExpressTuple, EVENT_EXPRESS_DATA, "all_events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(tau_trigger, ExpressTuple, EVENT_EXPRESS_DATA)
#undef VAR
#undef SUMMARY_DATA
