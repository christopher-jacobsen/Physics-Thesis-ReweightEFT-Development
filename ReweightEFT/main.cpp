//
//  main.cpp
//  ReweightEFT
//
//  Created by Christopher Jacobsen on 11/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "ReweightEFT.h"
#include "ModelCompare.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////

#define OBSFILL [](TH1D & h, double w, const HepMC::GenVertex & s)

static const ModelCompare::ObservableVector Observables1 =
{
  //{ "PTZ", "P_{T}(Z)",  300,     0, 1500, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       OBSFILL{ RootUtil::FillHistPT(  h, w, s, 24);     } },
  //{ "MWZ", "M(WZ)",     250,     0, 3750, "M(WZ) [GeV/c^2]",  "Events per 15 GeV/c^2",    OBSFILL{ RootUtil::FillHistMass(h, w, s, 24, 23); } },
  //{ "PTZ", "P_{T}(Z)",  100,     0,  500, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       OBSFILL{ RootUtil::FillHistPT(  h, w, s, 24);     } },
  //{ "MWZ", "M(WZ)",     100,     0, 1500, "M(WZ) [GeV/c^2]",  "Events per 15 GeV/c^2",    OBSFILL{ RootUtil::FillHistMass(h, w, s, 24, 23); } },
    { "PTZ", "P_{T}(Z)",  150,     0,  750, "P_{T}(Z) [GeV/c]", "Events per 5 GeV/c",       OBSFILL{ RootUtil::FillHistPT(  h, w, s, 24);     } },
    { "MWZ", "M(WZ)",     140,     0, 2800, "M(WZ) [GeV/c^2]",  "Events per 20 GeV/c^2",    OBSFILL{ RootUtil::FillHistMass(h, w, s, 24, 23); } },
    { "RAZ", "Y(Z)",      100,    -5,    5, "Y(Z)",             "Events per bin",           OBSFILL{ RootUtil::FillHistRap( h, w, s, 24);     } },
    { "ETZ", "#eta(Z)",   100,   -10,   10, "#eta(Z)",          "Events per bin",           OBSFILL{ RootUtil::FillHistEta( h, w, s, 24);     } },
    { "PHZ", "#phi(Z)",   100, -M_PI, M_PI, "#phi(Z)",          "Events per bin",           OBSFILL{ RootUtil::FillHistPhi( h, w, s, 24);     } },
};

////////////////////////////////////////////////////////////////////////////////

static const ModelCompare::ModelFileVector Models_1E4 =
{
    { "events/SM_220_1E4.hepmc2g.gz",           "SM_220",       "SM (2.2.0)",   18.2613, 0.154079, 10000 },
    { "weight/EFT_all_weight_1E4.hepmc2g.gz",   "EFT_all",      "EFT (all)",    42.9091, 0.379962, 10000 },
};

static const ModelCompare::ModelFileVector Models_1E6 =
{
    { "events/SM_220_1E6.hepmc2g.gz",           "SM_220",       "SM (2.2.0)",   18.5537, 0.0156025, 1000000 },
    { "weight/EFT_all_weight_1E6.hepmc2g.gz",   "EFT_all",      "EFT (all)",    42.8051, 0.0379205, 1000000 },
};

////////////////////////////////////////////////////////////////////////////////

static const RootUtil::CStringVector CoefNames_EFT_all =
{
    "F_0_0",
    "F_0_1_ocWWW",
    "F_0_2_ocW",
    "F_0_3_ocB",
    "F_1_1_ocWWW",
    "F_1_2_ocWWW_ocW",
    "F_1_3_ocWWW_ocB",
    "F_2_2_ocW",
    "F_2_3_ocW_ocB",
    "F_3_3_ocB",
};

static const ReweightEFT::ParamVector Params_EFT_all =
{
    { "ocWWW", 3E-5 },
    { "ocW",   5E-5 },
    { "ocB",   9E-4 },
};

static const ReweightEFT::ParamVector Params_EFT_SM =
{
    { "ocWWW", 0 },
    { "ocW",   0 },
    { "ocB",   0 },
};

////////////////////////////////////////////////////////////////////////////////
int main(void)
{
    //ReweightEFT::ReweightEFT( "reweight/EFT_all_to_self_1E4.root", Observables1, Models_1E4[1], Models_1E4[1], Params_EFT_all, Params_EFT_all, CoefNames_EFT_all );

    //ReweightEFT::ReweightEFT( "reweight/EFT_all_to_SM_1E4.root", Observables1, Models_1E4[1], Models_1E4[0], Params_EFT_all, Params_EFT_SM, CoefNames_EFT_all );
    ReweightEFT::ReweightEFT( "reweight/EFT_all_to_SM_1E6.root", Observables1, Models_1E6[1], Models_1E6[0], Params_EFT_all, Params_EFT_SM, CoefNames_EFT_all );

    LogMsgInfo( "Done." );
    return 0;
}
