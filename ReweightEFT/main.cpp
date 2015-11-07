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

using namespace RootUtil;
using namespace ModelCompare;
using namespace ReweightEFT;

////////////////////////////////////////////////////////////////////////////////
/*
static const ModelCompare::ObservableVector Observables1 =
{
    // phase-space observables

    { "PTZ", "P_{T}(Z)",  150,     0,  750, "P_{T}(Z) [GeV/c]",     "Events per 5 GeV/c",       GETOBS{ GetObs(s,v,c, GetObsPT,   24);     } },
    { "MWZ", "M(WZ)",     150,     0, 3000, "M(WZ) [GeV/c^{2}]",    "Events per 20 GeV/c^{2}",  GETOBS{ GetObs(s,v,c, GetObsMass, 24, 23); } },
    { "RAZ", "Y(Z)",      100,    -5,    5, "Y(Z)",                 "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsRap,  24);     } },
//  { "ETZ", "#eta(Z)",   100,   -10,   10, "#eta(Z)",              "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsEta,  24);     } },
//  { "PHZ", "#phi(Z)",   100, -M_PI, M_PI, "#phi(Z)",              "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsPhi,  24);     } },

    // optimal observables

    { "cWWW_O1",    "O_{1}(cWWW)",  1000,  -6E4,  1E4,  "O_{1}(cWWW) [GeV^{2}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_1_ocWWW"); } },
    { "cWWW_O2",    "O_{2}(cWWW)",  1000,   0,   14E11, "O_{2}(cWWW) [GeV^{4}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_1_1_ocWWW"); } },

    { "cW_O1",      "O_{1}(cW)",    1000, -12E5,  1E5,  "O_{1}(cW) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_2_ocW");   } },
    { "cW_O2",      "O_{2}(cW)",    1000,   0,   54E10, "O_{2}(cW) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_2_2_ocW");   } },

    { "cB_O1",      "O_{1}(cB)",    1000,  -1E3,  5E3,  "O_{1}(cB) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_3_ocB");   } },
    { "cB_O2",      "O_{2}(cB)",    1000,   0,   13E7,  "O_{2}(cB) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_3_3_ocB");   } },

    // mean optimal observables

//  { "s",          "#sqrt{s}",     100, 0, 3000,   "#sqrt{s} [GeV]", "Events per bin",             GETOBS{ GetObs(s,v,c, GetObsSqrtS); } },  // same as M(WZ)

    { "cWWW_O1vS",  "O_{1}(cWWW)",  100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cWWW) [GeV^{2}]", GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_1_ocWWW"); },    2, DefaultTProfileFactory },
    { "cWWW_O2vS",  "O_{2}(cWWW)",  100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cWWW) [GeV^{4}]", GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_1_1_ocWWW"); },    2, DefaultTProfileFactory },

    { "cW_O1vS",    "O_{1}(cW)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cW) [GeV^{2}]",   GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_2_ocW"); },      2, DefaultTProfileFactory },
    { "cW_O2vS",    "O_{2}(cW)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cW) [GeV^{4}]",   GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_2_2_ocW"); },      2, DefaultTProfileFactory },

    { "cB_O1vS",    "O_{1}(cB)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cB) [GeV^{2}]",   GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_3_ocB"); },      2, DefaultTProfileFactory },
    { "cB_O2vS",    "O_{2}(cB)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cB) [GeV^{4}]",   GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_3_3_ocB"); },      2, DefaultTProfileFactory },

    //

    { "cWWW_O1dSvS",    "O_{1}(cWWW)",  100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cWWW)/s",             GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_1_ocWWW", 1); },    2, DefaultTProfileFactory },
//  { "cWWW_O2dSvS",    "O_{2}(cWWW)",  100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cWWW)/s [GeV^{2}]",   GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_1_1_ocWWW", 1); },    2, DefaultTProfileFactory },
    { "cWWW_O2dS2vS",   "O_{2}(cWWW)",  100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cWWW)/s^{2}",         GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_1_1_ocWWW", 2); },    2, DefaultTProfileFactory },

    { "cW_O1dSvS",      "O_{1}(cW)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cW)/s",               GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_2_ocW",   1); },    2, DefaultTProfileFactory },
//  { "cW_O2dSvS",      "O_{2}(cW)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cW)/s [GeV^{2}]",     GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_2_2_ocW",   1); },    2, DefaultTProfileFactory },
    { "cW_O2dS2vS",     "O_{2}(cW)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cW)/s^{2}",           GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_2_2_ocW",   2); },    2, DefaultTProfileFactory },

    { "cB_O1dSvS",      "O_{1}(cB)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{1}(cB)/s",               GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_3_ocB",   1); },    2, DefaultTProfileFactory },
//  { "cB_O2dSvS",      "O_{2}(cB)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cB)/s [GeV^{2}]",     GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_3_3_ocB",   1); },    2, DefaultTProfileFactory },
    { "cB_O2dS2vS",     "O_{2}(cB)",    100, 0, 3000,   "#sqrt{s} [GeV]", "Mean O_{2}(cB)/s^{2}",           GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_3_3_ocB",   2); },    2, DefaultTProfileFactory },

};
*/
////////////////////////////////////////////////////////////////////////////////

// optimal bins
static const ModelCompare::ObservableVector Observables2 =
{

    // phase-space observables

    { "PTZ",        "P_{T}(Z)",      750,      0,    750,   "P_{T}(Z) [GeV/c]",     "Events per GeV/c",         GETOBS{ GetObs(s,v,c, GetObsPT,   24);     } },
    { "MWZ",        "M(WZ)",        1500,      0,   3000,   "M(WZ) [GeV/c^{2}]",    "Events per 2 GeV/c^{2}",   GETOBS{ GetObs(s,v,c, GetObsMass, 24, 23); } },
    { "RAZ",        "Y(Z)",          200,     -5,      5,   "Y(Z)",                 "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsRap,  24);     } },
//  { "ETZ",        "#eta(Z)",       100,    -10,     10,   "#eta(Z)",              "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsEta,  24);     } },
//  { "PHZ",        "#phi(Z)",       100,  -M_PI,   M_PI,   "#phi(Z)",              "Events per bin",           GETOBS{ GetObs(s,v,c, GetObsPhi,  24);     } },

    // optimal observables

    { "cWWW_O1",    "O_{1}(cWWW)",     2000,  -6E4,  1E4,   "O_{1}(cWWW) [GeV^{2}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_1_ocWWW"); } },
    { "cWWW_O2",    "O_{2}(cWWW)",    10000,   0,   14E11,  "O_{2}(cWWW) [GeV^{4}]", "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_1_1_ocWWW"); } },

    { "cW_O1",      "O_{1}(cW)",      10000, -12E5,  1E5,   "O_{1}(cW) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_2_ocW");   } },
    { "cW_O2",      "O_{2}(cW)",      10000,   0,   54E10,  "O_{2}(cW) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_2_2_ocW");   } },

    { "cB_O1",      "O_{1}(cB)",       1000,  -1E3,  5E3,   "O_{1}(cB) [GeV^{2}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_0_3_ocB");   } },
    { "cB_O2",      "O_{2}(cB)",      10000,   0,   13E7,   "O_{2}(cB) [GeV^{4}]",   "Events per bin",  GETOBS{ GetObs(s,v,c, GetObsOpt, "F_3_3_ocB");   } },

    // mean optimal observables

//  { "s",              "#sqrt{s}",     1500, 0, 3000,      "#sqrt{s} [GeV]", "Events per bin",                GETOBS{ GetObs(s,v,c, GetObsSqrtS); } },  // same as M(WZ)

    { "cWWW_O1vS",      "O_{1}(cWWW)",  1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cWWW) [GeV^{2}]",    GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_1_ocWWW"); },    2, DefaultTProfileFactory },
    { "cWWW_O2vS",      "O_{2}(cWWW)",  1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cWWW) [GeV^{4}]",    GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_1_1_ocWWW"); },    2, DefaultTProfileFactory },

    { "cW_O1vS",        "O_{1}(cW)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cW) [GeV^{2}]",      GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_2_ocW"); },      2, DefaultTProfileFactory },
    { "cW_O2vS",        "O_{2}(cW)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cW) [GeV^{4}]",      GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_2_2_ocW"); },      2, DefaultTProfileFactory },

    { "cB_O1vS",        "O_{1}(cB)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cB) [GeV^{2}]",      GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_0_3_ocB"); },      2, DefaultTProfileFactory },
    { "cB_O2vS",        "O_{2}(cB)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cB) [GeV^{4}]",      GETOBS{ GetObsOpt_vs_sqrtS(s,v,c, "F_3_3_ocB"); },      2, DefaultTProfileFactory },

    //

    { "cWWW_O1dSvS",    "O_{1}(cWWW)",  1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cWWW)/s",            GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_1_ocWWW", 1); },    2, DefaultTProfileFactory },
//  { "cWWW_O2dSvS",    "O_{2}(cWWW)",  1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cWWW)/s [GeV^{2}]",  GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_1_1_ocWWW", 1); },    2, DefaultTProfileFactory },
    { "cWWW_O2dS2vS",   "O_{2}(cWWW)",  1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cWWW)/s^{2}",        GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_1_1_ocWWW", 2); },    2, DefaultTProfileFactory },

    { "cW_O1dSvS",      "O_{1}(cW)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cW)/s",              GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_2_ocW",   1); },    2, DefaultTProfileFactory },
//  { "cW_O2dSvS",      "O_{2}(cW)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cW)/s [GeV^{2}]",    GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_2_2_ocW",   1); },    2, DefaultTProfileFactory },
    { "cW_O2dS2vS",     "O_{2}(cW)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cW)/s^{2}",          GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_2_2_ocW",   2); },    2, DefaultTProfileFactory },

    { "cB_O1dSvS",      "O_{1}(cB)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{1}(cB)/s",              GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_0_3_ocB",   1); },    2, DefaultTProfileFactory },
//  { "cB_O2dSvS",      "O_{2}(cB)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cB)/s [GeV^{2}]",    GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_3_3_ocB",   1); },    2, DefaultTProfileFactory },
    { "cB_O2dS2vS",     "O_{2}(cB)",    1500, 0, 3000,      "#sqrt{s} [GeV]", "Mean O_{2}(cB)/s^{2}",          GETOBS{ GetObsOptDivSN_vs_sqrtS(s,v,c, "F_3_3_ocB",   2); },    2, DefaultTProfileFactory },

};

////////////////////////////////////////////////////////////////////////////////

static const ModelCompare::ModelFileVector Models_1E4 =
{
    { "weight/SM_220_weight_1E4.hepmc2g.gz",    "SM_220",       "SM (2.2.0)",   18.2613, 0.154079, 10000 },
    { "weight/EFT_all_weight_1E4.hepmc2g.gz",   "EFT_all",      "EFT (all)",    42.9091, 0.379962, 10000 },
};

static const ModelCompare::ModelFileVector Models_1E6 =
{
  //{ "weight/SM_220_weight_1E6.hepmc2g.gz",        "SM_220",       "SM (2.2.0)",   18.5537, 0.0156025, 1000000 },
    { "weight/SM_UFO_220_weight_1E6.hepmc2g.gz",    "SM_UFO",       "SM-EFT",       18.5476, 0.0155949, 1000000 },
    { "weight/EFT_all_weight_1E6.hepmc2g.gz",       "EFT_all",      "EFT (all)",    42.8051, 0.0379205, 1000000 },
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

static const double Luminosity = 10.0;

////////////////////////////////////////////////////////////////////////////////
int main(void)
{
  //ReweightEFT::ReweightEFT( "reweight/EFT_all_to_self_1E4.root", Observables1, CoefNames_EFT_all, Models_1E4[1], Params_EFT_all, Models_1E4[1], Params_EFT_all );

//    ReweightEFT::ReweightEFT( "reweight/EFT_all_to_SM_1E4.root", Observables1, CoefNames_EFT_all,
//                              Models_1E4[0], Params_EFT_SM, Models_1E4[1], Params_EFT_all,
//                              Luminosity, "reweight/Cache_1E4.root" );

    ReweightEFT::ReweightEFT( "reweight/EFT_all_to_SM_1E6.root", Observables2, CoefNames_EFT_all,
                              Models_1E6[0], Params_EFT_SM, Models_1E6[1], Params_EFT_all,
                              Luminosity, "reweight/Cache_1E6.root" );

    LogMsgInfo( "\nDone." );
    return 0;
}
