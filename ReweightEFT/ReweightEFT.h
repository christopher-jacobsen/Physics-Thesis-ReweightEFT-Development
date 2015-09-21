//
//  ReweightEFT.h
//  ReweightEFT
//
//  Created by Christopher Jacobsen on 11/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#ifndef REWEIGHT_EFT_H
#define REWEIGHT_EFT_H

#include "ModelCompare.h"
#include "common.h"

////////////////////////////////////////////////////////////////////////////////

namespace ReweightEFT
{

////////////////////////////////////////////////////////////////////////////////

struct Parameter
{
    const char * name;
    double       value;
};

typedef std::vector<Parameter> ParamVector;

////////////////////////////////////////////////////////////////////////////////

double GetObsSqrtS( const HepMC::GenVertex & signal );
double GetObsOpt(   const HepMC::GenVertex & signal, const char * coefName );

void FillHistSqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal );
void FillHistOpt(   TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName );

void FillHistOpt_vs_sqrtS(     TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName );
void FillHistO2divS2_vs_sqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName );


void ReweightEFT( const char * outputFileName, const ModelCompare::ObservableVector & observables,
                  const ModelCompare::ModelFile & eventSource, const ModelCompare::ModelFile & eventTarget,
                  const ParamVector &             sourceParam, const ParamVector &             targetParam,
                  const RootUtil::CStringVector & coefNames );

////////////////////////////////////////////////////////////////////////////////

} // namespace ReweightEFT

////////////////////////////////////////////////////////////////////////////////

#endif // REWEIGHT_EFT_H
