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

void GetObsOpt_vs_sqrtS(      const HepMC::GenVertex & signal, double * values, size_t count, const char * coefName );
void GetObsOptDivSN_vs_sqrtS( const HepMC::GenVertex & signal, double * values, size_t count, const char * coefName, int order );

void FillHistSqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal );
void FillHistOpt(   TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName );

void FillHistOpt_vs_sqrtS(      TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName );
void FillHistOptDivSN_vs_sqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName, int order );

////////////////////////////////////////////////////////////////////////////////

void CalcEvalVector( const RootUtil::CStringVector & coefNames, const ParamVector & params, std::vector<double> & evals );

TH1D * ReweightHist( const TH1D & sourceData, const RootUtil::ConstTH1DVector & sourceCoefs,
                     const std::vector<double> & sourceEval, const std::vector<double> & targetEval,
                     const char * name = nullptr, const char * title = nullptr );

////////////////////////////////////////////////////////////////////////////////

void LoadReweightFiles( // inputs:
                        const ModelCompare::ObservableVector & observables,
                        const RootUtil::CStringVector & coefNames,
                        const ModelCompare::ModelFile & targetFile,
                        const ModelCompare::ModelFile & sourceFile, const ParamVector & sourceParam,
                        // outputs:
                        RootUtil::TH1DVector &              targetData,     // targetData[observable]
                        RootUtil::TH1DVector &              sourceData,     // sourceData[observable]
                        std::vector<RootUtil::TH1DVector> & sourceCoefs,    // sourceCoefs[observable][coefficient]
                        std::vector<double> &               sourceEval,     // sourceEval[coefficient]
                        RootUtil::TH1DVector &              rawTargetData,  // rawTargetData[observable]
                        RootUtil::TH1DVector &              rawSourceData,  // rawSourceData[observable]
                        // optionals
                        double luminosity = 1.0,
                        const char * cacheFilename = nullptr
                        );

void ReweightEFT( const char * outputFileName,
                  const ModelCompare::ObservableVector & observables,
                  const RootUtil::CStringVector & coefNames,
                  const ModelCompare::ModelFile & targetFile, const ParamVector & targetParam,
                  const ModelCompare::ModelFile & sourceFile, const ParamVector & sourceParam,
                  double luminosity = 1.0,
                  const char * cacheFilename = nullptr );

////////////////////////////////////////////////////////////////////////////////

} // namespace ReweightEFT

////////////////////////////////////////////////////////////////////////////////

#endif // REWEIGHT_EFT_H
