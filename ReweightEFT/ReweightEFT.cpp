//
//  ReweightEFT.cpp
//  ReweightEFT
//
//  Created by Christopher Jacobsen on 11/09/15.
//  Copyright (c) 2015 Christopher Jacobsen. All rights reserved.
//

#include "ReweightEFT.h"
#include "ModelCompare.h"
#include "RootUtil.h"

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>

// HepMC includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/WeightContainer.h>

////////////////////////////////////////////////////////////////////////////////

using namespace RootUtil;

////////////////////////////////////////////////////////////////////////////////

namespace ReweightEFT
{

////////////////////////////////////////////////////////////////////////////////
double GetObsSqrtS( const HepMC::GenVertex & signal )
{
    TLorentzVector total;

    auto itrPart = signal.particles_in_const_begin();
    auto endPart = signal.particles_in_const_end();
    for ( ; itrPart != endPart; ++itrPart)
    {
        const HepMC::GenParticle * pPart = *itrPart;

        TLorentzVector vec = ToLorentz( pPart->momentum() );

        total += vec;
    }

    return total.M();
}

////////////////////////////////////////////////////////////////////////////////
double GetObsOpt( const HepMC::GenVertex & signal, const char * coefName )
{
    const HepMC::GenEvent * pEvent = signal.parent_event();
    if (!pEvent)
        ThrowError( "Signal has no parent event" );

    const HepMC::WeightContainer & weights = pEvent->weights();

    if (!weights.has_key( coefName ))
        ThrowError( "Coefficient " + std::string(coefName) + " not found in event." );
    if (!weights.has_key( "F_0_0" ))
        ThrowError( "Coefficient F_0_0 not found in event." );

    double coefWeight = weights[coefName];
    double coefBase   = weights["F_0_0"];

    if (coefBase == 0)
        ThrowError( "Coefficient F_0_0 is zero" );

    double result = coefWeight / coefBase;
    return result;
}

////////////////////////////////////////////////////////////////////////////////
void FillHistSqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal )
{
    hist.Fill( GetObsSqrtS( signal ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistOpt( TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName )
{
    hist.Fill( GetObsOpt( signal, coefName ), weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistOpt_vs_sqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName )
{
    TProfile * pProfile = dynamic_cast<TProfile *>( &hist );
    if (!pProfile)
        ThrowError( "FillHistRelCoef_SqrtS requires a TProfile" );

    double opt    = GetObsOpt(   signal, coefName );
    double sqrt_s = GetObsSqrtS( signal );

    pProfile->Fill( sqrt_s, opt, weight );
}

////////////////////////////////////////////////////////////////////////////////
void FillHistO2divS2_vs_sqrtS( TH1D & hist, double weight, const HepMC::GenVertex & signal, const char * coefName )
{
    TProfile * pProfile = dynamic_cast<TProfile *>( &hist );
    if (!pProfile)
        ThrowError( "FillHistRelCoef_SqrtS requires a TProfile" );

    double opt    = GetObsOpt(   signal, coefName );
    double sqrt_s = GetObsSqrtS( signal );
    double s2     = sqrt_s * sqrt_s * sqrt_s * sqrt_s;

    if (s2 == 0.0)
        ThrowError( "s^2 is zero" );

    double y = opt / s2;

    pProfile->Fill( sqrt_s, y, weight );
}


////////////////////////////////////////////////////////////////////////////////
void LoadCoefHistData( // inputs:
                       const ModelCompare::ObservableVector & observables,
                       const CStringVector & coefNames, const std::vector<double> & coefEval,
                       const ModelCompare::ModelFile & eventFile,
                       // outputs:
                       TH1DVector & eventData, std::vector<TH1DVector> & coefHists )
{
    eventData.clear();
    coefHists.clear();

    if (observables.empty() || coefNames.empty() || coefEval.empty())
        return;

    if (coefNames.size() != coefEval.size())
        ThrowError( "Size mismatch in call to LoadCoefHistData" );

    // create histograms

    for ( const ModelCompare::Observable & obs : observables )
    {
        // create event data histograms
        {
            TH1D * pHist = obs.MakeHist( eventFile.modelName, eventFile.modelTitle );
            eventData.push_back( pHist );
        }

        // create coefficient histograms
        {
            TH1DVector coefData;

            for ( const char * coefName : coefNames )
            {
                TH1D * pHist = obs.MakeHist( eventFile.modelName, eventFile.modelTitle, coefName, coefName );
                coefData.push_back(pHist);
            }

            coefHists.push_back( coefData );
        }
    }

    // fill histograms

    auto FillFunc = [&](const HepMC::GenVertex & signal)
    {
        // get coefficent factors from event weights
        std::vector<double> coefFactors;
        {
            const HepMC::GenEvent * pEvent = signal.parent_event();
            if (!pEvent)
                ThrowError( "Signal has no parent event" );

            const HepMC::WeightContainer & weights = pEvent->weights();

            coefFactors.reserve( coefNames.size() );

            for ( const std::string name : coefNames )
            {
                if (!weights.has_key(name))
                    ThrowError( "Coefficient " + name + " not found in event." );
                double w = weights[name];
                coefFactors.push_back(w);
            }
        }

        // calculate evaluation ME
        double evalME = 0.0;
        {
            for (size_t i = 0; i < coefFactors.size(); ++i)
            {
                evalME += coefFactors[i] * coefEval[i];
            }

            if (evalME <= 0)
                ThrowError( "Evaluation ME is zero or negative." );
        }

        auto itrData    = eventData.cbegin();
        auto itrCoefObs = coefHists.cbegin();
        for ( const ModelCompare::Observable & obs : observables )
        {
            obs.fillFunction( **itrData++, 1, signal );

            auto itrCoef = (*itrCoefObs++).cbegin();
            for  ( double cf : coefFactors )
            {
                double w = cf / evalME;
                obs.fillFunction( **itrCoef++, w, signal );
            }
        }
    };

    LoadEvents( eventFile.fileName, FillFunc );
}

////////////////////////////////////////////////////////////////////////////////
void CalcEvalVector( const CStringVector & coefNames, const ParamVector & params, std::vector<double> & evals )
{
    evals.clear();

    std::map<std::string, double> paramEvals;
    {
        paramEvals["F_0_0"] = 1.0;

        for (size_t j = 0; j < params.size(); ++j)
        {
            std::string name = "F_0_" + std::to_string(j+1) + "_" + params[j].name;
            paramEvals[name] = params[j].value;
        }

        for (size_t i = 0; i < params.size(); ++i)
        {
            for (size_t j = i; j < params.size(); ++j)
            {
                std::string name = "F_" + std::to_string(i+1) + "_" + std::to_string(j+1) + "_" + std::string(params[i].name);
                if (i != j) name += "_" + std::string(params[j].name);
                paramEvals[name] = params[i].value * params[j].value;
            }
        }
    }

    if (coefNames.size() != paramEvals.size())
        ThrowError( "Size mismatch between parameters and coefficient names." );

    for ( const char * name : coefNames )
    {
        auto itrEval = paramEvals.find( name );
        if (itrEval == paramEvals.end())
            ThrowError( "Coefficient " + std::string(name) + " not found." );

        evals.push_back( itrEval->second );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightHist( const ConstTH1DVector & coefData, const std::vector<double> & coefEval,
                           const char * name = nullptr, const char * title = nullptr )
{
    TH1D * pHist = (TH1D *)coefData[0]->Clone();    // polymorphic clone
    pHist->SetDirectory( nullptr);                  // ensure not owned by any directory

    pHist->Scale( coefEval[0] );

    for (size_t c = 1; c < coefEval.size(); ++c)
    {
        pHist->Add( coefData[c], coefEval[c] );
    }

    // set name and title
    {
        std::string sName  = "RW_"       + std::string(pHist->GetName());
        std::string sTitle = "Reweight " + std::string(pHist->GetTitle());

        pHist->SetName(  name  && name[0]  ? name  : sName.c_str()  );
        pHist->SetTitle( title && title[0] ? title : sTitle.c_str() );
    }

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightFactor( const ConstTH1DVector & coefData, const std::vector<double> & evalNum, const std::vector<double> & evalDen,
                             bool bClearError = false,
                             const char * name = nullptr, const char * title = nullptr )
{
    TH1D * pHistNum = CreateReweightHist( coefData, evalNum );  // can be TH1D or TProfile
    TH1D * pHistDen = CreateReweightHist( coefData, evalDen );  // can be TH1D or TProfile

    pHistNum = ConvertTProfileToTH1D( pHistNum, true );  // also delete old pHistNum
    pHistDen = ConvertTProfileToTH1D( pHistDen, true );  // also delete old pHistDen

    pHistNum->Divide( pHistDen );

    delete pHistDen;

    if (bClearError)
    {
        Int_t nSize = pHistNum->GetSize();
        for (Int_t bin = 0; bin < nSize; ++bin)  // include under/overflow bins
        {
            pHistNum->SetBinError(bin, 0);
        }

        pHistNum->ResetStats();    // force recalculation of sumw2
    }

    // set name and title
    {
        std::string sName  = "RWF_" + std::string(coefData[0]->GetName());
        std::string sTitle = "Reweight factor " + std::string(coefData[0]->GetTitle());

        pHistNum->SetName(  name  && name[0]  ? name  : sName.c_str()  );
        pHistNum->SetTitle( title && title[0] ? title : sTitle.c_str() );
    }

    return pHistNum;
}

////////////////////////////////////////////////////////////////////////////////
void ApplyReweightFactor( TH1D & target, const TH1D & reweightFactor )
{
    // Multiply each bin in target by a corresponding reweight factor.
    // Do not use the reweight factor error to error propagate.
    // Handle both TH1D and TProfile.

    // Note: For TH1D, could use target.Multiply( reweightFactor) if errors in reweightFactor are zero.
    //       However, this does not work for TProfile, where Multiply is not supported.
    // Note: The following code is based upon TProfile::Scale and TH1D::Scale

    //LogMsgInfo( "------ reweightFactor (before) ------" );
    //reweightFactor.Print("all");

    //LogMsgInfo( "------ target (before) ------" );
    //LogMsgHistStats( target );
    //target.Print("all");

    Int_t nSize = reweightFactor.GetSize();  // includes under/overflow

    if ((nSize != target.GetSize()) || (nSize != target.GetSumw2N()))
        ThrowError( "ApplyReweightFactor: bins mismatch between reweight factor and target histograms." );

    Double_t * targetContent = target.GetArray();
    Double_t * targetSumw2   = target.GetSumw2()->GetArray();

    for (Int_t bin = 0; bin < nSize; ++bin)  // include under/overflow
    {
        Double_t factor = reweightFactor.GetBinContent(bin);

        targetContent[bin] *= factor;
        targetSumw2[bin]   *= factor * factor;
    }

    target.ResetStats();

    //LogMsgInfo( "------ target (after) ------" );
    //LogMsgHistStats( target );
    //target.Print("all");
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistEffectiveEntries( const TH1D & hist )
{
    LogMsgInfo( "%hs: entries = %g, eff. entries = %g, sum bins = %g",
        FMT_HS(hist.GetName()),
        FMT_F(hist.GetEntries()), FMT_F(hist.GetEffectiveEntries()),
        FMT_F(hist.GetSumOfWeights()) );
}

////////////////////////////////////////////////////////////////////////////////
void LogMsgHistEffectiveEntries( const ConstTH1DVector & hists )
{
    for (const TH1D * pHist : hists)
    {
        if (pHist)
            LogMsgHistEffectiveEntries( *pHist );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * ReweightHist( const TH1D & sourceData, const ConstTH1DVector & sourceCoefs,
                     const std::vector<double> & sourceEval, const std::vector<double> & targetEval,
                     const char * name = nullptr, const char * title = nullptr )
{
    // create target histogram

    TH1D * pTarget = (TH1D *)sourceData.Clone();   // polymorphic clone
    pTarget->SetDirectory( nullptr );              // ensure not owned by any directory

    // set name and title
    {
        std::string sName  = "RW_"       + std::string(pTarget->GetName());
        std::string sTitle = "Reweight " + std::string(pTarget->GetTitle());

        pTarget->SetName(  name  && name[0]  ? name  : sName.c_str()  );
        pTarget->SetTitle( title && title[0] ? title : sTitle.c_str() );
    }

    // get the reweight factor, with errors zeroed
    std::unique_ptr<TH1D> upReweightFactor( CreateReweightFactor( sourceCoefs, targetEval, sourceEval, true ) );

    ApplyReweightFactor( *pTarget, *upReweightFactor );

    LogMsgHistEffectiveEntries(*upReweightFactor);
    LogMsgHistEffectiveEntries(*pTarget);

    return pTarget;
}

////////////////////////////////////////////////////////////////////////////////
void LoadReweightFiles( // inputs:
                        const ModelCompare::ObservableVector & observables,
                        const CStringVector & coefNames,
                        const ModelCompare::ModelFile & targetFile,
                        const ModelCompare::ModelFile & sourceFile, const ParamVector & sourceParam,
                        // outputs:
                        TH1DVector &                targetData,     // targetData[observable]
                        TH1DVector &                sourceData,     // sourceData[observable]
                        std::vector<TH1DVector>  &  sourceCoefs,    // sourceCoefs[observable][coefficient]
                        std::vector<double> &       sourceEval      // sourceEval[coefficient]
                        )
{
    // load target data

    {
        std::vector<TH1DVector> histData;
        ModelCompare::LoadHistData( { targetFile }, observables, histData );

        targetData = histData[0];

        // scale to 1 fb^1
        ModelCompare::ScaleHistToLuminosity( 1.0, targetData, targetFile );

        LogMsgHistUnderOverflow(    ToConstTH1DVector(targetData) );
        LogMsgHistEffectiveEntries( ToConstTH1DVector(targetData) );
    }

    // load source data and reweight coefficients

    CalcEvalVector( coefNames, sourceParam, sourceEval );

    {
        LoadCoefHistData( observables, coefNames, sourceEval, sourceFile,  // inputs
                          sourceData, sourceCoefs );                       // outputs

        // process sourceData
        {
            ModelCompare::ScaleHistToLuminosity( 1, sourceData, sourceFile );  // scale to 1 fb^-1

            LogMsgHistUnderOverflow(    ToConstTH1DVector(sourceData) );
            LogMsgHistEffectiveEntries( ToConstTH1DVector(sourceData) );
        }

        // process sourceCoefs

        for (const auto & coefData : sourceCoefs)
        {
            ModelCompare::ScaleHistToLuminosity( 1, coefData, sourceFile );   // scale to 1 fb^-1

            LogMsgHistEffectiveEntries( ToConstTH1DVector(coefData) );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void ReweightEFT( const char * outputFileName,
                  const ModelCompare::ObservableVector & observables,
                  const CStringVector & coefNames,
                  const ModelCompare::ModelFile & targetFile, const ParamVector & targetParam,
                  const ModelCompare::ModelFile & sourceFile, const ParamVector & sourceParam )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);
    // enable automatic sumw2 for every histogram
    TH1::SetDefaultSumw2(kTRUE);

    // ------

    LogMsgInfo( "Output file: %hs", FMT_HS(outputFileName) );
    std::unique_ptr<TFile> upOutputFile( new TFile( outputFileName, "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(outputFileName) );
        ThrowError( std::invalid_argument( outputFileName ) );
    }

    // load input histograms

    TH1DVector              targetData;     // targetData[observable]
    TH1DVector              sourceData;     // sourceData[observable]
    std::vector<TH1DVector> sourceCoefs;    // sourceCoefs[observable][coefficient]
    std::vector<double>     sourceEval;     // sourceEval[coefficient]

    LoadReweightFiles( observables, coefNames, targetFile, sourceFile, sourceParam, // inputs
                       targetData, sourceData, sourceCoefs, sourceEval );           // outputs

    // write input hists

    WriteHists( upOutputFile.get(), targetData );   // output file takes ownership of histograms
    WriteHists( upOutputFile.get(), sourceData );   // output file takes ownership of histograms

    // uncomment to save coefficent hists to file
    //  for (const auto & coefData : sourceCoefs)
    //      WriteHists( upOutputFile.get(), coefData );  // output file takes ownership of histograms

    // reweight source to target

    std::vector<double> targetEval;
    CalcEvalVector( coefNames, targetParam, targetEval );

    for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
    {
        LogMsgInfo( "\n****** Reweigting %hs to %hs for %hs ******\n",
                    FMT_HS(sourceFile.modelTitle), FMT_HS(targetFile.modelTitle),
                    FMT_HS(observables[obsIndex].title) );

        TH1D * pReweightSource = ReweightHist( *sourceData[obsIndex], ToConstTH1DVector( sourceCoefs[obsIndex] ),
                                               sourceEval, targetEval );

        WriteHists( upOutputFile.get(), { pReweightSource } );  // output file takes ownership of histograms

        /*
        if (obsIndex == 0)
        {
            const TH1D & target   = *targetData[obsIndex];
            const TH1D & reweight = *pReweightSource;

            LogMsgInfo( "--------- target -----------" );
            for (Int_t i = 0; i <= target.GetNbinsX() + 1; ++i)
            {
                LogMsgInfo( "Bin %3i: %.15E ± %.15E",
                            FMT_I(i), FMT_F(target.GetBinContent(i)), FMT_F(target.GetBinError(i)) );
            }

            LogMsgInfo( "--------- reweight -----------" );
            for (Int_t i = 0; i <= reweight.GetNbinsX() + 1; ++i)
            {
                LogMsgInfo( "Bin %3i: %.15E ± %.15E (%.4G)",
                            FMT_I(i), FMT_F(reweight.GetBinContent(i)), FMT_F(reweight.GetBinError(i)),
                            FMT_F(reweight.GetBinError(i) / target.GetBinError(i)) );
            }
        }
        */

        std::string reweightName  = "RW_"       + std::string(sourceFile.modelName);
        std::string reweightTitle = "Reweight " + std::string(sourceFile.modelTitle);

        ModelCompare::ModelFileVector   models  = { targetFile,           sourceFile      };
        ConstTH1DVector                 obsData = { targetData[obsIndex], pReweightSource };
        TH1DVector                      obsComp;

        models[1].modelName  = reweightName.c_str();
        models[1].modelTitle = reweightTitle.c_str();

        ModelCompare::CalculateCompareHists( observables[obsIndex], obsData, obsComp, models, { kBlue, kRed } );

        WriteHists( upOutputFile.get(), obsComp );  // output file takes ownership of histograms

        {
            std::string figName  = "fig_" + std::string(obsComp[0]->GetName());
            std::string figTitle = obsComp[0]->GetTitle();

            ModelCompare::WriteCompareFigure( figName.c_str(), figTitle.c_str(), obsData, ToConstTH1DVector(obsComp), { kBlue, kRed } );
        }
    }

    upOutputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////

}
