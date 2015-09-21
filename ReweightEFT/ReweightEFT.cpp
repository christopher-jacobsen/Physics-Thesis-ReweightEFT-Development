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
void LoadCoefHistData( const ModelCompare::ModelFile & eventFile, const ModelCompare::ObservableVector & observables,
                       const CStringVector & coefNames, const std::vector<double> & coefEval,
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
void ScaleToLuminosity( double luminosity, const TH1DVector & hists, const ModelCompare::ModelFile & eventFile, bool bApplyCrossSectionError = false )
{
    double scale = luminosity * eventFile.crossSection * 1000 / eventFile.crossSectionEvents;

    for (TH1D * pHist : hists)
    {
        LogMsgInfo( "Scaling %hs with %g", FMT_HS(pHist->GetName()), FMT_F(scale) );

        pHist->Scale( scale );

        if (bApplyCrossSectionError)
        {
            double relError = eventFile.crossSectionError / eventFile.crossSection;

            for (Int_t bin = 0; bin <= pHist->GetNbinsX() + 1; ++bin)
            {
                Double_t binContent = pHist->GetBinContent(bin);
                Double_t addError   = binContent * relError;

                Double_t binError   = pHist->GetBinError(bin);
                Double_t newError   = std::sqrt( binError * binError + addError * addError );

                pHist->SetBinError( bin, newError );
            }

            pHist->ResetStats();  // force recalculation of sumw2
        }
    }
}

#if 0
////////////////////////////////////////////////////////////////////////////////
void TrimHistNearZero( TH1D & hist )
{
    //const Double_t epsilon = std::numeric_limits<Double_t>::epsilon();

    bool bTrimmed = false;

    Int_t nBins = hist.GetNbinsX();
    for (Int_t bin = 0; bin <= nBins + 1; ++bin)  // include under/overflow bins
    {
        Double_t value = std::abs( hist.GetBinContent(bin) );  // handle both positive & negative
        if ((value > 0) && (value < 1))
        {
            Double_t error       = hist.GetBinError(bin);
            Double_t upperSigma5 = value + 5 * error;
            Double_t lowerSigma5 = value - 5 * error;

            // if (value + 5 sigma < 1) && (value - 5sigma < 0))
            if ((upperSigma5 < 1) && (lowerSigma5 < 0))
            {
                // extremely unlikely this is a 1 and highly likely this is a zero

                LogMsgInfo( "%hs setting bin %i: %g (±%g) to 0", FMT_HS(hist.GetName()), FMT_I(bin), FMT_F(value), FMT_F(error) );

                hist.SetBinContent( bin, 0 );
                hist.SetBinError(   bin, 0 );
                bTrimmed = true;
            }
        }
    }

    if (bTrimmed)
        hist.ResetStats();
}

////////////////////////////////////////////////////////////////////////////////
void TrimHistNearZero( const TH1DVector & hists )
{
    for (TH1D * pHist : hists)
    {
        if (pHist)
            TrimHistNearZero( *pHist );
    }
}
#endif

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightHist( const ConstTH1DVector & coefData, const std::vector<double> & coefEval )
{
    TH1D * pHist = (TH1D *)coefData[0]->Clone();    // polymorphic clone
    pHist->SetDirectory( nullptr);                  // ensure not owned by any directory

    pHist->Scale( coefEval[0] );

    for (size_t c = 1; c < coefEval.size(); ++c)
    {
        pHist->Add( coefData[c], coefEval[c] );
    }

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
TH1D * ConvertTProfileToTH1D( TH1D * pProfile, bool bDeleteProfile )
{
    if (!pProfile->InheritsFrom(TProfile::Class()))
        return pProfile;

    TH1D * pHist = ((TProfile *)pProfile)->ProjectionX();

    if (bDeleteProfile)
        delete pProfile;

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightFactor( const ConstTH1DVector & coefData, const std::vector<double> & evalNum, const std::vector<double> & evalDen, bool bClearError = false )
{
    TH1D * pHistNum = CreateReweightHist( coefData, evalNum );  // can be TH1D or TProfile
    TH1D * pHistDen = CreateReweightHist( coefData, evalDen );  // can be TH1D or TProfile

    pHistNum = ConvertTProfileToTH1D( pHistNum, true );
    pHistDen = ConvertTProfileToTH1D( pHistDen, true );

    pHistNum->Divide( pHistDen );

    delete pHistDen;

    if (bClearError)
    {
        Int_t nBins = pHistNum->GetNbinsX();
        for (Int_t bin = 0; bin <= nBins + 1; ++bin)
        {
            pHistNum->SetBinError(bin, 0);
        }

        pHistNum->ResetStats();    // force recalculation of sumw2
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
void ReweightEFT( const char * outputFileName, const ModelCompare::ObservableVector & observables,
                  const ModelCompare::ModelFile & sourceFile, const ModelCompare::ModelFile & targetFile,
                  const ParamVector &             sourceParam, const ParamVector &            targetParam,
                  const CStringVector & coefNames )
{
    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);

    LogMsgInfo( "Output file: %hs", FMT_HS(outputFileName) );
    std::unique_ptr<TFile> upOutputFile( new TFile( outputFileName, "RECREATE" ) );
    if (upOutputFile->IsZombie() || !upOutputFile->IsOpen())    // IsZombie is true if constructor failed
    {
        LogMsgError( "Failed to create output file (%hs).", FMT_HS(outputFileName) );
        ThrowError( std::invalid_argument( outputFileName ) );
    }

    // load target data

    TH1DVector targetData; // targetData[observable]
    {
        std::vector<TH1DVector> histData;
        ModelCompare::LoadHistData( { targetFile }, observables, histData );

        targetData = histData[0];

        // scale to 1 fb^1
        ScaleToLuminosity( 1.0, targetData, targetFile );

        LogMsgHistUnderOverflow( ToConstTH1DVector(targetData) );
        LogMsgHistEffectiveEntries( ToConstTH1DVector(targetData) );

        WriteHists( upOutputFile.get(), targetData );  // output file takes ownership of histograms
    }

    // load source data and reweight coefficients

    std::vector<double> sourceEval;
    CalcEvalVector( coefNames, sourceParam, sourceEval );

    TH1DVector              sourceData;         // sourceData[observable]
    std::vector<TH1DVector> reweightCoefData;   // reweightCoefData[observable][coefficient]
    {
        LoadCoefHistData( sourceFile, observables, coefNames, sourceEval, sourceData, reweightCoefData );

        // process sourceData
        {
            ScaleToLuminosity( 1, sourceData, sourceFile );  // scale to 1 fb^-1

            LogMsgHistUnderOverflow( ToConstTH1DVector(sourceData) );
            LogMsgHistEffectiveEntries( ToConstTH1DVector(sourceData) );

            WriteHists( upOutputFile.get(), sourceData );
        }

        // process reweightCoefData

        for (const auto & coefData : reweightCoefData)
        {
            ScaleToLuminosity( 1, coefData, sourceFile );   // scale to 1 fb^-1

            LogMsgHistEffectiveEntries( ToConstTH1DVector(coefData) );

            // uncomment to save coefficent hists to file
            // WriteHists( upOutputFile.get(), coefData );
        }
    }

    /*
    // reweight source to source

    TH1DVector reweightSource;

    for ( const TH1DVector & coefData : reweightCoefData )
    {
        TH1D * pHist = CreateReweightHist( ToConstTH1DVector(coefData), sourceEval );
        reweightSource.push_back( pHist );
    }
    */

    // reweight source to target

    std::vector<double> targetEval;
    CalcEvalVector( coefNames, targetParam, targetEval );

    for (size_t obsIndex = 0; obsIndex < observables.size(); ++obsIndex)
    {
        const ConstTH1DVector & coefData = ToConstTH1DVector( reweightCoefData[obsIndex] );

        // get the reweight factor, with errors zeroed
        std::unique_ptr<TH1D> upReweightFactor( CreateReweightFactor( coefData, targetEval,  sourceEval, true ) );

        std::unique_ptr<TH1D> upReweightTarget( (TH1D *)sourceData[obsIndex]->Clone() );    // polymorphic clone
        upReweightTarget->SetDirectory( nullptr );                                          // ensure not owned by any directory

        ApplyReweightFactor( *upReweightTarget, *upReweightFactor );

        //TH1D * pOldReweight = CreateReweightHist( coefData, targetEval );

        LogMsgInfo( "reweightFactor:" );
        LogMsgHistEffectiveEntries(*upReweightFactor);

        LogMsgInfo( "reweightTarget:" );
        LogMsgHistEffectiveEntries(*upReweightTarget);

        /*
        if (obsIndex == 0)
        {
            TH1D *  pHist   = nullptr;
            Int_t   bin     = 221;

            pHist = sourceData[obsIndex];
            LogMsgInfo( "%hs bin %i: %g (±%g)", FMT_HS(pHist->GetName()), FMT_I(bin), FMT_F(pHist->GetBinContent(bin)), FMT_F(pHist->GetBinError(bin)) );

            pHist = upReweightFactor.get();
            LogMsgInfo( "%hs bin %i: %g (±%g)", FMT_HS(pHist->GetName()), FMT_I(bin), FMT_F(pHist->GetBinContent(bin)), FMT_F(pHist->GetBinError(bin)) );
        }
        */

        // TODO: write reweight target (needs name and title)

        //LogMsgInfo( "oldReweightTarget:" );
        //LogMsgHistEffectiveEntries(*pOldReweight);

        ConstTH1DVector                 obsData = { targetData[obsIndex], upReweightTarget.get() };
        TH1DVector                      obsComp;
        ModelCompare::ModelFileVector   models = { targetFile, sourceFile };

        /*
        if (obsIndex == 0)
        {
            const TH1D & target   = *obsData[0];
            const TH1D & reweight = *obsData[1];

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
