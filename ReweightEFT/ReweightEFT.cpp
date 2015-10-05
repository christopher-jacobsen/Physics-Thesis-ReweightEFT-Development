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
#include <TStyle.h>
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

// gain access to TProfile fBinEntries
struct MyProfile : public TProfile
{
    using TProfile::fBinEntries;
};

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

    /*
    if (pProfile->FindFixBin(sqrt_s) == 88)
    {
        //if (strcmp(pProfile->GetName(), "EFT_all_cWWW_O1vS_F_0_0") == 0)
        {
            LogMsgInfo( "%hs: w=%g y=%g", FMT_HS(pProfile->GetName()), FMT_F(weight), FMT_F(opt) );
        }
    }
    */

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
                pHist->Sumw2(kFALSE);  // do not use sumw2 for coefficient histograms
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
void HistReweightAdd( TH1D & hist, const TH1D & other, double cf )
{
    // This method should be used instead of TH1D::Add or TProfile::Add for adding
    // reweight coefficient histograms.

    // Neither TH1D::Add nor TProfile::Add change the number of effective entries,
    // when scaling other prior to adding the result to hist.
    // They essentially scale the weights within the sums of other.
    //
    // TH1D::Add:       sumw  *= cf             -> content *=  cf
    //                  sumw2 *= cf*cf          -> error   *= |cf|
    //                                          -> nEff    *=  1    [nEff = sumw^2/sumw2]
    //
    // TProfile::Add:   sumw     *=  cf         -> mean  *= sign(cf)
    //                  sumw2    *= |cf|        -> sigma *= 1
    //                  binEnt   *= |cf]        -> error *= 1
    //                  binSumw2 *= |cf|*|cf|   -> nEff  *= 1       [nEff = binEnt^2/binSumw2]
    //
    // For example, when adding a hist to itself with cf = -1, we get:
    //
    // TH1D::Add:       sumw   = 0          -> content  = 0
    //                  sumw2 *= 2          -> error   *= sqrt(2) != 0
    //                                      -> nEff     = 0
    //
    // TProfile::Add:   sumw      = 0       -> mean   = 0
    //                  sumw2    *= 2       -> sigma  = sqrt(sumw2/binEnt) != 0
    //                  binEnt   *= 2       -> error  = sigma/sqrt(nEff)   != 0
    //                  binSumw2 *= 2       -> nEff  *= 2
    //
    // The above is mathematically correct, if one wants to propagate errors.
    // However, we want to construct a final reweighted histogram from the
    // weighted sum of several other histograms, in such a way that the
    // reweighted result is equivalent to the result that would have been
    // achieved if filling the histogram in one pass.
    //
    // That is for when summing over event index k:
    //      sum(w0(k))*c0 + sum(w1(k))*c1 + ... = sum(w0(k)*c0 + w1(k)*c1 + ...)
    //
    // So we need to implement:
    // TH1D:        sumw  *= cf             -> content *= cf
    //              sumw2  = |sumw|         -> error    = sqrt(nEff)
    //                                      -> nEff     = |sumw|
    //
    // TProfile:    sumw     *= cf          -> mean  *= 1
    //              sumw2    *= cf          -> sigma *= 1
    //              binEnt   *= cf          -> error  = sigma/sqrt(nEff)
    //              binSumw2  = |binEnt|    -> nEff   = |binEnt|
    //
    // Note: that anyplace a weight is not used directly (squared or absolute)
    // that it is not possible to accomplish our goal. For example:
    //      sum(w0(k)^2)*c0 + sum(w1(k)^2)*c1 + ... != sum((w0(k)*c0 + w1(k)*c1 + ...)^2)
    //      sum(|w0(k)|)*c0 + sum(|w1(k)|)*c1 + ... = sum(|w0(k)*c0 + w1(k)*c1 + ...|)
    // Thus we cannot correct sumw2 for TH1D or binSumw2 for TProfile.
    // Fortunately we can assume that the resulting effective events is an average count,
    // and the error is then sqrt(nEff) for TH1D and sigma/sqrt(nEff) for TProfile.
    //
    // Using this technique, adding a hist to itself with cf = -1 we get:
    //
    // TH1D:        sumw  = 0           -> content = 0
    //              sumw2 = 0           -> error   = 0
    //                                  -> nEff    = 0
    //
    // TProfile:    sumw     = 0        -> mean  = 0
    //              sumw2    = 0        -> sigma = 0
    //              binEnt   = 0        -> error = 0
    //              binSumw2 = 0        -> nEff  = 0
    //
    // And this is what we want.

    if (hist.GetSize() != other.GetSize())
        ThrowError( "ReweightHistAdd: size mismatch in histograms" );

    if ((hist.GetSize() == 0) || (cf == 0))
        return;

          MyProfile * hp = hist .InheritsFrom(TProfile::Class()) ? static_cast<      MyProfile *>(&hist)  : nullptr;
    const MyProfile * op = other.InheritsFrom(TProfile::Class()) ? static_cast<const MyProfile *>(&other) : nullptr;

    if ((hp == nullptr) != (op == nullptr))
        ThrowError( "Cannot add TProfile to TH1D. Types must match" );

    Double_t * hW  = hist.GetArray();
    Double_t * hW2 = hist.GetSumw2()->GetArray();                       // TH1D: can be null
    Double_t * hB  = hp ? hp->fBinEntries.GetArray()    : nullptr;
    Double_t * hB2 = hp ? hp->GetBinSumw2()->GetArray() : nullptr;      // TProfile: can be null

    const Double_t * oW  = other.GetArray();
    const Double_t * oW2 = other.GetSumw2()->GetArray();                // TH1D: can be null
    const Double_t * oB  = op ? op->fBinEntries.GetArray() : nullptr;

    if (!hW || !oW || (hp && !(hW2 && oW2 && hB && oB)))
        ThrowError( "ReweightHistAdd: Internal error - unexpected buffer state." );

    const Int_t nSize = hist.GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)
    {
        hW[bin] += oW[bin] * cf;

        if (!hp)    // TH1D
        {
            if (hW2)
                hW2[bin] = std::abs( hW[bin] );     // -> error = sqrt(nEff)
        }
        else    // TProfile
        {
            hW2[bin] += oW2[bin] * cf;
             hB[bin] +=  oB[bin] * cf;

            if (hB2)
                hB2[bin] = std::abs( hB[bin] );     // -> error = sigma/sqrt(nEff)
        }
    }

    hist.ResetStats();
}

////////////////////////////////////////////////////////////////////////////////
void ValidateReweightHist( const TH1D & hist )
{
    // validate constructed reweight histogram
    // TH1D:     non-negative sumw and sumw2
    // TProfile: non-negative binEntries and binSumw2

    const Double_t * sumw  = nullptr;
    const Double_t * sumw2 = nullptr;

    if (!hist.InheritsFrom(TProfile::Class()))
    {
        sumw  = hist.GetArray();
        sumw2 = hist.GetSumw2()->GetArray();
    }
    else
    {
        sumw  = static_cast<const MyProfile &>(hist).fBinEntries.GetArray();
        sumw2 = static_cast<const MyProfile &>(hist).GetBinSumw2()->GetArray();
    }

    const Int_t nSize = hist.GetSize();
    for (Int_t bin = 0; bin < nSize; ++bin)
    {
        if (sumw[bin] < 0)
            ThrowError( "%s: negative sumw of %g in bin %i", FMT_HS(hist.GetName()), FMT_F(sumw[bin]), FMT_I(bin) );

        if (sumw2[bin] < 0)
            ThrowError( "%s: negative sumw2 of %g in bin %i", FMT_HS(hist.GetName()), FMT_F(sumw[bin]), FMT_I(bin) );
    }
}

////////////////////////////////////////////////////////////////////////////////
TH1D * CreateReweightHist( const ConstTH1DVector & coefData, const std::vector<double> & coefEval,
                           const char * name = nullptr, const char * title = nullptr )
{
    TH1D * pHist = (TH1D *)coefData[0]->Clone();    // polymorphic clone
    pHist->SetDirectory( nullptr);                  // ensure not owned by any directory

    // set name and title
    {
        std::string sName  = "RW_"       + std::string(pHist->GetName());
        std::string sTitle = "Reweight " + std::string(pHist->GetTitle());

        pHist->SetName(  name  && name[0]  ? name  : sName.c_str()  );
        pHist->SetTitle( title && title[0] ? title : sTitle.c_str() );
    }

    pHist->Reset("M");  // clear everything, including minimum/maximum

    // We desire sumw2 to be enabled for the reweighted result.
    // coefData does not normally have Sumw2 set, so neither does its clone.
    // We can enable sumw2 before or after the reweight calculation.
    // HistReweightAdd (used below) ensures the "correct" result when sumw is enabled.
    // Doing it before the calculation, saves on a ResetStats call,
    // which would be necessary if enabled afterwards.
    pHist->Sumw2(kTRUE);

    for (size_t c = 0; c < coefEval.size(); ++c)
    {
        HistReweightAdd( *pHist, *coefData[c], coefEval[c] );   // we cannot use Add, as it gives the incorrect result
    }

    // validate reweight histogram
    ValidateReweightHist(*pHist);

    //LogMsgHistEffectiveEntries(*pHist);

    return pHist;
}

////////////////////////////////////////////////////////////////////////////////
TH1D * ReweightHist( const TH1D & sourceData, const ConstTH1DVector & sourceCoefs,
                     const std::vector<double> & sourceEval, const std::vector<double> & targetEval,
                     const char * name /*= nullptr*/, const char * title /*= nullptr*/ )
{
    // set name and title
    std::string sName  = "RW_"       + std::string(sourceData.GetName());
    std::string sTitle = "Reweight " + std::string(sourceData.GetTitle());

    if (!name  || !name[0])  name  = sName.c_str();
    if (!title || !title[0]) title = sTitle.c_str();

    TH1D * pTarget = CreateReweightHist( sourceCoefs, targetEval, name, title );

    /*
    LogMsgInfo("--- RW ---");
    LogMsgHistEffectiveEntries(sourceData);
    LogMsgHistDump(sourceData);
    LogMsgInfo("----------");
    LogMsgHistEffectiveEntries(*pTarget);
    LogMsgHistDump(*pTarget);
    LogMsgInfo("----------");
    */

    //LogMsgHistDump(*pTarget);
    //LogMsgHistDump(*sourceCoefs[0]);

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
                        std::vector<double> &       sourceEval,     // sourceEval[coefficient]
                        TH1DVector &                rawTargetData,  // rawTargetData[observable]
                        TH1DVector &                rawSourceData   // rawSourceData[observable]
                        )
{
    // load target data

    {
        std::vector<TH1DVector> histData;
        ModelCompare::LoadHistData( { targetFile }, observables, histData );

        targetData = histData[0];

        // copy targetData to rawTargetData before modifications
        for (const TH1D * pHist : targetData)
        {
            TH1D * pRawHist = (TH1D *)pHist->Clone();
            pRawHist->SetDirectory( nullptr );
            pRawHist->SetName( ("raw_" + std::string(pRawHist->GetName())).c_str() );
            rawTargetData.push_back( pRawHist );
        }

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

        // copy sourceData to rawSourceData before modifications
        for (const TH1D * pHist : sourceData)
        {
            TH1D * pRawHist = (TH1D *)pHist->Clone();
            pRawHist->SetDirectory( nullptr );
            pRawHist->SetName( ("raw_" + std::string(pRawHist->GetName())).c_str() );
            rawSourceData.push_back( pRawHist );
        }

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
    LogMsgInfo( "Style: %hs", FMT_HS(gStyle->GetName()) );

    // disable automatic histogram addition to current directory
    TH1::AddDirectory(kFALSE);
    // enable automatic sumw2 for every histogram
    TH1::SetDefaultSumw2(kTRUE);

    // modify the global style
    gStyle->SetPaperSize( TStyle::kA4 );
    gStyle->SetTitleOffset( 1.3, "xyz" ); // increase title offsets a little more
    gStyle->SetPadTopMargin(   0.03 );
    gStyle->SetPadRightMargin( 0.03 );
    gStyle->SetOptTitle( kFALSE );

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
    TH1DVector              rawTargetData;  // rawTargetData[observable]
    TH1DVector              rawSourceData;  // rawSourceData[observable]

    LoadReweightFiles( observables, coefNames, targetFile, sourceFile, sourceParam, // inputs
                       targetData, sourceData, sourceCoefs, sourceEval,             // outputs
                       rawTargetData, rawSourceData );

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

        ModelCompare::ModelFileVector   models  = { targetFile,              sourceFile              };
        ConstTH1DVector                 obsData = { targetData[obsIndex],    pReweightSource         };
        ConstTH1DVector                 rawData = { rawTargetData[obsIndex], rawSourceData[obsIndex] };
        TH1DVector                      obsComp;

        models[1].modelName  = reweightName.c_str();
        models[1].modelTitle = reweightTitle.c_str();

        ModelCompare::CalculateCompareHists( observables[obsIndex], obsData, obsComp, models, { kBlue, kRed } );

        WriteHists( upOutputFile.get(), obsComp );  // output file takes ownership of histograms

        {
            std::string figName  = "fig_" + std::string(obsComp[0]->GetName());
            std::string figTitle = obsComp[0]->GetTitle();

            ModelCompare::WriteCompareFigure( figName.c_str(), figTitle.c_str(), obsData, ToConstTH1DVector(obsComp), { kBlue, kRed }, rawData );
        }
    }

    upOutputFile->Close();
}

////////////////////////////////////////////////////////////////////////////////

}
