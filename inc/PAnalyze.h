#ifndef _PAnalyze_h__
#define _PAnalyze_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TF1.h"
#include "TRandom3.h"

class	PAnalyze  : public PPhysics
{
private:
    TH1*    TaggerAccScal;
    TH1*    LiveTimeScal;

    TH1*    CorrTaggScal;
    TH1*    PolarizeScal;
    TH1*    Helicity;

    TH1*    Tagg_Tm;
    TH1*    Tagg_0;
    TH1*    Tagg_0_R;
    TH1*    Tagg_1;
    TH1*    Tagg_1_R;

    TH1*    Inc_Tm;
    TH2*    Inc_0;
    TH2*    Inc_0_R;
    TH2*    Inc_1;
    TH2*    Inc_1_R;

    TH1*    Pi0_IM_A;
    TH1*    Pi0_IM_E;
    TH1*    Pi0_IM_I;
    TH3*    Pi0_IM_CC;
    TH3*    Pi0_IM_CT;

    TH2*    Pi0_Sp;
    TH1*    Pi0_CA;

    TH2*    Pi0_Tm_NE;
    TH2*    Pi0_Tm_NI;
    TH2*    Pi0_Tm_CE;
    TH2*    Pi0_Tm_CI;
    TH2*    Pi0_Tm_WE;
    TH2*    Pi0_Tm_WI;
    TH2*    Pi0_Tm_TE;
    TH2*    Pi0_Tm_TI;

    TH3*    Pi0_OA;
    TH3*    Pi0_OA_R;
    TH3*    Pi0_OA_Cut;
    TH3*    Pi0_OA_Cut_R;

    TH3*    Pi0_MM_NE_0;
    TH3*    Pi0_MM_NE_0_R;
    TH3*    Pi0_MM_NE_1;
    TH3*    Pi0_MM_NE_1_R;
    TH3*    Pi0_Ph_NE_0;
    TH3*    Pi0_Ph_NE_0_R;
    TH3*    Pi0_Ph_NE_1;
    TH3*    Pi0_Ph_NE_1_R;

    TH3*    Pi0_MM_NI_0;
    TH3*    Pi0_MM_NI_0_R;
    TH3*    Pi0_MM_NI_1;
    TH3*    Pi0_MM_NI_1_R;
    TH3*    Pi0_Ph_NI_0;
    TH3*    Pi0_Ph_NI_0_R;
    TH3*    Pi0_Ph_NI_1;
    TH3*    Pi0_Ph_NI_1_R;

    TH3*    Pi0_MM_CE_0;
    TH3*    Pi0_MM_CE_0_R;
    TH3*    Pi0_MM_CE_1;
    TH3*    Pi0_MM_CE_1_R;
    TH3*    Pi0_Ph_CE_0;
    TH3*    Pi0_Ph_CE_0_R;
    TH3*    Pi0_Ph_CE_1;
    TH3*    Pi0_Ph_CE_1_R;

    TH3*    Pi0_MM_CI_0;
    TH3*    Pi0_MM_CI_0_R;
    TH3*    Pi0_MM_CI_1;
    TH3*    Pi0_MM_CI_1_R;
    TH3*    Pi0_Ph_CI_0;
    TH3*    Pi0_Ph_CI_0_R;
    TH3*    Pi0_Ph_CI_1;
    TH3*    Pi0_Ph_CI_1_R;

    TH3*    Pi0_MM_WE_0;
    TH3*    Pi0_MM_WE_0_R;
    TH3*    Pi0_MM_WE_1;
    TH3*    Pi0_MM_WE_1_R;
    TH3*    Pi0_Ph_WE_0;
    TH3*    Pi0_Ph_WE_0_R;
    TH3*    Pi0_Ph_WE_1;
    TH3*    Pi0_Ph_WE_1_R;

    TH3*    Pi0_MM_WI_0;
    TH3*    Pi0_MM_WI_0_R;
    TH3*    Pi0_MM_WI_1;
    TH3*    Pi0_MM_WI_1_R;
    TH3*    Pi0_Ph_WI_0;
    TH3*    Pi0_Ph_WI_0_R;
    TH3*    Pi0_Ph_WI_1;
    TH3*    Pi0_Ph_WI_1_R;

    TH3*    Pi0_MM_TE_0;
    TH3*    Pi0_MM_TE_0_R;
    TH3*    Pi0_MM_TE_1;
    TH3*    Pi0_MM_TE_1_R;
    TH3*    Pi0_Ph_TE_0;
    TH3*    Pi0_Ph_TE_0_R;
    TH3*    Pi0_Ph_TE_1;
    TH3*    Pi0_Ph_TE_1_R;

    TH3*    Pi0_MM_TI_0;
    TH3*    Pi0_MM_TI_0_R;
    TH3*    Pi0_MM_TI_1;
    TH3*    Pi0_MM_TI_1_R;
    TH3*    Pi0_Ph_TI_0;
    TH3*    Pi0_Ph_TI_0_R;
    TH3*    Pi0_Ph_TI_1;
    TH3*    Pi0_Ph_TI_1_R;

    TH3*    Pi0_Re_All;
    TH3*    Pi0_Re_All_R;
    TH3*    Pi0_Re_Det;
    TH3*    Pi0_Re_Det_R;
    TH3*    Pi0_Re_Dif;
    TH3*    Pi0_Re_Dif_R;
    TH3*    Pi0_Re_NoE;
    TH3*    Pi0_Re_NoE_R;

    TH2*    PiP_Sp;
    TH1*    PiP_CA;

    TH2*    PiP_Tm_NE;
    TH2*    PiP_Tm_NI;
    TH2*    PiP_Tm_CE;
    TH2*    PiP_Tm_CI;
    TH2*    PiP_Tm_WE;
    TH2*    PiP_Tm_WI;
    TH2*    PiP_Tm_TE;
    TH2*    PiP_Tm_TI;

    TH3*    PiP_OA;
    TH3*    PiP_OA_R;
    TH3*    PiP_OA_Cut;
    TH3*    PiP_OA_Cut_R;

    TH3*    PiP_MM_NE_0;
    TH3*    PiP_MM_NE_0_R;
    TH3*    PiP_MM_NE_1;
    TH3*    PiP_MM_NE_1_R;
    TH3*    PiP_Ph_NE_0;
    TH3*    PiP_Ph_NE_0_R;
    TH3*    PiP_Ph_NE_1;
    TH3*    PiP_Ph_NE_1_R;

    TH3*    PiP_MM_NI_0;
    TH3*    PiP_MM_NI_0_R;
    TH3*    PiP_MM_NI_1;
    TH3*    PiP_MM_NI_1_R;
    TH3*    PiP_Ph_NI_0;
    TH3*    PiP_Ph_NI_0_R;
    TH3*    PiP_Ph_NI_1;
    TH3*    PiP_Ph_NI_1_R;

    TH3*    PiP_MM_CE_0;
    TH3*    PiP_MM_CE_0_R;
    TH3*    PiP_MM_CE_1;
    TH3*    PiP_MM_CE_1_R;
    TH3*    PiP_Ph_CE_0;
    TH3*    PiP_Ph_CE_0_R;
    TH3*    PiP_Ph_CE_1;
    TH3*    PiP_Ph_CE_1_R;

    TH3*    PiP_MM_CI_0;
    TH3*    PiP_MM_CI_0_R;
    TH3*    PiP_MM_CI_1;
    TH3*    PiP_MM_CI_1_R;
    TH3*    PiP_Ph_CI_0;
    TH3*    PiP_Ph_CI_0_R;
    TH3*    PiP_Ph_CI_1;
    TH3*    PiP_Ph_CI_1_R;

    TH3*    PiP_MM_WE_0;
    TH3*    PiP_MM_WE_0_R;
    TH3*    PiP_MM_WE_1;
    TH3*    PiP_MM_WE_1_R;
    TH3*    PiP_Ph_WE_0;
    TH3*    PiP_Ph_WE_0_R;
    TH3*    PiP_Ph_WE_1;
    TH3*    PiP_Ph_WE_1_R;

    TH3*    PiP_MM_WI_0;
    TH3*    PiP_MM_WI_0_R;
    TH3*    PiP_MM_WI_1;
    TH3*    PiP_MM_WI_1_R;
    TH3*    PiP_Ph_WI_0;
    TH3*    PiP_Ph_WI_0_R;
    TH3*    PiP_Ph_WI_1;
    TH3*    PiP_Ph_WI_1_R;

    TH3*    PiP_MM_TE_0;
    TH3*    PiP_MM_TE_0_R;
    TH3*    PiP_MM_TE_1;
    TH3*    PiP_MM_TE_1_R;
    TH3*    PiP_Ph_TE_0;
    TH3*    PiP_Ph_TE_0_R;
    TH3*    PiP_Ph_TE_1;
    TH3*    PiP_Ph_TE_1_R;

    TH3*    PiP_MM_TI_0;
    TH3*    PiP_MM_TI_0_R;
    TH3*    PiP_MM_TI_1;
    TH3*    PiP_MM_TI_1_R;
    TH3*    PiP_Ph_TI_0;
    TH3*    PiP_Ph_TI_0_R;
    TH3*    PiP_Ph_TI_1;
    TH3*    PiP_Ph_TI_1_R;

    TH2*    Comp_Sp;
    TH1*    Comp_CA;

    TH2*    Comp_Tm_NE;
    TH2*    Comp_Tm_NI;
    TH2*    Comp_Tm_CE;
    TH2*    Comp_Tm_CI;
    TH2*    Comp_Tm_WE;
    TH2*    Comp_Tm_WI;
    TH2*    Comp_Tm_TE;
    TH2*    Comp_Tm_TI;

    TH3*    Comp_OA;
    TH3*    Comp_OA_R;
    TH3*    Comp_OA_Cut;
    TH3*    Comp_OA_Cut_R;

    TH3*    Comp_MM_NE_0;
    TH3*    Comp_MM_NE_0_R;
    TH3*    Comp_MM_NE_1;
    TH3*    Comp_MM_NE_1_R;
    TH3*    Comp_Ph_NE_0;
    TH3*    Comp_Ph_NE_0_R;
    TH3*    Comp_Ph_NE_1;
    TH3*    Comp_Ph_NE_1_R;

    TH3*    Comp_MM_NI_0;
    TH3*    Comp_MM_NI_0_R;
    TH3*    Comp_MM_NI_1;
    TH3*    Comp_MM_NI_1_R;
    TH3*    Comp_Ph_NI_0;
    TH3*    Comp_Ph_NI_0_R;
    TH3*    Comp_Ph_NI_1;
    TH3*    Comp_Ph_NI_1_R;

    TH3*    Comp_MM_CE_0;
    TH3*    Comp_MM_CE_0_R;
    TH3*    Comp_MM_CE_1;
    TH3*    Comp_MM_CE_1_R;
    TH3*    Comp_Ph_CE_0;
    TH3*    Comp_Ph_CE_0_R;
    TH3*    Comp_Ph_CE_1;
    TH3*    Comp_Ph_CE_1_R;

    TH3*    Comp_MM_CI_0;
    TH3*    Comp_MM_CI_0_R;
    TH3*    Comp_MM_CI_1;
    TH3*    Comp_MM_CI_1_R;
    TH3*    Comp_Ph_CI_0;
    TH3*    Comp_Ph_CI_0_R;
    TH3*    Comp_Ph_CI_1;
    TH3*    Comp_Ph_CI_1_R;

    TH3*    Comp_MM_WE_0;
    TH3*    Comp_MM_WE_0_R;
    TH3*    Comp_MM_WE_1;
    TH3*    Comp_MM_WE_1_R;
    TH3*    Comp_Ph_WE_0;
    TH3*    Comp_Ph_WE_0_R;
    TH3*    Comp_Ph_WE_1;
    TH3*    Comp_Ph_WE_1_R;

    TH3*    Comp_MM_WI_0;
    TH3*    Comp_MM_WI_0_R;
    TH3*    Comp_MM_WI_1;
    TH3*    Comp_MM_WI_1_R;
    TH3*    Comp_Ph_WI_0;
    TH3*    Comp_Ph_WI_0_R;
    TH3*    Comp_Ph_WI_1;
    TH3*    Comp_Ph_WI_1_R;

    TH3*    Comp_MM_TE_0;
    TH3*    Comp_MM_TE_0_R;
    TH3*    Comp_MM_TE_1;
    TH3*    Comp_MM_TE_1_R;
    TH3*    Comp_Ph_TE_0;
    TH3*    Comp_Ph_TE_0_R;
    TH3*    Comp_Ph_TE_1;
    TH3*    Comp_Ph_TE_1_R;

    TH3*    Comp_MM_TI_0;
    TH3*    Comp_MM_TI_0_R;
    TH3*    Comp_MM_TI_1;
    TH3*    Comp_MM_TI_1_R;
    TH3*    Comp_Ph_TI_0;
    TH3*    Comp_Ph_TI_0_R;
    TH3*    Comp_Ph_TI_1;
    TH3*    Comp_Ph_TI_1_R;

    TH3*    Comp_CS;
    TH3*    Comp_CS_MM;
    TH3*    Comp_CS_MM_R;
    TH3*    Reco_CS;
    TH3*    Reco_CS_MM;
    TH3*    Reco_CS_MM_R;

    Int_t   verbosity;
    Bool_t  match_charge;

    Double_t IMCut;
    Double_t MMLoC;
    Double_t MMHiC;
    Double_t OACut;
    Double_t ESCut;

    Bool_t  save_randoms;
    Bool_t  split_search;
    Bool_t  pure_mwpc;

    Double_t taps_eff;
    std::vector<Bool_t> ignoreTrack;

    Bool_t  cir_beam;
    Bool_t  lin_beam;
    Double_t beamPol;
    std::vector<Int_t> beamPolTime;
    std::vector<Double_t> beamPolMeas;

    Double_t targPol;
    std::vector<Int_t> targPolTime;
    std::vector<Double_t> targPolMeas;

    Bool_t  firstEvent;

protected:
    virtual Bool_t Start();
    virtual void ProcessEvent();
    virtual void ProcessScalerRead();
    virtual Bool_t Write();

public:
    PAnalyze();
    virtual ~PAnalyze();
    virtual Bool_t Init();
    Bool_t InitVerbosity();
    Bool_t InitMatchCharge();
    Bool_t InitInvariantMass();
    Bool_t InitMissingMass();
    Bool_t InitOpeningAngle();
    Bool_t InitEnergySum();
    Bool_t InitSaveRandoms();
    Bool_t InitSplitSearch();
    Bool_t InitPureMWPC();
    Bool_t InitTAPSEff();
    Bool_t InitBeamPol();
    Bool_t InitTargPol();
    TLorentzVector AdjustMass(TLorentzVector lv, Double_t mass);
    Double_t CalcCircBeamPol(Double_t E_e, Double_t P_e, Double_t E_g);

};
#endif
