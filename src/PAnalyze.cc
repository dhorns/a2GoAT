#include "PAnalyze.h"

PAnalyze::PAnalyze()
{
    CorrTaggScal = new TH1D("CorrTaggScal", "Livetime corrected tagger scalers;Tagger Channel;LT*Scaler", 352, 0, 352);
    PolarizeScal = new TH1D("PolarizeScal", "Livetime and polarization corrected tagger scalers;Tagger Channel;LT*Pt*Pg*Scaler", 352, 0, 352);

    Helicity = new TH1D("Helicity", "Helicity", 3, 0, 3);

    // Inclusive histograms
    Tagg_Tm = new TH1D("Tagg_Tm", "Tagger Time;t_{#gamma} (ns)", 1400, -700, 700);
    Tagg_0 = new TH1D("Tagg_0", "Tagger Hits;Tagger Channel", 352, 0, 352);
    Tagg_0_R = new TH1D("Tagg_0_R", "Tagger Hits;Tagger Channel", 352, 0, 352);
    Tagg_1 = new TH1D("Tagg_1", "Tagger Hits;Tagger Channel", 352, 0, 352);
    Tagg_1_R = new TH1D("Tagg_1_R", "Tagger Hits;Tagger Channel", 352, 0, 352);

    Inc_Tm = new TH1D("Inc_Tm", "Tagger - CB Time;t_{Tagger}-t_{CB} (ns)", 1400, -700, 700);
    Inc_0 = new TH2D("Inc_0", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);
    Inc_0_R = new TH2D("Inc_0_R", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);
    Inc_1 = new TH2D("Inc_1", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);
    Inc_1_R = new TH2D("Inc_1_R", "Total Inclusive Hits;Tagger Channel;Number of Tracks", 352, 0, 352, 10, 0, 10);

    // Pi0 histograms
    Pi0_IM_A = new TH1D("Pi0_IM_A", "Pi0 Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0_IM_E = new TH1D("Pi0_IM_E", "Pi0 Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0_IM_I = new TH1D("Pi0_IM_I", "Pi0 Invariant Mass;m_{#gamma#gamma} (MeV)", 400, 0, 400);
    Pi0_IM_CC = new TH3D("Pi0_IM_CC", "CB Invariant Mass;E_{#gamma_1};E_{#gamma_2};m_{#gamma#gamma} (MeV)", 80, 0, 800, 80, 0, 800, 200, 0, 200);
    Pi0_IM_CT = new TH3D("Pi0_IM_CT", "CB/TAPS Invariant Mass;E_{#gamma_{CB}};E_{#gamma_{TAPS}};m_{#gamma#gamma} (MeV)", 80, 0, 800, 80, 0, 800, 200, 0, 200);

    Pi0_Sp = new TH2D("Pi0_Sp", "Split OA vs Energy Ratio;E_{split}/E_{#gamma};#theta_{OA} (deg)", 200, 0, 1, 180, 0, 180);
    Pi0_CA = new TH1D("Pi0_CA", "Pi0 Recoil Coplanarity Angle;#phi_{#pi^{0}}-#phi_{p} (deg)", 360, 0, 360);

    Pi0_Tm_NE = new TH2D("Pi0_Tm_NE", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_NI = new TH2D("Pi0_Tm_NI", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_CE = new TH2D("Pi0_Tm_CE", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_CI = new TH2D("Pi0_Tm_CI", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_WE = new TH2D("Pi0_Tm_WE", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_WI = new TH2D("Pi0_Tm_WI", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_TE = new TH2D("Pi0_Tm_TE", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);
    Pi0_Tm_TI = new TH2D("Pi0_Tm_TI", "Pi0 Time;Tagger Channel;t_{#gamma}-t_{#pi^{0}} (ns)", 352, 0, 352, 1400, -700, 700);

    Pi0_OA = new TH3D("Pi0_OA", "Pi0 Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Pi0_OA_R = new TH3D("Pi0_OA_R", "Pi0 Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Pi0_OA_Cut = new TH3D("Pi0_OA_Cut", "Pi0 Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Pi0_OA_Cut_R = new TH3D("Pi0_OA_Cut_R", "Pi0 Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);

    Pi0_MM_NE_0 = new TH3D("Pi0_MM_NE_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NE_0_R = new TH3D("Pi0_MM_NE_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NE_1 = new TH3D("Pi0_MM_NE_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NE_1_R = new TH3D("Pi0_MM_NE_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_NE_0 = new TH3D("Pi0_Ph_NE_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NE_0_R = new TH3D("Pi0_Ph_NE_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NE_1 = new TH3D("Pi0_Ph_NE_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NE_1_R = new TH3D("Pi0_Ph_NE_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_NI_0 = new TH3D("Pi0_MM_NI_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NI_0_R = new TH3D("Pi0_MM_NI_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NI_1 = new TH3D("Pi0_MM_NI_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_NI_1_R = new TH3D("Pi0_MM_NI_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_NI_0 = new TH3D("Pi0_Ph_NI_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NI_0_R = new TH3D("Pi0_Ph_NI_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NI_1 = new TH3D("Pi0_Ph_NI_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_NI_1_R = new TH3D("Pi0_Ph_NI_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_CE_0 = new TH3D("Pi0_MM_CE_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CE_0_R = new TH3D("Pi0_MM_CE_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CE_1 = new TH3D("Pi0_MM_CE_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CE_1_R = new TH3D("Pi0_MM_CE_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_CE_0 = new TH3D("Pi0_Ph_CE_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CE_0_R = new TH3D("Pi0_Ph_CE_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CE_1 = new TH3D("Pi0_Ph_CE_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CE_1_R = new TH3D("Pi0_Ph_CE_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_CI_0 = new TH3D("Pi0_MM_CI_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CI_0_R = new TH3D("Pi0_MM_CI_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CI_1 = new TH3D("Pi0_MM_CI_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_CI_1_R = new TH3D("Pi0_MM_CI_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_CI_0 = new TH3D("Pi0_Ph_CI_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CI_0_R = new TH3D("Pi0_Ph_CI_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CI_1 = new TH3D("Pi0_Ph_CI_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_CI_1_R = new TH3D("Pi0_Ph_CI_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_WE_0 = new TH3D("Pi0_MM_WE_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WE_0_R = new TH3D("Pi0_MM_WE_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WE_1 = new TH3D("Pi0_MM_WE_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WE_1_R = new TH3D("Pi0_MM_WE_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_WE_0 = new TH3D("Pi0_Ph_WE_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WE_0_R = new TH3D("Pi0_Ph_WE_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WE_1 = new TH3D("Pi0_Ph_WE_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WE_1_R = new TH3D("Pi0_Ph_WE_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_WI_0 = new TH3D("Pi0_MM_WI_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WI_0_R = new TH3D("Pi0_MM_WI_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WI_1 = new TH3D("Pi0_MM_WI_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_WI_1_R = new TH3D("Pi0_MM_WI_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_WI_0 = new TH3D("Pi0_Ph_WI_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WI_0_R = new TH3D("Pi0_Ph_WI_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WI_1 = new TH3D("Pi0_Ph_WI_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_WI_1_R = new TH3D("Pi0_Ph_WI_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_TE_0 = new TH3D("Pi0_MM_TE_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TE_0_R = new TH3D("Pi0_MM_TE_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TE_1 = new TH3D("Pi0_MM_TE_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TE_1_R = new TH3D("Pi0_MM_TE_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_TE_0 = new TH3D("Pi0_Ph_TE_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TE_0_R = new TH3D("Pi0_Ph_TE_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TE_1 = new TH3D("Pi0_Ph_TE_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TE_1_R = new TH3D("Pi0_Ph_TE_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_MM_TI_0 = new TH3D("Pi0_MM_TI_0", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TI_0_R = new TH3D("Pi0_MM_TI_0_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TI_1 = new TH3D("Pi0_MM_TI_1", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_MM_TI_1_R = new TH3D("Pi0_MM_TI_1_R", "Pi0 Missing Mass;Tagger Channel;#theta_{#pi^{0}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Pi0_Ph_TI_0 = new TH3D("Pi0_Ph_TI_0", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TI_0_R = new TH3D("Pi0_Ph_TI_0_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TI_1 = new TH3D("Pi0_Ph_TI_1", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Pi0_Ph_TI_1_R = new TH3D("Pi0_Ph_TI_1_R", "Pi0 Phi Distribution;Tagger Channel;#theta_{#pi^{0}} (deg);#phi_{#pi^{0}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Pi0_Re_All = new TH3D("Pi0_Re_All", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);
    Pi0_Re_All_R = new TH3D("Pi0_Re_All_R", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);
    Pi0_Re_Det = new TH3D("Pi0_Re_Det", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);
    Pi0_Re_Det_R = new TH3D("Pi0_Re_Det_R", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);
    Pi0_Re_Dif = new TH3D("Pi0_Re_Dif", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);E_{miss}-E_{det} (MeV)", 352, 0, 352, 12, 0, 60, 200, -50, 150);
    Pi0_Re_Dif_R = new TH3D("Pi0_Re_Dif_R", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);E_{miss}-E_{det} (MeV)", 352, 0, 352, 12, 0, 60, 200, -50, 150);
    Pi0_Re_NoE = new TH3D("Pi0_Re_NoE", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);
    Pi0_Re_NoE_R = new TH3D("Pi0_Re_NoE_R", "Pi0 Recoil Detection;Tagger Channel;#theta_{miss} (deg);Recoil Energy (MeV)", 352, 0, 352, 12, 0, 60, 200, 0, 200);

    // Pi+ histograms
    PiP_Sp = new TH2D("PiP_Sp", "Split OA vs Energy Ratio;E_{split}/E_{#gamma};#theta_{OA} (deg)", 200, 0, 1, 180, 0, 180);
    PiP_CA = new TH1D("PiP_CA", "Pi+ Recoil Coplanarity Angle;#phi_{#pi^{+}}-#phi_{p} (deg)", 360, 0, 360);

    PiP_Tm_NE = new TH2D("PiP_Tm_NE", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_NI = new TH2D("PiP_Tm_NI", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_CE = new TH2D("PiP_Tm_CE", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_CI = new TH2D("PiP_Tm_CI", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_WE = new TH2D("PiP_Tm_WE", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_WI = new TH2D("PiP_Tm_WI", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_TE = new TH2D("PiP_Tm_TE", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);
    PiP_Tm_TI = new TH2D("PiP_Tm_TI", "Pi+ Time;Tagger Channel;t_{#gamma}-t_{#pi^{+}} (ns)", 352, 0, 352, 1400, -700, 700);

    PiP_OA = new TH3D("PiP_OA", "Pi+ Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    PiP_OA_R = new TH3D("PiP_OA_R", "Pi+ Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    PiP_OA_Cut = new TH3D("PiP_OA_Cut", "Pi+ Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    PiP_OA_Cut_R = new TH3D("PiP_OA_Cut_R", "Pi+ Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);

    PiP_MM_NE_0 = new TH3D("PiP_MM_NE_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NE_0_R = new TH3D("PiP_MM_NE_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NE_1 = new TH3D("PiP_MM_NE_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NE_1_R = new TH3D("PiP_MM_NE_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_NE_0 = new TH3D("PiP_Ph_NE_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NE_0_R = new TH3D("PiP_Ph_NE_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NE_1 = new TH3D("PiP_Ph_NE_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NE_1_R = new TH3D("PiP_Ph_NE_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_NI_0 = new TH3D("PiP_MM_NI_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NI_0_R = new TH3D("PiP_MM_NI_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NI_1 = new TH3D("PiP_MM_NI_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_NI_1_R = new TH3D("PiP_MM_NI_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_NI_0 = new TH3D("PiP_Ph_NI_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NI_0_R = new TH3D("PiP_Ph_NI_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NI_1 = new TH3D("PiP_Ph_NI_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_NI_1_R = new TH3D("PiP_Ph_NI_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_CE_0 = new TH3D("PiP_MM_CE_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CE_0_R = new TH3D("PiP_MM_CE_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CE_1 = new TH3D("PiP_MM_CE_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CE_1_R = new TH3D("PiP_MM_CE_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_CE_0 = new TH3D("PiP_Ph_CE_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CE_0_R = new TH3D("PiP_Ph_CE_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CE_1 = new TH3D("PiP_Ph_CE_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CE_1_R = new TH3D("PiP_Ph_CE_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_CI_0 = new TH3D("PiP_MM_CI_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CI_0_R = new TH3D("PiP_MM_CI_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CI_1 = new TH3D("PiP_MM_CI_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_CI_1_R = new TH3D("PiP_MM_CI_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_CI_0 = new TH3D("PiP_Ph_CI_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CI_0_R = new TH3D("PiP_Ph_CI_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CI_1 = new TH3D("PiP_Ph_CI_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_CI_1_R = new TH3D("PiP_Ph_CI_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_WE_0 = new TH3D("PiP_MM_WE_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WE_0_R = new TH3D("PiP_MM_WE_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WE_1 = new TH3D("PiP_MM_WE_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WE_1_R = new TH3D("PiP_MM_WE_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_WE_0 = new TH3D("PiP_Ph_WE_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WE_0_R = new TH3D("PiP_Ph_WE_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WE_1 = new TH3D("PiP_Ph_WE_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WE_1_R = new TH3D("PiP_Ph_WE_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_WI_0 = new TH3D("PiP_MM_WI_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WI_0_R = new TH3D("PiP_MM_WI_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WI_1 = new TH3D("PiP_MM_WI_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_WI_1_R = new TH3D("PiP_MM_WI_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_WI_0 = new TH3D("PiP_Ph_WI_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WI_0_R = new TH3D("PiP_Ph_WI_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WI_1 = new TH3D("PiP_Ph_WI_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_WI_1_R = new TH3D("PiP_Ph_WI_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_TE_0 = new TH3D("PiP_MM_TE_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TE_0_R = new TH3D("PiP_MM_TE_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TE_1 = new TH3D("PiP_MM_TE_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TE_1_R = new TH3D("PiP_MM_TE_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_TE_0 = new TH3D("PiP_Ph_TE_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TE_0_R = new TH3D("PiP_Ph_TE_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TE_1 = new TH3D("PiP_Ph_TE_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TE_1_R = new TH3D("PiP_Ph_TE_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    PiP_MM_TI_0 = new TH3D("PiP_MM_TI_0", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TI_0_R = new TH3D("PiP_MM_TI_0_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TI_1 = new TH3D("PiP_MM_TI_1", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_MM_TI_1_R = new TH3D("PiP_MM_TI_1_R", "Pi+ Missing Mass;Tagger Channel;#theta_{#pi^{+}} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    PiP_Ph_TI_0 = new TH3D("PiP_Ph_TI_0", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TI_0_R = new TH3D("PiP_Ph_TI_0_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TI_1 = new TH3D("PiP_Ph_TI_1", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    PiP_Ph_TI_1_R = new TH3D("PiP_Ph_TI_1_R", "Pi+ Phi Distribution;Tagger Channel;#theta_{#pi^{+}} (deg);#phi_{#pi^{+}} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    // Compton histograms
    Comp_Sp = new TH2D("Comp_Sp", "Split OA vs Energy Ratio;E_{split}/E_{#gamma};#theta_{OA} (deg)", 200, 0, 1, 180, 0, 180);
    Comp_CA = new TH1D("Comp_CA", "Compton Recoil Coplanarity Angle;#phi_{#gamma}-#phi_{p} (deg)", 360, 0, 360);

    Comp_Tm_NE = new TH2D("Comp_Tm_NE", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_NI = new TH2D("Comp_Tm_NI", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_CE = new TH2D("Comp_Tm_CE", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_CI = new TH2D("Comp_Tm_CI", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_WE = new TH2D("Comp_Tm_WE", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_WI = new TH2D("Comp_Tm_WI", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_TE = new TH2D("Comp_Tm_TE", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);
    Comp_Tm_TI = new TH2D("Comp_Tm_TI", "Compton Time;Tagger Channel;t_{#gamma}-t_{#gamma^{'}} (ns)", 352, 0, 352, 1400, -700, 700);

    Comp_OA = new TH3D("Comp_OA", "Compton Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Comp_OA_R = new TH3D("Comp_OA_R", "Compton Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Comp_OA_Cut = new TH3D("Comp_OA_Cut", "Compton Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);
    Comp_OA_Cut_R = new TH3D("Comp_OA_Cut_R", "Compton Recoil Opening Angle;Tagger Channel;#theta_{miss} (deg);Opening Angle (deg)", 352, 0, 352, 12, 0, 60, 90, 0, 180);

    Comp_MM_NE_0 = new TH3D("Comp_MM_NE_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NE_0_R = new TH3D("Comp_MM_NE_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NE_1 = new TH3D("Comp_MM_NE_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NE_1_R = new TH3D("Comp_MM_NE_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NE_0 = new TH3D("Comp_Ph_NE_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NE_0_R = new TH3D("Comp_Ph_NE_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NE_1 = new TH3D("Comp_Ph_NE_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NE_1_R = new TH3D("Comp_Ph_NE_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_NI_0 = new TH3D("Comp_MM_NI_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NI_0_R = new TH3D("Comp_MM_NI_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NI_1 = new TH3D("Comp_MM_NI_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_NI_1_R = new TH3D("Comp_MM_NI_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_NI_0 = new TH3D("Comp_Ph_NI_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NI_0_R = new TH3D("Comp_Ph_NI_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NI_1 = new TH3D("Comp_Ph_NI_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_NI_1_R = new TH3D("Comp_Ph_NI_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_CE_0 = new TH3D("Comp_MM_CE_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CE_0_R = new TH3D("Comp_MM_CE_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CE_1 = new TH3D("Comp_MM_CE_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CE_1_R = new TH3D("Comp_MM_CE_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_CE_0 = new TH3D("Comp_Ph_CE_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CE_0_R = new TH3D("Comp_Ph_CE_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CE_1 = new TH3D("Comp_Ph_CE_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CE_1_R = new TH3D("Comp_Ph_CE_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_CI_0 = new TH3D("Comp_MM_CI_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CI_0_R = new TH3D("Comp_MM_CI_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CI_1 = new TH3D("Comp_MM_CI_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_CI_1_R = new TH3D("Comp_MM_CI_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_CI_0 = new TH3D("Comp_Ph_CI_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CI_0_R = new TH3D("Comp_Ph_CI_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CI_1 = new TH3D("Comp_Ph_CI_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_CI_1_R = new TH3D("Comp_Ph_CI_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_WE_0 = new TH3D("Comp_MM_WE_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WE_0_R = new TH3D("Comp_MM_WE_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WE_1 = new TH3D("Comp_MM_WE_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WE_1_R = new TH3D("Comp_MM_WE_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_WE_0 = new TH3D("Comp_Ph_WE_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WE_0_R = new TH3D("Comp_Ph_WE_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WE_1 = new TH3D("Comp_Ph_WE_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WE_1_R = new TH3D("Comp_Ph_WE_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_WI_0 = new TH3D("Comp_MM_WI_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WI_0_R = new TH3D("Comp_MM_WI_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WI_1 = new TH3D("Comp_MM_WI_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_WI_1_R = new TH3D("Comp_MM_WI_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_WI_0 = new TH3D("Comp_Ph_WI_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WI_0_R = new TH3D("Comp_Ph_WI_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WI_1 = new TH3D("Comp_Ph_WI_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_WI_1_R = new TH3D("Comp_Ph_WI_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_TE_0 = new TH3D("Comp_MM_TE_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TE_0_R = new TH3D("Comp_MM_TE_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TE_1 = new TH3D("Comp_MM_TE_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TE_1_R = new TH3D("Comp_MM_TE_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_TE_0 = new TH3D("Comp_Ph_TE_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TE_0_R = new TH3D("Comp_Ph_TE_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TE_1 = new TH3D("Comp_Ph_TE_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TE_1_R = new TH3D("Comp_Ph_TE_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_MM_TI_0 = new TH3D("Comp_MM_TI_0", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TI_0_R = new TH3D("Comp_MM_TI_0_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TI_1 = new TH3D("Comp_MM_TI_1", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_MM_TI_1_R = new TH3D("Comp_MM_TI_1_R", "Compton Missing Mass;Tagger Channel;#theta_{#gamma} (deg);m_{miss}-m_{targ} (MeV)", 352, 0, 352, 36, 0, 180, 80, -80, 120);
    Comp_Ph_TI_0 = new TH3D("Comp_Ph_TI_0", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TI_0_R = new TH3D("Comp_Ph_TI_0_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TI_1 = new TH3D("Comp_Ph_TI_1", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);
    Comp_Ph_TI_1_R = new TH3D("Comp_Ph_TI_1_R", "Compton Phi Distribution;Tagger Channel;#theta_{#gamma} (deg);#phi_{#gamma} (deg)", 352, 0, 352, 36, 0, 180, 72, -180, 180);

    Comp_CS = new TH3D("Comp_CS", "Compton Cluster Size;E_{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Comp_CS_MM = new TH3D("Comp_CS_MM", "Compton Cluster Size Cut by MM;E_{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Comp_CS_MM_R = new TH3D("Comp_CS_MM_R", "Compton Cluster Size Cut by MM;E_{#gamma} (MeV);#theta_{#gamma} (deg);Cluster Size", 300, 0, 300, 36, 0, 180, 12, 0, 12);
    Reco_CS = new TH3D("Reco_CS", "Recoil Cluster Size;E_{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 18, 0, 90, 12, 0, 12);
    Reco_CS_MM = new TH3D("Reco_CS_MM", "Recoil Cluster Size Cut by MM;E_{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 18, 0, 90, 12, 0, 12);
    Reco_CS_MM_R = new TH3D("Reco_CS_MM_R", "Recoil Cluster Size Cut by MM;E_{p} (MeV);#theta_{p} (deg);Cluster Size", 300, 0, 300, 18, 0, 90, 12, 0, 12);

    verbosity = 0;
    match_charge = true;

    IMCut = 134.98;
    MMLoC = 800.00;
    MMHiC = 938.27;
    OACut = 180;
    ESCut = 0;

    save_randoms = false;
    split_search = false;
    pure_mwpc = true;

    taps_eff = 1;

    cir_beam = false;
    lin_beam = false;
    beamPol = 1;
    targPol = 1;

    firstEvent = true;
    gRandom->SetSeed(0);
}

PAnalyze::~PAnalyze()
{
}

Bool_t	PAnalyze::Init()
{
    cout << "Initialising physics analysis..." << endl;
    cout << "--------------------------------------------------" << endl << endl;

    if(!InitBackgroundCuts()) return kFALSE;
    if(!InitTargetMass()) return kFALSE;
    if(!InitVerbosity()) return kFALSE;
    if(!InitMatchCharge()) return kFALSE;
    if(!InitInvariantMass()) return kFALSE;
    if(!InitMissingMass()) return kFALSE;
    if(!InitOpeningAngle()) return kFALSE;
    if(!InitEnergySum()) return kFALSE;
    if(!InitSaveRandoms()) return kFALSE;
    if(!InitSplitSearch()) return kFALSE;
    if(!InitPureMWPC()) return kFALSE;
    if(!InitBeamPol()) return kFALSE;
    if(!InitTargPol()) return kFALSE;

    if(!PPhysics::Init()) return kFALSE;

    TaggerAccScal = GetScalerHist("TaggerAccScal");
    if(!TaggerAccScal)
    {
        cout << "No tagger scaler histogram available" << endl;
        return kFALSE;
    }
    LiveTimeScal = GetScalerHist("LiveTimeScal");
    if(!LiveTimeScal)
    {
        cout << "No live time histogram available" << endl;
        return kFALSE;
    }

    cout << "--------------------------------------------------" << endl;
    return kTRUE;
}

Bool_t	PAnalyze::Start()
{
    if(!IsAcquFile())
    {
        cout << "ERROR: Input File is not an Acqu file." << endl;
        return kFALSE;
    }
    SetAsGoATFile();

    TraverseValidEvents();

    //GoosyTagger(TaggerAccScal);
    //GoosyVuprom(TaggerAccScal);
    GoosyNewFPDRecabled(TaggerAccScal);

    return kTRUE;
}

void	PAnalyze::ProcessEvent()
{
    //////////////////////////////////////////////////
    // Read in polarization info on first event
    //////////////////////////////////////////////////
    if (firstEvent)
    {
        Int_t fileTime = GetSetupParameters()->GetTimeStamp();
        cout << endl << "Time Stamp for this run: " << fileTime << endl;
        Int_t timeDiff = fileTime;

        for (UInt_t i = 0; i<beamPolTime.size(); i++)
        {
            if (TMath::Abs(beamPolTime.at(i) - fileTime) < timeDiff)
            {
                beamPol = beamPolMeas.at(i);
                timeDiff = TMath::Abs(beamPolTime.at(i) - fileTime);
            }
        }
        for (UInt_t i = 0; i<targPolTime.size(); i+=2)
        {
            if ((targPolTime.at(i) < fileTime) && (fileTime < targPolTime.at(i+1)))
            {
                Double_t tau = ((targPolTime.at(i)-targPolTime.at(i+1))/TMath::Log(targPolMeas.at(i+1)/targPolMeas.at(i)));
                targPol = targPolMeas.at(i)*TMath::Exp((targPolTime.at(i)-fileTime)/tau);
                break;
            }
        }

        cout << "E Beam polarization for this run: " << 100*beamPol << "%" << endl;
        cout << "Target polarization for this run: " << 100*targPol << "%" << endl << endl;
        firstEvent = false;
    }

    if (GetEventNumber()%100000 == 0) cout << "Event " << GetEventNumber() << endl;

    //////////////////////////////////////////////////
    // Reject events below energy sum cut
    //////////////////////////////////////////////////
    if (GetTrigger()->GetEnergySum() < ESCut) return;

    //////////////////////////////////////////////////
    // Event weighting for MC
    //////////////////////////////////////////////////
    Float_t event_weight = 1;
    if (IsMCFile()) event_weight = GetTrigger()->GetMCWeight();

    //////////////////////////////////////////////////
    // Figure out 'helicity' state, reject errors
    //////////////////////////////////////////////////
    Bool_t b_hel = 0;
    if (cir_beam && !(IsMCFile()))
    {
        if(!(GetTrigger()->HasHelicity())) return;
        b_hel = GetTrigger()->GetHelicity();
    }
    for (Int_t i=0; i<(GetTrigger()->GetNErrors()); i++)
    {
        if (!(IsMCFile()) && (cir_beam || (GetTrigger()->GetErrorModuleID())[i] >= 0 || (GetTrigger()->GetErrorModuleIndex())[i] >= 0))
        {
            Helicity->Fill(2);
            return;
        }
    }
    Helicity->Fill(b_hel);

    //////////////////////////////////////////////////
    // Kinematic variables
    //////////////////////////////////////////////////
    Int_t n_accept, n_ignore, i_trk0, i_trk1, i_splt, i_part_sz, i_reco_sz, i_tagg_ch;
    Double_t d_part_tm, d_tagg_tm, d_aver_tm, d_subt_tm;
    Double_t d_tagg_en, d_trk0_en, d_trk1_en, d_part_en, d_reco_en, d_splt_en, d_miss_ma;
    Double_t d_part_th, d_part_ph, d_reco_th, d_miss_th, d_CA, d_OA, d_temp;
    TLorentzVector lv_trk0, lv_trk1, lv_splt, lv_part, lv_ptot, lv_beam, lv_miss, lv_piP;
    TLorentzVector lv_targ = GetTarget();
    TVector3 v_reco, v_splt;

    // Currently unused center of mass stuff
    //TLorentzVector lv_part_cm;
    //TVector3 v_lab_cm;
    //Double_t d_part_th_cm, d_part_ph_cm;

    Bool_t b_pi0, b_piP, b_comp, b_pi0_CC, b_pi0_CT, b_pi0_TC, b_NE, b_NI, b_CE, b_CI, b_TE, b_TI, b_WE, b_WI, b_cut_CA, b_cut_OA;
    b_pi0 = b_piP = b_comp = b_pi0_CC = b_pi0_CT = b_pi0_TC = b_NE = b_NI = b_CE = b_CI = b_TE = b_TI = b_WE = b_WI = b_cut_CA = b_cut_OA = false;

    n_accept = n_ignore = 0;
    i_trk0 = i_trk1 = i_splt = -1;
    Double_t d_min_IM = 0;

    ignoreTrack.clear();

    //////////////////////////////////////////////////
    // Check for tracks in MWPC and TAPS
    //////////////////////////////////////////////////
    for (Int_t i=0; i<(GetTracks()->GetNTracks()); i++)
    {
        if (((GetTracks()->GetClusterEnergy(i) == 0) && !pure_mwpc) ||
            ((GetTracks()->HasTAPS(i)) && (gRandom->Rndm() >= taps_eff)))
        {
            n_ignore++;
            ignoreTrack.push_back(true);
        }
        else
        {
            n_accept++;
            ignoreTrack.push_back(false);
        }
    }

    //////////////////////////////////////////////////
    // Initial pi0 stuff
    //////////////////////////////////////////////////
    for (Int_t i=0; i<(GetTracks()->GetNTracks()-1); i++)
    {
        d_trk0_en = GetTracks()->GetClusterEnergy(i);
        if (d_trk0_en == 0 || ignoreTrack.at(i)) continue;
        if (GetTracks()->IsCharged(i) && match_charge) continue;
        lv_trk0 = GetTracks()->GetVector(i);

        for (Int_t j=i+1; j<(GetTracks()->GetNTracks()); j++)
        {
            d_trk1_en = GetTracks()->GetClusterEnergy(j);
            if (d_trk1_en == 0 || ignoreTrack.at(j)) continue;
            if (GetTracks()->IsCharged(j) && match_charge) continue;
            lv_trk1 = GetTracks()->GetVector(j);

            lv_ptot = lv_trk0 + lv_trk1;
            if (TMath::Abs(134.98 - lv_ptot.M()) < TMath::Abs(134.98 - d_min_IM))
            {
                if ((GetTracks()->HasCB(i)) && (GetTracks()->HasCB(j)))
                {
                    d_part_tm = 0.5*(GetTracks()->GetTime(i) + GetTracks()->GetTime(j));
                    b_pi0_CC = true;
                    b_pi0_CT = b_pi0_TC = false;
                }
                else if (GetTracks()->HasCB(i))
                {
                    d_part_tm = (GetTracks()->GetTime(i));
                    b_pi0_CT = true;
                    b_pi0_CC = b_pi0_TC = false;
                }
                else if (GetTracks()->HasCB(j))
                {
                    d_part_tm = (GetTracks()->GetTime(j));
                    b_pi0_TC = true;
                    b_pi0_CC = b_pi0_CT = false;
                }
                else continue;
                i_trk0 = i;
                i_trk1 = j;
                i_splt = -1;
                lv_part = lv_trk0 + lv_trk1;                
                d_part_th = lv_part.Theta()*TMath::RadToDeg();
                d_part_ph = lv_part.Phi()*TMath::RadToDeg();
                d_min_IM = lv_part.M();

                b_NE = b_NI = false;
                if (n_accept == 2) b_NE = true;
                else b_NI = true;
            }

            if (!split_search) continue;
            //////////////////////////////////////////////////
            // Look for splits
            //////////////////////////////////////////////////
            for (Int_t k=0; k<(GetTracks()->GetNTracks()); k++)
            {
                d_splt_en = GetTracks()->GetClusterEnergy(k);
                if (k == i || k == j || d_splt_en == 0 || ignoreTrack.at(k)) continue;
                //if (GetTracks()->IsCharged(k) && match_charge) continue;
                lv_splt = GetTracks()->GetVector(k);
                v_splt = GetTracks()->GetUnitVector(k);

                Pi0_Sp->Fill(d_splt_en/d_trk0_en,(TMath::RadToDeg()*lv_trk0.Angle(v_splt)));
                Pi0_Sp->Fill(d_splt_en/d_trk1_en,(TMath::RadToDeg()*lv_trk1.Angle(v_splt)));

                if ((TMath::RadToDeg()*lv_trk0.Angle(v_splt)) < 40*(1-(2*d_splt_en/d_trk0_en)) ||
                        (TMath::RadToDeg()*lv_trk1.Angle(v_splt)) < 40*(1-(2*d_splt_en/d_trk1_en)))
                {
                    lv_ptot = lv_trk0 + lv_trk1 + lv_splt;
                    if (TMath::Abs(134.98 - lv_ptot.M()) < TMath::Abs(134.98 - d_min_IM))
                    {
                        if ((GetTracks()->HasCB(i)) && (GetTracks()->HasCB(j)))
                        {
                            d_part_tm = 0.5*(GetTracks()->GetTime(i) + GetTracks()->GetTime(j));
                            b_pi0_CC = true;
                            b_pi0_CT = b_pi0_TC = false;
                        }
                        else if (GetTracks()->HasCB(i))
                        {
                            d_part_tm = (GetTracks()->GetTime(i));
                            b_pi0_CT = true;
                            b_pi0_CC = b_pi0_TC = false;
                        }
                        else if (GetTracks()->HasCB(j))
                        {
                            d_part_tm = (GetTracks()->GetTime(j));
                            b_pi0_TC = true;
                            b_pi0_CC = b_pi0_CT = false;
                        }
                        else continue;
                        i_trk0 = i;
                        i_trk1 = j;
                        i_splt = k;
                        lv_part = lv_trk0 + lv_trk1 + lv_splt;
                        d_part_th = lv_part.Theta()*TMath::RadToDeg();
                        d_part_ph = lv_part.Phi()*TMath::RadToDeg();
                        //pi0_XX = lv_trk0 + lv_trk1;
                        d_min_IM = lv_part.M();

                        b_NE = b_NI = false;
                        if (n_accept == 3) b_NE = true;
                        else b_NI = true;
                    }
                }
            }
        }
    }

    //////////////////////////////////////////////////
    // Found a pi0 candidate
    //////////////////////////////////////////////////
    if (d_min_IM > 0)
    {
        b_pi0 = (TMath::Abs(134.98 - d_min_IM) < IMCut);

        Pi0_IM_A->Fill(d_min_IM, event_weight);
        if (b_NE) Pi0_IM_E->Fill(d_min_IM, event_weight);
        else if(b_NI) Pi0_IM_I->Fill(d_min_IM, event_weight);

        if (b_pi0_CC) Pi0_IM_CC->Fill(d_trk0_en, d_trk1_en, d_min_IM, event_weight);
        else if (b_pi0_CT) Pi0_IM_CT->Fill(d_trk0_en, d_trk1_en, d_min_IM, event_weight);
        else if (b_pi0_TC) Pi0_IM_CT->Fill(d_trk1_en, d_trk0_en, d_min_IM, event_weight);

        if (b_pi0)
        {
            //////////////////////////////////////////////////
            // Look for best recoil
            //////////////////////////////////////////////////
            d_CA = 360;

            for (Int_t i=0; i<(GetTracks()->GetNTracks()); i++)
            {
                if (i == i_trk0 || i == i_trk1 || i == i_splt || ignoreTrack.at(i)) continue;
                if (GetTracks()->IsNeutral(i) && match_charge) continue;
                d_temp = TMath::Abs(d_part_ph-GetTracks()->GetPhi(i));
                if (TMath::Abs(180.0-d_temp) < TMath::Abs(180.0-d_CA))
                {
                    d_CA = d_temp;
                    if (TMath::Abs(180.0-d_CA) < OACut)
                    {
                        b_cut_CA = true;
                        v_reco = GetTracks()->GetUnitVector(i);
                        d_reco_en = GetTracks()->GetClusterEnergy(i);
                        d_reco_th = GetTracks()->GetTheta(i);
                        i_reco_sz = GetTracks()->GetClusterSize(i);

                        b_NE = b_NI = b_CE = b_CI = b_WE = b_WI = b_TE = b_TI = false;
                        if (GetTracks()->HasCB(i) && d_reco_en > 0)
                        {
                            if (i_splt >= 0 && n_accept == 4) b_CE = true;
                            else if (i_splt < 0 && n_accept == 3) b_CE = true;
                            else b_CI = true;
                        }
                        else if (GetTracks()->HasCB(i))
                        {
                            if (i_splt >= 0 && n_accept == 4) b_WE = true;
                            else if (i_splt < 0 && n_accept == 3) b_WE = true;
                            else b_WI = true;
                        }
                        else
                        {
                            if (i_splt >= 0 && n_accept == 4) b_TE = true;
                            else if (i_splt < 0 && n_accept == 3) b_TE = true;
                            else b_TI = true;
                        }
                    }
                }
            }

            Pi0_CA->Fill(d_CA, event_weight);
        }
    }

    //////////////////////////////////////////////////
    // Initial pi+/compton stuff
    //////////////////////////////////////////////////
    if (!b_pi0)
    {
        i_trk0 = i_trk1 = i_splt = -1;
        d_part_en = 0;
        d_reco_en = 0;

        for (Int_t i=0; i<(GetTracks()->GetNTracks()); i++)
        {
            d_trk0_en = GetTracks()->GetClusterEnergy(i);
            if (d_trk0_en <= d_part_en || ignoreTrack.at(i)) continue;
            b_piP = ((GetTracks()->IsCharged(i) && match_charge) || !match_charge);
            b_comp = ((GetTracks()->IsNeutral(i) && match_charge) || !match_charge);
            i_trk0 = i;
            i_splt = -1;
            lv_part = GetTracks()->GetVector(i);
            d_part_tm = GetTracks()->GetTime(i);
            d_part_en = GetTracks()->GetClusterEnergy(i);
            d_part_th = GetTracks()->GetTheta(i);
            d_part_ph = GetTracks()->GetPhi(i);
            i_part_sz = GetTracks()->GetClusterSize(i);

            b_NE = b_NI = false;
            if (n_accept == 1) b_NE = true;
            else b_NI = true;

            if (!split_search) continue;
            //////////////////////////////////////////////////
            // Look for splits
            //////////////////////////////////////////////////
            for (Int_t j=0; j<(GetTracks()->GetNTracks()); j++)
            {
                d_splt_en = GetTracks()->GetClusterEnergy(j);
                if (j == i || d_splt_en == 0 || ignoreTrack.at(j)) continue;
                //if (GetTracks()->IsNeutral(j) && match_charge && b_piP) continue;
                //if (GetTracks()->IsCharged(j) && match_charge && b_comp) continue;
                lv_splt = GetTracks()->GetVector(j);
                v_splt = GetTracks()->GetUnitVector(j);

                if (b_piP) PiP_Sp->Fill(d_splt_en/d_part_en,(TMath::RadToDeg()*lv_part.Angle(v_splt)));
                else if (b_comp) Comp_Sp->Fill(d_splt_en/d_part_en,(TMath::RadToDeg()*lv_part.Angle(v_splt)));
                else
                {
                    PiP_Sp->Fill(d_splt_en/d_part_en,(TMath::RadToDeg()*lv_part.Angle(v_splt)));
                    Comp_Sp->Fill(d_splt_en/d_part_en,(TMath::RadToDeg()*lv_part.Angle(v_splt)));
                }

                if ((TMath::RadToDeg()*lv_part.Angle(v_splt)) < 40*(1-(2*d_splt_en/d_trk0_en)))
                {
                    i_trk0 = i;
                    i_splt = j;
                    lv_ptot = lv_part + lv_splt;
                    lv_part = AdjustMass(lv_ptot, 0);
                    d_part_en = lv_part.E();
                    d_part_th = lv_part.Theta()*TMath::RadToDeg();
                    d_part_ph = lv_part.Phi()*TMath::RadToDeg();
                    i_part_sz += GetTracks()->GetClusterSize(j);

                    b_NE = b_NI = false;
                    if (n_accept == 2) b_NE = true;
                    else b_NI = true;
                }
            }
        }

        //////////////////////////////////////////////////
        // Look for best recoil
        //////////////////////////////////////////////////
        d_CA = 360;

        for (Int_t i=0; i<(GetTracks()->GetNTracks()); i++)
        {
            if (i_trk0 < 0 || i == i_trk0 || i == i_splt || d_part_en == 0 || ignoreTrack.at(i)) continue;
            if (GetTracks()->IsCharged(i) && match_charge && b_piP) continue;
            if (GetTracks()->IsNeutral(i) && match_charge && b_comp) continue;
            d_temp = TMath::Abs(d_part_ph-GetTracks()->GetPhi(i));
            if (TMath::Abs(180.0-d_temp) < TMath::Abs(180.0-d_CA))
            {
                d_CA = d_temp;
                if (TMath::Abs(180.0-d_CA) < OACut)
                {
                    b_cut_CA = true;
                    v_reco = GetTracks()->GetUnitVector(i);
                    d_reco_en = GetTracks()->GetClusterEnergy(i);
                    d_reco_th = GetTracks()->GetTheta(i);
                    i_reco_sz = GetTracks()->GetClusterSize(i);

                    b_NE = b_NI = b_CE = b_CI = b_WE = b_WI = b_TE = b_TI = false;
                    if (GetTracks()->HasCB(i) && d_reco_en > 0)
                    {
                        if (i_splt >= 0 && n_accept == 3) b_CE = true;
                        else if (i_splt < 0 && n_accept == 2) b_CE = true;
                        else b_CI = true;
                    }
                    else if (GetTracks()->HasCB(i))
                    {
                        if (i_splt >= 0 && n_accept == 3) b_WE = true;
                        else if (i_splt < 0 && n_accept == 2) b_WE = true;
                        else b_WI = true;
                    }
                    else
                    {
                        if (i_splt >= 0 && n_accept == 3) b_TE = true;
                        else if (i_splt < 0 && n_accept == 2) b_TE = true;
                        else b_TI = true;
                    }
                }
            }
        }

        if (b_piP)
        {
            PiP_CA->Fill(d_CA, event_weight);

            lv_piP = AdjustMass(lv_part, 139.57);
        }
        if (b_comp)
        {
            Comp_CA->Fill(d_CA, event_weight);

            if (b_cut_CA)
            {
                if (d_reco_en > 0)
                {
                    Comp_CS->Fill(d_part_en, d_part_th, i_part_sz, event_weight);
                    Reco_CS->Fill(d_reco_en, d_reco_th, i_reco_sz, event_weight);
                }
            }
        }
    }

    //////////////////////////////////////////////////
    // Decide which state to fill
    //////////////////////////////////////////////////
    Bool_t b_fill_Ph_0, b_fill_Ph_1, b_fill_MM_0, b_fill_MM_1;
    b_fill_Ph_0 = b_fill_Ph_1 = b_fill_MM_0 = b_fill_MM_1 = false;
    if (lin_beam)
    {
        if (d_part_ph < -90 || (d_part_ph > 0 && d_part_ph < 90)) b_fill_MM_0 = true;
        else b_fill_MM_1 = true;
    }
    else
    {
        if(b_hel) b_fill_MM_1 = true;
        else b_fill_MM_0 = true;
    }

    //////////////////////////////////////////////////
    // Get average time for total inclusive
    //////////////////////////////////////////////////
    Int_t n_clus = 0;
    d_aver_tm = 0;
    for (Int_t i = 0; i < GetTracks()->GetNTracks(); i++)
    {
        if (GetTracks()->GetClusterEnergy(i) != 0)
        {
            d_aver_tm += GetTracks()->GetTime(i);
            n_clus++;
        }
    }
    d_aver_tm = d_aver_tm/n_clus;

    //////////////////////////////////////////////////
    // Tagger loop
    //////////////////////////////////////////////////
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        if (RejectTagged(i)) continue;

        i_tagg_ch = GetTagger()->GetTaggedChannel(i);
        d_tagg_tm = GetTagger()->GetTaggedTime(i);
        d_tagg_en = GetTagger()->GetTaggedEnergy(i);
        lv_beam = TLorentzVector(0., 0., d_tagg_en, d_tagg_en);

        //////////////////////////////////////////////////
        // Total inclusive stuff
        //////////////////////////////////////////////////
        Tagg_Tm->Fill(d_tagg_tm);
        if (b_hel)
        {
            if (GHistBGSub::IsPrompt(d_tagg_tm)) Tagg_1->Fill(i_tagg_ch);
            else if (GHistBGSub::IsRandom(d_tagg_tm)) Tagg_1_R->Fill(i_tagg_ch);
        }
        else
        {
            if (GHistBGSub::IsPrompt(d_tagg_tm)) Tagg_0->Fill(i_tagg_ch);
            else if (GHistBGSub::IsRandom(d_tagg_tm)) Tagg_0_R->Fill(i_tagg_ch);
        }

        d_subt_tm = d_tagg_tm - d_aver_tm;
        Inc_Tm->Fill(d_subt_tm);
        if (b_hel)
        {
            if (GHistBGSub::IsPrompt(d_subt_tm)) Inc_1->Fill(i_tagg_ch, GetTracks()->GetNTracks());
            else if (GHistBGSub::IsRandom(d_subt_tm)) Inc_1_R->Fill(i_tagg_ch, GetTracks()->GetNTracks());
        }
        else
        {
            if (GHistBGSub::IsPrompt(d_subt_tm)) Inc_0->Fill(i_tagg_ch, GetTracks()->GetNTracks());
            else if (GHistBGSub::IsRandom(d_subt_tm)) Inc_0_R->Fill(i_tagg_ch, GetTracks()->GetNTracks());
        }

        //////////////////////////////////////////////////
        // Good pi0 event
        //////////////////////////////////////////////////
        if (b_pi0)
        {
            d_subt_tm = d_tagg_tm - d_part_tm;
            lv_miss = lv_beam + lv_targ - lv_part;
            d_miss_ma = lv_miss.M()-lv_targ.M();
            d_miss_th = lv_miss.Theta()*TMath::RadToDeg();

            b_cut_OA = false;
            if (b_cut_CA)
            {
                d_OA = (TMath::RadToDeg()*lv_miss.Angle(v_reco));
                b_cut_OA = (d_OA < OACut);
            }

            //////////////////////////////////////////////////
            // Recoil detection stuff
            //////////////////////////////////////////////////
            if (GHistBGSub::IsPrompt(d_subt_tm))
            {
                //Pi0_MM->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                if (b_cut_CA) Pi0_OA->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    Pi0_Re_All->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M(), event_weight);
                    if (b_cut_CA) Pi0_OA_Cut->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                    if (b_cut_OA)
                    {
                        if (d_reco_en > 0)
                        {
                            Pi0_Re_Det->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M(), event_weight);
                            Pi0_Re_Dif->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M()-d_reco_en, event_weight);
                        }
                        else Pi0_Re_NoE->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M(), event_weight);
                    }
                }
            }
            else if (GHistBGSub::IsRandom(d_subt_tm))
            {
                //Pi0_MM_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                if (b_cut_CA) Pi0_OA_R->Fill(i_tagg_ch, d_miss_th, d_OA,event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    Pi0_Re_All_R->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M(), event_weight);
                    if (b_cut_CA) Pi0_OA_Cut_R->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                    if (b_cut_OA)
                    {
                        if (d_reco_en > 0)
                        {
                            Pi0_Re_Det_R->Fill(i_tagg_ch, d_miss_th, lv_miss.E()-lv_miss.M(), event_weight);
                            Pi0_Re_Dif_R->Fill(i_tagg_ch, d_miss_th,lv_miss.E()-lv_miss.M()-d_reco_en,  event_weight);
                        }
                        else Pi0_Re_NoE_R->Fill(i_tagg_ch, d_miss_th,lv_miss.E()-lv_miss.M(), event_weight);
                    }
                }
            }

            //////////////////////////////////////////////////
            // State selection
            //////////////////////////////////////////////////
            if (!b_hel) b_fill_Ph_0 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);
            else b_fill_Ph_1 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);

            // Exclusive event, without recoil detected
            if (b_NE)
            {
                Pi0_Tm_NE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_NE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_NE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_NE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_NE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_NE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_NE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_NE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_NE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, without recoil detected (add all of the now failed recoil ones in)
            else if (b_NI || !b_cut_OA)
            {
                Pi0_Tm_NI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_NI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_NI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_NI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_NI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_NI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_NI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_NI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_NI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // From here on the opening angle must have been satisfied
            // Exclusive event, with recoil detected in CB
            else if (b_CE)
            {
                Pi0_Tm_CE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_CE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_CE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_CE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_CE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_CE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_CE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_CE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_CE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in CB
            else if (b_CI)
            {
                Pi0_Tm_CI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_CI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_CI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_CI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_CI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_CI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_CI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_CI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_CI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in MWPC
            else if (b_WE)
            {
                Pi0_Tm_WE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_WE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_WE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_WE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_WE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_WE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_WE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_WE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_WE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in MWPC
            else if (b_WI)
            {
                Pi0_Tm_WI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_WI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_WI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_WI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_WI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_WI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_WI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_WI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_WI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in TAPS
            else if (b_TE)
            {
                Pi0_Tm_TE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_TE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_TE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_TE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_TE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_TE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_TE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_TE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_TE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in TAPS
            else if (b_TI)
            {
                Pi0_Tm_TI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_TI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_TI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_TI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_TI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Pi0_MM_TI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Pi0_MM_TI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Pi0_Ph_TI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Pi0_Ph_TI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
        }

        //////////////////////////////////////////////////
        // Good pi+ event
        //////////////////////////////////////////////////
        if (b_piP)
        {
            d_subt_tm = d_tagg_tm - d_part_tm;
            lv_miss = lv_beam + lv_targ - lv_piP;
            d_miss_ma = lv_miss.M()-lv_targ.M();
            d_miss_th = lv_miss.Theta()*TMath::RadToDeg();

            b_cut_OA = false;
            if (b_cut_CA)
            {
                d_OA = (TMath::RadToDeg()*lv_miss.Angle(v_reco));
                b_cut_OA = (d_OA < OACut);
            }

            //////////////////////////////////////////////////
            // Recoil detection stuff
            //////////////////////////////////////////////////
            if (GHistBGSub::IsPrompt(d_subt_tm))
            {
                if (b_cut_CA) PiP_OA->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    if (b_cut_CA) PiP_OA_Cut->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                }
            }
            else if (GHistBGSub::IsRandom(d_subt_tm))
            {
                if (b_cut_CA) PiP_OA_R->Fill(i_tagg_ch, d_miss_th, d_OA,event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    if (b_cut_CA) PiP_OA_Cut_R->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                }
            }

            //////////////////////////////////////////////////
            // State selection
            //////////////////////////////////////////////////
            if (!b_hel) b_fill_Ph_0 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);
            else b_fill_Ph_1 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);

            // Exclusive event, without recoil detection
            if (b_NE)
            {
                PiP_Tm_NE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_NE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_NE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_NE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_NE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_NE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_NE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_NE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_NE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, without recoil detection (add all of the now failed recoil ones in)
            else if (b_NI || !b_cut_OA)
            {
                PiP_Tm_NI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_NI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_NI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_NI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_NI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_NI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_NI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_NI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_NI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // From here on the opening angle must have been satisfied
            // Exclusive event, with recoil detected in CB
            else if (b_CE)
            {
                PiP_Tm_CE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_CE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_CE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_CE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_CE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_CE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_CE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_CE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_CE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in CB
            else if (b_CI)
            {
                PiP_Tm_CI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_CI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_CI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_CI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_CI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_CI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_CI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_CI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_CI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in MWPC
            else if (b_WE)
            {
                PiP_Tm_WE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_WE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_WE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_WE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_WE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_WE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_WE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_WE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_WE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in MWPC
            else if (b_WI)
            {
                PiP_Tm_WI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_WI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_WI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_WI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_WI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_WI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_WI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_WI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_WI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in TAPS
            else if (b_TE)
            {
                PiP_Tm_TE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_TE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_TE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_TE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_TE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_TE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_TE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_TE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_TE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in TAPS
            else if (b_TI)
            {
                PiP_Tm_TI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_TI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_TI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_TI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_TI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) PiP_MM_TI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) PiP_MM_TI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) PiP_Ph_TI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) PiP_Ph_TI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
        }

        //////////////////////////////////////////////////
        // Good Compton event
        //////////////////////////////////////////////////
        if(b_comp)
        {
            d_subt_tm = d_tagg_tm - d_part_tm;
            lv_miss = lv_beam + lv_targ - lv_part;
            d_miss_ma = lv_miss.M()-lv_targ.M();
            d_miss_th = lv_miss.Theta()*TMath::RadToDeg();

            /*
            lv_ptot = lv_beam + lv_targ;
            v_lab_cm = -lv_ptot.BoostVector();
            lv_part_cm = lv_part;
            lv_part_cm.Boost(v_lab_cm);
            d_part_th_cm = lv_part_cm.Theta()*TMath::RadToDeg();
            d_part_ph_cm = lv_part_cm.Phi()*TMath::RadToDeg();
            */

            b_cut_OA = false;
            if (b_cut_CA)
            {
                d_OA = (TMath::RadToDeg()*lv_miss.Angle(v_reco));
                b_cut_OA = (d_OA < OACut);
            }

            //////////////////////////////////////////////////
            // Recoil detection stuff
            //////////////////////////////////////////////////
            if (GHistBGSub::IsPrompt(d_subt_tm))
            {
                if (b_cut_CA) Comp_OA->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    if (b_cut_CA) Comp_OA_Cut->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                    if (b_cut_OA)
                    {
                        if (d_reco_en > 0)
                        {
                            Comp_CS_MM->Fill(d_part_en, d_part_th, i_part_sz, event_weight);
                            Reco_CS_MM->Fill(d_reco_en, d_reco_th, i_reco_sz, event_weight);
                        }
                    }
                }
            }
            else if (GHistBGSub::IsRandom(d_subt_tm))
            {
                if (b_cut_CA) Comp_OA_R->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                if (TMath::Abs(d_miss_ma) < 15)
                {
                    if (b_cut_CA) Comp_OA_Cut_R->Fill(i_tagg_ch, d_miss_th, d_OA, event_weight);
                    if (b_cut_OA)
                    {
                        if (d_reco_en > 0)
                        {
                            Comp_CS_MM_R->Fill(d_part_en, d_part_th, i_part_sz, event_weight);
                            Reco_CS_MM_R->Fill(d_reco_en, d_reco_th, i_reco_sz, event_weight);
                        }
                    }
                }
            }

            //////////////////////////////////////////////////
            // State selection
            //////////////////////////////////////////////////
            if (!b_hel) b_fill_Ph_0 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);
            else b_fill_Ph_1 = (lv_miss.M() >= MMLoC && lv_miss.M() < MMHiC);

            // Exclusive event, without recoil detection
            if (b_NE)
            {
                Comp_Tm_NE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_NE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_NE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_NE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_NE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_NE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_NE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_NE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_NE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, without recoil detection (add all of the now failed recoil ones in)
            else if (b_NI || !b_cut_OA)
            {
                Comp_Tm_NI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_NI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_NI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_NI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_NI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_NI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_NI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_NI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_NI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // From here on the opening angle must have been satisfied
            // Exclusive event, with recoil detected in CB
            else if (b_CE)
            {
                Comp_Tm_CE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_CE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_CE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_CE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_CE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_CE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_CE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_CE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_CE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in CB
            else if (b_CI)
            {
                Comp_Tm_CI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_CI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_CI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_CI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_CI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_CI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_CI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_CI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_CI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in MWPC
            else if (b_WE)
            {
                Comp_Tm_WE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_WE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_WE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_WE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_WE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_WE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_WE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_WE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_WE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in MWPC
            else if (b_WI)
            {
                Comp_Tm_WI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_WI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_WI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_WI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_WI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_WI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_WI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_WI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_WI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Exclusive event, with recoil detected in TAPS
            else if (b_TE)
            {
                Comp_Tm_TE->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_TE_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_TE_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_TE_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_TE_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_TE_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_TE_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_TE_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_TE_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
            // Inclusive event, with recoil detected in TAPS
            else if (b_TI)
            {
                Comp_Tm_TI->Fill(i_tagg_ch, d_subt_tm, event_weight);
                if (GHistBGSub::IsPrompt(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_TI_0->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_TI_1->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_TI_0->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_TI_1->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
                else if (GHistBGSub::IsRandom(d_subt_tm))
                {
                    if (b_fill_MM_0) Comp_MM_TI_0_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    else if (b_fill_MM_1) Comp_MM_TI_1_R->Fill(i_tagg_ch, d_part_th, d_miss_ma, event_weight);
                    if (b_fill_Ph_0) Comp_Ph_TI_0_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                    else if (b_fill_Ph_1) Comp_Ph_TI_1_R->Fill(i_tagg_ch, d_part_th, d_part_ph, event_weight);
                }
            }
        }
    }
}

//////////////////////////////////////////////////
// Set verbosity (not currently used)
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitVerbosity()
{
    Int_t sc1;
    string config = ReadConfig("Verbosity");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Setting verbosity: " << sc1 << endl << endl;
        verbosity = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Verbosity not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Set charge matching
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitMatchCharge()
{
    Int_t sc1;
    string config = ReadConfig("Match-Charge");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Setting charge matching: " << sc1 << endl << endl;
        match_charge = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Charge matching not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Invariant mass cut limits
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitInvariantMass()
{
    Double_t sc1;
    string config = ReadConfig("Invariant-Mass-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting invariant mass width: +/-" << sc1 << " MeV " << endl << endl;
        IMCut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Invariant mass cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Missing mass cut limits
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitMissingMass()
{
    Double_t sc1, sc2;
    string config = ReadConfig("Missing-Mass-Cut");
    if(sscanf( config.c_str(), "%lf%lf\n", &sc1, &sc2) == 2)
    {
        cout << "Setting missing mass cut: " << sc1 << " to " << sc2 << " MeV " << endl << endl;
        MMLoC = sc1;
        MMHiC = sc2;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Missing mass cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Opening angle cut limit
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitOpeningAngle()
{
    Double_t sc1;
    string config = ReadConfig("Opening-Angle-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting opening angle: " << sc1 << " deg " << endl << endl;
        OACut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Opening angle cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Energy sum threshold
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitEnergySum()
{
    Double_t sc1;
    string config = ReadConfig("Energy-Sum-Cut");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "Setting energy sum: " << sc1 << " MeV " << endl << endl;
        ESCut = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Energy sum cut not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Option to save separate histograms for randoms (and prompts)
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitSaveRandoms()
{
    Int_t sc1;
    string config = ReadConfig("Save-Randoms");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Save random hists: " << sc1 << endl << endl;
        save_randoms = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Save randoms not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Option to perform a simplistic split-off search
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitSplitSearch()
{
    Int_t sc1;
    string config = ReadConfig("Split-Search");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Search for splits: " << sc1 << endl << endl;
        split_search = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Split search not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Option to use pure MWPC tracks
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitPureMWPC()
{
    Int_t sc1;
    string config = ReadConfig("Pure-MWPC");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Use pure MWPC tracks: " << sc1 << endl << endl;
        pure_mwpc = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Pure MWPC tracks not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Option to adjust TAPS detection efficiency
//////////////////////////////////////////////////
Bool_t 	PAnalyze::InitTAPSEff()
{
    Double_t sc1;
    string config = ReadConfig("TAPS-Efficiency");
    if(sscanf( config.c_str(), "%lf\n", &sc1) == 1)
    {
        cout << "TAPS efficiency: " << sc1 << endl << endl;
        taps_eff = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "TAPS efficiency not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

//////////////////////////////////////////////////
// Configure beam polarization
//////////////////////////////////////////////////
Bool_t  PAnalyze::InitBeamPol()
{
    string config;
    Int_t instance = 0;
    Int_t time;
    Double_t meas;

    do
    {
        config = ReadConfig("Beam-Polarization",instance);

        if(sscanf( config.c_str(), "%d %lf\n", &time, &meas) == 2)
        {
            beamPolTime.push_back(time);
            beamPolMeas.push_back(meas);
            cout << "Beam polarization measurement: " << 100*meas << "% at " << time << endl;
            instance++;
        }
        else if(strcmp(config.c_str(), "nokey") != 0)
        {
            cout << "Beam polarization not set correctly" << endl;
            return kFALSE;
        }
    } while (strcmp(config.c_str(), "nokey") != 0);

    if(instance)
    {
        cir_beam = true;
        cout << endl;
    }

    Int_t sc1;
    config = ReadConfig("Linear-Pol");
    if(sscanf( config.c_str(), "%d\n", &sc1) == 1)
    {
        cout << "Use linear polarization: " << sc1 << endl << endl;
        lin_beam = sc1;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Linear polarization not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;
}

//////////////////////////////////////////////////
// Configure target polarization
//////////////////////////////////////////////////
Bool_t  PAnalyze::InitTargPol()
{
    string config;
    Int_t instance = 0;
    Int_t time1, time2;
    Double_t meas1, meas2;

    do
    {
        config = ReadConfig("Target-Polarization",instance);

        if(sscanf( config.c_str(), "%d %lf %d %lf\n", &time1, &meas1, &time2, &meas2) == 4)
        {
            targPolTime.push_back(time1);
            targPolMeas.push_back(meas1);
            targPolTime.push_back(time2);
            targPolMeas.push_back(meas2);
            cout << "Target polarization measurement: " << endl;
            cout << "\tInitial: " << 100*meas1 << "% at " << time1 << endl;
            cout << "\tFinal:   " << 100*meas2 << "% at " << time2 << endl;
            instance++;
        }
        else if(strcmp(config.c_str(), "nokey") != 0)
        {
            cout << "Target polarization not set correctly" << endl;
            return kFALSE;
        }
    } while (strcmp(config.c_str(), "nokey") != 0);

    if(instance) cout << endl;

    return kTRUE;
}

//////////////////////////////////////////////////
// Set a different mass for a lorentz vector
//////////////////////////////////////////////////
TLorentzVector  PAnalyze::AdjustMass(TLorentzVector lv, Double_t mass)
{
    Double_t th = lv.Theta();
    Double_t ph = lv.Phi();

    Double_t E 	= lv.E() - lv.M() + mass;
    Double_t P 	= TMath::Sqrt(E*E - mass*mass);
    Double_t Px = P* sin(th)*cos(ph);
    Double_t Py = P* sin(th)*sin(ph);
    Double_t Pz = P* cos(th);

    return TLorentzVector(Px, Py, Pz, E);
}

//////////////////////////////////////////////////
// Calculate the circular photon polarization
//////////////////////////////////////////////////
Double_t  PAnalyze::CalcCircBeamPol(Double_t E_e, Double_t P_e, Double_t E_g)
{
    Double_t P_g = P_e*(((4*E_g*E_e)-(E_g*E_g))/((4*E_e*E_e)-(4*E_g*E_e)+(3*E_g*E_g)));

    return P_g;
}

//////////////////////////////////////////////////
// Individual scaler read stuff
//////////////////////////////////////////////////
void	PAnalyze::ProcessScalerRead()
{
    if (lin_beam)
    {
        TH1D *TaggerPreScal = (TH1D*)TaggerAccScal->Clone("TaggerPreScal");
        TH1D *LiveTimePre = (TH1D*)LiveTimeScal->Clone("LiveTimePre");

        PPhysics::ProcessScalerRead();

        TH1D *TaggerCurScal = (TH1D*)TaggerAccScal->Clone("TaggerCurScal");
        TaggerCurScal->Add(TaggerPreScal,-1);
        TH1D *LiveTimeCur = (TH1D*)LiveTimeScal->Clone("LiveTimeCur");
        LiveTimeCur->Add(LiveTimePre,-1);

        //GoosyTagger(TaggerCurScal);
        //GoosyVuprom(TaggerCurScal);
        //GoosyNewFPD(TaggerCurScal);
        GoosyNewFPDRecabled(TaggerCurScal);

        Double_t livetime = 1;
        if (LiveTimeScal->GetBinContent(1) > 0) livetime = ((LiveTimeScal->GetBinContent(2))/(LiveTimeScal->GetBinContent(1)));
        Double_t counts;
        Double_t photonPol = 1;
        for (Int_t i=0; i<352; i++)
        {
            counts = TaggerCurScal->GetBinContent(i+1);
            photonPol = GetLinpol()->GetPolarizationDegree(i);
            CorrTaggScal->SetBinContent(i+1,((CorrTaggScal->GetBinContent(i+1))+(livetime*counts)));
            PolarizeScal->SetBinContent(i+1,((PolarizeScal->GetBinContent(i+1))+(livetime*targPol*photonPol*counts)));
        }

        delete TaggerPreScal;
        delete LiveTimePre;
        delete TaggerCurScal;
        delete LiveTimeCur;
    }
    else PPhysics::ProcessScalerRead();
}

//////////////////////////////////////////////////
// Process histograms, prompt/random subtraction
//////////////////////////////////////////////////
Bool_t	PAnalyze::Write()
{
    if (!lin_beam)
    {
        Double_t livetime = 1;
        if (LiveTimeScal->GetBinContent(1) > 0) livetime = ((LiveTimeScal->GetBinContent(2))/(LiveTimeScal->GetBinContent(1)));
        Double_t counts;
        Double_t photonPol = 1;
        for (Int_t i=0; i<352; i++)
        {
            counts = TaggerAccScal->GetBinContent(i+1);
            if(beamPol < 1) photonPol = CalcCircBeamPol(450,beamPol,GetSetupParameters()->GetTaggerPhotonEnergy(i));
            CorrTaggScal->SetBinContent(i+1,(livetime*counts));
            PolarizeScal->SetBinContent(i+1,(livetime*targPol*photonPol*counts));
        }
    }

    Double_t ratio = GHistBGSub::GetBackgroundSumbtractionFactor();

    if(save_randoms && !(IsMCFile()))
    {
        TH1D *Tagg_0_P = (TH1D*)Tagg_0->Clone("Tagg_0_P");
        TH1D *Tagg_1_P = (TH1D*)Tagg_1->Clone("Tagg_1_P");

        TH2D *Inc_0_P = (TH2D*)Inc_0->Clone("Inc_0_P");
        TH2D *Inc_1_P = (TH2D*)Inc_1->Clone("Inc_1_P");

        TH3D *Pi0_OA_P = (TH3D*)Pi0_OA->Clone("Pi0_OA_P");
        TH3D *Pi0_OA_Cut_P = (TH3D*)Pi0_OA_Cut->Clone("Pi0_OA_Cut_P");

        TH3D *Pi0_MM_NE_0_P = (TH3D*)Pi0_MM_NE_0->Clone("Pi0_MM_NE_0_P");
        TH3D *Pi0_MM_NE_1_P = (TH3D*)Pi0_MM_NE_1->Clone("Pi0_MM_NE_1_P");
        TH3D *Pi0_Ph_NE_0_P = (TH3D*)Pi0_Ph_NE_0->Clone("Pi0_Ph_NE_0_P");
        TH3D *Pi0_Ph_NE_1_P = (TH3D*)Pi0_Ph_NE_1->Clone("Pi0_Ph_NE_1_P");

        TH3D *Pi0_MM_NI_0_P = (TH3D*)Pi0_MM_NI_0->Clone("Pi0_MM_NI_0_P");
        TH3D *Pi0_MM_NI_1_P = (TH3D*)Pi0_MM_NI_1->Clone("Pi0_MM_NI_1_P");
        TH3D *Pi0_Ph_NI_0_P = (TH3D*)Pi0_Ph_NI_0->Clone("Pi0_Ph_NI_0_P");
        TH3D *Pi0_Ph_NI_1_P = (TH3D*)Pi0_Ph_NI_1->Clone("Pi0_Ph_NI_1_P");

        TH3D *Pi0_MM_CE_0_P = (TH3D*)Pi0_MM_CE_0->Clone("Pi0_MM_CE_0_P");
        TH3D *Pi0_MM_CE_1_P = (TH3D*)Pi0_MM_CE_1->Clone("Pi0_MM_CE_1_P");
        TH3D *Pi0_Ph_CE_0_P = (TH3D*)Pi0_Ph_CE_0->Clone("Pi0_Ph_CE_0_P");
        TH3D *Pi0_Ph_CE_1_P = (TH3D*)Pi0_Ph_CE_1->Clone("Pi0_Ph_CE_1_P");

        TH3D *Pi0_MM_CI_0_P = (TH3D*)Pi0_MM_CI_0->Clone("Pi0_MM_CI_0_P");
        TH3D *Pi0_MM_CI_1_P = (TH3D*)Pi0_MM_CI_1->Clone("Pi0_MM_CI_1_P");
        TH3D *Pi0_Ph_CI_0_P = (TH3D*)Pi0_Ph_CI_0->Clone("Pi0_Ph_CI_0_P");
        TH3D *Pi0_Ph_CI_1_P = (TH3D*)Pi0_Ph_CI_1->Clone("Pi0_Ph_CI_1_P");

        TH3D *Pi0_MM_WE_0_P = (TH3D*)Pi0_MM_WE_0->Clone("Pi0_MM_WE_0_P");
        TH3D *Pi0_MM_WE_1_P = (TH3D*)Pi0_MM_WE_1->Clone("Pi0_MM_WE_1_P");
        TH3D *Pi0_Ph_WE_0_P = (TH3D*)Pi0_Ph_WE_0->Clone("Pi0_Ph_WE_0_P");
        TH3D *Pi0_Ph_WE_1_P = (TH3D*)Pi0_Ph_WE_1->Clone("Pi0_Ph_WE_1_P");

        TH3D *Pi0_MM_WI_0_P = (TH3D*)Pi0_MM_WI_0->Clone("Pi0_MM_WI_0_P");
        TH3D *Pi0_MM_WI_1_P = (TH3D*)Pi0_MM_WI_1->Clone("Pi0_MM_WI_1_P");
        TH3D *Pi0_Ph_WI_0_P = (TH3D*)Pi0_Ph_WI_0->Clone("Pi0_Ph_WI_0_P");
        TH3D *Pi0_Ph_WI_1_P = (TH3D*)Pi0_Ph_WI_1->Clone("Pi0_Ph_WI_1_P");

        TH3D *Pi0_MM_TE_0_P = (TH3D*)Pi0_MM_TE_0->Clone("Pi0_MM_TE_0_P");
        TH3D *Pi0_MM_TE_1_P = (TH3D*)Pi0_MM_TE_1->Clone("Pi0_MM_TE_1_P");
        TH3D *Pi0_Ph_TE_0_P = (TH3D*)Pi0_Ph_TE_0->Clone("Pi0_Ph_TE_0_P");
        TH3D *Pi0_Ph_TE_1_P = (TH3D*)Pi0_Ph_TE_1->Clone("Pi0_Ph_TE_1_P");

        TH3D *Pi0_MM_TI_0_P = (TH3D*)Pi0_MM_TI_0->Clone("Pi0_MM_TI_0_P");
        TH3D *Pi0_MM_TI_1_P = (TH3D*)Pi0_MM_TI_1->Clone("Pi0_MM_TI_1_P");
        TH3D *Pi0_Ph_TI_0_P = (TH3D*)Pi0_Ph_TI_0->Clone("Pi0_Ph_TI_0_P");
        TH3D *Pi0_Ph_TI_1_P = (TH3D*)Pi0_Ph_TI_1->Clone("Pi0_Ph_TI_1_P");

        TH3D *Pi0_Re_All_P = (TH3D*)Pi0_Re_All->Clone("Pi0_Re_All_P");
        TH3D *Pi0_Re_Det_P = (TH3D*)Pi0_Re_Det->Clone("Pi0_Re_Det_P");
        TH3D *Pi0_Re_Dif_P = (TH3D*)Pi0_Re_Dif->Clone("Pi0_Re_Dif_P");
        TH3D *Pi0_Re_NoE_P = (TH3D*)Pi0_Re_NoE->Clone("Pi0_Re_NoE_P");

        TH3D *PiP_OA_P = (TH3D*)PiP_OA->Clone("PiP_OA_P");
        TH3D *PiP_OA_Cut_P = (TH3D*)PiP_OA_Cut->Clone("PiP_OA_Cut_P");

        TH3D *PiP_MM_NE_0_P = (TH3D*)PiP_MM_NE_0->Clone("PiP_MM_NE_0_P");
        TH3D *PiP_MM_NE_1_P = (TH3D*)PiP_MM_NE_1->Clone("PiP_MM_NE_1_P");
        TH3D *PiP_Ph_NE_0_P = (TH3D*)PiP_Ph_NE_0->Clone("PiP_Ph_NE_0_P");
        TH3D *PiP_Ph_NE_1_P = (TH3D*)PiP_Ph_NE_1->Clone("PiP_Ph_NE_1_P");

        TH3D *PiP_MM_NI_0_P = (TH3D*)PiP_MM_NI_0->Clone("PiP_MM_NI_0_P");
        TH3D *PiP_MM_NI_1_P = (TH3D*)PiP_MM_NI_1->Clone("PiP_MM_NI_1_P");
        TH3D *PiP_Ph_NI_0_P = (TH3D*)PiP_Ph_NI_0->Clone("PiP_Ph_NI_0_P");
        TH3D *PiP_Ph_NI_1_P = (TH3D*)PiP_Ph_NI_1->Clone("PiP_Ph_NI_1_P");

        TH3D *PiP_MM_CE_0_P = (TH3D*)PiP_MM_CE_0->Clone("PiP_MM_CE_0_P");
        TH3D *PiP_MM_CE_1_P = (TH3D*)PiP_MM_CE_1->Clone("PiP_MM_CE_1_P");
        TH3D *PiP_Ph_CE_0_P = (TH3D*)PiP_Ph_CE_0->Clone("PiP_Ph_CE_0_P");
        TH3D *PiP_Ph_CE_1_P = (TH3D*)PiP_Ph_CE_1->Clone("PiP_Ph_CE_1_P");

        TH3D *PiP_MM_CI_0_P = (TH3D*)PiP_MM_CI_0->Clone("PiP_MM_CI_0_P");
        TH3D *PiP_MM_CI_1_P = (TH3D*)PiP_MM_CI_1->Clone("PiP_MM_CI_1_P");
        TH3D *PiP_Ph_CI_0_P = (TH3D*)PiP_Ph_CI_0->Clone("PiP_Ph_CI_0_P");
        TH3D *PiP_Ph_CI_1_P = (TH3D*)PiP_Ph_CI_1->Clone("PiP_Ph_CI_1_P");

        TH3D *PiP_MM_WE_0_P = (TH3D*)PiP_MM_WE_0->Clone("PiP_MM_WE_0_P");
        TH3D *PiP_MM_WE_1_P = (TH3D*)PiP_MM_WE_1->Clone("PiP_MM_WE_1_P");
        TH3D *PiP_Ph_WE_0_P = (TH3D*)PiP_Ph_WE_0->Clone("PiP_Ph_WE_0_P");
        TH3D *PiP_Ph_WE_1_P = (TH3D*)PiP_Ph_WE_1->Clone("PiP_Ph_WE_1_P");

        TH3D *PiP_MM_WI_0_P = (TH3D*)PiP_MM_WI_0->Clone("PiP_MM_WI_0_P");
        TH3D *PiP_MM_WI_1_P = (TH3D*)PiP_MM_WI_1->Clone("PiP_MM_WI_1_P");
        TH3D *PiP_Ph_WI_0_P = (TH3D*)PiP_Ph_WI_0->Clone("PiP_Ph_WI_0_P");
        TH3D *PiP_Ph_WI_1_P = (TH3D*)PiP_Ph_WI_1->Clone("PiP_Ph_WI_1_P");

        TH3D *PiP_MM_TE_0_P = (TH3D*)PiP_MM_TE_0->Clone("PiP_MM_TE_0_P");
        TH3D *PiP_MM_TE_1_P = (TH3D*)PiP_MM_TE_1->Clone("PiP_MM_TE_1_P");
        TH3D *PiP_Ph_TE_0_P = (TH3D*)PiP_Ph_TE_0->Clone("PiP_Ph_TE_0_P");
        TH3D *PiP_Ph_TE_1_P = (TH3D*)PiP_Ph_TE_1->Clone("PiP_Ph_TE_1_P");

        TH3D *PiP_MM_TI_0_P = (TH3D*)PiP_MM_TI_0->Clone("PiP_MM_TI_0_P");
        TH3D *PiP_MM_TI_1_P = (TH3D*)PiP_MM_TI_1->Clone("PiP_MM_TI_1_P");
        TH3D *PiP_Ph_TI_0_P = (TH3D*)PiP_Ph_TI_0->Clone("PiP_Ph_TI_0_P");
        TH3D *PiP_Ph_TI_1_P = (TH3D*)PiP_Ph_TI_1->Clone("PiP_Ph_TI_1_P");

        TH3D *Comp_OA_P = (TH3D*)Comp_OA->Clone("Comp_OA_P");
        TH3D *Comp_OA_Cut_P = (TH3D*)Comp_OA_Cut->Clone("Comp_OA_Cut_P");

        TH3D *Comp_MM_NE_0_P = (TH3D*)Comp_MM_NE_0->Clone("Comp_MM_NE_0_P");
        TH3D *Comp_MM_NE_1_P = (TH3D*)Comp_MM_NE_1->Clone("Comp_MM_NE_1_P");
        TH3D *Comp_Ph_NE_0_P = (TH3D*)Comp_Ph_NE_0->Clone("Comp_Ph_NE_0_P");
        TH3D *Comp_Ph_NE_1_P = (TH3D*)Comp_Ph_NE_1->Clone("Comp_Ph_NE_1_P");

        TH3D *Comp_MM_NI_0_P = (TH3D*)Comp_MM_NI_0->Clone("Comp_MM_NI_0_P");
        TH3D *Comp_MM_NI_1_P = (TH3D*)Comp_MM_NI_1->Clone("Comp_MM_NI_1_P");
        TH3D *Comp_Ph_NI_0_P = (TH3D*)Comp_Ph_NI_0->Clone("Comp_Ph_NI_0_P");
        TH3D *Comp_Ph_NI_1_P = (TH3D*)Comp_Ph_NI_1->Clone("Comp_Ph_NI_1_P");

        TH3D *Comp_MM_CE_0_P = (TH3D*)Comp_MM_CE_0->Clone("Comp_MM_CE_0_P");
        TH3D *Comp_MM_CE_1_P = (TH3D*)Comp_MM_CE_1->Clone("Comp_MM_CE_1_P");
        TH3D *Comp_Ph_CE_0_P = (TH3D*)Comp_Ph_CE_0->Clone("Comp_Ph_CE_0_P");
        TH3D *Comp_Ph_CE_1_P = (TH3D*)Comp_Ph_CE_1->Clone("Comp_Ph_CE_1_P");

        TH3D *Comp_MM_CI_0_P = (TH3D*)Comp_MM_CI_0->Clone("Comp_MM_CI_0_P");
        TH3D *Comp_MM_CI_1_P = (TH3D*)Comp_MM_CI_1->Clone("Comp_MM_CI_1_P");
        TH3D *Comp_Ph_CI_0_P = (TH3D*)Comp_Ph_CI_0->Clone("Comp_Ph_CI_0_P");
        TH3D *Comp_Ph_CI_1_P = (TH3D*)Comp_Ph_CI_1->Clone("Comp_Ph_CI_1_P");

        TH3D *Comp_MM_WE_0_P = (TH3D*)Comp_MM_WE_0->Clone("Comp_MM_WE_0_P");
        TH3D *Comp_MM_WE_1_P = (TH3D*)Comp_MM_WE_1->Clone("Comp_MM_WE_1_P");
        TH3D *Comp_Ph_WE_0_P = (TH3D*)Comp_Ph_WE_0->Clone("Comp_Ph_WE_0_P");
        TH3D *Comp_Ph_WE_1_P = (TH3D*)Comp_Ph_WE_1->Clone("Comp_Ph_WE_1_P");

        TH3D *Comp_MM_WI_0_P = (TH3D*)Comp_MM_WI_0->Clone("Comp_MM_WI_0_P");
        TH3D *Comp_MM_WI_1_P = (TH3D*)Comp_MM_WI_1->Clone("Comp_MM_WI_1_P");
        TH3D *Comp_Ph_WI_0_P = (TH3D*)Comp_Ph_WI_0->Clone("Comp_Ph_WI_0_P");
        TH3D *Comp_Ph_WI_1_P = (TH3D*)Comp_Ph_WI_1->Clone("Comp_Ph_WI_1_P");

        TH3D *Comp_MM_TE_0_P = (TH3D*)Comp_MM_TE_0->Clone("Comp_MM_TE_0_P");
        TH3D *Comp_MM_TE_1_P = (TH3D*)Comp_MM_TE_1->Clone("Comp_MM_TE_1_P");
        TH3D *Comp_Ph_TE_0_P = (TH3D*)Comp_Ph_TE_0->Clone("Comp_Ph_TE_0_P");
        TH3D *Comp_Ph_TE_1_P = (TH3D*)Comp_Ph_TE_1->Clone("Comp_Ph_TE_1_P");

        TH3D *Comp_MM_TI_0_P = (TH3D*)Comp_MM_TI_0->Clone("Comp_MM_TI_0_P");
        TH3D *Comp_MM_TI_1_P = (TH3D*)Comp_MM_TI_1->Clone("Comp_MM_TI_1_P");
        TH3D *Comp_Ph_TI_0_P = (TH3D*)Comp_Ph_TI_0->Clone("Comp_Ph_TI_0_P");
        TH3D *Comp_Ph_TI_1_P = (TH3D*)Comp_Ph_TI_1->Clone("Comp_Ph_TI_1_P");

        TH3D *Comp_CS_MM_P = (TH3D*)Comp_CS_MM->Clone("Comp_CS_MM_P");
        TH3D *Reco_CS_MM_P = (TH3D*)Reco_CS_MM->Clone("Reco_CS_MM_P");

        Tagg_0_P->Write();
        Tagg_1_P->Write();

        Inc_0_P->Write();
        Inc_1_P->Write();

        Pi0_OA_P->Write();
        Pi0_OA_Cut_P->Write();

        Pi0_MM_NE_0_P->Write();
        Pi0_MM_NI_0_P->Write();
        Pi0_MM_CE_0_P->Write();
        Pi0_MM_CI_0_P->Write();
        Pi0_MM_WE_0_P->Write();
        Pi0_MM_WI_0_P->Write();
        Pi0_MM_TE_0_P->Write();
        Pi0_MM_TI_0_P->Write();

        Pi0_Re_All_P->Write();
        Pi0_Re_Det_P->Write();
        Pi0_Re_Dif_P->Write();
        Pi0_Re_NoE_P->Write();

        PiP_OA_P->Write();
        PiP_OA_Cut_P->Write();

        PiP_MM_NE_0_P->Write();
        PiP_MM_NI_0_P->Write();
        PiP_MM_CE_0_P->Write();
        PiP_MM_CI_0_P->Write();
        PiP_MM_WE_0_P->Write();
        PiP_MM_WI_0_P->Write();
        PiP_MM_TE_0_P->Write();
        PiP_MM_TI_0_P->Write();

        Comp_OA_P->Write();
        Comp_OA_Cut_P->Write();

        Comp_MM_NE_0_P->Write();
        Comp_MM_NI_0_P->Write();
        Comp_MM_CE_0_P->Write();
        Comp_MM_CI_0_P->Write();
        Comp_MM_WE_0_P->Write();
        Comp_MM_WI_0_P->Write();
        Comp_MM_TE_0_P->Write();
        Comp_MM_TI_0_P->Write();

        if(cir_beam || lin_beam)
        {
            Pi0_MM_NE_1_P->Write();
            Pi0_Ph_NE_0_P->Write();
            Pi0_Ph_NE_1_P->Write();

            Pi0_MM_NI_1_P->Write();
            Pi0_Ph_NI_0_P->Write();
            Pi0_Ph_NI_1_P->Write();

            Pi0_MM_CE_1_P->Write();
            Pi0_Ph_CE_0_P->Write();
            Pi0_Ph_CE_1_P->Write();

            Pi0_MM_CI_1_P->Write();
            Pi0_Ph_CI_0_P->Write();
            Pi0_Ph_CI_1_P->Write();

            Pi0_MM_WE_1_P->Write();
            Pi0_Ph_WE_0_P->Write();
            Pi0_Ph_WE_1_P->Write();

            Pi0_MM_WI_1_P->Write();
            Pi0_Ph_WI_0_P->Write();
            Pi0_Ph_WI_1_P->Write();

            Pi0_MM_TE_1_P->Write();
            Pi0_Ph_TE_0_P->Write();
            Pi0_Ph_TE_1_P->Write();

            Pi0_MM_TI_1_P->Write();
            Pi0_Ph_TI_0_P->Write();
            Pi0_Ph_TI_1_P->Write();

            PiP_MM_NE_1_P->Write();
            PiP_Ph_NE_0_P->Write();
            PiP_Ph_NE_1_P->Write();

            PiP_MM_NI_1_P->Write();
            PiP_Ph_NI_0_P->Write();
            PiP_Ph_NI_1_P->Write();

            PiP_MM_CE_1_P->Write();
            PiP_Ph_CE_0_P->Write();
            PiP_Ph_CE_1_P->Write();

            PiP_MM_CI_1_P->Write();
            PiP_Ph_CI_0_P->Write();
            PiP_Ph_CI_1_P->Write();

            PiP_MM_WE_1_P->Write();
            PiP_Ph_WE_0_P->Write();
            PiP_Ph_WE_1_P->Write();

            PiP_MM_WI_1_P->Write();
            PiP_Ph_WI_0_P->Write();
            PiP_Ph_WI_1_P->Write();

            PiP_MM_TE_1_P->Write();
            PiP_Ph_TE_0_P->Write();
            PiP_Ph_TE_1_P->Write();

            PiP_MM_TI_1_P->Write();
            PiP_Ph_TI_0_P->Write();
            PiP_Ph_TI_1_P->Write();

            Comp_MM_NE_1_P->Write();
            Comp_Ph_NE_0_P->Write();
            Comp_Ph_NE_1_P->Write();

            Comp_MM_NI_1_P->Write();
            Comp_Ph_NI_0_P->Write();
            Comp_Ph_NI_1_P->Write();

            Comp_MM_CE_1_P->Write();
            Comp_Ph_CE_0_P->Write();
            Comp_Ph_CE_1_P->Write();

            Comp_MM_CI_1_P->Write();
            Comp_Ph_CI_0_P->Write();
            Comp_Ph_CI_1_P->Write();

            Comp_MM_WE_1_P->Write();
            Comp_Ph_WE_0_P->Write();
            Comp_Ph_WE_1_P->Write();

            Comp_MM_WI_1_P->Write();
            Comp_Ph_WI_0_P->Write();
            Comp_Ph_WI_1_P->Write();

            Comp_MM_TE_1_P->Write();
            Comp_Ph_TE_0_P->Write();
            Comp_Ph_TE_1_P->Write();

            Comp_MM_TI_1_P->Write();
            Comp_Ph_TI_0_P->Write();
            Comp_Ph_TI_1_P->Write();
        }

        Comp_CS_MM_P->Write();
        Reco_CS_MM_P->Write();

        delete Tagg_0_P;
        delete Tagg_1_P;

        delete Inc_0_P;
        delete Inc_1_P;

        delete Pi0_OA_P;
        delete Pi0_OA_Cut_P;

        delete Pi0_MM_NE_0_P;
        delete Pi0_MM_NE_1_P;
        delete Pi0_Ph_NE_0_P;
        delete Pi0_Ph_NE_1_P;

        delete Pi0_MM_NI_0_P;
        delete Pi0_MM_NI_1_P;
        delete Pi0_Ph_NI_0_P;
        delete Pi0_Ph_NI_1_P;

        delete Pi0_MM_CE_0_P;
        delete Pi0_MM_CE_1_P;
        delete Pi0_Ph_CE_0_P;
        delete Pi0_Ph_CE_1_P;

        delete Pi0_MM_CI_0_P;
        delete Pi0_MM_CI_1_P;
        delete Pi0_Ph_CI_0_P;
        delete Pi0_Ph_CI_1_P;

        delete Pi0_MM_WE_0_P;
        delete Pi0_MM_WE_1_P;
        delete Pi0_Ph_WE_0_P;
        delete Pi0_Ph_WE_1_P;

        delete Pi0_MM_WI_0_P;
        delete Pi0_MM_WI_1_P;
        delete Pi0_Ph_WI_0_P;
        delete Pi0_Ph_WI_1_P;

        delete Pi0_MM_TE_0_P;
        delete Pi0_MM_TE_1_P;
        delete Pi0_Ph_TE_0_P;
        delete Pi0_Ph_TE_1_P;

        delete Pi0_MM_TI_0_P;
        delete Pi0_MM_TI_1_P;
        delete Pi0_Ph_TI_0_P;
        delete Pi0_Ph_TI_1_P;

        delete Pi0_Re_All_P;
        delete Pi0_Re_Det_P;
        delete Pi0_Re_Dif_P;
        delete Pi0_Re_NoE_P;

        delete PiP_OA_P;
        delete PiP_OA_Cut_P;

        delete PiP_MM_NE_0_P;
        delete PiP_MM_NE_1_P;
        delete PiP_Ph_NE_0_P;
        delete PiP_Ph_NE_1_P;

        delete PiP_MM_NI_0_P;
        delete PiP_MM_NI_1_P;
        delete PiP_Ph_NI_0_P;
        delete PiP_Ph_NI_1_P;

        delete PiP_MM_CE_0_P;
        delete PiP_MM_CE_1_P;
        delete PiP_Ph_CE_0_P;
        delete PiP_Ph_CE_1_P;

        delete PiP_MM_CI_0_P;
        delete PiP_MM_CI_1_P;
        delete PiP_Ph_CI_0_P;
        delete PiP_Ph_CI_1_P;

        delete PiP_MM_WE_0_P;
        delete PiP_MM_WE_1_P;
        delete PiP_Ph_WE_0_P;
        delete PiP_Ph_WE_1_P;

        delete PiP_MM_WI_0_P;
        delete PiP_MM_WI_1_P;
        delete PiP_Ph_WI_0_P;
        delete PiP_Ph_WI_1_P;

        delete PiP_MM_TE_0_P;
        delete PiP_MM_TE_1_P;
        delete PiP_Ph_TE_0_P;
        delete PiP_Ph_TE_1_P;

        delete PiP_MM_TI_0_P;
        delete PiP_MM_TI_1_P;
        delete PiP_Ph_TI_0_P;
        delete PiP_Ph_TI_1_P;

        delete Comp_OA_P;
        delete Comp_OA_Cut_P;

        delete Comp_MM_NE_0_P;
        delete Comp_MM_NE_1_P;
        delete Comp_Ph_NE_0_P;
        delete Comp_Ph_NE_1_P;

        delete Comp_MM_NI_0_P;
        delete Comp_MM_NI_1_P;
        delete Comp_Ph_NI_0_P;
        delete Comp_Ph_NI_1_P;

        delete Comp_MM_CE_0_P;
        delete Comp_MM_CE_1_P;
        delete Comp_Ph_CE_0_P;
        delete Comp_Ph_CE_1_P;

        delete Comp_MM_CI_0_P;
        delete Comp_MM_CI_1_P;
        delete Comp_Ph_CI_0_P;
        delete Comp_Ph_CI_1_P;

        delete Comp_MM_WE_0_P;
        delete Comp_MM_WE_1_P;
        delete Comp_Ph_WE_0_P;
        delete Comp_Ph_WE_1_P;

        delete Comp_MM_WI_0_P;
        delete Comp_MM_WI_1_P;
        delete Comp_Ph_WI_0_P;
        delete Comp_Ph_WI_1_P;

        delete Comp_MM_TE_0_P;
        delete Comp_MM_TE_1_P;
        delete Comp_Ph_TE_0_P;
        delete Comp_Ph_TE_1_P;

        delete Comp_MM_TI_0_P;
        delete Comp_MM_TI_1_P;
        delete Comp_Ph_TI_0_P;
        delete Comp_Ph_TI_1_P;

        delete Comp_CS_MM_P;
        delete Reco_CS_MM_P;

    }
    if(!(IsMCFile()))
    {
        Tagg_0->Sumw2();
        Tagg_0->Add(Tagg_0_R,-ratio);
        Tagg_1->Sumw2();
        Tagg_1->Add(Tagg_1_R,-ratio);

        Inc_0->Sumw2();
        Inc_0->Add(Inc_0_R,-ratio);
        Inc_1->Sumw2();
        Inc_1->Add(Inc_1_R,-ratio);

        Pi0_OA->Sumw2();
        Pi0_OA->Add(Pi0_OA_R,-ratio);
        Pi0_OA_Cut->Sumw2();
        Pi0_OA_Cut->Add(Pi0_OA_Cut_R,-ratio);

        Pi0_MM_NE_0->Sumw2();
        Pi0_MM_NE_0->Add(Pi0_MM_NE_0_R,-ratio);
        Pi0_MM_NE_1->Sumw2();
        Pi0_MM_NE_1->Add(Pi0_MM_NE_1_R,-ratio);
        Pi0_Ph_NE_0->Sumw2();
        Pi0_Ph_NE_0->Add(Pi0_Ph_NE_0_R,-ratio);
        Pi0_Ph_NE_1->Sumw2();
        Pi0_Ph_NE_1->Add(Pi0_Ph_NE_1_R,-ratio);

        Pi0_MM_NI_0->Sumw2();
        Pi0_MM_NI_0->Add(Pi0_MM_NI_0_R,-ratio);
        Pi0_MM_NI_1->Sumw2();
        Pi0_MM_NI_1->Add(Pi0_MM_NI_1_R,-ratio);
        Pi0_Ph_NI_0->Sumw2();
        Pi0_Ph_NI_0->Add(Pi0_Ph_NI_0_R,-ratio);
        Pi0_Ph_NI_1->Sumw2();
        Pi0_Ph_NI_1->Add(Pi0_Ph_NI_1_R,-ratio);

        Pi0_MM_CE_0->Sumw2();
        Pi0_MM_CE_0->Add(Pi0_MM_CE_0_R,-ratio);
        Pi0_MM_CE_1->Sumw2();
        Pi0_MM_CE_1->Add(Pi0_MM_CE_1_R,-ratio);
        Pi0_Ph_CE_0->Sumw2();
        Pi0_Ph_CE_0->Add(Pi0_Ph_CE_0_R,-ratio);
        Pi0_Ph_CE_1->Sumw2();
        Pi0_Ph_CE_1->Add(Pi0_Ph_CE_1_R,-ratio);

        Pi0_MM_CI_0->Sumw2();
        Pi0_MM_CI_0->Add(Pi0_MM_CI_0_R,-ratio);
        Pi0_MM_CI_1->Sumw2();
        Pi0_MM_CI_1->Add(Pi0_MM_CI_1_R,-ratio);
        Pi0_Ph_CI_0->Sumw2();
        Pi0_Ph_CI_0->Add(Pi0_Ph_CI_0_R,-ratio);
        Pi0_Ph_CI_1->Sumw2();
        Pi0_Ph_CI_1->Add(Pi0_Ph_CI_1_R,-ratio);

        Pi0_MM_WE_0->Sumw2();
        Pi0_MM_WE_0->Add(Pi0_MM_WE_0_R,-ratio);
        Pi0_MM_WE_1->Sumw2();
        Pi0_MM_WE_1->Add(Pi0_MM_WE_1_R,-ratio);
        Pi0_Ph_WE_0->Sumw2();
        Pi0_Ph_WE_0->Add(Pi0_Ph_WE_0_R,-ratio);
        Pi0_Ph_WE_1->Sumw2();
        Pi0_Ph_WE_1->Add(Pi0_Ph_WE_1_R,-ratio);

        Pi0_MM_WI_0->Sumw2();
        Pi0_MM_WI_0->Add(Pi0_MM_WI_0_R,-ratio);
        Pi0_MM_WI_1->Sumw2();
        Pi0_MM_WI_1->Add(Pi0_MM_WI_1_R,-ratio);
        Pi0_Ph_WI_0->Sumw2();
        Pi0_Ph_WI_0->Add(Pi0_Ph_WI_0_R,-ratio);
        Pi0_Ph_WI_1->Sumw2();
        Pi0_Ph_WI_1->Add(Pi0_Ph_WI_1_R,-ratio);

        Pi0_MM_TE_0->Sumw2();
        Pi0_MM_TE_0->Add(Pi0_MM_TE_0_R,-ratio);
        Pi0_MM_TE_1->Sumw2();
        Pi0_MM_TE_1->Add(Pi0_MM_TE_1_R,-ratio);
        Pi0_Ph_TE_0->Sumw2();
        Pi0_Ph_TE_0->Add(Pi0_Ph_TE_0_R,-ratio);
        Pi0_Ph_TE_1->Sumw2();
        Pi0_Ph_TE_1->Add(Pi0_Ph_TE_1_R,-ratio);

        Pi0_MM_TI_0->Sumw2();
        Pi0_MM_TI_0->Add(Pi0_MM_TI_0_R,-ratio);
        Pi0_MM_TI_1->Sumw2();
        Pi0_MM_TI_1->Add(Pi0_MM_TI_1_R,-ratio);
        Pi0_Ph_TI_0->Sumw2();
        Pi0_Ph_TI_0->Add(Pi0_Ph_TI_0_R,-ratio);
        Pi0_Ph_TI_1->Sumw2();
        Pi0_Ph_TI_1->Add(Pi0_Ph_TI_1_R,-ratio);

        Pi0_Re_All->Sumw2();
        Pi0_Re_All->Add(Pi0_Re_All_R,-ratio);
        Pi0_Re_Det->Sumw2();
        Pi0_Re_Det->Add(Pi0_Re_Det_R,-ratio);
        Pi0_Re_Dif->Sumw2();
        Pi0_Re_Dif->Add(Pi0_Re_Dif_R,-ratio);
        Pi0_Re_NoE->Sumw2();
        Pi0_Re_NoE->Add(Pi0_Re_NoE_R,-ratio);

        PiP_OA->Sumw2();
        PiP_OA->Add(PiP_OA_R,-ratio);
        PiP_OA_Cut->Sumw2();
        PiP_OA_Cut->Add(PiP_OA_Cut_R,-ratio);

        PiP_MM_NE_0->Sumw2();
        PiP_MM_NE_0->Add(PiP_MM_NE_0_R,-ratio);
        PiP_MM_NE_1->Sumw2();
        PiP_MM_NE_1->Add(PiP_MM_NE_1_R,-ratio);
        PiP_Ph_NE_0->Sumw2();
        PiP_Ph_NE_0->Add(PiP_Ph_NE_0_R,-ratio);
        PiP_Ph_NE_1->Sumw2();
        PiP_Ph_NE_1->Add(PiP_Ph_NE_1_R,-ratio);

        PiP_MM_NI_0->Sumw2();
        PiP_MM_NI_0->Add(PiP_MM_NI_0_R,-ratio);
        PiP_MM_NI_1->Sumw2();
        PiP_MM_NI_1->Add(PiP_MM_NI_1_R,-ratio);
        PiP_Ph_NI_0->Sumw2();
        PiP_Ph_NI_0->Add(PiP_Ph_NI_0_R,-ratio);
        PiP_Ph_NI_1->Sumw2();
        PiP_Ph_NI_1->Add(PiP_Ph_NI_1_R,-ratio);

        PiP_MM_CE_0->Sumw2();
        PiP_MM_CE_0->Add(PiP_MM_CE_0_R,-ratio);
        PiP_MM_CE_1->Sumw2();
        PiP_MM_CE_1->Add(PiP_MM_CE_1_R,-ratio);
        PiP_Ph_CE_0->Sumw2();
        PiP_Ph_CE_0->Add(PiP_Ph_CE_0_R,-ratio);
        PiP_Ph_CE_1->Sumw2();
        PiP_Ph_CE_1->Add(PiP_Ph_CE_1_R,-ratio);

        PiP_MM_CI_0->Sumw2();
        PiP_MM_CI_0->Add(PiP_MM_CI_0_R,-ratio);
        PiP_MM_CI_1->Sumw2();
        PiP_MM_CI_1->Add(PiP_MM_CI_1_R,-ratio);
        PiP_Ph_CI_0->Sumw2();
        PiP_Ph_CI_0->Add(PiP_Ph_CI_0_R,-ratio);
        PiP_Ph_CI_1->Sumw2();
        PiP_Ph_CI_1->Add(PiP_Ph_CI_1_R,-ratio);

        PiP_MM_WE_0->Sumw2();
        PiP_MM_WE_0->Add(PiP_MM_WE_0_R,-ratio);
        PiP_MM_WE_1->Sumw2();
        PiP_MM_WE_1->Add(PiP_MM_WE_1_R,-ratio);
        PiP_Ph_WE_0->Sumw2();
        PiP_Ph_WE_0->Add(PiP_Ph_WE_0_R,-ratio);
        PiP_Ph_WE_1->Sumw2();
        PiP_Ph_WE_1->Add(PiP_Ph_WE_1_R,-ratio);

        PiP_MM_WI_0->Sumw2();
        PiP_MM_WI_0->Add(PiP_MM_WI_0_R,-ratio);
        PiP_MM_WI_1->Sumw2();
        PiP_MM_WI_1->Add(PiP_MM_WI_1_R,-ratio);
        PiP_Ph_WI_0->Sumw2();
        PiP_Ph_WI_0->Add(PiP_Ph_WI_0_R,-ratio);
        PiP_Ph_WI_1->Sumw2();
        PiP_Ph_WI_1->Add(PiP_Ph_WI_1_R,-ratio);

        PiP_MM_TE_0->Sumw2();
        PiP_MM_TE_0->Add(PiP_MM_TE_0_R,-ratio);
        PiP_MM_TE_1->Sumw2();
        PiP_MM_TE_1->Add(PiP_MM_TE_1_R,-ratio);
        PiP_Ph_TE_0->Sumw2();
        PiP_Ph_TE_0->Add(PiP_Ph_TE_0_R,-ratio);
        PiP_Ph_TE_1->Sumw2();
        PiP_Ph_TE_1->Add(PiP_Ph_TE_1_R,-ratio);

        PiP_MM_TI_0->Sumw2();
        PiP_MM_TI_0->Add(PiP_MM_TI_0_R,-ratio);
        PiP_MM_TI_1->Sumw2();
        PiP_MM_TI_1->Add(PiP_MM_TI_1_R,-ratio);
        PiP_Ph_TI_0->Sumw2();
        PiP_Ph_TI_0->Add(PiP_Ph_TI_0_R,-ratio);
        PiP_Ph_TI_1->Sumw2();
        PiP_Ph_TI_1->Add(PiP_Ph_TI_1_R,-ratio);

        Comp_OA->Sumw2();
        Comp_OA->Add(Comp_OA_R,-ratio);
        Comp_OA_Cut->Sumw2();
        Comp_OA_Cut->Add(Comp_OA_Cut_R,-ratio);

        Comp_MM_NE_0->Sumw2();
        Comp_MM_NE_0->Add(Comp_MM_NE_0_R,-ratio);
        Comp_MM_NE_1->Sumw2();
        Comp_MM_NE_1->Add(Comp_MM_NE_1_R,-ratio);
        Comp_Ph_NE_0->Sumw2();
        Comp_Ph_NE_0->Add(Comp_Ph_NE_0_R,-ratio);
        Comp_Ph_NE_1->Sumw2();
        Comp_Ph_NE_1->Add(Comp_Ph_NE_1_R,-ratio);

        Comp_MM_NI_0->Sumw2();
        Comp_MM_NI_0->Add(Comp_MM_NI_0_R,-ratio);
        Comp_MM_NI_1->Sumw2();
        Comp_MM_NI_1->Add(Comp_MM_NI_1_R,-ratio);
        Comp_Ph_NI_0->Sumw2();
        Comp_Ph_NI_0->Add(Comp_Ph_NI_0_R,-ratio);
        Comp_Ph_NI_1->Sumw2();
        Comp_Ph_NI_1->Add(Comp_Ph_NI_1_R,-ratio);

        Comp_MM_CE_0->Sumw2();
        Comp_MM_CE_0->Add(Comp_MM_CE_0_R,-ratio);
        Comp_MM_CE_1->Sumw2();
        Comp_MM_CE_1->Add(Comp_MM_CE_1_R,-ratio);
        Comp_Ph_CE_0->Sumw2();
        Comp_Ph_CE_0->Add(Comp_Ph_CE_0_R,-ratio);
        Comp_Ph_CE_1->Sumw2();
        Comp_Ph_CE_1->Add(Comp_Ph_CE_1_R,-ratio);

        Comp_MM_CI_0->Sumw2();
        Comp_MM_CI_0->Add(Comp_MM_CI_0_R,-ratio);
        Comp_MM_CI_1->Sumw2();
        Comp_MM_CI_1->Add(Comp_MM_CI_1_R,-ratio);
        Comp_Ph_CI_0->Sumw2();
        Comp_Ph_CI_0->Add(Comp_Ph_CI_0_R,-ratio);
        Comp_Ph_CI_1->Sumw2();
        Comp_Ph_CI_1->Add(Comp_Ph_CI_1_R,-ratio);

        Comp_MM_WE_0->Sumw2();
        Comp_MM_WE_0->Add(Comp_MM_WE_0_R,-ratio);
        Comp_MM_WE_1->Sumw2();
        Comp_MM_WE_1->Add(Comp_MM_WE_1_R,-ratio);
        Comp_Ph_WE_0->Sumw2();
        Comp_Ph_WE_0->Add(Comp_Ph_WE_0_R,-ratio);
        Comp_Ph_WE_1->Sumw2();
        Comp_Ph_WE_1->Add(Comp_Ph_WE_1_R,-ratio);

        Comp_MM_WI_0->Sumw2();
        Comp_MM_WI_0->Add(Comp_MM_WI_0_R,-ratio);
        Comp_MM_WI_1->Sumw2();
        Comp_MM_WI_1->Add(Comp_MM_WI_1_R,-ratio);
        Comp_Ph_WI_0->Sumw2();
        Comp_Ph_WI_0->Add(Comp_Ph_WI_0_R,-ratio);
        Comp_Ph_WI_1->Sumw2();
        Comp_Ph_WI_1->Add(Comp_Ph_WI_1_R,-ratio);

        Comp_MM_TE_0->Sumw2();
        Comp_MM_TE_0->Add(Comp_MM_TE_0_R,-ratio);
        Comp_MM_TE_1->Sumw2();
        Comp_MM_TE_1->Add(Comp_MM_TE_1_R,-ratio);
        Comp_Ph_TE_0->Sumw2();
        Comp_Ph_TE_0->Add(Comp_Ph_TE_0_R,-ratio);
        Comp_Ph_TE_1->Sumw2();
        Comp_Ph_TE_1->Add(Comp_Ph_TE_1_R,-ratio);

        Comp_MM_TI_0->Sumw2();
        Comp_MM_TI_0->Add(Comp_MM_TI_0_R,-ratio);
        Comp_MM_TI_1->Sumw2();
        Comp_MM_TI_1->Add(Comp_MM_TI_1_R,-ratio);
        Comp_Ph_TI_0->Sumw2();
        Comp_Ph_TI_0->Add(Comp_Ph_TI_0_R,-ratio);
        Comp_Ph_TI_1->Sumw2();
        Comp_Ph_TI_1->Add(Comp_Ph_TI_1_R,-ratio);

        Comp_CS_MM->Sumw2();
        Comp_CS_MM->Add(Comp_CS_MM_R,-ratio);
        Reco_CS_MM->Sumw2();
        Reco_CS_MM->Add(Reco_CS_MM_R,-ratio);
    }
    if(save_randoms && !(IsMCFile()))
    {
        Tagg_0_R->Write();
        Tagg_1_R->Write();

        Inc_0_R->Write();
        Inc_1_R->Write();

        Pi0_OA_R->Write();
        Pi0_OA_Cut_R->Write();

        Pi0_MM_NE_0_R->Write();
        Pi0_MM_NI_0_R->Write();
        Pi0_MM_CE_0_R->Write();
        Pi0_MM_CI_0_R->Write();
        Pi0_MM_WE_0_R->Write();
        Pi0_MM_WI_0_R->Write();
        Pi0_MM_TE_0_R->Write();
        Pi0_MM_TI_0_R->Write();

        Pi0_Re_All_R->Write();
        Pi0_Re_Det_R->Write();
        Pi0_Re_Dif_R->Write();
        Pi0_Re_NoE_R->Write();

        PiP_OA_R->Write();
        PiP_OA_Cut_R->Write();

        PiP_MM_NE_0_R->Write();
        PiP_MM_NI_0_R->Write();
        PiP_MM_CE_0_R->Write();
        PiP_MM_CI_0_R->Write();
        PiP_MM_WE_0_R->Write();
        PiP_MM_WI_0_R->Write();
        PiP_MM_TE_0_R->Write();
        PiP_MM_TI_0_R->Write();

        Comp_OA_R->Write();
        Comp_OA_Cut_R->Write();

        Comp_MM_NE_0_R->Write();
        Comp_MM_NI_0_R->Write();
        Comp_MM_CE_0_R->Write();
        Comp_MM_CI_0_R->Write();
        Comp_MM_WE_0_R->Write();
        Comp_MM_WI_0_R->Write();
        Comp_MM_TE_0_R->Write();
        Comp_MM_TI_0_R->Write();

        Comp_CS_MM_R->Write();
        Reco_CS_MM_R->Write();

        if(cir_beam || lin_beam)
        {
            Pi0_MM_NE_1_R->Write();
            Pi0_Ph_NE_0_R->Write();
            Pi0_Ph_NE_1_R->Write();

            Pi0_MM_NI_1_R->Write();
            Pi0_Ph_NI_0_R->Write();
            Pi0_Ph_NI_1_R->Write();

            Pi0_MM_CE_1_R->Write();
            Pi0_Ph_CE_0_R->Write();
            Pi0_Ph_CE_1_R->Write();

            Pi0_MM_CI_1_R->Write();
            Pi0_Ph_CI_0_R->Write();
            Pi0_Ph_CI_1_R->Write();

            Pi0_MM_WE_1_R->Write();
            Pi0_Ph_WE_0_R->Write();
            Pi0_Ph_WE_1_R->Write();

            Pi0_MM_WI_1_R->Write();
            Pi0_Ph_WI_0_R->Write();
            Pi0_Ph_WI_1_R->Write();

            Pi0_MM_TE_1_R->Write();
            Pi0_Ph_TE_0_R->Write();
            Pi0_Ph_TE_1_R->Write();

            Pi0_MM_TI_1_R->Write();
            Pi0_Ph_TI_0_R->Write();
            Pi0_Ph_TI_1_R->Write();

            PiP_MM_NE_1_R->Write();
            PiP_Ph_NE_0_R->Write();
            PiP_Ph_NE_1_R->Write();

            PiP_MM_NI_1_R->Write();
            PiP_Ph_NI_0_R->Write();
            PiP_Ph_NI_1_R->Write();

            PiP_MM_CE_1_R->Write();
            PiP_Ph_CE_0_R->Write();
            PiP_Ph_CE_1_R->Write();

            PiP_MM_CI_1_R->Write();
            PiP_Ph_CI_0_R->Write();
            PiP_Ph_CI_1_R->Write();

            PiP_MM_WE_1_R->Write();
            PiP_Ph_WE_0_R->Write();
            PiP_Ph_WE_1_R->Write();

            PiP_MM_WI_1_R->Write();
            PiP_Ph_WI_0_R->Write();
            PiP_Ph_WI_1_R->Write();

            PiP_MM_TE_1_R->Write();
            PiP_Ph_TE_0_R->Write();
            PiP_Ph_TE_1_R->Write();

            PiP_MM_TI_1_R->Write();
            PiP_Ph_TI_0_R->Write();
            PiP_Ph_TI_1_R->Write();

            Comp_MM_NE_1_R->Write();
            Comp_Ph_NE_0_R->Write();
            Comp_Ph_NE_1_R->Write();

            Comp_MM_NI_1_R->Write();
            Comp_Ph_NI_0_R->Write();
            Comp_Ph_NI_1_R->Write();

            Comp_MM_CE_1_R->Write();
            Comp_Ph_CE_0_R->Write();
            Comp_Ph_CE_1_R->Write();

            Comp_MM_CI_1_R->Write();
            Comp_Ph_CI_0_R->Write();
            Comp_Ph_CI_1_R->Write();

            Comp_MM_WE_1_R->Write();
            Comp_Ph_WE_0_R->Write();
            Comp_Ph_WE_1_R->Write();

            Comp_MM_WI_1_R->Write();
            Comp_Ph_WI_0_R->Write();
            Comp_Ph_WI_1_R->Write();

            Comp_MM_TE_1_R->Write();
            Comp_Ph_TE_0_R->Write();
            Comp_Ph_TE_1_R->Write();

            Comp_MM_TI_1_R->Write();
            Comp_Ph_TI_0_R->Write();
            Comp_Ph_TI_1_R->Write();
        }
    }

    delete Tagg_0_R;
    delete Tagg_1_R;

    delete Inc_0_R;
    delete Inc_1_R;

    delete Pi0_OA_R;
    delete Pi0_OA_Cut_R;

    delete Pi0_MM_NE_0_R;
    delete Pi0_MM_NE_1_R;
    delete Pi0_Ph_NE_0_R;
    delete Pi0_Ph_NE_1_R;

    delete Pi0_MM_NI_0_R;
    delete Pi0_MM_NI_1_R;
    delete Pi0_Ph_NI_0_R;
    delete Pi0_Ph_NI_1_R;

    delete Pi0_MM_CE_0_R;
    delete Pi0_MM_CE_1_R;
    delete Pi0_Ph_CE_0_R;
    delete Pi0_Ph_CE_1_R;

    delete Pi0_MM_CI_0_R;
    delete Pi0_MM_CI_1_R;
    delete Pi0_Ph_CI_0_R;
    delete Pi0_Ph_CI_1_R;

    delete Pi0_MM_WE_0_R;
    delete Pi0_MM_WE_1_R;
    delete Pi0_Ph_WE_0_R;
    delete Pi0_Ph_WE_1_R;

    delete Pi0_MM_WI_0_R;
    delete Pi0_MM_WI_1_R;
    delete Pi0_Ph_WI_0_R;
    delete Pi0_Ph_WI_1_R;

    delete Pi0_MM_TE_0_R;
    delete Pi0_MM_TE_1_R;
    delete Pi0_Ph_TE_0_R;
    delete Pi0_Ph_TE_1_R;

    delete Pi0_MM_TI_0_R;
    delete Pi0_MM_TI_1_R;
    delete Pi0_Ph_TI_0_R;
    delete Pi0_Ph_TI_1_R;

    delete Pi0_Re_All_R;
    delete Pi0_Re_Det_R;
    delete Pi0_Re_Dif_R;
    delete Pi0_Re_NoE_R;

    delete PiP_OA_R;
    delete PiP_OA_Cut_R;

    delete PiP_MM_NE_0_R;
    delete PiP_MM_NE_1_R;
    delete PiP_Ph_NE_0_R;
    delete PiP_Ph_NE_1_R;

    delete PiP_MM_NI_0_R;
    delete PiP_MM_NI_1_R;
    delete PiP_Ph_NI_0_R;
    delete PiP_Ph_NI_1_R;

    delete PiP_MM_CE_0_R;
    delete PiP_MM_CE_1_R;
    delete PiP_Ph_CE_0_R;
    delete PiP_Ph_CE_1_R;

    delete PiP_MM_CI_0_R;
    delete PiP_MM_CI_1_R;
    delete PiP_Ph_CI_0_R;
    delete PiP_Ph_CI_1_R;

    delete PiP_MM_WE_0_R;
    delete PiP_MM_WE_1_R;
    delete PiP_Ph_WE_0_R;
    delete PiP_Ph_WE_1_R;

    delete PiP_MM_WI_0_R;
    delete PiP_MM_WI_1_R;
    delete PiP_Ph_WI_0_R;
    delete PiP_Ph_WI_1_R;

    delete PiP_MM_TE_0_R;
    delete PiP_MM_TE_1_R;
    delete PiP_Ph_TE_0_R;
    delete PiP_Ph_TE_1_R;

    delete PiP_MM_TI_0_R;
    delete PiP_MM_TI_1_R;
    delete PiP_Ph_TI_0_R;
    delete PiP_Ph_TI_1_R;

    delete Comp_OA_R;
    delete Comp_OA_Cut_R;

    delete Comp_MM_NE_0_R;
    delete Comp_MM_NE_1_R;
    delete Comp_Ph_NE_0_R;
    delete Comp_Ph_NE_1_R;

    delete Comp_MM_NI_0_R;
    delete Comp_MM_NI_1_R;
    delete Comp_Ph_NI_0_R;
    delete Comp_Ph_NI_1_R;

    delete Comp_MM_CE_0_R;
    delete Comp_MM_CE_1_R;
    delete Comp_Ph_CE_0_R;
    delete Comp_Ph_CE_1_R;

    delete Comp_MM_CI_0_R;
    delete Comp_MM_CI_1_R;
    delete Comp_Ph_CI_0_R;
    delete Comp_Ph_CI_1_R;

    delete Comp_MM_WE_0_R;
    delete Comp_MM_WE_1_R;
    delete Comp_Ph_WE_0_R;
    delete Comp_Ph_WE_1_R;

    delete Comp_MM_WI_0_R;
    delete Comp_MM_WI_1_R;
    delete Comp_Ph_WI_0_R;
    delete Comp_Ph_WI_1_R;

    delete Comp_MM_TE_0_R;
    delete Comp_MM_TE_1_R;
    delete Comp_Ph_TE_0_R;
    delete Comp_Ph_TE_1_R;

    delete Comp_MM_TI_0_R;
    delete Comp_MM_TI_1_R;
    delete Comp_Ph_TI_0_R;
    delete Comp_Ph_TI_1_R;

    delete Comp_CS_MM_R;
    delete Reco_CS_MM_R;

    if(!cir_beam && !lin_beam)
    {
        delete Pi0_MM_NE_1;
        delete Pi0_Ph_NE_0;
        delete Pi0_Ph_NE_1;

        delete Pi0_MM_NI_1;
        delete Pi0_Ph_NI_0;
        delete Pi0_Ph_NI_1;

        delete Pi0_MM_CE_1;
        delete Pi0_Ph_CE_0;
        delete Pi0_Ph_CE_1;

        delete Pi0_MM_CI_1;
        delete Pi0_Ph_CI_0;
        delete Pi0_Ph_CI_1;

        delete Pi0_MM_WE_1;
        delete Pi0_Ph_WE_0;
        delete Pi0_Ph_WE_1;

        delete Pi0_MM_WI_1;
        delete Pi0_Ph_WI_0;
        delete Pi0_Ph_WI_1;

        delete Pi0_MM_TE_1;
        delete Pi0_Ph_TE_0;
        delete Pi0_Ph_TE_1;

        delete Pi0_MM_TI_1;
        delete Pi0_Ph_TI_0;
        delete Pi0_Ph_TI_1;

        delete PiP_MM_NE_1;
        delete PiP_Ph_NE_0;
        delete PiP_Ph_NE_1;

        delete PiP_MM_NI_1;
        delete PiP_Ph_NI_0;
        delete PiP_Ph_NI_1;

        delete PiP_MM_CE_1;
        delete PiP_Ph_CE_0;
        delete PiP_Ph_CE_1;

        delete PiP_MM_CI_1;
        delete PiP_Ph_CI_0;
        delete PiP_Ph_CI_1;

        delete PiP_MM_WE_1;
        delete PiP_Ph_WE_0;
        delete PiP_Ph_WE_1;

        delete PiP_MM_WI_1;
        delete PiP_Ph_WI_0;
        delete PiP_Ph_WI_1;

        delete PiP_MM_TE_1;
        delete PiP_Ph_TE_0;
        delete PiP_Ph_TE_1;

        delete PiP_MM_TI_1;
        delete PiP_Ph_TI_0;
        delete PiP_Ph_TI_1;

        delete Comp_MM_NE_1;
        delete Comp_Ph_NE_0;
        delete Comp_Ph_NE_1;

        delete Comp_MM_NI_1;
        delete Comp_Ph_NI_0;
        delete Comp_Ph_NI_1;

        delete Comp_MM_CE_1;
        delete Comp_Ph_CE_0;
        delete Comp_Ph_CE_1;

        delete Comp_MM_CI_1;
        delete Comp_Ph_CI_0;
        delete Comp_Ph_CI_1;

        delete Comp_MM_WE_1;
        delete Comp_Ph_WE_0;
        delete Comp_Ph_WE_1;

        delete Comp_MM_WI_1;
        delete Comp_Ph_WI_0;
        delete Comp_Ph_WI_1;

        delete Comp_MM_TE_1;
        delete Comp_Ph_TE_0;
        delete Comp_Ph_TE_1;

        delete Comp_MM_TI_1;
        delete Comp_Ph_TI_0;
        delete Comp_Ph_TI_1;
    }

    // Write all GH1's and TObjects defined in this class
    if(!(GTreeManager::Write())) return false;

    return true;
}
