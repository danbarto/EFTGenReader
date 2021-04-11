# Wrapper to runGridpackValidation.C
# This code reads output tress produced by EFTLHEReader.cc and groups them by process+coeff, then processes them in runGridpackValidation.C

import subprocess
import datetime
import os
from utils import regex_match, getInfoByTag, getDirectories, groupByProcess, groupByCoefficient

ALL_INFO = [

    ######################################
    #### TOP19-001 Round 6 FP samples ####
    ######################################

    {# HanModelV4 FP R6 B2: ttWJet and ttZJet (and ttH but bad starting point for MC stats)
        'tag': 'ttXjet-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttllNuNuJetNoHiggs","ttlnuJet"], # Don't incluede ttH
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch2/",
    },
    {# HanModelV4 FP R6 B3: tHq
        'tag': 'tHq4f-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch3/",
    },
    {# HanModelV4 FP R6 B4: tZq
        'tag': 'tZq4f-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch4/",
    },
    {# HanModelV4 FP R6 B7: ttHJet
        'tag': 'ttHjet-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch7/",
    },

    ############################

    ### Various TOP-19-001 and pheno paper studies ###

    {# HanModelV4 ttHJet dedicated ctG=[-3,3] axis scan ####### AXIS SCAN #######
        # Apparently I did not specify basepath when I added this... would probably need to go searching for the tag name if we wanted to use it
        'tag': 'ttHJet_HanV4ctGAxisScan_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 ttHJet dedicated ctW=[-4,4] axis scan (for smeft comp, qed1, qcd2, dim6=2) ####### AXIS SCAN #######
        # Note: dim6=2 is NOT the same as NP=2, so can't compare this with NLO as was intended
        'tag': 'ttHJet_HanV4_cbW-AxisScan-withRwgt_smeftComp_QED1_QCD2_DIM62-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/",
    },
    {# HanModelV4 ttHJet xqcut and qcut scan (12 samples in all)
        # Apparently I did not specify basepath when I added this... would probably need to go searching for the tag name if we wanted to use it
        'tag': 'ttHJet_HanV4xqcutTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttHJet"],
        'c_wl': ["HanV4ttXjetxqcut10qCut19"],
        'r_wl': [],
    },
    {# HanModelV4 ttXJet comp with SMEFT (in proc card: QED=1, QCD=2, DIM6=2, and ttZ not ttll, ttW not ttlnu)
        # NOTE: This is a bad comparison! DIM6 should NOT be 2 for this comparison (NP=2 in smeft NLO is more or less DIM6=1, not DIM6=2 in dim6Top). Do not use!
        'tag': 'ttXJet_HanV4_semftComp_QED1_QCD2_DIM62-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        #'p_wl': ["ttHJetSMEFTcomp","ttWJetSMEFTcomp","ttZJetSMEFTcomp"],
        'p_wl': ["ttHJetSMEFTcomp"],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/",
    },
    {# HanModelV4 ttHJet, comp with NLO, so QED=1, QCD=2 (and DIM6=1) (Reza's NLO have NP=2, but this is NOT equivalent to DIM6=2!, so this DIM6=1 is ok)
     # NOTE: maching should be off, so this sample is actully NOT ok since matching was on! Do not use!
        'tag': 'ttHJet_HanV4_withRwgt_smeftComp_QED1_QCD2_DIM61-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# dim6Top May 19 2020 update version, ttH, ttHJet (for comp with HanV4)
        'tag': 'ttH-ttHJet_dim6Top-vMay2020-normChromoTrue-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# HanV4 ttW ttWJet ttZ ttZJet (i.e. NOT ttlnu!) for no order constraints and for QED=1,QCD=2 for comp to Reza's NLO ttW, ttZ
        # NOTE: The ttWJet and ttZJet here are actually ONLY Jet (i.e. only +1p, not 0+1p), so don't use them for the comparison
        # So in effect this is a ttW, ttZ sample with and without QED=1,QCCD=2
        'tag': 'ttW-ttWJet-ttZ-ttZJet_QED-QCD-order-tests-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True, # ForPheno (before moreStats-goodStartPt remake)
        'include': False,
        'p_wl': ["ttW","ttZ"], # Don't use ttWJet, ttZJet #
        'c_wl': ["HanV4withRwgt"], # Skipping the QED1,QCD2 ones, not good for NLO comp #
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# HanV4 ttWJet ttZJet (i.e. NOT ttlnu!) with no order constraints and with QED=1,QCD=2 
        # NOTE: The ttWJet and ttZJet should actually be 0+1p (unlike in "ttW-ttWJet-ttZ-ttZJet_QED-QCD-order-tests-GEN"), but matching was on
        # and it should not be (since there is no extra jet) so don't use the QED=1, QCD=2 ones
        # So in effect this is a ttWJet, ttZJet no order constraints sample
        'tag': 'ttWJet-ttZJet_HanV4-0plus1p-QCDQED-OrderChecks-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True, # ForPheno (before moreStats-goodStartPt remake)
        'include': False,
        'p_wl': [],
        'c_wl': ["HanV40plus1pwithRwgt"], # Don't use "HanV40plus1pQED1QCD2withRwgt" as it wrongly has matching on
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# HanV4 ttHJet, ttWJet, ttZJet all with QED=1, QCD=2, and matching off (since no extra parton), also ttH with QED=1, QCD=2 
        'tag': 'ttH-ttXJet_HanV4-0pttH-0plus1pttXJet-noMatching-QCDQED-OrderChecks-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttHJet","ttWJet","ttZJet","ttH"], # ttH is sort of annoyingly stuck in withe these ttXJet samples
        #'p_wl': ["ttHJet","ttWJet","ttZJet"], # ttH is sort of annoyingly stuck in withe these ttXJet samples
        'p_wl': ["ttHJet"], # ttH is sort of annoyingly stuck in withe these ttXJet samples
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# HanV4 ttX (ttH, ttW, ttZ), ttXJet (ttHJet, ttWJet, ttZJet) all with QED=1, QCD=3
        # Note: The QCD=3 _should_ be correct for comparing with the NLO (which have QCD=2, but also [QCD])
        'tag': 'ttX-ttXJet_HanV4-QED1-QCD3-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True, # ForPheno (before moreStats-goodStartPt remake)
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# ttV, ttVJet start point checks
        'tag': 'ttV-ttVJet_HanV4_QED1-and-noConstraints_startPtChecks-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttZ","ttZJet"],
        #'p_wl': ["ttW","ttWJet"],
        'p_wl': ["ttH","ttW","ttZ"],
        'c_wl': ["HanV4QED1QCD3startPtChecks"],
        #'c_wl': ["HanV4startPtChecks"],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    {# ttX, ttXJet QED2 start point checks
        'tag': 'ttX-ttXJet_HanV4-QED2-startPtChecks-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttH","ttHJet"],
        #'p_wl': ["ttW","ttWJet"],
        'p_wl': ["ttHJet","ttWJet","ttZJet"],
        'c_wl': [],
        #'r_wl': ["run1"],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
    ############################

    ### The final samples for the pheno paper ###
    {# HanV4 ttX (ttH, ttW, ttZ) ttXJet, all with and without QED1QCD3 constrints
        # Note: 300k events, and the ttV samples are at better starting points than the original, so the stats should be better
        # FOR PHENO PAPER 300k events
        'tag': 'ttX-ttXJet_HanV4_QED1-and-noConstraints_moreStats-goodStartPt-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': [],
        #'p_wl': ["ttH"],
        'p_wl': ["ttHJet"],
        #'p_wl': ["ttWJet"],
        #'p_wl': ["ttZJet"],
        #'p_wl': ["ttH","ttHJet"], 
        #'p_wl': ["ttW","ttWJet"], 
        #'p_wl': ["ttZ","ttZJet"], 
        #'p_wl': ["ttZ","ttZJet","ttW","ttWJet"], 
        #'p_wl': ["ttH","ttHJet","ttZ","ttZJet","ttW","ttWJet"], 
        #'c_wl': [],
        #'c_wl': ["HanV4goodStartPt"], # ttV: HanV4goodStartPt, HanV4QED1QCD3goodStartPt
        #'c_wl': ["HanV4lModel16DttllScanpoints","HanV4ModelNoJets16DttllScanpoints"], # This is for ttH without any cuts
        'c_wl': ["HanV4goodStartPt","HanV4lModel16DttllScanpoints","HanV4ModelNoJets16DttllScanpoints"], # All samples, no QED cuts
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },

    {# ttH start point check, just to make sure the "normal" start point is okay for ttH 0j (can't believe we havn't checked this before?)
        'tag': 'ttH_HanV4ttH0pStartPtDoubleCheck-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },

    #####################################################
    # New test samples made for the full run 2 analysis #
    #####################################################

    ### Samples for preliminary validation studies (for model and framework) ###

    {# ttHJet, with new (May 2020) dim6Top model, but still using old (pre updated) genproductions scripts. Wiht gs norm T and F
        # Note: Something seems to have gone wrong wiht these gridpacks (at least the gs True one), as it seems the integrate step ran very fast, and two of the 1d curves are flata. Possibly it ran into issues and skipped something.
        # Do not use!
        'tag': 'ttHJet_dim6TopMay20_testing-old-genprod-updated-model-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# ttHJet, with the new genprod framework and HanV4 as well as (May 2020) dim6Top model (with gs T and F)
        'tag': 'ttHJet_testUpdateGenprod-testModels-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': ["testUpdateGenprodHanV4"],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# ttHJet, with the old  genprod framework and HanV4 as well as (May 2020) dim6Top model (with gs T and F)
        'tag': 'ttHJet_testOldGenprod-testModels-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# ttXJet, tXq: All five signal processes from TOP-19-001, but with new genprod and new dim6Top
        'tag': 'ttXJet-tXq_testUpdateGenproddim6TopMay20GST-testAllProcs-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# ttbar and ttbar+j samples with three different starting points with new dim6Top with new genprod to check WC dependence
        'tag': 'ttbar-ttbarJet_newGenprodNewDim6Top_testBkgDependence-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': ["run2"],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# ttXJet, tXq: All five signal processes from TOP-19-001, but with new genprod and new dim6Top, BUT with MG 2.6.0
        'tag': 'ttXJet-tXq_testUpdateGenprodMG260dim6TopMay20GST-testAllProcs-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    {# tHq sample with old genprod framework and HanV4 for comparison to reference
        'tag': 'tHq4f_testOldGenprod-HanV4-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },
    # 1d samples
    {# ttHJet, 1d scans for cptb (since we saw a discrepency) and cpt (for reference)
        'tag': 'ttHJet_testGenprodVersions-dim6TopMat20GST-1dScans-cpt-cptb-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        #'c_wl': ["cptTestOldGenprodTestdim6TopMay20GSTAxisScan","cptbTestOldGenprodTestdim6TopMay20GSTAxisScan"],
        'c_wl': ["cptV5TestUpdateGenprodTestdim6TopMay20GSTAxisScan","cptV6TestUpdateGenprodTestdim6TopMay20GSTAxisScan","cptbV4TestUpdateGenprodTestdim6TopMay20GSTAxisScan","cptbV2TestUpdateGenprodTestdim6TopMay20GSTAxisScan"],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/PreliminaryStudies/",
    },

    ### UL samples for testing UL configs (summary tree step made with 10_6_8) ###

    { # Test UL 16 samples, all 5 processes
        'tag': 'ttXJet-tXq_GEN_ULCheckUL16-GEN',
        'grp_name': '',
        'version': 'v6',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # Test UL 16APV samples, all 5 processes
        'tag': 'ttXJet-tXq_GEN_ULCheckUL16APV-GEN',
        'grp_name': '',
        'version': 'v6',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # Test UL 17 samples, all 5 processes
        'tag': 'ttXJet-tXq_GEN_ULCheckUL17-GEN',
        'grp_name': '',
        'version': 'v6',
        #'include': True,
        'include': False,
        'p_wl': ["ttHJet","ttlnuJet","tllq4fNoSchanWNoHiggs0p","tHq4f"],
        #'p_wl': ["ttHJet"],
        #'p_wl': []
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # Test UL 18 samples, all 5 processes
        'tag': 'ttXJet-tXq_GEN_ULCheckUL17-GEN',
        'grp_name': '',
        'version': 'v6',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },

    ### Checking dim6 syntaxes (for the studies for the double insertion studies) ###

    { # (UL17) Checking the dim6^2 syntax for ttHJet
        'tag': 'ttHJet_dim6TopMay20GST-checkDIM6Syntaxes_UL17-GEN',
        'grp_name': '',
        'version': 'v8',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking the dim6^2 syntax for all the processes (except ttll, which did not finish in time)
        'tag': 'ttXJet-tXq_dim6TopMay20GST-and-ECO_checkDIM6SQSyntax_UL17-GEN',
        'grp_name': '',
        'version': 'v9',
        #'include': True,
        'include': False,
        #'p_wl': ["ttHJet"],
        'p_wl': [],
        #'c_wl': ["testDIM6SQdim6TopMay20GST"],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking the dim6^2 syntax and dim6=1 syntaxes at various starting points
        'tag': 'ttHJet_dim6TopMay20GST-StartPtChecks_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking the dim6^2 syntax and dim6=1 syntaxes at two starting points, drop ctG from gridpack
        'tag': 'ttHJet_dim6TopMay20GST-ctG0-testForNegativeWgts-StartPtChecks_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking the dim6^2 syntax and dim6=1 syntaxes at two starting points, with JUST ctG and ctp in the gridpack
        'tag': 'ttH-ttHJet_dim6TopMay20GST_JustctGctp-check-dim6syntaxes_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },

    ### Checking H and ll interference ###

    { # (UL17) Checking samples with H produced seperately to see if we are missing interferecne
        'tag': 'ttX-tXq_dim6TopMay20GST_testHiggsInterferenceWithLL_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v2', ## NOTE: v4 is actually latest version, but in this case do nto use (it's empty)
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking samples with H produced seperately to see if we are missing interferecne, v2 at different starting points
        'tag': 'ttX-tXq_dim6TopMay20GST_testV2HiggsInterferenceWithLL_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking samples with H produced seperately to see if we are missing interferecne, the two samples above were both processed together here, just to make sure no differences
        'tag': 'ttX-tXq_dim6TopMay20GST_testBothSetsOfSamplesTogetherHiggsInterferenceWithLL_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Checking samples with H produced seperately to see if we are missing interferecne, the two samples above were both processed together here, just to make sure no differences
        'tag': 'ttX-tXq_dim6TopMay20GST_testBothSetsOfSamplesTogetherHiggsInterferenceWithLL-rerunLHE_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },
    { # (UL17) Remake ttHTOll sample with mll cut to be consistent wiht ttll samples
        'tag': 'ttHTOll_testHiggsInterferenceWithllWithMMLL10dim6TopMay20GST_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ULChecks/",
    },

    ### Full Run2 Validation Checks ###

    { # (UL17) 22d ttHJet, ttlnuJet, tllq, tHq (no ttllJet since it took too long...)
        'tag': 'ttHJet-ttlnuJet-tllq-tHq_dim6TopMay20GST_all22WCsBaselineStartPtTOP19001_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },
    { # (UL17) 22d ttllJet and ttbarjet followup on the sample above
        'tag': 'ttllJet-ttbarJet_dim6TopMay20GST_all22WCsBaselineStartPtTOP19001_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },
    { # (UL17) dim6=0 for a SM reference
        'tag': 'ttXJet-tXq-ttbarJet_dim6TopMay20GST_all22WCsDim6Eq0_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },
    { # (UL17) ttHJet, ttlnuJet, ttbarJet, tllq, tHq: 22d samples at 7 different starting points
        'tag': 'ttHJet-ttlnuJet-ttbarJet-tllq-tHq_dim6TopMay20GST_all22WCsStartPtCheckV2_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },
    { # (UL17) ttHJet, ttlnuJet, ttllJet, ttbarJet, tllq, tHq: 22d samples at 7 different starting points
        'tag': 'ttHJet-ttlnuJet-ttllJet-ttbarJet-tllq-tHq_dim6TopMay20GST_all22WCsStartPtCheck_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },
    { # (UL17) ttHJet, ttlnuJet, ttllJet, ttbarJet qCut scans
        'tag': 'ttXJet_dim6TopMay20GST_run0StartPt_qCutScan_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': True,
        #'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },

    # Axis scans #
    { # (UL17) ttbar axis scan
        'tag': 'ttbarJet_dim6TopMay20GST_1dAxisScans-2heavy-2heavy2light_GENUL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/FullR2Studies/ValidationChecks/",
    },

    ######################################
    ### Pheno paper JHEP review checks ###
    ######################################

    { # (UL17) ttW, ttWJet: Checking 1d cbW and cptb
        'tag': 'ttW-ttWJet_cbW-cptb-1dChecks_dim6TopMay20GST_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttW"],
        'p_wl': ["ttWJet"],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/ForPhenoJhepReviewStudies/"
    },
    { # (UL17) ttWJet double check the qCut scans for the reviewer
        'tag': 'ttWJet_sampleForDoubleCheckingQcut_dim6TopMay20GST_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/ForPhenoJhepReviewStudies/"
    },
    { # (UL17) ttZJet double check the qCut scans for the reviewer
        'tag': 'ttZJet_sampleForDoubleCheckingQcut_dim6TopMay20GST_GEN_UL17-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/ForPhenoJhepReviewStudies/"
    },

]

# These are the tags to get the runs which should be used as reference points in the plotting
REF_TAGS = [
    # Axis scans (don't need to include as "True" in ALL_INFO)
    #'ttHJet_HanV4ctGAxisScan_analysisEtaCut-GEN'
    #'ttHJet_HanV4_cbW-AxisScan-withRwgt_smeftComp_QED1_QCD2_DIM62-GEN' # Note, the proc name is ttHJetSMEFTcomp, might need to take care of this in REF_PROCESS_MAP
    #'ttHJet_testGenprodVersions-dim6TopMat20GST-1dScans-cpt-cptb-GEN'
    #'ttbarJet_dim6TopMay20GST_1dAxisScans-2heavy-2heavy2light_GENUL17-GEN'
]

# Dictionary to map certain MG processes to another for use in finding reference samples
REF_PROCESS_MAP = {
    #'tllq4fMatchedNoHiggs': 'tllq',
    'ttHJet': 'ttH',
    'ttH': 'ttHJet',
    #'ttH': 'ttHJetSMEFTcomp'
    #'ttlnuJet': 'ttlnu',
    #'tHq4fMatched': 'tHq',
    #'ttllNuNuJetNoHiggs': 'ttll'
}

# These are not really used anymore
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2019_04_19/"
HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/"

def runByProcess(save_info_dir_path):

    cleanName = True
    #NOTE: The output name could be duplicated and overwrite a previous run
    all_grouped_file_dirs_dict = {}
    spacing = 0
    for idx,info in enumerate(ALL_INFO):
        if not info['include']:
            continue
        # For multiple base paths:
        if "basepath" in info.keys():
            path = os.path.join(info['basepath'],info['tag'],info['version'])
        else:
            path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
        grouped_dirs = groupByProcess(path,info['tag'],'',info['p_wl'],info['c_wl'],info['r_wl'])
        #print "\nThe current grouped dirs are:" , grouped_dirs , "\n"
        for tup,dirs in grouped_dirs.iteritems():
            if cleanName: # Note: This is also for grouping
                '''
                if "ttH" in tup[1]:
                    proc = "ttH"
                elif "ttlnu" in tup[1] or "ttW" in tup[1]:
                    proc = "ttW"
                elif "ttll" in tup[1] or "ttZ" in tup[1]:
                    proc = "ttZ"
                elif "tHq" in tup[1] or "tHTOllq" in tup[1]:
                    proc = tup[1]
                elif "tllq" in tup[1]:
                    proc = tup[1]
                elif "ttbar" in tup[1]:
                    #proc = tup[1]
                    proc = "ttbar"
                else:
                    print "\nError: Unknown process" , tup[1] , "exiting"
                    raise BaseException
            else:
                proc = tup[1]
            '''
            proc = tup[1]  # Do not do any grouping
            if not all_grouped_file_dirs_dict.has_key(proc):
                all_grouped_file_dirs_dict[proc] = []
            all_grouped_file_dirs_dict[proc].extend(dirs)

    # Run root macro once per process
    count = 0
    print "\nAll grouped dirs:" , all_grouped_file_dirs_dict
    for proc,fdirs in all_grouped_file_dirs_dict.iteritems():
        output_name = proc

        ref_dirs = []
        for rtag in REF_TAGS:
            print "tag!" , rtag
            info_list = getInfoByTag(ALL_INFO,rtag)
            for info in info_list:
                if "basepath" in info.keys():
                    ref_path = os.path.join(info['basepath'],info['tag'],info['version'])
                else:
                    ref_path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
                search_proc = proc
                print "proc: " , proc , REF_PROCESS_MAP.keys()
                if REF_PROCESS_MAP.has_key(proc):
                    #search_proc = REF_PROCESS_MAP[process]
                    search_proc = REF_PROCESS_MAP[proc]
                ref_dirs += getDirectories(ref_path,
                    p_wl=[search_proc],
                    #c_wl=[],
                    c_wl=info["c_wl"],
                    r_wl=info["r_wl"]
                )
        print "\nFile dirs:"
        for x in fdirs:
            print "\t{}".format(x)
        print "\nRef dirs:" , ref_dirs , "\n"

        dir_inputs = 'dirpaths.txt'
        ref_inputs = 'refpaths.txt'
        with open(dir_inputs,'w') as f:
            # Write the path directories to all relevant coeffs/runs
            for fd in sorted(fdirs):
                l = "%s\n" % (fd)
                f.write(l)
        with open(ref_inputs,'w') as f:
            for fd in sorted(ref_dirs):
                l = "%s\n" % (fd)
                f.write(l)

        print ""
        print "[%d/%d] %s (dirs %d, ref %d):" % (count+1,len(all_grouped_file_dirs_dict.keys()),output_name.ljust(spacing),len(fdirs),len(ref_dirs))
        #subprocess.check_call(['root','-b','-l','-q','runGridpackValidation.C+(\"%s\",\"%s\",\"%s\")' % (output_name,dir_inputs,ref_inputs)])
        subprocess.check_call(['root','-b','-l','-q','runGridpackValidation.C+(\"%s\",\"%s\",\"%s\",\"%s\")' % (output_name,save_info_dir_path,dir_inputs,ref_inputs)])
        count += 1

def runByCoeff(tags,runs):
    procs = ["ttH","ttll","ttlnu","tllq"]
    for idx,tag in enumerate(tags):
        print "[%d/%d] %s:" % (idx+1,len(tags),tag)
        fdirs = []  # fdirs can have runs from different processes/grp_tags/runs, so make sure we can handle that
        for idx,info in enumerate(ALL_INFO):
            if not info['include']:
                continue
            path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
            fdirs += getDirectories(path,
                p_wl=procs,
                c_wl=tag,
                r_wl=runs
            )
        dir_inputs = 'dirpaths.txt'
        with open(dir_inputs,'w') as f:
            for fd in sorted(fdirs):
                l = "%s\n" % (fd)
                f.write(l)
        subprocess.check_call(['root','-b','-l','-q','runLayeredPlots.C+(\"%s\")' % (dir_inputs)])
    return

# This function is very bad, we really should just write to the json automatically (not by hand).
# But I don't know how to do that from c++, so doing it manually in runGridpackValidation.C, so have to set it up manually here too.
def make_output_file_for_fit_info(save_file,action=None):
    timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')
    if action == "Begin":
        if os.path.isfile(save_file):
            print "\nERROR: File for saving fit info already exists: {f}. Exiting.\n".format(f=save_file)
            raise Exception
        with open (save_file,"w") as f:
            f.write("{\n")
            f.write("\t\"timestamp_start\" : \"{t}\",\n".format(t=timestamp_tag))
    elif action == "End":
        with open (save_file,"a") as f:
            f.write("\t\"timestamp_end\" : \"{t}\"\n".format(t=timestamp_tag))
            f.write("}")
    else:
        print "\nERROR: Unknown option {o}. Exiting.\n".format(o=action)

timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')

save_dir = "testing"
save_name = str(timestamp_tag)+".txt"

#save_dir = "fit_coeffs/start_pt_checks/sample_info/"
#save_name = "fit_info_dim6eq0.json"
#save_name = "fit_info_top19001_samples.json"
#save_name = "fit_info_all22WCsBaselineStartPtTOP19001.json"
#save_name = "fit_info_all22WCsStartPtCheck.json"
#save_name = "fit_info_all22WCsStartPtCheckV2.json"

save_file = os.path.join(save_dir,save_name)

make_output_file_for_fit_info(save_file,"Begin")
runByProcess(save_file)
make_output_file_for_fit_info(save_file,"End")
