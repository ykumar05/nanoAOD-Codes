//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 20 09:54:49 2019 by ROOT version 6.16/00
// from TTree Events/Events
// found on file: DYM50_tree1.root
//////////////////////////////////////////////////////////

/*
This header file contains a class called RHAna. 
(1) First the leaves are declared. This is a list of variables available to you.
(2) Then one branch per leaf is declared (these are lines such as TBranch *b_name)
(3) Then functions are declared. Focus on the User Added Functions first.
(4) Then variables are declared.
(5) An Init function which connects variables in code, to the input file.Input file is the skim.root and the variable in the code are defined as below
(6) ReadLimited function which if activated, only reads in few branches to speed
    up the code. 
 */




#ifndef RHAna_h
#define RHAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "TLorentzVector.h"

#include <vector>
#include <fstream>
#include <iostream>
// Header file for the classes stored in the TTree if any.
using namespace std;

// Header file for the classes stored in the TTree if any.

class RHAna : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   /* Float_t         btagWeight_CSVV2; */
   /* Float_t         btagWeight_CMVA; */
   /* Float_t         CaloMET_phi; */
   /* Float_t         CaloMET_pt; */
   /* Float_t         CaloMET_sumEt; */
   /* Float_t         ChsMET_phi; */
   /* Float_t         ChsMET_pt; */
   /* Float_t         ChsMET_sumEt; */
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[100];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[100];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[100];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[100];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[100];   //[nElectron]
   Float_t         Electron_dxy[100];   //[nElectron]
   Float_t         Electron_dxyErr[100];   //[nElectron]
   Float_t         Electron_dz[100];   //[nElectron]
   Float_t         Electron_dzErr[100];   //[nElectron]
   Float_t         Electron_eCorr[100];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[100];   //[nElectron]
   Float_t         Electron_energyErr[100];   //[nElectron]
   Float_t         Electron_eta[100];   //[nElectron]
   Float_t         Electron_hoe[100];   //[nElectron]
   Float_t         Electron_ip3d[100];   //[nElectron]
   Float_t         Electron_jetRelIso[100];   //[nElectron]
   Float_t         Electron_mass[100];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[100];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[100];   //[nElectron]
   /* Float_t         Electron_mvaFall17V1Iso[100];   //[nElectron] */
   /* Float_t         Electron_mvaFall17V1noIso[100];   //[nElectron] */
   /* Float_t         Electron_mvaFall17V2Iso[100];   //[nElectron] */
   /* Float_t         Electron_mvaFall17V2noIso[100];   //[nElectron] */
   /* Float_t         Electron_mvaSpring16GP[100];   //[nElectron] */
   /* Float_t         Electron_mvaSpring16HZZ[100];   //[nElectron] */
   Float_t         Electron_pfRelIso03_all[100];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[100];   //[nElectron]
   Float_t         Electron_phi[100];   //[nElectron]
   Float_t         Electron_pt[100];   //[nElectron]
   Float_t         Electron_r9[100];   //[nElectron]
   Float_t         Electron_sieie[100];   //[nElectron]
   Float_t         Electron_sip3d[100];   //[nElectron]
   Float_t         Electron_mvaTTH[100];   //[nElectron]
   Int_t           Electron_charge[100];   //[nElectron]
   Int_t           Electron_cutBased[100];   //[nElectron]
   /* Int_t           Electron_cutBased_Fall17_V1[100];   //[nElectron] */
   /* Int_t           Electron_cutBased_HLTPreSel[100];   //[nElectron] */
   /* Int_t           Electron_cutBased_Spring15[100];   //[nElectron] */
   /* Int_t           Electron_cutBased_Sum16[100];   //[nElectron] */
   Int_t           Electron_jetIdx[100];   //[nElectron]
   Int_t           Electron_pdgId[100];   //[nElectron]
   Int_t           Electron_photonIdx[100];   //[nElectron]
   Int_t           Electron_tightCharge[100];   //[nElectron]
   /* Int_t           Electron_vidNestedWPBitmap[100];   //[nElectron] */
   /* Int_t           Electron_vidNestedWPBitmapSpring15[100];   //[nElectron] */
   /* Int_t           Electron_vidNestedWPBitmapSum16[100];   //[nElectron] */
   Bool_t          Electron_convVeto[100];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[100];   //[nElectron]
   Bool_t          Electron_isPFcand[100];   //[nElectron]
   UChar_t         Electron_lostHits[100];   //[nElectron]
   /* Bool_t          Electron_mvaFall17V1Iso_WP80[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V1Iso_WP90[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V1Iso_WPL[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V1noIso_WP80[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V1noIso_WP90[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V1noIso_WPL[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2Iso_WP80[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2Iso_WP90[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2Iso_WPL[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2noIso_WP80[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2noIso_WP90[100];   //[nElectron] */
   /* Bool_t          Electron_mvaFall17V2noIso_WPL[100];   //[nElectron] */
   /* Bool_t          Electron_mvaSpring16GP_WP80[100];   //[nElectron] */
   /* Bool_t          Electron_mvaSpring16GP_WP90[100];   //[nElectron] */
   /* Bool_t          Electron_mvaSpring16HZZ_WPL[100];   //[nElectron] */
   UInt_t          nGenJet;
   Float_t         GenJet_eta[50];   //[nGenJet]
   Float_t         GenJet_mass[50];   //[nGenJet]
   Float_t         GenJet_phi[50];   //[nGenJet]
   Float_t         GenJet_pt[50];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[1000];   //[nGenPart]
   Float_t         GenPart_mass[1000];   //[nGenPart]
   Float_t         GenPart_phi[1000];   //[nGenPart]
   Float_t         GenPart_pt[1000];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[1000];   //[nGenPart]
   Int_t           GenPart_pdgId[1000];   //[nGenPart]
   Int_t           GenPart_status[1000];   //[nGenPart]
   Int_t           GenPart_statusFlags[1000];   //[nGenPart]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         genWeight;
   /* Float_t         LHEWeight_originalXWGTUP; */
   /* UInt_t          nLHEPdfWeight; */
   /* Float_t         LHEPdfWeight[500];   //[nLHEPdfWeight] */
   /* UInt_t          nLHEScaleWeight; */
   /* Float_t         LHEScaleWeight[50];   //[nLHEScaleWeight] */
   /* UInt_t          nPSWeight; */
   /* Float_t         PSWeight[20];   //[nPSWeight] */
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[50];   //[nIsoTrack]
   Float_t         IsoTrack_dz[50];   //[nIsoTrack]
   Float_t         IsoTrack_eta[50];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[50];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[50];   //[nIsoTrack]
   Float_t         IsoTrack_phi[50];   //[nIsoTrack]
   Float_t         IsoTrack_pt[50];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[50];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[50];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[50];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[50];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[50];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[50];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[50];   //[nIsoTrack]
   UInt_t          nJet;
   Float_t         Jet_area[100];   //[nJet]
   Float_t         Jet_btagCMVA[100];   //[nJet]
   Float_t         Jet_btagCSVV2[100];   //[nJet]
   Float_t         Jet_btagDeepB[100];   //[nJet]
   Float_t         Jet_btagDeepC[100];   //[nJet]
   Float_t         Jet_btagDeepFlavB[100];   //[nJet]
   Float_t         Jet_chEmEF[100];   //[nJet]
   Float_t         Jet_chHEF[100];   //[nJet]
   Float_t         Jet_eta[100];   //[nJet]
   Float_t         Jet_mass[100];   //[nJet]
   Float_t         Jet_muEF[100];   //[nJet]
   Float_t         Jet_neEmEF[100];   //[nJet]
   Float_t         Jet_neHEF[100];   //[nJet]
   Float_t         Jet_phi[100];   //[nJet]
   Float_t         Jet_pt[100];   //[nJet]
   Float_t         Jet_qgl[100];   //[nJet]
   Float_t         Jet_rawFactor[100];   //[nJet]
   Float_t         Jet_bRegCorr[100];   //[nJet]
   Float_t         Jet_bRegRes[100];   //[nJet]
   Int_t           Jet_electronIdx1[100];   //[nJet]
   Int_t           Jet_electronIdx2[100];   //[nJet]
   Int_t           Jet_jetId[100];   //[nJet]
   Int_t           Jet_muonIdx1[100];   //[nJet]
   Int_t           Jet_muonIdx2[100];   //[nJet]
   Int_t           Jet_nConstituents[100];   //[nJet]
   Int_t           Jet_nElectrons[100];   //[nJet]
   Int_t           Jet_nMuons[100];   //[nJet]
   Int_t           Jet_puId[100];   //[nJet]
   /* Float_t         GenMET_phi; */
   /* Float_t         GenMET_pt; */
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         METFixEE2017_phi;
   Float_t         METFixEE2017_pt;
   Float_t         MET_sumEt;
   UInt_t          nMuon;
   Float_t         Muon_dxy[100];   //[nMuon]
   Float_t         Muon_dxyErr[100];   //[nMuon]
   Float_t         Muon_dz[100];   //[nMuon]
   Float_t         Muon_dzErr[100];   //[nMuon]
   Float_t         Muon_eta[100];   //[nMuon]
   Float_t         Muon_ip3d[100];   //[nMuon]
   Float_t         Muon_jetRelIso[100];   //[nMuon]
   Float_t         Muon_mass[100];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[100];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[100];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[100];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[100];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[100];   //[nMuon]
   Float_t         Muon_phi[100];   //[nMuon]
   Float_t         Muon_pt[100];   //[nMuon]
   Float_t         Muon_ptErr[100];   //[nMuon]
   Float_t         Muon_segmentComp[100];   //[nMuon]
   Float_t         Muon_sip3d[100];   //[nMuon]
   Float_t         Muon_mvaTTH[100];   //[nMuon]
   Int_t           Muon_charge[100];   //[nMuon]
   Int_t           Muon_jetIdx[100];   //[nMuon]
   Int_t           Muon_nStations[100];   //[nMuon]
   Int_t           Muon_nTrackerLayers[100];   //[nMuon]
   Int_t           Muon_pdgId[100];   //[nMuon]
   Int_t           Muon_tightCharge[100];   //[nMuon]
   UChar_t         Muon_highPtId[100];   //[nMuon]
   Bool_t          Muon_inTimeMuon[100];   //[nMuon]
   Bool_t          Muon_isGlobal[100];   //[nMuon]
   Bool_t          Muon_isPFcand[100];   //[nMuon]
   Bool_t          Muon_isTracker[100];   //[nMuon]
   Bool_t          Muon_mediumId[100];   //[nMuon]
   Bool_t          Muon_mediumPromptId[100];   //[nMuon]
   UChar_t         Muon_miniIsoId[100];   //[nMuon]
   UChar_t         Muon_multiIsoId[100];   //[nMuon]
   UChar_t         Muon_mvaId[100];   //[nMuon]
   UChar_t         Muon_pfIsoId[100];   //[nMuon]
   Bool_t          Muon_softId[100];   //[nMuon]
   Bool_t          Muon_softMvaId[100];   //[nMuon]
   Bool_t          Muon_tightId[100];   //[nMuon]
   UChar_t         Muon_tkIsoId[100];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[100];   //[nMuon]
    UInt_t          nPhoton;// */
    Float_t         Photon_eCorr[100];   //[nPhoton] */
    Float_t         Photon_energyErr[100];   //[nPhoton] */
    Float_t         Photon_eta[100];   //[nPhoton] */
    Float_t         Photon_hoe[100];   //[nPhoton] */
    Float_t         Photon_mass[100];   //[nPhoton] */
    Float_t         Photon_mvaID[100];   //[nPhoton] */
    Float_t         Photon_mvaID17[100];   //[nPhoton] */
    Float_t         Photon_pfRelIso03_all[100];   //[nPhoton] */
    Float_t         Photon_pfRelIso03_chg[100];   //[nPhoton] */
    Float_t         Photon_phi[100];   //[nPhoton] */
    Float_t         Photon_pt[100];   //[nPhoton] */
    Float_t         Photon_r9[100];   //[nPhoton] */
    Float_t         Photon_sieie[100];   //[nPhoton] */
    Int_t           Photon_charge[100];   //[nPhoton] */
    Int_t           Photon_cutBased[100];   //[nPhoton] */
    Int_t           Photon_cutBased17Bitmap[100];   //[nPhoton] */
    Int_t           Photon_electronIdx[100];   //[nPhoton] */
    Int_t           Photon_jetIdx[100];   //[nPhoton] */
    Int_t           Photon_pdgId[100];   //[nPhoton] */
    Int_t           Photon_vidNestedWPBitmap[100];   //[nPhoton] */
    Bool_t          Photon_electronVeto[100];   //[nPhoton] */
    Bool_t          Photon_isScEtaEB[100];   //[nPhoton] */
    Bool_t          Photon_isScEtaEE[100];   //[nPhoton] */
    Bool_t          Photon_mvaID17_WP80[100];   //[nPhoton]
    Bool_t          Photon_mvaID17_WP90[100];   //[nPhoton]
    Bool_t          Photon_mvaID_WP80[100];   //[nPhoton]
    Bool_t          Photon_mvaID_WP90[100];   //[nPhoton]
    Bool_t          Photon_pixelSeed[100];   //[nPhoton] 
   Float_t         Pileup_nTrueInt;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nTau;
   Float_t         Tau_chargedIso[100];   //[nTau]
   Float_t         Tau_dxy[100];   //[nTau]
   Float_t         Tau_dz[100];   //[nTau]
   Float_t         Tau_eta[100];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[100];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[100];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[100];   //[nTau]
   Float_t         Tau_mass[100];   //[nTau]
   Float_t         Tau_neutralIso[100];   //[nTau]
   Float_t         Tau_phi[100];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[100];   //[nTau]
   Float_t         Tau_pt[100];   //[nTau]
   Float_t         Tau_puCorr[100];   //[nTau]
   Float_t         Tau_rawAntiEle[100];   //[nTau]
   Float_t         Tau_rawIso[100];   //[nTau]
   Float_t         Tau_rawIsodR03[100];   //[nTau]
   Float_t         Tau_rawMVAnewDM2017v2[100];   //[nTau]
   Float_t         Tau_rawMVAoldDM[100];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v1[100];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v2[100];   //[nTau]
   Float_t         Tau_rawMVAoldDMdR032017v2[100];   //[nTau]
   Int_t           Tau_charge[100];   //[nTau]
   Int_t           Tau_decayMode[100];   //[nTau]
   Int_t           Tau_jetIdx[100];   //[nTau]
   Int_t           Tau_rawAntiEleCat[100];   //[nTau]
   UChar_t         Tau_idAntiEle[100];   //[nTau]
   UChar_t         Tau_idAntiMu[100];   //[nTau]
   Bool_t          Tau_idDecayMode[100];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[100];   //[nTau]
   UChar_t         Tau_idDeepTau2017v2p1VSmu[100];
   UChar_t         Tau_idDeepTau2017v2p1VSe[100];
   UChar_t         Tau_idDeepTau2017v2p1VSjet[100];
   UChar_t         Tau_idMVAnewDM2017v2[100];   //[nTau]
   UChar_t         Tau_idMVAoldDM[100];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v1[100];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v2[100];   //[nTau]
   UChar_t         Tau_idMVAoldDMdR032017v2[100];   //[nTau]
   UInt_t          nGenVisTau;	
   Float_t         GenVisTau_eta[100];   //[nGenVisTau]
   Float_t         GenVisTau_mass[100];   //[nGenVisTau]
   Float_t         GenVisTau_phi[100];//[nGenVisTau]
   Float_t         GenVisTau_pt[100];   //[nGenVisTau]
   Int_t           GenVisTau_status[100];   //[nGenVisTau]
   Int_t           GenVisTau_charge[100];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[100];   //[nGenVisTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[500];   //[nTrigObj]
   Float_t         TrigObj_eta[500];   //[nTrigObj]
   Float_t         TrigObj_phi[500];   //[nTrigObj]
   Float_t         TrigObj_l1pt[500];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[500];   //[nTrigObj]
   Float_t         TrigObj_l2pt[500];   //[nTrigObj]
   Int_t           TrigObj_id[500];   //[nTrigObj]
   Int_t           TrigObj_l1iso[500];   //[nTrigObj]
   Int_t           TrigObj_l1charge[500];   //[nTrigObj]
   Int_t           TrigObj_filterBits[500];   //[nTrigObj]
   /* Int_t           genTtbarId; */
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[100];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[50];   //[nSV]
   Float_t         SV_dlenSig[50];   //[nSV]
   Float_t         SV_pAngle[50];   //[nSV]
   Int_t           Electron_genPartIdx[100];   //[nElectron]
   UChar_t         Electron_genPartFlav[100];   //[nElectron]
   /* Int_t           GenJetAK8_partonFlavour[50];   //[nGenJetAK8] */
   /* UChar_t         GenJetAK8_hadronFlavour[50];   //[nGenJetAK8] */
   Int_t           GenJet_partonFlavour[50];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[50];   //[nGenJet]
   Int_t           Jet_genJetIdx[100];   //[nJet]
   Int_t           Jet_hadronFlavour[100];   //[nJet]
   Int_t           Jet_partonFlavour[100];   //[nJet]
   Int_t           Muon_genPartIdx[100];   //[nMuon]
   UChar_t         Muon_genPartFlav[100];   //[nMuon]
   /* Int_t           Photon_genPartIdx[100];   //[nPhoton] */
   /* UChar_t         Photon_genPartFlav[100];   //[nPhoton] */
   /* Float_t         MET_fiducialGenPhi; */
   /* Float_t         MET_fiducialGenPt; */
   /* UChar_t         Electron_cleanmask[100];   //[nElectron] */
   /* UChar_t         Jet_cleanmask[100];   //[nJet] */
   /* UChar_t         Muon_cleanmask[100];   //[nMuon] */
   /* UChar_t         Photon_cleanmask[100];   //[nPhoton] */
   /* UChar_t         Tau_cleanmask[100];   //[nTau] */
   Float_t         SV_chi2[10];   //[nSV]
   Float_t         SV_eta[10];   //[nSV]
   Float_t         SV_mass[10];   //[nSV]
   Float_t         SV_ndof[10];   //[nSV]
   Float_t         SV_phi[10];   //[nSV]
   Float_t         SV_pt[10];   //[nSV]
   Float_t         SV_x[10];   //[nSV]
   Float_t         SV_y[10];   //[nSV]
   Float_t         SV_z[10];   //[nSV]
   Int_t           Tau_genPartIdx[100];   //[nTau]
   UChar_t         Tau_genPartFlav[100];   //[nTau]
   /* Bool_t          HLT_Ele27_WPLoose_Gsf; */
   /* Bool_t          HLT_Ele27_WPLoose_Gsf_WHbbBoost; */
   Bool_t          HLT_Ele27_WPTight_Gsf;
   /* Bool_t          HLT_Ele27_WPTight_Gsf_L1JetTauSeeded; */
   /* Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf; */
   /* Bool_t          HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1; */
   /* Bool_t          HLT_Ele27_eta2p1_WPTight_Gsf; */
   /* Bool_t          HLT_Ele30_WPTight_Gsf; */
   /* Bool_t          HLT_Ele30_eta2p1_WPLoose_Gsf; */
   /* Bool_t          HLT_Ele30_eta2p1_WPTight_Gsf; */
   Bool_t          HLT_Ele32_WPTight_Gsf;
   Bool_t          HLT_Ele32_WPTight_Gsf_L1DoubleEG;
   /* Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf; */
   /* Bool_t          HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1; */
   /* Bool_t          HLT_Ele32_eta2p1_WPTight_Gsf; */
   /* Bool_t          HLT_Ele35_WPLoose_Gsf; */
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu27;
   Bool_t          HLT_IsoTkMu24;
   /* Bool_t          HLT_IsoTkMu27; */
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   /* TBranch        *b_btagWeight_CSVV2;   //! */
   /* TBranch        *b_btagWeight_CMVA;   //! */
   /* TBranch        *b_CaloMET_phi;   //! */
   /* TBranch        *b_CaloMET_pt;   //! */
   /* TBranch        *b_CaloMET_sumEt;   //! */
   /* TBranch        *b_ChsMET_phi;   //! */
   /* TBranch        *b_ChsMET_pt;   //! */
   /* TBranch        *b_ChsMET_sumEt;   //! */
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   /* TBranch        *b_Electron_mvaFall17V1Iso;   //! */
   /* TBranch        *b_Electron_mvaFall17V1noIso;   //! */
   /* TBranch        *b_Electron_mvaFall17V2Iso;   //! */
   /* TBranch        *b_Electron_mvaFall17V2noIso;   //! */
   /* TBranch        *b_Electron_mvaSpring16GP;   //! */
   /* TBranch        *b_Electron_mvaSpring16HZZ;   //! */
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   /* TBranch        *b_Electron_cutBased_Fall17_V1;   //! */
   /* TBranch        *b_Electron_cutBased_HLTPreSel;   //! */
   /* TBranch        *b_Electron_cutBased_Spring15;   //! */
   /* TBranch        *b_Electron_cutBased_Sum16;   //! */
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   /* TBranch        *b_Electron_vidNestedWPBitmap;   //! */
   /* TBranch        *b_Electron_vidNestedWPBitmapSpring15;   //! */
   /* TBranch        *b_Electron_vidNestedWPBitmapSum16;   //! */
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_lostHits;   //!
   /* TBranch        *b_Electron_mvaFall17V1Iso_WP80;   //! */
   /* TBranch        *b_Electron_mvaFall17V1Iso_WP90;   //! */
   /* TBranch        *b_Electron_mvaFall17V1Iso_WPL;   //! */
   /* TBranch        *b_Electron_mvaFall17V1noIso_WP80;   //! */
   /* TBranch        *b_Electron_mvaFall17V1noIso_WP90;   //! */
   /* TBranch        *b_Electron_mvaFall17V1noIso_WPL;   //! */
   /* TBranch        *b_Electron_mvaFall17V2Iso_WP80;   //! */
   /* TBranch        *b_Electron_mvaFall17V2Iso_WP90;   //! */
   /* TBranch        *b_Electron_mvaFall17V2Iso_WPL;   //! */
   /* TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //! */
   /* TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //! */
   /* TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //! */
   /* TBranch        *b_Electron_mvaSpring16GP_WP80;   //! */
   /* TBranch        *b_Electron_mvaSpring16GP_WP90;   //! */
   /* TBranch        *b_Electron_mvaSpring16HZZ_WPL;   //! */
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_genWeight;   //!
   /* TBranch        *b_LHEWeight_originalXWGTUP;   //! */
   /* TBranch        *b_nLHEPdfWeight;   //! */
   /* TBranch        *b_LHEPdfWeight;   //! */
   /* TBranch        *b_nLHEScaleWeight;   //! */
   /* TBranch        *b_LHEScaleWeight;   //! */
   /* TBranch        *b_nPSWeight;   //! */
   /* TBranch        *b_PSWeight;   //! */
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepC;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_bRegCorr;   //!
   TBranch        *b_Jet_bRegRes;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_puId;   //!
   /* TBranch        *b_GenMET_phi;   //! */
   /* TBranch        *b_GenMET_pt;   //! */
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_METFixEE2017_phi;   //!
   TBranch        *b_METFixEE2017_pt;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
    TBranch        *b_nPhoton;   //! */
    TBranch        *b_Photon_eCorr;   //! */
    TBranch        *b_Photon_energyErr;   //! */
    TBranch        *b_Photon_eta;   //! */
    TBranch        *b_Photon_hoe;   //! */
    TBranch        *b_Photon_mass;   //! */
    TBranch        *b_Photon_mvaID;   //! */
    TBranch        *b_Photon_mvaID17;   //! */
    TBranch        *b_Photon_pfRelIso03_all;   //! */
    TBranch        *b_Photon_pfRelIso03_chg;   //! */
    TBranch        *b_Photon_phi;   //! */
    TBranch        *b_Photon_pt;   //! */
    TBranch        *b_Photon_r9;   //! */
    TBranch        *b_Photon_sieie;   //! */
    TBranch        *b_Photon_charge;   //! */
    TBranch        *b_Photon_cutBased;   //! */
    TBranch        *b_Photon_cutBased17Bitmap;   //! */
    TBranch        *b_Photon_electronIdx;   //! */
    TBranch        *b_Photon_jetIdx;   //! */
    TBranch        *b_Photon_pdgId;   //! */
    TBranch        *b_Photon_vidNestedWPBitmap;   //! */
    TBranch        *b_Photon_electronVeto;   //! */
    TBranch        *b_Photon_isScEtaEB;   //! */
    TBranch        *b_Photon_isScEtaEE;   //! */
    TBranch        *b_Photon_mvaID17_WP80;   //! */
    TBranch        *b_Photon_mvaID17_WP90;   //! */
    TBranch        *b_Photon_mvaID_WP80;   //! */
    TBranch        *b_Photon_mvaID_WP90;   //! */
    TBranch        *b_Photon_pixelSeed;   //! */
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_rawAntiEle;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_rawMVAnewDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDM;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v1;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_rawAntiEleCat;   //!
   TBranch        *b_Tau_idAntiEle;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayMode;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idMVAnewDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDM;   //!
  TBranch        *b_Tau_idDeepTau2017v2p1VSmu;  //!
  TBranch        *b_Tau_idDeepTau2017v2p1VSe;  //!
  TBranch        *b_Tau_idDeepTau2017v2p1VSjet;  //!    
   TBranch        *b_Tau_idMVAoldDM2017v1;   //!
   TBranch        *b_Tau_idMVAoldDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_nGenVisTau;	
   TBranch        *b_GenVisTau_eta;  
   TBranch        *b_GenVisTau_mass;  
   TBranch        *b_GenVisTau_phi;
   TBranch        *b_GenVisTau_pt;   
   TBranch        *b_GenVisTau_status;
   TBranch        *b_GenVisTau_charge; 
   TBranch        *b_GenVisTau_genPartIdxMother;
   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   /* TBranch        *b_genTtbarId;   //! */
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   /* TBranch        *b_GenJetAK8_partonFlavour;   //! */
   /* TBranch        *b_GenJetAK8_hadronFlavour;   //! */
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   /* TBranch        *b_Photon_genPartIdx;   //! */
   /* TBranch        *b_Photon_genPartFlav;   //! */
   /* TBranch        *b_MET_fiducialGenPhi;   //! */
   /* TBranch        *b_MET_fiducialGenPt;   //! */
   /* TBranch        *b_Electron_cleanmask;   //! */
   /* TBranch        *b_Jet_cleanmask;   //! */
   /* TBranch        *b_Muon_cleanmask;   //! */
   /* TBranch        *b_Photon_cleanmask;   //! */
   /* TBranch        *b_Tau_cleanmask;   //! */
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_Tau_genPartIdx;   //!
   TBranch        *b_Tau_genPartFlav;   //!
   /* TBranch        *b_HLT_Ele27_WPLoose_Gsf;   //! */
   /* TBranch        *b_HLT_Ele27_WPLoose_Gsf_WHbbBoost;   //! */
   TBranch        *b_HLT_Ele27_WPTight_Gsf;   //!
   /* TBranch        *b_HLT_Ele27_WPTight_Gsf_L1JetTauSeeded;   //! */
   /* TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf;   //! */
   /* TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //! */
   /* TBranch        *b_HLT_Ele27_eta2p1_WPTight_Gsf;   //! */
   /* TBranch        *b_HLT_Ele30_WPTight_Gsf;   //! */
   /* TBranch        *b_HLT_Ele30_eta2p1_WPLoose_Gsf;   //! */
   /* TBranch        *b_HLT_Ele30_eta2p1_WPTight_Gsf;   //! */
   TBranch        *b_HLT_Ele32_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_L1DoubleEG;   //!
   /* TBranch        *b_HLT_Ele32_eta2p1_WPLoose_Gsf;   //! */
   /* TBranch        *b_HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1;   //! */
   /* TBranch        *b_HLT_Ele32_eta2p1_WPTight_Gsf;   //! */
   /* TBranch        *b_HLT_Ele35_WPLoose_Gsf;   //! */
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoTkMu24;   //!
   /* TBranch        *b_HLT_IsoTkMu27;   //! */
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!

   RHAna(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~RHAna() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //User added functions

  int ReadLimited(int level, Long64_t entry); //This is an important function.
  // Now here are some functions to set some flags
  void SetHstFileName(const char *HstFileName){ _HstFileName = HstFileName;}
  void SetSumFileName(const char *SumFileName){ _SumFileName = SumFileName;}
  void SetSample(int sample){_sample=sample;}
  void SetVerbose(int verbose){ _verbosity = verbose; }
  void SetData(int data){_data=data;}
  void SetYear(int year){_year = year;}
  void SetEra(TString era){_era=era;}
  void BookHistograms();
  bool pass_tau_cuts(int opt, int tin);
  void Sort(int opt);
  float delta_phi(float phi1, float phi2);
  bool TaulepCleaning(TLorentzVector t);
  int GenMother(int ind, int mom_ind);
  
public:
  struct Hists {
    TH1F *llep;
    TH1F *muons[5];
    TH1F *electrons[5];
    TH1F *emu[5];
    TH1F *taus[5];
    TH1F *isotrack[5];
    TH1F *jets[5];
    TH1F *photons[5];
    TH1F *dphi[8];
    TH1F *dR[8];
    TH2F *two[10];
    TH1F *met[8];
    TH1F *ht;
    TH1F *Genmu[5];
    TH1F *Genele[5];
    TH1F *Genemu[5];
    TH1F *Gentau[5];
    TH1F *mom[3];
    TH1F *isoM[5];	
    TH1F *isoE[5];		
    TH1F *btag[5];	
    TH1F *lep4[5];
    TH1F *mT[5];
    TH1F *impact[10];
  };
  
  struct Jet{
    int id;
    TLorentzVector v;
    int ind;
  };

  struct Photon{
    int id;
    TLorentzVector v;
    int ind;
  };

  struct Lepton {
    TLorentzVector v;
    int id; int ind;
    float wt;
    int flavor; int charge;
    bool lepcleaning; bool taucleaning;
    int momid;
    int genmatch;
    int jetmatch;
  };
  
 protected:
  Hists h;

 private:
  //Global variables go here. Make them global only if necessary.
  TFile *_HstFile;
  const char *_HstFileName;
  const char *_SumFileName;
  int _verbosity,_exclude;
  int nEvtTotal,nEvtRan;
  int _data, _lep, _year, _sample;
  bool GoodEvt, GoodEvt2016, GoodEvt2017, GoodEvt2018;
  float metpt, metphi;
  TString _era;
  // Here are arrays in which we store our objects
  vector<Lepton> goodMu, llep, taus, loosetaus, goodEle, allEle, allTau, goodMu_pos, goodMu_neg, GenMu, GenEle, isotracks,gentaus; 
  vector<Jet> goodJet, allJet, bJet;
  vector<Photon> goodPhoton, allPhoton;
  
   ClassDef(RHAna,0);
};

#endif
//A collection of trees constitute a chain though they are technically same objects but a tree comprises many branches but they are structurally different objects.
#ifdef RHAna_cxx
void RHAna::Init(TTree *tree)
{
  /* cout<<"Inside Init()"<<endl; */
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   /* fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi); */
   /* fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt); */
   /* fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt); */
   /* fChain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi); */
   /* fChain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt); */
   /* fChain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt); */
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);//Connects the nelectron from skim.root to the variable nelectron and 	  									its branch			
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   /* fChain->SetBranchAddress("Electron_mvaFall17V1Iso", Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1noIso", Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso); */
   /* fChain->SetBranchAddress("Electron_mvaSpring16GP", Electron_mvaSpring16GP, &b_Electron_mvaSpring16GP); */
   /* fChain->SetBranchAddress("Electron_mvaSpring16HZZ", Electron_mvaSpring16HZZ, &b_Electron_mvaSpring16HZZ); */
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   /* fChain->SetBranchAddress("Electron_cutBased_Fall17_V1", Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1); */
   /* fChain->SetBranchAddress("Electron_cutBased_HLTPreSel", Electron_cutBased_HLTPreSel, &b_Electron_cutBased_HLTPreSel); */
   /* fChain->SetBranchAddress("Electron_cutBased_Spring15", Electron_cutBased_Spring15, &b_Electron_cutBased_Spring15); */
   /* fChain->SetBranchAddress("Electron_cutBased_Sum16", Electron_cutBased_Sum16, &b_Electron_cutBased_Sum16); */
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   /* fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap); */
   /* fChain->SetBranchAddress("Electron_vidNestedWPBitmapSpring15", Electron_vidNestedWPBitmapSpring15, &b_Electron_vidNestedWPBitmapSpring15); */
   /* fChain->SetBranchAddress("Electron_vidNestedWPBitmapSum16", Electron_vidNestedWPBitmapSum16, &b_Electron_vidNestedWPBitmapSum16); */
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   /* fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", Electron_mvaFall17V1Iso_WP80, &b_Electron_mvaFall17V1Iso_WP80); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", Electron_mvaFall17V1Iso_WP90, &b_Electron_mvaFall17V1Iso_WP90); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", Electron_mvaFall17V1Iso_WPL, &b_Electron_mvaFall17V1Iso_WPL); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90); */
   /* fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL); */
   /* fChain->SetBranchAddress("Electron_mvaSpring16GP_WP80", Electron_mvaSpring16GP_WP80, &b_Electron_mvaSpring16GP_WP80); */
   /* fChain->SetBranchAddress("Electron_mvaSpring16GP_WP90", Electron_mvaSpring16GP_WP90, &b_Electron_mvaSpring16GP_WP90); */
   /* fChain->SetBranchAddress("Electron_mvaSpring16HZZ_WPL", Electron_mvaSpring16HZZ_WPL, &b_Electron_mvaSpring16HZZ_WPL); */
   fChain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   fChain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   fChain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   fChain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   fChain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   fChain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   fChain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   fChain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   fChain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   fChain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   fChain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   fChain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   fChain->SetBranchAddress("Tau_rawMVAnewDM2017v2", Tau_rawMVAnewDM2017v2, &b_Tau_rawMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v1", Tau_rawMVAoldDM2017v1, &b_Tau_rawMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_rawMVAoldDM2017v2", Tau_rawMVAoldDM2017v2, &b_Tau_rawMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", Tau_rawMVAoldDMdR032017v2, &b_Tau_rawMVAoldDMdR032017v2);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSmu",Tau_idDeepTau2017v2p1VSmu, &b_Tau_idDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSe",Tau_idDeepTau2017v2p1VSe, &b_Tau_idDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("Tau_idDeepTau2017v2p1VSjet",Tau_idDeepTau2017v2p1VSjet, &b_Tau_idDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("Tau_idMVAnewDM2017v2", Tau_idMVAnewDM2017v2, &b_Tau_idMVAnewDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v1", Tau_idMVAoldDM2017v1, &b_Tau_idMVAoldDM2017v1);
   fChain->SetBranchAddress("Tau_idMVAoldDM2017v2", Tau_idMVAoldDM2017v2, &b_Tau_idMVAoldDM2017v2);
   fChain->SetBranchAddress("Tau_idMVAoldDMdR032017v2", Tau_idMVAoldDMdR032017v2, &b_Tau_idMVAoldDMdR032017v2);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   fChain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
    fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton); //*/
    fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr); //*/
    fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr); //*/
    fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta); //*/
    fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe); //*/
    fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass); //*/
    fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID); //*/
    // fChain->SetBranchAddress("Photon_mvaID17", Photon_mvaID17, &b_Photon_mvaID17);// */
    fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all); //*/
    fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg); //*/
      fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);// */
    fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt); //*/
    fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9); //*/
    fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie); //*/
    fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge); //*/
    // fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased); //*/
    // fChain->SetBranchAddress("Photon_cutBased17Bitmap", Photon_cutBased17Bitmap, &b_Photon_cutBased17Bitmap); //*/
      fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);// */
      fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);// */
      fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);// */
      fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);// */
      fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);// */
    fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB); //*/
    fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE); //*/
    //fChain->SetBranchAddress("Photon_mvaID17_WP80", Photon_mvaID17_WP80, &b_Photon_mvaID17_WP80); //*/
    //fChain->SetBranchAddress("Photon_mvaID17_WP90", Photon_mvaID17_WP90, &b_Photon_mvaID17_WP90);// */
    fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80); //*/
    fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90); //*/
      fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);// */
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   fChain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
  fChain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
  fChain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
  fChain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   /* fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask); */
   /* fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask); */
   /* fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask); */
   /* fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask); */
   /* fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask); */
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);

   //Following branches are not available in data
   if(_data==0){
     /* fChain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2); */
     /* fChain->SetBranchAddress("btagWeight_CMVA", &btagWeight_CMVA, &b_btagWeight_CMVA); */
     fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
     fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
     fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
     fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
     fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
     fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
     fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
     fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
     fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
     fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
     fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
     fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
     fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
     fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
     /* fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar); */
     /* fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF); */
     /* fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight); */
     /* fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1); */
     /* fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2); */
     /* fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1); */
     /* fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2); */
     /* fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1); */
     /* fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2); */
     /* fChain->SetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP, &b_LHEWeight_originalXWGTUP); */
     /* fChain->SetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight, &b_nLHEPdfWeight); */
     /* fChain->SetBranchAddress("LHEPdfWeight", LHEPdfWeight, &b_LHEPdfWeight); */
     /* fChain->SetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight, &b_nLHEScaleWeight); */
     /* fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight); */
     /* fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight); */
     /* fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight); */
     /* fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi); */
     /* fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt); */
     fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
     fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
     fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
     fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
     fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
     fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
     fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
     fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
     fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
     fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
     fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
     fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
     fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
     /* fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx); */
     /* fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav); */
     /* fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi); */
     /* fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt); */
     fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
     fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
     /* fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step); */
   }
   if(_year==2017){
     fChain->SetBranchAddress("METFixEE2017_phi", &METFixEE2017_phi, &b_METFixEE2017_phi);
     fChain->SetBranchAddress("METFixEE2017_pt", &METFixEE2017_pt, &b_METFixEE2017_pt);
     if(_lep==0 && _era!="B" && _era!="C" && _data==1)
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
     if(_data==1)
       fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   }
   if(_year==2016)
     fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
   if(_year==2018){
     fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf, &b_HLT_Ele32_WPTight_Gsf);
     fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b_HLT_Ele32_WPTight_Gsf_L1DoubleEG);
   }
}

Bool_t RHAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
int RHAna::ReadLimited(int level, Long64_t entry)
{
  /* cout<<"Inside ReadLimited()"<<endl; */
  if(level==0){ GetEntry(entry); return 1; }
  if(level==1){
    //Turn on only those branches that we need
    b_run->GetEntry(entry);
    b_luminosityBlock->GetEntry(entry);
    b_event->GetEntry(entry);
    b_nElectron->GetEntry(entry);
    b_Electron_deltaEtaSC->GetEntry(entry);
    b_Electron_dr03EcalRecHitSumEt->GetEntry(entry);
    b_Electron_dr03HcalDepth1TowerSumEt->GetEntry(entry);
    b_Electron_dr03TkSumPt->GetEntry(entry);
    b_Electron_dr03TkSumPtHEEP->GetEntry(entry);
    b_Electron_dxy->GetEntry(entry);
    b_Electron_dxyErr->GetEntry(entry);
    b_Electron_dz->GetEntry(entry);
    b_Electron_dzErr->GetEntry(entry);
    b_Electron_eCorr->GetEntry(entry);
    b_Electron_eInvMinusPInv->GetEntry(entry);
    b_Electron_energyErr->GetEntry(entry);
    b_Electron_eta->GetEntry(entry);
    b_Electron_hoe->GetEntry(entry);
    b_Electron_ip3d->GetEntry(entry);
    b_Electron_jetRelIso->GetEntry(entry);
    b_Electron_mass->GetEntry(entry);
    b_Electron_miniPFRelIso_all->GetEntry(entry);
    b_Electron_miniPFRelIso_chg->GetEntry(entry);
    /* b_Electron_mvaFall17V1Iso->GetEntry(entry); */
    /* b_Electron_mvaFall17V1noIso->GetEntry(entry); */
    /* b_Electron_mvaFall17V2Iso->GetEntry(entry); */
    /* b_Electron_mvaFall17V2noIso->GetEntry(entry); */
    /* b_Electron_mvaSpring16GP->GetEntry(entry); */
    /* b_Electron_mvaSpring16HZZ->GetEntry(entry); */
    b_Electron_pfRelIso03_all->GetEntry(entry);
    b_Electron_pfRelIso03_chg->GetEntry(entry);
    b_Electron_phi->GetEntry(entry);
    b_Electron_pt->GetEntry(entry);
    b_Electron_r9->GetEntry(entry);
    b_Electron_sieie->GetEntry(entry);
    b_Electron_sip3d->GetEntry(entry);
    b_Electron_mvaTTH->GetEntry(entry);
    b_Electron_charge->GetEntry(entry);
    b_Electron_cutBased->GetEntry(entry);
    /* b_Electron_cutBased_Fall17_V1->GetEntry(entry); */
    /* b_Electron_cutBased_HLTPreSel->GetEntry(entry); */
    /* b_Electron_cutBased_Spring15->GetEntry(entry); */
    /* b_Electron_cutBased_Sum16->GetEntry(entry); */
    b_Electron_jetIdx->GetEntry(entry);
    b_Electron_pdgId->GetEntry(entry);
    b_Electron_photonIdx->GetEntry(entry);
    b_Electron_tightCharge->GetEntry(entry);
    /* b_Electron_vidNestedWPBitmap->GetEntry(entry); */
    /* b_Electron_vidNestedWPBitmapSpring15->GetEntry(entry); */
    /* b_Electron_vidNestedWPBitmapSum16->GetEntry(entry); */
    b_Electron_convVeto->GetEntry(entry);
    b_Electron_cutBased_HEEP->GetEntry(entry);
    b_Electron_isPFcand->GetEntry(entry);
    b_Electron_lostHits->GetEntry(entry);
    /* b_Electron_mvaFall17V1Iso_WP80->GetEntry(entry); */
    /* b_Electron_mvaFall17V1Iso_WP90->GetEntry(entry); */
    /* b_Electron_mvaFall17V1Iso_WPL->GetEntry(entry); */
    /* b_Electron_mvaFall17V1noIso_WP80->GetEntry(entry); */
    /* b_Electron_mvaFall17V1noIso_WP90->GetEntry(entry); */
    /* b_Electron_mvaFall17V1noIso_WPL->GetEntry(entry); */
    /* b_Electron_mvaFall17V2Iso_WP80->GetEntry(entry); */
    /* b_Electron_mvaFall17V2Iso_WP90->GetEntry(entry); */
    /* b_Electron_mvaFall17V2Iso_WPL->GetEntry(entry); */
    /* b_Electron_mvaFall17V2noIso_WP80->GetEntry(entry); */
    /* b_Electron_mvaFall17V2noIso_WP90->GetEntry(entry); */
    /* b_Electron_mvaFall17V2noIso_WPL->GetEntry(entry); */
    /* b_Electron_mvaSpring16GP_WP80->GetEntry(entry); */
    /* b_Electron_mvaSpring16GP_WP90->GetEntry(entry); */
    /* b_Electron_mvaSpring16HZZ_WPL->GetEntry(entry); */
    // b_nIsoTrack->GetEntry(entry);
    // b_IsoTrack_dxy->GetEntry(entry);
    // b_IsoTrack_dz->GetEntry(entry);
    // b_IsoTrack_eta->GetEntry(entry);
    // b_IsoTrack_pfRelIso03_all->GetEntry(entry);
    // b_IsoTrack_pfRelIso03_chg->GetEntry(entry);
    // b_IsoTrack_phi->GetEntry(entry);
    // b_IsoTrack_pt->GetEntry(entry);
    // b_IsoTrack_miniPFRelIso_all->GetEntry(entry);
    // b_IsoTrack_miniPFRelIso_chg->GetEntry(entry);
    // b_IsoTrack_fromPV->GetEntry(entry);
    // b_IsoTrack_pdgId->GetEntry(entry);
    // b_IsoTrack_isHighPurityTrack->GetEntry(entry);
    // b_IsoTrack_isPFcand->GetEntry(entry);
    // b_IsoTrack_isFromLostTrack->GetEntry(entry);
    b_nJet->GetEntry(entry);
    b_Jet_area->GetEntry(entry);
    b_Jet_btagCMVA->GetEntry(entry);
    b_Jet_btagCSVV2->GetEntry(entry);
    b_Jet_btagDeepB->GetEntry(entry);
    b_Jet_btagDeepC->GetEntry(entry);
    b_Jet_btagDeepFlavB->GetEntry(entry);
    b_Jet_chEmEF->GetEntry(entry);
    b_Jet_chHEF->GetEntry(entry);
    b_Jet_eta->GetEntry(entry);
    b_Jet_mass->GetEntry(entry);
    b_Jet_muEF->GetEntry(entry);
    b_Jet_neEmEF->GetEntry(entry);
    b_Jet_neHEF->GetEntry(entry);
    b_Jet_phi->GetEntry(entry);
    b_Jet_pt->GetEntry(entry);
    b_Jet_qgl->GetEntry(entry);
    b_Jet_rawFactor->GetEntry(entry);
    b_Jet_bRegCorr->GetEntry(entry);
    b_Jet_bRegRes->GetEntry(entry);
    b_Jet_electronIdx1->GetEntry(entry);
    b_Jet_electronIdx2->GetEntry(entry);
    b_Jet_jetId->GetEntry(entry);
    b_Jet_muonIdx1->GetEntry(entry);
    b_Jet_muonIdx2->GetEntry(entry);
    b_Jet_nConstituents->GetEntry(entry);
    b_Jet_nElectrons->GetEntry(entry);
    b_Jet_nMuons->GetEntry(entry);
    b_Jet_puId->GetEntry(entry);
    b_MET_MetUnclustEnUpDeltaX->GetEntry(entry);
    b_MET_MetUnclustEnUpDeltaY->GetEntry(entry);
    b_MET_phi->GetEntry(entry);
    b_MET_pt->GetEntry(entry);
    b_MET_sumEt->GetEntry(entry);
    b_nMuon->GetEntry(entry);
    b_Muon_dxy->GetEntry(entry);
    b_Muon_dxyErr->GetEntry(entry);
    b_Muon_dz->GetEntry(entry);
    b_Muon_dzErr->GetEntry(entry);
    b_Muon_eta->GetEntry(entry);
    b_Muon_ip3d->GetEntry(entry);
    b_Muon_jetRelIso->GetEntry(entry);
    b_Muon_mass->GetEntry(entry);
    b_Muon_miniPFRelIso_all->GetEntry(entry);
    b_Muon_miniPFRelIso_chg->GetEntry(entry);
    b_Muon_pfRelIso03_all->GetEntry(entry);
    b_Muon_pfRelIso03_chg->GetEntry(entry);
    b_Muon_pfRelIso04_all->GetEntry(entry);
    b_Muon_phi->GetEntry(entry);
    b_Muon_pt->GetEntry(entry);
    b_Muon_ptErr->GetEntry(entry);
    b_Muon_segmentComp->GetEntry(entry);
    b_Muon_sip3d->GetEntry(entry);
    b_Muon_mvaTTH->GetEntry(entry);
    b_Muon_charge->GetEntry(entry);
    b_Muon_jetIdx->GetEntry(entry);
    b_Muon_nStations->GetEntry(entry);
    b_Muon_nTrackerLayers->GetEntry(entry);
    b_Muon_pdgId->GetEntry(entry);
    b_Muon_tightCharge->GetEntry(entry);
    b_Muon_highPtId->GetEntry(entry);
    b_Muon_inTimeMuon->GetEntry(entry);
    b_Muon_isGlobal->GetEntry(entry);
    b_Muon_isPFcand->GetEntry(entry);
    b_Muon_isTracker->GetEntry(entry);
    b_Muon_mediumId->GetEntry(entry);
    b_Muon_mediumPromptId->GetEntry(entry);
    b_Muon_miniIsoId->GetEntry(entry);
    b_Muon_multiIsoId->GetEntry(entry);
    b_Muon_mvaId->GetEntry(entry);
    b_Muon_pfIsoId->GetEntry(entry);
    b_Muon_softId->GetEntry(entry);
    b_Muon_softMvaId->GetEntry(entry);
    b_Muon_tightId->GetEntry(entry);
    b_Muon_tkIsoId->GetEntry(entry);
    b_Muon_triggerIdLoose->GetEntry(entry);
     b_nPhoton->GetEntry(entry);// */
     b_Photon_eCorr->GetEntry(entry);// */
     b_Photon_energyErr->GetEntry(entry);// */
     b_Photon_eta->GetEntry(entry); //*/
     b_Photon_hoe->GetEntry(entry);// */
     b_Photon_mass->GetEntry(entry); //*/
     b_Photon_mvaID->GetEntry(entry);// */
     b_Photon_mvaID17->GetEntry(entry);// */
     b_Photon_pfRelIso03_all->GetEntry(entry);// */
     b_Photon_pfRelIso03_chg->GetEntry(entry); //*/
     b_Photon_phi->GetEntry(entry); //*/
     b_Photon_pt->GetEntry(entry); //*/
     b_Photon_r9->GetEntry(entry); //*/
     b_Photon_sieie->GetEntry(entry);// */
     b_Photon_charge->GetEntry(entry); //*/
     b_Photon_cutBased->GetEntry(entry); //*/
     b_Photon_cutBased17Bitmap->GetEntry(entry);// */
     b_Photon_electronIdx->GetEntry(entry); //*/
     b_Photon_jetIdx->GetEntry(entry);// */
     b_Photon_pdgId->GetEntry(entry); //*/
     b_Photon_vidNestedWPBitmap->GetEntry(entry); //*/
     b_Photon_electronVeto->GetEntry(entry); //*/
     b_Photon_isScEtaEB->GetEntry(entry); //*/
     b_Photon_isScEtaEE->GetEntry(entry); //*/
     b_Photon_mvaID17_WP80->GetEntry(entry);// */
     b_Photon_mvaID17_WP90->GetEntry(entry); // */
     b_Photon_mvaID_WP80->GetEntry(entry); //*/
     b_Photon_mvaID_WP90->GetEntry(entry); //*/
     b_Photon_pixelSeed->GetEntry(entry);// */
    b_PuppiMET_phi->GetEntry(entry);
    b_PuppiMET_pt->GetEntry(entry);
    b_PuppiMET_sumEt->GetEntry(entry);
    b_RawMET_phi->GetEntry(entry);
    b_RawMET_pt->GetEntry(entry);
    b_RawMET_sumEt->GetEntry(entry);
    b_nTau->GetEntry(entry);
    b_Tau_chargedIso->GetEntry(entry);
    b_Tau_dxy->GetEntry(entry);
    b_Tau_dz->GetEntry(entry);
    b_Tau_eta->GetEntry(entry);
    b_Tau_leadTkDeltaEta->GetEntry(entry);
    b_Tau_leadTkDeltaPhi->GetEntry(entry);
    b_Tau_leadTkPtOverTauPt->GetEntry(entry);
    b_Tau_mass->GetEntry(entry);
    b_Tau_neutralIso->GetEntry(entry);
    b_Tau_phi->GetEntry(entry);
    b_Tau_photonsOutsideSignalCone->GetEntry(entry);
    b_Tau_pt->GetEntry(entry);
    b_Tau_puCorr->GetEntry(entry);
    b_Tau_rawAntiEle->GetEntry(entry);
    b_Tau_rawIso->GetEntry(entry);
    b_Tau_rawIsodR03->GetEntry(entry);
    // b_Tau_rawMVAnewDM2017v2->GetEntry(entry);
    // b_Tau_rawMVAoldDM->GetEntry(entry);
    // b_Tau_rawMVAoldDM2017v1->GetEntry(entry);
    // b_Tau_rawMVAoldDM2017v2->GetEntry(entry);
    // b_Tau_rawMVAoldDMdR032017v2->GetEntry(entry);
    b_Tau_charge->GetEntry(entry);
    b_Tau_decayMode->GetEntry(entry);
    b_Tau_jetIdx->GetEntry(entry);
    b_Tau_rawAntiEleCat->GetEntry(entry);
    b_Tau_idAntiEle->GetEntry(entry);
    b_Tau_idAntiMu->GetEntry(entry);
    b_Tau_idDecayMode->GetEntry(entry);
    b_Tau_idDecayModeNewDMs->GetEntry(entry);
    //    b_Tau_idMVAnewDM2017v2->GetEntry(entry);
    b_Tau_idMVAoldDM->GetEntry(entry);
    // b_Tau_idMVAoldDM2017v1->GetEntry(entry);
    // b_Tau_idMVAoldDM2017v2->GetEntry(entry);
    // b_Tau_idMVAoldDMdR032017v2->GetEntry(entry);
    b_Tau_idDeepTau2017v2p1VSmu->GetEntry(entry);
    b_Tau_idDeepTau2017v2p1VSe->GetEntry(entry);
    b_Tau_idDeepTau2017v2p1VSjet->GetEntry(entry);
    b_nTrigObj->GetEntry(entry);
    b_TrigObj_pt->GetEntry(entry);
    b_TrigObj_eta->GetEntry(entry);
    b_TrigObj_phi->GetEntry(entry);
    b_TrigObj_l1pt->GetEntry(entry);
    b_TrigObj_l1pt_2->GetEntry(entry);
    b_TrigObj_l2pt->GetEntry(entry);
    b_TrigObj_id->GetEntry(entry);
    b_TrigObj_l1iso->GetEntry(entry);
    b_TrigObj_l1charge->GetEntry(entry);
    b_TrigObj_filterBits->GetEntry(entry);
    b_nOtherPV->GetEntry(entry);
    b_OtherPV_z->GetEntry(entry);
    b_PV_ndof->GetEntry(entry);
    b_PV_x->GetEntry(entry);
    b_PV_y->GetEntry(entry);
    b_PV_z->GetEntry(entry);
    b_PV_chi2->GetEntry(entry);
    b_PV_score->GetEntry(entry);
    b_PV_npvs->GetEntry(entry);
    b_PV_npvsGood->GetEntry(entry);
    /* b_SV_chi2->GetEntry(entry); */
    /* b_SV_eta->GetEntry(entry); */
    /* b_SV_mass->GetEntry(entry); */
    /* b_SV_ndof->GetEntry(entry); */
    /* b_SV_phi->GetEntry(entry); */
    /* b_SV_pt->GetEntry(entry); */
    /* b_SV_x->GetEntry(entry); */
    /* b_SV_y->GetEntry(entry); */
    /* b_SV_z->GetEntry(entry); */
    b_HLT_IsoMu24->GetEntry(entry);
    b_HLT_IsoMu27->GetEntry(entry);
    b_HLT_Ele27_WPTight_Gsf->GetEntry(entry);
    /* b_Electron_cleanmask->GetEntry(entry); */
    /* b_Jet_cleanmask->GetEntry(entry); */
    /* b_Muon_cleanmask->GetEntry(entry); */
    /* b_Photon_cleanmask->GetEntry(entry); */
    /* b_Tau_cleanmask->GetEntry(entry); */
    /* b_nSV->GetEntry(entry); */
    /* b_SV_dlen->GetEntry(entry); */
    /* b_SV_dlenSig->GetEntry(entry); */
    /* b_SV_pAngle->GetEntry(entry); */
    /* cout<<"common branches are fine"<<endl;     */
    b_Flag_HBHENoiseFilter->GetEntry(entry);   //!
    b_Flag_HBHENoiseIsoFilter->GetEntry(entry);   //!
    b_Flag_CSCTightHaloFilter->GetEntry(entry);   //!
    b_Flag_CSCTightHaloTrkMuUnvetoFilter->GetEntry(entry);   //!
    b_Flag_CSCTightHalo2015Filter->GetEntry(entry);   //!
    b_Flag_globalTightHalo2016Filter->GetEntry(entry);   //!
    b_Flag_globalSuperTightHalo2016Filter->GetEntry(entry);   //!
    b_Flag_HcalStripHaloFilter->GetEntry(entry);   //!
    b_Flag_hcalLaserEventFilter->GetEntry(entry);   //!
    b_Flag_EcalDeadCellTriggerPrimitiveFilter->GetEntry(entry);   //!
    b_Flag_EcalDeadCellBoundaryEnergyFilter->GetEntry(entry);   //!
    b_Flag_goodVertices->GetEntry(entry);   //!
    b_Flag_eeBadScFilter->GetEntry(entry);   //!
    b_Flag_ecalLaserCorrFilter->GetEntry(entry);   //!
    b_Flag_trkPOGFilters->GetEntry(entry);   //!
    b_Flag_chargedHadronTrackResolutionFilter->GetEntry(entry);   //!
    b_Flag_muonBadTrackFilter->GetEntry(entry);   //!
    b_Flag_BadChargedCandidateFilter->GetEntry(entry);   //!
    b_Flag_BadPFMuonFilter->GetEntry(entry);   //!
    b_Flag_BadChargedCandidateSummer16Filter->GetEntry(entry);   //!
    b_Flag_BadPFMuonSummer16Filter->GetEntry(entry);   //!
    b_Flag_trkPOG_manystripclus53X->GetEntry(entry);   //!
    b_Flag_trkPOG_toomanystripclus53X->GetEntry(entry);   //!
    b_Flag_trkPOG_logErrorTooManyClusters->GetEntry(entry);   //!
    b_Flag_METFilters->GetEntry(entry);   //!
    
    //Enable following branches only for MC
    if(_data==0){
      /* b_btagWeight_CSVV2->GetEntry(entry); */
      /* b_btagWeight_CMVA->GetEntry(entry); */
      b_nGenJet->GetEntry(entry);
      b_GenJet_eta->GetEntry(entry);
      b_GenJet_mass->GetEntry(entry);
      b_GenJet_phi->GetEntry(entry);
      b_GenJet_pt->GetEntry(entry);
      b_nGenPart->GetEntry(entry);
      b_GenPart_eta->GetEntry(entry);
      b_GenPart_mass->GetEntry(entry);
      b_GenPart_phi->GetEntry(entry);
      b_GenPart_pt->GetEntry(entry);
      b_GenPart_genPartIdxMother->GetEntry(entry);
      b_GenPart_pdgId->GetEntry(entry);
      b_GenPart_status->GetEntry(entry);
      b_GenPart_statusFlags->GetEntry(entry);
      /* b_Generator_binvar->GetEntry(entry); */
      /* b_Generator_scalePDF->GetEntry(entry); */
      /* b_Generator_weight->GetEntry(entry); */
      /* b_Generator_x1->GetEntry(entry); */
      /* b_Generator_x2->GetEntry(entry); */
      /* b_Generator_xpdf1->GetEntry(entry); */
      /* b_Generator_xpdf2->GetEntry(entry); */
      /* b_Generator_id1->GetEntry(entry); */
      /* b_Generator_id2->GetEntry(entry); */
       b_nGenVisTau->GetEntry(entry); 
       b_GenVisTau_eta->GetEntry(entry); 
       b_GenVisTau_mass->GetEntry(entry); 
       b_GenVisTau_phi->GetEntry(entry); 
       b_GenVisTau_pt->GetEntry(entry); 
       b_GenVisTau_charge->GetEntry(entry); 
       b_GenVisTau_genPartIdxMother->GetEntry(entry); 
       b_GenVisTau_status->GetEntry(entry); 
      /* b_genWeight->GetEntry(entry); */
      /* b_LHEWeight_originalXWGTUP->GetEntry(entry); */
      /* b_nLHEPdfWeight->GetEntry(entry); */
      /* b_LHEPdfWeight->GetEntry(entry); */
      /* b_nLHEScaleWeight->GetEntry(entry); */
      /* b_LHEScaleWeight->GetEntry(entry); */
      /* b_nPSWeight->GetEntry(entry); */
      /* b_PSWeight->GetEntry(entry); */
      /* b_GenMET_phi->GetEntry(entry); */
      /* b_GenMET_pt->GetEntry(entry); */
      /* b_genTtbarId->GetEntry(entry); */
      b_Tau_genPartIdx->GetEntry(entry);
      b_Tau_genPartFlav->GetEntry(entry);
      b_Electron_genPartIdx->GetEntry(entry);
      b_Electron_genPartFlav->GetEntry(entry);
      /* b_GenJetAK8_partonFlavour->GetEntry(entry); */
      /* b_GenJetAK8_hadronFlavour->GetEntry(entry); */
      b_GenJet_partonFlavour->GetEntry(entry);
      b_GenJet_hadronFlavour->GetEntry(entry);
      b_Jet_genJetIdx->GetEntry(entry);
      b_Jet_hadronFlavour->GetEntry(entry);
      b_Jet_partonFlavour->GetEntry(entry);
      b_Muon_genPartIdx->GetEntry(entry);
      b_Muon_genPartFlav->GetEntry(entry);
      /* b_Photon_genPartIdx->GetEntry(entry); */
      /* b_Photon_genPartFlav->GetEntry(entry); */
      /* b_MET_fiducialGenPhi->GetEntry(entry); */
      /* b_MET_fiducialGenPt->GetEntry(entry); */
      b_Pileup_nTrueInt->GetEntry(entry);
      b_Pileup_nPU->GetEntry(entry);
      b_Pileup_sumEOOT->GetEntry(entry);
      b_Pileup_sumLOOT->GetEntry(entry);
    }

    if(_year==2017){
      b_METFixEE2017_phi->GetEntry(entry);
      b_METFixEE2017_pt->GetEntry(entry);
      if(_data==1 && _lep==0 && _era!="B" && _era!="C")
	b_HLT_Ele32_WPTight_Gsf->GetEntry(entry);
      if(_data==1)
	b_HLT_Ele32_WPTight_Gsf_L1DoubleEG->GetEntry(entry);
    }

    if(_year==2016)
      b_HLT_IsoTkMu24->GetEntry(entry);

    if(_year==2018){
      b_HLT_Ele32_WPTight_Gsf->GetEntry(entry);
      b_HLT_Ele32_WPTight_Gsf_L1DoubleEG->GetEntry(entry);
    }
    b_Flag_ecalBadCalibFilter->GetEntry(entry);   //!

    return 1;
  }

  return 0;
    
  
}


#endif // #ifdef RHAna_cxx
