#define RHAna_cxx

#include "RHAna.h"
#include <TH2.h>
#include <TStyle.h>
#include<fstream>
#include<iomanip>

// The begin functions are called before the process function
// Things that should happen before looping over all events
// go here.
void RHAna::Begin(TTree * /*tree*/)
{
   TString option = GetOption();
}
void RHAna::SlaveBegin(TTree * /*tree*/)   
{
   TString option = GetOption();
   nEvtTotal = 0;
   nEvtRan = 0;
   //Create the histogram file
   _HstFile = new TFile(_HstFileName,"recreate");
   //Call the function to book the histograms we declared in Hists.
   BookHistograms();

}

void RHAna::SlaveTerminate()
{
  //Write histograms and close histogram file
  _HstFile->Write();
  _HstFile->Close();
  //Output to screen.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total events = "<<nEvtTotal<<endl;
  //Open the text output file
  ofstream fout(_SumFileName);
  //Put text output in the summary file.
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total events  = "<<nEvtTotal<<endl;

}
void RHAna::Terminate()
{
  //   cout<<"Inside Terminate()"<<endl;
}



Bool_t RHAna::Process(Long64_t entry)
{

  // Decide to read in whole event or just some branches.
  // ReadLimited(0,..) reads in whole event
  // ReadLimited(1,..) reads in the branches defined in that function in header file.
  int readevent = ReadLimited(0,entry);
  if(readevent==0){ cout<<"Did not read in any branches.. quitting."<<endl; return kTRUE;}
  
  //Output processing information to screen based on verbosity level.
  if(_verbosity==0 && nEvtTotal%1000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;

  nEvtTotal++;

  // Some cleaning flags for running over different years.
  GoodEvt2018 = (_year==2018 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && (_data ? Flag_eeBadScFilter : 1) : 1);

  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;

  // If this event passes above filters, then we process it futher.
  if(GoodEvt){
  

    nEvtRan++;

    /*###############################################################################################
    ##################################      SELECTION PROCESS     ###################################
    ###############################################################################################*/

    
    llep.clear(); 

    
//############  MUONS ############
    
    int nmu=0;
    int nmu_pos=0;
    int nmu_neg=0;
    goodMu.clear();
    goodMu_pos.clear();
    goodMu_neg.clear();
    
    for(unsigned int i=0; i< (nMuon); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105);
      temp.id = -13*Muon_charge[i]; temp.ind = i;  temp.charge = Muon_charge[i];
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;

      h.isoM[0]->Fill(Muon_pfRelIso04_all[i]);

      if(passCuts){
	goodMu.push_back(temp);
	llep.push_back(temp);}
      
      if(passCuts && temp.charge == 1){ //Creating the mu+ array
	goodMu_pos.push_back(temp);}
      
      if(passCuts && temp.charge == -1){ //creating the mu- array
	goodMu_neg.push_back(temp);}}

    Sort(6); //for mu+
    Sort(7); //for mu-
    
    
//#########  ELECTRONS  ###############
    
   int nel=0;
    goodEle.clear();
    for(unsigned int i=0; i< (nElectron); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
      temp.id = -11*Electron_charge[i]; temp.ind = i; temp.charge = Electron_charge[i];
      
      bool isprompt = false;
      if(fabs(temp.v.Eta())<=1.479)
	if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	  isprompt = true;      
      if(fabs(temp.v.Eta())>1.479)
	if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	  isprompt = true;

      bool passCuts = temp.v.Pt()>2 && fabs(temp.v.Eta())<2.5 && Electron_cutBased[i]>2;
      passCuts = passCuts && isprompt;

      //dR cleaning for the i-th electron :
      bool ismuonclose = false;
      for(int j=0; j<(int)goodMu.size(); j++){
	float dR_mue =goodMu.at(j).v.DeltaR(temp.v);
	if(dR_mue<0.4)
	  ismuonclose = true;
      }

      h.isoE[0]->Fill(Electron_pfRelIso03_all[i]);
      
      if(passCuts && !ismuonclose){
	goodEle.push_back(temp);
	llep.push_back(temp);
	nel++;
      }
    }

    Sort(0); //for muons and light leptons
    Sort(2); // for electrons


//*************ISOLATED Tracks********************
    
    isotracks.clear();
    for(unsigned int i=0; i< (nIsoTrack); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(IsoTrack_pt[i],IsoTrack_eta[i],IsoTrack_phi[i],0);
      // bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.5 && fabs(IsoTrack_pdgId[i])==15;
      // passCuts= passCuts && fabs(IsoTrack_dxy[i])<0.2 && fabs(IsoTrack_dz[i])<0.1;
      // if(passCuts){
      //isotracks.push_back(temp);}}
      isotracks.push_back(temp);}
    h.isotrack[0]->Fill((int)isotracks.size());
    for(unsigned int i=0; i<isotracks.size();i++){
      h.isotrack[1]->Fill(isotracks.at(i).v.Pt());
      h.isotrack[2]->Fill(isotracks.at(i).v.Eta());
      h.isotrack[3]->Fill(isotracks.at(i).v.Phi());}
      
    
 //**************GenVisTaus*********************
    
    gentaus.clear();
    for(unsigned int i=0; i<(UInt_t)nGenVisTau; i++){
      Lepton temp; temp.v.SetPtEtaPhiM(GenVisTau_pt[i],GenVisTau_eta[i],GenVisTau_phi[i],1.77); 
      temp.id = -15*GenVisTau_charge[i]; temp.ind = i; temp.charge = GenVisTau_charge[i];
      bool passCuts = temp.v.Pt()>5 && fabs(temp.v.Eta())<2.5;
      if(passCuts){
	gentaus.push_back(temp);}}
    h.Gentau[0]->Fill((int)gentaus.size());
    for(unsigned int j=0; j<gentaus.size();j++){
      h.Gentau[1]->Fill(gentaus.at(j).v.Pt());
      h.Gentau[2]->Fill(gentaus.at(j).v.Eta());
      h.Gentau[3]->Fill(gentaus.at(j).v.Phi());}
    
    
 //###########  TAUS  ################
    
   int ntau =0;
    taus.clear();
    for(unsigned int i=0; i< (nTau); i++){
      if(Tau_idDecayMode[i]&&(Tau_decayMode[i]<3||Tau_decayMode[i]>9)){
	//Tau energy scale correction
	float tlv_corr = 1.;
	if(_year==2016){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.994;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]>9)  tlv_corr = 1;
	}     
	if(_year==2017){
	  if(Tau_decayMode[i]==0) tlv_corr = 1.007;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==10) tlv_corr = 1.001;
	  if(Tau_decayMode[i]==11) tlv_corr = 0.999;
	}
	if(_year==2018){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.987;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]==10) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==11) tlv_corr = 1;
	}
	Lepton temp; temp.v.SetPtEtaPhiM(Tau_pt[i],Tau_eta[i],Tau_phi[i],1.77);
	temp.v *= tlv_corr; //energy correction
	temp.id = -15*Tau_charge[i]; temp.ind = i; temp.charge = Tau_charge[i];
	temp.lepcleaning = TaulepCleaning(temp.v);

	bool passCuts = temp.v.Pt()>20 && fabs(temp.v.Eta()<2.3);
	passCuts = passCuts && Tau_idDecayModeNewDMs[i]==1 && (Tau_decayMode[i]<3 || Tau_decayMode[i]>9);
	passCuts = passCuts && fabs(Tau_dz[i])<0.2;
	passCuts = passCuts && Tau_idDeepTau2017v2p1VSe[i] >= 15 && Tau_idDeepTau2017v2p1VSmu[i] >= 3; //loose WP
	passCuts = passCuts && Tau_idDeepTau2017v2p1VSjet[i] >= 31; //medium WP

	//dR cleaning for tau :
	
	bool ismuonclose = false;
	for(int j=0; j<(int)goodMu.size(); j++){
	  float dR_mutau =goodMu.at(j).v.DeltaR(temp.v);
	  if(dR_mutau<0.4)
	    ismuonclose = true;
	}

	bool iselectronclose = false;
	for(int j=0; j<(int)goodEle.size(); j++){
	  float dR_etau =goodEle.at(j).v.DeltaR(temp.v);
	  if(dR_etau<0.4)
	    iselectronclose = true;
	    }
	
	if(passCuts && !ismuonclose && !iselectronclose){
	  taus.push_back(temp);
	  ntau++;
	}
      }
    }
    
    Sort(1);//for taus

    
    
//#############  JETS ################
  
    goodJet.clear();
    allJet.clear();
    bJet.clear();
    
    for(unsigned int i=0; i<(nJet); i++ ){
      Jet temp; temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]); 
      temp.ind = i; 
      bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<3.0;
      
      //dR cleaning for the i-th jet :
      
      bool ismuonclose = false;
      for(int j=0; j<(int)goodMu.size(); j++){
	float dR_mujet =goodMu.at(j).v.DeltaR(temp.v);
	if(dR_mujet<0.4){
	  ismuonclose = true;
	}
      }
      
      bool iselectronclose = false;
      for(int j=0; j<(int)goodEle.size(); j++){
	float dR_ejet =goodEle.at(j).v.DeltaR(temp.v);
	if(dR_ejet<0.4)
	  iselectronclose = true;
      }
      
       bool istauclose = false;
      for(int j=0; j<(int)taus.size(); j++){
	float dR_taujet =taus.at(j).v.DeltaR(temp.v);
	if(dR_taujet<0.4)
	  istauclose = true;
	  }
      
      if(passCuts){
	allJet.push_back(temp);
      }

      for(int j=0; j<(int)goodMu.size(); j++){
	if(passCuts && ismuonclose){
	  h.two[2]->Fill(temp.v.Pt(), goodMu.at(j).v.Pt());
	}
      }
      
      if(passCuts && !ismuonclose && !iselectronclose && !istauclose){
	goodJet.push_back(temp);
      }

      h.btag[0]->Fill(fabs(Jet_btagDeepB[i]));
      h.btag[1]->Fill(fabs(Jet_btagCSVV2[i]));

      bool passCuts_b = passCuts && !ismuonclose && !iselectronclose && !istauclose && Jet_btagDeepB[i]>0.8 && 0.8<Jet_btagCSVV2[i]<1.1;

      if(passCuts_b){
	bJet.push_back(temp);
      }

    }
    
    Sort(3); //For goodJet
    Sort(5); //for allJet
    
    
//##############  PHOTONS ####################
    
    int nphot=0;
    goodPhoton.clear();
    
    for(unsigned int i=0; i<(nPhoton); i++ ){
      //Prepare a temporary object
      Photon temp; temp.v.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],Photon_phi[i],0.000000); 
      temp.ind = i; 
      // this boolean passCuts is where all selections are applied
      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<3.0;
      // if this object passes cuts, then put in the array

      //dR Cleaning :

      bool ismuonclose = false;
      for(int j=0; j<(int)goodMu.size(); j++){
	float dR_muphot =goodMu.at(j).v.DeltaR(temp.v);
	if(dR_muphot<0.4){
	  ismuonclose = true;
	}
      }
      
      bool iselectronclose = false;
      for(int j=0; j<(int)goodEle.size(); j++){
	float dR_ephot =goodEle.at(j).v.DeltaR(temp.v);
	if(dR_ephot<0.4)
	  iselectronclose = true;
      }
      
      bool istauclose = false;
      for(int j=0; j<(int)taus.size(); j++){
	float dR_tauphot =taus.at(j).v.DeltaR(temp.v);
	if(dR_tauphot<0.4)
	  istauclose = true;
	  }

      bool isjetclose = false;
      for(int j=0; j<(int)goodJet.size(); j++){
	float dR_jetphot =goodJet.at(j).v.DeltaR(temp.v);
	if(dR_jetphot<0.4)
	  isjetclose = true;
	  }

      passCuts = passCuts && !ismuonclose && !iselectronclose && !istauclose && !isjetclose;
      
      if(passCuts){
	goodPhoton.push_back(temp);
       	nphot++;
      }
    }
    Sort(4); //for photons

//################  MET ####################
       
    // This is where missing Et for the event is defined.
    if(_year == 2017){
      metpt = METFixEE2017_pt;
      metphi = METFixEE2017_phi;
    }
    else{
      metpt = MET_pt;
      metphi = MET_phi;
    }

/*##########################   Adding GEN particles  #########################################*/

    int genmu =0;
    GenMu.clear();
    for(unsigned int i=0; i<(nGenPart); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
      temp.ind = i;

      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4;
      bool truemuon = false;
      if(fabs(GenPart_pdgId[i])==13 && GenPart_status[i]==1)
	truemuon = true;
      
      if(truemuon==true && passCuts){
	GenMu.push_back(temp);
	genmu++;
      }	
    }
    Sort(8); //Sorting GenMu
    h.Genmu[0]->Fill((int)GenMu.size());
    
    int genele =0;
    GenEle.clear();
    for(unsigned int i=0; i<(nGenPart); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
      temp.ind = i;

      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.5;
      bool true_electron = false;
      if(fabs(GenPart_pdgId[i])==11 && GenPart_status[i]==1)
	true_electron = true;
      
      if(true_electron && passCuts){
	GenEle.push_back(temp);
	genele++;
      }	
    }
    Sort(9); //Sorting GenEle
    h.Genele[0]->Fill((int)GenEle.size());
    
    ofstream file("Truth.txt");
    file<<"Index"<<setw(9)<<"PdgId"<<setw(11)<<"MotherId"<<setw(10)<<"Status"<<endl;
    for(unsigned int i=0;  i<nGenPart; i++){
      file<<i<<setw(10)<<GenPart_pdgId[i]<<setw(10)<<GenPart_genPartIdxMother[i]<<setw(11)<<GenPart_status[i];
      file<<endl;
    }
   
    file.close();
    
    //########################################################################################
    //################################     ANALYSIS    #######################################
    //########################################################################################
    
    h.llep->Fill((int)llep.size());

    //##### [0] dimuon events ############
    h.muons[0]->Fill((int)goodMu.size());
    if((int)goodMu.size()>1){
      h.muons[2]->Fill(goodMu.at(0).v.Pt());
      h.muons[3]->Fill(goodMu.at(1).v.Pt());
      h.muons[1]->Fill( (goodMu.at(0).v+goodMu.at(1).v).M());
      h.dphi[0]->Fill(goodMu.at(0).v.DeltaPhi(goodMu.at(1).v));
    }

    
    //##### [1] dielectron events ############
    h.electrons[0]->Fill((int)goodEle.size());
    if((int)goodEle.size()>1){
      h.electrons[2]->Fill(goodEle.at(0).v.Pt());
      h.electrons[3]->Fill(goodEle.at(1).v.Pt());
      h.electrons[1]->Fill((goodEle.at(0).v+goodEle.at(1).v).M());
      h.dphi[1]->Fill(goodEle.at(0).v.DeltaPhi(goodEle.at(1).v));
    }

    //##### [2] emu events ##############
    if((int)goodMu.size()>0 && (int)goodEle.size()>0){
      h.emu[0]->Fill((int)goodEle.size() + (int)goodMu.size());
      h.emu[1]->Fill((goodMu.at(0).v+goodEle.at(0).v).M() );
    }
    
    
    //#### [3] tau-tau events  #################
    h.taus[0]->Fill((int)taus.size());
    if((int)taus.size()>1){
      h.taus[2]->Fill(taus.at(0).v.Pt());
      h.taus[3]->Fill(taus.at(1).v.Pt());
      h.taus[1]->Fill((taus.at(0).v+taus.at(1).v).M());
    }
    
    //##### [4] jet-jet events  ##################
    h.jets[0]->Fill((int)goodJet.size());
    if((int)goodJet.size()>1){
      h.jets[2]->Fill(goodJet.at(0).v.Pt());
      h.jets[3]->Fill(goodJet.at(1).v.Pt());
      h.jets[1]->Fill( (goodJet.at(0).v+goodJet.at(1).v).M());
    }
    
    //##### [5] photon -photon events ##########
    h.photons[0]->Fill((int)goodPhoton.size());
    if((int)goodPhoton.size()>1){
      h.photons[2]->Fill(goodPhoton.at(0).v.Pt());
      h.photons[3]->Fill(goodPhoton.at(1).v.Pt());
      h.photons[1]->Fill((goodPhoton.at(0).v+goodPhoton.at(1).v).M());
    }

    
    // ####### 2D PLOTS ########################

    if((int)goodMu.size()>1 && (int)goodEle.size()>1){ 
      h.two[0]->Fill(goodMu.at(0).v.Pt(),goodEle.at(0).v.Pt());
      h.two[5]->Fill(((goodMu.at(0).v+goodMu.at(1).v).M()), ((goodEle.at(0).v+goodEle.at(1).v).M()));
    }

    //2D plot :
    for(int i=0; i<(int)goodMu.size(); i++){
      for(int j=0; j<(int)goodJet.size(); j++){
	h.two[1]->Fill(goodMu.at(i).v.Pt(), goodJet.at(j).v.Pt());
      }
    }
    
    //Positive and Negative Muons :
    for(int i=0; i<(int)goodMu_pos.size(); i++){
      for(int j=0; j<(int)goodMu_neg.size(); j++){
	h.two[3]->Fill(goodMu_pos.at(i).v.Pt(), goodMu_neg.at(j).v.Pt());
      }
    }
    
    

    //##########################  HT ###############################

    float ht=0.0;
    for(int i=0; i<(int)allJet.size(); i++){
      ht += allJet.at(i).v.Pt();
    }
    h.ht->Fill(ht);
    h.two[4]->Fill(ht, MET_pt);
    
    //################################  MET ANALYSIS : ########################################

    h.met[0]->Fill(MET_pt);
    h.met[1]->Fill(MET_phi);
    //plotting MET for mumu events :
    if((int)goodMu.size()>1){
      h.met[2]->Fill(MET_pt);
      h.met[5]->Fill(MET_phi);
    }
    //plotting MET for mumu +0 jet events :
    if((int)goodMu.size()>1 && (int)goodJet.size()==0){
      h.met[3]->Fill(MET_pt);
      h.met[6]->Fill(MET_phi);
    }
    //plotting MET for mumu + 1Jet events :
    if((int)goodMu.size()>1 && (int)goodJet.size()==1){
      h.met[4]->Fill(MET_pt);
      h.met[7]->Fill(MET_phi);
    }
    //plotting delta phi with METS
    for(int i=0; i<(int)goodMu.size(); i++){
      h.dphi[2]->Fill(delta_phi(goodMu.at(i).v.Phi(),MET_phi));
    }
    for(int i=0; i<(int)goodEle.size(); i++){
      h.dphi[3]->Fill(delta_phi(goodEle.at(i).v.Phi(),MET_phi));
    }
     
    //##############################   dR cleaning analysis  ####################################
    
    //Jets and Electrons
    float dRmin_JE =1000.0;
    for(int j=0; j<(int)goodEle.size();j++){
      for(int i=0; i<(int)goodJet.size(); i++){
	float  dR_value= goodEle.at(0).v.DeltaR(goodJet.at(i).v);
	if(dR_value<dRmin_JE){
	  dRmin_JE=dR_value;
	  h.dR[0]->Fill(dRmin_JE);
	}
      }
    }
   

    //Jets and Muons
    float dRmin_JM =1000.0;
    for(int j=0; j<(int)goodMu.size();j++){
      for(int i=0; i<(int)goodJet.size(); i++){
	float  dR_value= goodMu.at(0).v.DeltaR(goodJet.at(i).v);
	if(dR_value<dRmin_JM){
	  dRmin_JM=dR_value;
	  h.dR[1]->Fill(dRmin_JM);
	}
      }
    }

    //  GenMu and goodMu :
    float dRmin_matching_mu =1000.0;
    for(int j=0; j<(int)goodMu.size();j++){
      for(int i=0; i<(int)GenMu.size(); i++){
	float  dR_value= goodMu.at(j).v.DeltaR(GenMu.at(i).v);
	if(dR_value<dRmin_matching_mu){
	  dRmin_matching_mu=dR_value;
	  h.dR[2]->Fill(dRmin_matching_mu);}}}
    
    //GenTaus and Isolated tracks

    float dRmin_matching_track=1000.0;
    for(unsigned int i=0; i<isotracks.size(); i++){
      for(unsigned j=0; j<gentaus.size(); j++){
	float dR_value= isotracks.at(i).v.DeltaR(gentaus.at(j).v);
	if(dR_value<dRmin_matching_track){
	  dRmin_matching_track=dR_value;
	  h.dR[3]->Fill(dRmin_matching_track);}}}

    //##############################   GenPart analysis    ####################################
    
     //GenMu :

    if((int)GenMu.size()>1){  //dimuon events
      h.Genmu[1]->Fill(GenMu.at(0).v.Pt());
      h.Genmu[2]->Fill( (GenMu.at(0).v+GenMu.at(1).v).M() );
    }

    //GenEle :

    if((int)GenEle.size()>1){  //dielectron events
      h.Genele[1]->Fill(GenEle.at(0).v.Pt());
      h.Genele[2]->Fill( (GenEle.at(0).v+GenEle.at(1).v).M() );
    }

    //Gen e-mu events :
    if((int)GenMu.size()>0 && (int)GenEle.size()>0){
      h.Genemu[0]->Fill( (GenEle.at(0).v+GenMu.at(0).v).M() );
    }
    
    //##############################   b tagging  #######################################

    h.btag[2]->Fill((int)bJet.size());

    //############# Quad light lepton mass ###############
    
    if((int)goodMu.size()>3){
      h.lep4[0]->Fill( (goodMu.at(0).v+goodMu.at(1).v+goodMu.at(2).v+goodMu.at(3).v).M() );
    }
    
    if((int)goodEle.size()>3){
      h.lep4[1]->Fill( (goodEle.at(0).v+goodEle.at(1).v+goodEle.at(2).v+goodEle.at(3).v).M() );
    }

      
    //#############################  transverse mass calculation ####################################

    //for muons :

    if((int)goodMu.size()>0){
      float mT_muon;
      float dphi_muon_met = delta_phi(goodMu.at(0).v.Phi(), MET_phi);
      mT_muon = sqrt(2*goodMu.at(0).v.Pt()*MET_pt*(1-cos(dphi_muon_met)));
      h.mT[0]->Fill(mT_muon);
    }
    
    //################# finding matching GenPart and finding mother ##################################
    //Finding mother of the leading  muon
    
    if((int)goodMu.size()>0 && (int)GenMu.size()>0){
      //we are selecting the events which contain atleast one muon
      //first we calculate the minimum of dR with the GenMu particles
      float dRmin;
      dRmin=1000.0;
      int m=-1;
      for(int j=0; j<(int)GenMu.size(); j++){
	float dR = goodMu.at(0).v.DeltaR(GenMu.at(j).v); //for the leading muon
	if(dR<dRmin){
	  dRmin=dR;
	  m=j;
	}
      }
      //for the leading muon, dR_min has been found
      //now we find the mother of the leading muon:
      int mother;
      mother =GenMother(GenMu.at(m).ind, GenPart_genPartIdxMother[GenMu.at(m).ind]);
      h.mom[1]->Fill(mother);
      //mother of the leading muon is given by int mother.
      //This will add the mom_ID to the histogram.
      
      //dxy/dz/Isolation plotting
      
      if(fabs(mother)<30){
	h.impact[4]->Fill(Muon_dxy[goodMu.at(0).ind]);
	h.impact[5]->Fill(Muon_dz[goodMu.at(0).ind]);
	h.isoM[1]->Fill(Muon_pfRelIso03_all[goodMu.at(0).ind]);
      }
      if(fabs(mother)>100){
	h.impact[6]->Fill(Muon_dxy[goodMu.at(0).ind]);
	h.impact[7]->Fill(Muon_dz[goodMu.at(0).ind]);
	h.isoM[2]->Fill(Muon_pfRelIso03_all[goodMu.at(0).ind]);
      }
    }

    //Finding mother of the leading Electron
    
    if((int)goodEle.size()>0 && (int)GenEle.size()>0){
      //we are selecting the events which contain atleast one muon
      //first we calculate the minimum of dR with the GenMu particles
      float dRmin;
      dRmin=1000.0;
      int m=-1;
      for(int j=0; j<(int)GenEle.size(); j++){
	float dR = goodEle.at(0).v.DeltaR(GenEle.at(j).v); //for the leading muon
	if(dR<dRmin){
	  dRmin=dR;
	  m=j;
	}
      }
      //for the leading muon, dR_min has been found
      //now we find the mother of the leading muon:
      int mother;
      mother =GenMother(GenEle.at(m).ind, GenPart_genPartIdxMother[GenEle.at(m).ind]);
      h.mom[0]->Fill(mother);
      //mother of the leading muon is given by int mother.
      //This will add the mom_ID to the histogram.
      
      //dxy/dz/Isolation plotting
      
      if(fabs(mother)<30){
	h.impact[0]->Fill(Electron_dxy[goodEle.at(0).ind]);
	h.impact[1]->Fill(Electron_dz[goodEle.at(0).ind]);
	h.isoE[1]->Fill(Electron_pfRelIso03_all[goodEle.at(0).ind]);
      }
      if(fabs(mother)>100){
	h.impact[2]->Fill(Electron_dxy[goodEle.at(0).ind]);
	h.impact[3]->Fill(Electron_dz[goodEle.at(0).ind]);
	h.isoE[2]->Fill(Electron_pfRelIso03_all[goodEle.at(0).ind]);
      }
    }

  
    
  }
  return kTRUE;
}

//Functions :
// ====================
int RHAna::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
	ind = mom_ind;
	mom_ind = GenPart_genPartIdxMother[ind];
	p_id = GenPart_pdgId[ind];
	m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}
    // ====================


void RHAna::Sort(int opt)
{
  //Sort selected objects by pT (always descending).
  //Opt 0 is for sorting muons, option 2 for sorting electrons and option 1 for sorting taus.
  if(opt==0){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j)); }}
    for(int i=0; i<(int)llep.size()-1; i++){
      for(int j=i+1; j<(int)llep.size(); j++){
	if( llep[i].v.Pt() < llep[j].v.Pt() ) swap(llep.at(i),llep.at(j)); }}
  }
  if(opt==1){
    for(int i=0; i<(int)taus.size()-1; i++){
      for(int j=i+1; j<(int)taus.size(); j++){
	if( taus[i].v.Pt() < taus[j].v.Pt() ) swap(taus.at(i),taus.at(j));}}
  }
  if(opt==2){
    for(int i=0; i<(int)goodEle.size()-1; i++){
      for(int j=i+1; j<(int)goodEle.size(); j++){
	if( goodEle[i].v.Pt() < goodEle[j].v.Pt() ) swap(goodEle.at(i),goodEle.at(j));}}
  }
  if(opt==3){
    for(int i=0; i<(int)goodJet.size()-1; i++){
      for(int j=i+1; j<(int)goodJet.size(); j++){
	if( goodJet[i].v.Pt() < goodJet[j].v.Pt() ) swap(goodJet.at(i),goodJet.at(j));}}
  }
  if(opt==4){ //Photons
    for(int i=0; i<(int)goodPhoton.size()-1; i++){
      for(int j=i+1; j<(int)goodPhoton.size(); j++){
	if( goodPhoton[i].v.Pt() < goodPhoton[j].v.Pt() ) swap(goodPhoton.at(i),goodPhoton.at(j));}}
  }
  if(opt==5){ //Unfiltered jets
    for(int i=0; i<(int)allJet.size()-1; i++){
      for(int j=i+1; j<(int)allJet.size(); j++){
	if( allJet[i].v.Pt() < allJet[j].v.Pt() ) swap(allJet.at(i),allJet.at(j));}}
  }
  if(opt==6){
    for(int i=0; i<(int)goodMu_pos.size()-1; i++){
      for(int j=i+1; j<(int)goodMu_pos.size(); j++){
	if( goodMu_pos[i].v.Pt() < goodMu_pos[j].v.Pt() ) swap(goodMu_pos.at(i),goodMu_pos.at(j)); }}
  }
  if(opt==7){
    for(int i=0; i<(int)goodMu_neg.size()-1; i++){
      for(int j=i+1; j<(int)goodMu_neg.size(); j++){
	if( goodMu_neg[i].v.Pt() < goodMu_neg[j].v.Pt() ) swap(goodMu_neg.at(i),goodMu_neg.at(j)); }}
  }
  if(opt==8){
    for(int i=0; i<(int)GenMu.size()-1; i++){
      for(int j=i+1; j<(int)GenMu.size(); j++){
	if( GenMu[i].v.Pt() < GenMu[j].v.Pt() ) swap(GenMu.at(i),GenMu.at(j)); }}
  }
  if(opt==9){
    for(int i=0; i<(int)GenEle.size()-1; i++){
      for(int j=i+1; j<(int)GenEle.size(); j++){
	if( GenEle[i].v.Pt() < GenEle[j].v.Pt() ) swap(GenEle.at(i),GenEle.at(j)); }}
  }
  
  
}
float RHAna::delta_phi(float phi1, float phi2)
{
  //Calculate the correct deltaPhi=phi1-phi2
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}
bool RHAna::TaulepCleaning(TLorentzVector t)
{
  bool result=true;
  for(int i=0; i<(int)llep.size(); i++){
    if(t.DeltaR(llep[i].v)<0.4){
        result=false;
	break;
    }
  }
  return result;
} 
bool RHAna::pass_tau_cuts(int opt, int tin)
{
  bool result = false;
  result = Tau_idDecayModeNewDMs[tin]==1 && (Tau_decayMode[tin]<3 || Tau_decayMode[tin]>9);
  result = result &&  fabs(Tau_dz[tin])<0.2;
  result = result && Tau_idDeepTau2017v2p1VSe[tin] >= 15 && Tau_idDeepTau2017v2p1VSmu[tin] >= 3;
  //loose, loose
  if(opt==1) result = result && Tau_idDeepTau2017v2p1VSjet[tin] >= 7;//VLoose
  if(opt==2) result = result && Tau_idDeepTau2017v2p1VSjet[tin] >= 15;//Loose
  if(opt==3) result = result && Tau_idDeepTau2017v2p1VSjet[tin] >= 31;//Medium
  if(opt==4) result = result && Tau_idDeepTau2017v2p1VSjet[tin] >= 63;//Tight
  if(opt==5) result = result && Tau_idDeepTau2017v2p1VSjet[tin] >= 127;//VTight

  return result;

}
 
void RHAna::BookHistograms()
{
  
  //Isotrack Plots:
  h.isotrack[0]  = new TH1F("isotrack0","nIsoTrack",4,0,4);
  h.isotrack[1]  = new TH1F("isotrack1","Isotrack_Pt",100,0,100);
  h.isotrack[2]  = new TH1F("isotrack2","Isotrack_Eta",100,-4,4);
  h.isotrack[3]  = new TH1F("isotrack3","Isotrack_Phi",100,-3.3,3.3);

  //GenTau plots:
  h.Gentau[0]    = new TH1F("Gentau0","nGenTaus",4,0,4);
  h.Gentau[1]    = new TH1F("Gentau1","Gentau_Pt",200,0,150);
  h.Gentau[2]    = new TH1F("Gentau2","Gentau_Eta",100,-4,4);
  h.Gentau[3]    = new TH1F("Gentau3","Gentau_Phi",100,-3.3,3.3);
  
  //Number plots
  h.llep         = new TH1F("llep","NlightLep",10,0,10);
  h.muons[0]     = new TH1F("muons0", "Number of Muons", 15, 0, 15);
  h.electrons[0] = new TH1F("electrons0", "Number of Electrons", 15, 0, 15);
  h.emu[0]       = new TH1F("emu[0]", "No. of emu events", 15, 0, 15);
  h.taus[0]      = new TH1F("taus0", "Number of Taus", 15, 0, 15);
  h.jets[0]      = new TH1F("jets0", "Number of Jets", 15, 0, 15);
  h.photons[0]   = new TH1F("photons0", "Number of Photons", 15, 0, 15);
  
  //Mass plots
  h.muons[1]     = new TH1F("muons1", "Dimuon mass", 250, 0, 250);
  h.electrons[1] = new TH1F("electrons1", "Dielectron mass", 250, 0, 250);
  h.emu[1]       = new TH1F("emu1", "Electron_Muon_Mass", 250, 0, 250);
  h.taus[1]      = new TH1F("taus1", "DiTau mass", 250, 0, 250);
  h.jets[1]      = new TH1F("jets1", "Dijet Mass", 250, 0, 250);
  h.photons[1]   = new TH1F("photons1", "DiPhoton mass", 250, 0, 250);

  //Leading pT plots
  h.muons[2]     = new TH1F("muons2", "Leading Muon Pt", 100, 0, 100);
  h.electrons[2] = new TH1F("electrons2", "Leading Electron Pt", 100, 0, 100);
  h.taus[2]      = new TH1F("taus2", "Leading Tau Pt", 100, 0, 100);
  h.jets[2]      = new TH1F("jets2", "Leading Jet Pt", 100, 0, 100);
  h.photons[2]   = new TH1F("photons2", "Leading Photon Pt", 100, 0, 100); 
 
  h.ht = new TH1F("ht", "Sum of pT(alljets)", 1000, 0, 1000);

  //Sub-Leading pT plots

  h.muons[3]     = new TH1F("muons3", "Subleadingmuon Pt",100,0,100);
  h.electrons[3] = new TH1F("electrons3", "Subleading Electron Pt", 100, 0, 100);
  h.taus[3]      = new TH1F("taus3", "Subleading Tau Pt", 100, 0, 100);
  h.jets[3]      = new TH1F("jets3", "Subleading Jet Pt", 100, 0, 100);
  h.photons[3]   = new TH1F("photons3", "Subleading Photon Pt", 100, 0, 100); 

  //2D plots
  h.two[0]       = new TH2F("Two0", "LeadingEle_LeadingMuon_Pt", 200, 0, 200, 200, 0, 200);
  h.two[1]       = new TH2F("Two1", "Mu_Jet_Pt", 200, 0, 200, 200, 0, 200);
  h.two[2]       = new TH2F("Two2", "Mu_ Jet pT (#Delta R<0.4)", 200, 0, 200, 200, 0, 200);
  h.two[3]       = new TH2F("Two3", "mu+_mu-_ pT", 200, 0, 200, 200, 0, 200);
  h.two[4]       = new TH2F("Two4", "HT_MET_pT", 200, 0, 1000, 200, 0, 500);
  h.two[5]       = new TH2F("Two5", "Dimuon_Dielectron_mass", 200, 0, 200, 200, 0, 200);

  //Plotting METpT and METphi:
  h.met[0]       = new TH1F("met0", "Total MET pT", 200, 0, 200);
  h.met[2]       = new TH1F("met1", "met pT for dimuon events", 200, 0, 200);
  h.met[3]       = new TH1F("met2", "met pT for dimuon+0Jet events", 200, 0, 200);
  h.met[4]       = new TH1F("met3", "met pT for dimuon+1Jet events", 200, 0, 200);
  h.met[1]       = new TH1F("met4", "met phi for all events",160,-4,4);
  h.met[5]       = new TH1F("met5", "met phi for dimuon events",160,-4,4);
  h.met[6]       = new TH1F("met6", "met phi for dimuon+0 Jet events",160,-4,4);
  h.met[7]       = new TH1F("met7", "met phi for dimuon+1 Jet events",160,-4,4);

  //dR plots 
  h.dR[0]        = new TH1F("dR0","dR_Ele_Jet", 100, 0,10);
  h.dR[1]        = new TH1F("dR1","dR_Muon_Jet", 100, 0, 10);
  h.dR[2]        = new TH1F("dR2","dR_RecoMuon_GenMuon", 100, 0, 0.5);
  h.dR[3]        = new TH1F("dR3","dR_IsolatedTracks_GenTaus",100,0,2);
  
  //dphi plots :
  h.dphi[0]      = new TH1F("dphi0","dPhi_Muons",100,0,TMath::Pi());
  h.dphi[1]      = new TH1F("dphi1","dPhi_Electrons",100,0,TMath::Pi());
  h.dphi[2]      = new TH1F("dphi2", "dPhi_Muons_MET", 100, 0, TMath::Pi());
  h.dphi[3]      = new TH1F("dphi3", "dPhi_Electrons_MET", 100, 0, TMath::Pi());

  //GenMU Plots :
  h.Genmu[0]     = new TH1F("GenMu0", "NUmber of GenMuons", 10, 0, 10);
  h.Genmu[1]     = new TH1F("GenMu1", "GenMuons_leading pT", 100, 0, 100);
  h.Genmu[2]     = new TH1F("GenMu2", "Gen_Dimuon_Mass", 250, 0, 250);
  
  //GenEle Plots :
  h.Genele[0]    = new TH1F("GenEle0", "Nuber of GenElectrons", 10, 0, 10);
  h.Genele[1]    = new TH1F("GenEle1", "GenElectrons_leading pT", 100, 0, 100);
  h.Genele[2]    = new TH1F("GenEle2", "Gen_DilectronMass", 250, 0, 250);


  //GenEMu Plots :
  h.Genemu[0]    = new TH1F("Genemu0", "Gen_emuMass", 250, 0, 250);
  
  //isolation plots 
 
  h.isoE[0]       = new TH1F("isoE0", "Electron_pfRelIso03_all", 100, 0,30);
  h.isoE[1]       = new TH1F("isoE1", "Electron_pfRelIso04_all_MomId<30",100,0,0.4);
  h.isoE[2]       = new TH1F("isoE2", "Electron_pfRelIso04_all_MomId>30",100,0,0.4);
  h.isoM[0]       = new TH1F("isoM0", "Muon_pfRelIso04_all", 100, 0,30);
  h.isoM[1]       = new TH1F("isoM1", "Muon_pfRelIso04_all_MomId<30",100,0,0.4);
  h.isoM[2]       = new TH1F("isoM2", "Muon_pfRelIso04_all_MomId>30", 100, 0,0.4);
  
  //btagging
  h.btag[0]      = new TH1F("btag0","Jet_btagDeepB", 200, 0, 2.5);
  h.btag[1]      = new TH1F("btag1","Jet_btagCSVV2", 200, 0, 2.5);
  h.btag[2]      = new TH1F("btag2","Number of b Jets", 10, 0, 10);

  //Quad lepton masses
  h.lep4[0]      =  new TH1F("lep4_Muons", "QuadMuon Mass", 400, 0, 400);
  h.lep4[1]      =  new TH1F("lep4_Electrons","QuadElectron Mass", 400, 0, 400);

  //Transverse mass : Muon
  h.mT[0]        = new TH1F("mT0", "TransverseMass_Muon", 400, 0, 400);
  
  //Mother Id	
  h.mom[0]       = new TH1F("Mom0", "Electron_motherid", 1200, -600, 600);
  h.mom[1]       = new TH1F("Mom1", "Muon_motherid", 1200, -600, 600);

  //Impact Parameters
  h.impact[0]   = new TH1F("impact0","Electron_dxy_<30",100,-1,1);
  h.impact[1]   = new TH1F("impact1","Electron_dz_<30",100,-1,1);
  h.impact[2]   = new TH1F("impact2","Electron_dxy_>100",100,-1,1);
  h.impact[3]   = new TH1F("impact3","Electron_dz_>100",100,-1,1);
  h.impact[4]   = new TH1F("impact4","Muon_dxy_<30",100,-1,1);
  h.impact[5]   = new TH1F("impact5","Muon_dz_<30",100,-1,1);
  h.impact[6]   = new TH1F("impact6","Muon_dxy_>100",100,-1,1);
  h.impact[7]   = new TH1F("impact7","Muon_dz_>100",100,-1,1); 
}

