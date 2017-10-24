#define ATLAS_ss_leptons_cxx
// The class definition in ATLAS_multib.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("ATLAS_ss_leptons.C")
// root> T->Process("ATLAS_ss_leptons.C","some options")
// root> T->Process("ATLAS_ss_leptons.C+")
//

#include "ATLAS_ss_leptons.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>
#include <iostream>
using namespace std;
#include <fstream>
#include <vector>
#include <string>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TTree.h>
#include <TFile.h>


Int_t neventraw;
Double_t nevent,levent,levent1,levent2,levent3,levent4,levent5;
Double_t nSR3L1, nSR3L2, nSR0b1, nSR0b2, nSR1b, nSR3b, nSR1bDD, nSR3bDD, nSR1bGG, nchargeE, nchargeMu;
Bool_t use_kfac = 0;
Double_t kfactor, kgogo, kgosq, ksqsq, ksqsb;

Double_t HT, meff_incl, meff4j, mT, mT_min_b, deltaPhi4j_min;
//Bool_t human_output = 0;

//Bool_t isICHEP = 1;  //if 0 does the analyses on the 2015 dataset, if 1 uses the 2016 ICHEP release

void ATLAS_ss_leptons::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   //initialize the running total of the total number of events in the root file
   neventraw=0;
   nevent = 0;
   levent = 0;
   levent1 = 0;
   levent2 = 0;
   levent3 = 0;   
   levent4 = 0;   
   levent5 = 0;   
   //initialize the running total of number of events in each signal region
   nSR3L1=0;
   nSR3L2=0;
   nSR0b1=0;
   nSR0b2=0;
   nSR1b=0; 
   nSR3b=0;
   nSR1bDD=0;
   nSR3bDD=0;
   nSR1bGG=0;
   nchargeE=0;
   nchargeMu=0;
  
   // reading k-factors file (which just contains a triple kgogo, kgosq, ksqsq, ksbsq)
   if (use_kfac){
     ifstream kfile("kfactors.dat");
     kfile >> kgogo >> kgosq >> ksqsq >> ksqsb;
     kfile.close(); 
   }

}

void ATLAS_ss_leptons::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t ATLAS_ss_leptons::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ATLAS_ss_leptons::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  //get the next event in the root file
  fChain->GetTree()->GetEntry(entry);
  neventraw++;
  
  // reading K-factors from file - done once per root file, in Begin section
  kfactor=1; // initialize this in case anything goes wrong
  //identifying sparticles in hard process
  int nSUSY=0; // counts sparticles in event
  int ngo = 0;  // counts gluinos
  int nsq = 0;  // counts  squarks 
  int nsb = 0;  // counts anti squarks
  
  if ( use_kfac ) {
  for (int i = 0; i< Particle_size; i++){ // HERE: I could put 20 instead of Particle_size, it works anyway
    //if (nevent/kfactor <= 1 && i < 100){	    cout << " {"<< i <<", "<<Particle_Status[i]<<","<< Particle_PID[i] << "} "<<nSUSY<< ngo<< nsq<<endl;}
//    if (TMath::Abs(Particle_Status[i]) < 30 && TMath::Abs(Particle_Status[i]) > 20 ){ // hardest process has status 21-29 (Abs not needed)
      if (TMath::Abs(Particle_PID[i]) > 1000000){   // SUSY particles (antiparticels have minus sign)
        nSUSY++;
        if (nSUSY <= 2) {  // stops after the first two SUSY particles (core process)
          if (Particle_PID[i] == 1000021) ngo++;
          if ((Particle_PID[i] >= 1000001 && Particle_PID[i] <= 1000006) || (Particle_PID[i] >= 2000001 && Particle_PID[i] <= 2000006)) nsq++; // counts squarks
          if ((-Particle_PID[i] >= 1000001 && -Particle_PID[i] <= 1000006) || (-Particle_PID[i] >= 2000001 && -Particle_PID[i] <= 2000006)) nsb++; // counts anti-squarks
          //if (nevent/kfactor <= 1){
	  //cout << " {"<< i <<", "<<Particle_Status[i]<<","<< Particle_PID[i] << "} "<< nSUSY << ngo << nsq <<endl;;}
        }
      }
//    }
  }
  // assign kfactor to event
  if (ngo == 2) kfactor=kgogo;
  if ((ngo == 1 && nsq == 1) || (ngo == 1 && nsb ==1)  ) kfactor=kgosq; // gluino-(anti)squark production
  if (nsq == 1 && nsb == 1) kfactor=ksqsb;
  if ( nsq== 2  || nsb == 2 ) kfactor=ksqsq;  
  }
  //increment the number of events
  nevent+=kfactor;
  
  
  //dummy 4-vector useful shortly
  TLorentzVector temp;
  TLorentzVector temp1;
  TLorentzVector temp2;
  //create vectors of 4-vectors which will be populated by the particles in the event
  vector<TLorentzVector> jet;
  vector<TLorentzVector> bjet;
  vector<TLorentzVector> electron;
  vector<TLorentzVector> sigElectron;
  vector<TLorentzVector> sigElectron20;
  vector<TLorentzVector> ElectronP;
  vector<TLorentzVector> ElectronM;
  vector<TLorentzVector> muon;
  vector<TLorentzVector> muon20;
  vector<TLorentzVector> muonP;
  vector<TLorentzVector> muonM;
  TLorentzVector lept;
 
//  vector<TLorentzVector> tau;
//  vector<TLorentzVector> photon;

  //loop over the electrons in the event, fill the electron vector with a list of 4-vectors, ordered by p_T
  //electron candidates have p_T>10, |eta|<2.47
  for (int i = 0; i < LooseElectron_size; i++){
    temp.SetPtEtaPhiM(LooseElectron_PT[i],LooseElectron_Eta[i],LooseElectron_Phi[i],0.0);
    if(temp.Pt() > 10 && TMath::Abs(temp.Eta()) < 2.0){
      sigElectron.push_back(temp);
      if(temp.Pt() > 20) {
        sigElectron20.push_back(temp);;
        if(LooseElectron_Charge[i] == -1) {
          nchargeE++;
          ElectronM.push_back(temp);
        }
        if(LooseElectron_Charge[i] == +1) {
          ElectronP.push_back(temp);
        }
      }
      //cout << nchargeE << "    nchargeE" << endl;
      //cout << Electron_Charge[i] << "    chargeE" << endl;
      
    }
  }


 //loop over the electron candidates to get the signal electrons in the event, fill the signal electron vector with a list of 4-vectors, ordered by p_T
  //Signal electrons  have p_T>10, |eta|<2
  //for (int i = 0; i < (int)electron.size(); i++){
  //temp.SetPtEtaPhiM(electron[i].Pt(),electron[i].Eta(),electron[i].Phi,0.0);
  //if(TMath::Abs(temp.Eta()) < 2){
  //  sigElectron.push_back(temp);
  //}
  //}

  //loop over the muons in the event, fill the muon vector with a list of 4-vectors, ordered by p_T
  //muon candidates have p_T>10, |eta|<2.5
  for (int i = 0; i < LooseMuon_size; i++){
    temp.SetPtEtaPhiM(LooseMuon_PT[i],LooseMuon_Eta[i],LooseMuon_Phi[i],0.0);
    if(temp.Pt() > 10 && TMath::Abs(temp.Eta()) < 2.5){
      muon.push_back(temp);
      if(LooseMuon_Charge[i] == -1) nchargeMu++;
      //cout << nchargeMu << "     nchargeMu" << endl;
      if(temp.Pt() > 20) {
        muon20.push_back(temp);
        if(LooseMuon_Charge[i] == -1) {
          nchargeMu++;
          muonM.push_back(temp);
        }
        if(LooseMuon_Charge[i] == +1) {
          muonP.push_back(temp);
        }
      }
    }
  }
  //missing E_T 4-vector which will be populated shortly
  TLorentzVector missing_ET;
  //fill the missing ET 4-vector - add muons to it
  missing_ET.SetPtEtaPhiM(MissingET_MET[0],0.0,MissingET_Phi[0],0.0);
  double mux, muy; mux=0; muy=0;
  TLorentzVector sum = TLorentzVector(0,0,0,0);
  TLorentzVector sum1 = TLorentzVector(0,0,0,0);
  for(unsigned int i = 0; i < muon.size(); i++){
    mux += muon[i].Px();
    muy += muon[i].Py();
  }
  double newMETx = missing_ET.Px()-mux;
  double newMETy = missing_ET.Py()-muy;
  missing_ET.SetPxPyPzE(newMETx,newMETx,0.0,sqrt(newMETx*newMETx + newMETy*newMETy));

  //loop over the jets in the event, fill the jet vector with a list of 4-vectors, ordered by p_T
  //jets must have p_T>20, |eta|<2.8
  for (int i = 0; i < Jet_size; i++){
    temp.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
    if(temp.Pt() > 20 && TMath::Abs(temp.Eta()) < 2.8 /*&& Jet_TauTag[i] == 0*/){
	  jet.push_back(temp);
	  
	  if (Jet_BTag[i] && TMath::Abs(temp.Eta()) < 2.5){
	    bjet.push_back(temp);	  
	  }
	}
  }
  /*if (bjet.size()>=3){ levent++;
  //cout << "Event with 3+ b-jets: "<< neventraw <<endl;
  for (int i = 0; i< Particle_size; i++){ // HERE: I could put 20 instead of Particle_size, it works anyway
    //if (nevent/kfactor <= 1 && i < 100){	    cout << " {"<< i <<", "<<Particle_Status[i]<<","<< Particle_PID[i] << "} "<<nSUSY<< ngo<< nsq<<endl;}
    if (neventraw == 8907 ){//TMath::Abs(Particle_Status[i]) < 30 && TMath::Abs(Particle_Status[i]) > 20 ){ 
       cout << i <<", "<<Particle_Status[i]<<","<< Particle_PID[i] << endl;
    }
  }  //
  } */ 

  //calculate m_eff(incl) for Gtt: sum of pT of jets, letpons and ET
  meff_incl = 0;
  for (size_t i = 0; i < jet.size(); i++){
    meff_incl += jet[i].Pt();
  }
  for (size_t i = 0; i < sigElectron.size(); i++){
    meff_incl += sigElectron[i].Pt();
  }
  for (size_t i = 0; i < muon.size(); i++){
    meff_incl += muon[i].Pt();
  }
  //add MET
  meff_incl += missing_ET.Pt();
  
  if (neventraw == 28){
    cout << muonP.size() << ElectronP.size()  << muonM.size() << ElectronM.size(); 
  }
  
  if ( (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2 )){
    levent1++;
    if (bjet.size() >= 3  && bjet[0].Pt() > 20 && bjet[1].Pt() > 20 && bjet[2].Pt() > 20) {
       levent2++;
       if (jet.size() >=6 && jet[5].Pt() > 25 ){
         levent3++;
         if( missing_ET.Pt() > 150){
           levent4++;
           if (meff_incl > 600 )levent5++;
         }
       }
    }
  }
  
  //SR3L1 : N_Sig_lep>=3, N_bjets_20=0,Njets>=4, pT_jets>40, ET_miss>150, meff>=0    
  if ( (sigElectron.size()+muon.size()) >= 3 && bjet.size() == 0 && jet.size() >= 4 && jet[3].Pt() > 40 && missing_ET.Pt() > 150 ) nSR3L1+=kfactor;

 //SR3L2 : N_Sig_lep>=3, N_bjets_20=0,Njets>=4, pT_jets>40, ET_miss>200, meff>=1500    
  if ( (sigElectron.size() + muon.size()) >= 3 && bjet.size() == 0  && jet.size() >= 4 && jet[3].Pt() > 40 && missing_ET.Pt() > 200 && meff_incl > 1500 ) nSR3L2+=kfactor;

 //SR0b1 : N_Sig_lep>=2, N_bjets_20=0, Njets>=6, pT_jets>25, ET_miss>150, meff>=500    
  if ( (sigElectron20.size()+muon20.size()) >= 2 && (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2) && bjet.size() == 0  && jet.size() >= 6 && jet[5].Pt() > 25 && missing_ET.Pt() > 150 && meff_incl > 500 ) nSR0b1+=kfactor;

 //SR0b2 : N_Sig_lep>=2, N_bjets_20=0, Njets>=6, pT_jets>40, ET_miss>150, meff>=900    
  if ( (sigElectron20.size()+muon20.size()) >= 2 && (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2)&& bjet.size() == 0  && jet.size() >= 6 && jet[5].Pt() > 40 && missing_ET.Pt() > 150 && meff_incl > 900 ) nSR0b2+=kfactor;

 //SR1b : N_Sig_lep>=2, N_bjets_20=1, Njets>=6, pT_jets>25, ET_miss>200, meff>=650    
  if ( (sigElectron20.size()+muon20.size()) >= 2 && (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2) && bjet.size() >= 1 && bjet[0].Pt() > 20  && jet.size() >=6 && jet[5].Pt() > 25 && missing_ET.Pt() > 200 && meff_incl > 650 ) nSR1b+=kfactor;

 //SR3b : N_Sig_lep>=2, N_bjets_20=3, Njets>=6, pT_jets>25, ET_miss>150, meff>=600    
  if ( (sigElectron20.size()+muon20.size()) >= 2 && (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2) && bjet.size() >= 3  && bjet[0].Pt() > 20 && bjet[1].Pt() > 20 && bjet[2].Pt() > 20 && jet.size() >=6 && jet[5].Pt() > 25 && missing_ET.Pt() > 150 && meff_incl > 600 ) nSR3b+=kfactor;

 //SR1bDD : N_Sig_lep>=2, N_bjets_20=1, Njets>=4, pT_jets>50, ET_miss>= 0, meff>=1200 && >= 2 negatively charged leptons    
  if ( (sigElectron20.size()+muon20.size()) >= 2 && bjet.size() >=1  && bjet[0].Pt() > 20 && jet.size() >= 4 && jet[3].Pt() > 50 &&  meff_incl > 1200 && muonM.size()+ElectronM.size() >=2 )  nSR1bDD+=kfactor;


 //SR3bDD : N_Sig_lep>=2, N_bjets_20=3, Njets>=4, pT_jets>50, ET_miss>=0, meff>=100  && >= 2 negatively charged leptons     
  if ( (sigElectron20.size()+muon20.size()) >= 2 && bjet.size() >= 3   && bjet[0].Pt() > 20 && bjet[1].Pt() > 20 && bjet[2].Pt() > 20 && jet.size() >= 4 && jet[3].Pt() > 50  && meff_incl > 1000 && muonM.size()+ElectronM.size() >=2) nSR3bDD+=kfactor;

 //SR1bGG : N_Sig_lep>=2, N_bjets_20=1 , Njets>=6, pT_jets>50, ET_miss>= 0, meff>=1800   
  if ( (sigElectron20.size()+muon20.size()) >= 2 && (muonP.size()+ElectronP.size() >=2 || muonM.size()+ElectronM.size() >=2) && bjet.size() >= 1   && bjet[0].Pt() > 20  && jet.size() >=6 && jet[5].Pt() > 50  && meff_incl > 1800 ) nSR1bGG+=kfactor;
  
  return kTRUE;
}

void ATLAS_ss_leptons::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void ATLAS_ss_leptons::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.   
   
  //if(human_output){
  cout << neventraw << endl;
  cout << nevent << endl;
//  cout << levent << endl;
  cout << levent1/neventraw*607 << endl;
  cout << levent2/neventraw*607 << endl;
  cout << levent3/neventraw*607 << endl;
  cout << levent4/neventraw*607 << endl;
  cout << levent5/neventraw*607 << endl;
  
  //cout << "\t" << "Gbb-SR A" << "\t" << "Gbb-SR B" << "\t" << "Gbb-SR C" << endl;
  cout << "\t"  << 1.0 * nSR3L1/nevent;
  cout << "\t"  << 1.0 * nSR3L2/nevent;
  cout << "\t"  << 1.0 * nSR0b1/nevent;
  cout << "\t"  << 1.0 * nSR0b2/nevent;
  cout << "\t"  << 1.0 * nSR1b/nevent;
  cout << "\t"  << 1.0 * nSR3b/nevent;
  cout << "\t"  << 1.0 * nSR1bDD/nevent;
  cout << "\t"  << 1.0 * nSR3bDD/nevent;
  cout << "\t"  << 1.0 * nSR1bGG/nevent;
  cout<< endl;
  
  /*
  }
  else {
  cout << 1.0 * nbbA/nevent;
  cout << "\t"  << 1.0 * nbbB/nevent;
  if (!isICHEP) cout << "\t"  << 1.0 * nbbC/nevent;
  cout<< endl;
  */
  
  ofstream outfile;
  outfile.open("efficiencies_combined.dat", ios::app); 
  outfile << "ATLAS_ss_leptons"<<"\t";
  outfile << 1.0 * nSR3L1/nevent;
  outfile << "\t"  << 1.0 * nSR3L2/nevent;
  outfile << "\t"  << 1.0 * nSR0b1/nevent;
  outfile << "\t"  << 1.0 * nSR0b2/nevent;
  outfile << "\t"  << 1.0 * nSR1b/nevent;
  outfile << "\t"  << 1.0 * nSR3b/nevent;
  outfile << "\t"  << 1.0 * nSR1bDD/nevent;
  outfile << "\t"  << 1.0 * nSR3bDD/nevent;
  outfile << "\t"  << 1.0 * nSR1bGG/nevent;
  outfile << endl;
  outfile.close();
  
  
  //  }


}
