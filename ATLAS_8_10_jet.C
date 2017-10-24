#define ATLAS_8_10_jet_cxx
// The class definition in ATLAS_8_10_jet.h has been generated automatically
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
// root> T->Process("ATLAS_8_10_jet.C")
// root> T->Process("ATLAS_8_10_jet.C","some options")
// root> T->Process("ATLAS_8_10_jet.C+")
//

#include "ATLAS_8_10_jet.h"
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
Double_t nevent,levent,levent1,levent2,levent3;
Double_t n8j50MJ340,n9j50MJ340,n10j50MJ340;
Double_t n8j50MJ500,n9j50MJ500,n10j50MJ500;
Bool_t use_kfac = 0;
Double_t kfactor, kgogo, kgosq, ksqsq, ksqsb;

Double_t HT;
Double_t MJsum;
Bool_t human_output = 0;

   //checking
vector<int> nbjets;

void ATLAS_8_10_jet::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   //initialize the running total of the total number of events in the root file
   nevent = 0;
   neventraw=0;
   levent = 0;
   levent1 = 0;
   levent2 = 0;
   levent3 = 0;   
   //initialize the running total of number of events in each signal region
   n8j50MJ340 = 0;
   n9j50MJ340 = 0;
   n10j50MJ340 = 0;
   n8j50MJ500 = 0;
   n9j50MJ500 = 0;
   n10j50MJ500 = 0;
   
   // reading k-factors file (which just contains a triple kgogo, kgosq, ksqsq)
   if (use_kfac){
     ifstream kfile("kfactors.dat");
     kfile >> kgogo >> kgosq >> ksqsq >> ksqsb;
     kfile.close(); 
   }
   


}

void ATLAS_8_10_jet::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t ATLAS_8_10_jet::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ATLAS_8_10_jet::GetEntry() or TBranch::GetEntry()
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
//    if (TMath::Abs(Particle_Status[i]) < 30 && TMath::Abs(Particle_Status[i]) > 20 ){ // hardest process has status 21-29 (Abs not needed)
      if (TMath::Abs(Particle_PID[i]) > 1000000){   // SUSY particles (antiparticels have minus sign)
        nSUSY++;
        if (nSUSY <= 2) {  // stops after the first two SUSY particles (core process)
          if (Particle_PID[i] == 1000021) ngo++;
          if ((Particle_PID[i] >= 1000001 && Particle_PID[i] <= 1000006) || (Particle_PID[i] >= 2000001 && Particle_PID[i] <= 2000006)) nsq++; // counts squarks
          if ((-Particle_PID[i] >= 1000001 && -Particle_PID[i] <= 1000006) || (-Particle_PID[i] >= 2000001 && -Particle_PID[i] <= 2000006)) nsb++; // counts anti-squarks
          //cout << " {"<<i <<", "<<Particle_Status[i]<<","<< Particle_PID[i] << "} ";
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
  vector<TLorentzVector> fatjet;
  vector<TLorentzVector> electron;
  vector<TLorentzVector> muon;
//  vector<TLorentzVector> tau;
//  vector<TLorentzVector> photon;
  vector<TLorentzVector> bjet;
  vector<TLorentzVector> jet50;
  vector<TLorentzVector> jet80;
    

  //loop over the electrons in the event, fill the electron vector with a list of 4-vectors, ordered by p_T
  //electron candidates have p_T>10, |eta|<2.47
  for (int i = 0; i < LooseElectron_size; i++){
    temp.SetPtEtaPhiM(LooseElectron_PT[i],LooseElectron_Eta[i],LooseElectron_Phi[i],0.0);
    if(temp.Pt() > 10 && TMath::Abs(temp.Eta()) < 2.47){
      electron.push_back(temp);
    }
  }
  //loop over the muons in the event, fill the muon vector with a list of 4-vectors, ordered by p_T
  //muon candidates have p_T>10, |eta|<2.5
  for (int i = 0; i < LooseMuon_size; i++){
    temp.SetPtEtaPhiM(LooseMuon_PT[i],LooseMuon_Eta[i],LooseMuon_Phi[i],0.0);
    if(temp.Pt() > 10 && TMath::Abs(temp.Eta()) < 2.5){
      muon.push_back(temp);
    }
  }
  
  
  //loop over the jets in the event, fill the jet vector with a list of 4-vectors, ordered by p_T
  //jets must have p_T>20, |eta|<2.8
  for (int i = 0; i < Jet_size; i++){
    temp.SetPtEtaPhiM(Jet_PT[i],Jet_Eta[i],Jet_Phi[i],Jet_Mass[i]);
     if(temp.Pt() > 20 && TMath::Abs(temp.Eta()) < 2.8 ){    
  			jet.push_back(temp);
   			if ( (temp.Pt() > 50) && (TMath::Abs(temp.Eta())<2.0) ){
   				jet50.push_back(temp);
   			}
   			// B-tagging: 
   			if (Jet_BTag[i]){
   				bjet.push_back(temp);
   			}
    }
  }
  nbjets.push_back(bjet.size());
  
  //overlap leptons not near b-jets
  vector<int> el_rm;
  vector<TLorentzVector> clean_el;
  for (int i = 0; i < electron.size(); i++){
    for (int j = 0; j < bjet.size(); j++){
      if ( electron[i].DeltaR(bjet[j]) < 0.4){
        cout << "remove electron "<< i<<"!!"<<endl;
      }
      else { clean_el.push_back(electron[i]);
      }
    }
  }
  electron = clean_el;
  vector<int> mu_rm;
  vector<TLorentzVector> clean_mu;
  for (int i = 0; i < muon.size(); i++){
    for (int j = 0; j < bjet.size(); j++){
      if ( muon[i].DeltaR(bjet[j]) < 0.4){
        cout << "remove muon "<< i<<"!!"<<endl;
      }
      else {clean_mu.push_back(muon[i]);}
    }
  }
  muon = clean_mu;
  

// reclustered jets are R=0.4 jets re-clustered with anti-kt R=1.0 
// - similar to FatJET10, (untrimmed) but actually only the jets with pT>50 are used? is it very different?
  for (int i = 0; i < FatJet10_size; i++){
    //temp.SetPtEtaPhiM(FatJet10_Trimmed_PT[i],FatJet10_Trimmed_Eta[i],FatJet10_Trimmed_Phi[i],FatJet10_Trimmed_Mass[i]);
    temp1.SetPtEtaPhiM(FatJet10_PT[i],FatJet10_Eta[i],FatJet10_Phi[i],FatJet10_Mass[i]);
	if(temp1.Pt() > 100 && TMath::Abs(temp1.Eta()) < 1.5){
	  fatjet.push_back(temp1);
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

//calculate HT - only jets with pT>40 and eta<2.8
  	HT = 0;
	for (size_t i = 0; i < jet.size(); i++){
    	if( ( jet[i].Pt() > 40 ) && (TMath::Abs(jet[i].Eta())< 2.8 ) ){
    		//sum over scalar jet pT for all jets
    		HT += jet[i].Pt();
    	}
  	}
//calculate MJ - only jets with pT>40 and eta<2.8
	MJsum=0;
	for (size_t i = 0; i < fatjet.size(); i++){
    	MJsum += fatjet[i].M();
	}
// debugging
	if( ( electron.size()==0 ) && ( muon.size()==0 ) ){
		levent1++;
		if (jet50.size() >= 7) {
			levent2++;
			if ( missing_ET.Pt()/sqrt(HT) > 4 ){//cout << "Etmiss=" << missing_ET.Pt() << ", sqrt(HT)= "<<sqrt(HT) <<endl;
			   levent3++;
			}
		}
	}

// overall event cuts: ETmiss/Sqrt(HT) > 4, electrons or muons with pT>10 GeV are vetoed (Note: 10GeV cut unnecessary)
	if( ( missing_ET.Pt()/sqrt(HT) > 4 ) && ( electron.size()==0 ) && ( muon.size()==0 ) ){
		levent++;
		// SIGNAL REGIONS
		// 8-9-10j50 MJ340
		if(jet50.size() >= 8 && MJsum> 340){
		    //cout << "Etmiss=" << missing_ET.Pt() << ", sqrt(HT)= "<<sqrt(HT) <<endl;
			n8j50MJ340+=kfactor;	// 8j50
			if ( jet50.size() >= 9){
				n9j50MJ340+=kfactor; // 9j50
				if( jet50.size() >= 10) n10j50MJ340+=kfactor; // 10j50
			}
		}		
		if(jet50.size() >= 8 && MJsum> 500){
		    //cout << "Etmiss=" << missing_ET.Pt() << ", sqrt(HT)= "<<sqrt(HT) <<endl;
			n8j50MJ500+=kfactor;	// 8j50
			if ( jet50.size() >= 9){
				n9j50MJ500+=kfactor; // 9j50
				if( jet50.size() >= 10) n10j50MJ500+=kfactor; // 10j50
			}
		}		

	}
   return kTRUE;
}

void ATLAS_8_10_jet::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void ATLAS_8_10_jet::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


//cout << ngluino << endl;

if(human_output){
  cout << nevent << endl;
  cout << neventraw << endl;
  cout << levent << endl;
  cout << levent1 << endl;
  cout << levent2 << endl;
  cout << levent3 << endl;

}
else {
  cout << neventraw << endl;
  cout << nevent << endl;
  cout << "\t" << 1.0*n8j50MJ340/nevent;
  cout << "\t" << 1.0*n9j50MJ340/nevent;
  cout << "\t" << 1.0*n10j50MJ340/nevent;
  cout << "\t" << 1.0*n8j50MJ500/nevent;
  cout << "\t" << 1.0*n9j50MJ500/nevent;
  cout << "\t" << 1.0*n10j50MJ500/nevent << endl;

ofstream outfile;
  outfile.open("efficiencies_combined.dat", ofstream::app); 
  outfile << "ATLAS_8_10_jet";
  outfile << "\t"  << 1.0*n8j50MJ340/nevent;
  outfile << "\t" << 1.0*n9j50MJ340/nevent;
  outfile << "\t" << 1.0*n10j50MJ340/nevent;
  outfile << "\t" << 1.0*n8j50MJ500/nevent;
  outfile << "\t" << 1.0*n9j50MJ500/nevent; 
  outfile << "\t" << 1.0*n10j50MJ500/nevent << endl;
  outfile.close();
}
}  

