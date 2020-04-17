#include <cstdlib>
#include <iostream>
#include <chrono>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

#define clight 29.972 //cm/ns
#define DELTA_TIME 2 //ns
#define EBEAM 10.604
using namespace clas12;

TLorentzVector  Correct_Electron(TLorentzVector x){

  Double_t E_new, Px_el, Py_el, Pz_el;
  TLorentzVector el_new;

  E_new = x.E()-0.03689+0.1412*x.E()-0.04316*pow(x.E(),2)+0.007046*pow(x.E(),3)-0.0004055*pow(x.E(),4);

  Px_el = E_new*(x.Px()/x.Rho());
  Py_el = E_new*(x.Py()/x.Rho());
  Pz_el = E_new*(x.Pz()/x.Rho());

  el_new.SetXYZM(Px_el, Py_el, Pz_el, 0.000511);

  return el_new;
}







void test3(){
  HipoChain chain;
  vtp_ptr vtpp;
  chain.Add("/lustre19/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim2/skim2_*");

  gROOT->cd();
  TH1D *h_nElectrons=new TH1D("h_nElectrons","h_nElectrons",20,-0.5,19.5);
  TH1D *h_nPhotons=new TH1D("h_nPhotons","h_nPhotons",20,-0.5,19.5);
  TH1D *h_nProtons=new TH1D("h_nProtons","h_nProtons",20,-0.5,19.5);
  TH1D *h_dT=new TH1D("h_dT","h_dT",800,-50,50);
  TH1D *h_dTProt=new TH1D("h_dTProt","h_dTProt",800,-50,50);
  TH2D *h_dTvsE=new TH2D("h_dTvsE","h_dTvsE",800,-50,50,400,0,4);
  TH1D *h_nGood=new TH1D("h_nGood","h_nGood",5,-0.5,4.5);

  TH1D *h_M2G_all=new TH1D("h_M2G_all","h_M2G_all",500,0,1);
  TH1D *h_M2G_good=new TH1D("h_M2G_good","h_M2G_good",500,0,1);
  TH1D *h_M2GProt_good=new TH1D("h_M2GProt_good","h_M2GProt_good",500,0,1);
  TH1D *h_M2G_sel=new TH1D("h_M2G_sel","h_M2G_sel",500,0,1);
  TH1D *h_M2GProt_sel=new TH1D("h_M2GProt_sel","h_M2GProt_sel",500,0,1);

  TH1D *h_MM2_good=new TH1D("h_MM2_good","h_MM2_good",500,-.5,4.5);
  TH1D *h_MM2Prot_good=new TH1D("h_MM2Prot_good","h_MM2Prot_good",500,-.5,4.5);
  TH1D *h_MM2NewpeProt_good=new TH1D("h_MM2NewpeProt_good","h_MM2NewpeProt_good",500,-.5,4.5);
  TH1D *h_MM2_sel=new TH1D("h_MM2_sel","h_MM2_sel",500,-.5,4.5);
  TH1D *h_MM2Prot_sel=new TH1D("h_MM2Prot_sel","h_MM2Prot_sel",500,-.5,4.5);
  TH1D *h_MM2_selB=new TH1D("h_MM2_selB","h_MM2_selB",500,-.5,4.5);
  TH1D *h_MM2Prot_selB=new TH1D("h_MM2Prot_selB","h_MM2Prot_selB",500,-.5,4.5);

  TH1D *h_Trg_all=new TH1D("h_Trg_all","hTrg_all",32,-0.5,31.5);
  TH1D *h_Trg_good=new TH1D("h_Trg_good","hTrg_good",32,-0.5,31.5);
  TH1D *h_Trg_sel=new TH1D("h_Trg_sel","hTrg_sel",32,-0.5,31.5);
  TH1D *h_TrgProt_sel=new TH1D("h_TrgProt_sel","hTrgProt_sel",32,-0.5,31.5);

  TLorentzVector pbeam(0,0,EBEAM,EBEAM);
  TLorentzVector ptarget(0,0,0,0.9383);

  for(int ifile=0;ifile<chain.GetNFiles();ifile++){
    clas12reader c12(chain.GetFileName(ifile).Data());


    // c12.useFTBased();

    while(c12.next()==true){
      if(!c12.getDetParticles().empty()){

	vtpp=0;
	vtpp=c12.vtp();
	if (vtpp==0){
	  cout<<"NO VTP BANK?"<<endl;
	  continue;
	}

	auto particlesFT=c12.getByRegion(clas12::FT);
	auto nPhotons=0;
	auto nElectrons=0;
	auto idx=0;

	vector<int> idxPhotons;
	vector<int> idxElectrons;
	//Loop on all particles, count them
	for (auto particle : particlesFT){
	  switch (particle->getPid()){
	  case 11:
	    nElectrons++;
	    idxElectrons.push_back(idx);
	    break;
	  case 22:
	    nPhotons++;
	    idxPhotons.push_back(idx);
	  }
	  idx++;
	}
	h_nElectrons->Fill(nElectrons);
	h_nPhotons->Fill(nPhotons);


	auto particlesCD=c12.getByRegion(clas12::CD);
	auto nProtons=0;
	auto idp=0;

	vector<int> idpProtons;
	//Loop on all particles, count them
	for (auto particle : particlesCD){
	  switch (particle->getPid()){
	  case 2212:
	    nProtons++;
	    idpProtons.push_back(idp);
	  }
	  idp++;
	}
	h_nProtons->Fill(nProtons);

	//Ignore events with less than 1 e- and 2 g and 1 p(this is the exclusive case(see test3.c for inclusive))
	if (nElectrons<1) continue;
	if (nPhotons<2) continue;
	if (nProtons<1) continue;

	//Fill MM with ALL photon pairs
	for (int ig1=0;ig1<idxPhotons.size();ig1++){
	  for (int ig2=ig1+1;ig2<idxPhotons.size();ig2++){
	    TLorentzVector pg1,pg2;
	    pg1.SetXYZM(particlesFT[ig1]->par()->getPx(),particlesFT[ig1]->par()->getPy(),particlesFT[ig1]->par()->getPz(),0);
	    pg2.SetXYZM(particlesFT[ig2]->par()->getPx(),particlesFT[ig2]->par()->getPy(),particlesFT[ig2]->par()->getPz(),0);
	    h_M2G_all->Fill((pg1+pg2).M());
	  }
	}

	//check all trigger bits
	for (int ibit=0;ibit<32;ibit++){
	  if(c12.checkTriggerBit(ibit)) h_Trg_all->Fill(ibit);
	}

	//For each particle, construct a vector with the indexes of the other particles that are in time with this. Also include the particle itself
	vector<vector<int>> particleGroups;
	//Check all the relative times and group particles
	for (int ii=0;ii<particlesFT.size();ii++){
	  vector<int> inTimeWithThisParticle;
	  inTimeWithThisParticle.push_back(ii);
	  for (int jj=0;jj<particlesFT.size();jj++){
	    if (ii==jj) continue;
	    auto T1=particlesFT[ii]->ft(FTCAL)->getTime();
	    auto T2=particlesFT[jj]->ft(FTCAL)->getTime();
	    auto E1=particlesFT[ii]->ft(FTCAL)->getEnergy();
	    auto E2=particlesFT[jj]->ft(FTCAL)->getEnergy();
	    TVector3 p1(particlesFT[ii]->ft(FTCAL)->getX(),particlesFT[ii]->ft(FTCAL)->getY(),particlesFT[ii]->ft(FTCAL)->getZ());
	    TVector3 p2(particlesFT[jj]->ft(FTCAL)->getX(),particlesFT[jj]->ft(FTCAL)->getY(),particlesFT[jj]->ft(FTCAL)->getZ());
	    T1=T1-p1.Mag()/clight;
	    T2=T2-p2.Mag()/clight;
	    auto dT=T2-T1;

	    h_dT->Fill(dT);

	    if (fabs(dT)<DELTA_TIME){
	      inTimeWithThisParticle.push_back(jj);
	    }

	    //Fill dTvsE with the lower energy
	    (E2 < E1 ? h_dTvsE->Fill(dT,E2) : h_dTvsE->Fill(dT,E1));
	  }
	  particleGroups.push_back(inTimeWithThisParticle);
	}

	//For the moment, I did everything symmetrically, so that if particles A B C are in time, there are 3 groups.
	//Let's create one group only, keeping only the group IF the first particle has the minimum index in the group
	vector<vector<int>> particleGroupsReduced;
	for (int ii=0;ii<particleGroups.size();ii++){
	  auto idxFirst=particleGroups[ii][0];
	  if (idxFirst==*min_element(particleGroups[ii].begin(),particleGroups[ii].end())){
	    particleGroupsReduced.push_back(particleGroups[ii]);
	  }
	}

	/*
	  for (int ii=0;ii<particleGroups.size();ii++){
	  auto T=particlesFT[ii]->ft(FTCAL)->getTime();
	  TVector3 p(particlesFT[ii]->ft(FTCAL)->getX(),particlesFT[ii]->ft(FTCAL)->getY(),particlesFT[ii]->ft(FTCAL)->getZ());
	  T=T-p.Mag()/clight;
	  cout<<"GROUP FOR PARTICLE: "<<ii<<" id[ "<<particlesFT[ii]->getPid()<<"] Time: "<<T<<" : ";
	  for (auto idx : particleGroups[ii]) cout<<idx<<" ";
	  cout<<endl;
	  }


	  for (int ii=0;ii<particleGroupsReduced.size();ii++){
	  auto T=particlesFT[ii]->ft(FTCAL)->getTime();
	  TVector3 p(particlesFT[ii]->ft(FTCAL)->getX(),particlesFT[ii]->ft(FTCAL)->getY(),particlesFT[ii]->ft(FTCAL)->getZ());
	  T=T-p.Mag()/clight;
	  cout<<"GROUP-REDUCED FOR PARTICLE: "<<ii<<" id[ "<<particlesFT[ii]->getPid()<<"] Time: "<<T<<" : ";
	  for (auto idx : particleGroupsReduced[ii]) cout<<idx<<" ";
	  cout<<endl;
	  }*/


	//now start checking the groups
	//Let's check how many there are with 1 electron and two photons  //I should probably change this to check groups with more than 3 particles just in case there is some large constant background unrelated to the real event, something like low energy gamma or electron from radiated materials
	auto nGood=0;
	vector<vector<int>> particleGroupsGood;
	for (auto group : particleGroupsReduced){
	  if (group.size()>3{cout<<"Bananas="<<group.size <<endl;}
         if (group.size()==3){
           int nG=(particlesFT[group[0]]->getPid()==22)+(particlesFT[group[1]]->getPid()==22)+(particlesFT[group[2]]->getPid()==22);
           int nE=(particlesFT[group[0]]->getPid()==11)+(particlesFT[group[1]]->getPid()==11)+(particlesFT[group[2]]->getPid()==11);
           //if there are 2 photons and 1 electron, it is a good group. resort it so that the first particle is the e-
           if ((nG==2)&&(nE==1)){
             vector<int> goodGroup;
             if (particlesFT[group[0]]->getPid()==11){
               goodGroup.push_back(group[0]);
               goodGroup.push_back(group[1]);
               goodGroup.push_back(group[2]);
             }
             else if (particlesFT[group[1]]->getPid()==11){
               goodGroup.push_back(group[1]);
               goodGroup.push_back(group[2]);
               goodGroup.push_back(group[0]);
             }
             else{
               goodGroup.push_back(group[2]);
               goodGroup.push_back(group[0]);
               goodGroup.push_back(group[1]);
             }
             particleGroupsGood.push_back(goodGroup);
             nGood++;
           }
         }
       }
       h_nGood->Fill(nGood);

       if (nGood<1) continue;
       //check trigger bits good events
       for (int ibit=0;ibit<32;ibit++){
         if(c12.checkTriggerBit(ibit)) h_Trg_good->Fill(ibit);
       }

       //Fill MM with photon pairs from all the good groups.
       for (auto group : particleGroupsGood){
         TLorentzVector pg1,pg2,pe,Newpe;
         int ie=group[0];
         //The order is defined from before.
         int ig1=group[1];
         int ig2=group[2];
         pe.SetXYZM(particlesFT[ie]->par()->getPx(),particlesFT[ie]->par()->getPy(),particlesFT[ie]->par()->getPz(),0.511E-3);
         pg1.SetXYZM(particlesFT[ig1]->par()->getPx(),particlesFT[ig1]->par()->getPy(),particlesFT[ig1]->par()->getPz(),0);
         pg2.SetXYZM(particlesFT[ig2]->par()->getPx(),particlesFT[ig2]->par()->getPy(),particlesFT[ig2]->par()->getPz(),0);
	 Newpe = CorrectElectron(pe);
         auto Mgg=(pg1+pg2).M();
         h_M2G_good->Fill(Mgg);
         auto MM2=(pbeam+ptarget-pe-pg1-pg2);
         auto MM2Newpe=(pbeam+ptarget-Newpe-pg1-pg2);


         h_MM2_good->Fill(MM2.M2());
         h_MM2Newpe_good->Fill(MM2Newpe.M2());
         auto centerC=0.136;
         auto widthC=0.02;
         if (fabs(Mgg-centerC)<widthC){
           h_M2G_sel->Fill(Mgg);
           h_MM2_sel->Fill(MM2.M2());
	   h_MM2Newpe_sel->Fill(MM2Newpe.M2());
           //check trigger bits pi0 events
           for (int ibit=0;ibit<32;ibit++){
             if(c12.checkTriggerBit(ibit)) h_Trg_sel->Fill(ibit);
           }
         }
         if (fabs(Mgg-centerC+3*widthC/2)<widthC/2){
           h_MM2_selB->Fill(MM2.M2());
	   h_MM2Newpe_selB->Fill(MM2Newpe.M2());
         }
         if (fabs(Mgg-centerC-3*widthC/2)<widthC/2){
           h_MM2_selB->Fill(MM2.M2());
	   h_MM2Newpe_selB->Fill(MM2Newpe.M2());
         }
       }


//Start dealing with proton timing here too! Call it goodProt and bad Prot

//loop for protons timing etc, then within that a second loop for good groups and to compare. Take only the first group initially
	//loop round all particles in the CD(should be mostly protons since the requirement of a proton applied above)
       for (int ii=0;ii<particlesCD.size();ii++){

	 auto T1=particlesCD[ii]->sci(CTOF)->getTime();
	 auto E1=particlesCD[ii]->sci(CTOF)->getEnergy();
	 TVector3 p1(particlesCD[ii]->sci(CTOF)->getX(),particlesCD[ii]->sci(CTOF)->getY(),particlesCD[ii]->sci(CTOF)->getZ());//use getpathlength
	 auto B1=particlesCD[ii]->par()->getFTBBeta();
	 //auto B1=particles[ii]->par()->getBeta();
	 T1=T1-p1.Mag()/(B1*clight); //need beta
	 auto ProtPid=particlesCD[ii]->par()->getPid();
	 if(ProtPid==2212){
	   TLorentzVector pprot;
	   pprot.SetXYZM(particlesCD[ii]->par()->getPx(),particlesCD[ii]->par()->getPy(),particlesCD[ii]->par()->getPz(),0.93828);
	   for (auto group : particleGroupsGood){
	     TLorentzVector pg1,pg2,pe,Newpe;
	     int ie=group[0];
	     //The order is defined from before.
	     int ig1=group[1];
	     int ig2=group[2];
	     e.SetXYZM(particlesFT[ie]->par()->getPx(),particlesFT[ie]->par()->getPy(),particlesFT[ie]->par()->getPz(),0.511E-3);
	     pg1.SetXYZM(particlesFT[ig1]->par()->getPx(),particlesFT[ig1]->par()->getPy(),particlesFT[ig1]->par()->getPz(),0);
	     pg2.SetXYZM(particlesFT[ig2]->par()->getPx(),particlesFT[ig2]->par()->getPy(),particlesFT[ig2]->par()->getPz(),0);

	     auto T2=particlesFT[ie]->ft(FTCAL)->getTime();
	     auto T3=particlesFT[ig1]->ft(FTCAL)->getTime();
	     auto T4=particlesFT[ig2]->ft(FTCAL)->getTime();
	     auto E2=particlesFT[ie]->ft(FTCAL)->getEnergy();
	     auto E3=particlesFT[ig1]->ft(FTCAL)->getEnergy();
	     auto E4=particlesFT[ig2]->ft(FTCAL)->getEnergy();
	     TVector3 p2(particlesFT[ie]->ft(FTCAL)->getX(),particlesFT[ie]->ft(FTCAL)->getY(),particlesFT[ie]->ft(FTCAL)->getZ());
	     TVector3 p3(particlesFT[ig1]->ft(FTCAL)->getX(),particlesFT[ig1]->ft(FTCAL)->getY(),particlesFT[ig1]->ft(FTCAL)->getZ());
	     TVector3 p4(particlesFT[ig2]->ft(FTCAL)->getX(),particlesFT[ig2]->ft(FTCAL)->getY(),particlesFT[ig2]->ft(FTCAL)->getZ());
	     T2=T2-p2.Mag()/clight;
	     T3=T3-p3.Mag()/clight;
	     T4=T4-p4.Mag()/clight;
	     //auto dT=T3-T2;
	     if(E2>E3 && E2>E4){
	       auto dTProt=T1-T2;			
	     }	
	     else if(E3>E4){
	       auto dTProt=T1-T3;
	     }

	     else{
	       auto dTProt=T1-T4;
	     }
	     h_dTProt->Fill(dTProt); 				//Create all these histograms above
	     if (fabs(dTProt)<DELTA_TIME){
              
	       Newpe = CorrectElectron(pe);
	       auto MggProt=(pg1+pg2).M();
	       h_M2GProt_good->Fill(MggProt);
	       auto MM2Prot=(pbeam+ptarget-pe-pg1-pg2);
	       auto MM2NewpeProt=(pbeam+ptarget-Newpe-pg1-pg2);
	       h_MM2Prot_good->Fill(MM2Prot.M2());
	       h_MM2NewpeProt_good->Fill(MM2NewpeProt.M2());
	       auto centerC=0.136;
	       auto widthC=0.02;
	       if (fabs(MggProt-centerC)<widthC){
		 h_M2GProt_sel->Fill(MggProt);
		 h_MM2Prot_sel->Fill(MM2Prot.M2());
		 //check trigger bits pi0 events
		 for (int ibit=0;ibit<32;ibit++){
		   if(c12.checkTriggerBit(ibit)) h_TrgProt_sel->Fill(ibit);
		 }
	       }
	       if (fabs(MggProt-centerC+3*widthC/2)<widthC/2){
		 h_MM2Prot_selB->Fill(MM2Prot.M2());
		 h_MM2NewpeProt_selB->Fill(MM2NewpeProt.M2());
	       }
	       if (fabs(MggProt-centerC-3*widthC/2)<widthC/2){
		 h_MM2Prot_selB->Fill(MM2Prot.M2());
		 h_MM2NewpeProt_selB->Fill(MM2NewpeProt.M2());
	       }

       
	     }			   


				   
	   }

	 }
	  
       }


	}
      }
    }

 std::cout<<"DONE LOOP"<<std::endl;
 TFile *fout = new TFile("TEstout.root","recreate");

 gROOT->cd();
 TIter next(gDirectory->GetList());
 TObject* object = 0;
 fout->cd();
 while (object = next()) {
   if (object->InheritsFrom(TH1::Class())) {
     fout->cd();
     object->Write();
     std::cout << "Wrote " << object->GetName() << std::endl;
   }
 }                                                         
 fout->Close();
}



