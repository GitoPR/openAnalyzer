// -*- C++ -*-
//
// Package:    CICADA/openAnalyzer
// Class:      openAnalyzer
//
/**\class openAnalyzer openAnalyzer.cc CICADA/openAnalyzer/plugins/openAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Wed, 10 May 2023 15:15:24 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "CICADA/openAnalyzer/plugins/openAnalyzer.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::PFRecHitCollection;
using std::vector;


openAnalyzer::openAnalyzer(const edm::ParameterSet& iConfig)
  : pfToken_(consumes<vector <reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("pfcands"))) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  
  folder = tfs_-> mkdir("data");

  
  
  regionTree = folder.make<TTree>("EfficiencyTree", "Efficiency Tree");
  regionTree->Branch("run",        &run,     "run/I");
  regionTree->Branch("lumi",       &lumi,    "lumi/I");
  regionTree->Branch("event",      &event,   "event/I");
  regionTree->Branch("nvtx",       &nvtx,    "nvtx/D");
  regionTree->Branch("vRegionEt",  &vRegionEt  );
  regionTree->Branch("vRegionEta", &vRegionEta );
  regionTree->Branch("vRegionPhi", &vRegionPhi );
  regionTree->Branch("vRegionEG",  &vRegionEG  );
  regionTree->Branch("vRegionTau", &vRegionTau );
  

  regionTree->Branch("vRegionCal", &vRegionCal );
  regionTree->Branch("vRegionEcal", &vRegionEcal );
  regionTree->Branch("vRegionHcal", &vRegionHcal );
  
  nEvents       = folder.make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );

  regionHitEta  = folder.make<TH1F>( "regionHit_eta"  , "eta", 16, 1, 16. );
  regionHitPhi  = folder.make<TH1F>( "regionHit_phi"  , "phi", 16, 1, 16. );
  regionTotal   = folder.make<TH1F>( "regionHit_total"  , "fullmap", 16, 1, 16. );

  regionEta     = folder.make<TH1F>( "region_eta"  , "eta", 22, 1, 22. );
  regionPhi     = folder.make<TH1F>( "region_phi"  , "phi", 72, 1, 72. );
  regionPt      = folder.make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );
 
  regionCal     = folder.make<TH1F>("region_cal", "Cal", 10, 0 ,100);  
  regionHcal    = folder.make<TH1F>("region_hcal","Ecal", 10, 0, 100 ) ;
  regionEcal    = folder.make<TH1F>("region_ecal"," Hcal", 10 , 0 , 100) ; 
 
regionEtaFine   = folder.make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
regionPhiFine   = folder.make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );

  


}

openAnalyzer::~openAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
/* In this analyzer we'd like to convert pfcandidate raweta/rawphi to tower eta/tower phi
 */

void openAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
 
  edm::Handle<vector <reco::PFCandidate> > pfCands;
  if(!iEvent.getByToken(pfToken_, pfCands))
    std::cout<<"ERROR GETTING THE PFCANDS"<<std::endl;
  std::cout<<"size of pfCands: "<<pfCands->size()<<std::endl;

  //First, create the regions
  tRegion allRegions[22][18];

  for(int ieta = 0; ieta < 21; ieta++){
    for(int iphi = 0; iphi < 17; iphi++){
      allRegions[ieta][iphi].pt = 0;
      allRegions[ieta][iphi].iEta = ieta;
      allRegions[ieta][iphi].iPhi = iphi;
      allRegions[ieta][iphi].eta = convertRCTEtaCentral(ieta);
      allRegions[ieta][iphi].phi = convertRCTPhiCentral(iphi);
      allRegions[ieta][iphi].phiUp = convertRCTPhiUpperBound(iphi);
      allRegions[ieta][iphi].phiLow = convertRCTPhiLowerBound(iphi);
      allRegions[ieta][iphi].etaLeft = convertRCTEtaLeftBound(ieta);
      allRegions[ieta][iphi].etaRight = convertRCTEtaRightBound(ieta);
      allRegions[ieta][iphi].hcalEnergy = 0; 
      allRegions[ieta][iphi].ecalEnergy = 0;
      allRegions[ieta][iphi].calEnergy = 0;
      

                                                                                                              



 // COMPLETE ME: add in the left and right values for each eta and phi value

    }
  }
  


  //clear the region vectors
  vRegionEt.clear();
  vRegionEta.clear();
  vRegionPhi.clear();
  vRegionCal.clear();
  vRegionEcal.clear();
  vRegionHcal.clear();
  //Tau and EG are not needed currently
  vRegionTau.clear();
  vRegionEG.clear();


  int i = 0;
  for (const auto& pfCand : iEvent.get(pfToken_)) {
    // get HCAL and ECAL energy
    float hcalEnergy = pfCand.rawHcalEnergy();
    float ecalEnergy = pfCand.rawEcalEnergy();
    //Figure out where the energy was deposited
    int towerEta = convertEtaToTowerEta(pfCand.positionAtECALEntrance().eta());
    int towerPhi = convertPhiToTowerPhi(pfCand.positionAtECALEntrance().phi()); 
    //find the appropriate region and add in region energy to the tRegion
                                                                                                              
    std::cout<< "ecal energy:" << ecalEnergy<< std::endl;                                                                                                        
    std::cout<< "hcal energy:" << hcalEnergy<< std::endl;    

    //Add that to the correct tRegion in the AllRegions collection
	
	int eta = pfCand.positionAtECALEntrance().eta();
	int phi = pfCand.positionAtECALEntrance().phi();
	int Region_ieta = 0; 
	int Region_iphi = 0; 

   

	for(int ieta = 0; ieta < 22; ieta ++ ) {

	  if (eta < convertRCTEtaRightBound(ieta) && eta > convertRCTEtaLeftBound(ieta)){
	    
	    Region_ieta = ieta; 
	       break; 
	  }

	}
	  for(int iphi = 0 ; iphi < 18 ; iphi++) {

	    if ( phi > convertRCTPhiLowerBound(iphi) && phi < convertRCTPhiUpperBound(iphi)){
	      
	      Region_iphi = iphi; 
	      break;
	      
	    }
      }    

	  allRegions[Region_ieta][Region_iphi].hcalEnergy = hcalEnergy;
	  allRegions[Region_ieta][Region_iphi].ecalEnergy = ecalEnergy;
  
	  allRegions[Region_ieta][Region_iphi].calEnergy = hcalEnergy + ecalEnergy; 

	  /*
	  std::cout<< "ecal energy:" << allRegions[Region_ieta][Region_iphi].ecalEnergy<< std::endl; 
	  std::cout<< "hcal energy:" << allRegions[Region_ieta][Region_iphi].hcalEnergy<< std::endl;
	  std::cout<< "cal energy:" << allRegions[Region_ieta][Region_iphi].calEnergy<<std::endl;
	  */

    ///////ignore for now
    //if(i<20){
    //std::cout<<"tower eta, tower phi: "<<towerEta<<", "<<towerPhi<<std::endl;
    //i++;
    //}
    // Save the following lines for now
    //EcalEBTriggerPrimitiveDigi tempDigi;
    //tempDigi.id().eta() = towerEta;
    //std::cout<<"digi towereta: "<<tempDigi.id().eta()<<std::endl;
    //std::cout<<"pfCand HCAL Energy: "<< pfCand.rawHcalEnergy()<<std::endl;
  }


  for(auto region: allRegions){
    
    //std::cout<<"region eta,phi: "<<eta<<" , "<<phi<<std::endl;
    vRegionEt.push_back(region->pt);
    vRegionEta.push_back(region->eta);
    vRegionPhi.push_back(region->phi);
    vRegionCal.push_back(region->calEnergy); 
    vRegionEcal.push_back(region->ecalEnergy); 
    vRegionHcal.push_back(region->hcalEnergy);

    //we don't have this for now
    //vRegionEG.push_back(isEgammaLike);
    //vRegionEG.push_back(isTauLike);
     
    regionPt->Fill(region->pt);
    regionEta->Fill(region->eta);
    regionPhi->Fill(region->phi);
    regionCal->Fill(region->calEnergy);
    regionEcal->Fill(region->hcalEnergy);
    regionHcal->Fill(region->ecalEnergy);  

}
  regionTree->Fill();

  //I think this might be the proble, cause the Tree branches dont include vRegionCal and so on// 

 


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void openAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void openAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void openAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(openAnalyzer);
