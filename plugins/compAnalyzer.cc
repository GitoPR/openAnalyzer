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
// test comment
 
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CICADA/openAnalyzer/plugins/compAnalyzer.h" 

// To calculate the position at ECAL 
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "TLorentzVector.h"  
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::PFRecHitCollection;
using std::vector;


compAnalyzer::compAnalyzer(const edm::ParameterSet& iConfig)
  : pfToken_(consumes<vector <pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfcands")))
    //    regionSource_(consumes<vector<L1CaloRegion> >(iConfig.getParameter<edm::InputTag>("L1CaloRegion")))



 {
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
 
  //initialize the pfcand regionTree. 

  regionTree = tfs_->make<TTree>("EfficiencyTree", "Efficiency Tree");
  regionTree->Branch("run",        &run,     "run/I");
  regionTree->Branch("lumi",       &lumi,    "lumi/I");
  regionTree->Branch("event",      &event,   "event/I");
  regionTree->Branch("nvtx",       &nvtx,    "nvtx/D");
  regionTree->Branch("vRegionEt",  &vRegionEt  );
  regionTree->Branch("vRegionEta", &vRegionEta );
  regionTree->Branch("vRegionPhi", &vRegionPhi );
  regionTree->Branch("vRegionEG",  &vRegionEG  );
  regionTree->Branch("vRegionTau", &vRegionTau );
  regionTree->Branch("vRegionInputPhi", &vRegionInputPhi) ; 
  regionTree->Branch("vRegionInputEta", &vRegionInputEta ) ; 

  regionTree->Branch("vRegionCal", &vRegionCal );
  regionTree->Branch("vRegionEcal", &vRegionEcal );
  regionTree->Branch("vRegionHcal", &vRegionHcal );
  
  regionTree -> Branch ("vcompRegionEta", &vcompRegionEta);
  regionTree -> Branch ("vcompRegionPhi" , &vcompRegionPhi);
  regionTree -> Branch ("vcompRegionEt", &vcompRegionEt);
  


  nEvents       = tfs_->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );

  regionHitEta  = tfs_->make<TH1F>( "regionHit_eta"  , "eta", 16, 1, 16. );
  regionHitPhi  = tfs_->make<TH1F>( "regionHit_phi"  , "phi", 16, 1, 16. );
  regionTotal   = tfs_->make<TH1F>( "regionHit_total"  , "fullmap", 16, 1, 16. );

  regionEta     = tfs_->make<TH1F>( "region_eta"  , "eta", 22, 1, 22. );
  regionPhi     = tfs_->make<TH1F>( "region_phi"  , "phi", 72, 1, 72. );
  regionPt      = tfs_->make<TH1F>( "region_pt"  , "pt", 100, 0, 100. );
  regionCompEta =  tfs_->make<TH1F>( "region_compEta"  , "compEta", 100, 0, 100. );  
  regionCompPhi =  tfs_->make<TH1F>( "region_compPhi"  , "compPhi", 100, 0, 100. ); 
  regionCompEt =  tfs_->make<TH1F>( "region_compEt"  , "compEt", 100, 0, 100. );  


  regionCal     = tfs_->make<TH1F>("region_cal", "Cal", 10, 0 ,10);  
  regionHcal    = tfs_->make<TH1F>("region_hcal","Hcal", 10, 0, 10 ) ;
  regionEcal    = tfs_->make<TH1F>("region_ecal"," Ecal", 10 , 0 , 10) ; 
 
  
  regionEtaFine   = tfs_->make<TH1F>( "region_eta_Fine"  , "eta", 88, 1, 88. );
  regionPhiFine   = tfs_->make<TH1F>( "region_phi_Fine"  , "phi", 72, 1, 72. );

  regionInputEta     = tfs_->make<TH1F>( "region_ieta"  , "ieta", 22, 1, 22. ); 
  regionInputPhi     = tfs_->make<TH1F>( "region_iphi"  , "iphi", 22, 1, 22. ); 

  l1CaloHist  = tfs_->make<TH2F>("l1CaloHist", "Eta/Phi Et distribution", 22, 0, 22 , 18,0, 18 );
  pfcCalHist = tfs_->make<TH2F>("ParticleCandidateHist", "Eta vs. Phi Cal distribution;eta;phi", 22,0,22,18,0,18); 
  pfcEcalHist = tfs_->make<TH2F>("pfcEcalHist", "Eta vs. Phi Ecal distribution;eta;phi", 20,-3,3,20,-3.14,3.14);
  pfcHcalHist = tfs_->make<TH2F>("pfcHcalHist", "Eta vs. Phi Hcal distribution;eta;phi", 20,-3,3,20,-3.14,3.14);
} 

compAnalyzer::~compAnalyzer() {
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

void compAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
 
    run = iEvent.id().run();
    lumi = iEvent.id().luminosityBlock();
    event = iEvent.id().event();

     
  //gotta push back the run, lumi, event information into the corresponding ttree vector. 

  //edm::Handle<vector <pat::PackedCandidate> > pfCands;
  edm::Handle<pat::PackedCandidateCollection> pfCands;
  if(!iEvent.getByToken(pfToken_, pfCands))
    std::cout<<"ERROR GETTING THE PFCANDS"<<std::endl;
  std::cout<<"size of pfCands: "<<pfCands->size()<<std::endl;

  /*  edm::Handle<vector <L1CaloRegion> > regions;
  if(!iEvent.getByToken(regionSource_, regions))
    std::cout <<"Error GETTING regions" << std:: endl; 
    std::cout <<"size of regions" << regions -> size() << std::endl; */ 


  //First, create the regions
  tRegion allRegions[22][18];

  for(int ieta = 0; ieta < 22; ieta++){
    for(int iphi = 0; iphi < 18; iphi++){
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
      
      //      std::cout<< "iphi" << allRegions[ieta][iphi].phi << std::endl;
    }
  }

  //region Vectors for pfCand analyzer


  //clear the region vectors
  vRegionEt.clear();
  vRegionEta.clear();
  vRegionPhi.clear();
  vRegionCal.clear();
  vRegionEcal.clear();
  vRegionHcal.clear();
  vRegionInputEta.clear();
  vRegionInputPhi.clear();
  //Tau and EG are not needed currently
  //vRegionTau.clear();
  // vRegionEG.clear();
  //  propagatedPF.clear() ; 
  // std::cout << "check1"  << std::endl;
  for (const auto& pf : *pfCands.product()) {

     std::cout << "pf:" << pf  << std::endl;

    //    reco::CandidatePtr packedcand = pf.sourceCandidatePtr(k);
    //    const reco::PFCandidate* pf_Cand = dynamic_cast<const reco::PFCandidate*>(&(*packedcand));
    
    // Construct Particle object 
    RawParticle particle(pf.p4());  	 
    particle.setVertex(pf.vertex().x(), pf.vertex().y(), pf.vertex().z(), 0); 
    particle.setMass(pf.mass()); 


    // Construct BaseParticlePropagator object and call to propagate to Ecal Entrance. 
      float field_z = 4.; 
    BaseParticlePropagator prop(particle, 0. ,0. , field_z);  // Do the last three inputs of the function describe the physics under which the particle propagates? Why does it stay the same?
    prop.propagateToEcalEntrance(); 

    TLorentzVector corrPF;

    if (prop.getSuccess() != 0){

      GlobalPoint ecal_pos(prop.particle().vertex().x(), prop.particle().vertex().y() , prop.particle().vertex().z());   

      corrPF.SetPtEtaPhiM(prop.particle().Pt(), 
				   ecal_pos.eta(), 
				   ecal_pos.phi(), 
				   prop.particle().mass()); 

    /*  std::cout << "Corrected PF pt/eta/phi"
		<< pf.pt() << ", " << pf.eta() << ", " << pf.phi()
		<< "to" 
		<< corrPF.Pt() << ", " << corrPF.Eta() << ", " << corrPF.Phi()
		<< "position at ecal  eta/phi" 
		<< ecal_pos.eta()<< ", " << ecal_pos.phi() 
		<< std::endl; */ 
 
	// I think I might not need this last vector but I'll keep it and ask about it.      
	// propagatedPF.push_back(corrPF);  
	} 

   
    // get HCAL and ECAL energy
    //After a few difficulties with trying to access the energies this is what we think works. 
    float calEnergy  = pf.caloFraction()* pf.energy(); 
    float hcalEnergy = calEnergy * pf.hcalFraction()  ; 
    float ecalEnergy = calEnergy - hcalEnergy;  
    std::cout << "pf.caloFractio: " << pf.caloFraction() << std::endl; 
    std::cout << "pf.energy: " << pf.energy() << std::endl; 
    std::cout <<"pf.hcalFraction: " << pf.hcalFraction() << std::endl; 
    //     std::cout << "pf" << pf<< std::endl;
    //     std::cout << "pfc calEnergy:" << pf.caloFraction()* pf.energy() << std::endl;
    //       std::cout << "pfc ecalEnergy:" << ecalEnergy << std::endl; 
    //std::cout << "pfc hcalEnergy:" << hcalEnergy << std::endl; 


    //find the appropriate region and add in region energy to the tRegion
    //We're grabbing these from the corrected/calculated propagated particles
    
    float eta = corrPF.Eta(); 
    float phi = corrPF.Phi(); 
    //       count  = count + 1 ;
    //std::cout << "pfc eta:" << eta << std::endl;
    // std::cout << "pfc phi:" << phi << std::endl;
    //       std::cout << "\n" << std::endl;
    
    int Region_ieta = 0; 
    int Region_iphi = 0; 

    // Iterate through regions, find the bounds for given ieta/iphi then store region number. 
    for(int ieta = 0; ieta < 22; ieta ++ ) {
      if (eta > convertRCTEtaLeftBound(ieta) && eta < convertRCTEtaRightBound(ieta)){
	Region_ieta = ieta; 
	
      }
    }
    
    for(int iphi = 0 ; iphi < 18 ; iphi++) {
      if ( phi > convertRCTPhiLowerBound(iphi) && phi < convertRCTPhiUpperBound(iphi)){
	Region_iphi = iphi; 
	//	               std::cout<<"Region Iphi: "<< Region_iphi<<std::endl;
      }
    }    
    // Exception to account for the special case of region 9. 
    if(phi > 2.967059728 || phi < -2.967059728){
      Region_iphi = 9;
    }
    
    
    
    //        std::cout<<"Region Iphi: "<< Region_iphi<<std::endl;                                                                                   
    //	std::cout<<"Region ieta: "<< Region_ieta<<std::endl;  
    
    
    
    // Assign the hcal, ecal , and cal energies to the appropriate, just found region 
    allRegions[Region_ieta][Region_iphi].hcalEnergy += hcalEnergy;
    allRegions[Region_ieta][Region_iphi].ecalEnergy += ecalEnergy;
    allRegions[Region_ieta][Region_iphi].calEnergy += calEnergy;
    

    //Some Checks
    /*  std::cout<< "ecal energy:" << allRegions[Region_ieta][Region_iphi].ecalEnergy<< std::endl; 
	std::cout<< "hcal energy:" << allRegions[Region_ieta][Region_iphi].hcalEnergy<< std::endl;
	std::cout<< "cal energy:" << allRegions[Region_ieta][Region_iphi].calEnergy<<std::endl; */ 
    

    ///////ignore for now
    //if(i<20){    //std::cout<<"tower eta, tower phi: "<<towerEta<<", "<<towerPhi<<std::endl;
 //i++;
    //}
    // Save the following lines for now
    //EcalEBTriggerPrimitiveDigi tempDigi;
    //tempDigi.id().eta() = towerEta;
    //std::cout<<"digi towereta: "<<tempDigi.id().eta()<<std::endl;
    //std::cout<<"pfCand HCAL Energy: "<< pfCand.rawHcalEnergy()<<std::endl;

      }

  /*for(auto region: allRegions){

    
    //     std::cout<<"region eta,phi: "<<region->eta<<" , "<<region->phi<<std::endl;


    
    for(int ieta = 0; ieta < 22; ieta++){
      for(int iphi = 0; iphi < 18; iphi++){

    std::cout<<"phi: "<<allRegions[ieta][iphi].phi<<std::endl;
    std::cout<<"hcal: "<<allRegions[ieta][iphi].hcalEnergy<<std::endl;
      }
      } 

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
      regionTree->Fill(); */ 

  for(int ieta = 0; ieta < 22; ieta++){                                                                                                                                                                   
    for(int iphi = 0; iphi < 18; iphi++){ 

      //             std::cout<<"hcal: "<<allRegions[ieta][iphi].hcalEnergy<<std::endl; 

      vRegionEt.push_back(allRegions[ieta][iphi].pt);                                               
      vRegionEta.push_back(allRegions[ieta][iphi].eta);                                                  
      vRegionPhi.push_back(allRegions[ieta][iphi].phi);                                                
      vRegionCal.push_back(allRegions[ieta][iphi].calEnergy);   
      vRegionEcal.push_back(allRegions[ieta][iphi].ecalEnergy);                                           
      vRegionHcal.push_back(allRegions[ieta][iphi].hcalEnergy); 
      vRegionInputEta.push_back(allRegions[ieta][iphi].iEta);
      vRegionInputPhi.push_back(allRegions[ieta][iphi].iPhi); 
      
      pfcCalHist->Fill(allRegions[ieta][iphi].iEta,allRegions[ieta][iphi].iPhi,allRegions[ieta][iphi].calEnergy); 
      pfcEcalHist->Fill(allRegions[ieta][iphi].eta,allRegions[ieta][iphi].phi,allRegions[ieta][iphi].ecalEnergy);
      pfcHcalHist->Fill(allRegions[ieta][iphi].eta,allRegions[ieta][iphi].phi,allRegions[ieta][iphi].hcalEnergy);

      regionPt->Fill(allRegions[ieta][iphi].pt);                                                        
      regionEta->Fill(allRegions[ieta][iphi].eta);             
      regionPhi->Fill(allRegions[ieta][iphi].phi);
      //std::cout << "allRegions[ieta][iphi].calEnergy : " << allRegions[ieta][iphi].calEnergy << std::endl;
      regionCal->Fill(allRegions[ieta][iphi].calEnergy);  
      
      regionEcal->Fill(allRegions[ieta][iphi].ecalEnergy);                                  
      regionHcal->Fill(allRegions[ieta][iphi].hcalEnergy);
      regionInputEta->Fill(allRegions[ieta][iphi].iEta); 
      regionInputPhi->Fill(allRegions[ieta][iphi].iPhi); 

      //std::cout<<"cal: "<<allRegions[ieta][iphi].calEnergy<<std::endl;
      //std::cout << "iEta: " << allRegions[ieta][iphi].iEta << std:: endl; 
      //std::cout<< "iPhi: " << allRegions[ieta][iphi].iPhi << std:: endl;

      
 }

}

       regionTree->Fill();

       //       for(auto vi = std::begin(vRegionCal); vi != std::end(vRegionCal) ; vi++) {

       //double value = *vi; 
	   
       //  std::cout << "value :" << value << std::endl;   
       
       //  }
       /*
       uint32_t count = 0  ; 

      for (vector<L1CaloRegion>::const_iterator testRegion = regions->begin();  testRegion != regions -> end(); ++testRegion){

	
	count = count +1 ; 
	
	uint32_t test_et = testRegion -> et(); 
	uint32_t test_rEta =  testRegion -> id().ieta(); 
	uint32_t test_rPhi = testRegion -> id().iphi();  


	std::cout << "test_et :  " << test_et << std::endl;	
	std::cout << "test_reta :  " << test_rEta << std::endl; 
	std::cout << "test_rPhi:  " << test_rPhi << std::endl;
	std::cout << "\n" << std::endl; 
	std::cout << "count :" << count << std::endl; 
	std::cout << "\n" << std::endl;

	vcompRegionEta.push_back(test_rEta); 
	vcompRegionPhi.push_back(test_rPhi);
	vcompRegionEt.push_back(test_et);
	l1CaloHist-> Fill(test_rEta,test_rPhi, test_et);
 
	
  }

 */ 


      /*
      for(auto vi = std::begin(vcompRegionEta); vi != std::end(vcompRegionEta) ; vi++) {                                                                                                                                                                                      

      double value = *vi;                                                                                                                                                                                                                                                    

      std::cout << "vcompRegionEta :" << value  << std::endl;                                                                                                                                                                                                                        

        }    
      */ 


       //   regionTree->Fill(); 
       

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void compAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void compAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void compAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(compAnalyzer);
