// this is a test for github
#ifndef openAnalyzer_H
#define openAnalyzer_H

// system include files
#include <memory>
#include <unistd.h>
#include <iostream>
#include <fstream>

#include <iostream>
#include <fstream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// GCT and RCT data formats
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "L1Trigger/L1TCaloLayer1/src/L1UCTCollections.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"


#include <memory>
#include <math.h>
#include <vector>
#include <list>

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"
class openAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit openAnalyzer(const edm::ParameterSet&);
  ~openAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------

  struct tTower{
    int iEta;
    int iPhi;
    int etaRight;
    int etaLeft;
    int phiUp;
    int phiLow;
  };

  struct tRegion{
    float pt;
    int iEta;
    int iPhi;
    float eta;
    float phi;
    float etaRight;
    float etaLeft;
    float phiUp;
    float phiLow;
    float hcalEnergy;
    float ecalEnergy;
    float calEnergy; 
  };



  edm::Service<TFileService> tfs_;

  std::ofstream file0, file1, file10;
  
  std::vector<float> vRegionEt;
  std::vector<float> vRegionEta;
  std::vector<float> vRegionPhi;
  std::vector<float> vRegionTau;
  std::vector<float> vRegionEG;
  std::vector<float> vRegionCal;
  std::vector<float> vRegionEcal;
  std::vector<float> vRegionHcal; 
  std::vector<float> vRegionInputEta; 
  std::vector<float> vRegionInputPhi; 


  TH1F* nEvents;

  TH1F* regionEta;
  TH1F* regionPhi;
  TH1F* regionPt;
  TH1F* regionEtaFine;
  TH1F* regionPhiFine;
  TH1F* regionCal; 
  TH1F* regionHcal;
  TH1F* regionEcal;
  TH1F* regionTotal;
  TH1F* regionInputEta; 
  TH1F* regionInputPhi; 

  TH1F* regionHitEta;
  TH1F* regionHitPhi;
  TTree* regionTree;
  TFileDirectory folder;

  int run, lumi, event;
  float nvtx;

  int convertEtaToTowerEta(double inputEta) {
    const double tpgEtaValues[27] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.848,
      1.956, // Last two HE bins are 0.432 and 0.828 wide
      2.064,
      2.172,
      2.379,
      2.586,
      2.793,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };


    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	//std::cout<<"found to be true"<<std::endl;
	//int tpgEta = n;
	//Positive eta is >28
	//negative eta is 0 to 27
	if(inputEta>0){
	  //std::cout<<"returning input eta >0 so + 28"<<std::endl;
	  return n + 28;}
	else{
	  //std::cout<<"returning input eta <0 so n"<<std::endl;
	  return n;}
	break;
      }
    }
        std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
  int convertPhiToTowerPhi(double inputPhi){
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  }

  /* This function takes as input the value of the 
     RCT Phi Region and then returns the radian decial phi value
   */
  
  float convertRCTPhiCentral(uint32_t inputPhi) {
    const double regionPhiValues[18] = {
					0.000, //0th region
					0.349, //1st region
					0.698, //2nd region

					1.047, //.... 
					1.396,
					1.744,
					2.093,
					2.442,
					2.791,
					-3.14159,
					-2.791,
					-2.442,
					-2.093,
					-1.744,
					-1.396,
					-1.047,
					-0.698,
					-0.349
    };
    return regionPhiValues[inputPhi];
  };

    /* This function takes as input the value of the 
       RCT Phi Region and then returns the radian decial phi value
   */

  //FINISHME
   float convertRCTPhiUpperBound(uint32_t inputPhi) {
   const double regionPhiValues[18] = {


					0.174, //0th region upper bound
					0.523, //1st region upper bound
					0.872, //2nd region upper bound
					1.221, //3rd region upper bound
					1.570, //4
					1.919, //5 
					2.268, //6
					2.617, //7
					2.967, //8
					-2.967,//9
 					-2.617,//10
					-2.268,//11
					-1.919, //12
					-1.570, //13 
					-1.221, //14
					-0.872, //15
					-0.523, //16 
					-0.174 //17 


   };


   return regionPhiValues[inputPhi];

   };

   
  //FINISHME

     
        float convertRCTPhiLowerBound(uint32_t inputPhi) {
   const double regionPhiValues[18] = {
					-0.174, //0
					0.174, //1th region 
					0.523, //2st region 
					0.872, //3nd region 
					1.221, //4rd region 
					1.570, //5
					1.919, //6
					2.268, //7
					2.617, //8 
					2.967, //9 
					-2.967, //10
					-2.617, //11
					-2.268, //12
					-1.919,//13 
					-1.570, //14
					-1.221,//15
					-0.872,//16
					-0.523 //17 
   };
   return regionPhiValues[inputPhi];
 };


  

  float convertRCTEtaCentral(uint32_t inputEta) {
    const double regionEtaValues[22] = {
      -4.75, //eta 0
      -4.25, //eta 1
      -3.75, //eta 2
      -3.25, //eta 3
      -2.5,  //eta 4
      -1.93, //eta 5
      -1.566, 
      -1.218,
      -0.87,
      -0.522,
      -0.174, //eta 10
      0.174,  // eta 11
      0.522, // 12 
      0.87, //13
      1.218, //14
      1.566, //15
      1.93, //16
      2.5, //17 
      3.25, //18 
      3.75, //19 
      4.25, //20 
      4.75 //21 
    };
    return regionEtaValues[inputEta];
  };

  float convertRCTEtaRightBound (uint32_t inputEta) {
    const double regionEtaValues[22] = {
      -4.5, //1 
      -4, //2
      -3.5, //3
      -3, //4
      -2.172, //5 
      -1.740, //6 
      -1.392, //7 
      -1.044, //8
      -0.696, //9 
      -0.348, //10
      0,//11
      0.348,//12
      0.696, //13
      1.044, //14
      1.392, //15
     1.74,
     2.172,
     3.00,
     3.50,
     4.00,
     4.50,
     5.00
    };

    return regionEtaValues[inputEta]; 
  };


  float convertRCTEtaLeftBound (uint32_t inputEta) {
    const double regionEtaValues[22] ={
      -5,  //1
      -4.5, //2 
      -4, //3 
      -3.5,// 4
      -3, //5 
      -2.172, //6 
      -1.740,//7 
      -1.392,//8 
      -1.044,// 9
      -0.696, //10
      -0.348, //11
      0, //12
      0.348, //13
      0.696, //14 
      1.044,
      1.392,
      1.74,
      2.172,
      3,
      3.5,
      4,
      4.5
    }; 

      return regionEtaValues[inputEta];	

  };



  
  edm::EDGetTokenT< vector <reco::PFCandidate> > pfToken_;  //used to select what ecal rechits to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
#endif
