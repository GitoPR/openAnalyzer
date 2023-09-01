#include <cstring>
#include <string>

void compAnalyzer_PlottingScript () {

  // Open the file
  std::unique_ptr<TFile> myFile(TFile::Open ("/nfs_scratch/jorgeeh/compAnalyzer_MC2018.root"));


  //Get the file 
  TTree* efficiencyTree = (TTree*) myFile -> Get("demo/EfficiencyTree"); 


  TH1F *script_Ecal = new TH1F("vRegionEcal", "Ecal_Energy", 50, 0 , 3 );  
  //  TH1F *script_compRegionEt = new TH1F ("vcompRegionEt", "Transverse_Energy", 50, 0 ,3);
  // TH1F *script_compRegionPhi = new TH1F ("vcompRegionPhi" , "Phi" , 50 , 0, 3);
  //  TH1F *script_compRegionEta = new TH1F ("vcompRegionEta" , "Eta", 50 , 0 , 3);  


  // Draw one of the branches 
  efficiencyTree -> Draw("vRegionEcal >> script_RegionEcal");
  //  efficiencyTree -> Draw("vcompRegionPhi" , "script_compRegionPhi" );
  // efficiencyTree -> Draw("vcompRegionEta", "script_compRegionEta");
  //  efficiencyTree -> Draw("vcompRegionEt","script_compRegionEt")  


  // Create TCanvas
  TCanvas* TCan = new TCanvas("TCan", "" , 100, 20, 800, 600) ;  

  //Generally we need this line before we draw anything
  TCan -> cd(); 

  //Draw the histogram 
  script_Ecal -> SetMarkerColor(0); 
  script_Ecal -> SetLineWidth(1); 
  script_Ecal -> SetFillStyle(1001); 



  script_Ecal -> Draw("Ecal"); 
  TCan -> SaveAs ("../plot/compAnalyzer_testplot.pdf"); 

    //Clean up 
    delete TCan; 

}  
