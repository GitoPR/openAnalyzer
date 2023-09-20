#include <cstring>
#include <string>

void compAnalyzer_PlottingScript () {

  // Open the file
  std::unique_ptr<TFile> myFile(TFile::Open ("/nfs_scratch/jorgeeh/test_compAnalyzer.root"));


  //Get the file 
  TTree* efficiencyTree = (TTree*) myFile -> Get("demo/EfficiencyTree"); 
  
  

  TH1F *Cal = new TH1F("Cal", "Cal_Energy", 20, 0.0 ,20.0 );  
  //  TH1F *script_compRegionEt = new TH1F ("vcompRegionEt", "Transverse_Energy", 50, 0 ,3);
  // TH1F *script_compRegionPhi = new TH1F ("vcompRegionPhi" , "Phi" , 50 , 0, 3);
  //  TH1F *script_compRegionEta = new TH1F ("vcompRegionEta" , "Eta", 50 , 0 , 3);  

  std::vector<float>* vCal = nullptr;
  

  efficiencyTree->SetBranchAddress("vRegionCal",  &vCal);

  // Draw one of the branches 



  Int_t nentries = (Int_t)efficiencyTree->GetEntries();
  //  std::cout << "number of entries :" << nentries << std::endl;

  for(Int_t i = 1; i<nentries; i++){
    // std::cout << "Entry :" << i << std::endl;
    efficiencyTree -> GetEntry(i); 
    //    std::cout << "check to see if the code got here" << std::endl;
  for(auto vi = vCal->begin(); vi != vCal->end(); vi++) { 

  	double value = *vi ; 	
	// 		std::cout << "value :" << value << std::endl;
	
   	if (!std::isnan(value)){
	  if (value > 0){
	    Cal->Fill(value);
	  
	  }  
  	  //std::cout << "Filling Th1F with " << value << std::endl;
   	} 
    }
  }
  //  efficiencyTree -> Draw("vcompRegionPhi" , "script_compRegionPhi" );
  // efficiencyTree -> Draw("vcompRegionEta", "script_compRegionEta");
  //  efficiencyTree -> Draw("vcompRegionEt","script_compRegionEt")  


  // Create TCanvas
  TCanvas* TCan = new TCanvas("TCan", "" , 100, 20, 800, 600) ;  

  //Generally we need this line before we draw anything
  TCan->cd(); 

  //Draw the histogram 
  Cal->SetMarkerColor(0); 
  Cal->SetLineWidth(1); 
  Cal->SetFillStyle(1000); 



  Cal->Draw("HIST"); 
  TCan->SaveAs("./compAnalyzer_testplot3.pdf"); 

    //Clean up 
    delete TCan; 

}  
