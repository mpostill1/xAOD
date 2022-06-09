//Intialise functions
#include "xAODAnaHelpers/SlidingWindow.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include <cstdio>
#include "TMath.h"
#include "xAODAnaHelpers/JetContainer.h"
#include "xAODJet/JetContainer.h"
#include <vector>
#include <iostream>
std::vector<std::vector<double>> SlidingWindow::SlidingWindowExecutetopo(const xAOD::JetContainer* jetcont, const xAOD::CaloClusterContainer* clustercont){
 



//define variables
float pTower;

//Create a sliding $\eta$ window object of type !SlidingWindow::Window scanning $-4.9 < \eta < 4.9$ in 98 equal-distant steps of 0.1 with a window of width 0.8;


 SlidingWindow::Window etaWindows("MyTowerIntegrator fwd",98,-4.9,4.9,0.8);

//Loop on the towers and fill the etaWindows object with the towers like
 double ftarea(0.1*0.1);
 //double area = 0.1*0.1;          // assuming 0.1 x 0.1 tower grid
 // for (const xAOD::CaloCluster* pTower :*clustercont) {
 //if ( pTower->e() > 0. ) { etaWindows.addObject(*pTower,area); }
 //}
for ( const auto* cft : *clustercont ) { 
    if ( cft->e(xAOD::CaloCluster::UNCALIBRATED) > 0. ) { etaWindows.addObject(*cft,ftarea); }
 }
//Plot $\rho(\eta)$ by retrieving the median values:

 std::vector<double> medians(etaWindows.findMedians());
 SlidingWindow::Window::binning_t bins = etaWindows.binning();
 for ( size_t i(0); i<bins.size(); ++i ) {
   double eta(bins.at(i).binCenter());
   double rho(medians.at(etaWindows.binIndex(eta)));
    // ... plot ...
 }
	      
//FInd local $\rho$ for jets

 std::vector<double> medianpostill(etaWindows.findMedians());
 for (const xAOD::Jet* pJet :*jetcont ) {
   double rho(medianpostill.at(etaWindows.binIndex(pJet->eta())));
   double ptcorr=0;
   double jet_pt=pJet->pt();
   std::cout << "Check on pjet "<< pJet <<std::endl;
   std::cout << "Check on pjet eta "<< pJet->eta() <<std::endl;
   std::cout << "Check on etaWindows.binIndex(pJet->eta()) "<< etaWindows.binIndex(pJet->eta()) <<std::endl;
   std::cout << "Check on medianpostill.at(etaWindows.binIndex(pJet->eta()))) "<< medianpostill.at(etaWindows.binIndex(pJet->eta())) <<std::endl;
   double active_area=pJet->getAttribute<double>("ActiveArea");
   ptcorr=jet_pt-rho*active_area;
   std::cout <<"event density="<<rho<<std::endl;
    // ... keep/remove jet etc. ...
  }


//clear this vector if there's some junk in memory
  // final vector to return - this and anything associated to it won't be deleted!
 std::vector<std::vector<double>> myReturnVector;
 myReturnVector.clear();
 myReturnVector.shrink_to_fit();
 for (const xAOD::Jet* pJet :*jetcont ) {
   double eta=pJet->eta();
   double rho(medianpostill.at(etaWindows.binIndex(pJet->eta())));
   //clear this vector if there's some junk in memory
   std::vector<double> myEntryVec; // this object will be assigned to something in scope (our return vector) so will be saved from deletion
   myEntryVec.clear();
   if(std::fabs(eta)>2.5) {   
   
     myEntryVec.shrink_to_fit(); // not 100% sure on syntax on this - check it

     myEntryVec.push_back(eta); // push back appends to end of empty list, so remember the ordering
     myEntryVec.push_back(rho);
     //now push the entry vector to the return. Now the object in memory is assigned to the ReturnVector which now exists outside of the for loop.
     myReturnVector.push_back(myEntryVec);
   }

 }

 // check everything reads out:
 for(std::vector<double> Entry: myReturnVector){
   std::cout<<"Eta "<<Entry[0]<<" Rho "<<Entry[1]<<std::endl;
 }
 return myReturnVector;
}


