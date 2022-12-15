// The manual on fastjet.fr is very useful (click on the manual tab and look for specific items in the document).
// I also find the doxygen page very useful to look up interfaces and implementation details (click on the doxygen
// tab and use the search for quick access). 
// 
//
#include "fastjet/PseudoJet.hh"            // FastJet EDM (Sect. 3.1, p. 11ff)
//#include "xAODAnaHelpers/fastjet.h"
#include "fastjet/ClusterSequenceArea.hh"  // clustering jet with area (p. 38ff)
#include "fastjet/AreaDefinition.hh"       // area definition          (p. 38ff)
#include "fastjet/JetDefinition.hh"        // jet definition           (p. 21ff)
#include "xAODAnaHelpers/fastjetmp.h"
////////////////////////////////////////////////////////////////////////
// Define a class that links the fastjet::PseudoJet with xAOD objects //
////////////////////////////////////////////////////////////////////////

// I suggest to use a templated class here so that you can easily link pFlow or CaloCluster
// or whatever and retrieve it without type casting. Note that the base class
// fastjet::PseudoJet::UserInfoBase does not require any purely virtual interface to be implemented.
//////////////////////////////////////////////////////////////
// Example function to get areas from four-momentum type objects //
///////////////////////////////////////////////////////////////////

// Example for topo-clusters:
//
// Two-step approach: first cluster inclusive (no pt- or rapidity selection) jets and then get cluster areas.
//
// Assumes that the xAOD::CaloClusterContainer has been picked up elsewhere.


// Retrieve all topo-clusters (pflows etc.) clustered into jets. The return PseudoJets each have a link to the original xAOD object
// and an area.  


///////////////////////////////////
// How to use with SlidingWindow //
///////////////////////////////////

//bool yourImplementation() {
  // pick up topo-clustres and forward towers from somewhere...
  //xAOD::CaloClusterContainer* pClusCont = ...; // topo-clusters
  //xAOD::CaloClusterContainer* pFwdTCont = ...; // forward towers
  // get clusters with area
  //std::vector<fastjet::PseudoJet> cpjets = getConstituentsWithAreas(*pClusCont);
  // fill sliding window (with topo-clusters and topo-towers)
  //SlidingWindow::Window rapWindow("Rho",98,-4.9,4.9,0.8);
  //double ftarea(0.1*0.1);
  //double rapmax(2.5); 
  // -- topo-clusters within |y|<2.5
  //for ( const auto& cpj : cpjets ) { 
//if ( std::abs(cpj.rap()) < rapmax ) {
//    const xAOD::CaloCluster* cpc = getLinkedObject<xAOD::CaloCluster>(cpj); 
//    if ( cpc != nullptr ) { rapWindow.addObject(*cpc,cpj.area()); }
///  }
//  }
  // -- forward towers (non-overlapping with topo-clusters)
//  for ( const auto* cft : *pFwdTCont ) { 
//    if ( cft->e(xAOD::CaloCluster::UNCALIBRATED) > 0. ) { rapWindow.addObject(*cft,ftarea); }

  // ... get profiles etc...
  

