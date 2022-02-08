#include "xAODAnaHelpers/ClusterHists.h"
#include <xAODAnaHelpers/ClusterHistsAlgo.h>
#include <math.h>

ANA_MSG_SOURCE(msgClusterHists, "ClusterHists")

ClusterHists :: ClusterHists (std::string name, std::string detailStr) :
  HistogramManager(name, detailStr)
{
}

ClusterHists :: ~ClusterHists () {}

StatusCode ClusterHists::initialize() {
  

  // These plots are always made
  m_ccl_n   = book(m_name, "n", "cluster multiplicity", 80, 0, 800);
  m_ccl_e   = book(m_name, "e", "cluster e [GeV]", 100, -5, 15);
  m_ccl_eta = book(m_name, "eta", "cluster #eta", 80, -4, 4);
  m_ccl_phi = book(m_name, "phi", "cluster #phi", 120, -TMath::Pi(), TMath::Pi());
  m_ccl_rho = book(m_name, "rho", "cluster #rho",100,0,70000);
  
  // 2D plots
  m_ccl_eta_vs_phi = book(m_name, "eta_vs_phi", "cluster #phi", 120, -TMath::Pi(), TMath::Pi(), "cluster #eta", 80, -4, 4);
  m_ccl_e_vs_eta   = book(m_name, "e_vs_eta", "cluster #eta", 80, -5, 5, "cluster e [GeV]", 100, -5, 15);
  m_ccl_e_vs_phi   = book(m_name, "e_vs_phi", "cluster #phi", 120, -TMath::Pi(), TMath::Pi(), "cluster e [GeV]", 100, -5, 15);
  m_ccl_eta_vs_rho = book(m_name, "eta_vs_rho", "cluster #eta", 80, -5, 5, "cluster #rho ", 100, 0,70000);
  // if worker is passed to the class add histograms to the output

  return StatusCode::SUCCESS;
  
}
 
StatusCode ClusterHists::RhoEtaHistFill(float eventweight,std::vector<std::vector<double>> RhoEtaEntries)
  {
    for (std::vector<double>Entry: RhoEtaEntries)
      {
      m_ccl_eta_vs_rho->Fill(Entry[0],Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }
 
StatusCode ClusterHists::RhoHistFill(float eventweight,std::vector<std::vector<double>> RhoEntries)
  {
    for (std::vector<double>Entry: RhoEntries)
      {
      m_ccl_rho->Fill(Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }
 


StatusCode ClusterHists::execute( const xAOD::CaloClusterContainer* ccls, float eventWeight ) {
  using namespace msgClusterHists;
  xAOD::CaloClusterContainer::const_iterator ccl_itr = ccls->begin();
  xAOD::CaloClusterContainer::const_iterator ccl_end = ccls->end();
  for( ; ccl_itr != ccl_end; ++ccl_itr ) {
    ANA_CHECK( this->execute( (*ccl_itr), eventWeight ));
  }

  m_ccl_n -> Fill( ccls->size(), eventWeight );

  return StatusCode::SUCCESS;
}

StatusCode ClusterHists::execute( const xAOD::CaloCluster* ccl, float eventWeight ) {

  //basic
  float cclE   = ccl->e()/1e3;
  float cclEta = ccl->eta();
  float cclPhi = ccl->phi();
  //float cclrho = ccl->rhoEM();

  m_ccl_e          -> Fill( cclE,   eventWeight );
  m_ccl_eta        -> Fill( cclEta, eventWeight );
  m_ccl_phi        -> Fill( cclPhi, eventWeight );
  //m_ccl_rho        -> Fill( cclrho, eventWeight );
  // 2D plots
  m_ccl_eta_vs_phi -> Fill( cclPhi, cclEta,  eventWeight );
  m_ccl_e_vs_eta   -> Fill( cclEta, cclE,    eventWeight );
  m_ccl_e_vs_phi   -> Fill( cclPhi, cclE,    eventWeight );
  //m_ccl_eta_vs_rho -> fill( cclrho, cclEta,  eventWeight );
  return StatusCode::SUCCESS;




}








