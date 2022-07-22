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
  m_ccl_n   = book(m_name, "n", "tower multiplicity", 80, 0, 800);
  m_ccl_e   = book(m_name, "e", "tower e [GeV]", 100, -5, 15);
  m_ccl_eta = book(m_name, "eta", "tower/cluster #eta", 80, -5, 5);
  m_ccl_eta->SetYTitle ("Number of towers/clusters");
  m_ccl_eta->SetStats(0);
  m_ccl_phi = book(m_name, "phi", "tower #phi", 120, -TMath::Pi(), TMath::Pi());
  m_ccl_pt = book(m_name,"pt", "tower/cluster pt",100,0,1000);
  // m_ccl_jet = book(m_name,"njets","njets",100,0,1000)
  //m_ccl_rhoEM = book(m_name, "rhoEM", "cluster #rho (MeV/Area)",100,0,70000); 
  m_ccl_rho_central = book(m_name, "rho", "tower #rho_{median} (MeV/Area)",100,0,70000);
  m_ccl_rho_fwd = book(m_name, "rho", "tower #rho_{median} (MeV/Area)",100,0,70000);
  
  // 2D plots
  m_ccl_eta_vs_phi = book(m_name, "eta_vs_phi", "tower #phi", 120, -TMath::Pi(), TMath::Pi(), "cluster #eta", 80, -4, 4);
  m_ccl_e_vs_eta   = book(m_name, "e_vs_eta", "tower #eta", 80, -5, 5, "tower e [GeV]", 100, -5, 15);
  m_ccl_e_vs_phi   = book(m_name, "e_vs_phi", "tower #phi", 120, -TMath::Pi(), TMath::Pi(), "cluster e [GeV]", 100, -5, 15);
  m_ccl_eta_vs_pt  = book(m_name, "eta_vs_pt", "tower/cluster #eta", 80, -5, 5, "tower p_{t} (MeV) ", 100, 0,1000);
  //m_ccl_eta_vs_pt->GetZaxis()->SetRangeUser(1,2000);
  m_ccl_eta_vs_rho_central = book(m_name, "eta_vs_rho_central", "SlidingWindow #eta", 80, -5, 5, "#rho_{median} (MeV/Area) ", 50, 0,70000);
  m_ccl_eta_vs_rho_fwd = book(m_name, "eta_vs_rho_fwd", "SlidingWindow #eta", 80, -5, 5, "#rho_{median} (MeV/Area) ", 50, 0,70000);
  m_ccl_eta_vs_rho_combined = book(m_name, "eta_vs_rho_combined", "SlidingWindow #eta", 80, -5, 5, "#rho_{median} (MeV/Area) ", 50, 0,70000);
  m_ccl_eta_vs_rho_combined_containers = book(m_name, "eta_vs_rho_combined_containers", "SlidingWindow #eta", 80, -5, 5, "#rho_{median} (MeV/Area) ", 50, 0,70000);
  m_ccl_eta_vs_rho_combined_containers_binned_mean = book(m_name, "eta_vs_rho_combined_containers", "SlidingWindow #eta", 10, -5, 5, "#binned mean rho_{median} (MeV/Area) ", 50, 0,70000);
  
  //m_ccl_eta_vs_jet = book((m_name, "eta_vs_njets", "cluster #eta", 80, -5, 5,m_name,"njets",100,0,1000)
  // if worker is passed to the class add histograms to the output

  return StatusCode::SUCCESS;
  
}




StatusCode ClusterHists::RhoEtaHistFill_combined_containers(float eventweight,std::vector<std::vector<double>> RhoEtaEntryCombined_containers)
{//etabining=[-5,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5];//create eta bin
    for (std::vector<double>Entry: RhoEtaEntryCombined_containers)
      {
	m_ccl_eta_vs_rho_combined_containers->Fill(Entry[0],Entry[1]);



	 //for (etabinning){//scan across eta bin
	   //if (etabining[i]<Entry[0]<etabining[i+1]){//find which bin eta falls into 
	     //rho[i].add?=Entry[1];//store rho value 
	     //}
  
	   //}
	 
	 //} 

    //for 

      //need mean of rho here
      //m_ccl_eta_vs_rho_combined_containers_binned_mean->Fill(eta,mean);//fill mean values
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
      }
    return StatusCode::SUCCESS;
}




StatusCode ClusterHists::RhoEtaHistFill_combined(float eventweight,std::vector<std::vector<double>> RhoEtaEntrycentral,std::vector<std::vector<double>> RhoEtaEntryfwd)
  {
    for (std::vector<double>Entry: RhoEtaEntrycentral)
      {
	if (std::abs(Entry[0])<2.5) {m_ccl_eta_vs_rho_combined->Fill(Entry[0],Entry[1]);} 
      } 
      for (std::vector<double>Entry: RhoEtaEntryfwd)
      {
	if (std::abs(Entry[0])>2.5) {m_ccl_eta_vs_rho_combined->Fill(Entry[0],Entry[1]);} 
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }


StatusCode ClusterHists::RhoEtaHistFill_central(float eventweight,std::vector<std::vector<double>> RhoEtaEntries)
  {
    for (std::vector<double>Entry: RhoEtaEntries)
      {
      m_ccl_eta_vs_rho_central->Fill(Entry[0],Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }

StatusCode ClusterHists::RhoEtaHistFill_fwd(float eventweight,std::vector<std::vector<double>> RhoEtaEntries)
  {
    for (std::vector<double>Entry: RhoEtaEntries)
      {
      m_ccl_eta_vs_rho_fwd->Fill(Entry[0],Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }
 
StatusCode ClusterHists::RhoHistFill_central(float eventweight,std::vector<std::vector<double>> RhoEntries)
  {
    for (std::vector<double>Entry: RhoEntries)
      {
      m_ccl_rho_central->Fill(Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }
 
StatusCode ClusterHists::RhoHistFill_fwd(float eventweight,std::vector<std::vector<double>> RhoEntries)
  {
    for (std::vector<double>Entry: RhoEntries)
      {
      m_ccl_rho_central->Fill(Entry[1]); // you need to know what axes are and what these correspond to
      } 
    //since we have to return a StatusCode, return the obvious one (SUCCESS)
    return StatusCode::SUCCESS;
  }

//StatusCode ClusterHists::execute( const xAOD::CaloClusterContainer* ccls, float eventWeight ) {
StatusCode ClusterHists::execute(const xAOD::IParticleContainer* Ipars, float eventWeight ) {
  using namespace msgClusterHists;
  xAOD::IParticleContainer::const_iterator ccl_itr = Ipars->begin();
  xAOD::IParticleContainer::const_iterator ccl_end = Ipars->end();
  for( ; ccl_itr != ccl_end; ++ccl_itr ) {
    ANA_CHECK( this->execute( (*ccl_itr), eventWeight ));
  }

  m_ccl_n -> Fill( Ipars->size(), eventWeight );

  return StatusCode::SUCCESS;
}

//StatusCode ClusterHists::execute( const xAOD::CaloCluster* ccl, float eventWeight ) {
StatusCode ClusterHists::execute(const xAOD::IParticle* Ipar, float eventWeight ) {

  //basic
  float cclE   = Ipar->e()/1e3;
  float cclEta = Ipar->eta();
  float cclPhi = Ipar->phi();
  float cclpt  = Ipar->pt();
  //float cclnjets = ccl->njets();
  //float cclrhoEM = ccl->rhoEM();

  m_ccl_e          -> Fill( cclE,   eventWeight );
  m_ccl_eta        -> Fill( cclEta, eventWeight );
  m_ccl_phi        -> Fill( cclPhi, eventWeight );
  m_ccl_pt         -> Fill( cclpt, eventWeight  );
  // m_ccl_rhoEM        -> Fill( cclrhoEM, eventWeight );
  // 2D plots
  m_ccl_eta_vs_phi -> Fill( cclPhi, cclEta,  eventWeight );
  m_ccl_e_vs_eta   -> Fill( cclEta, cclE,    eventWeight );
  m_ccl_e_vs_phi   -> Fill( cclPhi, cclE,    eventWeight );
  m_ccl_eta_vs_pt  -> Fill( cclEta, cclpt,   eventWeight );
  //m_ccl_eta_vs_rho -> Fill( cclEta, m_ccl_rho,  eventWeight );
  return StatusCode::SUCCESS;




}








