#include "xAODAnaHelpers/ClusterHists.h"
#include <xAODAnaHelpers/ClusterHistsAlgo.h>
#include <math.h>
#include "TGraph.h"
#include "TCanvas.h"
#include <cmath>
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
  m_ccl_rho_original = book(m_name, "rho", "original #rho_{median} (MeV/Area)",100,0,70000);
  m_ccl_nomralised = book(m_name,"rho_normalised","ratio of binned mean rho_{median}", 48, -5,5);
  
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
  m_ccl_eta_vs_rho_combined_containers_binned_mean = book(m_name, "eta_vs_rho_combined_containers_average", "SlidingWindow #eta", 101, -5, 5, "binned mean rho_{median} (MeV/Area) ", 100, 0,70000);
  m_ccl_eta_vs_rho_combined_containers_binned_mean_comparsion = book(m_name, "eta_vs_rho_combined_containers_average_ratio", "SlidingWindow #eta",101, -5, 5, "ratio of binned mean rho_{median}", 100, 0,1.5);
  m_graph_normalisedratio=book("NormalisedRatio","jet #eta","Ratio",101);
  //m_ccl_eta_vs_jet = book((m_name, "eta_vs_njets", "cluster #eta", 80, -5, 5,m_name,"njets",100,0,1000)
  // if worker is passed to the class add histograms to the output

  return StatusCode::SUCCESS;
  
}


StatusCode ClusterHists::RhoEtaHistFill_combined_containers_average(std::vector<std::vector<std::vector<double>>> RhoEtaEntryCombined_containers,std::vector<float> eventweight_total){
  //for (std::vector<std::vector<double>>containers_entry: RhoEtaEntryCombined_containers){
    //float eventweight = containers_entry[0];}
  std::cout << "Average function start " <<std::endl;
  std::vector<std::vector<double>> sortedreturnvector;
  sortedreturnvector.clear();
  sortedreturnvector.shrink_to_fit();
  double rhoajusted;
  double k; 
  std::list<double> etabining;
  std::list<double> etabining_small_bins;
  std::vector<double> AverageEntryvec;
  double rhoaverage;
  double rhosummed;
  double rhomax;
  std::vector<double> MyEntryvectorsorted;
  etabining = {-5,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};//create eta bin
  etabining_small_bins={-5,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
  
  for (double i:etabining_small_bins){//scan across eta bin
	MyEntryvectorsorted.clear();
	MyEntryvectorsorted.shrink_to_fit();
	MyEntryvectorsorted.push_back(i);
	k=0;
	std::cout << "i value "<< i <<std::endl;
	for (std::vector<std::vector<double>> top_Entry: RhoEtaEntryCombined_containers){
	  // if (i<Entry[0]<i+0.5)
	 
	  for (std::vector<double> Entry: top_Entry){
	  if( (i<Entry[0])&&(Entry[0]<i+0.1)){//find which bin eta falls into 
	    std::cout << "eta value "<< Entry[0] <<std::endl;
	    rhoajusted = (Entry[1])*(eventweight_total[k]);
	    MyEntryvectorsorted.push_back(rhoajusted);
	  }
	  }	
	  k=++k;
	     }
	sortedreturnvector.push_back(MyEntryvectorsorted);
      }
      
      
      std::vector<std::vector<double>> myReturnVectoraverage2;
      myReturnVectoraverage2.clear();
      myReturnVectoraverage2.shrink_to_fit();
      for (std::vector<double>Entry: sortedreturnvector){
	rhomax = Entry.size();
	rhosummed =0;
	for (int i = 1; i<rhomax; i++){
	  rhosummed += Entry[i];}
	AverageEntryvec.clear();
	AverageEntryvec.shrink_to_fit();
	rhoaverage = rhosummed/(rhomax-1);
	//std::cout << "eta value "<< Entry[0] <<std::endl;
	//std::cout << "rhoaverage value "<< rhoaverage <<std::endl;	
	std::cout << "eta assigment "<< Entry[0] <<std::endl;
	AverageEntryvec.push_back(Entry[0]);
	std::cout << "rhoaverage "<< rhoaverage <<std::endl;
	AverageEntryvec.push_back(rhoaverage);
	myReturnVectoraverage2.push_back(AverageEntryvec);
	
      }
      for (std::vector<double>Entry: myReturnVectoraverage2)
	{
	  std::cout << "eta fill "<< Entry[0] <<std::endl;
	  std::cout << "rho fill "<< Entry[1] <<std::endl;
	  m_ccl_eta_vs_rho_combined_containers_binned_mean->Fill(Entry[0],Entry[1]);
	}
      
      
      std::cout << "Average function end " <<std::endl;
      return StatusCode::SUCCESS;
}



StatusCode ClusterHists::RhoEtaHistFill_combined_containers_average_comparsion(std::vector<std::vector<std::vector<double>>> RhoEtaEntryCombined_containers,std::vector<float> eventweight_total){
    //for (std::vector<std::vector<double>>containers_entry: RhoEtaEntryCombined_containers){
    //float eventweight = containers_entry[0];}
  std::cout << "Average function start " <<std::endl;
  std::vector<std::vector<double>> sortedreturnvector;
  sortedreturnvector.clear();
  sortedreturnvector.shrink_to_fit();
  double rhoajusted;
  double k; 
  std::list<double> etabining;
  std::list<double> etabining_small_bins;
  std::list<double> etabining_medium_bins;
  std::vector<double> AverageEntryvec;
  double rhoaverage =0;
  double rhosummed;
  double rhomax;
  double counter;
  double rho_summed_foravg;
  double RMSrho_partial;
  double RMSrho;
  std::vector<double> MyEntryvectorsorted;
  etabining = {-5,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};//create eta bin
  etabining_small_bins={-5,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.0,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0};
  etabining_medium_bins={-5,-4.8,-4.6,-4.4,-4.2,-4.0,-3.8,-3.6,-3.4,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6};
  
  for (double i:etabining_medium_bins){//scan across eta bin
	MyEntryvectorsorted.clear();
	MyEntryvectorsorted.shrink_to_fit();
	MyEntryvectorsorted.push_back(i);
	k=0;
	std::cout << "i value "<< i <<std::endl;
	for (std::vector<std::vector<double>> top_Entry: RhoEtaEntryCombined_containers){
	  // if (i<Entry[0]<i+0.5)
	  for (std::vector<double> Entry: top_Entry){
	    if( (Entry[0]<2.0)&&(Entry[0]>-2.0)){
	      rho_summed_foravg +=Entry[1];
	      counter +=1;
	    }
	  if( (i<Entry[0])&&(Entry[0]<i+0.19999)){//find which bin eta falls into 
	    std::cout << "eta value "<< Entry[0] <<std::endl;
	    rhoajusted = (Entry[1]);
	    MyEntryvectorsorted.push_back(rhoajusted);
	  }
	  }	
	  k=++k;
	     }
	sortedreturnvector.push_back(MyEntryvectorsorted);
      }
  
      
      std::vector<std::vector<double>> myReturnVectoraverage2;
      myReturnVectoraverage2.clear();
      myReturnVectoraverage2.shrink_to_fit();
      double rho_central = rho_summed_foravg/counter;
      std::vector<double> normalised_error;
      double N;
      normalised_error.clear();
      normalised_error.shrink_to_fit();
      for (std::vector<double>Entry: sortedreturnvector){
	rhomax = Entry.size();
	rhosummed =0;
	RMSrho_partial =0;
	N=0;
	for (int i = 1; i<rhomax; i++){
	  rhosummed += Entry[i];
	  RMSrho_partial += pow(Entry[i]/rho_central,2);
	  N+=1;
	}
	AverageEntryvec.clear();
	AverageEntryvec.shrink_to_fit();
	rhoaverage = rhosummed/(rhomax-1);
	//std::cout << "eta value "<< Entry[0] <<std::endl;
	//std::cout << "rhoaverage value "<< rhoaverage <<std::endl;	
	std::cout << "eta assigment "<< Entry[0] <<std::endl;
	AverageEntryvec.push_back(Entry[0]);
	std::cout << "rhoaverage_comparison "<< rhoaverage <<std::endl;
	AverageEntryvec.push_back(rhoaverage/rho_central);
	myReturnVectoraverage2.push_back(AverageEntryvec);
	RMSrho=sqrt(RMSrho_partial)/sqrt(N);
	normalised_error.push_back(RMSrho/sqrt(N));
	}
      
       
      auto Arr_size= myReturnVectoraverage2.size(); // cast to int to shut up compiler warns, really doesn't matter.
      //Check if its Double or Float type
      Double_t  xTgraph[Arr_size],yTgraph[Arr_size];
      Int_t i=0;
      for (std::vector<double>Entry: myReturnVectoraverage2)
	{
	  std::cout << "eta fill "<< Entry[0] <<std::endl;
	  std::cout << "rho fill "<< Entry[1] <<std::endl;

	  xTgraph[i] = Entry[0];
	  yTgraph[i] = Entry[1];
    	 
	  m_graph_normalisedratio->SetPoint( Arr_size,xTgraph[i],yTgraph[i]);
	  m_ccl_eta_vs_rho_combined_containers_binned_mean_comparsion->Fill(Entry[0],Entry[1]);
	  //m_ccl_nomralised->Fill(Entry[0],Entry[1]);
	  //m_ccl_nomralised->SetBinError(Entry[0],Entry[1],normalised_error[i]);	
	  i+=1;
	}
      for (unsigned int j=0; j<normalised_error.size(); j++){
      m_ccl_nomralised->SetBinContent(j,myReturnVectoraverage2[j][1]);
      m_ccl_nomralised->SetBinError(j,normalised_error[j]);}
      
      //m_ccl_nomralised->Sumw2();

	
      //TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,2,600,400);
      //m_graph_normalisedratio->SetPointX(Arr_size,xTgraph);
      //m_graph_normalisedratio->SetPointY(Arr_size,yTgraph);
      
      //m_graph_normalisedratio->SetPoint( Arr_size,xTgraph,yTgraph);	
      // auto m_graph_normalisedratio = (TGraph*) new TGraph(Arr_size,xTgraph,yTgraph);
      //m_graph_normalisedratio->SetTitle("Tnomralisedrho");
      //m_graph_normalisedratio->GetXaxis()->SetTitle("jet #eta");
      //m_graph_normalisedratio->GetYaxis()->SetTitle("Normalised binned mean average (median rho");	
      //m_graph_normalisedratio->Fit("pol2","","",3, 5); 
      //m_graph_normalisedratio->Print("V");
      //m_graph_nomralisedratio->SaveAS("Tnormalisedrho.txt")
      //m_graph_normalisedratio->Draw("AC*");
      //m_graph_normalisedratio->Write();
      std::cout << "Average function end " <<std::endl;
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
StatusCode ClusterHists::RhoHistFill_orig(float eventweight,float rho_orig)
  {
    m_ccl_rho_central->Fill(rho_orig); // you need to know what axes are and what these correspond to 
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








