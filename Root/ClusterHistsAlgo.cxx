#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODJet/JetContainer.h"
#include "xAODPFlow/FlowElementContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include <xAODAnaHelpers/ClusterHistsAlgo.h>
#include <xAODAnaHelpers/SlidingWindow.h>
#include <xAODAnaHelpers/ClusterHists.h>
#include <math.h>
#include <xAODJet/JetContainer.h>
#include <xAODCore/AuxContainerBase.h>
#include <iostream>
// this is needed to distribute the algorithm to the workers
ClassImp(ClusterHistsAlgo)

ClusterHistsAlgo :: ClusterHistsAlgo () :
    Algorithm("ClusterHistsAlgo")
{
}

EL::StatusCode ClusterHistsAlgo :: setupJob (EL::Job& job)
{
  job.useXAOD();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("ClusterHistsAlgo").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ClusterHistsAlgo :: histInitialize ()
{

  ANA_MSG_INFO( m_name );
  ANA_CHECK( xAH::Algorithm::algInitialize());
  // needed here and not in initalize since this is called first
  if( m_inContainerName.empty() || m_detailStr.empty() ){
    ANA_MSG_ERROR( "One or more required configuration values are empty");
    return EL::StatusCode::FAILURE;
  }

  // declare class and add histograms to output
  m_plots = new ClusterHists(m_name, m_detailStr);
  //rho_plots = new; 
  ANA_CHECK( m_plots -> initialize());
  //rho_plots->m_ccl_eta_vs_rho();
  m_plots -> record( wk() );
// this is needed to distribute the algorithm to the workers
  //ClassImp(ClusterHistsAlgo);

  //ClusterHistsAlgo :: ClusterHistsAlgo () :
  // Algorithm("ClusterHistsAlgo")
  //{
  //}
  return EL::StatusCode::SUCCESS;
  //EL::StatusCode ClusterHistsAlgo :: setupJob (EL::Job& job)
  //{
  //job.useXAOD();

  // let's initialize the algorithm to use the xAODRootAccess package
  //xAOD::Init("ClusterHistsAlgo").ignore(); // call before opening first file

  //return EL::StatusCode::SUCCESS;
}

//EL::StatusCode ClusterHistsAlgo :: histInitialize ()
//{

//ANA_MSG_INFO( m_name );
//ANA_CHECK( xAH::Algorithm::algInitialize());
 // needed here and not in initalize since this is called first
//if( m_inContainerName.empty() || m_detailStr.empty() ){
//ANA_MSG_ERROR( "One or more required configuration values are empty");
//return EL::StatusCode::FAILURE;
// }


  // declare class and add histograms to output
  //m_plots = new ClusterHists(m_name, m_detailStr);
  //rho_plots = new; 
// ANA_CHECK( m_plots -> initialize());
  //rho_plots->m_ccl_eta_vs_rho();
// m_plots -> record( wk() );
  
//return EL::StatusCode::SUCCESS;


EL::StatusCode ClusterHistsAlgo :: fileExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode ClusterHistsAlgo :: changeInput (bool /*firstFile*/) { return EL::StatusCode::SUCCESS; }

EL::StatusCode ClusterHistsAlgo :: initialize ()
{
  ANA_MSG_INFO( "ClusterHistsAlgo");
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ClusterHistsAlgo :: execute ()
{


  
  
  const xAOD::EventInfo* eventInfo(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store, msg()) );


  float eventWeight(1);
  if( eventInfo->isAvailable< float >( "mcEventWeight" ) ) {
    eventWeight = eventInfo->auxdecor< float >( "mcEventWeight" );
  }
  //define new container that combines both central and forward topo towers. 
  //xAOD::IParticleContainer Combined_Container(SG::VIEW_ELEMENTS);
  const xAOD::CaloClusterContainer* ccls(nullptr);
  const xAOD::CaloClusterContainer* Fwd(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(ccls, m_inContainerName, m_event, m_store, msg()) );
  ANA_CHECK( m_plots->execute( ccls, eventWeight ));
  ANA_CHECK( HelperFunctions::retrieve(Fwd,  m_inContainerName, m_event, m_store, msg()) );
  ANA_CHECK( m_plots->execute( Fwd, eventWeight ));  

  






   // Retrieve input jet container
  const xAOD::FlowElementContainer* CFE_container=nullptr;
  const xAOD::FlowElementContainer* NFE_container=nullptr;
  
  ANA_CHECK( HelperFunctions::retrieve(CFE_container,  "JetETMissChargedParticleFlowObjects", m_event, m_store, msg()) );
  ANA_CHECK( HelperFunctions::retrieve(NFE_container,  "JetETMissNeutralParticleFlowObjects", m_event, m_store, msg()) );

  // Create the new container and its auxiliary store.
  //const xAOD::IParticleContainer* central_combined=nullptr;
  auto central_combined = std::make_unique<xAOD::IParticleContainer>();
  auto central_combinedAux = std::make_unique<xAOD::AuxContainerBase>();
  //const xAOD::AuxContainerBase central_combinedAux=nullptr;
  //std::unique_ptr<const xAOD::AuxContainerBase> central_combinedAux= std::make_unique<xAOD::AuxContainerBase>();
    central_combined->setStore (central_combinedAux.get()); //< Connect the two
  for (const xAOD::FlowElement* CFE_Element : *CFE_container) {
    // Copy this CFE to the output container:
    xAOD::FlowElement* CFE = new xAOD::FlowElement();
    central_combined->push_back(CFE); // CFE acquires the central_combined auxstore
    *CFE=*CFE_Element; // copies auxdata from one auxstore to the other
  }
  
  for (const xAOD::FlowElement* NFE_Element : *NFE_container) {
    // Copy this NFE to the output container:
    xAOD::FlowElement* NFE= new xAOD::FlowElement();
    central_combined->push_back(NFE); // jet acquires the goodJets auxstore
    *NFE = *NFE_Element; // copies auxdata from one auxstore to the other
  }
  
  
  std::cout << "works up to this point";

  ANA_CHECK( m_plots->execute( central_combined.get(), eventWeight ));
  const xAOD::JetContainer* AntiKt4EMPFlowJets(nullptr);
  ANA_CHECK(HelperFunctions::retrieve(AntiKt4EMPFlowJets,"AntiKt4EMPFlowJets",m_event,m_store, msg() ) );
  //ANA_CHECK( m_plots->execute( AntiKt4EMPFlowJets, eventWeight ));
  
 
  //central_combined=CFE_container+NFE_container
  //Combined_Container= central_combined+Fwd
  std::vector<std::vector<double>> vec_central=SlidingWindow::SlidingWindowExecute(AntiKt4EMPFlowJets, central_combined.get());
  std::vector<std::vector<double>> vec_forward=SlidingWindow::SlidingWindowExecutetopo(AntiKt4EMPFlowJets, Fwd);
  //std::vector<std::vector<double>> vec_forward=SlidingWindow::SlidingWindowExecutetCombined(AntiKt4EMPFlowJets, Combined_Container);
  //std::vector<std::vector<double>> vec_combined=SlidingWindow::SlidingWindowExecuteFill(AntiKt4EMPFlowJets, CMB);
  //ANA_CHECK( m_plots->execute(SlidingWindow::myReturnVector, eventWeight));
  //return EL::StatusCode::SUCCESS;
  ANA_CHECK(m_plots->RhoHistFill_central(eventWeight, vec_central) ) ;
  ANA_CHECK(m_plots->RhoEtaHistFill_central(eventWeight, vec_central) ) ;
  ANA_CHECK(m_plots->RhoHistFill_fwd(eventWeight, vec_forward) ) ;
  ANA_CHECK(m_plots->RhoEtaHistFill_fwd(eventWeight, vec_forward) ) ;
  //ANA_CHECK(m_plots->RhoHistFill(eventWeight, vec_combined) ) ;
  //ANA_CHECK(m_plots->RhoEtaHistFill(eventWeight, vec_combined) ) ;

  // Record the objects into the event store
  ANA_CHECK (evtStore()->record (central_combined.release(), "central_combined"));
  ANA_CHECK (evtStore()->record (central_combinedAux.release(), "central_combinedAux"));


  return EL::StatusCode::SUCCESS;
}


EL::StatusCode ClusterHistsAlgo :: postExecute () { return EL::StatusCode::SUCCESS; }
EL::StatusCode ClusterHistsAlgo :: finalize () { return EL::StatusCode::SUCCESS; }
EL::StatusCode ClusterHistsAlgo :: histFinalize ()
{
  // clean up memory
  if(m_plots) delete m_plots;
  ANA_CHECK( xAH::Algorithm::algFinalize());
  return EL::StatusCode::SUCCESS;

}
