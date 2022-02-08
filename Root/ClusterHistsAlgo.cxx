#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODJet/JetContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include <xAODAnaHelpers/ClusterHistsAlgo.h>
#include <xAODAnaHelpers/SlidingWindow.h>
#include <xAODAnaHelpers/ClusterHists.h>
#include <math.h>
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
 
  const xAOD::CaloClusterContainer* ccls(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(ccls, m_inContainerName, m_event, m_store, msg()) );
  ANA_CHECK( m_plots->execute( ccls, eventWeight ));
  const xAOD::JetContainer* AntiKt4EMPFlowJets(nullptr);
  ANA_CHECK(HelperFunctions::retrieve(AntiKt4EMPFlowJets,"AntiKt4EMPFlowJets",m_event,m_store, msg() ) );
  std::vector<std::vector<double>> vec_test=SlidingWindow::SlidingWindowExecute(AntiKt4EMPFlowJets, ccls);
  //ANA_CHECK( m_plots->execute(SlidingWindow::myReturnVector, eventWeight));
  //return EL::StatusCode::SUCCESS;
  ANA_CHECK(m_plots->RhoEtaHistFill(eventWeight, vec_test) ) ;
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
