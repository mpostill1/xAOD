#ifndef xAODAnaHelpers_ClusterHistsAlgo_H
#define xAODAnaHelpers_ClusterHistsAlgo_H

#include <xAODAnaHelpers/ClusterHists.h>

// algorithm wrapper
#include "xAODAnaHelpers/Algorithm.h"
//std::vector<std::vector<std::vector<double>>>  vec_combined_total;
class ClusterHistsAlgo : public xAH::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  std::string m_inContainerName = "";
  // configuration variables
  std::string m_detailStr = "";
  //xAOD::IParticleContainer central_combined(nullptr); 
  //std::unique_ptr<DataVector<xAOD::IParticle> central_combined;
  std::vector<std::vector<std::vector<double>>>  vec_combined_total;
  std::vector<std::vector<std::vector<double>>>  vec_normalised_total;
  std::vector<float>  vec_eventweight_total;
private:
  ClusterHists* m_plots = nullptr; //!
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  // this is a standard constructor
  ClusterHistsAlgo ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  /// @cond
  // this is needed to distribute the algorithm to the workers
  ClassDef(ClusterHistsAlgo, 1);
  /// @endcond
  
};

#endif
