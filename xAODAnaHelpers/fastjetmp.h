// -*- c++ -*-
#ifndef FASTJETMP_H
#define FASTJETMP_H
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "TH1D.h"
#include "TProfile.h"
#include "xAODAnaHelpers/JetContainer.h"
#include "xAODJet/JetContainer.h"
//#include <boost/tuple/tuple.hpp>"
//#include <map>
#include <string>
#include <vector>

#include "fastjet/PseudoJet.hh"

#include "fastjet/PseudoJet.hh"            // FastJet EDM (Sect. 3.1, p. 11ff)
//#include "xAODAnaHelpers/fastjet.h"
#include "fastjet/ClusterSequenceArea.hh"  // clustering jet with area (p. 38ff)
#include "fastjet/AreaDefinition.hh"       // area definition          (p. 38ff)
#include "fastjet/JetDefinition.hh"        // jet definition           (p. 21ff)
#include "xAODAnaHelpers/fastjetmp.h"

template<class OBJTYPE> class ObjectLink : virtual public fastjet::PseudoJet::UserInfoBase {
private:
  const OBJTYPE* _linkedObj = { nullptr };
public:
  ObjectLink() { }
  ObjectLink(const OBJTYPE* pObj) : _linkedObj(pObj) { }
  virtual ~ObjectLink() { };
  const OBJTYPE* operator*()  const { return _linkedObj; }
  const OBJTYPE* operator->() const { return _linkedObj; } 
  const OBJTYPE* cptr()       const { return _linkedObj; }
};

//template<class OBJTYPE> const OBJTYPE* getLinkedObj(const fastjet::PseudoJet& cjet);
/////








// function retrieving linked object (return NULL in case of wrong object type or no valid link)
template<class OBJTYPE> const OBJTYPE* getLinkedObj(const fastjet::PseudoJet& cjet) { 
  // check object type at return
  const OBJTYPE*  myobject = cjet.has_user_info<OBJTYPE>() ? &cjet.user_info<OBJTYPE>() :  nullptr;
// check if nullptr
  std::cout << "Check on if cpc is a null pointer on not number 2 electric bugalo, is this a null pointer? 1 yes, 0 no  "<< (myobject==nullptr) <<std::endl; 
  return myobject;
  //return cjet.has_user_info<OBJTYPE>() ? &cjet.user_info<OBJTYPE>() :  nullptr;
}

//std::vector<fastjet::PseudoJet> cpjets;
//std::vector<fastjet::PseudoJet> getJetsWithAreas(const xAOD::IParticleContainer& Iparcont);





//inline std::vector<fastjet::PseudoJet> getJetsWithAreas(const xAOD::IParticleContainer& Iparcont) {
inline std::vector<fastjet::PseudoJet> getJets(const xAOD::IParticleContainer& Iparcont) {
  // generate PseudoJets on EM scale
  //std::vector<fastjet::PseudoJet>cpjets.reserve(Iparcont.size());
  std::vector<fastjet::PseudoJet> cpjets; cpjets.reserve(Iparcont.size());
  for ( auto cptr : Iparcont) {
    // consider only positive energy clusters
    if ( cptr->e() > 0. ) { 
      // add PseudoJet to list - note that mass = 0 for topo-clusters 
      cpjets.push_back(fastjet::PtYPhiM(cptr->pt(),cptr->eta(),cptr->phi(),0.));
      // store the link to the original signal (topo-cluster) object (the link object is held as shared pointer) 
      cpjets.back().set_user_info(new ObjectLink<xAOD::IParticle>(cptr));
    }
    //const xAOD::IParticle* cpcc = getLinkedObj<xAOD::IParticle>(cpjets);
    //std::cout << "Check on cpc withni getJets "<< (cpcc==nullptr) <<std::endl;
  }
  return cpjets;
    }
inline fastjet::ClusterSequenceArea getJetsWithAreas(std::vector<fastjet::PseudoJet> cpjets) {
  // set up jet definition
  double jradius(0.4);
  fastjet::JetDefinition  jdef(fastjet::kt_algorithm,jradius);
  // area definition (here Voronoi with usual factor 0.9 for the radius of the circle intersecting with Voronoi cells (do not change)
  fastjet::AreaDefinition adef(fastjet::voronoi_area,fastjet::VoronoiAreaSpec(0.9));
  // cluster the jets
  //
  // Note the sequence: you need to first retrieve the jets from the cluster sequence and then tell the cluster sequence to only
  // destroy itself when the last PseudoJet associated with it is going out of scope. This means you can happily cope the jet 
  // container as many times as you like with its content staying connected with a valid cluster sequence! This is important 
  // because only the cluster sequnce stores important jet features, including the constituents adn the areas!
  //fastjet::ClusterSequenceArea cseq(cpjets,jdef,adef); 
  fastjet::ClusterSequenceArea cs(cpjets,jdef,adef);  // get your cluster sequence for a jets with area
  //std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
  //std::cout << "jets="<<jets.constituents().begin(); <<std::endl;
  //std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets(); // all jets without cuts
  //cs.delete_self_when_unused();       // this lets the cluster sequence "live" until the last object referencing it goes away
  //cseq.delete_self_when_unused();  // this tells the cluster algorithm to destroy itself only after the last reference to it is out of scope
  //
  //return jets; 
  return cs;
  //return cs.inclusive_jets();
}






inline std::vector<fastjet::PseudoJet> getConstituentsWithAreas(const xAOD::IParticleContainer& Iparcont,fastjet::ClusterSequenceArea cs, std::vector<fastjet::PseudoJet> jets) { 
  // get jets from all clusters
  //std::vector<fastjet::PseudoJet> jets(getJetsWithAreas(Iparcont));
  //fastjet::ClusterSequenceArea cs(getJetsWithAreas(Iparcont));
  //std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
  // retrieve constituents
  std::vector<fastjet::PseudoJet> cpjets; cpjets.reserve(Iparcont.size());
  //for ( const auto& jj : jets ) { cpjets.insert(cpjets.end(),jj.constituents().begin(),jj.constituents().end()); }
  //
  for ( const auto& jj : jets ) {std::vector<fastjet::PseudoJet> constits = jj.constituents();
    for (auto c : constits) cpjets.push_back(c);}
  return cpjets; 
}












#endif 
