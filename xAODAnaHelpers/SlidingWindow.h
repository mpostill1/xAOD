// -*- c++ -*-
#ifndef SLIDINGWINDOW_H
#define SLIDINGWINDOW_H
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "TH1D.h"
#include "TProfile.h"
#include "xAODAnaHelpers/JetContainer.h"
#include "xAODJet/JetContainer.h"
//#include <boost/tuple/tuple.hpp>"
//#include <map>
#include <string>
#include <vector>

#include "fastjet/PseudoJet.hh"            // FastJet EDM (Sect. 3.1, p. 11ff)
//#include "xAODAnaHelpers/fastjet.h"
#include "fastjet/ClusterSequenceArea.hh"  // clustering jet with area (p. 38ff)
#include "fastjet/AreaDefinition.hh"       // area definition          (p. 38ff)
#include "fastjet/JetDefinition.hh"        // jet definition           (p. 21ff)
#include "xAODAnaHelpers/fastjetmp.h"

namespace SlidingWindow
{
  std::vector<std::vector<double>> SlidingWindowExecute(const xAOD::JetContainer*  jetcont,const xAOD::IParticleContainer*  clustercont);
  std::vector<std::vector<double>> SlidingWindowExecutetopo(const xAOD::JetContainer*  jetcont,const xAOD::CaloClusterContainer*  clustercont);
  //std::vector<std::vector<double>> myReturnVector;
  fastjet::ClusterSequenceArea cs(std::vector<fastjet::PseudoJet> cpjets,fastjet::JetDefinition jdef,fastjet::AreaDefinition adef);
  //std::vector<fastjet::PseudoJet> cpjet;
  class Bin
  {
  public:
  Bin() : m_left(0.),m_right(0.),m_center(0.),m_area(0.) { }
  Bin(double center,double area,double left,double right) : m_left(left),m_right(right),m_center(center),m_area(area) { }
    virtual ~Bin() { }

    virtual double binCenter() const { return m_center; }
    virtual double binArea()   const { return m_area;   }
    virtual double binLower()  const { return m_left;   }
    virtual double binUpper()  const { return m_right;  }

    virtual bool   inBin(double value) const 
    { return value >= binLower() && value < binUpper(); }

    bool operator<(const Bin& bin) const
    { return this->binCenter() < bin.binCenter(); }

  private:

    double m_left;
    double m_right;
    double m_center;
    double m_area;
  };

  class Window
  {
  public:

    typedef std::vector<double>             data_t;
    typedef std::vector<data_t>             collector_t;
    typedef std::vector<fastjet::PseudoJet> object_t;
    typedef std::vector<object_t>           object_collector_t;

    typedef Bin                    bin_t;
    typedef std::vector<bin_t>     binning_t;

    Window();
    Window(const std::string& pName,size_t nBins,double min,double max,double window,bool storeObjects=false);
    virtual ~Window();

    bool add(const object_t& mom)
    {
      if ( m_isDead ) return false;
      // find central bin and distribute
      bool isOk(true);
      for ( auto pm : mom ) { isOk = isOk && this->integrate(pm); }
      return isOk;
    }

    template<class OBJ>
    bool addObject(const OBJ& obj,double area=0.)
    {
      if ( m_isDead ) return false;
      // find central bin and distribute
      fastjet::PseudoJet fj(fastjet::PtYPhiM(obj.pt(),obj.phi(),obj.m()));
      return integrate(fj,area);
    }//obj.p4().Rap()

    virtual void   reset();

    virtual size_t binIndex(double value) const;
    virtual int    rootBinIndex(double value) const;
    virtual int    rootBinIndex(size_t index) const;

    bool findAllBins(double value,std::vector<size_t>& indices) const;

    const binning_t& binning() const { return m_bins; }

    void  Sliding_window_rho(xAOD::JetContainer* jetcont, xAOD::CaloClusterContainer* clustercont);
    double windowSize() const { return m_window; }
    double min() const { return m_min; }
    double max() const { return m_max; }
    double adjustedMin() const { return m_adjustMin; }
    double adjustedMax() const { return m_adjustMax; }

    TH1D*     areaProfile()   { return h_areaProfile; }
    TProfile* eventProfile()  { return p_eventProfile; }
    TH1D*     spreadProfile() { return h_spreadProfile; }

    std::vector<double> findMedians(); 

    const collector_t&        rhoCollection()    const { return m_collector; }
    const collector_t&        areaCollection()   const { return m_areas;     }
    const object_collector_t& objectCollection() const { return m_objects;   }

  protected:
    // stores
    binning_t          m_bins;
    collector_t        m_collector;
    collector_t        m_areas;
    object_collector_t m_objects;

    // requested binning and window size
    size_t m_nBins;
    double m_min;
    double m_max;
    double m_window;
    double m_binWidth;
    // adjusted binning after centering on mean value in range
    size_t m_adjustBins;
    size_t m_adjustBinRange;
    double m_adjustMin;
    double m_adjustMax;
    double m_adjustBinWidth;

    // process steering
    bool m_firstCall;
    bool m_isDead;
    bool m_storeObjects;
    static size_t m_invalidIndex;


    TH1D*     h_spreadProfile;
    TH1D*     h_areaProfile;
    TH1D*     h_medianProfile;
    TProfile* p_eventProfile;

    // helpers
    double area(int rootIdx) const;
    double center(int rootIdx) const;
    double lowEdge(int rootIdx) const;
    double upEdge(int rootIdx) const;
    const Bin& bin(int rootIdx) const;

    bool   integrate(const fastjet::PseudoJet& pMom);
    bool   integrate(const fastjet::PseudoJet& pMom,double area);
    void   setupCollector();
    void   resetCollector();
    double median(data_t& values);
    double median(object_t& jet);
    bool   isComplete(const fastjet::PseudoJet& jet) const;
    double rhoFromJet(const fastjet::PseudoJet& jet) const;

    struct _sortObjects 
    { 
      bool operator()(const fastjet::PseudoJet& p0,const fastjet::PseudoJet& p1) 
      { 
	if ( ( p0.has_valid_cs() && p0.area() > 0. ) && ( p1.has_valid_cs() && p1.area() > 0. ) ) { 
	  return p0.perp()/p0.area() < p1.perp()/p1.area(); 
	} else {
	  return false;
	}
      }
    };
  };
}

inline double SlidingWindow::Window::area(int rootIdx) const
{ return rootIdx > 0 && rootIdx <= (int)m_bins.size() 
    ? m_bins.at(rootIdx-1).binArea() : 0.; }

inline double SlidingWindow::Window::center(int rootIdx) const
{ return rootIdx > 0 && rootIdx <= (int)m_bins.size() 
    ? m_bins.at(rootIdx-1).binCenter() : 0.; }

inline double SlidingWindow::Window::lowEdge(int rootIdx) const
{ return rootIdx > 0 && rootIdx <= (int)m_bins.size() 
    ? m_bins.at(rootIdx-1).binLower() : 0.; }

inline double SlidingWindow::Window::upEdge(int rootIdx) const
{ return rootIdx > 0 && rootIdx <= (int)m_bins.size() 
    ? m_bins.at(rootIdx-1).binUpper() : 0.; }

inline const SlidingWindow::Bin& SlidingWindow::Window::bin(int rootIdx) const
{ return rootIdx > 0 && rootIdx <= (int)m_bins.size() 
    ? m_bins.at(rootIdx-1) : *m_bins.end(); }

inline int SlidingWindow::Window::rootBinIndex(double value) const
{ return this->rootBinIndex(this->binIndex(value)); }

inline int SlidingWindow::Window::rootBinIndex(size_t index) const
{ return index != m_invalidIndex ? index+1 : -1; }


inline bool SlidingWindow::Window::isComplete(const fastjet::PseudoJet& jet) const
{ return jet.has_valid_cs() && jet.area() > 0.; } 
inline double SlidingWindow::Window::rhoFromJet(const fastjet::PseudoJet& jet) const
{std::cout << "Check on jet.perp() "<< jet.perp() <<std::endl; 
  std::cout << "Check on jet.area() "<< jet.area() <<std::endl; 
return this->isComplete(jet) ? jet.perp()/jet.area() : 0.;} 



#endif
