
#include "xAODAnaHelpers/SlidingWindow.h"

#include "TMath.h"

//#include "TStopwatch.h"

#include <cstdio>

size_t SlidingWindow::Window::m_invalidIndex = size_t(-1);

SlidingWindow::Window::Window()
  : m_nBins(0)
  , m_min(0.)
  , m_max(0.)
  , m_window(0.)
  , m_binWidth(0.)
  , m_adjustBins(0)
  , m_adjustBinRange(0)
  , m_adjustMin(0.)
  , m_adjustMax(0.)
  , m_adjustBinWidth(0.)
  , m_firstCall(true)
  , m_isDead(true)
  , m_storeObjects(false)
  , h_spreadProfile((TH1D*)0)
  , h_areaProfile((TH1D*)0)
  , h_medianProfile((TH1D*)0)
  , p_eventProfile((TProfile*)0)
{ }

SlidingWindow::Window::Window(const std::string& pName,size_t nBins,double xMin,double xMax,double window,bool storeObjects)
  : m_nBins(nBins)
  , m_min(xMin)
  , m_max(xMax)
  , m_window(window)
  , m_binWidth(0.)
  , m_adjustBins(0)
  , m_adjustBinRange(0)
  , m_adjustMin(0.)
  , m_adjustMax(0.)
  , m_adjustBinWidth(0.)
  , m_firstCall(true)
  , m_isDead(true)
  , m_storeObjects(storeObjects)
  , h_spreadProfile((TH1D*)0)
  , h_areaProfile((TH1D*)0)
  , h_medianProfile((TH1D*)0)
  , p_eventProfile((TProfile*)0)
{
  static double _twoPi = TMath::TwoPi();
  // invalid bin options
  if ( m_nBins == 0 || ( m_min >= m_max ) ) {
    printf("SlidingWindow::Window::Window():%s ERROR Invalid bin specifications #bins/min/max = %i/%4.1f/%4.1f\n",
	   pName.c_str(),(int)m_nBins,m_min,m_max);
    delete this;
  }
  // invalid window size
  if ( m_window <= 0. ) { 
    printf("SlidingWindow::Window::Window():%s ERROR Invalid window size %4.1f\n",pName.c_str(),m_window);
    delete this;
  }

  // setup
  m_binWidth = (m_max-m_min)/((double)m_nBins); 
  // book event profile
  m_adjustBins = m_nBins % 2 == 0 ? m_nBins + 1 : m_nBins; ///////////////////////
  double fx(m_binWidth/2.);                                // Re-adjusting bins //
  double ctr((m_min+m_max)/2.);                            ///////////////////////  
  size_t ir((m_adjustBins-1)/2); 
  m_adjustMin      = (ctr-fx)-ir*m_binWidth;
  m_adjustMax      = (ctr+fx)+ir*m_binWidth;
  m_adjustBinWidth = (m_adjustMax-m_adjustMin)/((double)m_adjustBins); 
  p_eventProfile = new TProfile(pName.c_str(),pName.c_str(),m_adjustBins,m_adjustMin,m_adjustMax);
  if ( p_eventProfile == 0 ) return;
  // store area profile
  std::string hName = "WindowArea_"+pName;
  h_areaProfile = new TH1D(hName.c_str(),hName.c_str(),m_adjustBins,m_adjustMin,m_adjustMax);
  // bin spread
  hName = "WindowSpread_"+pName;
  h_spreadProfile = new TH1D(hName.c_str(),hName.c_str(),m_adjustBins,m_adjustMin,m_adjustMax);
  // median profile
  hName = "Median_"+pName;
  h_medianProfile = new TH1D(hName.c_str(),hName.c_str(),m_adjustBins,m_adjustMin,m_adjustMax);

  // adjusted parameters and descriptors
  if ( h_areaProfile == 0 || h_spreadProfile == 0 || h_medianProfile == 0 ) return;

  // calculate areas and binning
  double xm(m_adjustMin+m_adjustBinWidth/2.);
  for ( size_t i(0); i<m_adjustBins; ++i )
    {
      double mn = TMath::Max(xm-m_window/2.,m_adjustMin);
      double mx = TMath::Min(xm+m_window/2.,m_adjustMax);
      double ar((mx-mn)*_twoPi);
      m_bins.push_back(Bin(xm,ar,mn,mx));
      xm += m_adjustBinWidth;
    }

  // range of window in bin index space (one-sided)
  m_adjustBinRange = (size_t)((m_window/2.)/m_adjustBinWidth);

  // plot areas and set up collector
  for ( auto bin : m_bins ) { h_areaProfile->Fill(bin.binCenter(),bin.binArea()); }
  this->setupCollector();
  // set up internal store
  m_isDead = false;

  printf("SlidingWindow::Window::Window():%s INFO Specifications\n",pName.c_str());
  printf("SlidingWindow::Window::Window():%s INFO Number of bins .......... %i\n",pName.c_str(),(int)m_adjustBins-1);
  printf("SlidingWindow::Window::Window():%s INFO Range ................... [%5.2f,%5.2f]\n",pName.c_str(),m_adjustMin,m_adjustMax);
  printf("SlidingWindow::Window::Window():%s INFO Bin width ............... %5.2f\n",pName.c_str(),m_adjustBinWidth);
  printf("SlidingWindow::Window::Window():%s INFO Window size ............. %5.2f\n",pName.c_str(),m_window);
  printf("SlidingWindow::Window::Window():%s INFO Window index range ...... +/-%i\n",pName.c_str(),(int)m_adjustBinRange);
  int ii(0);
  int mb((int)m_adjustBinRange);
  //  int ms((int)m_adjustBins-1);
  for ( auto bin : m_bins ) { 
    double loRange(TMath::Max(bin.binCenter()-mb*m_adjustBinWidth,m_adjustMin));
    double hiRange(TMath::Min(bin.binCenter()+mb*m_adjustBinWidth,m_adjustMax));
    printf("SlidingWindow::Window::Window():%s INFO Bin %3i left/center/right %5.2f/%5.2f/%5.2f range [%5.2f,%5.2f] (area) %5.2f\n",
	   pName.c_str(),ii++,bin.binLower(),bin.binCenter(),bin.binUpper(),loRange,hiRange,bin.binArea());
  }
}

SlidingWindow::Window::~Window()
{ }

void SlidingWindow::Window::reset()
{
  if ( m_isDead ) return;
  //
  if ( m_firstCall ) { m_firstCall = false; return; }

  for ( size_t i(0); i<m_collector.size(); ++i  ) { 
    for ( auto rhomed : m_collector.at(i) ) { p_eventProfile->Fill(m_bins.at(i).binCenter(),rhomed); }
  }

  // 
  h_medianProfile->Reset();

  // reset event summation
  this->resetCollector();
}

bool SlidingWindow::Window::integrate(const fastjet::PseudoJet& pMom)
{
  // retrieve area
  double area = this->isComplete(pMom) ? pMom.area() : 0.; 
  // 
  return this->integrate(pMom,area);
}

bool SlidingWindow::Window::integrate(const fastjet::PseudoJet& pMom,double area)
{
  // exclude outside boundaries 
  double eta(pMom.pseudorapidity());
  if ( eta < m_min || eta > m_max ) return false;
  // find index range
  double rhoLocal = area > 0. ? pMom.pt()/area : 0.;
  int idx(this->binIndex(eta));
  //  printf("SidingWindow::Window::integrate() INFO eta value/index %5.2f/%3i\n",eta,idx);
  int ilo(TMath::Max(idx-(int)m_adjustBinRange,0));
  int ihi(TMath::Min(idx+m_adjustBinRange,m_adjustBins-1));
  // store data for central bin only
  m_areas[idx].push_back(area); 
  if ( m_storeObjects ) { m_objects[idx].push_back(pMom); }
  // store data in windows
  for ( int i(ilo); i<=ihi; ++i ) { m_collector[i].push_back(rhoLocal); }
  return true;
}

void SlidingWindow::Window::setupCollector()
{
  static size_t _reserve = 10000;
  m_collector.resize(m_bins.size(),data_t());
  m_areas.resize(m_bins.size(),data_t());
  if ( !m_storeObjects ) { 
    for ( size_t i(0); i<m_bins.size(); ++i ) { m_collector.at(i).reserve(_reserve); m_areas.at(i).reserve(_reserve); }
    printf("SlidingWindow::setupCollector [%s] - collector/area cache size %i/%i, reserved %i\n",p_eventProfile->GetName(),(int)m_collector.size(),(int)m_areas.size(),(int)_reserve);
  } else {
    m_objects.resize(m_bins.size(),object_t());
    for ( size_t i(0); i<m_bins.size(); ++i ) { m_collector.at(i).reserve(_reserve); m_areas.at(i).reserve(_reserve); m_objects.at(i).reserve(_reserve); }
    printf("SlidingWindow::setupCollector [%s] - collector/area/object cache size %i/%i/%i, reserved %i\n",p_eventProfile->GetName(),
	   (int)m_collector.size(),(int)m_areas.size(),(int)m_objects.size(),(int)_reserve);
  }
}

void SlidingWindow::Window::resetCollector()
{ 
  if ( m_storeObjects ) { 
    for ( size_t i(0); i<m_bins.size(); ++i ) { m_collector[i].clear(); m_areas[i].clear(); m_objects[i].clear(); } 
  } else {
    for ( size_t i(0); i<m_bins.size(); ++i ) { m_collector[i].clear(); m_areas[i].clear(); } 
  }
}

size_t SlidingWindow::Window::binIndex(double value) const
{ return ( value < m_adjustMin || value >= m_adjustMax ) ? m_invalidIndex : (size_t)((value-m_adjustMin)/m_adjustBinWidth); }

bool SlidingWindow::Window::findAllBins(double value,std::vector<size_t>& indices) const
{
  indices.clear();
  int idx(this->binIndex(value));
  if ( idx < 0 || idx >= (int)m_adjustBins ) { return false; }
  int ilo(TMath::Max(idx-(int)m_adjustBinRange,0));
  int ihi(TMath::Min(idx+m_adjustBinRange,m_adjustBins-1));
  for ( int i(ilo); i<=ihi; ++i ) { indices.push_back(i); } 
  return true; 
}

std::vector<double> SlidingWindow::Window::findMedians()
{
  data_t medians; medians.resize(m_collector.size(),0.);
  if ( !m_storeObjects ) { 
    for ( size_t i(0); i<m_collector.size(); ++i ) { medians[i] = this->median(m_collector.at(i)); }
  } else {
    for ( size_t i(0); i<m_objects.size(); ++i ) { medians[i] = this->median(m_objects.at(i)); }
  }
  // fill distribution if needed
  if ( h_medianProfile->GetEntries() == 0 ) { 
    for ( auto bin : m_bins ) { 
      h_medianProfile->Fill(bin.binCenter(),medians.at(this->binIndex(bin.binCenter())));
    } 
  }

  return medians;
}

double SlidingWindow::Window::median(data_t& values)                 //////////////////////
{std::cout << "Check on call "<< values.size()<<std::endl;                                                                    // Sorts input data //
  if ( !values.empty() ) {                                           //////////////////////
    if ( values.size() == 1 ) { return values.front(); }
    std::sort(values.begin(),values.end());
    size_t n(values.size());
    return n % 2 != 0 ? values.at((n-1)/2) : (values.at(n/2-1)+values.at(n/2))/2.;
  } else { return 0.; }
}

double SlidingWindow::Window::median(object_t& objects)              //////////////////////   
{                                                                     // Sorts input data //
  static _sortObjects _sorter;                                       ////////////////////// 
  std::cout << "Check on rho "<< rhoFromJet(objects.at(objects.size())) <<std::endl;
  if ( !objects.empty() ) {                                          
    if ( objects.size() == 1 ) { return this->rhoFromJet(objects.front()); }
    std::sort(objects.begin(),objects.end(),_sorter);
    if ( this->rhoFromJet(objects.front()) >= this->rhoFromJet(objects.back()) ) { printf("SlidingWindow::Window::median() WARNING BUG wrong sort order"); }
    size_t n(objects.size());
    if ( n % 2 != 0 ) { 
      return this->rhoFromJet(objects.at((n-1)/2)); 
      //std::cout << "Check on rho "<< rhoFromJet(objects.at(n)) <<std::endl;
    } else {
      return (this->rhoFromJet(objects.at((n-1)/2))+this->rhoFromJet(objects.at(n/2)))/2.;
    }
  } else { return 0.; }
}


//std::vector<double> SlidingWindow::Window::findrho()
//{
//data_t medians; medians.resize(m_collector.size(),0.);
// if ( !m_storeObjects ) { 
//  for ( size_t i(0); i<m_collector.size(); ++i ) { medians[i] = this->workoutrho(m_collector.at(i)); }
// } else {
//   for ( size_t i(0); i<m_objects.size(); ++i ) { medians[i] = this->workoutrho(m_objects.at(i)); }
//}
  // fill distribution if needed
// if ( h_medianProfile->GetEntries() == 0 ) { 
//   for ( auto bin : m_bins ) { 
//    h_medianProfile->Fill(bin.binCenter(),workoutrho.at(this->binIndex(bin.binCenter())));
//   } 
// }

// return medians;
//}



//double SlidingWindow::Window::workoutrho(object_t& objects)              //////////////////////   
//{                                                                    // Sorts input data //
// static _sortObjects _sorter;                                       //////////////////////                                                                   
// if ( !objects.empty() ) {                                          
//  if ( objects.size() == 1 ) { return this->rhoFromJet(objects.front()); }
//  std::sort(objects.begin(),objects.end(),_sorter);
//  if ( this->rhoFromJet(objects.front()) >= this->rhoFromJet(objects.back()) ) { printf("SlidingWindow::Window::median() WARNING BUG wrong sort order"); }
//  size_t n(objects.size());
///  if ( n % 2 != 0 ) { 
//   return this->rhoFromJet(objects); 
//  } else {
//    return (this->rhoFromJet(objects);
//  }
//} else { return 0.; }
//}
