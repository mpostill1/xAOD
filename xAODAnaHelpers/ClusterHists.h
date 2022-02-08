#ifndef xAODAnaHelpers_ClusterHists_H
#define xAODAnaHelpers_ClusterHists_H

#include "xAODAnaHelpers/HistogramManager.h"
#include "xAODCaloEvent/CaloClusterContainer.h"

//std::vector<std::vector<double>> RhoEtaEntries;
//std::vector<double> entry;
//std::vector<std::vector<double>> vec_test;

ANA_MSG_HEADER(msgClusterHists)

class ClusterHists : public HistogramManager
{
  public:
    ClusterHists(std::string name, std::string detailStr );
    ~ClusterHists();

    StatusCode initialize();
    StatusCode execute( const xAOD::CaloClusterContainer* ccls, float eventWeight );
    StatusCode execute( const xAOD::CaloCluster* ccl, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload
    StatusCode RhoEtaHistFill(float eventweight, std::vector<std::vector<double>> RhoEtaEntries );
    StatusCode RhoHistFill(float eventweight, std::vector<std::vector<double>> RhoEntries );
  protected:
    // bools to control which histograms are filled
    bool m_fillDebugging;        //!

  private:
    // Histograms
    TH1F* m_ccl_n; //!
    TH1F* m_ccl_e; //!
    TH1F* m_ccl_eta; //!
    TH1F* m_ccl_phi; //!
    TH1F* m_ccl_rho; //!
    TH2F* m_ccl_eta_vs_phi; //!
    TH2F* m_ccl_e_vs_eta; //!
    TH2F* m_ccl_e_vs_phi; //!
    TH2F* m_ccl_eta_vs_rho; //!
    
};


#endif
