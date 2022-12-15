#ifndef xAODAnaHelpers_ClusterHists_H
#define xAODAnaHelpers_ClusterHists_H
#include "TGraph.h"
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
    //StatusCode execute( const xAOD::CaloClusterContainer* ccls, float eventWeight );
    //StatusCode execute( const xAOD::CaloCluster* ccl, float eventWeight );
    StatusCode execute(const xAOD::IParticleContainer* Ipars, float eventWeight );
    StatusCode execute(const xAOD::IParticle* Ipar, float eventWeight );
    using HistogramManager::book; // make other overloaded versions of book() to show up in subclass
    using HistogramManager::execute; // overload
    StatusCode RhoEtaHistFill_central(float eventweight, std::vector<std::vector<double>> RhoEtaEntries );
    StatusCode RhoHistFill_orig(float eventweight,float rho_orig);
    StatusCode RhoHistFill_central(float eventweight, std::vector<std::vector<double>> RhoEntries );
    StatusCode RhoEtaHistFill_fwd(float eventweight, std::vector<std::vector<double>> RhoEtaEntries );
    StatusCode RhoHistFill_fwd(float eventweight, std::vector<std::vector<double>> RhoEntries );
    StatusCode RhoEtaHistFill_combined(float eventweight,std::vector<std::vector<double>> RhoEtaEntrycentral,std::vector<std::vector<double>> RhoEtaEntryfwd);
    StatusCode RhoEtaHistFill_combined_containers(float eventweight,std::vector<std::vector<double>> RhoEtaEntryCombined_containers);
    StatusCode RhoEtaHistFill_combined_containers_average(std::vector<std::vector<std::vector<double>>> RhoEtaEntryCombined_containers, std::vector<float>  eventweight_total);
    StatusCode RhoEtaHistFill_combined_containers_average_comparsion(std::vector<std::vector<std::vector<double>>> RhoEtaEntryCombined_containers, std::vector<float>  eventweight_total);
  protected:
    // bools to control which histograms are filled
    bool m_fillDebugging;        //!

  private:
    // Histograms
    TH1F* m_ccl_n; //!
    TH1F* m_ccl_e; //!
    TH1F* m_ccl_eta; //!
    TH1F* m_ccl_phi; //!
    TH1F* m_ccl_pt; //!
    TH1F* m_ccl_rho_central; //!
    TH1F* m_ccl_rho_fwd; //!
    TH1F* m_ccl_rho_original; //!
    TH1F* m_ccl_rhoEM; //!
    TH1F* m_ccl_nomralised; //!
    TH2F* m_ccl_eta_vs_phi; //!
    TH2F* m_ccl_e_vs_eta; //!
    TH2F* m_ccl_e_vs_phi; //!
    TH2F* m_ccl_eta_vs_rho_central; //!
    TH2F* m_ccl_eta_vs_rho_fwd; //!
    TH2F* m_ccl_eta_vs_pt; //!
    TH2F* m_ccl_eta_vs_rho_combined; //!
    TH2F* m_ccl_eta_vs_rho_combined_containers; //!
    TH2F* m_ccl_eta_vs_rho_combined_containers_binned_mean; //!
    TH2F* m_ccl_eta_vs_rho_combined_containers_binned_mean_comparsion; //!
    TGraph*  m_graph_normalisedratio; //!

    
};


#endif
