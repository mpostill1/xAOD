#include "xAODAnaHelpers/PhotonContainer.h"

#include <iostream>

using namespace xAH;
using std::vector;
using std::string;

PhotonContainer::PhotonContainer(const std::string& name, const std::string& detailStr, float units, bool mc)
  : ParticleContainer(name, detailStr, units, mc, true)
{


  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    m_ptcone20                    = new std::vector<float> ();
    //m_ptvarcone20                 = new std::vector<float> ();
    m_topoetcone20                = new std::vector<float> ();
    m_isIsolated_Cone20           = new std::vector<int>   ();
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    m_ptcone30                    = new std::vector<float> ();
    //m_ptvarcone30                 = new std::vector<float> ();
    m_topoetcone30                = new std::vector<float> ();
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    m_ptcone40                    = new std::vector<float> ();
    //m_ptvarcone40                 = new std::vector<float> ();
    m_topoetcone40                = new std::vector<float> ();
    m_isIsolated_Cone40CaloOnly   = new std::vector<int>   ();
    m_isIsolated_Cone40           = new std::vector<int>   ();
  }

      // PID
  if(m_infoSwitch.m_PID){
    m_n_IsLoose  = 0;
    m_n_IsMedium = 0;
    m_n_IsTight  = 0;

    m_IsLoose    = new std::vector<int>   ();
    m_IsMedium   = new std::vector<int>   ();
    m_IsTight    = new std::vector<int>   ();
  }

  if(m_infoSwitch.m_purity){
      //Purity
      m_radhad1    = new std::vector<float> ();
      m_radhad     = new std::vector<float> ();
      m_e277	   = new std::vector<float> ();
      m_reta	   = new std::vector<float> ();
      m_rphi	   = new std::vector<float> ();
      m_weta2      = new std::vector<float> ();
      m_f1	   = new std::vector<float> ();
      m_wtot	   = new std::vector<float> ();
      m_deltae     = new std::vector<float> ();
      m_eratio     = new std::vector<float> ();
      //std::vector<float> m_w1
  }

  if(m_infoSwitch.m_effSF && m_mc){
    m_LooseEffSF =new std::vector<float>();
    m_MediumEffSF=new std::vector<float>();
    m_TightEffSF =new std::vector<float>();

    m_LooseEffSF_Error =new std::vector<float>();
    m_MediumEffSF_Error=new std::vector<float>();
    m_TightEffSF_Error =new std::vector<float>();
  }

  if(m_infoSwitch.m_trigger){
    m_trigMatched=new std::vector<std::vector<std::string> >();
  }
}

PhotonContainer::~PhotonContainer()
{

  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    delete m_ptcone20;
    //delete m_ptvarcone20;
    delete m_topoetcone20;
    delete m_isIsolated_Cone20;
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    delete m_ptcone30;
    //delete m_ptvarcone30;
    delete m_topoetcone30;
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    delete m_ptcone40;
    //delete m_ptvarcone40;
    delete m_topoetcone40;
    delete m_isIsolated_Cone40;
    delete m_isIsolated_Cone40CaloOnly;
  }

  // PID
  if(m_infoSwitch.m_PID){
    delete m_IsLoose;
    delete m_IsMedium;
    delete m_IsTight;
  }

  if(m_infoSwitch.m_purity){
    delete m_radhad1;
    delete m_radhad ;
    delete m_e277	;
    delete m_reta	;
    delete m_rphi	;
    delete m_weta2	;
    delete m_f1	;
    delete m_wtot	;
    delete m_deltae;
    delete m_eratio;
    //std::vector<float> m_w1
  }

  if(m_infoSwitch.m_effSF && m_mc){
    delete m_LooseEffSF;
    delete m_MediumEffSF;
    delete m_TightEffSF;

    delete m_LooseEffSF_Error;
    delete m_MediumEffSF_Error;
    delete m_TightEffSF_Error;
  }

  if(m_infoSwitch.m_trigger){
    delete m_trigMatched;
  }
}

void PhotonContainer::setTree(TTree *tree)
{
  //
  // Connect branches
  ParticleContainer::setTree(tree);

  tree->SetBranchStatus  ("nph" , 1);
  tree->SetBranchAddress ("nph" , &m_n);

  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    connectBranch<float>(tree, "ptcone20",                  &m_ptcone20                  );
    //connectBranch<float>(tree, "ptvarcone20",               &m_ptvarcone20               );
    connectBranch<float>(tree, "topoetcone20",              &m_topoetcone20              );
    connectBranch<int>  (tree, "isIsolated_Cone20",         &m_isIsolated_Cone20         );
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    connectBranch<float>(tree, "ptcone30",                  &m_ptcone30                  );
    //connectBranch<float>(tree, "ptvarcone30",               &m_ptvarcone30               );
    connectBranch<float>(tree, "topoetcone30",              &m_topoetcone30              );
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    connectBranch<float>(tree, "ptcone40",                  &m_ptcone40                  );
    //connectBranch<float>(tree, "ptvarcone40",               &m_ptvarcone40               );
    connectBranch<float>(tree, "topoetcone40",              &m_topoetcone40              );
    connectBranch<int>  (tree, "isIsolated_Cone40CaloOnly", &m_isIsolated_Cone40CaloOnly );
    connectBranch<int>  (tree, "isIsolated_Cone40",         &m_isIsolated_Cone40         );
  }

  // PID
  if(m_infoSwitch.m_PID){
    tree->SetBranchStatus (("n"+m_name+"_IsLoose").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsLoose").c_str(),      &m_n_IsLoose);
    connectBranch<int>(tree,  "IsLoose"  , &m_IsLoose );

    tree->SetBranchStatus (("n"+m_name+"_IsMedium").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsMedium").c_str(),      &m_n_IsMedium);
    connectBranch<int>(tree,  "IsMedium" , &m_IsMedium);

    tree->SetBranchStatus (("n"+m_name+"_IsTight").c_str(),     1);
    tree->SetBranchAddress(("n"+m_name+"_IsTight").c_str(),      &m_n_IsTight);
    connectBranch<int>(tree,  "IsTight"  , &m_IsTight );
  }


  if(m_infoSwitch.m_purity){
    connectBranch<float>(tree,"radhad1", &m_radhad1);
    connectBranch<float>(tree,"radhad" , &m_radhad );
    connectBranch<float>(tree,"e277"   , &m_e277   );
    connectBranch<float>(tree,"reta"   , &m_reta   );
    connectBranch<float>(tree,"rphi"   , &m_rphi   );
    connectBranch<float>(tree,"weta2"  , &m_weta2  );
    connectBranch<float>(tree,"f1"     , &m_f1     );
    connectBranch<float>(tree,"wtot"   , &m_wtot   );
    connectBranch<float>(tree,"deltae" , &m_deltae );
    connectBranch<float>(tree,"eratio" , &m_eratio );
  }

  if(m_infoSwitch.m_effSF && m_mc)
    {
      connectBranch<float>(tree, "LooseEffSF", &m_LooseEffSF);
      connectBranch<float>(tree, "MediumEffSF",&m_MediumEffSF);
      connectBranch<float>(tree, "TightEffSF", &m_TightEffSF);

      connectBranch<float>(tree, "LooseEffSF_Error", &m_LooseEffSF_Error);
      connectBranch<float>(tree, "MediumEffSF_Error",&m_MediumEffSF_Error);
      connectBranch<float>(tree, "TightEffSF_Error", &m_TightEffSF_Error);
    }

  if(m_infoSwitch.m_trigger)
    {
      connectBranch<std::vector<std::string> >(tree, "trigMatched", &m_trigMatched);
    }

}

void PhotonContainer::updateParticle(uint idx, Photon& photon)
{
  ParticleContainer::updateParticle(idx,photon);

  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    photon.ptcone20 =                   m_ptcone20                  ->at(idx);
    photon.topoetcone20 =               m_topoetcone20              ->at(idx);
    //photon.ptvarcone20 =                m_ptvarcone20               ->at(idx);
    photon.isIsolated_Cone20 =          m_isIsolated_Cone20         ->at(idx);
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    photon.ptcone30 =                   m_ptcone30                  ->at(idx);
    photon.topoetcone30 =               m_topoetcone30              ->at(idx);
    //photon.ptvarcone30 =                m_ptvarcone30               ->at(idx);
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    photon.ptcone40 =                   m_ptcone40                  ->at(idx);
    //photon.ptvarcone40 =                m_ptvarcone40               ->at(idx);
    photon.topoetcone40 =               m_topoetcone40              ->at(idx);
    photon.isIsolated_Cone40CaloOnly =  m_isIsolated_Cone40CaloOnly ->at(idx);
    photon.isIsolated_Cone40 =          m_isIsolated_Cone40         ->at(idx);
  }

  // PID
  if(m_infoSwitch.m_PID){
    photon.IsLoose =   m_IsLoose ->at(idx);
    photon.IsMedium =  m_IsMedium->at(idx);
    photon.IsTight =   m_IsTight ->at(idx);
  }


  if(m_infoSwitch.m_purity){
    photon.radhad1 = m_radhad1->at(idx);
    photon.radhad =  m_radhad ->at(idx);
    photon.e277 =    m_e277   ->at(idx);
    photon.reta =    m_reta   ->at(idx);
    photon.rphi =    m_rphi   ->at(idx);
    photon.weta2 =   m_weta2  ->at(idx);
    photon.f1 =      m_f1     ->at(idx);
    photon.wtot =    m_wtot   ->at(idx);
    photon.deltae =  m_deltae ->at(idx);
    photon.eratio =  m_eratio ->at(idx);
  }

  if(m_infoSwitch.m_effSF && m_mc){
    photon.LooseEffSF =m_LooseEffSF ->at(idx);
    photon.MediumEffSF=m_MediumEffSF->at(idx);
    photon.TightEffSF =m_TightEffSF ->at(idx);

    photon.LooseEffSF_Error =m_LooseEffSF_Error ->at(idx);
    photon.MediumEffSF_Error=m_MediumEffSF_Error->at(idx);
    photon.TightEffSF_Error =m_TightEffSF_Error ->at(idx);
  }

  if(m_infoSwitch.m_trigger){
    photon.trigMatched =m_trigMatched->at(idx);
  }

}


void PhotonContainer::setBranches(TTree *tree)
{
  ParticleContainer::setBranches(tree);


  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    setBranch<float>(tree, "ptcone20",                  m_ptcone20                  );
    //setBranch<float>(tree, "ptvarcone20",               m_ptvarcone20               );
    setBranch<float>(tree, "topoetcone20",              m_topoetcone20              );
    setBranch<int>  (tree, "isIsolated_Cone20",         m_isIsolated_Cone20         );
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    setBranch<float>(tree, "ptcone30",                  m_ptcone30                  );
    //setBranch<float>(tree, "ptvarcone30",               m_ptvarcone30               );
    setBranch<float>(tree, "topoetcone30",              m_topoetcone30              );
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    setBranch<float>(tree, "ptcone40",                  m_ptcone40                  );
    //setBranch<float>(tree, "ptvarcone40",               m_ptvarcone40               );
    setBranch<float>(tree, "topoetcone40",              m_topoetcone40              );
    setBranch<int>  (tree, "isIsolated_Cone40CaloOnly", m_isIsolated_Cone40CaloOnly );
    setBranch<int>  (tree, "isIsolated_Cone40",         m_isIsolated_Cone40         );
  }

  // PID
  if(m_infoSwitch.m_PID){
    tree->Branch(("n"+m_name+"_IsLoose").c_str(),      &m_n_IsLoose);
    setBranch<int>(tree,  "IsLoose"  , m_IsLoose );

    tree->Branch(("n"+m_name+"_IsMedium").c_str(),      &m_n_IsMedium);
    setBranch<int>(tree,  "IsMedium" , m_IsMedium);

    tree->Branch(("n"+m_name+"_IsTight").c_str(),      &m_n_IsTight);
    setBranch<int>(tree,  "IsTight"  , m_IsTight );
  }

  // purity
  if(m_infoSwitch.m_purity){
    setBranch<float>(tree,"radhad1", m_radhad1);
    setBranch<float>(tree,"radhad" , m_radhad );
    setBranch<float>(tree,"e277"   , m_e277   );
    setBranch<float>(tree,"reta"   , m_reta   );
    setBranch<float>(tree,"rphi"   , m_rphi   );
    setBranch<float>(tree,"weta2"  , m_weta2  );
    setBranch<float>(tree,"f1"     , m_f1     );
    setBranch<float>(tree,"wtot"   , m_wtot   );
    setBranch<float>(tree,"deltae" , m_deltae );
    setBranch<float>(tree,"eratio" , m_eratio );
  }

  // effSF
  if(m_infoSwitch.m_effSF && m_mc){
    setBranch<float>(tree, "LooseEffSF" , m_LooseEffSF);
    setBranch<float>(tree, "MediumEffSF", m_MediumEffSF);
    setBranch<float>(tree, "TightEffSF" , m_TightEffSF);

    setBranch<float>(tree, "LooseEffSF_Error" , m_LooseEffSF_Error);
    setBranch<float>(tree, "MediumEffSF_Error", m_MediumEffSF_Error);
    setBranch<float>(tree, "TightEffSF_Error" , m_TightEffSF_Error);
  }

  // trigger
  if(m_infoSwitch.m_trigger){
    setBranch<std::vector<std::string> >(tree, "trigMatched", m_trigMatched);
  }

  return;
}



void PhotonContainer::clear()
{

  ParticleContainer::clear();

  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    m_ptcone20		  -> clear() ;
    //m_ptvarcone20	  -> clear() ;
    m_topoetcone20	  -> clear() ;
    m_isIsolated_Cone20	  -> clear() ;
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    m_ptcone30		  -> clear() ;
    //m_ptvarcone30	  -> clear() ;
    m_topoetcone30	  -> clear() ;
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    m_ptcone40		  -> clear() ;
    //m_ptvarcone40	  -> clear() ;
    m_topoetcone40        -> clear();
    m_isIsolated_Cone40CaloOnly-> clear();
    m_isIsolated_Cone40	  -> clear() ;
  }

  // PID
  if(m_infoSwitch.m_PID){
    m_n_IsLoose = 0;
    m_IsLoose -> clear();

    m_n_IsMedium = 0;
    m_IsMedium-> clear();

    m_n_IsTight = 0;
    m_IsTight -> clear();
  }

  // purity
  if(m_infoSwitch.m_purity){
    m_radhad1-> clear();
    m_radhad -> clear();
    m_e277   -> clear()	;
    m_reta   -> clear()	;
    m_rphi   -> clear()	;
    m_weta2  -> clear()	;
    m_f1     -> clear()	;
    m_wtot   -> clear()	;
    m_deltae -> clear();
    m_eratio -> clear();
    //std::vector<float> m_w1
  }

  // effSF
  if(m_infoSwitch.m_effSF && m_mc){
    m_LooseEffSF ->clear();
    m_MediumEffSF->clear();
    m_TightEffSF ->clear();

    m_LooseEffSF_Error ->clear();
    m_MediumEffSF_Error->clear();
    m_TightEffSF_Error ->clear();
  }

  // trigger
  if(m_infoSwitch.m_trigger){
    m_trigMatched->clear();
  }

}


void PhotonContainer::FillPhoton( const xAOD::Photon* photon ){
  return FillPhoton(static_cast<const xAOD::IParticle*>(photon));
}

void PhotonContainer::FillPhoton( const xAOD::IParticle* particle )
{

  ParticleContainer::FillParticle(particle);

  const xAOD::Photon* photon=dynamic_cast<const xAOD::Photon*>(particle);


  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "20") != m_infoSwitch.m_isoCones.end()){
    m_ptcone20     -> push_back( photon->isolation( xAOD::Iso::ptcone20    ) / m_units  );
    //m_ptvarcone20  -> push_back( photon->isolation( xAOD::Iso::ptvarcone20 ) / m_units  );
    m_topoetcone20 -> push_back( photon->isolation( xAOD::Iso::topoetcone20) / m_units  );
    static SG::AuxElement::Accessor<char> isIsoCone20Acc            ("isIsolated_FixedCutLoose");
    safeFill<char, int, xAOD::Photon>(photon, isIsoCone20Acc, m_isIsolated_Cone20, -1);
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "30") != m_infoSwitch.m_isoCones.end()){
    m_ptcone30     -> push_back( photon->isolation( xAOD::Iso::ptcone30    ) / m_units  );
    //m_ptvarcone30  -> push_back( photon->isolation( xAOD::Iso::ptvarcone30 ) / m_units  );
    m_topoetcone30 -> push_back( photon->isolation( xAOD::Iso::topoetcone30) / m_units  );
  }
  if(m_infoSwitch.m_isolation || std::find(m_infoSwitch.m_isoCones.begin(), m_infoSwitch.m_isoCones.end(), "40") != m_infoSwitch.m_isoCones.end()){
    m_ptcone40     -> push_back( photon->isolation( xAOD::Iso::ptcone40    ) / m_units  );
    //m_ptvarcone40  -> push_back( photon->isolation( xAOD::Iso::ptvarcone40 ) / m_units  );
    m_topoetcone40 -> push_back( photon->isolation( xAOD::Iso::topoetcone40) / m_units  );

    static SG::AuxElement::Accessor<char> isIsoCone40CaloOnlyAcc    ("isIsolated_FixedCutTightCaloOnly");
    safeFill<char, int, xAOD::Photon>(photon, isIsoCone40CaloOnlyAcc, m_isIsolated_Cone40CaloOnly, -1);

    static SG::AuxElement::Accessor<char> isIsoCone40Acc            ("isIsolated_FixedCutTight");
    safeFill<char, int, xAOD::Photon>(photon, isIsoCone40Acc, m_isIsolated_Cone40, -1);
  }

  if ( m_infoSwitch.m_PID ) {

    static SG::AuxElement::Accessor<bool> phLooseAcc  ("PhotonID_Loose");
    safeFill<bool, int, xAOD::Photon>(photon, phLooseAcc, m_IsLoose, -1);

    static SG::AuxElement::Accessor<bool> phMediumAcc ("PhotonID_Medium");
    safeFill<bool, int, xAOD::Photon>(photon, phMediumAcc, m_IsMedium, -1);

    static SG::AuxElement::Accessor<bool> phTightAcc  ("PhotonID_Tight");
    safeFill<bool, int, xAOD::Photon>(photon, phTightAcc, m_IsTight, -1);

  }

  if (m_infoSwitch.m_purity) {
    static SG::AuxElement::Accessor<float> radhad1  ("Rhad1"  );
    static SG::AuxElement::Accessor<float> radhad   ("Rhad"   );
    static SG::AuxElement::Accessor<float> e277     ("e277"   );
    static SG::AuxElement::Accessor<float> reta     ("Reta"   );
    static SG::AuxElement::Accessor<float> rphi     ("Rphi"   );
    static SG::AuxElement::Accessor<float> weta2    ("weta2"  );
    static SG::AuxElement::Accessor<float> f1       ("f1"     );
    static SG::AuxElement::Accessor<float> wtot     ("wtots1" );
    //static SG::AuxElement::Accessor<float> w1       ("w1"     );
    static SG::AuxElement::Accessor<float> deltae   ("DeltaE" );
    static SG::AuxElement::Accessor<float> eratio   ("Eratio" );

    m_radhad1  -> push_back( radhad1(*photon) );
    m_radhad   -> push_back( radhad (*photon) );
    m_e277     -> push_back( e277   (*photon) );
    m_reta     -> push_back( reta   (*photon) );
    m_rphi     -> push_back( rphi   (*photon) );
    m_weta2    -> push_back( weta2  (*photon) );
    m_f1       -> push_back( f1     (*photon) );
    m_wtot     -> push_back( wtot   (*photon) );
    m_deltae   -> push_back( deltae (*photon) );
    m_eratio   -> push_back( eratio (*photon) );
  }

  if (m_infoSwitch.m_effSF && m_mc) {
    static SG::AuxElement::Accessor<float> PhotonID_Tight_EffSF  ("PhotonID_Tight_EffSF"  );
    static SG::AuxElement::Accessor<float> PhotonID_Medium_EffSF ("PhotonID_Medium_EffSF" );
    static SG::AuxElement::Accessor<float> PhotonID_Loose_EffSF  ("PhotonID_Loose_EffSF"  );

    static SG::AuxElement::Accessor<float> PhotonID_Tight_EffSF_Error  ("PhotonID_Tight_EffSF_Error" );
    static SG::AuxElement::Accessor<float> PhotonID_Medium_EffSF_Error ("PhotonID_Medium_EffSF_Error" );
    static SG::AuxElement::Accessor<float> PhotonID_Loose_EffSF_Error  ("PhotonID_Loose_EffSF_Error" );


    m_TightEffSF  ->push_back( PhotonID_Tight_EffSF (*photon) );
    m_MediumEffSF ->push_back( PhotonID_Medium_EffSF(*photon) );
    m_LooseEffSF  ->push_back( PhotonID_Loose_EffSF (*photon) );

    m_TightEffSF_Error  ->push_back( PhotonID_Tight_EffSF_Error (*photon) );
    m_MediumEffSF_Error ->push_back( PhotonID_Medium_EffSF_Error(*photon) );
    m_LooseEffSF_Error  ->push_back( PhotonID_Loose_EffSF_Error (*photon) );
  }

  if (m_infoSwitch.m_trigger) {
    static SG::AuxElement::Accessor< std::vector< std::string> > trigMatched("trigMatched");

    m_trigMatched ->push_back( trigMatched(*photon) );
  }

  return;
}
