/*
   Package:    PsiTrakRootupler
   Class:      PsiTrakRootupler
 
   Description: make rootuple of J/psi Track combination

   Original Author:  Alberto Sanchez Hernandez
   Created:  based on Alessandro Degano RootupleChic

*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

//
// class declaration
//

class PsiTrakRootupler : public edm::EDAnalyzer {
   public:
      explicit PsiTrakRootupler(const edm::ParameterSet&);
      ~PsiTrakRootupler();

      bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  std::string file_name;
  edm::InputTag oniat_cand_Label;
  edm::InputTag oniat_rf_cand_Label;
  edm::InputTag primaryVertices_Label;
  edm::InputTag triggerResults_Label;
  int  oniat_pdgid_, onia_pdgid_, trak_pdgid_;
  bool isMC_,OnlyBest_;

  UInt_t run,event,numPrimaryVertices, trigger;  

  TLorentzVector oniat_p4;
  TLorentzVector psi_p4;
  TLorentzVector track_p4;
  TLorentzVector muonp_p4;
  TLorentzVector muonn_p4;

  TLorentzVector oniat_rf_p4;
  TLorentzVector psi_rf_p4;
  TLorentzVector track_rf_p4;

  Int_t    oniat_charge, psi_triggerMatch, track_nvsh, track_nvph;
  Double_t oniat_vProb,  oniat_vChi2, oniat_cosAlpha, oniat_ctauPV, oniat_ctauErrPV;
  Double_t track_d0, track_d0Err, track_dz, track_dxy;
  Double_t psi_vProb, psi_vChi2, psi_DCA, psi_ctauPV, psi_ctauErrPV, psi_cosAlpha, psi_nSigma;

  TLorentzVector gen_oniat_p4;
  Int_t          gen_oniat_pdgId;
  TLorentzVector gen_psi_p4;
  Int_t          gen_onia_pdgId;
  TLorentzVector gen_track_p4;
  TLorentzVector gen_muonp_p4;
  TLorentzVector gen_muonn_p4;

  TTree* oniat_tree;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
static const double pi0_mass    =  0.1349766;
static const Double_t psi1SMass =  3.09691;
static const Double_t psi2SMass =  3.68610;
static const Double_t ups1SMass =  9.46030;

// 2011 par
//static const double Y_sig_par_A = 0.058;
//static const double Y_sig_par_B = 0.047;
//static const double Y_sig_par_C = 0.22;

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

//
// constructors and destructor
//
PsiTrakRootupler::PsiTrakRootupler(const edm::ParameterSet& iConfig):
	oniat_cand_Label(iConfig.getParameter<edm::InputTag>("oniat_cand")),
        oniat_rf_cand_Label(iConfig.getParameter<edm::InputTag>("oniat_rf_cand")),
	primaryVertices_Label(iConfig.getParameter<edm::InputTag>("primaryVertices")),
	triggerResults_Label(iConfig.getParameter<edm::InputTag>("TriggerResults")),
        oniat_pdgid_(iConfig.getParameter<uint32_t>("oniat_pdgid")),
        onia_pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
        trak_pdgid_(iConfig.getParameter<uint32_t>("trak_pdgid")),
	isMC_(iConfig.getParameter<bool>("isMC")),
        OnlyBest_(iConfig.getParameter<bool>("OnlyBest"))
{
	edm::Service<TFileService> fs;
        oniat_tree = fs->make<TTree>("OniaTrakTree","Tree of Onia and Track");

        oniat_tree->Branch("run",                &run,                "run/I");
        oniat_tree->Branch("event",              &event,              "event/I");
        oniat_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
        oniat_tree->Branch("trigger",            &trigger,            "trigger/I");

        oniat_tree->Branch("oniat_p4",   "TLorentzVector", &oniat_p4);
        oniat_tree->Branch("track_p4",   "TLorentzVector", &track_p4);
        oniat_tree->Branch("psi_p4",     "TLorentzVector", &psi_p4);
        oniat_tree->Branch("muonp_p4",   "TLorentzVector", &muonp_p4);
        oniat_tree->Branch("muonn_p4",   "TLorentzVector", &muonn_p4);

        oniat_tree->Branch("oniat_rf_p4", "TLorentzVector", &oniat_rf_p4);
        oniat_tree->Branch("psi_rf_p4",   "TLorentzVector", &psi_rf_p4);
        oniat_tree->Branch("track_rf_p4", "TLorentzVector", &track_rf_p4);

        oniat_tree->Branch("psi_vProb",        &psi_vProb,        "psi_vProb/D");
        oniat_tree->Branch("psi_vNChi2",       &psi_vChi2,        "psi_vNChi2/D");
        oniat_tree->Branch("psi_DCA",          &psi_DCA,          "psi_DCA/D");
        oniat_tree->Branch("psi_ctauPV",       &psi_ctauPV,       "psi_ctauPV/D");
        oniat_tree->Branch("psi_ctauErrPV",    &psi_ctauErrPV,    "psi_ctauErrPV/D");
        oniat_tree->Branch("psi_cosAlpha",     &psi_cosAlpha,     "psi_cosAlpha/D");
        oniat_tree->Branch("psi_nSigma",       &psi_nSigma,       "psi_nSigma/D");
        oniat_tree->Branch("psi_triggerMatch", &psi_triggerMatch, "psi_triggerMatch/I");

        oniat_tree->Branch("oniat_vProb",      &oniat_vProb,        "oniat_vProb/D");
        oniat_tree->Branch("oniat_vChi2",      &oniat_vChi2,        "oniat_vChi2/D");
        oniat_tree->Branch("oniat_cosAlpha",   &oniat_cosAlpha,     "oniat_cosAlpha/D");
        oniat_tree->Branch("oniat_ctauPV",     &oniat_ctauPV,       "oniat_ctauPV/D");
        oniat_tree->Branch("oniat_ctauErrPV",  &oniat_ctauErrPV,    "oniat_ctauErrPV/D");
        oniat_tree->Branch("oniat_charge",     &oniat_charge,       "oniat_charge/I");

        oniat_tree->Branch("track_d0",    &track_d0,    "track_d0/D");
        oniat_tree->Branch("track_d0Err", &track_d0Err, "track_d0Err/D");
        oniat_tree->Branch("track_dz",    &track_dz,    "track_dz/D");
        oniat_tree->Branch("track_dxy",   &track_dxy,   "track_dxy/D");
        oniat_tree->Branch("track_nvsh",  &track_nvsh,  "track_nvsh/I");
        oniat_tree->Branch("track_nvph",  &track_nvph,  "track_nvph/I");

	if(isMC_)
	  {
            oniat_tree->Branch("gen_oniat_pdgId", &gen_oniat_pdgId, "gen_oniat_pdgId/I");
            oniat_tree->Branch("gen_onia_pdgId",  &gen_onia_pdgId,  "gen_onia_pdgId/I");
	    oniat_tree->Branch("gen_oniat_p4",    "TLorentzVector", &gen_oniat_p4);
	    oniat_tree->Branch("gen_psi_p4",      "TLorentzVector", &gen_psi_p4);
	    oniat_tree->Branch("gen_track_p4",    "TLorentzVector", &gen_track_p4);
            oniat_tree->Branch("gen_muonp_p4",    "TLorentzVector", &gen_muonp_p4);
            oniat_tree->Branch("gen_muonn_p4",    "TLorentzVector", &gen_muonn_p4);
	  }
}

PsiTrakRootupler::~PsiTrakRootupler() {}

//
// member functions
//

bool PsiTrakRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


// ------------ method called for each event  ------------
void PsiTrakRootupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
//  using namespace edm;
  using namespace std;

  edm::Handle<std::vector<pat::CompositeCandidate>> oniat_cand_handle;
  iEvent.getByLabel(oniat_cand_Label, oniat_cand_handle);

  edm::Handle<std::vector<pat::CompositeCandidate>> oniat_rf_cand_handle;
  iEvent.getByLabel(oniat_rf_cand_Label, oniat_rf_cand_handle);

  edm::Handle<std::vector<reco::Vertex >> primaryVertices_handle;
  iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByLabel( triggerResults_Label , triggerResults_handle);

  numPrimaryVertices = primaryVertices_handle->size();
  run = iEvent.id().run();
  event = iEvent.id().event();

  edm::Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByLabel("prunedGenParticles",pruned);
  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByLabel("packedGenParticles",packed);

  if (isMC_ && packed.isValid() && pruned.isValid()) {
     gen_oniat_pdgId  = 0;
     gen_onia_pdgId = 0;
     gen_oniat_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     int foundit   = 0;
     for (size_t i=0; i<pruned->size(); i++) {
         int p_id = (*pruned)[i].pdgId();
         const reco::Candidate *aonia = &(*pruned)[i];
         if ( (abs(p_id) ==  oniat_pdgid_) && aonia->status() == 2) {
            gen_oniat_p4.SetPtEtaPhiM(aonia->pt(),aonia->eta(),aonia->phi(),aonia->mass());
            gen_oniat_pdgId = p_id;
            foundit++;
            for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
               const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
               const reco::Candidate * d = &(*packed)[j];
               if ( motherInPrunedCollection != nullptr && (abs(d->pdgId()) == trak_pdgid_)  && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_track_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonn_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonp_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( foundit == 4 ) break;
            }
            if ( foundit == 4 ) {
               gen_psi_p4 = gen_muonn_p4 + gen_muonp_p4;   // this should take into account FSR
               gen_onia_pdgId = onia_pdgid_;
               break;
            } else {
               foundit = 0;
               gen_oniat_pdgId = 0;
               gen_onia_pdgId  = 0;
               gen_oniat_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            }
         }  // if ( p_id
     } // for (size

//... sanity check
     if (!gen_oniat_pdgId) std::cout << "PsiTrakRootupler: didn't find the given decay " << run << "," << event << std::endl;
  } // end if isMC

	// grab Trigger information
	// save it in variable trigger, trigger is an int between 0 and 7, in binary it is:
	// (pass 10)(pass 8)(pass 0)
	// ex. 7 = pass 0, 8 and 10
	// ex. 6 = pass 8, 10
        // ex. 1 = pass 0

  trigger = 0;
  if ( triggerResults_handle.isValid() ) {
    const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	    
    vector <unsigned int> bits_0;
    for ( int version = 1; version<2; version ++ ) {
      stringstream ss;
      ss<<"HLT_Dimuon13_PsiPrime_v"<<version;
      bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
    }  

    vector <unsigned int> bits_1;
    for ( int version = 1; version<2; version ++ ) {
      stringstream ss;
      ss<<"HLT_Dimuon13_Upsilon_v"<<version;
      bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
    }

    vector <unsigned int> bits_2;
    for ( int version = 1; version<2; version ++ ) {
      stringstream ss;
      ss<<"HLT_Dimuon20_Jpsi_v"<<version;
      bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss.str()).label().c_str()));
    }

    for (unsigned int i=0; i<bits_0.size(); i++) {
      unsigned int bit = bits_0[i];
      if ( bit < triggerResults_handle->size() ){
	if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
          trigger += 1;
          break;
        }
      }
    }

    for (unsigned int i=0; i<bits_1.size(); i++) {
      unsigned int bit = bits_1[i];
      if ( bit < triggerResults_handle->size() ){
	if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
          trigger += 2;
          break;
        }
      }
    }

    for (unsigned int i=0; i<bits_2.size(); i++) {
      unsigned int bit = bits_2[i];
      if ( bit < triggerResults_handle->size() ){
	if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
          trigger += 4;
          break;
        }
      }
    }
  } else {
    std::cout<< "No valid trigger information found " << run << "," << event <<std::endl;
  }
	
// grabbing oniat information
  if (!oniat_cand_handle.isValid()) std::cout<< "No oniat information " << run << "," << event <<std::endl;
  if (!oniat_rf_cand_handle.isValid()) std::cout<< "No oniat_rf information " << run << "," << event <<std::endl;
// get rf information. Notice we are just keeping combinations with succesfull vertex fit
  if (oniat_rf_cand_handle.isValid() && oniat_cand_handle.isValid()) {
    pat::CompositeCandidate oniat_rf_cand, oniat_cand, *onia_cand;
    for (unsigned int i=0; i< oniat_rf_cand_handle->size(); i++){
      oniat_rf_cand   = oniat_rf_cand_handle->at(i);
      int   bindx     = oniat_rf_cand.userInt("bIndex");
      oniat_vProb     = oniat_rf_cand.userFloat("vProb");
      oniat_vChi2     = oniat_rf_cand.userFloat("vChi2");
      oniat_cosAlpha  = oniat_rf_cand.userFloat("cosAlpha");
      oniat_ctauPV    = oniat_rf_cand.userFloat("ctauPV");
      oniat_ctauErrPV = oniat_rf_cand.userFloat("ctauErrPV"); 
      oniat_charge    = oniat_rf_cand.charge();
      if (oniat_charge != oniat_rf_cand.daughter("trak")->charge()) 
         std::cout << "Charge mismatch " << run << "," << event << std::endl;
      oniat_rf_p4.SetPtEtaPhiM(oniat_rf_cand.pt(),oniat_rf_cand.eta(),oniat_rf_cand.phi(),oniat_rf_cand.mass());
      psi_rf_p4.SetPtEtaPhiM(oniat_rf_cand.daughter("onia")->pt(),oniat_rf_cand.daughter("onia")->eta(),
                              oniat_rf_cand.daughter("onia")->phi(),oniat_rf_cand.daughter("onia")->mass());
      track_rf_p4.SetPtEtaPhiM(oniat_rf_cand.daughter("trak")->pt(),oniat_rf_cand.daughter("trak")->eta(),
                              oniat_rf_cand.daughter("trak")->phi(),oniat_rf_cand.daughter("trak")->mass());
      if (bindx<0 || bindx>(int) oniat_cand_handle->size()) {
        std::cout << "Incorrect index for oniat combination " << run << "," << event <<"," << bindx << std::endl;
        continue; 
      }
      oniat_cand = oniat_cand_handle->at(bindx);
      onia_cand = dynamic_cast <pat::CompositeCandidate *>(oniat_cand.daughter("onia"));

      track_d0      = oniat_cand.userFloat("trak_d0");
      track_d0Err   = oniat_cand.userFloat("trak_d0err");
      track_dz      = oniat_cand.userFloat("trak_dz");
      track_dxy     = oniat_cand.userFloat("trak_dxy");
      track_nvsh    = oniat_cand.userInt("trak_nvsh");
      track_nvph    = oniat_cand.userInt("trak_nvph");

      psi_vProb        = onia_cand->userFloat("vProb");
      psi_vChi2        = onia_cand->userFloat("vNChi2");
      psi_DCA          = onia_cand->userFloat("DCA");
      psi_ctauPV       = onia_cand->userFloat("ppdlPV");
      psi_ctauErrPV    = onia_cand->userFloat("ppdlErrPV");
      psi_cosAlpha     = onia_cand->userFloat("cosAlpha");
      psi_triggerMatch = 0;

      reco::Candidate::LorentzVector vP = onia_cand->daughter("muon1")->p4();
      reco::Candidate::LorentzVector vM = onia_cand->daughter("muon2")->p4();
      if (onia_cand->daughter("muon1")->charge() < 0) {
         vP = onia_cand->daughter("muon2")->p4();
         vM = onia_cand->daughter("muon1")->p4();
      }

      double kmass = 0.4936770;
      oniat_p4.SetPtEtaPhiM(oniat_cand.pt(),oniat_cand.eta(),oniat_cand.phi(),oniat_cand.mass());
      psi_p4.SetPtEtaPhiM(onia_cand->pt(),onia_cand->eta(),onia_cand->phi(),onia_cand->mass());
      track_p4.SetPtEtaPhiM(oniat_cand.daughter("trak")->pt(),oniat_cand.daughter("trak")->eta(),
                           oniat_cand.daughter("trak")->phi(), kmass);
      muonp_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
      muonn_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

      double sigma = Y_sig_par_A + Y_sig_par_B*pow(fabs(psi_p4.Rapidity()),2) + Y_sig_par_C*pow(fabs(psi_p4.Rapidity()),3);
      sigma /= 1000.;
      sigma *= psi1SMass/ups1SMass;
      psi_nSigma = fabs(psi_p4.M() - psi1SMass) / sigma;

      oniat_tree->Fill();
      if (OnlyBest_) break;  // oniat candidates are sorted by vProb
    }
  }
 
}

// ------------ method called once each job just before starting event loop  ------------
void PsiTrakRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void PsiTrakRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void PsiTrakRootupler::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void PsiTrakRootupler::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void PsiTrakRootupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void PsiTrakRootupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PsiTrakRootupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PsiTrakRootupler);
