/**
   \file
   Declaration of OniaTrakTrakProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __OniaTrakTrakProducer_h_
#define __OniaTrakTrakProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <TLorentzVector.h>
#include <vector>

/**
   Create a HF candidate by mathing Onia(chi,psi,etc.) and a track (K, pi, etc.)
 */

class OniaTrakTrakProducer : public edm::EDProducer {

 public:
  explicit OniaTrakTrakProducer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);
  
  virtual void endJob();

  edm::InputTag OniaCollection_;
  edm::InputTag TrakCollection_;
  std::vector<double> OniaMassCuts_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> OniaTrakTrakMassCuts_;
  std::vector<double> MassTraks_;
  bool OnlyBest_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeOniaTTCandidate(const pat::CompositeCandidate& onia, 
						    const pat::CompositeCandidate& tt);
  const pat::CompositeCandidate makeTTCandidate(const pat::GenericParticle& trak1,
                                                const pat::GenericParticle& trak2);
  int candidates;
  int nevents;
  int nonia;
  int nreco;
};

#endif // __OniaTrakTrakProducer_h_
