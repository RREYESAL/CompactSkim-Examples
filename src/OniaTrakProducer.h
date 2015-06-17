/**
   \file
   Declaration of OniaTrakProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __OniaTrakProducer_h_
#define __OniaTrakProducer_h_

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

class OniaTrakProducer : public edm::EDProducer {

 public:
  explicit OniaTrakProducer(const edm::ParameterSet& ps);
 
 private:

  virtual void produce(edm::Event& event, const edm::EventSetup& esetup);
  
  virtual void endJob();

  edm::InputTag OniaCollection_;
  edm::InputTag TrakCollection_;
  std::vector<double> OniaMassCuts_;
  std::vector<double> OniaTrakMassCuts_;
  bool OnlyBest_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  const pat::CompositeCandidate makeOniaTCandidate(const pat::CompositeCandidate& onia, 
						   const pat::GenericParticle& trak);
  int candidates;
  int nevents;
  int nonia;
  int nreco;
};

#endif // __OniaTrakProducer_h_
