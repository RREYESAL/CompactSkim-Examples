#include <CompactSkim/Examples/src/OniaTrakProducer.h>

OniaTrakProducer::OniaTrakProducer(const edm::ParameterSet& ps):
  OniaCollection_(ps.getParameter<edm::InputTag>("Onia")),
  TrakCollection_(ps.getParameter<edm::InputTag>("Trak")),
  OniaMassCuts_(ps.getParameter<std::vector<double>>("OniaMassCuts")),
  OniaTrakMassCuts_(ps.getParameter<std::vector<double>>("OniaTrakMassCuts")),
  OnlyBest_(ps.getParameter<bool>("OnlyBest"))
{
  produces<pat::CompositeCandidateCollection>("OniaTrakCandidates");
  candidates = 0;
  nevents = 0;
  nonia = 0;
  nreco = 0;
}
 
void OniaTrakProducer::produce(edm::Event& event, const edm::EventSetup& esetup){

  std::auto_ptr<pat::CompositeCandidateCollection> OniaTCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> onia;
  event.getByLabel(OniaCollection_,onia);

  edm::Handle<std::vector<pat::GenericParticle> > trak;
  event.getByLabel(TrakCollection_,trak);

  uint ncombo = 0;
  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];
  float OniaTrakMassMax_ = OniaTrakMassCuts_[1];
  float OniaTrakMassMin_ = OniaTrakMassCuts_[0];

  int ionia = -1;
// Note: Dimuon cand are sorted by decreasing vertex probability then first chi is associated with "best" dimuon 
  for (pat::CompositeCandidateCollection::const_iterator oniaCand = onia->begin(); oniaCand != onia->end(); ++oniaCand){
     ionia++;
// If J/psi use reco mass otherwise use mQ
     float oniaM = oniaCand->mass();
     if ( oniaM < OniaMassMax_  && oniaM > OniaMassMin_ ) {

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon1"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(oniaCand->daughter("muon2"));
       
// loop on track candidates, make OniaT candidate
       for (std::vector<pat::GenericParticle>::const_iterator trakCand = trak->begin(); trakCand!= trak->end(); ++trakCand){
         if ( IsTheSame(*trakCand,*pmu1) || IsTheSame(*trakCand,*pmu2) ) continue;
	 pat::CompositeCandidate OniaTCand = makeOniaTCandidate(*oniaCand, *trakCand);
         if ( OniaTCand.mass() < OniaTrakMassMax_ && OniaTCand.mass() > OniaTrakMassMin_) {
           OniaTCand.addUserInt("oIndex",ionia);
           float d0 = trakCand->track()->d0();
           float d0err = trakCand->track()->d0Error();
           float dz = trakCand->track()->dz(oniaCand->vertex());
           float dxy = trakCand->track()->dxy(oniaCand->vertex());
           int   nvstriphits = trakCand->track()->hitPattern().numberOfValidStripHits();
           int   nvpixelhits = trakCand->track()->hitPattern().numberOfValidPixelHits();
           OniaTCand.addUserFloat("trak_d0",d0);
           OniaTCand.addUserFloat("trak_d0err",d0err);
           OniaTCand.addUserFloat("trak_dz",dz);
           OniaTCand.addUserFloat("trak_dxy",dxy);
           OniaTCand.addUserInt("trak_nvsh",nvstriphits);
           OniaTCand.addUserInt("trak_nvph",nvpixelhits);
	   OniaTCandColl->push_back(OniaTCand);
	   candidates++;
           ncombo++;
         }
       }

     }
     if (OnlyBest_) break;
  }
  if ( ncombo != OniaTCandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != OniaT ("<<OniaTCandColl->size()<<")"<< std::endl;
  if ( onia->size() > 0 )  nonia++;
  if ( ncombo > 0 ) nreco++; 
  event.put(OniaTCandColl,"OniaTrakCandidates"); 
  nevents++;
}

void OniaTrakProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "OniaTrak Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with Onia candidates " << nonia << std::endl;
  std::cout << "Events with OniaTrak candidates " << nreco << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " OniaTrak candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

bool OniaTrakProducer::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
 
const pat::CompositeCandidate 
OniaTrakProducer::makeOniaTCandidate(const pat::CompositeCandidate& onia, 
				  const pat::GenericParticle& trak){

  pat::CompositeCandidate OniaTCand;
  OniaTCand.addDaughter(onia,"onia");
  OniaTCand.addDaughter(trak,"trak");
  OniaTCand.setVertex(onia.vertex());
  OniaTCand.setCharge(trak.charge());

  double m_kaon = 0.4936770;
  math::XYZVector mom_kaon = trak.momentum();
  double e_kaon = sqrt(m_kaon*m_kaon + mom_kaon.Mag2());
  math::XYZTLorentzVector p4_kaon = math::XYZTLorentzVector(mom_kaon.X(),mom_kaon.Y(),mom_kaon.Z(),e_kaon);
  reco::Candidate::LorentzVector vOniaT = onia.p4() + p4_kaon;
  OniaTCand.setP4(vOniaT);

  return OniaTCand;

}

reco::Candidate::LorentzVector OniaTrakProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(OniaTrakProducer);
