// -*- C++ -*-
//
// Package:    ​RecoBTag/​SecondaryVertex
// Class:      DeepNNTagInfoProducer
//
/**\class DeepNNTagInfoProducer DeepNNTagInfoProducer.cc ​RecoBTag/DeepFlavour/plugins/DeepNNTagInfoProducer.cc
 *
 * Description: EDProducer that produces collection of DeepNNTagInfos
 *
 * Implementation:
 *    A collection of CandIPTagInfo and CandSecondaryVertexTagInfo and a CombinedSVComputer ESHandle is taken as input and a collection of DeepNNTagInfos
 *    is produced as output.
 */
//
// Original Author:  Mauro Verzetti (U. Rochester)
// Added templated functionality : Praveen C Tiwari & Jyothsna Komaragiri (IISc)
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "RecoBTag/SecondaryVertex/interface/CombinedSVComputer.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

//$$
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//$$

#include <map>

////////////////////////////////////////////////////////////
//
// TemplatedDeepNNTagInfoProducer
//
// class declaration
//
template <typename IPTag, typename SVTag>
class TemplatedDeepNNTagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit TemplatedDeepNNTagInfoProducer(const edm::ParameterSet&);
  ~TemplatedDeepNNTagInfoProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override {}

  // ----------member data ---------------------------
  const edm::EDGetTokenT<std::vector<SVTag> > svSrc_;
  CombinedSVComputer computer_;

//$$
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
//$$
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
template <typename IPTag, typename SVTag>
TemplatedDeepNNTagInfoProducer<IPTag, SVTag>::TemplatedDeepNNTagInfoProducer(const edm::ParameterSet& iConfig)
    : svSrc_(consumes<std::vector<SVTag> >(iConfig.getParameter<edm::InputTag>("svTagInfos"))),
      computer_(iConfig.getParameter<edm::ParameterSet>("computer")) {
  produces<std::vector<reco::ShallowTagInfo> >();

//$$
  puToken_  = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
  vtxToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
//$$  vtxToken_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices4D"));
  jetToken_ = consumes<edm::View<pat::Jet>>(edm::InputTag("slimmedJetsPuppi"));
//$$
}

template <typename IPTag, typename SVTag>
TemplatedDeepNNTagInfoProducer<IPTag, SVTag>::~TemplatedDeepNNTagInfoProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
template <typename IPTag, typename SVTag>
void TemplatedDeepNNTagInfoProducer<IPTag, SVTag>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get input TagInfos
  edm::Handle<std::vector<SVTag> > svTagInfos;
  iEvent.getByToken(svSrc_, svTagInfos);

  // create the output collection
  auto tagInfos = std::make_unique<std::vector<reco::ShallowTagInfo> >();

//$$
  // PU density
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
//   if (vertices->empty()) return; // skip the event if no PV found

  edm::Handle<std::vector <PileupSummaryInfo> > PUInfo;
  iEvent.getByToken(puToken_, PUInfo);

  auto npv_0_z = vertices->at(0).z();

  float PUrho = 0;
  std::vector<PileupSummaryInfo>::const_iterator ipu;
  for (ipu = PUInfo->begin(); ipu != PUInfo->end(); ++ipu) {
      if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0

      for (unsigned int i=0; i<ipu->getPU_zpositions().size(); ++i) {
          auto PU_z = (ipu->getPU_zpositions())[i];
          if ( std::abs(PU_z - npv_0_z) < 1) PUrho++;
      }
  }	
  PUrho /= 20.;
  
  // Event time
  bool PV4D = false; // if event time is taken from PV4D

  const reco::Vertex  *pv;
  pv = &(*vertices->begin());

// all these are equivalent:
//   float PVtime      = vertices->at(0).t();
//   float PVtimeError = vertices->at(0).tError();
//   float PVtime      = (*vertices)[0].t();
//   float PVtimeError = (*vertices)[0].tError();
  float PVtime      = (pv)[0].t();
  float PVtimeError = (pv)[0].tError();

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken_, jets);

  std::vector<size_t> indices(jets->size());
  for (size_t i=0; i<jets->size(); i++) indices.at(i)=i;

  edm::View<pat::Jet>::const_iterator jetIter;

  float event_time	 = 0;
  float event_timeWeight = 0;
  float event_timeNtk	 = 0;

  if ( !PV4D ) { 
    for (size_t j=0; j<indices.size(); j++) {
      size_t jetidx=indices.at(j);
      jetIter = jets->begin()+jetidx;
      const pat::Jet& jet = *jetIter;
      for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++) {
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
      if ( !PackedCandidate ) continue;
      if ( PackedCandidate->charge() == 0 ) continue;
        auto track = PackedCandidate->bestTrack();
      if ( !track ) continue;
        float track_dxy = track->dxy(pv->position());
        float track_dz  = track->dz(pv->position());
        float track_time      = track->t0();
        float track_timeError = track->covt0t0();
        float track_pt    = track->pt();
        float time_weight = track_pt * track_pt;
        if ( track_timeError > 0. && abs(track_time) < 1 
             && abs(track_dxy) < 0.05 && abs(track_dz) < 0.10 ) {
          event_timeNtk	   += 1;
          event_timeWeight += time_weight;
          event_time	   += track_time * time_weight;
        }
      }
    }
    if ( event_timeNtk > 0 ) event_time /= event_timeWeight;
    else		     event_time = -1;
  }

  else {
    if ( PVtimeError > 0. ) {
      event_timeNtk = 1;
      event_time = PVtime;
    }
    else {
      event_timeNtk = 0;
      event_time = -1;
    }
  }

//   std::cout << std::endl;
//   std::cout << " in DeepNNTagInfoProducer " << PUrho << " " << PVtimeError << " " << PVtime << std::endl;
//$$

  // loop over TagInfos
  for (auto iterTI = svTagInfos->begin(); iterTI != svTagInfos->end(); ++iterTI) {
    // get TagInfos
    const SVTag& svTagInfo = *(iterTI);
    const IPTag& ipInfo = *(iterTI->trackIPTagInfoRef().get());

    reco::TaggingVariableList vars = computer_(ipInfo, svTagInfo);
    std::vector<float> tagValList = vars.getList(reco::btau::trackEtaRel, false);
    vars.insert(reco::btau::jetNTracksEtaRel, tagValList.size());
    tagValList = vars.getList(reco::btau::trackSip2dSig, false);
    vars.insert(reco::btau::jetNSelectedTracks, tagValList.size());

//$$$$$$    vars.finalize();  //fix the TaggingVariableList, nothing should be added/removed

    //If not SV found set it to 0, not to non-existent
    if (!vars.checkTag(reco::btau::jetNSecondaryVertices))
      vars.insert(reco::btau::jetNSecondaryVertices, 0);
    if (!vars.checkTag(reco::btau::vertexNTracks))
      vars.insert(reco::btau::vertexNTracks, 0);
    
//$$
    // Jet time
    float jet_time	 = 0;
    float jet_timeWeight = 0;
    float jet_timeNtk	 = 0;
    
    // loop on jets
    bool ping = false;
    for (size_t j=0; j<indices.size(); j++) {
      size_t jetidx=indices.at(j);
      jetIter = jets->begin()+jetidx;
      const pat::Jet& jet = *jetIter;

      float deta = TMath::Abs(jet.eta() - iterTI->jet()->eta());
      float dphi = TMath::Abs(jet.phi() - iterTI->jet()->phi());
      if (dphi > 3.141593 ) dphi -= 2.*3.141593;
      if (deta < 0.01 && dphi < 0.01) {
        ping = true;

        for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++) {
          const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
        if ( !PackedCandidate ) continue;
        if ( PackedCandidate->charge() == 0 ) continue;
          auto track = PackedCandidate->bestTrack();
        if ( !track ) continue;
          float track_time      = track->t0();
          float track_timeError = track->covt0t0();
          float track_pt    = track->pt();
          float time_weight = track_pt * track_pt;
          if ( track_timeError > 0. && abs(track_time) < 1 ) {
            jet_timeNtk += 1;
            jet_timeWeight += time_weight;
            jet_time += track_time * time_weight;
          }
        }  // end loop on tracks in jets

        if ( jet_timeNtk > 0 && event_timeNtk > 0 && jet_timeWeight > 0 ) {
          jet_time = jet_time / jet_timeWeight - event_time;
          jet_time = TMath::Abs(jet_time);
        }
        else jet_time = -1; 
      }
    } // end loop on jets
    if ( !ping ) jet_time = -1; 

    // Jet Vertex time and Vertex time
    float jet_vertex_time	= 0;
    float jet_vertex_timeWeight = 0;
    float jet_vertex_timeNtk	= 0;

    float vertex_timeNtk    = 0;
    float vertex_timeWeight = 0;
    float vertex_time       = 0;

    int  nSV = -2;
    if constexpr(std::is_same_v<SVTag, reco::CandSecondaryVertexTagInfo>) {
      nSV = svTagInfo.nVertices();
      if ( nSV > 0 && ping ) {

        // loop on SV
        for ( unsigned int isv = 0; isv < svTagInfo.nVertices(); ++isv ) {
          const reco::VertexCompositePtrCandidate &vertex = svTagInfo.secondaryVertex(isv);

          // loop on jets
          for (size_t j=0; j<indices.size(); j++) {
            size_t jetidx=indices.at(j);
            jetIter = jets->begin()+jetidx;
            const pat::Jet& jet = *jetIter;

            float deta = TMath::Abs(jet.eta() - iterTI->jet()->eta());
            float dphi = TMath::Abs(jet.phi() - iterTI->jet()->phi());
            if (dphi > 3.141593 ) dphi -= 2.*3.141593;
            if (deta < 0.01 && dphi < 0.01) {

              for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++) {
                const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
              if ( !PackedCandidate ) continue;
              if ( PackedCandidate->charge() == 0 ) continue;
                auto track = PackedCandidate->bestTrack();
              if ( !track ) continue;
                float track_time   = track->t0();
                float track_timeError = track->covt0t0();
                float track_pt	= track->pt();
                float time_weight = track_pt * track_pt;

                if ( track_timeError > 0. && abs(track_time) < 1 ) {
                  bool SVtrack = false;
                  const std::vector<reco::CandidatePtr> &trks = vertex.daughterPtrVector();
                  for ( std::vector<reco::CandidatePtr>::const_iterator trk = trks.begin(); trk != trks.end(); ++trk ) {
                  if ( (*trk)->charge() != track->charge() ) continue;
                    float dpt  = TMath::Abs((*trk)->pt() / track->pt() - 1.);
                    float deta = TMath::Abs((*trk)->eta() - track->eta());
                    float dphi = TMath::Abs((*trk)->phi() - track->phi());
                    if ( dphi > 3.141593 ) dphi -= 2.*3.141593;
                    if ( dpt < 0.01 && deta < 0.01 && dphi < 0.01 ) SVtrack = true;
                  }
                  if ( SVtrack ) {
                    jet_vertex_timeNtk     += 1;
                    jet_vertex_timeWeight  += time_weight;
                    jet_vertex_time        += track_time * time_weight;
                  }
                  if ( SVtrack && isv == 0 ) {  // only one vertex for DeepCSV
                    vertex_timeNtk    += 1;
                    vertex_timeWeight += time_weight;
                    vertex_time	      += track_time * time_weight;
                  }
                }
              }  // end loop on tracks in jets
            }
          } // end loop on jets
        } // end loop on SVs 

        if ( jet_vertex_timeNtk > 0 && event_timeNtk > 0 && jet_vertex_timeWeight > 0 ) {
          jet_vertex_time = jet_vertex_time / jet_vertex_timeWeight - event_time;
          jet_vertex_time = TMath::Abs(jet_vertex_time);
        }
        else jet_vertex_time = -1; 
        if ( vertex_timeNtk > 0 && event_timeNtk > 0 && vertex_timeWeight > 0 ) {
          vertex_time = vertex_time / vertex_timeWeight - event_time;
          vertex_time = TMath::Abs(vertex_time);
        }
        else vertex_time = -1; 
      }
    }
    if ( !ping || nSV <= 0 ) {
      jet_vertex_time     = -1;
      vertex_time         = -1;
    }

    vars.insert(reco::btau::puDensity, PUrho);
    vars.insert(reco::btau::eventTime, event_time);
    vars.insert(reco::btau::jetTime, jet_time);
    vars.insert(reco::btau::jetVertexTime, jet_vertex_time);
    vars.insert(reco::btau::svTime, vertex_time);

    vars.finalize(); //fix the TaggingVariableList, nothing should be added/removed
    
//   std::cout << "    jet pt eta "  << iterTI->jet()->pt() << " " << iterTI->jet()->eta() 
//             << " PUrho " << PUrho 
// 	    << " Time evt jet jvx vtx " << event_time << " " << jet_time << " " << jet_vertex_time << " " << vertex_time 
// 	    << " nSV " << nSV 
// 	    << std::endl;
//$$

    tagInfos->emplace_back(vars, svTagInfo.jet());
  }

  // put the output in the event
  iEvent.put(std::move(tagInfos));
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <typename IPTag, typename SVTag>
void TemplatedDeepNNTagInfoProducer<IPTag, SVTag>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
//For PFJets
typedef TemplatedDeepNNTagInfoProducer<reco::CandIPTagInfo, reco::CandSecondaryVertexTagInfo> DeepNNTagInfoProducer;
//For CaloJets
typedef TemplatedDeepNNTagInfoProducer<reco::TrackIPTagInfo, reco::SecondaryVertexTagInfo> TrackDeepNNTagInfoProducer;
DEFINE_FWK_MODULE(DeepNNTagInfoProducer);
DEFINE_FWK_MODULE(TrackDeepNNTagInfoProducer);
