#ifndef EFTGENREADER_GENREADER_SELECTIONANALYZER_h
#define EFTGENREADER_GENREADER_SELECTIONANALYZER_h

#include <cstdlib>
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <boost/any.hpp>

#include <iostream>
#include <algorithm>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TRandom3.h>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TVector.h"
#include "TLorentzVector.h"

// Framework
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PyBind11ProcessDesc.h"
#include "DataFormats/Common/interface/Handle.h"

// Physics
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// end includes
// -----------------------------------------------

class EFTMaodHists: public edm::one::EDAnalyzer<>
{
    private:
        // EDAnalyzer-specific:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

        void tree_add_branches();
        void initialize_variables();

        //void dumpParticles(const reco::GenParticleCollection& particles);
        //void dumpJets(const std::vector<reco::GenJet>& jets);
        const reco::Candidate* GetGenMotherNoFsr(const reco::Candidate* p);
        //std::pair<const reco::Candidate*, const reco::Candidate*> GetGenDaughterNoFsr(const reco::Candidate* p);

        //std::vector<reco::GenJet>   MakeJetPtEtaCuts(const std::vector<reco::GenJet>& gen_jets, double pt_cut, double eta_cut);
        //reco::GenParticleCollection MakePtEtaCuts(const reco::GenParticleCollection& gen_particles, double pt_cut, double eta_cut);
        //reco::GenParticleCollection GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id);
        //reco::GenParticleCollection GetGenElectrons(const reco::GenParticleCollection& gen_particles);
        //reco::GenParticleCollection GetGenMuons(const reco::GenParticleCollection& gen_particles);
        reco::GenParticleCollection GetGenLeptons(const reco::GenParticleCollection& gen_particles);
        std::vector<reco::GenJet> GetGenJets(const std::vector<reco::GenJet>& inputs);
        //std::vector<pat::Electron> MakePtEtaCutsPatElectrons(const std::vector<pat::Electron>& pat_electrons, double pt_cut, double eta_cut);
        //std::vector<pat::Muon> MakePtEtaCutsPatMuons(const std::vector<pat::Muon>& pat_muon, double pt_cut, double eta_cut);
        //std::vector<pat::Jet> MakePatJetPtEtaCuts(const std::vector<pat::Jet>& pat_jet, double pt_cut, double eta_cut);
        //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> getSumTLV(reco::GenParticleCollection col);
        double getdPhi(reco::GenParticle p1, reco::GenParticle p2);
        //double getdR(reco::GenParticle p1, reco::GenParticle p2);
        //double getInvMass(reco::GenParticle p1, reco::GenParticle p2);

    public:
        explicit EFTMaodHists(const edm::ParameterSet&);
        ~EFTMaodHists();

        template <typename T> edm::Handle<T> get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token);

        template <typename T1, typename T2> std::vector<T1> CleanGenJets(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize);
        template <typename T1, typename T2> double getdR(T1 p1, T2 p2);
        template <typename T> std::vector<T> MakePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta);
        template <typename T> std::vector<T> rmParticleType(const std::vector<T>& input, std::vector<int> rmid);
        template <typename T> std::vector<T> keepParticleTypes(const std::vector<T>& input, std::vector<int> keep_id);

        std::ofstream fout;
        FILE * ffout;

        std::string sampleName;
        int sampleNumber;
        int eventcount;

        // runtime configurable parameters
        bool iseft;
        bool debug;
        double min_pt_jet;
        double min_pt_lep;
        double max_eta_jet;
        double max_eta_lep;

        edm::EDGetTokenT<LHEEventProduct> lheInfo_token_;
        edm::EDGetTokenT<GenEventInfoProduct> genInfo_token_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticles_token_;  // reco::GenParticlesCollection is an alias for std::vector<reco::GenParticle>>
        edm::EDGetTokenT<std::vector<reco::GenJet> > genJets_token_;

        edm::EDGetTokenT<std::vector<pat::Electron> > patElectrons_token_;
        edm::EDGetTokenT<std::vector<pat::Muon> > patMuons_token_;
        edm::EDGetTokenT<std::vector<pat::Jet> > patJets_token_;

        std::unordered_map<std::string,WCFit*> lep_fits_dict;
        std::vector<std::string> lep_category_names_vect; // keys in the dict (but do not have to be names of the fit tags)

        // declare the tree
        TTree * summaryTree;

        // tree branches
        double originalXWGTUP_intree;
        int eventnum_intree;
        int lumiBlock_intree;
        int runNumber_intree;

        WCFit eft_wgt_intree;
        double sm_wgt_intree;

        // Misc. counters
        int total_ls;
        int total_events;
        double total_orig_xsec;
        double total_sm_xsec;

        TH1EFT* h_maodgen_j0_pt;
        TH1EFT* h_maodgen_j0_eta;
        TH1EFT* h_maodgen_l0_pt;
        TH1EFT* h_maodgen_njetsclean;
        TH1D*   h_2l2j_counts;

};

void EFTMaodHists::tree_add_branches()
{
    summaryTree->Branch("originalXWGTUP",&originalXWGTUP_intree);
    summaryTree->Branch("eventnum",&eventnum_intree);
    summaryTree->Branch("lumiBlock",&lumiBlock_intree);
    summaryTree->Branch("runNumber",&runNumber_intree);

    summaryTree->Branch("eft_wgt",&eft_wgt_intree);
    summaryTree->Branch("sm_wgt",&sm_wgt_intree);
    //summaryTree->Branch("n_jets",&n_jet_intree);
    //summaryTree->Branch("n_bjets",&n_bjet_intree);
}

void EFTMaodHists::initialize_variables()
{
    originalXWGTUP_intree = -99;
    eventnum_intree = -1;
    lumiBlock_intree = -1;
    runNumber_intree = -1;
}

reco::GenParticleCollection EFTMaodHists::GetGenLeptons(const reco::GenParticleCollection& gen_particles) {
    reco::GenParticleCollection gen_leptons;
    bool is_lepton;
    bool is_neutrino;
    for (size_t i = 0; i < gen_particles.size(); i++) {
        const reco::GenParticle& p = gen_particles.at(i);
        const reco::Candidate* p_mom_noFSR = GetGenMotherNoFsr(&p); // This will walk up the chain and so doesn't get direct mothers
        const reco::Candidate* p_gmom_noFSR = GetGenMotherNoFsr(p_mom_noFSR);
        const reco::Candidate* p_mom = p.mother();  // The direct mother
        int id = p.pdgId();
        int gmom_noFSR_id = p_gmom_noFSR->pdgId();
        is_lepton = (abs(id) == 11 || abs(id) == 13 || abs(id) == 15);
        is_neutrino = (abs(id) == 12 || abs(id) == 14 || abs(id) == 16);

        if (!is_lepton && !is_neutrino) {
            continue;
        }

        int mom_id = id;    // If no mother, set to own id
        if (p_mom) mom_id = p_mom->pdgId();

        bool is_fromEWBoson = (abs(mom_id) >= 22 && abs(mom_id) <= 25);
        bool is_fromTopSystem = (abs(mom_id) == 24 && abs(gmom_noFSR_id) == 6);

        bool is_hard_process = p.isHardProcess();
        bool is_fromHardGluon = (abs(mom_id) == 21 && p_mom->status() == 21);   // status == 21 corresponds to incoming hard particle
        bool is_fromHardQCD = (abs(mom_id) >= 1 && abs(mom_id) <= 5 && p_mom->status() == 21);

        if (is_fromEWBoson) {
            gen_leptons.push_back(p);
        } else if (is_hard_process && is_fromHardGluon) {
            gen_leptons.push_back(p);
        } else if (is_hard_process && is_fromHardQCD) {
            gen_leptons.push_back(p);
        }
    }
    std::sort(gen_leptons.begin(),gen_leptons.end(), [] (reco::GenParticle a, reco::GenParticle b) { return a.p4().Pt() > b.p4().Pt();});
    return gen_leptons;
}

// This is ill-defined for particles with multiple mothers
const reco::Candidate* EFTMaodHists::GetGenMotherNoFsr(const reco::Candidate* p) {
    if (p->numberOfMothers()) {
        const reco::Candidate* mom = p->mother(0);
        if (mom->pdgId() != p->pdgId()) {
            return mom;
        } else {
            return GetGenMotherNoFsr(mom);
        }
    } else {
        return p;
    }
}

std::vector<reco::GenJet> EFTMaodHists::GetGenJets(const std::vector<reco::GenJet>& inputs) {
    std::vector<reco::GenJet> ret;
    for (size_t i = 0; i < inputs.size(); i++) {
        const reco::GenJet& j = inputs.at(i);
        if (j.p4().Pt() < min_pt_jet) {
            continue;
        } else if (max_eta_jet > 0.0 && fabs(j.eta()) >= max_eta_jet) {
            continue;
        }
        ret.push_back(j);
    }
    std::sort(ret.begin(),ret.end(), [] (reco::GenJet a, reco::GenJet b) { return a.p4().Pt() > b.p4().Pt();});
    return ret;
}

// Clean objects based on how close they are to other objects
template <typename T1, typename T2>
std::vector<T1> EFTMaodHists::CleanGenJets(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize) {
    std::vector<T1> cleaned_obj1_vector;
    bool isClean;
    for (size_t i = 0; i < obj1_vector.size(); i++){
        isClean = true;
        const T1& obj1 = obj1_vector.at(i);
        for (size_t j = 0; j < obj2_vector.size(); j++){
            const T2& obj2 = obj2_vector.at(j);
            if (getdR(obj1,obj2) <= coneSize){
                isClean = false;
            }
        }
        if (isClean) {
            cleaned_obj1_vector.push_back(obj1);
        }
    }
    return cleaned_obj1_vector;
}

// Get dR between two objects
template <typename T1, typename T2>
double EFTMaodHists::getdR(T1 p1, T2 p2) {
    double dR = (p1.p4().Eta() - p2.p4().Eta())*(p1.p4().Eta() - p2.p4().Eta());
    dR += (getdPhi(p1,p2)*getdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

double EFTMaodHists::getdPhi(reco::GenParticle p1, reco::GenParticle p2) {
    double dPhi = p2.p4().Phi() - p1.p4().Phi();
    double pi = 2.0*asin(1.0);
    if (dPhi>=pi) dPhi = dPhi - 2.0*pi;     // from TVector2
    if (dPhi<(-pi)) dPhi = dPhi + 2.0*pi;   // from TVector2
    return dPhi;
}

// Template function for making pt and eta cuts
template <typename T>
std::vector<T> EFTMaodHists::MakePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta) {
    std::vector<T> ret;
    for (size_t i = 0; i < input.size(); i++) {
        const T& p = input.at(i);
        if (p.p4().Pt() < min_pt) {
            continue;
        } else if (max_eta > 0.0 && fabs(p.eta()) >= max_eta) {
            continue;
        }
        ret.push_back(p);
    }
    std::sort(ret.begin(),ret.end(), [] (T a, T b) { return a.p4().Pt() > b.p4().Pt();});
    return ret;
}

// Template function for selecting particles with a certain ID
template <typename T>
std::vector<T> EFTMaodHists::keepParticleTypes(const std::vector<T>& input, std::vector<int> keep_ids) {
    std::vector<T> ret;
    for (size_t i = 0; i < input.size(); i++) {
        const T& p = input.at(i);
        int p_id = p.pdgId();
        bool keep_p = false;
        for (auto keep_id: keep_ids) {
            if (abs(p_id) == abs(keep_id)){
                keep_p = true;
            }
        }
        if (keep_p) {
            ret.push_back(p);
            //std::cout << "Keeping " << p_id << std::endl;
        } else {
            //std::cout << "Skipping: " << p_id << std::endl;
        }
    }
    return ret;
}

// Template function for getting rid of particles with a certain ID
template <typename T>
std::vector<T> EFTMaodHists::rmParticleType(const std::vector<T>& input, std::vector<int> rmids) {
    std::vector<T> ret;
    for (size_t i = 0; i < input.size(); i++) {
        const T& p = input.at(i);
        int p_id = p.pdgId();
        bool skip_p = false;
        for (auto rmid: rmids) {
            if (abs(p_id) == abs(rmid)){
                //std::cout << "Skipping this particle, it's on the skip list: " << p_id << std::endl;
                skip_p = true;
            }
        }
        if (not skip_p){
            ret.push_back(p);
        }
    }
    return ret;
}

#endif
