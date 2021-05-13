#include "EFTMaodHists.h"

EFTMaodHists::EFTMaodHists(const edm::ParameterSet& iConfig)
{
    debug = iConfig.getParameter<bool> ("debug");
    iseft = iConfig.getParameter<bool> ("iseft");

    min_pt_jet  = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep  = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");
    
    lheInfo_token_      = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfo"));
    genInfo_token_      = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GENInfo"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"));
    genJets_token_      = consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"));

    patElectrons_token_ = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("PatElectrons"));
    patMuons_token_     = consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("PatMuons"));
    patJets_token_      = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("PatJets"));
}

EFTMaodHists::~EFTMaodHists(){}

void EFTMaodHists::beginJob()
{
    total_ls = 0;
    total_events = 0;
    total_orig_xsec = 0.;
    total_sm_xsec = 0.;

    edm::Service<TFileService> newfs;

    h_maodgen_j0pt             = newfs->make<TH1EFT>("h_maodgen_j0pt","h_maodgen_j0pt",30,0,700);
    h_maodgen_njetsclean       = newfs->make<TH1EFT>("h_maodgen_njetsclean","h_maodgen_njetsclean",12,0,12);
    h_maodgen_nchlepsnotau     = newfs->make<TH1EFT>("h_maodgen_nchlepsnotau","h_maodgen_nchlepsnotau",8,0,8);
    h_maodgen_sumalllepptnotau = newfs->make<TH1EFT>("h_maodgen_sumalllepptnotau","h_maodgen_sumalllepptnotau",10,0,600);
    h_maodgen_ht               = newfs->make<TH1EFT>("h_maodgen_ht","h_maodgen_ht",20,0,1500);

    h_maodpat_j0pt             = newfs->make<TH1EFT>("h_maodpat_j0pt","h_maodpat_j0pt",30,0,700);
    h_maodpat_njetsclean       = newfs->make<TH1EFT>("h_maodpat_njetsclean","h_maodpat_njetsclean",12,0,12);
    h_maodpat_nchlepsnotau     = newfs->make<TH1EFT>("h_maodpat_nchlepsnotau","h_maodpat_nchlepsnotau",8,0,8);
    h_maodpat_sumalllepptnotau = newfs->make<TH1EFT>("h_maodpat_sumalllepptnotau","h_maodpat_sumalllepptnotau",10,0,600);
    h_maodpat_ht               = newfs->make<TH1EFT>("h_maodpat_ht","h_maodpat_ht",20,0,1500);

}

void EFTMaodHists::endJob()
{
    std::string delim = " | ";

    std::cout << "Total events processed: " << total_events << std::endl;
    std::cout << "Total LS: " << total_ls << std::endl;
    std::cout << "                " 
        << std::setw(13) << "xsec / LS" << delim
        << std::setw(13) << "xsec / Events" << delim
        << std::setw(13) << "xsec / 500"
        << std::endl;
    std::cout << "SM Xsec (pb):   "
        << std::setw(13) << total_sm_xsec / double(total_ls) << delim
        << std::setw(13) << total_sm_xsec / double(total_events) << delim
        << std::setw(13) << total_sm_xsec / 500.
        << std::endl;
    std::cout << "Orig Xsec (pb): "
        << std::setw(13) << total_orig_xsec / double(total_ls) << delim
        << std::setw(13) << total_orig_xsec / double(total_events) << delim
        << std::setw(13) << total_orig_xsec / 500.
        << std::endl;

    //WCPoint* rwgt_pt = new WCPoint("rwgt",0.0);

}

void EFTMaodHists::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    eventcount++;

    // set tree vars to default values
    initialize_variables();

    // Gen
    edm::Handle<LHEEventProduct>             LHEInfo;
    edm::Handle<reco::GenParticleCollection> prunedParticles;
    edm::Handle<std::vector<reco::GenJet>>   genJets;
    event.getByToken(lheInfo_token_,LHEInfo);
    event.getByToken(genParticles_token_,prunedParticles);
    event.getByToken(genJets_token_,genJets);

    // Reco
    edm::Handle<std::vector<pat::Electron>> patElectrons;
    edm::Handle<std::vector<pat::Muon>>     patMuons;
    edm::Handle<std::vector<pat::Jet>>      patJets;
    event.getByToken(patElectrons_token_,patElectrons);
    event.getByToken(patMuons_token_,patMuons);
    event.getByToken(patJets_token_,patJets);

    // Gen objects, baseline cuts and cleaning
    reco::GenParticleCollection gen_leptons  = GetGenLeptons(*prunedParticles);
    gen_leptons                              = MakePtEtaCuts(gen_leptons,15.0,2.5);
    std::vector<reco::GenJet> gen_jets       = MakePtEtaCuts(*genJets,30.0,2.5); // Pt, eta cuts
    std::vector<reco::GenJet> gen_jets_clean = CleanGenJets(gen_jets,gen_leptons,0.4);

    // Reco objects, baseline cuts and cleaning
    std::vector<pat::Electron> pat_electrons = *patElectrons;
    std::vector<pat::Muon> pat_muons         = *patMuons;
    pat_electrons                            = MakePtEtaCuts(pat_electrons,15,2.5); 
    pat_muons                                = MakePtEtaCuts(pat_muons,15,2.5);
    std::vector<pat::Jet> pat_jets           = MakePtEtaCuts(*patJets,30.0,2.5); // Pt, eta cuts
    std::vector<pat::Jet> pat_jets_clean     = CleanGenJets(pat_jets,pat_electrons,0.4); // Clean jets against e
    pat_jets_clean                           = CleanGenJets(pat_jets,pat_muons,0.4);     // Clean jets against mu

    std::cout << " \nStart event:\n" << std::endl;

    // Get WC fit
    originalXWGTUP_intree = LHEInfo->originalXWGTUP();  // original cross-section
    double sm_wgt = 0.;
    std::vector<WCPoint> wc_pts;
    if (iseft) {// Add EFT weights 
        for (auto wgt_info: LHEInfo->weights()) {
            auto LHEwgtstr = std::string(wgt_info.id);
            std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
            if (foundstr!=std::string::npos) {
                WCPoint wc_pt(wgt_info.id,wgt_info.wgt);
                wc_pts.push_back(wc_pt);
                if (wc_pt.isSMPoint()) {
                    sm_wgt = wgt_info.wgt;
                }
            }
        }
    } else {
        sm_wgt = originalXWGTUP_intree;
        WCPoint wc_pt("smpt",sm_wgt);
        wc_pts.push_back(wc_pt);
    }
    WCFit eft_fit(wc_pts,"");


    double gen_j0_pt = 0;
    double gen_ht = 0;
    double gen_leppt = 0;
    double nchlepsnotau = 0;

    double pat_j0_pt = 0;
    double pat_ht = 0;
    double pat_leppt = 0;
    double pat_mpt = 0;
    double pat_ept = 0;

    reco::GenParticleCollection gen_e_mu = rmParticleType(gen_leptons,{15,12,14,16});
    //gen_leptons = rmParticleType(gen_leptons,15,12,14,16);

    // GEN
    for (size_t i = 0; i < gen_jets_clean.size(); i++) {
        reco::GenJet j = gen_jets_clean.at(i);
        gen_ht = gen_ht + j.p4().Pt();
        //std::cout << "\ngen j pt: " << j.p4().Pt() << std::endl;
        if (i == 0){
            gen_j0_pt = j.p4().Pt();
        }
    }
    for (auto l: gen_e_mu){
    //for (auto l: gen_leptons){
        /*
        if (l.charge() != 0){           // Skip nu
            if (abs(l.pdgId()) != 15){  // Skip tau
                gen_leppt = gen_leppt + l.p4().Pt();
                nchlepsnotau = nchlepsnotau + 1;
            }
            else {
                std::cout << "Skipping tau!" << std::endl;
            }
        }
        else {
            std::cout << "Skipping nu!" << std::endl;
        }
        */
        gen_leppt = gen_leppt + l.p4().Pt();
        nchlepsnotau = nchlepsnotau + 1;
        if ( (abs(l.pdgId()) != 11) and (abs(l.pdgId()) != 13) ) {
            std::cout << "\nError: Not a mu or an e, not sure how this happened. Go check the logic. ID: " << l.pdgId() << "\n\n" << std::endl;
            throw std::runtime_error("");
        }
    }

    //if ( (gen_jets_clean.size()>=2) and 
    h_maodgen_j0pt->Fill(gen_j0_pt,1.0,eft_fit);
    h_maodgen_ht->Fill(gen_ht,1.0,eft_fit);
    h_maodgen_sumalllepptnotau->Fill(gen_leppt,1.0,eft_fit);

    h_maodgen_nchlepsnotau->Fill(nchlepsnotau,1.0,eft_fit);
    h_maodgen_njetsclean->Fill(gen_jets_clean.size(),1.0,eft_fit);

    // PAT
    for (size_t i = 0; i < pat_jets_clean.size(); i++) {
        pat::Jet j = pat_jets_clean.at(i);
        pat_ht = pat_ht + j.p4().Pt();
        //std::cout << "\npat j pt: " << j.p4().Pt() << std::endl;
        if (i == 0){
            pat_j0_pt = j.p4().Pt();
        }
    }
    for (auto e: pat_electrons){
        pat_leppt = pat_leppt + e.p4().Pt();
        pat_ept = pat_ept + e.p4().Pt();
        //std::cout << "e: " << e.pdgId() << std::endl;
    }
    for (auto m: pat_muons){
        pat_leppt = pat_leppt + m.p4().Pt();
        pat_mpt = pat_mpt + m.p4().Pt();
        //std::cout << "m: " << m.pdgId() << std::endl;
    }
    //std::cout << "comp to sum: " << pat_ept + pat_mpt << " " << pat_leppt << std::endl;
    h_maodpat_j0pt->Fill(pat_j0_pt,1.0,eft_fit);
    h_maodpat_ht->Fill(pat_ht,1.0,eft_fit);
    h_maodpat_sumalllepptnotau->Fill(pat_leppt,1.0,eft_fit);

    h_maodpat_njetsclean->Fill(pat_jets_clean.size(),1.0,eft_fit);
    h_maodpat_nchlepsnotau->Fill(pat_electrons.size()+pat_muons.size(),1.0,eft_fit);

    std::cout << "max j pt: " << gen_j0_pt << " " << pat_j0_pt << std::endl;
    //std::cout << "nleps getgenleps: " << GetGenLeptons(pat_electrons).size()+GetGenLeptons(pat_muons).size() << std::endl;
    std::cout << " " << std::endl;

}

void EFTMaodHists::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTMaodHists::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    total_events += eventcount;
}

void EFTMaodHists::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup)
{
    total_ls += 1;
}

void EFTMaodHists::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTMaodHists);

