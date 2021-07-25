// root -l -b -q exploreTH1EFT.C

void exploreTH1EFT(TString hist_path){

    // Open the file
    //TString hist_path = "/afs/crc.nd.edu/user/k/kmohrman/EFT/root_files/EFTMaodHists_output_tree.root"; // Total events processed: 1610
    TFile* maod_hist = TFile::Open(hist_path);
    maod_hist->Print();

    // Make the rwgt point
    std::string top19001_hi_str = "rwgt_cQei_4.59_cQl3i_8.97_cQlMi_4.99_cQq11_23.0_cQq13_4.0_cQq81_-13.0_cQq83_6.0_cbW_4.95_cpQ3_3.48_cpQM_21.65_cpt_12.31_cptb_12.63_ctG_1.18_ctW_2.87_ctZ_3.15_ctei_4.86_ctlSi_6.52_ctlTi_0.84_ctli_4.82_ctp_44.26_ctq1_-19.0_ctq8_10.0";
    std::string all_2_str = "rwgt_cQei_2.0_cQl3i_2.0_cQlMi_2.0_cQq11_2.0_cQq13_2.0_cQq81_2.0_cQq83_2.0_cbW_2.0_cpQ3_2.0_cpQM_2.0_cpt_2.0_cptb_2.0_ctG_2.0_ctW_2.0_ctZ_2.0_ctei_2.0_ctlSi_2.0_ctlTi_2.0_ctli_2.0_ctp_2.0_ctq1_2.0_ctq8_2.0";
    WCPoint* rwgt_pt = new WCPoint(all_2_str);
    //WCPoint* rwgt_pt = new WCPoint("rwgt_ctG_100.0");
    //WCPoint* rwgt_pt = new WCPoint("rwgt_cpt_100.0");
    WCPoint* sm_pt = new WCPoint("smpt");

    // Loop though the hists and plot them
    TDirectory* td = maod_hist->GetDirectory("EFTMaodHists");
    TIter next(td->GetListOfKeys());
    TKey* key; TLegend* leg;
    while ((key = (TKey*)next())) {


        TCanvas *c1 = new TCanvas("c1","",1000,1000);

        TString s = key->GetName();
        std::cout << "Hist name: " << s << std::endl;

        TH1EFT* h = (TH1EFT*)td->Get(key->GetName());

        // Check if it's a TH1EFT
        bool is_TH1EFT = false;
        if (h->IsA()->InheritsFrom(TH1EFT::Class())){
            is_TH1EFT = true;
        }

        // Reweight
        if (is_TH1EFT) {
            std::cout << "Reweighting hist" << std::endl;
            for (Int_t bin_idx = 0; bin_idx <= h->GetNbinsX()+1; bin_idx++) {
                double wcfit_bin_val_rwgt = h->GetBinFit(bin_idx).evalPoint(sm_pt);
                double wcfit_bin_err_rwgt = h->GetBinFit(bin_idx).evalPointError(sm_pt);
                h->SetBinContent(bin_idx,wcfit_bin_val_rwgt);
                h->SetBinError(bin_idx,wcfit_bin_err_rwgt);
                WCFit fit = h->GetBinFit(bin_idx);
                std::vector<std::string> names = fit.getNames(); 
                std::cout << "The WC names:" << std::endl;
                for (auto name: names){
                    std::cout << "\t" << name << std::endl;
                }
            }
        }

        return;
        h->SetLineColor(kRed);
        h->Draw("E SAME");
        //h->Draw("SAME,HIST");
        c1->Print(s+".png","png");

        delete c1;

    }

    delete sm_pt;
    delete rwgt_pt;

}
