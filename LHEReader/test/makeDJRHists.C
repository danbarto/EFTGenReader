#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TChain.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "EFTGenReader/EFTHelperUtilities/interface/Stopwatch.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"

void printProgress(int current_index, int total_entries, int interval=20) {
    if (current_index % max(int(total_entries*interval/100.),interval) == 0) {
        float fraction = 100.*current_index/total_entries;
        std::cout << int(fraction) << " % processed " << std::endl;
    }
}

void add_wcpts_to_dict(std::string wc_name , double min , double max , double delta , std::map<string,std::pair<string,double>> &rwgt_dict){
    double nsteps = (max-min)/delta;
    double wc_val;
    for (int i = 0; i <= nsteps; i++) {
        wc_val = min + i*delta;
        string dict_key = wc_name + "_" + std::to_string(wc_val);
        rwgt_dict[dict_key] = std::make_pair(wc_name,wc_val);
    }
}

void setcanvas(TCanvas *c1, TPad **pad) {
    c1->SetLeftMargin(0.0);
    c1->SetTopMargin(0.00);
    c1->SetRightMargin(0.00);
    c1->SetBottomMargin(0.0);

    pad[0] = new TPad("pad0","pad",0,0.5,0.5,1.0);
    pad[1] = new TPad("pad1","pad",0.5,0.5,1.0,1.0);
    pad[2] = new TPad("pad2","pad",0,0,0.5,0.5);
    pad[3] = new TPad("pad3","pad",0.5,0.0,1.0,0.5);
    for(int k=0; k<4; k++) {
        pad[k]->Draw();
    }
    return;
}

void setlegend(TLegend *legend, TH1D *hall, TH1D *hmult0, TH1D *hmult1, TH1D *hmult2, TH1D *hmult3) {
    legend->SetTextSize(0.055);
    legend->SetBorderSize(0);
    legend->SetLineColor(0);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);

    legend->AddEntry(hmult0,"0 partons");
    legend->AddEntry(hmult1,"1 parton");
    legend->AddEntry(hall,"Total");
    //legend->AddEntry(hmult2,"2 partons");
    //legend->AddEntry(hmult3,"3 partons");
    return;
}

std::vector<TH1EFT*> makeTH1EFTs(const char *name, TChain *tree, int djr_idx, string rwgt_str, int nbins, double xlow, double xhigh) {

    //std::string sample_type = "central_sample";
    std::string sample_type = "EFT_sample";

    std::string cat = "no_cuts";
    //std::string cat = "2lep";
    //std::string cat = "3lep";

    // Set max number of events for event loop
    //int max_events = -1; // Run over all events
    //int max_events = 30000; // Debug
    //int max_events = 200000;
    int max_events = 100;

    Stopwatch sw;
    std::set<int> unique_runs;

    TH1EFT *hall   = new TH1EFT(TString::Format("hall_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult0 = new TH1EFT(TString::Format("hmult0_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult1 = new TH1EFT(TString::Format("hmult1_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult2 = new TH1EFT(TString::Format("hmult2_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult3 = new TH1EFT(TString::Format("hmult3_%s",name),"",nbins,xlow,xhigh);

    int nMEpartons_intree;
    int nMEpartonsFiltered_intree;
    double genWgt_intree;
    std::vector<double> *djrvalues_intree = 0;
    WCFit* wcFit_intree = 0;
    double origxsec_intree;
    int lumiBlock_intree;

    double genLep_pt1_intree;
    double genLep_pt2_intree;
    double genLep_pt3_intree;
    double genJet_pt1_intree;
    double genJet_pt2_intree;
    double genJet_pt3_intree;
    double genJet_pt4_intree;

    tree->SetBranchAddress("nMEpartons",&nMEpartons_intree);
    tree->SetBranchAddress("nMEpartonsFiltered",&nMEpartonsFiltered_intree);
    tree->SetBranchAddress("genWgt",&genWgt_intree);
    tree->SetBranchAddress("DJRValues",&djrvalues_intree);
    tree->SetBranchAddress("wcFit",&wcFit_intree);
    tree->SetBranchAddress("originalXWGTUP",&origxsec_intree);
    tree->SetBranchAddress("lumiBlock",&lumiBlock_intree);

    tree->SetBranchAddress("genLep_pt1",&genLep_pt1_intree);
    tree->SetBranchAddress("genLep_pt2",&genLep_pt2_intree);
    tree->SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
    tree->SetBranchAddress("genJet_pt1",&genJet_pt1_intree);
    tree->SetBranchAddress("genJet_pt2",&genJet_pt2_intree);
    tree->SetBranchAddress("genJet_pt3",&genJet_pt3_intree);
    tree->SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

    WCFit inclusive_fit; 
    WCPoint* wc_pt = new WCPoint(rwgt_str,0);

    // Event loop:
    int tot_events = tree->GetEntries();
    sw.start("Full loop");
    int n_tot = 0;
    int n_count = 0;
    for (int i = 0; i < tot_events; i++) {
        n_tot = i;
        if (i == max_events){
            break;
        }
        sw.start("Event loop");
        tree->GetEntry(i);
        unique_runs.insert(lumiBlock_intree);
        if (sample_type == "EFT_sample"){
            if (djr_idx >= djrvalues_intree->size()) {
                continue;
            }
        }

        if (cat != "no_cuts"){
            // Cut on pt (this applies to both 2lep and 3lep catagories)
            if (genLep_pt1_intree < 25.0) {
                continue; // All events need at least 1 lepton with pt greater than 25
            } else if (genLep_pt2_intree < 15.0) {
                continue; // All events need at least 2 leptons (wiht pt greater than 15)
            }
            if (cat == "2lep") {
                // 2 leptons + 4 jets
                if (genLep_pt3_intree > 10.0) {
                    continue; // We want exactly 2 leptons with pt greater than 10 (skip if 3 leptons in event)
                } else if (genJet_pt4_intree < 30.0) {
                    continue; // We want at least 4 jets with pt greater than 30
                }
            } else if (cat == "3lep") { 
                // We won't split up 3 and 4 lep catagories
                if (genLep_pt3_intree < 10) {
                    continue; // We want at least 3 leptons with pt greater than 10
                } else if (genJet_pt2_intree < 30.0 ) {
                    continue; // We want at least 2 jets with pt greater than 30
                }
            } 
        }

        if (sample_type == "central_sample") { // For any sample w/o EFT rwgting
            //std::cout << origxsec_intree << std::endl;
            std::vector<WCPoint> tmp_pts;
            WCPoint tmp_SM_pt("smpt",origxsec_intree);
            tmp_pts.push_back(tmp_SM_pt);
            WCFit tmp_fit(tmp_pts,"");
            inclusive_fit.addFit(tmp_fit);
        }

        if (sample_type == "EFT_sample"){
            double djr_val = log10(djrvalues_intree->at(djr_idx)); // Remember to un comment when not doing central sampel!!!
            inclusive_fit.addFit(*wcFit_intree);

            printProgress(i,tot_events,10);

            double eft_weight = wcFit_intree->evalPoint(wc_pt);

            bool mult0 = (nMEpartonsFiltered_intree == 0);
            bool mult1 = (nMEpartonsFiltered_intree == 1);
            bool mult2 = (nMEpartonsFiltered_intree == 2);
            bool mult3 = (nMEpartonsFiltered_intree == 3);

            sw.start("Hist fill");
            
            wcFit_intree->scale(genWgt_intree);
            if (genWgt_intree != 0) {
                hall->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                if (mult0) {
                    hmult0->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult1) {
                    hmult1->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult2) {
                    hmult2->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult3) {
                    hmult3->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                }
            }
            sw.lap("Hist fill");
            sw.lap("Event loop");
        }
        n_count = n_count + 1;

    }
    sw.stop("Full loop");
    sw.readAllTimers(true,"");
    sw.readAllTimers(false,"");

    double norm_factor;
    if (sample_type == "EFT_sample"){
        norm_factor = unique_runs.size(); // number of lumi blocks
    } else if (sample_type == "central_sample"){
        norm_factor = max_events; // this one is correct?
        //norm_factor = n_count;
    }
    //int nlumi_blocks = unique_runs.size();
    // Divide by number of groups of 500 events (e.g. 1000000 events / 500 = 200 blocks)
    /*
    double incl_xsec = inclusive_fit.evalPoint(wc_pt) / nlumi_blocks; 
    double incl_xsec_err = inclusive_fit.evalPointError(wc_pt) / nlumi_blocks;
    */
    double incl_xsec = inclusive_fit.evalPoint(wc_pt) / norm_factor; 
    double incl_xsec_err = inclusive_fit.evalPointError(wc_pt) / norm_factor;

    std::cout << "\nXsec: " << incl_xsec << " +- " << incl_xsec_err 
              //<< " pb (events: " << n_count << "/" << max_events << " -> " << (float)n_count/max_events 
              << " pb (events: " << n_count << "/" << n_tot << " -> " << (float)n_count/n_tot 
              << ", type: " << sample_type 
              << ", cat: " << cat 
              << ", rwgt_str: \"" << rwgt_str << "\")\n"
              << std::endl;

    vector<TH1EFT*> hist_list; 
    hist_list.push_back(hall);
    hist_list.push_back(hmult0);
    hist_list.push_back(hmult1);
    hist_list.push_back(hmult2);
    hist_list.push_back(hmult3);

    delete wc_pt;
    return hist_list;

}

void makeplot(WCPoint* wc_pt,  const char* xlabel, std::vector<TH1EFT*> hist_list) {

    TH1EFT* hall   = hist_list[0];
    TH1EFT* hmult0 = hist_list[1];
    TH1EFT* hmult1 = hist_list[2];
    TH1EFT* hmult2 = hist_list[3];
    TH1EFT* hmult3 = hist_list[4];

    //WCPoint* wc_pt = new WCPoint(rwgt_str,0);
    for (Int_t bin_idx = 0; bin_idx <= hall->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hall->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hall->GetBinFit(bin_idx).evalPointError(wc_pt);
        hall->SetBinContent(bin_idx,wcfit_bin_val);
        hall->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult0->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult0->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult0->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult0->SetBinContent(bin_idx,wcfit_bin_val);
        hmult0->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult1->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult1->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult1->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult1->SetBinContent(bin_idx,wcfit_bin_val);
        hmult1->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult2->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult2->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult2->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult2->SetBinContent(bin_idx,wcfit_bin_val);
        hmult2->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult3->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult3->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult3->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult3->SetBinContent(bin_idx,wcfit_bin_val);
        hmult3->SetBinError(bin_idx,wcfit_bin_err);
    }

    double max_y = hall->GetBinContent(hall->GetMaximumBin());
    hall->GetYaxis()->SetRangeUser(0.01,50*(hall->GetBinContent(hall->GetMaximumBin()))); // Do not use GetMaximum() here since it will always be whatever we set it to first

    //hall->SetLineColor(921);
    //hmult0->SetLineColor(600);
    //hmult1->SetLineColor(629);
    hall->SetLineColor(14);
    hmult0->SetLineColor(4);
    hmult1->SetLineColor(2);
    hmult2->SetLineColor(419);
    hmult3->SetLineColor(810);

    hall->SetLineWidth(4);
    hmult0->SetLineWidth(2);
    hmult1->SetLineWidth(2);
    hmult0->SetLineStyle(2);
    hmult1->SetLineStyle(2);
    hmult2->SetLineStyle(2);
    hmult3->SetLineStyle(2);

    hall->GetXaxis()->SetTitle(xlabel);

    // Adjust sizes
    hall->GetXaxis()->SetTitleSize(0.08);
    hall->GetXaxis()->SetTitleOffset(0.8);
    hall->GetXaxis()->SetLabelSize(0.055);
    hall->GetYaxis()->SetLabelSize(0.055);
    hall->GetXaxis()->SetLabelOffset(0.009);
    hall->GetYaxis()->SetLabelOffset(0.009);

    //TLegend *legend=new TLegend(0.67,0.87-4*0.06,0.87,0.87); // left,top,right,bottom
    //TLegend *legend=new TLegend(0.14,0.87-3*0.06,0.32,0.87); // left,top,right,bottom
    TLegend *legend=new TLegend(0.14,0.89,0.32,0.9-0.07*3); // left,top,right,bottom
    setlegend(legend, hall, hmult0, hmult1, hmult2, hmult3);

    hall->Draw("EHIST");
    hmult0->Draw("EHISTSAME");
    hmult1->Draw("EHISTSAME");
    //hmult2->Draw("EHISTSAME");
    //hmult3->Draw("EHISTSAME");
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);

    legend->Draw();
    return;
}

void makeDJRHists(const TString & infile_spec, const TString & outfile, bool basePtSM, bool basePtRefPt, const TString & proc_name) {

    // Configure the scan
    bool IndividualWCScan = true;  // Fill the dictionary if we are doing a scan for the ctG plus the 2light2heavy WCs
    bool WCscan      = true;       // Fill the dictionary if we are testing all WC at their ref values

    if ( basePtSM == 1 and basePtRefPt == 1 ) {
        std::cout << "\nERROR: Cannot set both basePtSM and basePtRefPt to true - please choose one base point to vary the WC with respect to\n" << std::endl;
        return; 
    } else if ( basePtSM == 0 and basePtRefPt == 0 ){
        std::cout << "\nERROR: No base point chosen - please choose one base point to vary the WC with respect to\n" << std::endl;
        return; 
    } else if  ( basePtSM ){
        std::cout << "\nBase point selected: SM point\n" << std::endl;
    } else if ( basePtRefPt ) {
        std::cout << "\nBase point selected: Ref point\n" << std::endl;
    }

    // Get the tree
    TH1::SetDefaultSumw2();
    TChain *tree = new TChain("EFTLHEReader/summaryTree");
    std::ifstream infiles(infile_spec);
    TString fn;
    while (infiles >> fn) {
        tree->Add(fn);
    }

    int nbins = 50.;
    double djrmin = -0.5;
    double djrmax = 3.;


    // Make the TH1EFts:
    std::vector<TH1EFT*> hist_list0 = makeTH1EFTs("djr0",tree,0,"",nbins,djrmin,djrmax);
    ////std::vector<TH1EFT> hist_list0 =  makeTH1EFTs("djr0",tree,0,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list0 = makeTH1EFTs("djr0",tree,0,"rwgt_ctG_250.0",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list1 = makeTH1EFTs("djr1",tree,1,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list2 = makeTH1EFTs("djr2",tree,2,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list3 = makeTH1EFTs("djr3",tree,3,"",nbins,djrmin,djrmax);

    // Return here when just finding xsec
    //return;

    // Declare the dictionary (for reweighting) and incude some specific wc points
    std::map<string,std::pair<string,double> > rwgt_dict; 
    rwgt_dict["SM"] = std::make_pair("",0);
    rwgt_dict["RefPt"] = std::make_pair("",0);
    rwgt_dict["LitPt"] = std::make_pair("",0);
    rwgt_dict["arxiv1901-9wcPt"] = std::make_pair("",0);
    rwgt_dict["top19001hiPt"] = std::make_pair("",0);
    rwgt_dict["ttHJet30percPt"] = std::make_pair("",0);
    rwgt_dict["ttlnuJet30percPt"] = std::make_pair("",0);
    rwgt_dict["ttllJet30percPt"] = std::make_pair("",0);



    if ( IndividualWCScan ) {
        add_wcpts_to_dict("ctG",-4,4,2,rwgt_dict);
        add_wcpts_to_dict("cQq13",-4,4,2,rwgt_dict);
        add_wcpts_to_dict("cQq83",-4,4,2,rwgt_dict);
        add_wcpts_to_dict("cQq11",-15,15,5,rwgt_dict);
        add_wcpts_to_dict("cQq81",-15,15,5,rwgt_dict);
        add_wcpts_to_dict("ctq1",-15,15,5,rwgt_dict);
        add_wcpts_to_dict("ctq8",-12,12,4,rwgt_dict);
    }

    if ( WCscan ) {
        std::vector<string> wc_names_vect{ "ctW", "ctp", "cpQM", "ctei", "ctli", "cQei", "ctZ", "cQlMi", "cQl3i", "ctG", "ctlTi", "cbW", "cpQ3", "cptb", "cpt", "ctlSi" , "cQq13", "cQq83", "cQq11", "ctq1" , "cQq81", "ctq8"};
        std::vector<double> wc_vals_vect{   -1.8, -60  , -3.5  , -5.0  , 5.0   , -5.0  , 4.0  , 5.0    , -10.0  , 0.4  , 1.0    , 3.1  , 5.8   , -27.0 , 18   , -7.0    , 1.3    , 1.6    , 7.4    , 7.5    , 7.8    , 4.1   };
        //std::vector<double> ref_pt_vals_vect{ -8.303849, 64.337172, 45.883907, 24.328689, 24.43011, 23.757944, -6.093077, 23.951426, 21.540499, -3.609446, 21.809598, 49.595354, -51.106621, 136.133729, -43.552406, -20.005026}; // These are the values of the WC at the ref point
        //std::vector<double> ref_pt_vals_vect{-4.0 , 41.0 , 29.0 , 7.0 , 8.0 , 7.0 , -4.0 , 7.0 , -8.0 , -2.0 , -2.0 , -5.0 , -10.0 , -18.0 , -25.0 , -9.0}; // These numbers are from table 18 of the AN
        for (int i = 0; i < wc_names_vect.size(); i++) {
            //cout << wc_names_vect.at(i) << " " << ref_pt_vals_vect.at(i) << endl;
            if (basePtSM ) {
                string dict_key = wc_names_vect.at(i) + "_" + std::to_string(wc_vals_vect.at(i));
                rwgt_dict[dict_key] = std::make_pair(wc_names_vect.at(i),wc_vals_vect.at(i));
            } else if (basePtRefPt ) {
                string dict_key = wc_names_vect.at(i) + "_0";
                rwgt_dict[dict_key] = std::make_pair(wc_names_vect.at(i),0);
            }
        }
    }

    // Print the dictionary: 
    std::cout << "\nRwgt dictionary:" << std::endl;
    for (auto it =rwgt_dict.begin(); it !=rwgt_dict.end(); it++ ) {
        std::pair<string,double> p = it->second;
        std::cout << "    " << it->first << ": " << p.first << " " << p.second << std::endl ;
    }
    std::cout << "\n" << std::endl;

    // Set the WCPoint points: 
    string sm_pt_str = "";
    string ref_pt_str = "rwgt_ctW_-8.303849_ctp_64.337172_cpQM_45.883907_ctei_24.328689_ctli_24.43011_cQei_23.757944_ctZ_-6.093077_cQlMi_23.951426_cQl3i_21.540499_ctG_-3.609446_ctlTi_21.809598_cbW_49.595354_cpQ3_-51.106621_cptb_136.133729_cpt_-43.552406_ctlSi_-20.005026";

    string top19001_hi_str = "rwgt_ctp_44.26_cpQM_21.65_ctW_2.87_ctZ_3.15_ctG_1.18_cbW_4.95_cpQ3_3.48_cptb_12.63_cpt_12.31_cQl3i_8.97_cQlMi_4.99_cQei_4.59_ctli_4.82_ctei_4.86_ctlSi_6.52_ctlTi_0.84";
    std::string ttHJet_30perc_str   = "rwgt_cQei_100.0_cQl3i_100.0_cQlMi_100.0_cQq11_-1.25_cQq13_1.21_cQq81_2.17_cQq83_2.82_cbW_12.71_cpQ3_10.29_cpQM_-19.97_cpt_22.64_cptb_28.63_ctG_0.24_ctW_1.81_ctZ_-2.23_ctei_100.0_ctlSi_100.0_ctlTi_100.0_ctli_100.0_ctp_-2.29_ctq1_1.24_ctq8_2.03";
    std::string ttlnuJet_30perc_str = "rwgt_cQei_100.0_cQl3i_32.13_cQlMi_56.97_cQq11_-0.66_cQq13_0.44_cQq81_0.98_cQq83_-1.0_cbW_100.0_cpQ3_7.13_cpQM_-12.98_cpt_16.12_cptb_100.0_ctG_0.94_ctW_2.14_ctZ_-12.66_ctei_100.0_ctlSi_100.0_ctlTi_20.18_ctli_100.0_ctp_-89.0_ctq1_0.67_ctq8_0.89";
    std::string ttllJet_30perc_str  = "rwgt_cQei_-12.81_cQl3i_6.51_cQlMi_-10.06_cQq11_-0.94_cQq13_0.92_cQq81_1.7_cQq83_-2.19_cbW_12.09_cpQ3_15.15_cpQM_-2.86_cpt_3.94_cptb_42.52_ctG_0.57_ctW_2.6_ctZ_1.64_ctei_-12.88_ctlSi_-17.55_ctlTi_2.46_ctli_-9.74_ctp_100.0_ctq1_1.13_ctq8_1.87";

    // Construct the point from limits from the literature
    std::string two_heavy_lims = "ctG_0.4_ctW_-1.8_cbW_3.1_ctZ_4.0_cptb_-27_cpQ3_5.8_cpQM_-3.5_cpt_18_ctp_-60";     // From arxiv 1901
    std::string two_quark_two_lep_lims = "ctei_-5.0_ctli_5.0_cQei_-5.0_cQlMi_5.0_cQl3i_-10.0_ctlTi_1.0_ctlSi_-7.0"; // From TOP-19-001
    std::string two_heavy_two_light_lims = "cQq13_1.3_cQq83_1.6_cQq11_7.4_cQq81_7.8_ctq1_7.5_ctq8_4.1";             // From arxiv 1901
    std::string lit_str = "rwgt_" + two_heavy_lims + "_" + two_heavy_two_light_lims + "_" + two_quark_two_lep_lims;

    WCPoint* sm_pt = new WCPoint(sm_pt_str);
    WCPoint* ref_pt = new WCPoint(ref_pt_str); 
    WCPoint* lit_pt = new WCPoint(lit_str);
    WCPoint* arxiv1901_9wc_pt   = new WCPoint("rwgt_"+two_heavy_lims);
    WCPoint* top19001hi_pt      = new WCPoint(top19001_hi_str);
    WCPoint* ttHJet_30perc_pt   = new WCPoint(ttHJet_30perc_str);
    WCPoint* ttlnuJet_30perc_pt = new WCPoint(ttlnuJet_30perc_str);
    WCPoint* ttllJet_30perc_pt  = new WCPoint(ttllJet_30perc_str);
    WCPoint* tmp_pt;

    // Set the base point: This is the point our scan changes values **with respect to** (either SM point, or ref point)
    WCPoint* base_pt;
    if ( basePtSM ){
        base_pt = sm_pt;
    } else if ( basePtRefPt ){
        base_pt = ref_pt;
    } 

    // Loop over the dictionary and rwgt the hists: 
    for (auto i = rwgt_dict.begin(); i != rwgt_dict.end(); i++ ) {

        // Append rwgt info to png name:
        TString outfile_TString = (class TString)outfile;
        TString outfile_TString_pdf = (class TString)outfile;
        TString rwgt_string_key = i->first;
        outfile_TString.Append(rwgt_string_key);
        outfile_TString_pdf.Append(rwgt_string_key);
        if ( basePtSM ){
            outfile_TString.Append("_basePtSM");
            outfile_TString_pdf.Append("_basePtSM");
        } else if (basePtRefPt ){
            outfile_TString.Append("_basePtRefPt");
            outfile_TString_pdf.Append("_basePtRefPt");
        }
        outfile_TString.Append(".png");
        outfile_TString_pdf.Append(".pdf");

        tmp_pt = base_pt;
        //double wc_val = 0;
        double orig_val;

        std::pair<string,double> rwgt_pair = i->second; // rwgt_pair = (wc name, wc val), note wc name is empty for the 22d points

        // Set the rwgt pt that gets passed to makeplot:
        if (rwgt_string_key == "SM") {
            orig_val = sm_pt->getStrength(rwgt_pair.first);
            tmp_pt = sm_pt;
            std::cout << "    SM point" << std::endl;
        }

        else if (rwgt_string_key == "RefPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = ref_pt;
            std::cout << "    Ref point" << std::endl;

        } else if (rwgt_string_key == "LitPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = lit_pt;
            std::cout << "    Point with all 22 WCs set to values from the lit" << std::endl;

        } else if (rwgt_string_key == "arxiv1901-9wcPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = arxiv1901_9wc_pt;
            std::cout << "    Point used for pheno paper" << std::endl;

        } else if (rwgt_string_key == "top19001hiPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = top19001hi_pt;
            std::cout << "    TOP19001hi point" << std::endl;

        } else if (rwgt_string_key == "ttHJet30percPt"){
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = ttHJet_30perc_pt;
            std::cout << "    ttHJet30percPt point" << std::endl;

        } else if (rwgt_string_key == "ttlnuJet30percPt"){
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = ttlnuJet_30perc_pt;
            std::cout << "    ttlnuJet30percPt point" << std::endl;

        } else if (rwgt_string_key == "ttllJet30percPt"){
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = ttllJet_30perc_pt;
            std::cout << "    ttllJet30percPt point" << std::endl;

        } else {
            orig_val = base_pt->getStrength(rwgt_pair.first); 
            tmp_pt = base_pt;
            tmp_pt->setStrength(rwgt_pair.first,rwgt_pair.second);       
            std::cout << "    " << rwgt_string_key << ", " "rwgt point: " << rwgt_pair.first << ", " << rwgt_pair.second << std::endl;
        }

        //tmp_pt->setStrength(rwgt_pair.first,wc_val); // Set the base object to a new value

        std::cout << "\tOrig: " << orig_val << std::endl;
        std::cout << "\tNew: " << tmp_pt->getStrength(rwgt_pair.first) << std::endl;

        //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
        TCanvas *c1 = new TCanvas("c1", "c1", 800, 700);

        // For just the first plot (0->1): 
        c1->SetLeftMargin(0.0);
        c1->SetTopMargin(0.00);
        c1->SetRightMargin(0.00);
        c1->SetBottomMargin(0.00);
        TPad* pad[1];
        pad[0] = new TPad("pad0","pad",0,0,1.0,1.0);
        pad[0]->Draw();
        pad[0]->cd();
        pad[0]->SetBottomMargin(.17);
        //makeplot(tmp_pt,"DJR 0->1",hist_list0);
        makeplot(tmp_pt,"DJR 0#rightarrow1",hist_list0);

        // For all 4 plots: 
        //TPad *pad[4];
        //setcanvas(c1,pad);
        //pad[0]->cd();
        //makeplot(tmp_pt,"DJR 0->1",hist_list0);
        //pad[1]->cd();
        //makeplot(tmp_pt,"DJR 1->2",hist_list1);
        //pad[2]->cd();
        //makeplot(tmp_pt,"DJR 2->3",hist_list2);
        //pad[3]->cd();
        //makeplot(tmp_pt,"DJR 3->4",hist_list3);

        c1->Print(outfile_TString);
        c1->Print(outfile_TString_pdf);
        delete c1;

        // Reset the base object to the original state:
        tmp_pt->setStrength(rwgt_pair.first,orig_val);

        std::cout << "\tReset: " << "tmp: " << tmp_pt->getStrength(rwgt_pair.first) << ", " << "base: " << base_pt->getStrength(rwgt_pair.first) << std::endl;
    }

    return;

}


