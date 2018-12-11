#Developed by L.Dawson, A.Minotti, F.Xie
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

void alphaLengthPDF(){

  using namespace RooFit;
  using namespace RooStats;


  string foils="/Users/fang/SuperNEMO/Analysis/RadonBKG/SimuRecData/0Gauss_1E6/Bulk/sensitivity.root";
  string surface="/Users/fang/SuperNEMO/Analysis/RadonBKG/SimuRecData/0Gauss_1E6/Surf/sensitivity.root";
  string wires="/Users/fang/SuperNEMO/Analysis/RadonBKG/SimuRecData/0Gauss_1E6/wire/sensitivity.root";

  TFile *file1=new TFile(foils.c_str());
  TFile *file2=new TFile(surface.c_str());
  TFile *file3=new TFile(wires.c_str());

  TTree *treeFoil=(TTree*)file1->Get("Sensitivity");
  TTree *treeSurf=(TTree*)file2->Get("Sensitivity");
  TTree *treeWire=(TTree*)file3->Get("Sensitivity");

  int sampleBulk=treeFoil->GetEntries();
  int sampleSurf=treeSurf->GetEntries();
  int sampleWire=treeWire->GetEntries();

  gStyle->SetOptStat(1111);

  TH1F *h1 = new TH1F("h1", "h1", 100, 0, 500);
  TH1F *h2 = new TH1F("h2", "h2", 100, 0, 500);
  TH1F *h3 = new TH1F("h3", "h3", 100, 0, 500);
  TH1F *h_total = new TH1F("h_total", "h_total", 100, 0, 500);

  h1->SetTitle("Alpha track lengths Bi214, 3E5 1e1a field wires,10 #mus, xy 40 cm");
  h1->GetYaxis()->SetTitle("Events s^{-1} / 5mm"); //("Number of events");
  h1->GetXaxis()->SetTitle("Alpha length (mm)");
  h_total->SetTitle("Reference activity - tracker selection, 25 Gauss");
  h_total->GetYaxis()->SetTitle("Events s^{-1} / 5mm");
  h_total->GetXaxis()->SetTitle("Alpha length (mm)");

  h1->SetLineColor(kRed);
  h1->SetFillColor(kRed);
  h1->SetFillStyle(3003);
  h2->SetLineColor(kGreen+2);
  h2->SetFillColor(kGreen+2);
  h2->SetFillStyle(3003);
  h3->SetLineColor(kBlue);
  h3->SetFillColor(kBlue);
  h3->SetFillStyle(3003);
  h_total->SetLineColor(kBlack);

  bool topology_1e1alpha=0;
  vector<bool> *alphas_from_foil=0;
  vector<bool> *electron_hits_main_wall=0;
  vector<int> *electron_charges=0;
  vector<bool> *electronsFromFoil=0;
  double alpha_track_length=0;

  const bool verbose = false;

  int exposure=0;
  int exposure_in_seconds=0;
  int number_of_entries_bulk=0;
  int number_of_entries_surf=0;
  int number_of_entries_wire=0;
  int number_exp_bulk=0;
  int number_exp_surf=0;
  int number_exp_wire=0;
  int total_number_exp=0;
  double efficiency_bulk=0;
  double efficiency_surf=0;
  double efficiency_wire=0;
  double old_fraction_bulk=0;
  double old_fraction_surf=0;
  double old_fraction_wire=0;
  int total_number_entries=0;

  double activity_tracker=2.28E-3; // 2.28E-3 with flushing
  double activity_bulk=15.4E-3;//70E-6; //70 µBq target, 10 µBq/kg
  double activity_surf=0.1685E-3;//activity_tracker*0.078; //0.178 mBq //3.37mBq/20
//  double activity_wire=42.3E-3;//activity_tracker*0.922; // 2.10 mBq
  double activity_wire=2.1E-3;//activity_tracker*0.922; // 2.10 mBq

  double h1_scale=0;
  double h2_scale=0;
  double h3_scale=0;

  treeFoil->SetBranchAddress("reco.topology_1e1alpha",&topology_1e1alpha);
  treeFoil->SetBranchAddress("reco.electron_hits_mainwall", &electron_hits_main_wall);
  treeFoil->SetBranchAddress("reco.electron_charges", &electron_charges);
  treeFoil->SetBranchAddress("reco.electrons_from_foil", &electronsFromFoil);
  treeFoil->SetBranchAddress("reco.alpha_track_length", &alpha_track_length);
  treeFoil->SetBranchAddress("reco.alphas_from_foil", &alphas_from_foil);
  treeSurf->SetBranchAddress("reco.topology_1e1alpha",&topology_1e1alpha);
  treeSurf->SetBranchAddress("reco.electron_hits_mainwall", &electron_hits_main_wall);
  treeSurf->SetBranchAddress("reco.electron_charges", &electron_charges);
  treeSurf->SetBranchAddress("reco.electrons_from_foil", &electronsFromFoil);
  treeSurf->SetBranchAddress("reco.alpha_track_length", &alpha_track_length);
  treeSurf->SetBranchAddress("reco.alphas_from_foil", &alphas_from_foil);
  treeWire->SetBranchAddress("reco.topology_1e1alpha",&topology_1e1alpha);
  treeWire->SetBranchAddress("reco.electron_hits_mainwall", &electron_hits_main_wall);
  treeWire->SetBranchAddress("reco.electron_charges", &electron_charges);
  treeWire->SetBranchAddress("reco.electrons_from_foil", &electronsFromFoil);
  treeWire->SetBranchAddress("reco.alpha_track_length", &alpha_track_length);
  treeWire->SetBranchAddress("reco.alphas_from_foil", &alphas_from_foil);

  for(int entry=0; entry<sampleBulk; entry++){
    treeFoil->GetEntry(entry);
    if(topology_1e1alpha){
      if(electronsFromFoil->at(0)==0 && alphas_from_foil->at(0)==0){
          h1->Fill(alpha_track_length);
      }
    }
  }



  for(int entry=0; entry<sampleSurf; entry++){
    treeSurf->GetEntry(entry);
    if(topology_1e1alpha){
      if(electronsFromFoil->at(0)==0 && alphas_from_foil->at(0)==0){
          h2->Fill(alpha_track_length);
      }
    }
  }
  for(int entry=0; entry<sampleWire; entry++){
    treeWire->GetEntry(entry);
    if(topology_1e1alpha){
      if(electronsFromFoil->at(0)==0 && alphas_from_foil->at(0)==0){
            h3->Fill(alpha_track_length);
      }
    }
  }

  //normalise the histograms, scaling by activity
  number_of_entries_bulk = h1->GetEntries();
  cout<<"Number of entries bulk: "<<number_of_entries_bulk<<endl;
  number_of_entries_surf = h2->GetEntries();
  cout<<"Number of entries surf: "<<number_of_entries_surf<<endl;
  number_of_entries_wire = h3->GetEntries();
  cout<<"Number of entries wire: "<<number_of_entries_wire<<endl;

  h1_scale = (activity_bulk/sampleBulk);
  h2_scale = (activity_surf/sampleSurf);
  h3_scale = (activity_wire/sampleWire);

  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();

  h1->Scale(h1_scale);
  h2->Scale(h2_scale);
  h3->Scale(h3_scale);

  h_total->Draw("hist");
  h3->Draw("hist same");
  h2->Draw("hist same");
  h1->Draw("hist same");

  TLegend *leg = new TLegend(0.1293878,0.5787037,0.295102,0.9074074);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(h_total,"Combined events","l");
  leg->AddEntry(h1,"{}^{214}Bi Bulk","f");
  leg->AddEntry(h2,"{}^{214}Bi Surface","f");
  leg->AddEntry(h3,"{}^{214}Bi Tracker","f");
  leg->Draw("same");

  exposure=80; // 60 days
  exposure_in_seconds=exposure*24*60*60;

  double expoure_sqrt;
  expoure_sqrt=sqrt(exposure);

  double exposure_rsd;
  exposure_rsd=expoure_sqrt/exposure;

  //Sum the normalised histograms of the individual generators
  h_total->Add(h1);
  h_total->Add(h2);
  h_total->Add(h3);

  //Scale by the exposure which here is 60 days
  h_total->Scale(exposure_in_seconds);
  h1->Scale(exposure_in_seconds);
  h2->Scale(exposure_in_seconds);
  h3->Scale(exposure_in_seconds);

  //Calculate the old fractions by finding the expected number of entries
  // N_exp = efficiency*activity*exposure
  total_number_entries=h_total->GetEntries();
  efficiency_bulk=(double(number_of_entries_bulk)/double(sampleBulk));
  cout<<"Bulk efficiency: "<<efficiency_bulk<<endl;
  efficiency_surf=(double(number_of_entries_surf)/double(sampleSurf));
  cout<<"Surface efficiency: "<<efficiency_surf<<endl;
  efficiency_wire=(double(number_of_entries_wire)/double(sampleWire));
  cout<<"Tracker efficiency: "<<efficiency_wire<<endl;

  number_exp_bulk=efficiency_bulk*activity_bulk*exposure_in_seconds;
  number_exp_surf=efficiency_surf*activity_surf*exposure_in_seconds;
  number_exp_wire=efficiency_wire*activity_wire*exposure_in_seconds;
  total_number_exp=number_exp_bulk+number_exp_surf+number_exp_wire;

  const int n_pseudo = 100;
  const unsigned int poly_order = 0;

  double fit_limit_bulk[2];
  double fit_limit_surf[2];
  double fit_limit_wire[2];

  const double fit_limit_coeffs[2] = {1.0-exposure_rsd,1.0+exposure_rsd};

//  int bin_number = (1-fit_limit_coeffs[0])*50;
  int bin_number= 50;

  fit_limit_bulk[0] = activity_bulk*fit_limit_coeffs[0];
  fit_limit_bulk[1] = activity_bulk*fit_limit_coeffs[1];
  fit_limit_surf[0] = activity_surf*fit_limit_coeffs[0];
  fit_limit_surf[1] = activity_surf*fit_limit_coeffs[1];
  fit_limit_wire[0] = activity_wire*fit_limit_coeffs[0];
  fit_limit_wire[1] = activity_wire*fit_limit_coeffs[1];

  cout<<"fit_limit_wire_lower_bound is: " << fit_limit_wire[0]<< endl;
  cout<<"fit_limit_wire_higher_bound is: " << fit_limit_wire[1]<< endl;
  cout<<"real activity wire is: " << activity_wire<< endl;

  TH1D *h_fitted_activity_bulk;
  TH1D *h_fitted_activity_surf;
  TH1D *h_fitted_activity_wire;

  h_fitted_activity_bulk = new TH1D("h_fitted_activity_bulk","h_fitted_activity_bulk",
      bin_number,fit_limit_bulk[0],fit_limit_bulk[1]);
  h_fitted_activity_surf = new TH1D("h_fitted_activity_surf","h_fitted_activity_surf",
      bin_number,fit_limit_surf[0],fit_limit_surf[1]);
  h_fitted_activity_wire = new TH1D("h_fitted_activity_wire","h_fitted_activity_wire",
      bin_number,fit_limit_wire[0],fit_limit_wire[1]);

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(1);

  RooWorkspace w("w");

  RooRealVar track_length_obs("track_length_obs","track_length_obs",0,500);

  RooRealVar *fitted_activity_bulk;
  RooRealVar *fitted_activity_surf;
  RooRealVar *fitted_activity_wire;

  fitted_activity_bulk = new RooRealVar("fitted_activity_bulk","fitted_activity_bulk",
      (exposure_in_seconds*efficiency_bulk)*fit_limit_bulk[0],
      (exposure_in_seconds*efficiency_bulk)*fit_limit_bulk[1]);

  fitted_activity_surf = new RooRealVar("fitted_activity_surf","fitted_activity_surf",
      (exposure_in_seconds*efficiency_surf)*fit_limit_surf[0],
      (exposure_in_seconds*efficiency_surf)*fit_limit_surf[1]);

  fitted_activity_wire = new RooRealVar("fitted_activity_wire","fitted_activity_wire",
      (exposure_in_seconds*efficiency_wire)*fit_limit_wire[0],
      (exposure_in_seconds*efficiency_wire)*fit_limit_wire[1]);

  RooDataHist datahist("datahist","datahist",track_length_obs, h_total);
  RooHistPdf pdf("pdf","pdf",track_length_obs,datahist,poly_order);

  RooDataHist *datahist_bulk;
  RooDataHist *datahist_surf;
  RooDataHist *datahist_wire;

  datahist_bulk = new RooDataHist("datahist_bulk","datahist_bulk",track_length_obs,h1);
  datahist_surf = new RooDataHist("datahist_surf","datahist_surf",track_length_obs,h2);
  datahist_wire = new RooDataHist("datahist_wire","datahist_wire",track_length_obs,h3);

  RooHistPdf *pdf_bulk;
  RooHistPdf *pdf_surf;
  RooHistPdf *pdf_wire;

  pdf_bulk = new RooHistPdf("pdf_bulk","pdf_bulk",track_length_obs, *datahist_bulk,poly_order);
  pdf_surf = new RooHistPdf("pdf_surf","pdf_surf",track_length_obs, *datahist_surf,poly_order);
  pdf_wire = new RooHistPdf("pdf_wire","pdf_wire",track_length_obs, *datahist_wire,poly_order);


  RooAddPdf model("model","model",RooArgSet(*pdf_bulk,*pdf_surf,*pdf_wire),
      RooArgList(*fitted_activity_bulk,*fitted_activity_surf,*fitted_activity_wire));
//  RooAddPdf model("model","model",RooArgSet(*pdf_wire),RooArgList(*fitted_activity_wire));

  TCanvas *c_mock[n_pseudo];
  TLegend *l_mock[n_pseudo];

  TH1 *h_total_data[n_pseudo];
  TH1 *h_total_mock[n_pseudo];
  TH1 *h1_mock[n_pseudo];
  TH1 *h2_mock[n_pseudo];
  TH1 *h3_mock[n_pseudo];

  for(int i_pseudo=0;i_pseudo<n_pseudo;i_pseudo++) {


  RooDataSet *data = pdf.generate(track_length_obs, h_total->Integral(), Extended(kTRUE));

  data->SetName("data");
  model.fitTo(*data,PrintLevel(-1));

    debugging test start
      RooPlot * pl4 = track_length_obs.frame();
      data->plotOn(pl4,MarkerStyle(kFullSquare));
      pl4->Draw("same");

      TLegend *legData = new TLegend(0.1293878,0.5787037,0.295102,0.9074074);
      legData->SetFillColor(0);
      legData->SetLineColor(0);
      legData->AddEntry(pl4,"mock data","EP");
      legData->Draw("same");
    


  h_fitted_activity_bulk->Fill(fitted_activity_bulk->getValV()/(exposure_in_seconds*efficiency_bulk) );
  h_fitted_activity_surf->Fill(fitted_activity_surf->getValV()/(exposure_in_seconds*efficiency_surf) );
  h_fitted_activity_wire->Fill(fitted_activity_wire->getValV()/(exposure_in_seconds*efficiency_wire) );

  cout<< "fitted_activity_wire is "<< fitted_activity_wire->getValV()/(exposure_in_seconds*efficiency_wire) << endl;

if(verbose) {

    c_mock[i_pseudo] = new TCanvas("c_mock","c_mock",1961,344,700,502);
    l_mock[i_pseudo] = new TLegend(0.5,0.65,0.89,0.89);

    h_total_data[i_pseudo] = data->createHistogram("track_length_obs",100);
    cout << "INTEGRAL: " << h_total_data[i_pseudo]->Integral() << endl;
    h_total_data[i_pseudo]->GetXaxis()->SetTitle("#alpha track length (mm)");
    h_total_data[i_pseudo]->SetLineWidth(2);
    h_total_data[i_pseudo]->SetLineColor(1);
    h_total_data[i_pseudo]->SetFillStyle(0);
    h_total_data[i_pseudo]->Draw();
    h_total_data[i_pseudo]->SetTitle("mock data (60 days exposure)");
    l_mock[i_pseudo]->AddEntry(h_total_data[i_pseudo],h_total_data[i_pseudo]->GetTitle(),"lep");
    h_total_data[i_pseudo]->SetTitle("");

    h_total_mock[i_pseudo] = model.createHistogram("track_length_obs",100);
    h_total_mock[i_pseudo]->Scale(5);
    h_total_mock[i_pseudo]->SetLineWidth(2);
    h_total_mock[i_pseudo]->SetLineColor(16);
    h_total_mock[i_pseudo]->SetFillColor(16);
    h_total_mock[i_pseudo]->SetFillStyle(3003);
    h_total_mock[i_pseudo]->Draw("histsames");
    h_total_mock[i_pseudo]->SetTitle("global PDF fitted to mock data");
    l_mock[i_pseudo]->AddEntry(h_total_mock[i_pseudo],h_total_mock[i_pseudo]->GetTitle(),"l");

    h3_mock[i_pseudo] = pdf_wire->createHistogram("track_length_obs",100);
    h3_mock[i_pseudo]->Scale(fitted_activity_wire->getValV());
    h3_mock[i_pseudo]->SetLineWidth(2);
    h3_mock[i_pseudo]->SetLineColor(kRed);
    h3_mock[i_pseudo]->SetFillColor(kRed);
    h3_mock[i_pseudo]->SetFillStyle(3003);
    h3_mock[i_pseudo]->Draw("histsames");
    h3_mock[i_pseudo]->SetTitle("pdf_wire_mock");
    l_mock[i_pseudo]->AddEntry(h3_mock[i_pseudo],h3_mock[i_pseudo]->GetTitle(),"l");
    l_mock[i_pseudo]->Draw();

    }
  }

 TCanvas *c_wire;

    c_wire = new TCanvas("c_wire","c_wire",1961,344,700,502);

    h_fitted_activity_wire->SetLineWidth(2);
    h_fitted_activity_wire->SetLineColor(kGray);
    h_fitted_activity_wire->Draw();

h_fitted_activity_wire->SetTitle("h_fitted_activity_wire0");

    TF1 *myfit = new TF1("myfit","gaus");
    
   h_fitted_activity_wire->Fit(myfit,"","",fit_limit_wire[0],fit_limit_wire[1]); // option "LE" give succeddful fit

    cout << "FITTED HISTOGRAM OF MOCK-DATA ACTIVITIES WITH GAUSSIAN FUNCTION" << endl;
    cout << "GAUSSIAN MEAN: " << myfit->GetParameter(1) << endl;
    cout << "GAUSSIAN SIGMA: " << myfit->GetParameter(2) << endl;
    cout << "ACCURACY OF ACTIVITY MEASUREMENT AFTER " << exposure << " DAYS: " << (myfit->GetParameter(2)/myfit->GetParameter(1))*100. << "%" <<endl;

   c_wire->Write();

 }
