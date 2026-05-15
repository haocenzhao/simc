// select_accepted_DIS_dat.C
// Usage:
//   root -l -b -q select_accepted_DIS_dat.C
//
// Function:
//   1) Read accepted event_id from ../worksim/dis_hms_e.root
//      and ../worksim/dis_shms_e.root
//   2) Select corresponding events from the original DIS .dat file
//   3) Write new .dat files with exactly the same original line format
//   4) Make kinematic plots and save them into two PDF files
//   5) For Q2, |q| and W, overlay the values stored in the accepted ROOT file
//
// Output:
//   HMS_10.2GeV_DIS_recp_select.dat
//   SHMS_10.2GeV_DIS_recp_select.dat
//   HMS_E10.2GeV.pdf
//   SHMS_E10.2GeV.pdf

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLegend.h>
#include <TColor.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

static bool ReadAcceptedIDs(const std::string& acceptedRootFile,
                            const std::string& treeName,
                            const std::string& eventIDBranch,
                            std::set<Long64_t>& acceptedIDs) {
  acceptedIDs.clear();

  TFile fAcc(acceptedRootFile.c_str(), "READ");
  if (fAcc.IsZombie()) {
    std::cerr << "ERROR: cannot open accepted ROOT file: "
              << acceptedRootFile << "\n";
    return false;
  }

  TTree* h10 = (TTree*)fAcc.Get(treeName.c_str());
  if (!h10) {
    std::cerr << "ERROR: cannot find TTree '" << treeName
              << "' in file: " << acceptedRootFile << "\n";
    return false;
  }

  Float_t event_id_f = -1.0;

  h10->SetBranchStatus("*", 0);
  h10->SetBranchStatus(eventIDBranch.c_str(), 1);
  h10->SetBranchAddress(eventIDBranch.c_str(), &event_id_f);

  const Long64_t NaccTree = h10->GetEntries();

  for (Long64_t i = 0; i < NaccTree; i++) {
    h10->GetEntry(i);
    Long64_t eid = (Long64_t)std::llround(event_id_f);
    acceptedIDs.insert(eid);
  }

  std::cout << "\nRead " << NaccTree
            << " entries from " << acceptedRootFile << "\n";
  std::cout << "Number of unique accepted event_id = "
            << acceptedIDs.size() << "\n";

  if (acceptedIDs.empty()) {
    std::cerr << "ERROR: accepted event_id set is empty for "
              << acceptedRootFile << "\n";
    return false;
  }

  return true;
}

static bool FillRootKinematicsHistograms(const std::string& acceptedRootFile,
                                         const std::string& treeName,
                                         TH1D* h_Q2_root,
                                         TH1D* h_q_root,
                                         TH1D* h_W_root) {
  TFile fAcc(acceptedRootFile.c_str(), "READ");
  if (fAcc.IsZombie()) {
    std::cerr << "ERROR: cannot open accepted ROOT file for kinematics: "
              << acceptedRootFile << "\n";
    return false;
  }

  TTree* h10 = (TTree*)fAcc.Get(treeName.c_str());
  if (!h10) {
    std::cerr << "ERROR: cannot find TTree '" << treeName
              << "' in file: " << acceptedRootFile << "\n";
    return false;
  }

  Float_t Q2_root = 0.0;
  Float_t q_root  = 0.0;
  Float_t W_root  = 0.0;

  h10->SetBranchStatus("*", 0);
  h10->SetBranchStatus("Q2", 1);
  h10->SetBranchStatus("q",  1);
  h10->SetBranchStatus("W",  1);

  h10->SetBranchAddress("Q2", &Q2_root);
  h10->SetBranchAddress("q",  &q_root);
  h10->SetBranchAddress("W",  &W_root);

  const Long64_t N = h10->GetEntries();

  for (Long64_t i = 0; i < N; i++) {
    h10->GetEntry(i);

    h_Q2_root->Fill(Q2_root);
    h_q_root->Fill(q_root);

    if (W_root > 0.0) {
      h_W_root->Fill(W_root);
    }
  }

  std::cout << "Filled ROOT-file Q2/q/W histograms from "
            << N << " accepted entries in " << acceptedRootFile << "\n";

  return true;
}

static void SetCanvasMargin(TCanvas* c, bool hasColorBar = false) {
  c->SetLeftMargin(0.12);
  c->SetRightMargin(hasColorBar ? 0.14 : 0.05);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);
}

static void Draw1DToPDF(TH1D* h,
                        const std::string& figPdf,
                        const std::string& canvasName,
                        const std::string& canvasTitle) {
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 650);
  SetCanvasMargin(c, false);
  h->SetLineWidth(2);
  h->SetFillStyle(0);
  h->Draw("HIST");
  c->Print(figPdf.c_str());
}

static void DrawOverlay1DToPDF(TH1D* h_dat,
                               TH1D* h_root,
                               const std::string& figPdf,
                               const std::string& canvasName,
                               const std::string& canvasTitle,
                               const std::string& xTitle) {
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 650);
  SetCanvasMargin(c, false);

  h_dat->SetLineColor(kBlue);
  h_dat->SetLineWidth(2);
  h_dat->SetFillStyle(0);

  h_root->SetLineColor(kRed);
  h_root->SetLineWidth(2);
  h_root->SetFillStyle(0);

  h_dat->GetXaxis()->SetTitle(xTitle.c_str());
  h_root->GetXaxis()->SetTitle(xTitle.c_str());

  const double max_dat  = h_dat->GetMaximum();
  const double max_root = h_root->GetMaximum();
  const double ymax = std::max(max_dat, max_root);

  h_dat->SetMaximum(ymax * 1.20);
  h_dat->SetMinimum(0.0);

  h_dat->Draw("HIST");
  h_root->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_dat,  "calculated from selected dat", "l");
  leg->AddEntry(h_root, "stored in accepted ROOT", "l");
  leg->Draw();

  c->Print(figPdf.c_str());
}

static void Draw2DToPDF(TH2D* h,
                        const std::string& figPdf,
                        const std::string& canvasName,
                        const std::string& canvasTitle) {
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasTitle.c_str(), 800, 650);
  SetCanvasMargin(c, true);
  h->Draw("COLZ");
  c->Print(figPdf.c_str());
}

static void ProcessOneSpectrometer(const std::string& specName,
                                   const std::string& acceptedRootFile,
                                   const std::string& inDatFile,
                                   const std::string& outDatFile,
                                   const std::string& figPdf,
                                   const std::string& treeName,
                                   const std::string& eventIDBranch) {
  const double Ebeam_MeV = 10200.0;
  const double Mp_GeV = 0.9382720813;

  std::set<Long64_t> acceptedIDs;

  if (!ReadAcceptedIDs(acceptedRootFile, treeName, eventIDBranch, acceptedIDs)) {
    return;
  }

  std::ifstream in(inDatFile.c_str());
  if (!in) {
    std::cerr << "ERROR: cannot open input dat file: "
              << inDatFile << "\n";
    return;
  }

  std::ofstream out(outDatFile.c_str());
  if (!out) {
    std::cerr << "ERROR: cannot open output dat file: "
              << outDatFile << "\n";
    return;
  }

  std::map<int, Long64_t> pidCount;

  const std::string tag = specName + " accepted";

  // ============================================================
  // Histograms
  // ============================================================

  TH2D* h_theta_phi_rec = new TH2D(
      Form("h_theta_phi_rec_%s", specName.c_str()),
      Form("%s recoil particle #theta-#phi;#phi [deg];#theta [deg]", tag.c_str()),
      180, -90.0, 90.0,
      120,  80.0, 200.0);

  TH1D* h_p_rec = new TH1D(
      Form("h_p_rec_%s", specName.c_str()),
      Form("%s recoil particle momentum magnitude;p [MeV/c];Counts", tag.c_str()),
      100, 0.0, 1000.0);

  TH2D* h_pz_vs_p_rec = new TH2D(
      Form("h_pz_vs_p_rec_%s", specName.c_str()),
      Form("%s recoil particle p_{z} vs p_{#perp};p_{#perp} [MeV/c];p_{z} [MeV/c]", tag.c_str()),
      110, 0.0, 1000.0,
      100, -800.0, 800.0);

  TH2D* h_p_vs_theta_rec = new TH2D(
      Form("h_p_vs_theta_rec_%s", specName.c_str()),
      Form("%s recoil particle p vs #theta;#theta [deg];p [MeV/c]", tag.c_str()),
      120, 80.0, 200.0,
      100, 0.0, 1000.0);

  TH2D* h_theta_phi_e = new TH2D(
      Form("h_theta_phi_e_%s", specName.c_str()),
      Form("%s scattered electron #theta-#phi;#phi [deg];#theta [deg]", tag.c_str()),
      190, -190.0, 190.0,
      100,   0.0, 200.0);

  TH1D* h_p_e = new TH1D(
      Form("h_p_e_%s", specName.c_str()),
      Form("%s scattered electron momentum magnitude;p [MeV/c];Counts", tag.c_str()),
      120, -1000.0, 11000.0);

  TH2D* h_pz_vs_p_e = new TH2D(
      Form("h_pz_vs_p_e_%s", specName.c_str()),
      Form("%s scattered electron p_{z} vs p_{#perp};p_{#perp} [MeV/c];p_{z} [MeV/c]", tag.c_str()),
      125, 0.0, 2500.0,
      200, -100.0, 10000.0);

  TH2D* h_p_vs_theta_e = new TH2D(
      Form("h_p_vs_theta_e_%s", specName.c_str()),
      Form("%s scattered electron p vs #theta;#theta [deg];p [MeV/c]", tag.c_str()),
      85, 0.0, 190.0,
      100, 0.0, 10000.0);

  // Q2/q/W calculated from selected dat
  TH1D* h_Q2 = new TH1D(
      Form("h_Q2_%s", specName.c_str()),
      Form("%s Q^{2} comparison;Q^{2} [GeV^{2}];Counts", tag.c_str()),
      70, 1.0, 8.0);

  TH1D* h_q = new TH1D(
      Form("h_q_%s", specName.c_str()),
      Form("%s |#vec{q}| comparison;|#vec{q}| [GeV/c];Counts", tag.c_str()),
      60, 4.0, 10.0);

  TH1D* h_W = new TH1D(
      Form("h_W_%s", specName.c_str()),
      Form("%s W comparison;W [GeV];Counts", tag.c_str()),
      50, 1.0, 6.0);

  // Q2/q/W directly read from accepted ROOT file
  TH1D* h_Q2_root = new TH1D(
      Form("h_Q2_root_%s", specName.c_str()),
      Form("%s Q^{2} comparison;Q^{2} [GeV^{2}];Counts", tag.c_str()),
      70, 1.0, 8.0);

  TH1D* h_q_root = new TH1D(
      Form("h_q_root_%s", specName.c_str()),
      Form("%s |#vec{q}| comparison;|#vec{q}| [GeV/c];Counts", tag.c_str()),
      60, 4.0, 10.0);

  TH1D* h_W_root = new TH1D(
      Form("h_W_root_%s", specName.c_str()),
      Form("%s W comparison;W [GeV];Counts", tag.c_str()),
      50, 1.0, 6.0);

  FillRootKinematicsHistograms(acceptedRootFile,
                               treeName,
                               h_Q2_root,
                               h_q_root,
                               h_W_root);

  // ============================================================
  // Read dat file
  // ============================================================

  std::string line;
  Long64_t nInputDatEvents = 0;
  Long64_t nSelectedEvents = 0;

  while (std::getline(in, line)) {
    if (line.empty()) continue;

    if (line[0] == '#') {
      out << line << "\n";
      continue;
    }

    std::istringstream iss(line);

    Long64_t evtID = -1;
    double pe0 = 0.0, pe1 = 0.0, pe2 = 0.0;
    double pr0 = 0.0, pr1 = 0.0, pr2 = 0.0;
    int pidrec = 0;
    double zcm = 0.0;
    double weight = 1.0;

    iss >> evtID
        >> pe0 >> pe1 >> pe2
        >> pr0 >> pr1 >> pr2
        >> pidrec >> zcm >> weight;

    if (!iss) {
      std::cerr << "WARNING: skip malformed line:\n"
                << line << "\n";
      continue;
    }

    nInputDatEvents++;

    if (acceptedIDs.find(evtID) == acceptedIDs.end()) {
      continue;
    }

    // Keep exactly the same original dat line format.
    out << line << "\n";

    nSelectedEvents++;
    pidCount[pidrec]++;

    // ==========================================================
    // Recoil-particle kinematics
    // ==========================================================

    const double p_rec  = std::sqrt(pr0*pr0 + pr1*pr1 + pr2*pr2);
    const double pt_rec = std::sqrt(pr0*pr0 + pr1*pr1);

    double theta_rec_deg = 0.0;
    if (p_rec > 0.0) {
      double cth = pr2 / p_rec;
      if (cth >  1.0) cth =  1.0;
      if (cth < -1.0) cth = -1.0;
      theta_rec_deg = std::acos(cth) * 180.0 / TMath::Pi();
    }

    const double phi_rec_deg =
        std::atan2(pr1, pr0) * 180.0 / TMath::Pi();

    h_theta_phi_rec->Fill(phi_rec_deg, theta_rec_deg);
    h_p_rec->Fill(p_rec);
    h_pz_vs_p_rec->Fill(pt_rec, pr2);
    h_p_vs_theta_rec->Fill(theta_rec_deg, p_rec);

    // ==========================================================
    // Scattered-electron kinematics
    // ==========================================================

    const double p_e  = std::sqrt(pe0*pe0 + pe1*pe1 + pe2*pe2);
    const double pt_e = std::sqrt(pe0*pe0 + pe1*pe1);

    double theta_e_deg = 0.0;
    if (p_e > 0.0) {
      double cth = pe2 / p_e;
      if (cth >  1.0) cth =  1.0;
      if (cth < -1.0) cth = -1.0;
      theta_e_deg = std::acos(cth) * 180.0 / TMath::Pi();
    }

    const double phi_e_deg =
        std::atan2(pe1, pe0) * 180.0 / TMath::Pi();

    h_theta_phi_e->Fill(phi_e_deg, theta_e_deg);
    h_p_e->Fill(p_e);
    h_pz_vs_p_e->Fill(pt_e, pe2);
    h_p_vs_theta_e->Fill(theta_e_deg, p_e);

    // ==========================================================
    // Inclusive DIS kinematics calculated from selected dat
    // ==========================================================

    const double Eeprime_MeV = p_e;
    const double nu_MeV = Ebeam_MeV - Eeprime_MeV;

    const double qx_MeV = -pe0;
    const double qy_MeV = -pe1;
    const double qz_MeV = Ebeam_MeV - pe2;

    const double q_mag_MeV =
        std::sqrt(qx_MeV*qx_MeV + qy_MeV*qy_MeV + qz_MeV*qz_MeV);

    const double Q2_MeV2 =
        q_mag_MeV*q_mag_MeV - nu_MeV*nu_MeV;

    const double Q2_GeV2 = Q2_MeV2 / 1.0e6;
    const double q_mag_GeV = q_mag_MeV / 1000.0;
    const double nu_GeV = nu_MeV / 1000.0;

    const double W2_GeV2 =
        Mp_GeV*Mp_GeV + 2.0*Mp_GeV*nu_GeV - Q2_GeV2;

    double W_GeV = -999.0;
    if (W2_GeV2 > 0.0) {
      W_GeV = std::sqrt(W2_GeV2);
    }

    h_Q2->Fill(Q2_GeV2);
    h_q->Fill(q_mag_GeV);

    if (W_GeV > 0.0) {
      h_W->Fill(W_GeV);
    }
  }

  std::cout << "\n========== " << specName << " ==========\n";
  std::cout << "Input dat events checked: " << nInputDatEvents << "\n";
  std::cout << "Selected accepted events from dat:  " << nSelectedEvents << "\n";
  std::cout << "Accepted unique event_id from ROOT: " << acceptedIDs.size() << "\n";
  std::cout << "Difference ROOT - selected dat:     "
            << (Long64_t)acceptedIDs.size() - nSelectedEvents << "\n";
  std::cout << "Wrote selected dat file:   " << outDatFile << "\n";
  
  if (nSelectedEvents == 0) {
    std::cerr << "WARNING: no matched events found for " << specName << ".\n"
              << "Please check whether event_id in h10 matches "
              << "the first column EvtID in the dat file.\n";
    return;
  }

  // ============================================================
  // pidrec plot
  // ============================================================

  std::vector<int> pdgs;
  pdgs.reserve(pidCount.size());
  for (auto &kv : pidCount) pdgs.push_back(kv.first);
  std::sort(pdgs.begin(), pdgs.end());

  TCanvas* c_pid = new TCanvas(
      Form("c_pidrec_%s", specName.c_str()),
      Form("%s pidrec", specName.c_str()),
      700, 500);

  c_pid->SetLeftMargin(0.12);
  c_pid->SetRightMargin(0.04);
  c_pid->SetBottomMargin(0.20);
  c_pid->SetTopMargin(0.08);

  TH1F* h_pid = new TH1F(
      Form("h_pidrec_%s", specName.c_str()),
      Form("%s pidrec distribution;particle species (PDG);counts", tag.c_str()),
      (int)pdgs.size(), 0.5, (double)pdgs.size() + 0.5);

  Long64_t ymax = 0;
  for (int i = 0; i < (int)pdgs.size(); i++) {
    const int pdg = pdgs[i];
    const Long64_t cnt = pidCount[pdg];

    h_pid->SetBinContent(i + 1, (double)cnt);

    if (pdg == 2212) {
      h_pid->GetXaxis()->SetBinLabel(i + 1, "proton (2212)");
    } else if (pdg == 2112) {
      h_pid->GetXaxis()->SetBinLabel(i + 1, "neutron (2112)");
    } else {
      h_pid->GetXaxis()->SetBinLabel(i + 1, Form("%d", pdg));
    }

    if (cnt > ymax) ymax = cnt;
  }

  h_pid->SetMinimum(0.0);
  h_pid->SetMaximum((double)ymax * 1.20);
  h_pid->GetXaxis()->LabelsOption("v");
  h_pid->GetXaxis()->SetLabelSize(0.055);
  h_pid->GetYaxis()->SetLabelSize(0.05);
  h_pid->GetXaxis()->SetTitleSize(0.055);
  h_pid->GetYaxis()->SetTitleSize(0.06);
  h_pid->SetLineWidth(2);
  h_pid->SetFillStyle(0);
  h_pid->Draw("HIST");

  TLatex lat;
  lat.SetTextAlign(21);
  lat.SetTextSize(0.04);
  for (int i = 1; i <= h_pid->GetNbinsX(); i++) {
    const double y = h_pid->GetBinContent(i);
    lat.DrawLatex(i, y + 0.02 * (double)ymax,
                  Form("%lld", (Long64_t)y));
  }

  // ============================================================
  // Multi-page PDF output
  // ============================================================

  c_pid->Print((figPdf + "[").c_str());
  c_pid->Print(figPdf.c_str());

  Draw2DToPDF(h_theta_phi_rec,
              figPdf,
              Form("c_theta_phi_rec_%s", specName.c_str()),
              Form("%s theta_phi_rec", specName.c_str()));

  Draw1DToPDF(h_p_rec,
              figPdf,
              Form("c_p_rec_%s", specName.c_str()),
              Form("%s p_magnitude_rec", specName.c_str()));

  Draw2DToPDF(h_pz_vs_p_rec,
              figPdf,
              Form("c_pz_vs_p_rec_%s", specName.c_str()),
              Form("%s pz_vs_p_perp_rec", specName.c_str()));

  Draw2DToPDF(h_p_vs_theta_rec,
              figPdf,
              Form("c_p_vs_theta_rec_%s", specName.c_str()),
              Form("%s p_vs_theta_rec", specName.c_str()));

  Draw2DToPDF(h_theta_phi_e,
              figPdf,
              Form("c_theta_phi_e_%s", specName.c_str()),
              Form("%s theta_phi_e", specName.c_str()));

  Draw1DToPDF(h_p_e,
              figPdf,
              Form("c_p_e_%s", specName.c_str()),
              Form("%s p_magnitude_e", specName.c_str()));

  Draw2DToPDF(h_pz_vs_p_e,
              figPdf,
              Form("c_pz_vs_p_e_%s", specName.c_str()),
              Form("%s pz_vs_p_perp_e", specName.c_str()));

  Draw2DToPDF(h_p_vs_theta_e,
              figPdf,
              Form("c_p_vs_theta_e_%s", specName.c_str()),
              Form("%s p_vs_theta_e", specName.c_str()));

  DrawOverlay1DToPDF(h_Q2,
                     h_Q2_root,
                     figPdf,
                     Form("c_Q2_compare_%s", specName.c_str()),
                     Form("%s Q2 comparison", specName.c_str()),
                     "Q^{2} [GeV^{2}]");

  DrawOverlay1DToPDF(h_q,
                     h_q_root,
                     figPdf,
                     Form("c_q_compare_%s", specName.c_str()),
                     Form("%s q comparison", specName.c_str()),
                     "|#vec{q}| [GeV/c]");

  TCanvas* c_W =
      new TCanvas(Form("c_W_compare_%s", specName.c_str()),
                  Form("%s W comparison", specName.c_str()), 800, 650);
  SetCanvasMargin(c_W, false);

  h_W->SetLineColor(kBlue);
  h_W->SetLineWidth(2);
  h_W->SetFillStyle(0);

  h_W_root->SetLineColor(kRed);
  h_W_root->SetLineWidth(2);
  h_W_root->SetFillStyle(0);

  const double ymax_W = std::max(h_W->GetMaximum(), h_W_root->GetMaximum());
  h_W->SetMaximum(ymax_W * 1.20);
  h_W->SetMinimum(0.0);

  h_W->Draw("HIST");
  h_W_root->Draw("HIST SAME");

  TLegend* leg_W = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg_W->SetBorderSize(0);
  leg_W->SetFillStyle(0);
  leg_W->AddEntry(h_W,      "calculated from selected dat", "l");
  leg_W->AddEntry(h_W_root, "stored in accepted ROOT", "l");
  leg_W->Draw();

  c_W->Print(figPdf.c_str());
  c_W->Print((figPdf + "]").c_str());

  std::cout << "\nSaved plots to: " << figPdf << "\n";

  std::cout << "\n" << specName << " accepted PIDREC counts:\n";
  for (auto &kv : pidCount) {
    std::cout << "  " << kv.first << " : " << kv.second << "\n";
  }

  std::cout << "\nDone with " << specName << ".\n";
}

void select_accepted_DIS_dat() {
  const std::string treeName = "h10";
  const std::string eventIDBranch = "event_id";

  const std::string inDatFile =
      "test_band_tagged_full_PRC_E10.2GeV_recp.dat";

  const std::string hmsRootFile  = "../worksim/dis_hms_e.root";
  const std::string shmsRootFile = "../worksim/dis_shms_e.root";

  const std::string hmsOutDatFile  = "HMS_10.2GeV_DIS_recp_select.dat";
  const std::string shmsOutDatFile = "SHMS_10.2GeV_DIS_recp_select.dat";

  const std::string hmsFigPdf  = "HMS_E10.2GeV.pdf";
  const std::string shmsFigPdf = "SHMS_E10.2GeV.pdf";

  gStyle->SetOptStat(0);

  ProcessOneSpectrometer("HMS",
                         hmsRootFile,
                         inDatFile,
                         hmsOutDatFile,
                         hmsFigPdf,
                         treeName,
                         eventIDBranch);

  ProcessOneSpectrometer("SHMS",
                         shmsRootFile,
                         inDatFile,
                         shmsOutDatFile,
                         shmsFigPdf,
                         treeName,
                         eventIDBranch);

  std::cout << "\nAll done.\n";
}