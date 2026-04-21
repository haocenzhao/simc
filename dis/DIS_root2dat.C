// DIS_root2dat.C
// Usage: root -l -b -q DIS_root2dat.C
//
// ROOT branches units (given):
//   pe[3]   : GeV/c
//   prec[3] : GeV/c
//   zvtx    : cm
//   pidrec  : PDG code (unitless)
//   weight  : unitless (event weight)
//
// Output is converted to SIMC units (simulate.inc):
//   momentum: MeV/c
//   length  : cm
//
// Added features:
//   1) theta-phi TH2F for recoil particle prec and scattered electron pe
//   2) momentum magnitude TH1F for prec and pe
//   3) pz vs p_perp TH2F for prec and pe
//   4) p vs theta TH2F for prec and pe
//   5) all plots saved into the existing pdf as a multi-page PDF

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

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

static std::string BaseNoExt(const std::string& path) {
  std::string name = gSystem->BaseName(path.c_str());  // drop directory
  auto pos = name.find_last_of('.');
  if (pos != std::string::npos) name = name.substr(0, pos);
  return name;
}

void DIS_root2dat() {
  const std::string inFile =
    "/w/hallc-scshelf2102/c-lad/haocen/Deuteron_DIS/GEMC/generator/10.2/tagged_full/PRC/test_band_tagged_full_PRC_E10.2GeV_recp.root";
  const std::string base    = BaseNoExt(inFile);
  const std::string outFile = base + ".dat";              // written to CWD
  const std::string figPng  = base + "_pidrec.png";
  const std::string figPdf  = base + "_pidrec.pdf";

  // ============================================================
  // Unit conversion knobs: set to 1 if already in SIMC units
  //   ROOT:  pe/prec in GeV/c  -> SIMC: MeV/c  => pscale = 1000
  //   ROOT:  zvtx in cm        -> SIMC: cm     => xscale = 1
  // ============================================================
  const double pscale = 1000.0;  // (GeV/c -> MeV/c). set 1.0 if already MeV/c
  const double xscale = 1.0;     // (cm -> cm).     set 1.0 if already cm

  TFile f(inFile.c_str(), "READ");
  if (f.IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inFile << "\n";
    return;
  }

  TTree* T = (TTree*)f.Get("T");
  if (!T) {
    std::cerr << "ERROR: cannot find TTree 'T' in file: " << inFile << "\n";
    return;
  }

  // Branch variables (ROOT units)
  double pe[3]   = {0, 0, 0};
  double prec[3] = {0, 0, 0};
  int    pidrec  = 0;
  double zvtx    = 0.0;
  double weight  = 1.0;

  // Enable only needed branches
  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("pe", 1);
  T->SetBranchStatus("prec", 1);
  T->SetBranchStatus("pidrec", 1);
  T->SetBranchStatus("zvtx", 1);
  T->SetBranchStatus("weight", 1);

  // Bind branches
  T->SetBranchAddress("pe", pe);
  T->SetBranchAddress("prec", prec);
  T->SetBranchAddress("pidrec", &pidrec);
  T->SetBranchAddress("zvtx", &zvtx);
  T->SetBranchAddress("weight", &weight);

  std::ofstream out(outFile.c_str());
  if (!out) {
    std::cerr << "ERROR: cannot open output file: " << outFile << "\n";
    return;
  }
  out.setf(std::ios::scientific);
  out.precision(10);

  // Header with explicit output units (SIMC)
  out << "# EvtID "
      << "pe_px[MeV/c] pe_py[MeV/c] pe_pz[MeV/c] "
      << "prec_px[MeV/c] prec_py[MeV/c] prec_pz[MeV/c] "
      << "pidrec zvtx[cm] weight\n";

  // Count pidrec distribution
  std::map<int, Long64_t> pidCount;

  // ------------------------------------------------------------
  // Histograms for recoil particle prec
  // ------------------------------------------------------------
  gStyle->SetOptStat(0);

  TH2D* h_theta_phi_rec = new TH2D(
      "h_theta_phi_rec",
      "Recoil particle #theta-#phi;#phi [deg];#theta [deg]",
      180, -90.0, 90.0,
      120,  80.0, 200.0);

  TH1D* h_p_rec = new TH1D(
      "h_p_rec",
      "Recoil particle momentum magnitude;p [MeV/c];Counts",
      100, 0.0, 1000.0);

  TH2D* h_pz_vs_p_rec = new TH2D(
      "h_pz_vs_p_rec",
      "Recoil particle p_{z} vs p_{#perp};p_{#perp} [MeV/c];p_{z} [MeV/c]",
      110, -100.0, 1000.0,
      100, -800.0, 800.0);

  TH2D* h_p_vs_theta_rec = new TH2D(
      "h_p_vs_theta_rec",
      "Recoil particle p vs #theta;#theta [deg];p [MeV/c]",
      120, 80.0, 200.0,
      100, 0.0, 1000.0);

  // ------------------------------------------------------------
  // Histograms for scattered electron pe
  // ------------------------------------------------------------
  TH2D* h_theta_phi_e = new TH2D(
      "h_theta_phi_e",
      "Scattered electron #theta-#phi;#phi [deg];#theta [deg]",
      190, -190.0, 190.0,
      100,   0.0, 200.0);

  TH1D* h_p_e = new TH1D(
      "h_p_e",
      "Scattered electron momentum magnitude;p [MeV/c];Counts",
      120, -1000.0, 11000.0);

  TH2D* h_pz_vs_p_e = new TH2D(
      "h_pz_vs_p_e",
      "Scattered electron p_{z} vs p_{#perp};p_{#perp} [MeV/c];p_{z} [MeV/c]",
      125, 0.0, 2500.0,
      200, -100.0, 10000.0);

  TH2D* h_p_vs_theta_e = new TH2D(
      "h_p_vs_theta_e",
      "Scattered electron p vs #theta;#theta [deg];p [MeV/c]",
      85, 0.0, 190.0,
      100, 0.0, 10000.0);

  const Long64_t N = T->GetEntries();
  for (Long64_t i = 0; i < N; i++) {
    T->GetEntry(i);

    pidCount[pidrec]++;

    // Convert to SIMC units
    const double pe0 = pe[0]   * pscale;
    const double pe1 = pe[1]   * pscale;
    const double pe2 = pe[2]   * pscale;

    const double pr0 = prec[0] * pscale;
    const double pr1 = prec[1] * pscale;
    const double pr2 = prec[2] * pscale;

    const double zcm = zvtx * xscale;

    // write dat  (unchanged)
    out << i << " "
        << pe0 << " " << pe1 << " " << pe2 << " "
        << pr0 << " " << pr1 << " " << pr2 << " "
        << pidrec << " "
        << zcm << " "
        << weight << "\n";

    // ----------------------------
    // Recoil-particle kinematics
    // ----------------------------
    const double p_rec  = std::sqrt(pr0*pr0 + pr1*pr1 + pr2*pr2);
    const double pt_rec = std::sqrt(pr0*pr0 + pr1*pr1);

    double theta_rec_deg = 0.0;
    if (p_rec > 0.0) {
      double cth = pr2 / p_rec;
      if (cth >  1.0) cth =  1.0;
      if (cth < -1.0) cth = -1.0;
      theta_rec_deg = std::acos(cth) * 180.0 / TMath::Pi();
    }

    const double phi_rec_deg = std::atan2(pr1, pr0) * 180.0 / TMath::Pi();

    h_theta_phi_rec->Fill(phi_rec_deg, theta_rec_deg, weight);
    h_p_rec->Fill(p_rec, weight);
    h_pz_vs_p_rec->Fill(pt_rec, pr2, weight);
    h_p_vs_theta_rec->Fill(theta_rec_deg, p_rec, weight);

    // ----------------------------
    // Scattered-electron kinematics
    // ----------------------------
    const double p_e  = std::sqrt(pe0*pe0 + pe1*pe1 + pe2*pe2);
    const double pt_e = std::sqrt(pe0*pe0 + pe1*pe1);

    double theta_e_deg = 0.0;
    if (p_e > 0.0) {
      double cth = pe2 / p_e;
      if (cth >  1.0) cth =  1.0;
      if (cth < -1.0) cth = -1.0;
      theta_e_deg = std::acos(cth) * 180.0 / TMath::Pi();
    }

    const double phi_e_deg = std::atan2(pe1, pe0) * 180.0 / TMath::Pi();

    h_theta_phi_e->Fill(phi_e_deg, theta_e_deg, weight);
    h_p_e->Fill(p_e, weight);
    h_pz_vs_p_e->Fill(pt_e, pe2, weight);
    h_p_vs_theta_e->Fill(theta_e_deg, p_e, weight);
  }

  std::cout << "Wrote " << N << " events to " << outFile << "\n";
  std::cout << "Unit scales used: pscale=" << pscale
            << " (to MeV/c), xscale=" << xscale << " (to cm)\n";

  // ============================
  // Make pidrec plot
  // ============================
  if (pidCount.empty()) {
    std::cerr << "WARNING: pidrec distribution is empty, skip plotting.\n";
    return;
  }

  std::vector<int> pdgs;
  pdgs.reserve(pidCount.size());
  for (auto &kv : pidCount) pdgs.push_back(kv.first);
  std::sort(pdgs.begin(), pdgs.end());

  TCanvas* c_pid = new TCanvas("c_pidrec", "pidrec", 700, 500);
  c_pid->SetLeftMargin(0.12);
  c_pid->SetRightMargin(0.04);
  c_pid->SetBottomMargin(0.20);
  c_pid->SetTopMargin(0.08);

  TH1F* h_pid = new TH1F("h_pidrec", "pidrec distribution;particle species (PDG);counts",
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
    lat.DrawLatex(i, y + 0.02 * (double)ymax, Form("%lld", (Long64_t)y));
  }

  c_pid->SaveAs(figPng.c_str());

  // ============================
  // Multi-page PDF output
  // ============================
  c_pid->Print((figPdf + "[").c_str());
  c_pid->Print(figPdf.c_str());

  // Recoil particle: theta-phi
  TCanvas* c_theta_phi_rec = new TCanvas("c_theta_phi_rec", "theta_phi_rec", 800, 650);
  c_theta_phi_rec->SetLeftMargin(0.12);
  c_theta_phi_rec->SetRightMargin(0.14);
  c_theta_phi_rec->SetBottomMargin(0.12);
  c_theta_phi_rec->SetTopMargin(0.08);
  h_theta_phi_rec->Draw("COLZ");
  c_theta_phi_rec->Print(figPdf.c_str());

  // Recoil particle: p
  TCanvas* c_p_rec = new TCanvas("c_p_rec", "p_magnitude_rec", 800, 650);
  c_p_rec->SetLeftMargin(0.12);
  c_p_rec->SetRightMargin(0.05);
  c_p_rec->SetBottomMargin(0.12);
  c_p_rec->SetTopMargin(0.08);
  h_p_rec->SetLineWidth(2);
  h_p_rec->SetFillStyle(0);
  h_p_rec->Draw("HIST");
  c_p_rec->Print(figPdf.c_str());

  // Recoil particle: pz vs p_perp
  TCanvas* c_pz_vs_p_rec = new TCanvas("c_pz_vs_p_rec", "pz_vs_p_perp_rec", 800, 650);
  c_pz_vs_p_rec->SetLeftMargin(0.12);
  c_pz_vs_p_rec->SetRightMargin(0.14);
  c_pz_vs_p_rec->SetBottomMargin(0.12);
  c_pz_vs_p_rec->SetTopMargin(0.08);
  h_pz_vs_p_rec->Draw("COLZ");
  c_pz_vs_p_rec->Print(figPdf.c_str());

  // Recoil particle: p vs theta
  TCanvas* c_p_vs_theta_rec = new TCanvas("c_p_vs_theta_rec", "p_vs_theta_rec", 800, 650);
  c_p_vs_theta_rec->SetLeftMargin(0.12);
  c_p_vs_theta_rec->SetRightMargin(0.14);
  c_p_vs_theta_rec->SetBottomMargin(0.12);
  c_p_vs_theta_rec->SetTopMargin(0.08);
  h_p_vs_theta_rec->Draw("COLZ");
  c_p_vs_theta_rec->Print(figPdf.c_str());

  // Scattered electron: theta-phi
  TCanvas* c_theta_phi_e = new TCanvas("c_theta_phi_e", "theta_phi_e", 800, 650);
  c_theta_phi_e->SetLeftMargin(0.12);
  c_theta_phi_e->SetRightMargin(0.14);
  c_theta_phi_e->SetBottomMargin(0.12);
  c_theta_phi_e->SetTopMargin(0.08);
  h_theta_phi_e->Draw("COLZ");
  c_theta_phi_e->Print(figPdf.c_str());

  // Scattered electron: p
  TCanvas* c_p_e = new TCanvas("c_p_e", "p_magnitude_e", 800, 650);
  c_p_e->SetLeftMargin(0.12);
  c_p_e->SetRightMargin(0.05);
  c_p_e->SetBottomMargin(0.12);
  c_p_e->SetTopMargin(0.08);
  h_p_e->SetLineWidth(2);
  h_p_e->SetFillStyle(0);
  h_p_e->Draw("HIST");
  c_p_e->Print(figPdf.c_str());

  // Scattered electron: pz vs p_perp
  TCanvas* c_pz_vs_p_e = new TCanvas("c_pz_vs_p_e", "pz_vs_p_perp_e", 800, 650);
  c_pz_vs_p_e->SetLeftMargin(0.12);
  c_pz_vs_p_e->SetRightMargin(0.14);
  c_pz_vs_p_e->SetBottomMargin(0.12);
  c_pz_vs_p_e->SetTopMargin(0.08);
  h_pz_vs_p_e->Draw("COLZ");
  c_pz_vs_p_e->Print(figPdf.c_str());

  // Scattered electron: p vs theta
  TCanvas* c_p_vs_theta_e = new TCanvas("c_p_vs_theta_e", "p_vs_theta_e", 800, 650);
  c_p_vs_theta_e->SetLeftMargin(0.12);
  c_p_vs_theta_e->SetRightMargin(0.14);
  c_p_vs_theta_e->SetBottomMargin(0.12);
  c_p_vs_theta_e->SetTopMargin(0.08);
  h_p_vs_theta_e->Draw("COLZ");
  c_p_vs_theta_e->Print(figPdf.c_str());

  c_p_vs_theta_e->Print((figPdf + "]").c_str());

  std::cout << "Saved plots to " << figPdf << "\n";
  std::cout << "Saved pidrec png to " << figPng << "\n";

  std::cout << "\nPIDREC counts:\n";
  for (auto &kv : pidCount) {
    std::cout << "  " << kv.first << " : " << kv.second << "\n";
  }
}