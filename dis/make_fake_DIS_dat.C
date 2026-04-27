#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

std::string FormatForName(double x) {
  // Integer if very close to integer; otherwise compact decimal text.
  const double rx = std::round(x);
  std::ostringstream os;
  if (std::fabs(x - rx) < 1e-9) {
    os << static_cast<long long>(rx);
  } else {
    os << std::fixed << std::setprecision(2) << x;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    for (char &c : s) {
      if (c == '.') c = 'p';
      if (c == '-') c = 'm';
    }
    return s;
  }
  return os.str();
}

void SphericalToCartesian(double p, double theta_deg, double phi_deg,
                          double &px, double &py, double &pz) {
  const double theta = theta_deg * TMath::DegToRad();
  const double phi   = phi_deg   * TMath::DegToRad();

  px = p * std::sin(theta) * std::cos(phi);
  py = p * std::sin(theta) * std::sin(phi);
  pz = p * std::cos(theta);
}

} // namespace

void make_fake_DIS_dat(
    Long64_t NEvents        = 2000,
    double theta_center_deg = 10.0,
    double delta_deg        = 1.0,
    double Ee_MeV           = 4400.0,
    double phi_e_deg        = 0.0,
    unsigned long seed      = 12345UL) {

  // ------------------------------------------------------------
  // Fixed beam and fixed recoil-proton settings
  // ------------------------------------------------------------
  const double Ebeam_MeV = 10200.0;
  const double me_MeV    = 0.51099895;

  // User-requested fixed recoil proton
  const double pRec_MeV      = 100.0;
  const double theta_rec_deg = 140.0;
  const double phi_rec_deg   = 0.0;

  const int    pidrec  = 2212;
  const double zvtx_cm = 0.0;
  const double weight  = 1.0;

  if (NEvents <= 0) {
    std::cerr << "ERROR: NEvents must be > 0.\n";
    return;
  }
  if (delta_deg < 0.0) {
    std::cerr << "ERROR: delta_deg must be >= 0.\n";
    return;
  }
  if (Ee_MeV <= me_MeV) {
    std::cerr << "ERROR: Ee_MeV must be larger than electron mass.\n";
    return;
  }
  if (Ee_MeV > Ebeam_MeV) {
    std::cerr << "WARNING: Ee_MeV is larger than beam energy. "
              << "The file will still be generated as requested.\n";
  }

  // Scattered electron momentum magnitude from fixed energy
  const double pe_MeV = std::sqrt(Ee_MeV * Ee_MeV - me_MeV * me_MeV);

  // Fixed recoil proton momentum components
  double pr0 = 0.0, pr1 = 0.0, pr2 = 0.0;
  SphericalToCartesian(pRec_MeV, theta_rec_deg, phi_rec_deg, pr0, pr1, pr2);

  const std::string outFile =
      std::string("fakeDIS_") + FormatForName(theta_center_deg) +
      std::string("deg_") + FormatForName(Ee_MeV) +
      std::string("MeV.dat");

  std::ofstream out(outFile.c_str());
  if (!out) {
    std::cerr << "ERROR: cannot open output file: " << outFile << "\n";
    return;
  }

  // Keep the same output style as your converter code
  out.setf(std::ios::scientific);
  out.precision(10);

  out << "# EvtID "
      << "pe_px[MeV/c] pe_py[MeV/c] pe_pz[MeV/c] "
      << "prec_px[MeV/c] prec_py[MeV/c] prec_pz[MeV/c] "
      << "pidrec zvtx[cm] weight\n";

  TRandom3 rng(seed);

  for (Long64_t i = 0; i < NEvents; ++i) {
    const double theta_e_deg = rng.Uniform(theta_center_deg - delta_deg,
                                           theta_center_deg + delta_deg);

    double pe0 = 0.0, pe1 = 0.0, pe2 = 0.0;
    SphericalToCartesian(pe_MeV, theta_e_deg, phi_e_deg, pe0, pe1, pe2);

    // Write exactly the same column order as your converter
    out << i << " "
        << pe0 << " " << pe1 << " " << pe2 << " "
        << pr0 << " " << pr1 << " " << pr2 << " "
        << pidrec << " "
        << zvtx_cm << " "
        << weight << "\n";
  }

  out.close();

  std::cout << "Wrote " << NEvents << " events to " << outFile << "\n";
  std::cout << "Beam energy fixed at " << Ebeam_MeV << " MeV\n";
  std::cout << "Electron energy fixed at " << Ee_MeV << " MeV\n";
  std::cout << "Electron theta sampled uniformly in ["
            << (theta_center_deg - delta_deg) << ", "
            << (theta_center_deg + delta_deg) << "] deg\n";
  std::cout << "Electron phi fixed at " << phi_e_deg << " deg\n";
  std::cout << "Recoil proton fixed at p = " << pRec_MeV
            << " MeV/c, theta = " << theta_rec_deg
            << " deg, phi = " << phi_rec_deg << " deg\n";
}