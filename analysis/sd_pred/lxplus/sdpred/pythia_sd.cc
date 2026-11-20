// pythia_sd.cc — Pythia 8, poziom cząstek: mezony D0/D+/D*+ + proton wiodący (ξ, t) po obu stronach + ξ z cząstek centralnych.
// Zadanie #8 planu prac z 17.11: przewidywania dla definicji fiducjalnej pomiaru D + tag AFP (0,030 ≤ ξ ≤ 0,100).
// Na bazie Q_lxplus/truth/pythia_dd.cc (ten sam wybór D, fromB, mode) + łatka z DIFF_model/RAPORT.md §5.
//
// Użycie: pythia_sd <nEvents> <seed> <out.tsv> <plik.cmnd> [plik2.cmnd ...]     (ustawienia WYŁĄCZNIE z plików cmnd, czytane po kolei)
// Wyjście TSV — tylko zdarzenia z ≥ 1 D (|pdg| 421/411/413, ostatnia kopia) o pT > 2 GeV, |η| < 2,5:
//   E ev code hd nch lead xiA tA xiC tC xcalA xcalC xtrkA xtrkC xtk1A xtk1C
//   D ev pdg pt eta y phi fromB mode
//  code = info.code(); hd = isHardDiffractiveA + 2·isHardDiffractiveB (tylko Diffraction:doHard)
//  strona A = +z, C = −z. xi = 1 − |pz|/p_wiązki protonu końcowego o NAJMNIEJSZYM ξ po danej stronie (−1 = brak protonu), t = (p_wiązki − p')² [GeV²] (< 0).
//  xcal = Σ(E ± pz)/√s po widocznych cząstkach końcowych |η| < 4,9 (A: +, C: −) — odpowiednik ξ_cal; xtrk: naładowane |η| < 2,5, pT > 0,5; xtk1: pT > 0,1.
//  nch/lead: naładowane pT > 0,5, |η| < 2,5 (jak pythia_dd).
// JSON (<out.tsv>.json): σ_gen i liczniki per proces dla WSZYSTKICH zdarzeń: n, protony w oknie 0,03–0,10 (A, C), histogram ξ (0–0,2, 40 binów, A+C).
#include "Pythia8/Pythia.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>
using namespace Pythia8;

static const double PTSTORE = 2.0, ETASTORE = 2.5, XILO = 0.030, XIHI = 0.100;

static bool isB(int id) { int a = std::abs(id); return (a / 100) % 10 == 5 || (a / 1000) % 10 == 5; }
static void finals(const Event& ev, int i, std::vector<int>& out, int depth = 0) {
  if (depth > 20) return;
  if (ev[i].isFinal()) { out.push_back(ev[i].idAbs()); return; }
  for (int k : ev[i].daughterList()) finals(ev, k, out, depth + 1);
}
static int modeOf(const Event& ev, int i) {
  const Particle& p = ev[i]; int a = p.idAbs();
  std::vector<int> fs; finals(ev, i, fs); std::sort(fs.begin(), fs.end());
  if (a == 421) return (fs == std::vector<int>{211, 321}) ? 1 : 0;
  if (a == 411) return (fs == std::vector<int>{211, 211, 321}) ? 1 : 0;
  if (a == 413) {
    std::vector<int> d; for (int k : p.daughterList()) d.push_back(ev[k].idAbs()); std::sort(d.begin(), d.end());
    if (d != std::vector<int>{211, 421}) return 0;
    return (fs == std::vector<int>{211, 211, 321}) ? 1 : 0;
  }
  return 0;
}
struct Cnt { long n = 0, tagA = 0, tagC = 0, tagAC = 0; std::vector<long> hxi = std::vector<long>(40, 0); };

int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "użycie: %s nEvents seed out.tsv a.cmnd [b.cmnd ...]\n", argv[0]); return 1; }
  long nEv = std::atol(argv[1]); int seed = std::atoi(argv[2]); std::string out = argv[3];
  Pythia pythia("", false);
  std::string cm;
  for (int i = 4; i < argc; ++i) { if (!pythia.readFile(argv[i])) { fprintf(stderr, "BŁĄD cmnd: %s\n", argv[i]); return 3; } cm += std::string(i > 4 ? "," : "") + argv[i]; }
  for (std::string s : {std::string("Random:setSeed = on"), "Random:seed = " + std::to_string(seed % 900000000 + 1),
                        std::string("Next:numberCount = 0"), std::string("Next:numberShowEvent = 0"), std::string("Next:numberShowInfo = 0"),
                        std::string("Next:numberShowProcess = 0"), std::string("Init:showChangedSettings = on"), std::string("Init:showChangedParticleData = off")})
    if (!pythia.readString(s)) return 3;
  if (!pythia.init()) { fprintf(stderr, "BŁĄD init\n"); return 4; }
  const double eCM = pythia.info.eCM();
  const bool doHard = pythia.settings.flag("Diffraction:doHard");

  std::ofstream f(out);
  f << "# pythia_sd seed=" << seed << " cmnd=" << cm << " eCM=" << eCM << "\n";
  std::map<int, Cnt> cnt;
  long nOk = 0, nFail = 0, nWithD = 0;
  char buf[512];
  for (long iEv = 0; iEv < nEv; ++iEv) {
    if (!pythia.next()) { if (++nFail > 10000) break; continue; }
    ++nOk;
    const Event& ev = pythia.event;
    int code = pythia.info.code();
    int hd = doHard ? (int)pythia.info.isHardDiffractiveA() + 2 * (int)pythia.info.isHardDiffractiveB() : 0;
    const Vec4 pA = ev[1].p(), pC = ev[2].p();   // wiązki (A: +z)
    const double pb = pA.pAbs();
    int nch = 0; double lead = 0., xiA = -1, tA = 0, xiC = -1, tC = 0;
    double xc[2] = {0, 0}, xt[2] = {0, 0}, x1[2] = {0, 0};
    std::vector<int> ds;
    for (int i = 0; i < ev.size(); ++i) {
      const Particle& p = ev[i];
      int a = p.idAbs();
      if ((a == 421 || a == 411 || a == 413) && p.iBotCopyId() == i && p.pT() > PTSTORE && std::abs(p.eta()) < ETASTORE) ds.push_back(i);
      if (!p.isFinal()) continue;
      if (p.id() == 2212) {
        double xi = 1. - std::abs(p.pz()) / pb;
        if (p.pz() > 0 && (xiA < 0 || xi < xiA)) { xiA = xi; tA = (pA - p.p()).m2Calc(); }
        if (p.pz() < 0 && (xiC < 0 || xi < xiC)) { xiC = xi; tC = (pC - p.p()).m2Calc(); }
      }
      if (!p.isVisible()) continue;
      double ae = std::abs(p.eta());
      if (ae < 4.9) { xc[0] += p.e() + p.pz(); xc[1] += p.e() - p.pz(); }
      if (p.isCharged() && ae < 2.5) {
        if (p.pT() > 0.5) { ++nch; lead = std::max(lead, p.pT()); xt[0] += p.e() + p.pz(); xt[1] += p.e() - p.pz(); }
        if (p.pT() > 0.1) { x1[0] += p.e() + p.pz(); x1[1] += p.e() - p.pz(); }
      }
    }
    Cnt& c = cnt[code * 10 + hd];
    c.n++;
    bool wA = xiA >= XILO && xiA <= XIHI, wC = xiC >= XILO && xiC <= XIHI;
    c.tagA += wA; c.tagC += wC; c.tagAC += (wA && wC);
    if (xiA >= 0 && xiA < 0.2) c.hxi[(int)(xiA / 0.005)]++;
    if (xiC >= 0 && xiC < 0.2) c.hxi[(int)(xiC / 0.005)]++;
    if (ds.empty()) continue;
    ++nWithD;
    snprintf(buf, sizeof buf, "E\t%ld\t%d\t%d\t%d\t%.3f\t%.6g\t%.5g\t%.6g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n", iEv, code, hd, nch, lead,
             xiA, tA, xiC, tC, xc[0] / eCM, xc[1] / eCM, xt[0] / eCM, xt[1] / eCM, x1[0] / eCM, x1[1] / eCM);
    f << buf;
    for (int i : ds) {
      const Particle& p = ev[i];
      bool fromB = false;
      std::vector<int> stack = p.motherList(); std::vector<char> seen(ev.size(), 0); int guard = 0;
      while (!stack.empty() && guard++ < 2000) {
        int m = stack.back(); stack.pop_back();
        if (m <= 0 || seen[m]) continue;
        seen[m] = 1;
        if (isB(ev[m].id())) { fromB = true; break; }
        for (int mm : ev[m].motherList()) stack.push_back(mm);
      }
      snprintf(buf, sizeof buf, "D\t%ld\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n", iEv, p.id(), p.pT(), p.eta(), p.y(), p.phi(), (int)fromB, modeOf(ev, i));
      f << buf;
    }
  }
  f.close();
  pythia.stat();
  std::ofstream js(out + ".json");
  js << "{\n \"seed\": " << seed << ", \"cmnd\": \"" << cm << "\", \"eCM\": " << eCM << ", \"doHard\": " << (doHard ? 1 : 0)
     << ", \"sampleType\": " << pythia.settings.mode("Diffraction:sampleType") << ", \"PomFlux\": " << pythia.settings.mode("SigmaDiffractive:PomFlux") << ",\n";
  js << " \"n_requested\": " << nEv << ", \"n_ok\": " << nOk << ", \"n_fail\": " << nFail << ", \"n_with_D\": " << nWithD << ",\n";
  js << " \"n_tried\": " << pythia.info.nTried() << ", \"n_selected\": " << pythia.info.nSelected() << ", \"n_accepted\": " << pythia.info.nAccepted() << ",\n";
  js << " \"sigma_gen_mb_total\": " << pythia.info.sigmaGen() << ", \"sigma_err_mb_total\": " << pythia.info.sigmaErr() << ",\n";
  js << " \"xi_window\": [" << XILO << ", " << XIHI << "], \"proc\": {";
  bool first = true;
  for (auto& kv : cnt) {
    int code = kv.first / 10, hd = kv.first % 10; const Cnt& c = kv.second;
    js << (first ? "" : ", ") << "\n  \"" << code << "_" << hd << "\": {\"n\": " << c.n << ", \"sigma_mb\": " << pythia.info.sigmaGen(code)
       << ", \"tagA\": " << c.tagA << ", \"tagC\": " << c.tagC << ", \"tagAC\": " << c.tagAC << ", \"hxi\": [";
    for (size_t i = 0; i < c.hxi.size(); ++i) js << (i ? "," : "") << c.hxi[i];
    js << "]}";
    first = false;
  }
  js << "},\n \"pythia_version\": " << PYTHIA_VERSION << "\n}\n";
  js.close();
  return 0;
}
