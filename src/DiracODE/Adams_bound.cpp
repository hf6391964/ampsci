#include "DiracODE/Adams_bound.hpp"
#include "DiracODE/Adams_coefs.hpp"
#include "DiracODE/DiracODE.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

// Unfortunately, seems too messy to fix this
#pragma GCC diagnostic ignored "-Wsign-conversion"

/*
Program to solve single-electron bound-state Dirac problem for a (given)
local, central potential.
Based on method presented in book by W. Johnson.
Employs the Adams-Moulton method.
boundState is the main routine that is called from elsewhere.
All other functions called by boundState.
*/

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

namespace DiracODE {

using namespace Adams;
//******************************************************************************
void boundState(DiracSpinor &psi, const double en0,
                const std::vector<double> &v, const std::vector<double> &H_mag,
                const double alpha, int log_dele, const DiracSpinor *const VxFa,
                double zion)
// Solves local, spherical bound state dirac equation using Adams-Moulton
// method. Based on method presented in book by W. R. Johnson:
//   W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)
// I have added a few extensions to this method. In particular, I integrate past
// the classical turning point. (See below)
//
// See also:
//  * https://en.wikipedia.org/wiki/Linear_multistep_method
//  * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
//    Solving ordinary differential equations I: Nonstiff problems (2nd ed.),
//    Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
//  * Quarteroni, Alfio; Sacco, Riccardo; Saleri, Fausto (2000),
//    Matematica Numerica, Springer Verlag, ISBN 978-88-470-0077-3.
//  * http://mathworld.wolfram.com/AdamsMethod.html
//
// Rough description of method:
// 1. Start with initial 'guess' of energy
// 2. Find "practical infinity" (psi~0), and Classical turning point (e=v)
// 3. Performs 'inward' integration (Adams Moulton).
//    Integrates inwards from the practical infinity inwards to d_ctp points
//    past the classical turning point (ctp).
// 4. Performs 'outward' integration (Adams Moulton). Integrates from 0
//    outwards to d_ctp points past the ctp.
// 5. Matches the two functions around ctp, by re-scaling the 'inward' solution
//    (for f). Uses a weighted mean, with weights given by distance from ctp.
// 6. Checks the number of nodes the wf has. If too many or too few nodes, makes
//    a large change to the energy and tries again (from step 2).
//    If the correct number of nodes, uses perturbation theory to make minor
//    corrections to the energy to 'zoom in' (matching the in/out solution for
//    g), then re-starts from step 2.
// Continues until this energy adjustment falls below a prescribed threshold.
//
// Orbitals defined:
//   psi := (1/r) {f O_k, ig O_(-k)}
//
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Convergance goal. Default: 1e-14
  const double eps_goal = std::pow(10, -std::abs(log_dele));

  DEBUG(if (!(std::abs(psi.k) <= psi.n && psi.k != psi.n)) {
    std::cerr << "\nFail96 in Adams: bad state " << psi.symbol() << "\n";
    return;
  })
  DEBUG(std::cerr << "Start: " << psi.symbol() << ", en=" << psi.en << "\n";)

  const auto &rgrid = *psi.rgrid;

  // orbital should have (n-l-1) nodes:
  const int required_nodes = psi.n - psi.l() - 1;
  bool correct_nodes = false;
  TrackEnGuess sofar; // track higest/lowest energy guesses etc.

  // Start eigenvalue iterations:
  double t_en = en0; // set initial energy guss
  int t_pinf = 0;
  double t_eps = 1.0;
  double t_eps_prev = 1.0;
  double anorm = 1.0;
  int t_its = 1;
  for (; t_its < Param::max_its; ++t_its) {
    t_pinf = Adams::findPracticalInfinity(t_en, v, rgrid.r, Param::cALR);
    const int ctp =
        Adams::findClassicalTurningPoint(t_en, v, t_pinf, Param::d_ctp);

    // Find solution (f,g) to DE for given energy:
    // Also stores dg (gout-gin) for PT [used for PT to find better e]
    std::vector<double> dg(2 * Param::d_ctp + 1);
    Adams::trialDiracSolution(psi.f, psi.g, dg, t_en, psi.k, v, H_mag, rgrid,
                              ctp, Param::d_ctp, t_pinf, alpha, VxFa, zion);

    const int counted_nodes = Adams::countNodes(psi.f, t_pinf);

    // If correct number of nodes, use PT to make minor energy adjustment.
    // Otherwise, make large adjustmunt until correct # of nodes
    const double en_old = t_en;
    if (counted_nodes == required_nodes) {
      correct_nodes = true;
      anorm = Adams::calcNorm(psi.f, psi.g, rgrid.drdu, rgrid.du, t_pinf);
      t_en = Adams::smallEnergyChangePT(en_old, anorm, psi.f, dg, ctp,
                                        Param::d_ctp, alpha, sofar);
    } else {
      correct_nodes = false;
      const bool toomany_nodes =
          (counted_nodes > required_nodes) ? true : false;
      Adams::largeEnergyChange(&t_en, &sofar, Param::lfrac_de, toomany_nodes);
    }
    t_eps = fabs((t_en - en_old) / en_old);

    DEBUG(std::cerr << " :: it=" << t_its << " nodes:" << counted_nodes << "/"
                    << required_nodes << " new_en = " << t_en
                    << " delta=" << t_eps * t_en << " eps=" << t_eps << "\n";
          std::cin.get();)

    auto getting_worse = (t_its > 10 && t_eps >= 1.2 * t_eps_prev &&
                          correct_nodes && t_eps < 1.e-5);
    auto converged = (t_eps < eps_goal && correct_nodes);
    if (converged || getting_worse)
      break;
    t_eps_prev = t_eps;
  } // END itterations

  // If we never got correct # of nodes, never calc'd norm constant.
  // This is rare - means a failure. But hopefully, failure will go away after
  // a few more HF iterations..If we don't norm wf, HF will fail.
  if (!correct_nodes) {
    anorm = Adams::calcNorm(psi.f, psi.g, rgrid.drdu, rgrid.du, t_pinf);
    DEBUG(std::cerr << "\nFAIL-148: wrong nodes:"
                    << Adams::countNodes(psi.f, t_pinf) << "/" << required_nodes
                    << " for " << psi.symbol() << "\n";)
  }

  // store energy etc.
  psi.en = t_en;
  psi.eps = t_eps;
  psi.pinf = (std::size_t)t_pinf;
  psi.its = t_its;

  // Explicitely set 'tail' to zero (we may be re-using orbital)
  for (auto i = psi.pinf; i < rgrid.num_points; i++) {
    psi.f[i] = 0;
    psi.g[i] = 0;
  }
  // normalises the orbital (and zero's after pinf)
  const double an = 1. / std::sqrt(anorm);
  for (auto i = 0ul; i < psi.pinf; i++) {
    psi.f[i] = an * psi.f[i];
    psi.g[i] = an * psi.g[i];
  }

  return;
}

//******************************************************************************
void regularAtOrigin(DiracSpinor &phi, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto &gr = phi.rgrid;
  if (en != 0)
    phi.en = en;
  const auto pinf = Adams::findPracticalInfinity(phi.en, v, gr->r, Param::cALR);
  Adams::DiracMatrix Hd(*gr, v, phi.k, phi.en, alpha, H_mag);
  Adams::outwardAM(phi.f, phi.g, Hd, pinf - 1);
  // for safety: make sure zerod! (I may re-use existing orbitals!)
  for (std::size_t i = pinf; i < gr->num_points; i++) {
    phi.f[i] = 0;
    phi.g[i] = 0;
  }
  phi.pinf = pinf;
}

//******************************************************************************
void regularAtInfinity(DiracSpinor &phi, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_mag, const double alpha) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto &gr = phi.rgrid;
  if (en < 0)
    phi.en = en;
  const auto pinf = Adams::findPracticalInfinity(phi.en, v, gr->r, Param::cALR);
  Adams::DiracMatrix Hd(*gr, v, phi.k, phi.en, alpha, H_mag);
  Adams::inwardAM(phi.f, phi.g, Hd, 0, pinf - 1);
  for (std::size_t i = pinf; i < gr->num_points; i++) {
    phi.f[i] = 0;
    phi.g[i] = 0;
  }
  phi.pinf = pinf;
}

//******************************************************************************
//******************************************************************************
namespace Adams {

//******************************************************************************
void largeEnergyChange(double *en, TrackEnGuess *sofar_ptr, double frac_de,
                       bool toomany_nodes)
// wf did not have correct number of nodes. Make a large energy adjustment
// toomany_nodes=true means there were too many nodes
{
  auto &sofar = *sofar_ptr; // for ease of typing only
  auto etemp = *en;
  if (toomany_nodes) {
    ++sofar.count_toomany;
    if ((sofar.count_toomany == 1) || (*en < sofar.high_en))
      sofar.high_en = *en;
    etemp *= (1.0 + frac_de);
    if ((sofar.count_toofew != 0) && (etemp < sofar.low_en))
      etemp = 0.5 * (sofar.high_en + sofar.low_en);
  } else {
    ++sofar.count_toofew;
    if ((sofar.count_toofew == 1) || (*en > sofar.low_en))
      sofar.low_en = *en;
    etemp *= (1.0 - frac_de);
    if ((sofar.count_toomany != 0) && (etemp > sofar.high_en))
      etemp = 0.5 * (sofar.high_en + sofar.low_en);
  }
  *en = etemp;
}

//******************************************************************************
double calcNorm(const std::vector<double> &f, const std::vector<double> &g,
                const std::vector<double> &drdu, const double du,
                const int pinf) {
  const auto anormF = NumCalc::integrate(1.0, 0, pinf, f, f, drdu);
  const auto anormG = NumCalc::integrate(1.0, 0, pinf, g, g, drdu);
  return (anormF + anormG) * du;
}

//******************************************************************************
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, const int ctp,
                           const int d_ctp, const double alpha,
                           const TrackEnGuess &sofar)
// delta E = c*f(r)*[g_out(r)-g_in(r)] - evaluate at ctp
// nb: wf not yet normalised (anorm is input param)!
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  double p_del_q = f[ctp] * dg[d_ctp];
  double denom = 1.0;
  // weighted average around ctp:
  for (int i = 1; i <= d_ctp; i++) {
    const auto w = Param::weight(i);
    p_del_q +=
        0.5 * (f[ctp + i] * dg[d_ctp + i] + f[ctp - i] * dg[d_ctp - i]) * w;
    denom += w;
  }

  const double de = p_del_q / (alpha * anorm * denom);
  double new_en = en + de;

  if ((sofar.count_toofew != 0) && (new_en < sofar.low_en)) {
    new_en = 0.5 * (en + sofar.low_en);
  } else if ((sofar.count_toomany != 0) && (new_en > sofar.high_en)) {
    new_en = 0.5 * (en + sofar.high_en);
  } else if (new_en > 0) {
    // This only happens v. rarely. nodes correct, but P.T. gives silly result!
    new_en = (de > 0) ? 0.9 * en : 1.1 * en;
  }

  return new_en;
}

//******************************************************************************
int findPracticalInfinity(const double en, const std::vector<double> &v,
                          const std::vector<double> &r, const double alr)
// Find the practical infinity 'pinf'
// Step backwards from the last point (num_points-1) until
// (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
// XXX Note: unsafe, and a little slow. Would be better to use
// std::lower_bound.. but would need lambda or something for r^2 part?
{
  auto pinf = r.size() - 5;
  while ((en - v[pinf - 1]) * r[pinf - 1] * r[pinf - 1] + alr < 0) {
    --pinf;
    DEBUG(if (pinf <= 1) {
      std::cerr << "Fail290 in WF: pinf underflowed?\n";
      std::cin.get();
    })
  }
  return (int)pinf + 5;
}

//******************************************************************************
int findClassicalTurningPoint(const double en, const std::vector<double> &v,
                              const int pinf, const int d_ctp)
// Finds classical turning point 'ctp'
// Enforced to be between (0+ctp) and (pinf-ctp)
//  V(r) > E        [nb: both V and E are <0]
{
  auto low =
      std::lower_bound(v.begin() + d_ctp + 1, v.begin() + pinf - d_ctp, en);
  return (int)(low - v.begin()) - 1;
}

//******************************************************************************
int countNodes(const std::vector<double> &f, const int pinf)
// Just counts the number of times orbital (f) changes sign
{
  double sp = f[1];
  int counted_nodes = 0;
  for (int i = 2; i < pinf; ++i) {
    if (sp * f[i] < 0) { // ok if f[i]==0 ?? (-0<0?) XXX
      ++counted_nodes;
      sp = f[i];
    }
  }
  return counted_nodes;
}

//******************************************************************************
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v,
                        const std::vector<double> &H_mag, const Grid &gr,
                        const int ctp, const int d_ctp, const int pinf,
                        const double alpha, const DiracSpinor *const VxFa,
                        double zion)
// Performs inward (from pinf) and outward (from r0) integrations for given
// energy. Intergated in/out towards ctp +/- d_ctp [class. turn. point]
// Then, joins solutions, including weighted meshing around ctp +/ d_ctp
// Also: stores dg [the difference: (gout-gin)], which is used for PT
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  DiracMatrix Hd(gr, v, ka, en, alpha, H_mag, VxFa, zion);
  outwardAM(f, g, Hd, ctp + d_ctp);
  std::vector<double> f_in(gr.num_points), g_in(gr.num_points);
  inwardAM(f_in, g_in, Hd, ctp - d_ctp, pinf - 1);
  joinInOutSolutions(f, g, dg, f_in, g_in, ctp, d_ctp, pinf);
}

//******************************************************************************
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, const int ctp,
                        const int d_ctp, const int pinf) {

  // Find the re-scaling factor (for inward soln)
  // Average of points ctp+/-d_ctp [using a weighted average]
  double rescale = f[ctp] / f_in[ctp];
  double denom = 1;
  for (int i = 1; i <= d_ctp; i++) {
    auto w = Param::weight(i);
    rescale +=
        0.5 * (f[ctp + i] / f_in[ctp + i] + f[ctp - i] / f_in[ctp - i]) * w;
    denom += w;
  }
  rescale /= denom;

  // store difference between in/out solutions (for g) - after re-scaling
  // Used later for P.T.
  for (auto i = 0; i < (int)dg.size(); i++) {
    dg[i] = g[ctp - d_ctp + i] - g_in[ctp - d_ctp + i] * rescale;
  }

  // Join the in and outward solutions. "Meshed" around ctp +/- d_ctp
  for (int i = ctp - d_ctp; i <= ctp + d_ctp; i++) {
    //"Mesh" in the intermediate region, using weighted av.
    auto b =
        (i - ctp > 0) ? 1. - Param::weight(i - ctp) : Param::weight(i - ctp);
    auto a = 1. - b;
    f[i] = a * f[i] + b * f_in[i] * rescale;
    g[i] = a * g[i] + b * g_in[i] * rescale;
  }
  for (int i = ctp + d_ctp + 1; i < pinf; i++) {
    f[i] = f_in[i] * rescale;
    g[i] = g_in[i] * rescale;
  }
}

//******************************************************************************
void outwardAM(std::vector<double> &f, std::vector<double> &g,
               const DiracMatrix &Hd, const int nf)
// Program to start the OUTWARD integration.
// Starts from 0, and uses an expansion(?) to go to (num_loops*AMO).
// Then, it then call ADAMS-MOULTON, to finish
// (from num_loops*AMO+1 to nf = ctp+d_ctp)
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  const auto &r = Hd.pgr->r;
  const auto &drduor = Hd.pgr->drduor;
  const auto du = Hd.pgr->du;
  const auto &v = *(Hd.v);
  const auto alpha = Hd.alpha;
  const auto ka = Hd.k;
  const double az0 = -1 * v[0] * r[0] * alpha;
  const auto ka2 = (double)(ka * ka);
  const double ga0 = std::sqrt(ka2 - az0 * az0);

  // initial wf values
  // P(r) = r^gamma u(r) // f(r) = P(r)
  // Q(r) = r^gamma v(r) // g(r) = -Q(r)
  auto u0 = 1.0;
  auto v0 = (ka > 0) ? -(ga0 + ka) / az0 : az0 / (ga0 - ka);
  f[0] = std::pow(r[0], ga0) * u0;
  g[0] = -std::pow(r[0], ga0) * v0;

  // loop through and find first Param::num_loops*AMO points of wf
  for (int ln = 0; ln < Param::num_loops; ln++) {
    const int i0 = ln * Param::AMO + 1;

    // defines/populates em expansion coeficients (then inverts)
    std::array<double, Param::AMO> coefa, coefb, coefc, coefd;
    std::array<double, Param::AMO> ga;
    LinAlg::SqMatrix em(Param::AMO);
    const auto oid_du = Param::AMcoef.OId * du;
    for (int i = 0; i < Param::AMO; i++) {
      const auto az = -v[i + i0] * r[i + i0] * alpha;
      ga[i] = std::sqrt(ka * ka - az * az);
      const auto dror = drduor[i + i0];
      coefa[i] = oid_du * (Hd.a(i + i0) - ga[i] * dror);
      coefb[i] = oid_du * Hd.b(i + i0);
      coefc[i] = oid_du * Hd.c(i + i0);
      coefd[i] = oid_du * (Hd.d(i + i0) - ga[i] * dror);
      for (int j = 0; j < Param::AMO; j++) {
        em[i][j] = Param::AMcoef.OIe[i][j];
      }
      em[i][i] -= coefd[i];
    }
    // from here on, em is the inverted matrix
    em.invert();

    // defines/populates fm, s coefs
    std::array<double, Param::AMO> s;
    LinAlg::SqMatrix fm(Param::AMO);
    for (int i = 0; i < Param::AMO; i++) {
      s[i] = -Param::AMcoef.OIa[i] * u0;
      for (int j = 0; j < Param::AMO; j++) {
        fm[i][j] = Param::AMcoef.OIe[i][j] - coefb[i] * em[i][j] * coefc[j];
        s[i] -= coefb[i] * em[i][j] * Param::AMcoef.OIa[j] * v0;
      }
      fm[i][i] -= coefa[i];
    }
    // from here on, fm is the inverted matrix
    fm.invert();

    // writes u(r) in terms of coefs and the inverse of fm
    // P(r) = r^gamma u(r)
    std::array<double, Param::AMO> us;
    for (int i = 0; i < Param::AMO; i++) {
      us[i] = 0;
      for (int j = 0; j < Param::AMO; j++) {
        us[i] += fm[i][j] * s[j];
      }
    }

    // writes v(r) in terms of coefs + u(r)
    // Q(r) = r^gamma v(r)
    std::array<double, Param::AMO> vs;
    for (int i = 0; i < Param::AMO; i++) {
      vs[i] = 0;
      for (int j = 0; j < Param::AMO; j++) {
        vs[i] += em[i][j] * (coefc[j] * us[j] - Param::AMcoef.OIa[j] * v0);
      }
    }

    // writes wavefunction: P= r^gamma u(r) etc..
    for (int i = 0; i < Param::AMO; i++) {
      const auto r_ga = std::pow(r[i + i0], ga[i]);
      f[i + i0] = r_ga * us[i];
      g[i + i0] = -r_ga * vs[i];
    }

    // re-sets 'starting point' for next ln
    u0 = us.back();
    v0 = vs.back();

  } // END loop through outint [num_loops]

  // Call adamsmoulton to finish integration from (num_loops*AMO) to ctp+d_ctp
  const auto na = Param::num_loops * Param::AMO + 1;
  if (nf > na)
    adamsMoulton(f, g, Hd, na, nf);

  return;
}

//******************************************************************
void inwardAM(std::vector<double> &f, std::vector<double> &g,
              const DiracMatrix &Hd, const int nf, const int pinf)
// Program to start the INWARD integration.
// Starts from Pinf, and uses an expansion [WKB approx] to go to (pinf-AMO)
// i.e., gets last AMO points
// Then, it then call ADAMS-MOULTON, to finish (from num_loops*AMO+1
//   to nf = ctp-d_ctp)
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // short-cuts
  const auto alpha = Hd.alpha;
  const auto ka = Hd.k;
  const auto en = Hd.en;
  const auto &r = Hd.pgr->r;
  const auto &v = *(Hd.v);

  const auto alpha2 = alpha * alpha;
  const auto cc = 1. / alpha;
  const auto c2 = 1. / alpha2;
  const auto ka2 = (double)(ka * ka);

  const auto lambda = std::sqrt(-en * (2. + en * alpha2));
  const auto zeta = Hd.VxFa ? Hd.zion : -v[pinf] * r[pinf]; // XXX ?
  const auto zeta2 = zeta * zeta;
  const auto sigma = (1.0 + en * alpha2) * (zeta / lambda);
  const auto Ren = en + c2; // total relativistic energy

  // Generates the expansion coeficients for asymptotic wf up to order nx
  std::array<double, Param::nx> bx;
  std::array<double, Param::nx> ax;
  bx[0] = (ka + (zeta / lambda)) * (0.5 * alpha);
  for (int i = 0; i < Param::nx; i++) {
    ax[i] = (ka + (i + 1 - sigma) * Ren * alpha2 - zeta * lambda * alpha2) *
            bx[i] * cc / ((i + 1) * lambda);
    if (i < (Param::nx - 1))
      bx[i + 1] =
          (ka2 - std::pow((double(i + 1) - sigma), 2) - zeta2 * alpha2) *
          bx[i] / (2 * (i + 1) * lambda);
  }

  // Generates last `AMO' points for P and Q [actually AMO+1?]
  const double f1 = std::sqrt(1. + en * alpha2 * 0.5);
  const double f2 = std::sqrt(-en * 0.5) * alpha;
  for (int i = pinf; i >= (pinf - Param::AMO); i--) {
    const double rfac = std::pow(r[i], sigma) * exp(-lambda * r[i]);
    double ps = 1.;
    double qs = 0.;
    double rk = 1.;
    double xe = 1.;
    for (int k = 0; k < Param::nx; k++) {
      // this will loop until a) converge, or b) k=Param::nx
      rk *= r[i];
      ps += (ax[k] / rk);
      qs += (bx[k] / rk);
      xe = std::fmax(std::fabs(ax[k] / ps), std::fabs(bx[k] / qs)) / rk;
      if (xe < Param::nx_eps) {
        break;
      }
    }
    DEBUG(if (xe > 1.e-3) std::cerr
              << "WARNING: Asymp. expansion in inwardAM didn't converge: " << i
              << " " << xe << "\n";)
    f[i] = rfac * (f1 * ps + f2 * qs);
    g[i] = rfac * (f1 * qs - f2 * ps);
  }

  if ((pinf - Param::AMO - 1) >= nf)
    adamsMoulton(f, g, Hd, pinf - Param::AMO - 1, nf);
}

//******************************************************************************
void adamsMoulton(std::vector<double> &f, std::vector<double> &g,
                  const DiracMatrix &Hd, const int ni, const int nf)
// program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
//   * ni is starting (initial) point for integration
//   * nf is end (final) point for integration (nf=ctp+/-d_ctp)
{
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto nosteps = std::abs(nf - ni) + 1; // number of integration steps
  const auto inc = (nf > ni) ? 1 : -1;        //'increment' for integration
  if (nf == ni) {
    // This ever happen?
    return;
  }

  const auto du = Hd.pgr->du;
  // create arrays for wf derivatives
  // const auto num_points = Hd.pgr->num_points;
  // std::vector<double> df(num_points), dg(num_points);
  std::array<double, Param::AMO> df, dg;
  std::array<double, Param::AMO> am;
  // int k1 = ni - inc * Param::AMO; // nb: k1 is iterated
  for (auto i = 0, k1 = ni - inc * Param::AMO; i < Param::AMO; i++, k1 += inc) {
    df[i] = (Hd.a(k1) * f[k1] - Hd.b(k1) * g[k1]);
    dg[i] = (Hd.d(k1) * g[k1] - Hd.c(k1) * f[k1]);
    if (Hd.VxFa) {
      df[i] += -Hd.VxFa->g[k1] * Hd.pgr->drdu[k1]; // * du;
      dg[i] += Hd.VxFa->f[k1] * Hd.pgr->drdu[k1];  //* du;
    }
    am[i] = du * inc * Param::AMcoef.AMd * Param::AMcoef.AMa[i];
    // k1 += inc;
  }

  // integrates the function from ni to the c.t.p
  const double a0 = inc * du * Param::AMcoef.AMd * Param::AMcoef.AMaa;
  // int k2 = ni;
  for (int i = 0, k2 = ni; i < nosteps; i++, k2 += inc) {
    const double a = Hd.a(k2);
    const double b = Hd.b(k2);
    const double c = Hd.c(k2);
    const double d = Hd.d(k2);
    const double det_inv = 1.0 / (1.0 - a0 * a0 * (b * c - a * d));
    double sf = f[k2 - inc];
    double sg = g[k2 - inc];
    for (int l = 0; l < Param::AMO; l++) {
      sf += am[l] * df[l];
      sg += am[l] * dg[l];
    }
    f[k2] = (sf - a0 * (b * sg + d * sf)) * det_inv;
    g[k2] = (sg - a0 * (c * sf + a * sg)) * det_inv;

    if (k2 == nf)
      break;

    // loads next 'first' k values:
    for (int l = 0; l < (Param::AMO - 1); l++) {
      df[l] = df[l + 1];
      dg[l] = dg[l + 1];
    }
    // loads next 'first' deriv's:
    df[Param::AMO - 1] = a * f[k2] - b * g[k2];
    dg[Param::AMO - 1] = d * g[k2] - c * f[k2];
    if (Hd.VxFa) {
      df[Param::AMO - 1] += -Hd.VxFa->g[k2] * Hd.pgr->drdu[k2]; // * du;
      dg[Param::AMO - 1] += Hd.VxFa->f[k2] * Hd.pgr->drdu[k2];  // * du;
    }
    // k2 += inc;
  }

} // END adamsmoulton

//******************************************************************************
//******************************************************************************
DiracMatrix::DiracMatrix(const Grid &in_grid, const std::vector<double> &in_v,
                         const int in_k, const double in_en,
                         const double in_alpha,
                         const std::vector<double> &in_Hmag,
                         const DiracSpinor *const in_VxFa, double in_zion)
    : pgr(&in_grid),
      v(&in_v),
      Hmag(in_Hmag.empty() ? nullptr : &in_Hmag),
      VxFa(in_VxFa),
      zion(in_zion),
      k(in_k),
      en(in_en),
      alpha(in_alpha),
      c2(1.0 / in_alpha / in_alpha) {}

// update a and d for off-diag additional potential (magnetic form-fac, QED)
double DiracMatrix::a(std::size_t i) const {
  const auto h_mag = (Hmag == nullptr) ? 0.0 : (*Hmag)[i];
  return (double(-k)) * pgr->drduor[i] + h_mag * pgr->drdu[i];
  // return (double(-k)) * pgr->drduor[i];
}
double DiracMatrix::b(std::size_t i) const {
  return -alpha * (en + 2.0 * c2 - (*v)[i]) * pgr->drdu[i];
}
double DiracMatrix::c(std::size_t i) const {
  return alpha * (en - (*v)[i]) * pgr->drdu[i];
}
double DiracMatrix::d(std::size_t i) const {
  const auto h_mag = (Hmag == nullptr) ? 0.0 : (*Hmag)[i];
  // return double(k) * pgr->drduor[i];
  return double(k) * pgr->drduor[i] - h_mag * pgr->drdu[i];
}

} // namespace Adams
} // namespace DiracODE
