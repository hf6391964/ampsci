#include "DoubleCorePol.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "HF/MixedStates.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>

namespace HF {

// Inherit from ExternalField??
//******************************************************************************
DoubleCorePol::DoubleCorePol(const ExternalField *const dV1,
                             const ExternalField *const dV2)
    : m_dV1(dV1),
      m_dV2(dV2),
      m_rank(dV1->m_rank + dV2->m_rank),
      m_pi(dV1->m_pi * dV2->m_pi),
      p_core(dV2->p_core),
      m_vl(dV1->m_vl),     // for now, const Vl (same each l)
      m_Hmag(dV1->m_Hmag), //(same each l)
      m_alpha(dV1->m_alpha),
      m_imag(dV1->m_imag != dV2->m_imag),
      p_VBr(dV1->p_VBr)
//
{
  initialise_dPsi();
}

//******************************************************************************
void DoubleCorePol::initialise_dPsi() {
  // Initialise dPsi vectors, accounting for selection rules
  const bool print = false;
  m_X.resize(p_core->size());
  for (auto ic = 0u; ic < p_core->size(); ic++) {
    const auto &Fc = (*p_core)[ic];
    const auto pi_ch = Fc.parity() * m_pi;
    const auto tj_c = Fc.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    if (print)
      std::cout << "|" << Fc.symbol() << ">  -->  ";
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Angular::parity_l(l_minus) * pi_ch;
      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      const auto kappa = Angular::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, Fc.rgrid);
      m_X[ic].back().pinf = Fc.pinf;
      if (print)
        std::cout << "|" << m_X[ic].back().symbol() << "> + ";
    }
    if (print)
      std::cout << "\n";
  }
  m_Y = m_X;
}

//******************************************************************************
void DoubleCorePol::clear_dPsi() {
  // re-set p0/pinf? no need.
  for (auto &mx : m_X) {
    for (auto &m : mx) {
      m *= 0.0;
    }
  }
  m_Y = m_X;
}

//******************************************************************************
const std::vector<DiracSpinor> &DoubleCorePol::get_dPsis(const DiracSpinor &Fc,
                                                         dPsiType XorY) const {
  const auto index = static_cast<std::size_t>(
      std::find(p_core->cbegin(), p_core->cend(), Fc) - p_core->cbegin());
  //  Better way?
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//******************************************************************************
const DiracSpinor &DoubleCorePol::get_dPsi_x(const DiracSpinor &Fc,
                                             dPsiType XorY,
                                             const int kappa_x) const {
  const auto &dPsis = get_dPsis(Fc, XorY);
  auto match_kappa_x = [=](const auto &Fa) { return Fa.k == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//******************************************************************************
double DoubleCorePol::rme3js(const int twoja, const int twojb, int two_mb,
                             int two_q) const
// XXX Correct??
// rme3js = (-1)^{ja-ma} (ja, k, jb,\ -ma, q, mb)
{
  auto two_ma = two_mb - two_q; // -ma + mb + q = 0;
  // sig = (-1)^(ja - ma)
  auto sig = ((twoja - two_ma) / 2) % 2 == 0 ? 1 : -1;
  return sig *
         Angular::threej_2(twoja, 2 * m_rank, twojb, -two_ma, two_q, two_mb);
}

//******************************************************************************
void DoubleCorePol::solve_TDHFcore(const double omega, const int max_its,
                                   const bool print) {

  // XXX Note: ONLY this routine (and constructor) is different from regular
  // TDHF method. Spo, inherit from that, and over-write?
  // ALSO: Most of this routine is exactly the same. Break the 'rhs' into a
  // sepperate function, only overwrite that one
  // Maybe: Also put into ExternalField dir, call TDHF, RPAD, DCP???
  // --> All overload! CorePolarisation (or RPA) TDHartreeFock,

  const double converge_targ = 1.0e-9;
  const auto damper = rampedDamp(0.95, 0.5, 1, 20);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  const auto imag = m_imag;

  auto eps = 0.0;
  double ceiling_eps = 1.0;
  int worse_count = 0;
  double extra_damp = 0.0;
  int it = 1;
  if (print) {
    printf("TDHF (w=%.3f): .. \r", omega);
    std::cout << std::flush;
  }
  for (; it <= max_its; it++) {
    eps = 0.0;
    const auto a_damp = (it == 1) ? 0.0 : damper(it) + extra_damp;

    // eps for solveMixedState - doesn't need to be small!
    const auto eps_ms = (it == 1) ? 1.0e-8 : 1.0e-3;

    auto tmp_X = m_X;
    auto tmp_Y = m_Y;
#pragma omp parallel for
    for (auto ic = 0ul; ic < p_core->size(); ic++) {
      const auto &Fc = (*p_core)[ic];
      auto eps_c = 0.0;

      for (auto j = 0ul; j < tmp_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        const auto &X1s = m_dV1->get_dPsis(Fc, dPsiType::X);
        const auto &X2s = m_dV2->get_dPsis(Fc, dPsiType::X);
        // XXX Check selection rules!
        // XXX Three-j symbol required?? X2 and Fc different?
        auto rhs = dV_rhs(Xx.k, Fc, false);
        for (const auto &X1 : X1s)
          // selection rule!
          rhs += m_dV2->dV_rhs(Xx.k, X1, false);
        for (const auto &X2 : X2s)
          rhs += m_dV1->dV_rhs(Xx.k, X2, false);

        if (Xx.k == Fc.k && !imag) {
          const auto de = Fc * rhs;
          rhs -= de * Fc;
        }
        HF::solveMixedState(Xx, Fc, omega, m_vl, m_alpha, *p_core, rhs, eps_ms,
                            nullptr, p_VBr, m_Hmag);
        Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        if (delta > eps_c)
          eps_c = delta;
      }
      if (!staticQ) {
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          auto &Yx = tmp_Y[ic][j];
          const auto &oldY = m_Y[ic][j];
          const auto &Y1s = m_dV1->get_dPsis(Fc, dPsiType::Y);
          const auto &Y2s = m_dV2->get_dPsis(Fc, dPsiType::Y);
          // XXX Check selection rules!
          // XXX Three-j symbol required?? X2 and Fc different?
          auto rhs = dV_rhs(Yx.k, Fc, true);
          for (const auto &Y1 : Y1s)
            // selection rule!
            rhs += m_dV2->dV_rhs(Yx.k, Y1, true);
          for (const auto &Y2 : Y2s)
            rhs += m_dV1->dV_rhs(Yx.k, Y2, true);
          if (Yx.k == Fc.k && !imag) {
            const auto de = Fc * rhs;
            rhs -= de * Fc;
          }
          HF::solveMixedState(Yx, Fc, omega, m_vl, m_alpha, *p_core, rhs,
                              eps_ms, nullptr, p_VBr, m_Hmag);
          Yx = a_damp * oldY + (1.0 - a_damp) * Yx;
        }
      } else {
        const auto s = imag ? -1 : 1;
        for (auto j = 0ul; j < tmp_Y[ic].size(); j++) {
          tmp_Y[ic][j] = s * tmp_X[ic][j];
        }
      }
#pragma omp critical(compare)
      if (eps_c > eps)
        eps = eps_c;
    }
    m_X = tmp_X;
    m_Y = tmp_Y;

    if (print && it > 0 && it % 15 == 0) {
      printf("TDHF (w=%.3f): %2i %.1e \r", omega, it, eps);
      std::cout << std::flush;
    }

    if (it > 15 && eps > 1.1 * ceiling_eps) {
      ++worse_count;
      extra_damp = (it % 2) ? 0.5 : 0.25;
    } else if (eps < ceiling_eps) {
      worse_count = 0;
    }
    ceiling_eps = std::min(eps, ceiling_eps);

    if ((it > 1 && eps < converge_targ) || worse_count > 3)
      break;
  }
  if (print) {
    printf("TDHF (w=%.3f): %2i %.1e\n", omega, it, eps);
  }
  m_core_eps = eps;
}

//******************************************************************************
// does it matter if a or b is in the core?
double DoubleCorePol::dV(const DiracSpinor &Fn, const DiracSpinor &Fm,
                         bool conj, const DiracSpinor *const Fexcl) const {
  auto s = conj && m_imag ? -1 : 1; // careful, not always needed
  return s * Fn * dV_rhs(Fn.k, Fm, conj, Fexcl);
}

double DoubleCorePol::dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  auto conj = Fm.en > Fn.en;
  return dV(Fn, Fm, conj);
}

//******************************************************************************
DiracSpinor DoubleCorePol::dV_rhs(const int kappa_n, const DiracSpinor &Fm,
                                  bool conj,
                                  const DiracSpinor *const Fexcl) const {

  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = !conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_rank;
  const auto tkp1 = double(2 * k + 1);

  const auto tjn = Angular::twoj_k(kappa_n);
  const auto tjm = Fm.twoj();
  const auto Ckala = Angular::Ck_kk(k, kappa_n, Fm.k);

  auto dVFm = DiracSpinor(0, kappa_n, Fm.rgrid);
  dVFm.pinf = Fm.pinf;

#pragma omp parallel for
  for (auto ib = 0ul; ib < p_core->size(); ib++) {
    const auto &Fb = (*p_core)[ib];
    const auto tjb = Fb.twoj();
    const auto &X_betas = get_dPsis(Fb, ChiType);
    const auto &Y_betas = get_dPsis(Fb, EtaType);
    auto dVFm_c = DiracSpinor(0, kappa_n, Fm.rgrid);
    dVFm_c.pinf = Fm.pinf;

    // only for testing: exclude certain (core) states from dV sum
    if (Fexcl && (*Fexcl) == Fb)
      continue;

    for (auto ibeta = 0ul; ibeta < X_betas.size(); ++ibeta) {
      const auto &X_beta = X_betas[ibeta];
      const auto &Y_beta = Y_betas[ibeta];
      const auto tjbeta = X_beta.twoj();

      // Direct part:
      const auto Ckbeb = Angular::Ck_kk(k, X_beta.k, Fb.k);
      if (Ckala != 0 && Ckbeb != 0) {
        dVFm_c += (Ckala * Ckbeb / tkp1) *
                  Coulomb::Rkv_bcd(kappa_n, Fb, Fm, X_beta + Y_beta, k);
      }

      const auto s = Angular::evenQ_2(tjbeta - tjm) ? 1 : -1;

      // exchange part (X):
      const auto l_min_X =
          std::max(std::abs(tjn - tjbeta), std::abs(tjm - tjb)) / 2;
      const auto l_max_X = std::min((tjn + tjbeta), (tjm + tjb)) / 2;
      for (int l = l_min_X; l <= l_max_X; ++l) {
        const auto sixj = Angular::sixj_2(tjm, tjn, 2 * k, tjbeta, tjb, 2 * l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Angular::evenQ(k + l) ? 1 : -1;
        const auto Ckba = Angular::Ck_kk(l, Fm.k, Fb.k);
        const auto Ckalbe = Angular::Ck_kk(l, kappa_n, X_beta.k);
        if (Ckba == 0 || Ckalbe == 0)
          continue;
        dVFm_c +=
            (s * m1kpl * Ckba * Ckalbe * sixj) *
            Coulomb::Rkv_bcd(kappa_n, Fm, X_beta, Fb, l); // xxx use  Ymb ?
      }

      // exchange part (Y):
      const auto l_min_Y =
          std::max(std::abs(tjn - tjb), std::abs(tjm - tjbeta)) / 2;
      const auto l_max_Y = std::min((tjn + tjb), (tjm + tjbeta)) / 2;
      for (int l = l_min_Y; l <= l_max_Y; ++l) {
        const auto sixj = Angular::sixj_2(tjm, tjn, 2 * k, tjb, tjbeta, 2 * l);
        if (sixj == 0)
          continue;
        const auto m1kpl = Angular::evenQ(k + l) ? 1 : -1;
        const auto Ckbea = Angular::Ck_kk(l, Fm.k, Y_beta.k);
        const auto Ckbal = Angular::Ck_kk(l, kappa_n, Fb.k);
        if (Ckbea == 0 || Ckbal == 0)
          continue;
        dVFm_c += (s * m1kpl * Ckbea * Ckbal * sixj) *
                  Coulomb::Rkv_bcd(kappa_n, Fm, Fb, Y_beta, l);
      }

      // Breit part:
      if (p_VBr) {
        // Note: Not perfectly symmetric for E1 - some issue??
        // But, differences is extremely small, so maybe just numerics?*
        // Assymetry enters below number of digits VD presents..
        // nb: Agrees perfectly w/ Vladimir for E1, E2, and PNC
        dVFm_c += p_VBr->dVbrD_Fa(kappa_n, k, Fm, Fb, X_beta, Y_beta);
        dVFm_c += p_VBr->dVbrX_Fa(kappa_n, k, Fm, Fb, X_beta, Y_beta);
      }
    }
#pragma omp critical(sum_X_core)
    { dVFm += dVFm_c; }
  }

  return dVFm;
}

} // namespace HF
