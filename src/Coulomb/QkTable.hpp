#pragma once
#include "Angular/Angular_tables.hpp"
// #include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
class Grid;
class DiracSpinor;
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace Coulomb {

class QkTable {
  using Real = float; // make template
  using Nk = AtomData::nkappa;
  using Nkabcd = std::array<Nk, 4>;

  // map:
  std::map<Nkabcd, std::vector<Real>> m_Rmap_k{};
  // store vector of pointers into map, so that I can random access the map
  std::vector<std::pair<const Nkabcd, std::vector<Real>> *> m_pQk{};
  Angular::Ck_ab m_Ck{};
  Angular::SixJ m_6j{};

  const bool b_is_d;
  const bool a_is_c;
  const bool c_is_d; // required for exchange
  const bool all_same;

public:
  QkTable(const std::vector<DiracSpinor> *as,
          const std::vector<DiracSpinor> *bs = nullptr,
          const std::vector<DiracSpinor> *cs = nullptr,
          const std::vector<DiracSpinor> *ds = nullptr)
      : b_is_d(!bs || bs == ds),
        a_is_c(!cs || as == cs),
        c_is_d(!ds || cs == ds),
        all_same(b_is_d && a_is_c && (!bs || as == bs)) {

    assert(as != nullptr);
    if (!bs)
      bs = as;
    if (!cs)
      cs = bs;
    if (!ds)
      ds = cs;

    // Find maximum (2j), and fill the 3j and 6j symbol tables:
    const auto mtja = DiracSpinor::max_tj(*as);
    const auto mtjb = DiracSpinor::max_tj(*bs);
    const auto mtjc = DiracSpinor::max_tj(*cs);
    const auto mtjd = DiracSpinor::max_tj(*ds);
    const auto max_tj = std::max({mtja, mtjb, mtjc, mtjd});
    m_Ck.fill_maxK_twojmax(max_tj, max_tj);
    m_6j.fill(max_tj, max_tj);

    // First loop through and 'size' the arrays, without calculating anything.
    // Since, we can't do this part in parallel!
    std::cout << "Creating Qk lookup table Q_abcd:\n";
    if (all_same) {
      std::cout << "a=b=c=d = " << DiracSpinor::state_config(*as) << "\n";
    } else {
      if (a_is_c)
        std::cout << "a=c, ";
      if (b_is_d)
        std::cout << "b=d, ";
      if (c_is_d)
        std::cout << "c=d, ";
      std::cout << "a = " << DiracSpinor::state_config(*as) << "\n";
      std::cout << "b = " << DiracSpinor::state_config(*bs) << "\n";
      std::cout << "c = " << DiracSpinor::state_config(*cs) << "\n";
      std::cout << "d = " << DiracSpinor::state_config(*ds) << "\n";
    }

    // Creat the map: (but don't fill it yet)
    // Creating must be in seriel - parallelise the Q calculation below
    {
      IO::ChronoTimer tt("Create");
      for (const auto &Fa : *as) {
        for (const auto &Fb : *bs) {
          if (all_same && Fb < Fa)
            continue;
          for (const auto &Fc : *cs) {
            if (a_is_c && Fc < Fa)
              continue;
            for (const auto &Fd : *ds) {
              if ((all_same && Fd < Fa) || (b_is_d && Fd < Fb))
                continue;
              // Not all combinations will allow non-zero results: skip those
              const auto [kmin, kmax] = minmaxk(Fa.k, Fb.k, Fc.k, Fd.k);
              if (kmax < kmin)
                continue;
              // permute puts in to a<{b,c,d} and b<d order
              const auto abcd = permute(Fa.nk(), Fb.nk(), Fc.nk(), Fd.nk());
              // adds if new perumation; does nothing if not
              m_Rmap_k[abcd];
            }
          }
        }
      }
    }

    // store vector of pointers into map, so that I can random access the map
    // ANY change to the map may invalidate these pointers (but changes to the
    // _values_ in the map are fine)
    // Random access into the map allows us to loop over the map in parallel
    m_pQk.reserve(m_Rmap_k.size());
    for (auto &pair : m_Rmap_k) {
      m_pQk.emplace_back(&pair);
    }

    // Fill the Table:
    YkTable ykbd(as->front().rgrid, bs, ds);
    {
      IO::ChronoTimer tt("fill");
#pragma omp parallel for
      for (auto i = 0ul; i < m_Rmap_k.size(); ++i) {
        auto &[key, Rk] = *m_pQk[i];
        const auto &[a, b, c, d] = key;
        // Look up in the orbitals
        const auto &Fa = *find_orb(a, as, bs, cs, ds);
        const auto &Fb = *find_orb(b, as, bs, cs, ds);
        const auto &Fc = *find_orb(c, as, bs, cs, ds);
        const auto &Fd = *find_orb(d, as, bs, cs, ds);

        assert(Rk.empty());
        const auto [kmin, kmax] = minmaxk(key);
        if (kmax >= kmin) {
          const auto num_ks = std::size_t((kmax - kmin) / 2) + 1;
          Rk.reserve(num_ks); // re-size?
        }
        for (int k = kmin; k <= kmax; k += 2) {
          // depending on {a,b,c,d}, may have y_bd or y_ac
          // Further, may be {bd} or {db} -- flaw in YkTable
          auto ybd = ykbd.ptr_yk_ab(k, Fb, Fd);
          if (ybd == nullptr)
            ybd = ykbd.ptr_yk_ab(k, Fd, Fb);
          if (ybd != nullptr) {
            Rk.push_back(static_cast<Real>(
                Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k, *ybd, m_Ck)));
          } else {
            auto yac = ykbd.ptr_yk_ab(k, Fa, Fc);
            if (yac == nullptr)
              yac = ykbd.ptr_yk_ab(k, Fc, Fa);
            assert(yac != nullptr);
            Rk.push_back(static_cast<Real>(
                Coulomb::Qk_abcd(Fb, Fa, Fd, Fc, k, *yac, m_Ck)));
          }
          assert(Rk.back() != Real(0.0));
        }
      }
    }
    std::cout << m_Rmap_k.size() << " Q_abcd integrals calculated (x k)\n";
  }

  //****************************************************************************
  static const DiracSpinor *find_orb(const Nk &nk,
                                     const std::vector<DiracSpinor> *as,
                                     const std::vector<DiracSpinor> *bs,
                                     const std::vector<DiracSpinor> *cs,
                                     const std::vector<DiracSpinor> *ds) {
    // look in {a}
    auto pFnk = std::find(cbegin(*as), cend(*as), nk);
    if (pFnk != cend(*as))
      return &*pFnk;
    // If not in a, look in b
    pFnk = std::find(cbegin(*bs), cend(*bs), nk);
    if (pFnk != cend(*bs))
      return &*pFnk;
    // if not in b, look in c
    pFnk = std::find(cbegin(*cs), cend(*cs), nk);
    if (pFnk != cend(*cs))
      return &*pFnk;
    // if not in c, must be in d
    pFnk = std::find(cbegin(*ds), cend(*ds), nk);
    if (pFnk != cend(*ds))
      return &*pFnk;
    // Could not find???
    assert(false);
    return nullptr;
  }

  //****************************************************************************
  static std::pair<int, int> minmaxk(const Nkabcd &abcd) {
    return minmaxk(abcd[0].k, abcd[1].k, abcd[2].k, abcd[3].k);
  }
  static std::pair<int, int> minmaxk(int ka, int kb, int kc, int kd) {
    // min/max allowed k for (a,k,c), (b,k,d)
    const auto tja = Angular::twoj_k(ka);
    const auto tjb = Angular::twoj_k(kb);
    const auto tjc = Angular::twoj_k(kc);
    const auto tjd = Angular::twoj_k(kd);
    const auto minK = Angular::min_k_CkCk(tja, tjb, tjc, tjd);
    const auto maxK = Angular::max_k_CkCk(tja, tjb, tjc, tjd);

    // Only every second k is non-zero. Only store those
    // nb: true for Q, not R technically
    const auto la = Angular::l_k(ka);
    const auto lb = Angular::l_k(kb);
    const auto lc = Angular::l_k(kc);
    const auto ld = Angular::l_k(kd);
    // parity: Check if minK is a valid k for both:
    const auto p1 = Angular::parity(la, lc, minK) == 1;
    const auto p2 = Angular::parity(lb, ld, minK) == 1;

    // XXX Also kmax?

    if (p1 && p2)
      return {minK, maxK};

    if (!p1 && !p2)
      return {minK + 1, maxK};

    // If one is and one isn't OK, then there are NO k's for which this will
    // work!
    return {0, -1};
  }

  //****************************************************************************
  // Returns the "smallest first" Coulomb symmetry permutation of
  // {a,b,c,d}
  std::array<Nk, 4> permute(const Nk &a, const Nk &b, const Nk &c,
                            const Nk &d) {
    // put smallest first
    const auto min = std::min({a, b, c, d});
    if (min == a) {
      // options are abcd, and adcb
      return (b < d) ? std::array{a, b, c, d} : std::array{a, d, c, b};
    } else if (min == b) {
      // options are badc, and bcda
      return (a < c) ? std::array{b, a, d, c} : std::array{b, c, d, a};
    } else if (min == c) {
      // options are cbad, and cdab
      return (b < d) ? std::array{c, b, a, d} : std::array{c, d, a, b};
    } else if (min == d) {
      // options are dabc, and dcba
      return (a < c) ? std::array{d, a, b, c} : std::array{d, c, b, a};
    }
    assert(false);
    // Check - what happens when 'min' is not unique?
    // Must still be fine
  }
};

} // namespace Coulomb
