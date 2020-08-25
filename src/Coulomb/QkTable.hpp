#pragma once
#include "Angular/Angular_tables.hpp"
// #include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <utility>
#include <vector>
class Grid;
class DiracSpinor;
#include "Coulomb/YkTable.hpp"
#include "Physics/AtomData.hpp"

#include "Wavefunction/DiracSpinor.hpp"

namespace Coulomb {

class QkTable {
  using Real = double; // make template
  using Nk = AtomData::nkappa;
  using Nkabcd = std::array<Nk, 4>;

  // const std::vector<DiracSpinor> *const as, bs, cs, ds;

  std::vector<std::vector<Real>> m_Rabcd_k{};
  std::vector<Nkabcd> m_keys{};
  // m_Rabcd_k should be always same size as m_keys!

public:
  QkTable(const std::vector<DiracSpinor> *as,
          const std::vector<DiracSpinor> *bs = nullptr,
          const std::vector<DiracSpinor> *cs = nullptr,
          const std::vector<DiracSpinor> *ds = nullptr) {
    assert(as != nullptr);
    if (!bs)
      bs = as;
    if (!cs)
      cs = bs;
    if (!ds)
      ds = cs;

    YkTable ykbd(as->front().rgrid, bs, ds); // store this?

    // First loop through and 'size' the arrays, without calculating anything.
    // Since, we can't do this part in parallel!
    for (const auto &Fa : *as) {
      for (const auto &Fb : *bs) {
        for (const auto &Fc : *cs) {
          for (const auto &Fd : *ds) {
            // 1. Check if already calc'd
            const auto abcd = permute(Fa.nk(), Fb.nk(), Fc.nk(), Fd.nk());
            const auto pRk = find_Rk(abcd);
            if (pRk != nullptr)
              continue;
            // If so, continue.
            // If not, add new key
            m_keys.push_back(abcd);
            // and R vector place
            auto rk = m_Rabcd_k.emplace_back();
            const auto [kmin, kmax] = minmaxk(Fa.k, Fb.k, Fc.k, Fd.k);
            if (kmax < kmin)
              continue;
            const auto num_ks = std::size_t((kmax - kmin) / 2) + 1;
            rk.reserve(num_ks); // re-size?
          }
        }
      }
    }

    assert(m_Rabcd_k.size() == m_keys.size());

    /*
Fails here:

    5p_1/2,5p_3/2,5p_1/2,5p_3/2, 2
    k=0 2.82843
    k=1 0
    k=2 0
    k=3 0
    k=4 -0
    k=5 0

    */

    for (std::size_t i = 0; i < m_Rabcd_k.size(); ++i) {
      const auto &abcd = m_keys[i];
      auto &Rk = m_Rabcd_k[i];
      assert(Rk.empty());
      // 2. Find min/max k + store
      // 3. Calculate R for each k
      const auto [kmin, kmax] = minmaxk(abcd);
      for (auto j = 0ul; j < 4; ++j)
        std::cout << abcd[j].symbol() << ",";
      for (int k = kmin; k <= kmax; k += 2) {
        // Calculate R for each k!
        std::cout << " " << k;
      }
      std::cout << "\n";
      for (int k = 0; k <= kmax + 3; k++) {
        std::cout << "k=" << k << " "
                  << Angular::Ck_kk(k, abcd[0].k, abcd[2].k) *
                         Angular::Ck_kk(k, abcd[1].k, abcd[3].k)
                  << "\n";
      }
    }
  }

  //****************************************************************************
  static std::pair<int, int> minmaxk(const Nkabcd &abcd) {
    return minmaxk(abcd[0].k, abcd[1].k, abcd[2].k, abcd[3].k);
  }
  static std::pair<int, int> minmaxk(int ka, int kb, int kc, int kd) {
    // min of (a,c), (b,d)
    const auto tja = Angular::twoj_k(ka);
    const auto tjb = Angular::twoj_k(kb);
    const auto tjc = Angular::twoj_k(kc);
    const auto tjd = Angular::twoj_k(kd);
    const auto minK = Angular::min_lambda_tj(tja, tjb, tjc, tjd);
    const auto maxK = Angular::max_lambda_tj(tja, tjb, tjc, tjd);

    // Only every second k is non-zero. Only store those
    // nb: true for Q, not R technically
    const auto la = Angular::l_k(ka);
    const auto lb = Angular::l_k(kb);
    const auto lc = Angular::l_k(kc);
    const auto ld = Angular::l_k(kd);
    // parity: Check if minK is a valid k for both:
    const auto p1 = Angular::parity(la, lc, minK) == 1;
    const auto p2 = Angular::parity(lb, ld, minK) == 1;

    if (p1 && p2)
      return {minK, maxK};

    if (!p1 && !p2)
      return {minK + 1, maxK};

    // If one is and one isn't OK, then there are NO k's for which this will
    // work!
    return {0, -1};
  }

  // //****************************************************************************
  // Real get_Rk(int k, const Nk &a, const Nk &b, const Nk &c, const Nk &d) {
  //   auto [kmin, Rp] = lookup_R(permute(a, b, c, d));
  // }
  //
  // //****************************************************************************
  // std::pair<int, std::vector<Real> *> lookup_R(const Nkabcd &abcd);

  //****************************************************************************
  std::vector<Real> *find_Rk(const Nkabcd &abcd) {

    auto compare = [&abcd](const Nkabcd &el) {
      return abcd[0] == el[0] && abcd[1] == el[1] && abcd[2] == el[2] &&
             abcd[3] == el[3];
    };

    auto it = std::find_if(cbegin(m_keys), cend(m_keys), compare);
    if (it == cend(m_keys))
      return nullptr;
    const auto index = std::size_t(std::distance(cbegin(m_keys), it));
    assert(index < m_Rabcd_k.size());
    return &(m_Rabcd_k[index]);
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
    // XXX Check - what happens when 'min' is not unique?
    // Must still be fine
  }
}; // namespace Coulomb

} // namespace Coulomb
