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
  using Real = double; // make template
  using Nk = AtomData::nkappa;
  using Nkabcd = std::array<Nk, 4>;

  /*
  #include "Coulomb/QkTable.hpp"
    Coulomb::QkTable(&core, &excited, &excited, &core);
    // Coulomb::QkTable(&excited, &core, &excited, &core);
    // Coulomb::QkTable(&wf.basis);
  */
  std::vector<std::vector<Real>> m_Rabcd_k{};
  std::vector<Nkabcd> m_keys{};
  std::map<Nkabcd, std::vector<Real>> m_Rmap_k{};
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

    m_keys.reserve(144789);
    m_Rabcd_k.reserve(144789);

    {
      IO::ChronoTimer t1("perm");
      IO::ChronoTimer t2("find");
      IO::ChronoTimer t3("rest");
      IO::ChronoTimer t4("tot");
      t1.stop();
      t2.stop();
      t3.stop();

      // First loop through and 'size' the arrays, without calculating anything.
      // Since, we can't do this part in parallel!
      for (const auto &Fa : *as) {
        for (const auto &Fb : *bs) {
          for (const auto &Fc : *cs) {
            for (const auto &Fd : *ds) {
              // 1. Check if already calc'd
              t1.start();
              const auto abcd = permute(Fa.nk(), Fb.nk(), Fc.nk(), Fd.nk());
              t1.stop();
              t2.start();
              // const auto pRk = find_Rk(abcd);
              const auto pRk = m_Rmap_k.find(abcd);
              t2.stop();
              if (pRk != m_Rmap_k.end())
                continue;
              // If so, continue.
              // If not, add new key and R vector place
              t3.start();
              // m_keys.push_back(abcd);
              // auto rk = m_Rabcd_k.emplace_back();
              m_Rmap_k.insert({abcd, {}});
              t3.stop();
            }
          }
        }
      }
    }

    YkTable ykbd(as->front().rgrid, bs, ds);
    // note: below, not guarenteed we won't need y_{ac} !!

    for (auto &[key, Rk] : m_Rmap_k) {
      // const auto &abcd = m_keys[i];
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
        auto ybd = ykbd.ptr_yk_ab(k, Fb, Fd);
        if (ybd == nullptr)
          ybd = ykbd.ptr_yk_ab(k, Fd, Fb);
        // depends if we calculated y_{ac} or y_{bd}
        if (ybd != nullptr) {
          Rk.push_back(Coulomb::Rk_abcd(Fa, Fc, *ybd));
        } else {
          auto yac = ykbd.ptr_yk_ab(k, Fa, Fc);
          if (yac == nullptr)
            yac = ykbd.ptr_yk_ab(k, Fc, Fa);
          assert(yac != nullptr);
          Rk.push_back(Coulomb::Rk_abcd(Fb, Fd, *yac));
        }
        auto cc = Angular::Ck_kk(k, a.k, c.k) * Angular::Ck_kk(k, b.k, d.k);
        std::cout << Fa.shortSymbol() << "," << Fb.shortSymbol() << ","
                  << Fc.shortSymbol() << "," << Fd.shortSymbol() << " " << k
                  << " " << cc * Rk.back() << "\n";
        if (cc == 0.0) {
          std::cout << " *** \n";
          std::cin.get();
        }
      }
    }
    std::cout << m_Rmap_k.size() << "\n";
    std::cin.get();
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
    // min of (a,c), (b,d)
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
