#pragma once
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/QkTable.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <string>

namespace UnitTest {
bool QkTable(std::ostream &obuff) {
  bool pass = true;

  Wavefunction wf({1000, 1.0e-6, 100.0, 10.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
  wf.formBasis({"8spdfghi", 30, 9, 1.0e-4, 1.0e-6, 40.0, false});

  // Split basis into core/excited
  std::vector<DiracSpinor> core, excited;
  for (const auto &Fb : wf.basis) {
    if (wf.isInCore(Fb.n, Fb.k)) {
      core.push_back(Fb);
    } else {
      excited.push_back(Fb);
    }
  }

  const auto max_tj = DiracSpinor::max_tj(wf.basis);

  std::mt19937 gen(0.0); // seeded w/ constant; same each run
  using uniform_int = std::uniform_int_distribution<std::size_t>;
  uniform_int r_b_index(0, wf.basis.size() - 1);
  uniform_int r_c_index(0, core.size() - 1);
  uniform_int r_e_index(0, excited.size() - 1);

  {
    Coulomb::QkTable Qk(&wf.basis);
    int tries = 0;
    double worst_q = 0.0;
    std::string worst_abcdk{};
    {
      while (++tries < 1000000) {
        const auto &Fa = wf.basis[r_b_index(gen)];
        const auto &Fb = wf.basis[r_b_index(gen)];
        const auto &Fc = wf.basis[r_b_index(gen)];
        const auto &Fd = wf.basis[r_b_index(gen)];
        for (int k = 0; k <= max_tj; ++k) {
          const auto qk1 = Qk.get_Qk(Fa, Fb, Fc, Fd, k);
          const auto qk2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
          if (qk1 == 0.0 && qk2 == 0.0)
            continue;
          const auto del = std::abs(qk1 - qk2);
          if (del > worst_q) {
            worst_q = del;
            worst_abcdk = Fa.shortSymbol() + Fb.shortSymbol() +
                          Fc.shortSymbol() + Fd.shortSymbol() + "_" +
                          std::to_string(k);
          }
        }
      }
    }
    std::cout << worst_abcdk << " " << worst_q << "\n";
    pass &= qip::check_value(&obuff, "Qk_abcd " + worst_abcdk, worst_q, 0.0,
                             1.0e-7);
  }

  {
    Coulomb::QkTable Qk(&core, &excited, &core, &excited);
    int tries = 0;
    double worst_q = 0.0;
    std::string worst_abcdk{};
    {
      while (++tries < 1000000) {
        const auto &Fa = core[r_c_index(gen)];
        const auto &Fb = excited[r_e_index(gen)];
        const auto &Fc = core[r_c_index(gen)];
        const auto &Fd = excited[r_e_index(gen)];
        for (int k = 0; k <= max_tj; ++k) {
          const auto qk1 = Qk.get_Qk(Fa, Fb, Fc, Fd, k);
          const auto qk2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
          if (qk1 == 0.0 && qk2 == 0.0)
            continue;
          const auto del = std::abs(qk1 - qk2);
          if (del > worst_q) {
            worst_q = del;
            worst_abcdk = Fa.shortSymbol() + Fb.shortSymbol() +
                          Fc.shortSymbol() + Fd.shortSymbol() + "_" +
                          std::to_string(k);
          }
        }
      }
    }
    std::cout << worst_abcdk << " " << worst_q << "\n";
    pass &= qip::check_value(&obuff, "Qk_cece " + worst_abcdk, worst_q, 0.0,
                             1.0e-7);
  }

  {
    Coulomb::QkTable Qk(&excited, &excited, &core, &core);
    int tries = 0;
    double worst_q = 0.0;
    std::string worst_abcdk{};
    {
      while (++tries < 1000000) {
        const auto &Fa = excited[r_e_index(gen)];
        const auto &Fb = excited[r_e_index(gen)];
        const auto &Fc = core[r_c_index(gen)];
        const auto &Fd = core[r_c_index(gen)];
        for (int k = 0; k <= max_tj; ++k) {
          const auto qk1 = Qk.get_Qk(Fa, Fb, Fc, Fd, k);
          const auto qk2 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
          if (qk1 == 0.0 && qk2 == 0.0)
            continue;
          const auto del = std::abs(qk1 - qk2);
          if (del > worst_q) {
            worst_q = del;
            worst_abcdk = Fa.shortSymbol() + Fb.shortSymbol() +
                          Fc.shortSymbol() + Fd.shortSymbol() + "_" +
                          std::to_string(k);
          }
        }
      }
    }
    std::cout << worst_abcdk << " " << worst_q << "\n";
    pass &= qip::check_value(&obuff, "Qk_eecc " + worst_abcdk, worst_q, 0.0,
                             1.0e-7);
  }

  return pass;
}

} // namespace UnitTest
