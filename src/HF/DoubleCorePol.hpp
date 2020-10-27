#pragma once
#include "ExternalField.hpp"
#include <string>
#include <vector>

namespace HF {

class DoubleCorePol {
public:
  DoubleCorePol(const ExternalField *const dV1, const ExternalField *const dV2);

private:
  std::vector<std::vector<DiracSpinor>> m_X{};
  std::vector<std::vector<DiracSpinor>> m_Y{};

  // const DiracOperator::TensorOperator *const m_h; //??
  const ExternalField *const m_dV1;
  const ExternalField *const m_dV2;
  const int m_rank;
  const int m_pi;
  const std::vector<DiracSpinor> *const p_core;
  const std::vector<double> m_vl;
  const std::vector<double> m_Hmag;
  const double m_alpha;
  const bool m_imag;
  const Breit *const p_VBr;
  double m_core_eps = 1.0;
  // Angular::SixJ m_6j; // used?

public:
  void solve_TDHFcore(const double omega, int max_its = 200,
                      const bool print = true);
  //
  // //! Returns eps (convergance) of last solve_TDHFcore run
  // double get_eps() const { return m_core_eps; }

  //! @brief Clears the dPsi orbitals (sets to zero)
  void clear_dPsi();

  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb, bool conj,
            const DiracSpinor *const Fexcl = nullptr) const;
  double dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  DiracSpinor dV_rhs(const int kappa_n, const DiracSpinor &Fm,
                     bool conj = false,
                     const DiracSpinor *const Fexcl = nullptr) const;

  const std::vector<DiracSpinor> &get_dPsis(const DiracSpinor &Fc,
                                            dPsiType XorY) const;
  const DiracSpinor &get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                const int kappa_x) const;

  double rme3js(const int twoja, const int twojb, int two_mb,
                int two_q = 0) const;

private:
  void initialise_dPsi();

public:
  DoubleCorePol &operator=(const DoubleCorePol &) = delete;
  DoubleCorePol(const DoubleCorePol &) = default;
  ~DoubleCorePol() = default;
};

} // namespace HF
