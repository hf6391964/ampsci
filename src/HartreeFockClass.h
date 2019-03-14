#pragma once
#include <string>
#include <vector>
class ElectronOrbitals;
class DiracSpinor;
class Grid;

/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.

Note: can run without exchange (Hartree) - for tests only.
(Still, non-local potential, 'vex' different for each core state)

XXX Have option to give a list of valence states!
Can solve them to some degree in parallel
Requires re-writing the valence part (a little)

  Definitions:
  v^k_ab(r)   := Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
  rho(r')     := fa(r')*fb(r') + ga(r')gb(r')
  Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2 * parity(la+lb+k)
  vex[a]      := [v_ex*psi_a](r) *(psi_a/psi_a^2) (approx exchange)

*/

class HartreeFock {

public:
  HartreeFock(ElectronOrbitals &wf, const std::string &in_core,
              double eps_HF = 1.e-8, bool in_ExcludeExchange = false);

  void solveValence(int n, int kappa);

  double calculateCoreEnergy() const;

private:
  ElectronOrbitals *const p_wf;
  const Grid *const p_rgrid;

  double m_eps_HF = 1.e-8;

  const int MAX_HART_ITS = 64;
  const bool m_excludeExchange; // for testing

  // const int m_ngp;
  size_t m_num_core_states;
  std::vector<int> twoj_list;
  std::vector<int> kappa_index_list;
  int m_max_kappa_index_so_far;

  // The "localised"/approximate HF potential:
  std::vector<std::vector<double>> vex;

  // Store underlying arrays. These are 'double private'
  // nb: these are 'highly non-rectangular'
  std::vector<std::vector<std::vector<std::vector<double>>>> m_arr_v_abk_r;
  std::vector<std::vector<std::vector<double>>> m_arr_Lambda_nmk;

private:
  void hartree_fock_core();
  void starting_approx_core(const std::string &in_core);

  void form_core_Lambda_abk();
  void extend_Lambda_abk(int kappa_a);
  double get_Lambda_abk(size_t a, size_t b, int k) const;

  void initialise_m_arr_v_abk_r_core();
  void extend_m_arr_v_abk_r_valence(int kappa_a);

  int index_from_kappa(int ka) const;
  int twoj_from_index(int i) const;
  int kappa_from_index(int i);
  int l_from_index(int i) const;

  void form_vabk_core();
  void form_vbb0();
  void calculate_v_abk(const DiracSpinor &phi_a, const DiracSpinor &phi_b,
                       int k, std::vector<double> &vabk) const;

  const std::vector<double> &get_v_aa0(size_t a) const;
  const std::vector<std::vector<double>> &get_v_abk(size_t a, size_t b) const;

  void form_vabk_valence(size_t w);

  void form_vdir(std::vector<double> &vdir, bool re_scale = false);

  void form_approx_vex_core(std::vector<std::vector<double>> &vex);

  void form_approx_vex_a(size_t a, std::vector<double> &vex_a);
};
