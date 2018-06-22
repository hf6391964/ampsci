#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include "ContinuumOrbitals.h"
#include "HF_hartree.h"
#include "WIG_369j.h"
#include <gsl/gsl_sf_bessel.h>
// #include "dmAbsorptionFunctions.h"


double Ckk(int ka, int kb)
{
  double ja = 0.5*ATI_twoj_k(ka);
  double jb = 0.5*ATI_twoj_k(kb);
  if(fabs(ja-jb)>1) return 0; //already checked, but for safety
  int la = ATI_l_k(ka);
  int lb = ATI_l_k(kb);
  if(abs(la-lb)!=1) return 0;
  double tjs = WIG_3j(ja,jb,1.,-0.5,0.5,0.);
  return (2*ja+1.)*(2*jb+1.)*pow(tjs,2);
}



//******************************************************************************
int main(void){

  struct timeval start, end;
  gettimeofday(&start, NULL);

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha;// for non-relativistic approx
  double hart_del;//HART convergance

  std::vector<std::string> str_core; //States for the core

  //Green potential parameters
  int Gf;
  double Gh,Gd;

  // Max anglular momentums
  int n_min=0, n_max=10;
  int l_min=0, l_max=0;
  double hw_min, hw_max; //hbar-omega - incident energy (or m_A).
  int N_hw;

  std::string label; //label for output file

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmAbsorption.in"); //XXX set up so can use other inputs!
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    while(true){
      std::string str;
      ifs >> str;
      if(str=="0"||str=="."||str=="|"||str=="!") break;
      str_core.push_back(str);
    }
    getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;         getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;            getline(ifs,jnk);
    ifs >> hart_del;                  getline(ifs,jnk);
    ifs >> n_min >> n_max;            getline(ifs,jnk);
    ifs >> l_min >> l_max;            getline(ifs,jnk);
    ifs >> hw_min >> hw_max >> N_hw;  getline(ifs,jnk);
    ifs >> varalpha;                  getline(ifs,jnk);
    ifs >> label;                     getline(ifs,jnk);
    ifs.close();
  }

  //convert hw from eV to a.u.
  hw_min/=HARTREE_EV;
  hw_max/=HARTREE_EV;
  //allow easier single energy
  if(N_hw==1) hw_max=hw_min;
  //don't let the max l be too large (for S.B. look-up table)
  if(l_max>=n_max) l_max = n_max-1;

  //default Hartree convergance goal:
  if(hart_del==0) hart_del=1.e-6;

  //alpha can't be zero:
  if(varalpha==0) varalpha=1.e-25;

  //Look-up atomic number, Z, and also A
  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI_a[Z]; //if none given, get default A

  printf("\nRunning for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");


  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n"
  , wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  // Check if 'h' is small enough for oscillating region:
  double h_target = (M_PI/15)/sqrt(2.*hw_max);
  if(wf.h>2*h_target){
    std::cout<<"\nWARNING 101: Grid not dense enough for contimuum state with "
      <<"ec="<<hw_max<<" (h="<<wf.h<<", need h<"<<h_target<<")\n";
    std::cout<<"Program will continue, but continuum wfs may be bad.\n\n";
  }

  //If non-zero A is given, use spherical nucleus.
  if(A>0) wf.sphericalNucleus();


  if(str_core.size()!=0){
    if(Gf!=0){
      if(Gh==0) PRM_defaultGreen(Z,Gh,Gd);
      printf("Using Green potential: H=%.4f  d=%.4f\n",Gh,Gd);
    }
    else printf("Using Hartree potential (converge to %.0e)\n",hart_del);
    //Determine which states are in the core:
    int core_ok = wf.determineCore(str_core);
    if(core_ok==2){
      std::cout<<"Problem with core: ";
      for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
      std::cout<<"\n";
      return 1;
    }
  }

  //Solve Hartree potential for core.
  //NOTE have option for H-like !! XXX
  if(Gf==0){
    HF_hartreeCore(wf,hart_del);
  }else{
    //Fill the electron part of the (local/direct) potential
    for(int i=0; i<wf.ngp; i++){
      wf.vdir.push_back(PRM_green(Z,wf.r[i],Gh,Gd));
    }
    // Solve Dirac equation for each (bound) core state:
    wf.solveInitialCore();
  }

  //Find max 'n' in the core (used for valence energy guess)
  int max_core_n=0;
  for(int i=0; i<wf.num_core; i++)
    if(wf.nlist[i]>max_core_n) max_core_n = wf.nlist[i];

  //Calculate the valence (and excited) states
  for(int n=n_min; n<=n_max; n++){
    for(int l=l_min; l<=l_max; l++){ //loop over l
      if(l+1>n) continue;
      for(int tk=0; tk<2; tk++){ //loop over k (ie j=l +/- 1/2)
        int k;
        if(tk==0) k=l;      //j = l - 1/2
        else      k=-(l+1); //j = l + 1/2
        if(k==0) continue;  // no j = l - 1/2 for l=0
        if(wf.isInCore(n,k)) continue; //skip states already in the core
        //This energy guess works very well for Cs, Fr etc.
        //Works poorly (but still converges) for light atoms
        int dn=n-max_core_n;
        double neff=1.+dn;
        double x=1;
        if(max_core_n<4) x=0.25;
        if(l==1) neff+=0.5*x;
        if(l==2) neff+=2.*pow(x,0.5);
        if(l>=3) neff+=4.*x;
        double en_a = -0.5/pow(neff,2);
        wf.solveLocalDirac(n,k,en_a);
      }
    }
  }

  //XXX make en_guess(n,l,nmaxcore)
  //XXX make "solveValenceStates"

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);



  //Output energies:
  printf("\n n  l_j      k   inf  it    eps      En (au)      En (eV)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    if((int)m==wf.num_core)
      std::cout<<" ========= Valence: ======\n"
      << " n  l_j      k   inf  it    eps      En (au)      En (eV)\n";
    int i = sort_list[m];
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = ATI_twoj_k(k);
    int l = ATI_l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %2s_%2i/2 %3i  %4.0f %3i  %5.0e  %11.5f %12.4f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*HARTREE_EV);
  }

  //Continuum wavefunction object
  ContinuumOrbitals cntm(wf);

  std::vector<std::string> state_list;
  std::vector<double> Ink_list;
  std::vector<int> kap_list;






  double factor = (4./3.)*pow(M_PI,2)*ALPHA;

  // (1) Calculate each rad_int, as a function of energy
  // (2) There is one for each state, summed over k' - just sum, too hard o/wise
  // (3) If summing, need include the coef!
  // (4) Only need a few (2?) l' for each l ! don't waste time solving cntm!
  // (5) But, for each core & each hw, need re-solve cntm => order matters??


  double ecmax = 15.;

  std::vector< std::vector<double> > K_hw_nk;

  for(int ihw=0; ihw<N_hw; ihw++){
    std::vector<double> K_nk;
    double x  = double(ihw)/(N_hw-1.);
    double hw = hw_min*pow(hw_max/hw_min,x);
    for(size_t is=0; is<wf.nlist.size(); is++){
      //Only calculate for required states..
      int k = wf.klist[is];
      int l = ATI_l_k(k);
      int n = wf.nlist[is];
      int twoj = ATI_twoj_k(k);
      if(n<n_min || n>n_max || l<l_min || l>l_max) continue;
      std::string state=std::to_string(n)+ATI_l(l)
        +"_{"+std::to_string(twoj)+"/2}";
      state_list.push_back(state);
      Ink_list.push_back(-wf.en[is]);
      kap_list.push_back(k);
      double ef = hw + wf.en[is];
      if(ef<=0 || ef>ecmax){
        K_nk.push_back(0);
        continue;
      }
      cntm.clear(); //clear previous
      int lc1 = l-1;
      int lc2 = l+1;
      cntm.solveLocalContinuum(ef,lc1,lc1); //XXX note: solve for 1 too many..
      cntm.solveLocalContinuum(ef,lc2,lc2);
      double tmp_K=0;
      for(size_t isc=0; isc<cntm.klist.size(); isc++){
        int kc = cntm.klist[isc];
        int lc = ATI_l_k(kc);
        int twojc = ATI_twoj_k(kc);
        if(abs(twojc-twoj)>2) continue; //dont need check l, calcd only right 1s
        double ckk = Ckk(kc,k);
        if(ckk==0) continue;
        double tmp_Rint=0;
        for(int ir=0; ir<wf.ngp; ir++){
          tmp_Rint += (wf.p[is][ir]*cntm.p[isc][ir]
                      + wf.q[is][ir]*cntm.q[isc][ir])
                      *wf.r[ir]*wf.drdt[ir];
        }
        tmp_Rint *= wf.h;
        tmp_K += pow(tmp_Rint,2)*ckk;
        //std::cout<<n<<k<<" "<<hw<<" "<<ef<<" "<<tmp_K<<"\n";
      }
      K_nk.push_back(tmp_K);
    }
    K_hw_nk.push_back(K_nk);
  }






  // //nb: only makes sense to sum the core!!
  // std::vector<float> s_core(N_hw);
  // for(int ihw=0; ihw<N_hw; ihw++){
  //   for(size_t is=0; is<state_list.size(); is++) s_core[ihw]+=s_nk_hw[is][ihw];
  // }
  //
  // //Write out to text file (in gnuplot friendly form)
  // fname = "sigma-"+Z_str+"_"+label+".txt";
  // ofile.open(fname);
  // ofile<<"hw/au hw/eV Total ";
  // for(size_t is=0; is<state_list.size(); is++){
  //   int k=kap_list[is];
  //   if(k!=-1 && k<0) continue;
  //   ofile<<state_list[is]<<" ";
  // }
  // ofile<<"\n";
  // for(int ihw=0; ihw<N_hw; ihw++){
  //   double x=ihw/(N_hw-1.);
  //   double hw = hw_min*pow(hw_max/hw_min,x);
  //   ofile<<hw<<" "<<hw*HARTREE_EV<<" ";
  //   ofile<<s_hw[ihw]<<" ";
  //   for(size_t is=0; is<state_list.size(); is++){
  //     int k=kap_list[is];
  //     if(k!=-1 && k<0) continue;
  //     ofile<<s_nk_hw[is][ihw]<<" ";
  //   }
  //   ofile<<"\n";
  // }
  // ofile.close();



  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.1f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.1f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}




// //pre-calculate the spherical Bessel function look-up table for efficiency
// //XXX Only for the Born approximation
//
// std::vector< std::vector< std::vector<float> > > jLqr;
// std::cout<<std::endl;
// int max_L = l_max+1; //+1 for g_nk case
// int min_L = l_min-1;
// if(min_L<0) min_L=0;
// jLqr.resize(max_L-min_L+1, std::vector< std::vector<float> >
//   (N_hw, std::vector<float>(wf.ngp)));
// for(int L=min_L; L<=max_L; L++){
//   std::cout<<"\rCalculating spherical Bessel look-up table for L="
//     <<L<<"/"<<max_L<<" .. "<<std::flush;
//   #pragma omp parallel for
//   for(int ik=0; ik<N_hw; ik++){
//     double x=ik/(N_hw-1.);
//     double ef = hw_min*pow(hw_max/hw_min,x); //note: not hw!!
//     double k = sqrt(2.*ef);
//     for(int ir=0; ir<wf.ngp; ir++){
//       jLqr[L-min_L][ik][ir] = gsl_sf_bessel_jl(L, k*wf.r[ir]);
//     }
//   }
// }
// std::cout<<"done\n";
//
//
// std::cout<<"Calculating Chi_nk^2(k) ..";
// std::vector< std::vector<float> > chi2_nk_k;
// std::vector<std::string> state_list;
// std::vector<double> Ink_list;
// std::vector<int> kap_list;
// for(size_t is=0; is<wf.nlist.size(); is++){
//   //Only calculate for required states..
//   int k = wf.klist[is];
//   int l = ATI_l_k(k);
//   int n = wf.nlist[is];
//   int twoj = ATI_twoj_k(k);
//   if(n<n_min || l<l_min || l>l_max) continue;
//   if(n>n_max) continue;
//   std::string state=std::to_string(n)+ATI_l(l)+"_{"+std::to_string(twoj)+"/2}";
//   state_list.push_back(state);
//   Ink_list.push_back(-wf.en[is]);
//   kap_list.push_back(k);
//   //calculate chi_nk^2(k) for given nk
//   std::vector<float> chi2_k;
//   for(int ik=0; ik<N_hw; ik++){
//     //perform integral over r:
//     double fint=0, gint=0;
//     int ltil = ATI_l_k(-k);
//     for(int ir=0; ir<wf.ngp; ir++){
//       fint += wf.p[is][ir]*wf.r[ir]*jLqr[l-min_L][ik][ir]*wf.drdt[ir];
//       gint += wf.q[is][ir]*wf.r[ir]*jLqr[ltil-min_L][ik][ir]*wf.drdt[ir];
//     }
//     double tmp_chi = (pow(fint,2)+pow(gint,2))*pow(wf.h,2);
//     chi2_k.push_back(tmp_chi);
//   }
//   chi2_nk_k.push_back(chi2_k);
// }
// std::cout<<".. done!\n";
// double q_to_keV = (HARTREE_EV*CLIGHT)/1.e3;
//
// // XXX do for core + 'not' core ?? nah, prob don't need
// //Write out to text file (in gnuplot friendly form)
// std::ofstream ofile;
// std::string fname = "chi2-"+Z_str+"_"+label+".txt";
// ofile.open(fname);
// ofile<<"k/au k/eV en_f/eV ";
// for(size_t is=0; is<state_list.size(); is++)
//   ofile<<state_list[is]<<" ";
// ofile<<"\n";
// for(int ik=0; ik<N_hw; ik++){
//   double x=ik/(N_hw-1.);
//   double ef = hw_min*pow(hw_max/hw_min,x); //note: not hw!!
//   double k = sqrt(2.*ef);
//   ofile<<k<<" "<<k*q_to_keV<<" "<<ef*HARTREE_EV<<" ";
//   for(size_t is=0; is<state_list.size(); is++)
//     ofile<<chi2_nk_k[is][ik]<<" ";
//   ofile<<"\n";
// }
// ofile.close();
//
// //XXX if "is_in_core" XXX
// double xfill = 0.5; //XXX fr core, hard-coded!
//
// double factor = (16./3.)*M_PI*ALPHA;
//
// std::vector< std::vector<float> > s_nk_hw;
// s_nk_hw.resize(state_list.size(), std::vector<float>(N_hw)); //initialise
// for(int ihw=0; ihw<N_hw; ihw++){
//   double x=ihw/(N_hw-1.);
//   double hw = hw_min*pow(hw_max/hw_min,x);
//   for(size_t is=0; is<state_list.size(); is++){
//     double ef = hw - Ink_list[is];
//     if(ef<=0) continue;
//     double dj = (N_hw-1.)*log(ef/hw_min)/log(hw_max/hw_min);
//     int j = int(dj);
//     double B = dj-j;
//     double A = 1-B; //A of j, B of (j+1)
//     if(j<0 || j>= N_hw-1) std::cout<<"ERROR 308: j"<<j<<"\n";
//     double k = sqrt(2*ef);
//     double tmp_chi = A*chi2_nk_k[is][j] + B*chi2_nk_k[is][j+1];
//     double coef = factor * fabs(kap_list[is]) * xfill;
//     double sig = coef * (pow(k,3)/hw) * tmp_chi;
//     s_nk_hw[is][ihw]=sig;
//   }
// }
