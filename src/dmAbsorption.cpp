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

//Declare angular coeficient. [See Phys.Rev.D 93, 115037 (2016).]
double CLkk(int L, int ka, int kb);

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
  int min_n, max_n;
  int min l, max_l;
  int min lc, max lc;

  std::string label; //label for output file

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmAbsorption.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    while(true){
      std::string str;
      ifs >> str;
      if(str=="."||str=="|"||str=="!") break;
      str_core.push_back(str);
    }
    getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> hart_del;              getline(ifs,jnk);
    ifs >> varalpha;              getline(ifs,jnk);
    ifs >> label;                 getline(ifs,jnk);
    ifs.close();
  }

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
  if(Gf!=0) printf("Using Green potential: H=%.4f  d=%.4f\n",Gh,Gd);
  else printf("Using Hartree potential (converge to %.0e)\n",hart_del);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n"
  , wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  // Check if 'h' is small enough for oscillating region:
  //XXX check this - dE max??
  double h_target = (M_PI/15)/sqrt(2.*demax);
  if(wf.h>2*h_target){
    std::cout<<"\nWARNING 101: Grid not dense enough for contimuum state with "
      <<"ec="<<demax<<" (h="<<wf.h<<", need h<"<<h_target<<")\n";
    std::cout<<"Program will continue, but continuum wfs may be bad.\n\n";
  }

  //If non-zero A is given, use spherical nucleus.
  if(A>0) wf.sphericalNucleus();

  //Determine which states are in the core:
  int core_ok = wf.determineCore(str_core);
  if(core_ok==2){
    std::cout<<"Problem with core: ";
    for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
    std::cout<<"\n";
    return 1;
  }

  //Solve Hartree potential for core.
  //NOTE have option for H-like !!
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

  //XXX solve for the valence states here!

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n n l_j    k Rinf its    eps     En (au)     En (/cm)    En (eV)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = ATI_twoj_k(k);
    int l = ATI_l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*HARTREE_ICM, eni*HARTREE_EV);
  }

  //Continuum wavefunction object
  ContinuumOrbitals cntm(wf);

  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AK; //float ok?
  std::vector<float> qlst(qsteps);
  std::vector<float> dElst;
  std::vector<std::string> nklst;

  //pre-calculate the spherical Bessel function look-up table for efficiency
  //XXX Only for the Born approximation
  std::cout<<std::endl;
  std::vector< std::vector< std::vector<float> > > jLqr_f;
  jLqr_f.resize(max_L+1, std::vector< std::vector<float> >
    (qsteps, std::vector<float>(wf.ngp)));
  for(int L=0; L<=max_L; L++){
    std::cout<<"\rCalculating spherical Bessel look-up table for L="
    <<L<<"/"<<max_L<<" .. "<<std::flush;
    #pragma omp parallel for
    for(int iq=0; iq<qsteps; iq++){
      double x=iq/(qsteps-1.);
      double q = qmin*pow(qmax/qmin,x);
      for(int ir=0; ir<wf.ngp; ir++){
        jLqr_f[L][iq][ir] = gsl_sf_bessel_jl(L, q*wf.r[ir]);
      }
    }
  }
  std::cout<<"done\n";



//XXX calc cross-sections here!











  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.1f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.1f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}



//******************************************************************************
double CLkk(int L, int ka, int kb)
/*
Calculates the angular coeficient (averaged over all m)
B. M. Roberts, V. A. Dzuba, V. V. Flambaum, M. Pospelov, and Y. V. Stadnik,
Phys. Rev. D 93, 115037 (2016). [arXiv:1604.04559]
*/
{
  int two_ja = ATI_twoj_k(ka);
  int two_jb = ATI_twoj_k(kb);
  double ja = 0.5*two_ja;
  double jb = 0.5*two_jb;
  int la = ATI_l_k(ka);
  int lb = ATI_l_k(kb);

  double tjB = WIG_3j(jb,L,ja,-0.5,0,0.5);
  if(fabs(tjB)==0) return 0;
  double B = 1./pow(tjB,2);

  //(-1)^(ja etc) -> calc sign
  int s1 = -1;
  if((two_ja+two_jb-2*(la+lb))%4==0) s1=1;

  double tj1 = WIG_3j(lb,la,L,0,0,0);
  double A = (1./4)*s1*(2*L+1)*pow(tj1,2);
  double X = s1*(two_ja+1)*(two_jb+1)*pow(tj1,2);
  double tj2 = WIG_3j(lb,la,L,-1,1,0);
  double Y = 8*sqrt(la*(la+1)*lb*(lb+1))*tj1*tj2;
  double Z = -4*(ka+1)*(kb+1)*pow(tj2,2);

  return (A*B)*(X+Y+Z);
}
