#include "XEM_Constant.h"
extern "C"
{
  //  SUBROUTINE DEUT_ELASTIC(E0,PF_E,TSCAT,POL_BEAM,DEUT_SIG,DEUT_ASYM,
  //  #                            PN_HI,PT_HI,PL_HI,PN_HD,PT_HD,PL_HD); 
  void hyd_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*); 
  void deut_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*);
  void trit_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*); 
  void he3_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*);
  void he4_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*);
  void c12_elastic_(const double*,const double*, const double*, const double*, double*, double*,double*,double*, double*, double*, double*, double*);
  //       SUBROUTINE F1F2IN09(Z, A, QSQ, Wsq, F1, F2, rc)                       
  // void f1f2in09_(const double*,const double*, const double*, const double*, double*, double*, double*);
  //void f1f2qe09_(const double*,const double*, const double*, const double*, double*, double*);

}



inline double gCal_Elastic(const int kA, const int kZ, const double kEb, const double kDeg)         
{
/*
//--------------------------------------------------------------------
// calculate elastic cross section of H,2H,3H,3He,4He,12C
// fortran subroutines from MCEEP package
// build-in coulomb correction
// -Shujie, 2018.04
//--------------------------------------------------------------------
*/
  double dummy=0,xs=-1;
  

  
  if(kZ==1 && kA==1)           hyd_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else if(kZ==1 && kA==2)      deut_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else if(kZ==1 && kA==3)      trit_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else if(kZ==2 && kA==3)      he3_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else if(kZ==2 && kA==4)      he4_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else if(kZ==6 && kA==12)     c12_elastic_(&kEb,&dummy,&kDeg,&dummy,&xs,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  else xs=-1;
  //  cout<<"elastic:  "<<kA<<" "<<kZ<<"  "<<kEb<<"  "<<kDeg<<"  "<<xs<<endl;
  return xs;
}



