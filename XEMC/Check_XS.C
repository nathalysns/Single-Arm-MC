////////////////////////////////////////////////
// Generate Cross Section Look-up Table
//  --Zhihong Ye, 07/03/2012
////////////////////////////////////////////////
/*Include{{{*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <errno.h>
#include <sstream>
//#include <iterator>
//from root
#include <TMath.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TMatrix.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
/*}}}*/
using namespace::std;
using namespace::TMath;
#include "./SRC/XEMC.h"
//#include "/work/halla/e08014/disk1/yez/XEMC/SRC/XEMC.h"
//int wall = 1;//wall=0,1,2 target, entrance window, exit window
char* gTarget;
Double_t E0=-1.0;   //GeV
Double_t Ep=-1.0;   //GeV
Double_t Theta=-1.0; //Degree
Bool_t bAsk=kFALSE;
int getargs(int argc, char**argv);


int main(int argc,char** argv){
  string Head, Line;
  stringstream ss;

  int gerr = getargs(argc,argv);
  
  // input: grid of kinematics 
  char* id = argv[1]; // kin id 
  TString filename = Form("table_%s.inp",id);    
  TString Targets[7]={"Dummy","D2","H3","He3","Al_u","Al_d","C12"};
  
  int tarid = atoi(argv[2]);
  int Z, A;	
  double E0, Ep, Theta,Q2,xbj,sig;	
  
  
  ofstream outfile,outfile1;
  TString output = Form("extern.out.%s.%dx1",id,tarid);
  outfile.open(output.Data(), ios::trunc | ios::out);
  outfile << "   A "  << "    Z    "<< "XBJ" << "         " << "Q^2(GeV)" << "     " <<"deg"<< "     "<<"XS_QE"<<"      "<<"XS_DIS"<<"      "<<"XS_Born      XS_Rad"<<'\n';
  outfile.close();
  //TString output1= Form("/home/shujie/jlab/f1f217/%s_%d.xemc",id,tarid);
  TString output1 = Form("table.out.%s.%d",id,tarid);
  outfile1.open(output1.Data(), ios::trunc | ios::out);
  // outfile1 << "   A      Z     E0"  << "      "<< "Ep" << "      " << "theta   "<<"XS_Born   "<<"born/rad\n";
  outfile1.close();
  Double_t  cs_qe=0, cs_dis=0, cs_rad=0,cs_born=0,cs_el=0;
 
  //if (tarid>0) N=1;
  int i = 0;
  

  //  for (i=0;i<N;i++){
  ifstream infile;
  TString Target = Targets[tarid];

  //Basic input file, mostly important for radiated cross sections
  //  cout<< Target<<" target "<<Target.Data()<<endl;
  TString Target_Input = Form("input/%s_Input.dat",Target.Data());
  
  //TString filename = Form("/home/shujie/jlab/f1f217/%s.out",id);
  infile.open (filename.Data());
  // cout<<"open input file "<<filename.Data()<<endl;
  for(int qq=0;qq<6;qq++){
    getline(infile,Head); // get the header of the data.
  }

  // can't see why we need to call this 2017.8.29
  //	  XEM_TGT* XEMTarget= new XEM_TGT();   
  // XEM_TGT* XEMWall= new XEM_TGT();   
	 
  //	  XEMTarget->GetValueAZ(A,Z);
  // delete XEMTarget;//delete XEMWall;

  outfile1.open(output1.Data(), ios::out | ios::app);
  outfile.open(output.Data(), ios::out | ios::app);

  int cnt=0;
  while(!infile.eof()) // To get all the lines.
    {
      cnt++;
      ss.clear ();
      ss.str ("");
  
      getline(infile,Line);                                          
      ss << Line;  
      ss >> E0>>Ep>>Theta;
      Q2 = 2*E0*Ep*(1-cos(Theta/180.*3.14159265));

      xbj = Q2/2/0.938/(E0-Ep);
                                       	              	                       
      if( Q2<=25){

	cs_qe=0, cs_dis=0, cs_rad=0,cs_born=0;

	/*Set Target{*/
	if(tarid == 1) {
	  A = 2; Z = 1;}
	else if(tarid == 2) {
	  A = 3; Z = 1;}
	else if(tarid == 3) {
	  A = 3; Z = 2;}
	else if(tarid == 6) { 
	  A = 12; Z = 6;}
	else if(tarid==4||tarid==5||tarid==0) { 
	  A = 27; Z = 13;}
	//	else if(Target == "Ca40") {
	//		A = 40; Z = 20;}
	//	else if(Target == "Ca48") {
	//		A = 48; Z = 20;}
	//	else if(Target == "Dummy") { 
	//		A = 27; Z = 13;}
	//	else{
	//		cerr<<"I don't understand the Target!"<<endl;}
	//	/*}}}*/
	//
	//

	 
	//  double y = gGet_Y(E0, Ep,Theta*3.1415926/180., XEMTarget);
	//double y_wall = gGet_Y(E0, Ep,Theta*3.1415926/180., XEMWall);
	  
	
     
	//Define a event to calculate radiated cross section
	XEMC* Event = new XEMC(); 

	Event->Init(Target_Input.Data());
	Int_t err = -1;
	//	  cout<<A<<"  "<<Z<<endl;
	err = Event->Process(E0,Ep,Theta,A,Z,0.0);
	if(err>=0){
	  cs_qe   = Event->XS_QE();
	  cs_dis  = Event->XS_DIS();

	  cs_rad  = Event->XS_Rad();
	 cs_born = Event->XS_Born();
	 // cs_el = Event->XS_Elastic(); // fm2/sr
	}
	delete Event;
	//	   cout<<A<<"  "<<Z<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	////	cerr << Form("For Ep=%f, Theta=%f, xbj=%f, Q2=%f", Ep, Theta, xbj, Q2)<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	//	cerr << Form("@@@ XS_QE = %e, XS_DIS = %e, XS_Born = %e, XS_Rad = %e nb/sr/MeV, XS_Elastic= %e fm2/sr", cs_qe, cs_dis, cs_born,cs_rad, cs_el)<<endl;
	//	cerr <<"------------------------------------------------------------------------"<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	if(cnt % 100 ==0){

	  outfile.close();	//	   cout<<A<<"  "<<Z<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	cerr << Form("For Ep=%f, Theta=%f, xbj=%f, Q2=%f", Ep, Theta, xbj, Q2)<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
	cerr << Form("@@@ XS_QE = %e, XS_DIS = %e, XS_Born = %e, XS_Rad = %e nb/sr/MeV, XS_Elastic= %e fm2/sr", cs_qe, cs_dis, cs_born,cs_rad, cs_el)<<endl;
	cerr <<"------------------------------------------------------------------------"<<endl;
	//	  cerr <<"------------------------------------------------------------------------"<<endl;
;
	  outfile1.close();
	  outfile1.open(output1.Data(), ios::out | ios::app);
	  outfile.open(output.Data(), ios::out | ios::app);

	}
	
	outfile << setw(5)<<A << " "<<setw(5)<<A<< setw(10)<< xbj << "  "<< setw(10)<< Q2 << " "<< setw(5)<< Theta <<"   "<< setw(10)<<cs_qe<<"   "<< setw(10)<<cs_dis<<"   "<<setw(10)<<cs_born<< "   "<<setw(10)<< cs_rad<<'\n';

	  

	//	  if((cs_born/cs_rad)>1111) cs_rad=cs_born;
	if(cs_rad>-1e-34)
	//  outfile1 <<setw(5)<< A<<"  "<<setw(5)<<Z<<"  "<<setw(5)<<E0 << " "<<setw(10)<< Ep << "  "<<setw(5)<< Theta << "  "<<setw(10)<<cs_born*1000<< "   "<<setw(10)<< cs_born/cs_rad<< "  "<<xbj<< "  "<<Q2 <<'\n'; // nb/sr/GeV
	  outfile1 <<setw(5)<< A<<"  "<<setw(5)<<Z<<"  "<<setw(5)<<E0 << " "<<setw(10)<< Ep << "  "<<setw(5)<< Theta << "  "<<setw(10)<<cs_born*1000<< "   "<<setw(10)<< cs_born/cs_rad<< "  "<<xbj<< "  "<<Q2 <<'\n'; // nb/sr/GeV
	 
      }
 
    }
  outfile.close();
  outfile1.close();
}

/*int getargs(int argc,char** argv){{{*/
int getargs(int argc,char** argv){
  char* argptr;
  bool noStop;  
  int err=-1;

  for(int i=0;i<argc;i++){
    argptr = argv[i];
    if(*argptr=='-'){
      argptr++;
      switch (*argptr){
      case 'E':
	if(*(++argptr))
	  E0 = atof(argptr);
	err=0;
	break;
      case 'P':
	if(*(++argptr))
	  Ep = atof(argptr);
	err=0;
	break;
      case 'A':
	if(*(++argptr))
	  Theta = atof(argptr);
	err=0;
	break;
      case 'T':
	if(*(++argptr))
	  gTarget = argptr;
	err=0;
	break;
      case 'h':
	cerr <<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-E[GeV]  Beam Energy" <<endl;
	cerr <<"-P[GeV]  Scattering Momentum" <<endl;
	cerr <<"-A[Deg]  Scattering Angle"<<endl;
	cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48,...]  Target Name" <<endl;
	cerr <<endl<<"================================================"<<endl;
	noStop=false;
	err=-1;
	goto OUT;
	break;
      default:
	cerr <<"Unrecognized argument: " << argptr << endl;
	cerr <<endl<<"================================================"<<endl;
	cerr <<"Option: " <<endl;
	cerr <<"-E[GeV]  Beam Energy" <<endl;
	cerr <<"-P[GeV]  Scattering Momentum" <<endl;
	cerr <<"-A[Deg]  Scattering Angle"<<endl;
	cerr <<"-T[H2,He3,He4,C12,Ca40,Ca48,...]  Target Name" <<endl;
	cerr <<"             Zhihong Ye 10/04/2009" <<endl;
	cerr <<"================================================"<<endl;
	err=-1;
	break;
      }

    }
    noStop=true;
  OUT:if(!noStop){break;}
    continue;
  }
  return err; 
}
/**/
