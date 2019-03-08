#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <functional>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <utility>

/**
   This  code  match the  ellipsoid  in  the  simulation identified  with
   different ellipsoidal overdensities.   It start from the virial  haloes and see
   what are the haloes in the other catalogues that match.

   vesion MiEllipsoid-1.0
   created 26 - 03 - 2015 (cgiocoli@gmail.com)
   
 */

int nlist;

static const char fINT[] = "%i";

template <class T> std:: string conv (T &val, const char *fact){
  char VAL[20]; sprintf (VAL, fact, val);
  return std:: string(VAL);
}

/* with c++11 (lambda)
  template<class Vals>
  void sortingPermutation(const Vals& values, std::vector<int>& v){
  int size = values.size(); 
  v.clear(); v.reserve(size);
  for(int i=0; i < size; ++i)
  v.push_back(i);
  
  std::sort(v.begin(), v.end(), [&values](int a, int b) -> bool { 
  return values[a] < values[b];
  });
  }
*/

template<class T>
struct CmpPairs{
  CmpPairs(const std::vector<T> &v): v_(v) {}
  std::vector<T> v_;
  bool operator()(int a, int b){ return v_[a] < v_[b]; }
};

template<class T>
CmpPairs<T> CreateCmpPairs(const std::vector<T> & v) { return CmpPairs<T>(v); }

// coordinates have to be in [0,1] and valid for periodic BC
void makeLlist(std:: vector<double> x, std:: vector<double> y, std:: vector<double> z,double dl, int *llist, int ***ihead){
  int n = x.size();
  for(int i=0;i<nlist;i++){
    for(int j=0;j<nlist;j++){
      for(int k=0;k<nlist;k++){
	ihead[i][j][k] = -1;
      }
    }
  }
  for(int i=0;i<n;i++) llist[i] = -1;

  int j,k,l;
  for(int i=0;i<n;i++){
    x[i] = fmod(x[i],1.);
    y[i] = fmod(y[i],1.);
    z[i] = fmod(z[i],1.);
    
    j = int(x[i]/dl);
    k = int(y[i]/dl); 
    l = int(z[i]/dl);
    
    if(j>=0 && j<nlist && k>=0 && k<nlist && l>=0 && l<nlist){
      llist[i] = ihead[j][k][l];
      ihead[j][k][l] = i;
    }
  }
}

bool fexists(const std::string& filename) {
  std:: ifstream ifile(filename.c_str());
  return ifile;
}

int main(){
  std:: string pathvir,path2000rhoC,path1000rhoC,path500rhoC,path200rhoC,path200rhoB;
  std:: string epathvir,epath2000rhoC,epath1000rhoC,epath500rhoC,epath200rhoC,epath200rhoB;
  std:: string sbut;
  std:: ifstream filinput;
  std:: string fin = "INPUT_E";
  double boxl;
  int snap;
  bool f = fexists(fin);
  if(f){
    filinput.open(fin.c_str());
    filinput >> sbut;
    filinput >> pathvir;
    filinput >> sbut;
    filinput >> path2000rhoC;
    filinput >> sbut;
    filinput >> path1000rhoC;
    filinput >> sbut;
    filinput >> path500rhoC;
    filinput >> sbut;
    filinput >> path200rhoC;
    filinput >> sbut;
    filinput >> path200rhoB;
    filinput >> sbut;
    filinput >> epathvir;
    filinput >> sbut;
    filinput >> epath2000rhoC;
    filinput >> sbut;
    filinput >> epath1000rhoC;
    filinput >> sbut;
    filinput >> epath500rhoC;
    filinput >> sbut;
    filinput >> epath200rhoC;
    filinput >> sbut;
    filinput >> epath200rhoB;
    filinput >> sbut;
    filinput >> snap;
    filinput >> sbut;
    filinput >> boxl;
    filinput >> sbut;
    filinput >> nlist;
    filinput.close();
  }else{
    std:: cout << " INPUT_E file does not exsist I will STOP here !!! " << std:: endl;
    exit(1);
  }

  // since all the units will be in [0,1]
  double dl = 1./double(nlist);

  std:: cout << " INPUT_E file read " << std:: endl;
  std:: string filvir,fil2000rhoC,fil1000rhoC,fil500rhoC,fil200rhoC,fil200rhoB;
  std:: string efilvir,efil2000rhoC,efil1000rhoC,efil500rhoC,efil200rhoC,efil200rhoB;
  
  std:: string ssnap;

  if(snap<10) ssnap = "00"+conv(snap,fINT);
  else if(snap>=10 && snap<100) ssnap = "0"+conv(snap,fINT);
  else ssnap = conv(snap,fINT);

  // set the filenames
  filvir = pathvir + "/SO." + ssnap + ".gv";
  fil2000rhoC = path2000rhoC + "/SO." + ssnap + ".gv";
  fil1000rhoC = path1000rhoC + "/SO." + ssnap + ".gv";
  fil500rhoC = path500rhoC + "/SO." + ssnap + ".gv";
  fil200rhoC = path200rhoC + "/SO." + ssnap + ".gv";
  fil200rhoB = path200rhoB + "/SO." + ssnap + ".gv";

  // set the filenames
  efilvir = epathvir + "/EO." + ssnap + ".gv";
  efil2000rhoC = epath2000rhoC + "/EO." + ssnap + ".gv";
  efil1000rhoC = epath1000rhoC + "/EO." + ssnap + ".gv";
  efil500rhoC = epath500rhoC + "/EO." + ssnap + ".gv";
  efil200rhoC = epath200rhoC + "/EO." + ssnap + ".gv";
  efil200rhoB = epath200rhoB + "/EO." + ssnap + ".gv";
  
  f = fexists(filvir);
  if(f){
    std:: cout << " reading the following files : " << std:: endl;
    std:: cout << filvir << std:: endl;
  }
  else{
    std:: cout << filvir << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil2000rhoC);
  if(f){
    std:: cout << fil2000rhoC << std:: endl;
  }else{
    std:: cout << fil2000rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }
  
  f = fexists(fil1000rhoC);
  if(f){
    std:: cout << fil1000rhoC << std:: endl;
  }else{
    std:: cout << fil1000rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil500rhoC);
  if(f){
    std:: cout << fil500rhoC << std:: endl;
  }else{
    std:: cout << fil500rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil200rhoC);
  if(f){
    std:: cout << fil200rhoC << std:: endl;
  }else{
    std:: cout << fil200rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil200rhoB);
  if(f){
    std:: cout << fil200rhoB << std:: endl;
  }else{
    std:: cout << fil200rhoB << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(efilvir);
  if(f){
    std:: cout << " reading the following files : " << std:: endl;
    std:: cout << efilvir << std:: endl;
  }
  else{
    std:: cout << efilvir << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(efil2000rhoC);
  if(f){
    std:: cout << efil2000rhoC << std:: endl;
  }else{
    std:: cout << efil2000rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }
  
  f = fexists(efil1000rhoC);
  if(f){
    std:: cout << efil1000rhoC << std:: endl;
  }else{
    std:: cout << efil1000rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(efil500rhoC);
  if(f){
    std:: cout << efil500rhoC << std:: endl;
  }else{
    std:: cout << efil500rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(efil200rhoC);
  if(f){
    std:: cout << efil200rhoC << std:: endl;
  }else{
    std:: cout << efil200rhoC << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(efil200rhoB);
  if(f){
    std:: cout << efil200rhoB << std:: endl;
  }else{
    std:: cout << efil200rhoB << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }
 
  /////
  std:: ifstream ifilvir,ifil2000rhoC,ifil1000rhoC,ifil500rhoC,ifil200rhoC,ifil200rhoB;
  std:: ifstream iefilvir,iefil2000rhoC,iefil1000rhoC,iefil500rhoC,iefil200rhoC,iefil200rhoB;

  double msb,m,r,x,y,z,vx,vy,vz,vdisp,etot,spin,rbut,a,b,c;
  int idh,n,n2,nsb,ibut,idmb,ide;

  // open virial catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> mvir,rvir,xcmvir,ycmvir,zcmvir;
  std:: vector<int> idhvir,nvir0;  
  ifilvir.open(filvir.c_str());
  while(ifilvir >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idhvir.push_back(idh);
    mvir.push_back(m);
    nvir0.push_back(n);
    rvir.push_back(r/boxl/1.e+3);
    xcmvir.push_back(x/boxl);
    ycmvir.push_back(y/boxl);
    zcmvir.push_back(z/boxl);
  }
  ifilvir.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> mvirell,avirell,bvirell,cvirell,xcmevir,ycmevir,zcmevir,revir;
  std:: vector<int> idevir,idhevir,nvirell,nvir;  
  iefilvir.open(efilvir.c_str());
  while(iefilvir >> ide >> idh >> n >> n2 >> m >> a >> b >> c >> x >> y >> z){
    
    idevir.push_back(ide);
    idhevir.push_back(idh);
    mvirell.push_back(m);
    nvirell.push_back(n);
    // point at the virial radius of the corresponding SO halo!
    revir.push_back(rvir[idh-1]);
    nvir.push_back(nvir0[idh-1]);

    xcmevir.push_back(x/boxl);
    ycmevir.push_back(y/boxl);
    zcmevir.push_back(z/boxl);
    
    avirell.push_back(a);
    bvirell.push_back(b);
    cvirell.push_back(c);

  }
  iefilvir.close();  
  
  std:: cout << " 1 " << std:: endl;
  
  // open 2000rhoC catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m2000rhoC,r2000rhoC,xcm2000rhoC,ycm2000rhoC,zcm2000rhoC;
  std:: vector<int> idh2000rhoC;  
  ifil2000rhoC.open(fil2000rhoC.c_str());
  while(ifil2000rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh2000rhoC.push_back(idh);
    m2000rhoC.push_back(m);
    r2000rhoC.push_back(r/boxl/1.e+3);
    xcm2000rhoC.push_back(x/boxl);
    ycm2000rhoC.push_back(y/boxl);
    zcm2000rhoC.push_back(z/boxl);
  }
  ifil2000rhoC.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m2000rhoCell,a2000rhoCell,b2000rhoCell,c2000rhoCell,xcme2000rhoC,ycme2000rhoC,zcme2000rhoC,re2000rhoC;
  std:: vector<int> ide2000rhoC,idhe2000rhoC;  
  iefil2000rhoC.open(efil2000rhoC.c_str());
  while(iefil2000rhoC >> ide >> idh >> n >> n >> m >> a >> b >> c >> x >> y >> z){
    
    ide2000rhoC.push_back(ide);
    idhe2000rhoC.push_back(idh);
    m2000rhoCell.push_back(m);
    // point at the 2000rhoC radius of the corresponding SO halo!
    re2000rhoC.push_back(r2000rhoC[idh-1]);
    xcme2000rhoC.push_back(x/boxl);
    ycme2000rhoC.push_back(y/boxl);
    zcme2000rhoC.push_back(z/boxl);
    
    a2000rhoCell.push_back(a);
    b2000rhoCell.push_back(b);
    c2000rhoCell.push_back(c);

  }
  iefil2000rhoC.close();  

  int ***ihead2000rhoC;
  int *llist2000rhoC;
  // Allocate memory
  llist2000rhoC = new int[xcme2000rhoC.size()];
  ihead2000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead2000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead2000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcme2000rhoC,ycme2000rhoC,zcme2000rhoC,dl,llist2000rhoC,ihead2000rhoC);

  std:: cout << " 2 " << std:: endl;

  // open 1000rhoC catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m1000rhoC,r1000rhoC,xcm1000rhoC,ycm1000rhoC,zcm1000rhoC;
  std:: vector<int> idh1000rhoC;  
  ifil1000rhoC.open(fil1000rhoC.c_str());
  while(ifil1000rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh1000rhoC.push_back(idh);
    m1000rhoC.push_back(m);
    r1000rhoC.push_back(r/boxl/1.e+3);
    xcm1000rhoC.push_back(x/boxl);
    ycm1000rhoC.push_back(y/boxl);
    zcm1000rhoC.push_back(z/boxl);
  }
  ifil1000rhoC.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m1000rhoCell,a1000rhoCell,b1000rhoCell,c1000rhoCell,xcme1000rhoC,ycme1000rhoC,zcme1000rhoC,re1000rhoC;
  std:: vector<int> ide1000rhoC,idhe1000rhoC;  
  iefil1000rhoC.open(efil1000rhoC.c_str());
  while(iefil1000rhoC >> ide >> idh >> n >> n >> m >> a >> b >> c >> x >> y >> z){
    
    ide1000rhoC.push_back(ide);
    idhe1000rhoC.push_back(idh);
    m1000rhoCell.push_back(m);
    // point at the 1000rhoC radius of the corresponding SO halo!
    re1000rhoC.push_back(r1000rhoC[idh-1]);
    xcme1000rhoC.push_back(x/boxl);
    ycme1000rhoC.push_back(y/boxl);
    zcme1000rhoC.push_back(z/boxl);
    
    a1000rhoCell.push_back(a);
    b1000rhoCell.push_back(b);
    c1000rhoCell.push_back(c);

  }
  iefil1000rhoC.close();  

  int ***ihead1000rhoC;
  int *llist1000rhoC;
  // Allocate memory
  llist1000rhoC = new int[xcme1000rhoC.size()];
  ihead1000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead1000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead1000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcme1000rhoC,ycme1000rhoC,zcme1000rhoC,dl,llist1000rhoC,ihead1000rhoC);
  std:: cout << " 3 " << std:: endl;

  // open 500rhoC catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m500rhoC,r500rhoC,xcm500rhoC,ycm500rhoC,zcm500rhoC;
  std:: vector<int> idh500rhoC;  
  ifil500rhoC.open(fil500rhoC.c_str());
  while(ifil500rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh500rhoC.push_back(idh);
    m500rhoC.push_back(m);
    r500rhoC.push_back(r/boxl/1.e+3);
    xcm500rhoC.push_back(x/boxl);
    ycm500rhoC.push_back(y/boxl);
    zcm500rhoC.push_back(z/boxl);
  }
  ifil500rhoC.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m500rhoCell,a500rhoCell,b500rhoCell,c500rhoCell,xcme500rhoC,ycme500rhoC,zcme500rhoC,re500rhoC;
  std:: vector<int> ide500rhoC,idhe500rhoC;  
  iefil500rhoC.open(efil500rhoC.c_str());
  while(iefil500rhoC >> ide >> idh >> n >> n >> m >> a >> b >> c >> x >> y >> z){
    
    ide500rhoC.push_back(ide);
    idhe500rhoC.push_back(idh);
    m500rhoCell.push_back(m);
    // point at the 500rhoC radius of the corresponding SO halo!
    re500rhoC.push_back(r500rhoC[idh-1]);
    xcme500rhoC.push_back(x/boxl);
    ycme500rhoC.push_back(y/boxl);
    zcme500rhoC.push_back(z/boxl);
    
    a500rhoCell.push_back(a);
    b500rhoCell.push_back(b);
    c500rhoCell.push_back(c);

  }
  iefil500rhoC.close();  

  int ***ihead500rhoC;
  int *llist500rhoC;
  // Allocate memory
  llist500rhoC = new int[xcme500rhoC.size()];
  ihead500rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead500rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead500rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcme500rhoC,ycme500rhoC,zcme500rhoC,dl,llist500rhoC,ihead500rhoC);
  std:: cout << " 4 " << std:: endl;

  // open 200rhoC catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m200rhoC,r200rhoC,xcm200rhoC,ycm200rhoC,zcm200rhoC;
  std:: vector<int> idh200rhoC;  
  ifil200rhoC.open(fil200rhoC.c_str());
  while(ifil200rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh200rhoC.push_back(idh);
    m200rhoC.push_back(m);
    r200rhoC.push_back(r/boxl/1.e+3);
    xcm200rhoC.push_back(x/boxl);
    ycm200rhoC.push_back(y/boxl);
    zcm200rhoC.push_back(z/boxl);
  }
  ifil200rhoC.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m200rhoCell,a200rhoCell,b200rhoCell,c200rhoCell,xcme200rhoC,ycme200rhoC,zcme200rhoC,re200rhoC;
  std:: vector<int> ide200rhoC,idhe200rhoC;  
  iefil200rhoC.open(efil200rhoC.c_str());
  while(iefil200rhoC >> ide >> idh >> n >> n >> m >> a >> b >> c >> x >> y >> z){
    
    ide200rhoC.push_back(ide);
    idhe200rhoC.push_back(idh);
    m200rhoCell.push_back(m);
    // point at the 200rhoC radius of the corresponding SO halo!
    re200rhoC.push_back(r200rhoC[idh-1]);
    xcme200rhoC.push_back(x/boxl);
    ycme200rhoC.push_back(y/boxl);
    zcme200rhoC.push_back(z/boxl);
    
    a200rhoCell.push_back(a);
    b200rhoCell.push_back(b);
    c200rhoCell.push_back(c);

  }
  iefil200rhoC.close();  

  int ***ihead200rhoC;
  int *llist200rhoC;
  // Allocate memory
  llist200rhoC = new int[xcme200rhoC.size()];
  ihead200rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead200rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead200rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcme200rhoC,ycme200rhoC,zcme200rhoC,dl,llist200rhoC,ihead200rhoC);
  std:: cout << " 5 " << std:: endl;

  // open 200rhob catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m200rhoB,r200rhoB,xcm200rhoB,ycm200rhoB,zcm200rhoB;
  std:: vector<int> idh200rhoB;  
  ifil200rhoB.open(fil200rhoB.c_str());
  while(ifil200rhoB >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh200rhoB.push_back(idh);
    m200rhoB.push_back(m);
    r200rhoB.push_back(r/boxl/1.e+3);
    xcm200rhoB.push_back(x/boxl);
    ycm200rhoB.push_back(y/boxl);
    zcm200rhoB.push_back(z/boxl);
  }
  ifil200rhoB.close();  

  // open ellipsoid catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m200rhoBell,a200rhoBell,b200rhoBell,c200rhoBell,xcme200rhoB,ycme200rhoB,zcme200rhoB,re200rhoB;
  std:: vector<int> ide200rhoB,idhe200rhoB;  
  iefil200rhoB.open(efil200rhoB.c_str());
  while(iefil200rhoB >> ide >> idh >> n >> n >> m >> a >> b >> c >> x >> y >> z){
    
    ide200rhoB.push_back(ide);
    idhe200rhoB.push_back(idh);
    m200rhoBell.push_back(m);
    // point at the 200rhoB radius of the corresponding SO halo!
    re200rhoB.push_back(r200rhoB[idh-1]);
    xcme200rhoB.push_back(x/boxl);
    ycme200rhoB.push_back(y/boxl);
    zcme200rhoB.push_back(z/boxl);
    
    a200rhoBell.push_back(a);
    b200rhoBell.push_back(b);
    c200rhoBell.push_back(c);

  }
  iefil200rhoB.close();  

  int ***ihead200rhoB;
  int *llist200rhoB;
  // Allocate memory
  llist200rhoB = new int[xcme200rhoB.size()];
  ihead200rhoB = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead200rhoB[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead200rhoB[i][j] = new int[nlist];
  }
  makeLlist(xcme200rhoB,ycme200rhoB,zcme200rhoB,dl,llist200rhoB,ihead200rhoB);
  std:: cout << " 6 " << std:: endl;

  // number of haloes in the EO catalogues
  int nevir = idevir.size();
  int ne2000rhoC = ide2000rhoC.size();
  int ne1000rhoC = ide1000rhoC.size();
  int ne500rhoC = ide500rhoC.size();
  int ne200rhoC = ide200rhoC.size();
  int ne200rhoB = ide200rhoB.size();


  std:: vector<int> idMie2000rhoC(nevir);
  std:: vector<int> dMie2000rhoC(nevir);
  std:: vector<int> idMie1000rhoC(nevir);
  std:: vector<int> dMie1000rhoC(nevir);
  std:: vector<int> idMie500rhoC(nevir);
  std:: vector<int> dMie500rhoC(nevir);
  std:: vector<int> idMie200rhoC(nevir);
  std:: vector<int> dMie200rhoC(nevir);
  std:: vector<int> idMie200rhoB(nevir);
  std:: vector<int> dMie200rhoB(nevir);

  std:: string filout;
  filout = "MatchingEllipsoids." + ssnap + ".gv";
  std:: ofstream ofilout;
  ofilout.open(filout.c_str());

  std:: string filout2;
  filout2 = "MatchingEllipsoids.Prop." + ssnap + ".gv";
  std:: ofstream ofilout2;
  ofilout2.open(filout2.c_str());

  for(int i=0;i<nevir;i++){
    // for each halo compute the distance of the ellipsoid in the other catalogue
    std:: vector<double> d2000rhoC;
    std:: vector<int> epos2000rhoC,eid2000rhoC;
    std:: vector<double> d1000rhoC;
    std:: vector<int> epos1000rhoC,eid1000rhoC;
    std:: vector<double> d500rhoC;
    std:: vector<int> epos500rhoC,eid500rhoC;
    std:: vector<double> d200rhoC;
    std:: vector<int> epos200rhoC,eid200rhoC;
    std:: vector<double> d200rhoB;
    std:: vector<int> epos200rhoB,eid200rhoB;
    int sele2000rhoC=0;
    int sele1000rhoC=0;
    int sele500rhoC=0;
    int sele200rhoC=0;
    int sele200rhoB=0;
    double xi = xcmevir[i];
    double yi = ycmevir[i];
    double zi = zcmevir[i];
    // cell containing the ellipsoid
    int jx = int(xi/dl);
    int jy = int(yi/dl);
    int jz = int(zi/dl);
    // number of cells contaning x times the virial radius of the ellipsoid
    int jr = int(1.5*revir[i]/dl)+1;

    int jx1 = jx-jr;
    int jy1 = jy-jr;
    int jz1 = jz-jr;
    int jx2 = jx+jr;
    int jy2 = jy+jr;
    int jz2 = jz+jr;

    std:: cout << " searching for neighbours of ellipsoid vir " << i << " in other catalogues " << std:: endl;
    for(int l=jz1;l<=jz2;l++){
      int li = l;
      // wrapping for periodic BC
      if(li<0) li = nlist + li;
      if(li>=nlist) li = li - nlist;
      for(int k=jy1;k<=jy2;k++){
	int ki = k;
	// wrapping for periodic BC
	if(ki<0) ki = nlist + ki;
	if(ki>=nlist) ki = ki - nlist;
	for(int j=jx1;j<=jx2;j++){
	  int ji = j;
	  // wrapping for periodic BC
	  if(ji<0) ji = nlist + ji;
	  if(ji>=nlist) ji = ji - nlist;

	  // head ellipsoid for this cell in 2000rhoC
	  int mm = ihead2000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcme2000rhoC[mm] - xi;
	    double dy = ycme2000rhoC[mm] - yi;
	    double dz = zcme2000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*revir[i]){
	      d2000rhoC.push_back(dist);
	      epos2000rhoC.push_back(sele2000rhoC);
	      eid2000rhoC.push_back(mm);
	      sele2000rhoC++;	      
	    }
	    //...take next ellipsoid and repeat
	    mm = llist2000rhoC[mm];
	  }

	  // head ellipsoid for this cell in 1000rhoC
	  mm = ihead1000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcme1000rhoC[mm] - xi;
	    double dy = ycme1000rhoC[mm] - yi;
	    double dz = zcme1000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*revir[i]){
	      d1000rhoC.push_back(dist);
	      epos1000rhoC.push_back(sele1000rhoC);
	      eid1000rhoC.push_back(mm);
	      sele1000rhoC++;	      
	    }
	    //...take next ellipsoid and repeat
	    mm = llist1000rhoC[mm];
	  }

	  // head ellipsoid for this cell in 500rhoC
	  mm = ihead500rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcme500rhoC[mm] - xi;
	    double dy = ycme500rhoC[mm] - yi;
	    double dz = zcme500rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*revir[i]){
	      d500rhoC.push_back(dist);
	      epos500rhoC.push_back(sele500rhoC);
	      eid500rhoC.push_back(mm);
	      sele500rhoC++;	      
	    }
	    //...take next ellipsoid and repeat
	    mm = llist500rhoC[mm];
	  }

	  // head ellipsoid for this cell in 200rhoC
	  mm = ihead200rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcme200rhoC[mm] - xi;
	    double dy = ycme200rhoC[mm] - yi;
	    double dz = zcme200rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*revir[i]){
	      d200rhoC.push_back(dist);
	      epos200rhoC.push_back(sele200rhoC);
	      eid200rhoC.push_back(mm);
	      sele200rhoC++;	      
	    }
	    //...take next ellipsoid and repeat
	    mm = llist200rhoC[mm];
	  }

	  // head ellipsoid for this cell in 200rhoB
	  mm = ihead200rhoB[ji][ki][li];
	  while(mm>-1){
	    double dx = xcme200rhoB[mm] - xi;
	    double dy = ycme200rhoB[mm] - yi;
	    double dz = zcme200rhoB[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*revir[i]){
	      d200rhoB.push_back(dist);
	      epos200rhoB.push_back(sele200rhoB);
	      eid200rhoB.push_back(mm);
	      sele200rhoB++;	      
	    }
	    //...take next ellipsoid and repeat
	    mm = llist200rhoB[mm];
	  }
	}
      }
    }
    int iw2000rhoC,iwh2000rhoC;
    double mw2000rhoC,rw2000rhoC,dr2000rhoC,aw2000rhoC,bw2000rhoC,cw2000rhoC;
    if(sele2000rhoC>0){
      std::sort(epos2000rhoC.begin(),epos2000rhoC.end(), CreateCmpPairs(d2000rhoC));
      idMie2000rhoC[i] = eid2000rhoC[epos2000rhoC[0]];
      dMie2000rhoC[i] = d2000rhoC[epos2000rhoC[0]];
      iw2000rhoC=idMie2000rhoC[i]+1;
      mw2000rhoC=m2000rhoCell[idMie2000rhoC[i]];
      iwh2000rhoC=idhe2000rhoC[idMie2000rhoC[i]];
      rw2000rhoC=re2000rhoC[idMie2000rhoC[i]]*boxl;
      dr2000rhoC=d2000rhoC[epos2000rhoC[0]]*boxl;
      aw2000rhoC=a2000rhoCell[idMie2000rhoC[i]];
      bw2000rhoC=b2000rhoCell[idMie2000rhoC[i]];
      cw2000rhoC=c2000rhoCell[idMie2000rhoC[i]];
    }
    else{
      iw2000rhoC=0;
      mw2000rhoC=0.;
      iwh2000rhoC=0;
      rw2000rhoC=0.;
      dr2000rhoC=0.;
      aw2000rhoC=0;
      bw2000rhoC=0;
      cw2000rhoC=0;
    }
    int iw1000rhoC,iwh1000rhoC;
    double mw1000rhoC,rw1000rhoC,dr1000rhoC,aw1000rhoC,bw1000rhoC,cw1000rhoC;
    if(sele1000rhoC>0){
      std::sort(epos1000rhoC.begin(),epos1000rhoC.end(), CreateCmpPairs(d1000rhoC));
      idMie1000rhoC[i] = eid1000rhoC[epos1000rhoC[0]];
      dMie1000rhoC[i] = d1000rhoC[epos1000rhoC[0]];
      iw1000rhoC=idMie1000rhoC[i]+1;
      mw1000rhoC=m1000rhoCell[idMie1000rhoC[i]];
      iwh1000rhoC=idhe1000rhoC[idMie1000rhoC[i]];
      rw1000rhoC=re1000rhoC[idMie1000rhoC[i]]*boxl;
      dr1000rhoC=d1000rhoC[epos1000rhoC[0]]*boxl;
      aw1000rhoC=a1000rhoCell[idMie1000rhoC[i]];
      bw1000rhoC=b1000rhoCell[idMie1000rhoC[i]];
      cw1000rhoC=c1000rhoCell[idMie1000rhoC[i]];
    }     
    else{
      iw1000rhoC=0;
      mw1000rhoC=0.;
      iwh1000rhoC=0;
      rw1000rhoC=0.;
      dr1000rhoC=0.;
      aw1000rhoC=0;
      bw1000rhoC=0;
      cw1000rhoC=0;
    }
    int iw500rhoC,iwh500rhoC;
    double mw500rhoC,rw500rhoC,dr500rhoC,aw500rhoC,bw500rhoC,cw500rhoC;
    if(sele500rhoC>0){
      std::sort(epos500rhoC.begin(),epos500rhoC.end(), CreateCmpPairs(d500rhoC));
      idMie500rhoC[i] = eid500rhoC[epos500rhoC[0]];
      dMie500rhoC[i] = d500rhoC[epos500rhoC[0]];
      iw500rhoC=idMie500rhoC[i]+1;
      mw500rhoC=m500rhoCell[idMie500rhoC[i]];
      iwh500rhoC=idhe500rhoC[idMie500rhoC[i]];
      rw500rhoC=re500rhoC[idMie500rhoC[i]]*boxl;
      dr500rhoC=d500rhoC[epos500rhoC[0]]*boxl;
      aw500rhoC=a500rhoCell[idMie500rhoC[i]];
      bw500rhoC=b500rhoCell[idMie500rhoC[i]];
      cw500rhoC=c500rhoCell[idMie500rhoC[i]];
    }    
    else{
      iw500rhoC=0;
      mw500rhoC=0.;
      iwh500rhoC=0;
      rw500rhoC=0.;
      dr500rhoC=0.;
      aw500rhoC=0;
      bw500rhoC=0;
      cw500rhoC=0;
    } 
    int iw200rhoC,iwh200rhoC;
    double mw200rhoC,rw200rhoC,dr200rhoC,aw200rhoC,bw200rhoC,cw200rhoC;
    if(sele200rhoC>0){
      std::sort(epos200rhoC.begin(),epos200rhoC.end(), CreateCmpPairs(d200rhoC));
      idMie200rhoC[i] = eid200rhoC[epos200rhoC[0]];
      dMie200rhoC[i] = d200rhoC[epos200rhoC[0]];
      iw200rhoC=idMie200rhoC[i]+1;
      mw200rhoC=m200rhoCell[idMie200rhoC[i]];
      iwh200rhoC=idhe200rhoC[idMie200rhoC[i]];
      rw200rhoC=re200rhoC[idMie200rhoC[i]]*boxl;
      dr200rhoC=d200rhoC[epos200rhoC[0]]*boxl;
      aw200rhoC=a200rhoCell[idMie200rhoC[i]];
      bw200rhoC=b200rhoCell[idMie200rhoC[i]];
      cw200rhoC=c200rhoCell[idMie200rhoC[i]];
    }     
    else{
      iw200rhoC=0;
      mw200rhoC=0.;
      iwh200rhoC=0;
      rw200rhoC=0.;
      dr200rhoC=0.;
      aw200rhoC=0;
      bw200rhoC=0;
      cw200rhoC=0;
    } 
    int iw200rhoB,iwh200rhoB;
    double mw200rhoB,rw200rhoB,dr200rhoB,aw200rhoB,bw200rhoB,cw200rhoB;
    if(sele200rhoB>0){
      std::sort(epos200rhoB.begin(),epos200rhoB.end(), CreateCmpPairs(d200rhoB));
      idMie200rhoB[i] = eid200rhoB[epos200rhoB[0]];
      dMie200rhoB[i] = d200rhoB[epos200rhoB[0]];
      iw200rhoB=idMie200rhoB[i]+1;
      mw200rhoB=m200rhoBell[idMie200rhoB[i]];
      iwh200rhoB=idhe200rhoB[idMie200rhoB[i]];
      rw200rhoB=re200rhoB[idMie200rhoB[i]]*boxl;
      dr200rhoB=d200rhoB[epos200rhoB[0]]*boxl;
      aw200rhoB=a200rhoBell[idMie200rhoB[i]];
      bw200rhoB=b200rhoBell[idMie200rhoB[i]];
      cw200rhoB=c200rhoBell[idMie200rhoB[i]];
    }    
    else{
      iw200rhoB=0;
      mw200rhoB=0.;
      iwh200rhoB=0;
      rw200rhoB=0.;
      dr200rhoB=0.;
      aw200rhoB=0;
      bw200rhoB=0;
      cw200rhoB=0;
    }  
    //ofilout << std:: setw(7) << i+1 << "   " << std:: scientific << mvirell[i] << "  " << std:: fixed << std:: setprecision(6) << revir[i]*boxl << "   " 
    ofilout << std:: setw(7) << i+1 << "   " << std:: setw(9) << nvir[i] << "  " << std:: setw(9) << nvirell[i] << "   " 
	    << std:: fixed << std:: setprecision(6) << revir[i]*boxl << "   " 
	    << std:: setw(7) << iw2000rhoC << "   " << std:: scientific  << mw2000rhoC << "   " << std:: fixed << std:: setprecision(6)  << rw2000rhoC << "  " << dr2000rhoC << "   " 
	    << std:: setw(7) << iw1000rhoC << "   " << std:: scientific  << mw1000rhoC << "   " << std:: fixed << std:: setprecision(6) << rw1000rhoC << "  " << dr1000rhoC << "   " 
	    << std:: setw(7) << iw500rhoC << "   " << std:: scientific << mw500rhoC << "   " << std:: fixed << std:: setprecision(6) << rw500rhoC << "  " << dr500rhoC << "   " 
	    << std:: setw(7) << iw200rhoC << "   " << std:: scientific << mw200rhoC << "   " << std:: fixed << std:: setprecision(6) << rw200rhoC << "  " << dr200rhoC << "   " 
	    << std:: setw(7) << iw200rhoB << "   " << std:: scientific  << mw200rhoB << "   " << std:: fixed << std:: setprecision(6) << rw200rhoB << "  " << dr200rhoB << std:: endl;
    
    ofilout2 << std:: setw(7) << idhevir[i] << "  " << std:: fixed << std:: setprecision(6) << avirell[i] << "  " << std:: fixed << std:: setprecision(6) << bvirell[i] << "   " << std:: fixed << std:: setprecision(6) << cvirell[i] << "  " 
	     << std:: setw(7) << idhe2000rhoC[i] << "  " << std:: fixed << std:: setprecision(6) << aw2000rhoC << "  " << std:: fixed << std:: setprecision(6) << bw2000rhoC << "   " << std:: fixed << std:: setprecision(6) << cw2000rhoC << "  " 
	     << std:: setw(7) << idhe1000rhoC[i] << "  " << std:: fixed << std:: setprecision(6) << aw1000rhoC << "  " << std:: fixed << std:: setprecision(6) << bw1000rhoC << "   " << std:: fixed << std:: setprecision(6) << cw1000rhoC << "  " 
	     << std:: setw(7) << idhe500rhoC[i] << "  " << std:: fixed << std:: setprecision(6) << aw500rhoC << "  " << std:: fixed << std:: setprecision(6) << bw500rhoC << "   " << std:: fixed << std:: setprecision(6) << cw500rhoC << "  " 
	     << std:: setw(7) << idhe200rhoC[i] << "  " << std:: fixed << std:: setprecision(6) << aw200rhoC << "  " << std:: fixed << std:: setprecision(6) << bw200rhoC << "   " << std:: fixed << std:: setprecision(6) << cw200rhoC << "  " 
	     << std:: setw(7) << idhe200rhoB[i] << "  " << std:: fixed << std:: setprecision(6) << aw200rhoB << "  " << std:: fixed << std:: setprecision(6) << bw200rhoB << "   " << std:: fixed << std:: setprecision(6) << cw200rhoB << std:: endl;
      

  }
  ofilout.close();
  ofilout2.close();
}
