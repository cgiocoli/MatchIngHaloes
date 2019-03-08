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
#include <sys/stat.h>
#include <unistd.h>

/**
   This  code  match the  haloes  in  the  simulation identified  with
   different overdensities.   It start from the virial  haloes and see
   what are the haloes in the other catalogues that match.

   vesion MiHaloes-1.5
   created  24 - 09 - 2014 (cgiocoli@gmail.com)
   modified 26 - 06 - 2017 (cgiocoli@gmail.com) to read also 5000C and 50000C cats, for C. Nipoti
   
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
  return ifile.good();
}

int main(){
  std:: string pathvir,path2000rhoC,path1000rhoC,path500rhoC,path200rhoC,path200rhoB;
  std:: string path5000rhoC, path50000rhoC;
  std:: string sbut;
  std:: ifstream filinput;
  std:: string fin = "INPUTmih.ini";
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
    filinput >> path5000rhoC;
    filinput >> sbut;
    filinput >> path50000rhoC;    
    filinput >> sbut;
    filinput >> snap;
    filinput >> sbut;
    filinput >> boxl;
    filinput >> sbut;
    filinput >> nlist;
    filinput.close();
  }else{
    std:: cout << " INPUT file does not exsist I will STOP here !!! " << std:: endl;
    exit(1);
  }

  // since all the units will be in [0,1]
  double dl = 1./double(nlist);

  std:: cout << " INPUT file read " << std:: endl;
  std:: string filvir,fil2000rhoC,fil1000rhoC,fil500rhoC,fil200rhoC,fil200rhoB;
  std:: string fil5000rhoC, fil50000rhoC;
  std:: string ssnap;

  if(snap<10) ssnap = "00"+conv(snap,fINT);
  else if(snap>=10 && snap<100) ssnap = "0"+conv(snap,fINT);
  else ssnap = conv(snap,fINT);

  // set the filenames
  filvir = pathvir + "/SO." + ssnap + ".gv";
  fil50000rhoC = path50000rhoC + "/SO." + ssnap + ".gv";  
  fil5000rhoC = path5000rhoC + "/SO." + ssnap + ".gv";  
  fil2000rhoC = path2000rhoC + "/SO." + ssnap + ".gv";
  fil1000rhoC = path1000rhoC + "/SO." + ssnap + ".gv";
  fil500rhoC = path500rhoC + "/SO." + ssnap + ".gv";
  fil200rhoC = path200rhoC + "/SO." + ssnap + ".gv";
  fil200rhoB = path200rhoB + "/SO." + ssnap + ".gv";
  
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
  
  std:: ifstream ifilvir,ifil2000rhoC,ifil1000rhoC,ifil500rhoC,ifil200rhoC,ifil200rhoB;
  std:: ifstream ifil5000rhoC, ifil50000rhoC;
  
  double msb,m,r,x,y,z,vx,vy,vz,vdisp,etot,spin,rbut;
  int idh,n,nsb,ibut,idmb;

  // open virial catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> mvir,rvir,xcmvir,ycmvir,zcmvir;
  std:: vector<int> idhvir;  
  ifilvir.open(filvir.c_str());
  while(ifilvir >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idhvir.push_back(idh);
    mvir.push_back(m);
    rvir.push_back(r/boxl/1.e+3);
    xcmvir.push_back(x/boxl);
    ycmvir.push_back(y/boxl);
    zcmvir.push_back(z/boxl);
  }
  ifilvir.close();

  std:: cout << " 50000rhoC " << std:: endl;
  
  // open 50000rhoc catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m50000rhoC,r50000rhoC,xcm50000rhoC,ycm50000rhoC,zcm50000rhoC;
  std:: vector<int> idh50000rhoC;  
  ifil50000rhoC.open(fil50000rhoC.c_str());
  while(ifil50000rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh50000rhoC.push_back(idh);
    m50000rhoC.push_back(m);
    r50000rhoC.push_back(r/boxl/1.e+3);
    xcm50000rhoC.push_back(x/boxl);
    ycm50000rhoC.push_back(y/boxl);
    zcm50000rhoC.push_back(z/boxl);
  }
  ifil50000rhoC.close();  

  int ***ihead50000rhoC;
  int *llist50000rhoC;
  // Allocate memory
  llist50000rhoC = new int[xcm50000rhoC.size()];
  ihead50000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead50000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead50000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm50000rhoC,ycm50000rhoC,zcm50000rhoC,dl,llist50000rhoC,ihead50000rhoC);  

  std:: cout << " 5000rhoC " << std:: endl;
  
  // open 5000rhoc catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m5000rhoC,r5000rhoC,xcm5000rhoC,ycm5000rhoC,zcm5000rhoC;
  std:: vector<int> idh5000rhoC;  
  ifil5000rhoC.open(fil5000rhoC.c_str());
  while(ifil5000rhoC >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh5000rhoC.push_back(idh);
    m5000rhoC.push_back(m);
    r5000rhoC.push_back(r/boxl/1.e+3);
    xcm5000rhoC.push_back(x/boxl);
    ycm5000rhoC.push_back(y/boxl);
    zcm5000rhoC.push_back(z/boxl);
  }
  ifil5000rhoC.close();  

  int ***ihead5000rhoC;
  int *llist5000rhoC;
  // Allocate memory
  llist5000rhoC = new int[xcm5000rhoC.size()];
  ihead5000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead5000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead5000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm5000rhoC,ycm5000rhoC,zcm5000rhoC,dl,llist5000rhoC,ihead5000rhoC);  
  
  std:: cout << " 2000rhoC " << std:: endl;
  
  // open 2000rhoc catalogue and read ID, Mass, Radius and CM pos
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

  int ***ihead2000rhoC;
  int *llist2000rhoC;
  // Allocate memory
  llist2000rhoC = new int[xcm2000rhoC.size()];
  ihead2000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead2000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead2000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm2000rhoC,ycm2000rhoC,zcm2000rhoC,dl,llist2000rhoC,ihead2000rhoC);

  std:: cout << " 1000rhoC " << std:: endl;

  // open 1000rhoc catalogue and read ID, Mass, Radius and CM pos
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

  int ***ihead1000rhoC;
  int *llist1000rhoC;
  // Allocate memory
  llist1000rhoC = new int[xcm1000rhoC.size()];
  ihead1000rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead1000rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead1000rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm1000rhoC,ycm1000rhoC,zcm1000rhoC,dl,llist1000rhoC,ihead1000rhoC);
  
  std:: cout << " 500rhoC " << std:: endl;

  // open 500rhoc catalogue and read ID, Mass, Radius and CM pos
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

  int ***ihead500rhoC;
  int *llist500rhoC;
  // Allocate memory
  llist500rhoC = new int[xcm500rhoC.size()];
  ihead500rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead500rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead500rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm500rhoC,ycm500rhoC,zcm500rhoC,dl,llist500rhoC,ihead500rhoC);
  
  std:: cout << " 200rhoC " << std:: endl;

  // open 200rhoc catalogue and read ID, Mass, Radius and CM pos
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

  int ***ihead200rhoC;
  int *llist200rhoC;
  // Allocate memory
  llist200rhoC = new int[xcm200rhoC.size()];
  ihead200rhoC = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead200rhoC[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead200rhoC[i][j] = new int[nlist];
  }
  makeLlist(xcm200rhoC,ycm200rhoC,zcm200rhoC,dl,llist200rhoC,ihead200rhoC);
  
  std:: cout << " 200rhoB " << std:: endl;

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

  int ***ihead200rhoB;
  int *llist200rhoB;
  // Allocate memory
  llist200rhoB = new int[xcm200rhoB.size()];
  ihead200rhoB = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead200rhoB[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead200rhoB[i][j] = new int[nlist];
  }
  makeLlist(xcm200rhoB,ycm200rhoB,zcm200rhoB,dl,llist200rhoB,ihead200rhoB);
  
  std:: cout << " all linked lists done " << std:: endl;

  // number of haloes in the catalogues
  int nhvir = idhvir.size();
  int nh50000rhoC = idh50000rhoC.size();
  int nh5000rhoC = idh5000rhoC.size();
  int nh2000rhoC = idh2000rhoC.size();  
  int nh1000rhoC = idh1000rhoC.size();
  int nh500rhoC = idh500rhoC.size();
  int nh200rhoC = idh200rhoC.size();
  int nh200rhoB = idh200rhoB.size();

  std:: vector<int> idMih50000rhoC(nhvir);
  std:: vector<int> dMih50000rhoC(nhvir);
  std:: vector<int> idMih5000rhoC(nhvir);
  std:: vector<int> dMih5000rhoC(nhvir);
  std:: vector<int> idMih2000rhoC(nhvir);
  std:: vector<int> dMih2000rhoC(nhvir);
  std:: vector<int> idMih1000rhoC(nhvir);
  std:: vector<int> dMih1000rhoC(nhvir);
  std:: vector<int> idMih500rhoC(nhvir);
  std:: vector<int> dMih500rhoC(nhvir);
  std:: vector<int> idMih200rhoC(nhvir);
  std:: vector<int> dMih200rhoC(nhvir);
  std:: vector<int> idMih200rhoB(nhvir);
  std:: vector<int> dMih200rhoB(nhvir);

  std:: string filout;
  filout = "MatchingHaloes-1.5_." + ssnap + ".gv";
  std:: ofstream ofilout;
  ofilout.open(filout.c_str());

  for(int i=0;i<nhvir;i++){
    // for each halo compute the distance of the haloes in the other catalogue
    std:: vector<double> d50000rhoC;
    std:: vector<int> hpos50000rhoC,hid50000rhoC;
    std:: vector<double> d5000rhoC;
    std:: vector<int> hpos5000rhoC,hid5000rhoC;
    std:: vector<double> d2000rhoC;
    std:: vector<int> hpos2000rhoC,hid2000rhoC;
    std:: vector<double> d1000rhoC;
    std:: vector<int> hpos1000rhoC,hid1000rhoC;
    std:: vector<double> d500rhoC;
    std:: vector<int> hpos500rhoC,hid500rhoC;
    std:: vector<double> d200rhoC;
    std:: vector<int> hpos200rhoC,hid200rhoC;
    std:: vector<double> d200rhoB;
    std:: vector<int> hpos200rhoB,hid200rhoB;
    int selh50000rhoC=0;
    int selh5000rhoC=0;
    int selh2000rhoC=0;    
    int selh1000rhoC=0;
    int selh500rhoC=0;
    int selh200rhoC=0;
    int selh200rhoB=0;
    double xi = xcmvir[i];
    double yi = ycmvir[i];
    double zi = zcmvir[i];
    // cell containing the halo
    int jx = int(xi/dl);
    int jy = int(yi/dl);
    int jz = int(zi/dl);
    // number of cells contaning x times the virial radius of the halo
    int jr = int(1.5*rvir[i]/dl)+1;

    int jx1 = jx-jr;
    int jy1 = jy-jr;
    int jz1 = jz-jr;
    int jx2 = jx+jr;
    int jy2 = jy+jr;
    int jz2 = jz+jr;

    std:: cout << " searching for neighbours of halo vir " << i << " in other catalogues " << std:: endl;
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

	  // head halo for this cell in 50000rhoC
	  int mm = ihead50000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm50000rhoC[mm] - xi;
	    double dy = ycm50000rhoC[mm] - yi;
	    double dz = zcm50000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d50000rhoC.push_back(dist);
	      hpos50000rhoC.push_back(selh50000rhoC);
	      hid50000rhoC.push_back(mm);
	      selh50000rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist50000rhoC[mm];
	  }

	  // head halo for this cell in 5000rhoC
	  mm = ihead5000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm5000rhoC[mm] - xi;
	    double dy = ycm5000rhoC[mm] - yi;
	    double dz = zcm5000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d5000rhoC.push_back(dist);
	      hpos5000rhoC.push_back(selh5000rhoC);
	      hid5000rhoC.push_back(mm);
	      selh5000rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist5000rhoC[mm];
	  }

	  // head halo for this cell in 2000rhoC
	  mm = ihead2000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm2000rhoC[mm] - xi;
	    double dy = ycm2000rhoC[mm] - yi;
	    double dz = zcm2000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d2000rhoC.push_back(dist);
	      hpos2000rhoC.push_back(selh2000rhoC);
	      hid2000rhoC.push_back(mm);
	      selh2000rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist2000rhoC[mm];
	  }
	  
	  // head halo for this cell in 1000rhoC
	  mm = ihead1000rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm1000rhoC[mm] - xi;
	    double dy = ycm1000rhoC[mm] - yi;
	    double dz = zcm1000rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d1000rhoC.push_back(dist);
	      hpos1000rhoC.push_back(selh1000rhoC);
	      hid1000rhoC.push_back(mm);
	      selh1000rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist1000rhoC[mm];
	  }

	  // head halo for this cell in 500rhoC
	  mm = ihead500rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm500rhoC[mm] - xi;
	    double dy = ycm500rhoC[mm] - yi;
	    double dz = zcm500rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d500rhoC.push_back(dist);
	      hpos500rhoC.push_back(selh500rhoC);
	      hid500rhoC.push_back(mm);
	      selh500rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist500rhoC[mm];
	  }

	  // head halo for this cell in 200rhoC
	  mm = ihead200rhoC[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm200rhoC[mm] - xi;
	    double dy = ycm200rhoC[mm] - yi;
	    double dz = zcm200rhoC[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d200rhoC.push_back(dist);
	      hpos200rhoC.push_back(selh200rhoC);
	      hid200rhoC.push_back(mm);
	      selh200rhoC++;	      
	    }
	    //...take next halo and repeat
	    mm = llist200rhoC[mm];
	  }

	  // head halo for this cell in 200rhoB
	  mm = ihead200rhoB[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm200rhoB[mm] - xi;
	    double dy = ycm200rhoB[mm] - yi;
	    double dz = zcm200rhoB[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d200rhoB.push_back(dist);
	      hpos200rhoB.push_back(selh200rhoB);
	      hid200rhoB.push_back(mm);
	      selh200rhoB++;	      
	    }
	    //...take next halo and repeat
	    mm = llist200rhoB[mm];
	  }
	}
      }
    }
    
    int iw50000rhoC;
    double mw50000rhoC,rw50000rhoC,dr50000rhoC;
    if(selh50000rhoC>0){
      std::sort(hpos50000rhoC.begin(),hpos50000rhoC.end(), CreateCmpPairs(d50000rhoC));
      idMih50000rhoC[i] = hid50000rhoC[hpos50000rhoC[0]];
      dMih50000rhoC[i] = d50000rhoC[hpos50000rhoC[0]];
      iw50000rhoC=idMih50000rhoC[i]+1;
      mw50000rhoC=m50000rhoC[idMih50000rhoC[i]];
      rw50000rhoC=r50000rhoC[idMih50000rhoC[i]]*boxl;
      dr50000rhoC=d50000rhoC[hpos50000rhoC[0]]*boxl;
    }
    else{
      iw50000rhoC=0;
      mw50000rhoC=0.;
      rw50000rhoC=0.;
      dr50000rhoC=0.;
    }

    int iw5000rhoC;
    double mw5000rhoC,rw5000rhoC,dr5000rhoC;
    if(selh5000rhoC>0){
      std::sort(hpos5000rhoC.begin(),hpos5000rhoC.end(), CreateCmpPairs(d5000rhoC));
      idMih5000rhoC[i] = hid5000rhoC[hpos5000rhoC[0]];
      dMih5000rhoC[i] = d5000rhoC[hpos5000rhoC[0]];
      iw5000rhoC=idMih5000rhoC[i]+1;
      mw5000rhoC=m5000rhoC[idMih5000rhoC[i]];
      rw5000rhoC=r5000rhoC[idMih5000rhoC[i]]*boxl;
      dr5000rhoC=d5000rhoC[hpos5000rhoC[0]]*boxl;
    }
    else{
      iw5000rhoC=0;
      mw5000rhoC=0.;
      rw5000rhoC=0.;
      dr5000rhoC=0.;
    }

    int iw2000rhoC;
    double mw2000rhoC,rw2000rhoC,dr2000rhoC;
    if(selh2000rhoC>0){
      std::sort(hpos2000rhoC.begin(),hpos2000rhoC.end(), CreateCmpPairs(d2000rhoC));
      idMih2000rhoC[i] = hid2000rhoC[hpos2000rhoC[0]];
      dMih2000rhoC[i] = d2000rhoC[hpos2000rhoC[0]];
      iw2000rhoC=idMih2000rhoC[i]+1;
      mw2000rhoC=m2000rhoC[idMih2000rhoC[i]];
      rw2000rhoC=r2000rhoC[idMih2000rhoC[i]]*boxl;
      dr2000rhoC=d2000rhoC[hpos2000rhoC[0]]*boxl;
    }
    else{
      iw2000rhoC=0;
      mw2000rhoC=0.;
      rw2000rhoC=0.;
      dr2000rhoC=0.;
    }
    
    int iw1000rhoC;
    double mw1000rhoC,rw1000rhoC,dr1000rhoC;
    if(selh1000rhoC>0){
      std::sort(hpos1000rhoC.begin(),hpos1000rhoC.end(), CreateCmpPairs(d1000rhoC));
      idMih1000rhoC[i] = hid1000rhoC[hpos1000rhoC[0]];
      dMih1000rhoC[i] = d1000rhoC[hpos1000rhoC[0]];
      iw1000rhoC=idMih1000rhoC[i]+1;
      mw1000rhoC=m1000rhoC[idMih1000rhoC[i]];
      rw1000rhoC=r1000rhoC[idMih1000rhoC[i]]*boxl;
      dr1000rhoC=d1000rhoC[hpos1000rhoC[0]]*boxl;
    }     
    else{
      iw1000rhoC=0;
      mw1000rhoC=0.;
      rw1000rhoC=0.;
      dr1000rhoC=0.;
    }
    
    int iw500rhoC;
    double mw500rhoC,rw500rhoC,dr500rhoC;
    if(selh500rhoC>0){
      std::sort(hpos500rhoC.begin(),hpos500rhoC.end(), CreateCmpPairs(d500rhoC));
      idMih500rhoC[i] = hid500rhoC[hpos500rhoC[0]];
      dMih500rhoC[i] = d500rhoC[hpos500rhoC[0]];
      iw500rhoC=idMih500rhoC[i]+1;
      mw500rhoC=m500rhoC[idMih500rhoC[i]];
      rw500rhoC=r500rhoC[idMih500rhoC[i]]*boxl;
      dr500rhoC=d500rhoC[hpos500rhoC[0]]*boxl;
    }    
    else{
      iw500rhoC=0;
      mw500rhoC=0.;
      rw500rhoC=0.;
      dr500rhoC=0.;
    }
    
    int iw200rhoC;
    double mw200rhoC,rw200rhoC,dr200rhoC;
    if(selh200rhoC>0){
      std::sort(hpos200rhoC.begin(),hpos200rhoC.end(), CreateCmpPairs(d200rhoC));
      idMih200rhoC[i] = hid200rhoC[hpos200rhoC[0]];
      dMih200rhoC[i] = d200rhoC[hpos200rhoC[0]];
      iw200rhoC=idMih200rhoC[i]+1;
      mw200rhoC=m200rhoC[idMih200rhoC[i]];
      rw200rhoC=r200rhoC[idMih200rhoC[i]]*boxl;
      dr200rhoC=d200rhoC[hpos200rhoC[0]]*boxl;
    }     
    else{
      iw200rhoC=0;
      mw200rhoC=0.;
      rw200rhoC=0.;
      dr200rhoC=0.;
    }
    
    int iw200rhoB;
    double mw200rhoB,rw200rhoB,dr200rhoB;
    if(selh200rhoB>0){
      std::sort(hpos200rhoB.begin(),hpos200rhoB.end(), CreateCmpPairs(d200rhoB));
      idMih200rhoB[i] = hid200rhoB[hpos200rhoB[0]];
      dMih200rhoB[i] = d200rhoB[hpos200rhoB[0]];
      iw200rhoB=idMih200rhoB[i]+1;
      mw200rhoB=m200rhoB[idMih200rhoB[i]];
      rw200rhoB=r200rhoB[idMih200rhoB[i]]*boxl;
      dr200rhoB=d200rhoB[hpos200rhoB[0]]*boxl;
    }    
    else{
      iw200rhoB=0;
      mw200rhoB=0.;
      rw200rhoB=0.;
      dr200rhoB=0.;
    }
    
    ofilout << std:: setw(7) << i+1 << "   " << std:: scientific << mvir[i] << "  " << std:: fixed
	    << std:: setprecision(6) << rvir[i]*boxl << "   "
      	    << std:: setw(7) << iw50000rhoC << "   " << std:: scientific  << mw50000rhoC << "   "
	    << std:: fixed << std:: setprecision(6)  << rw50000rhoC << "  " << dr50000rhoC << "   "
	    << std:: setw(7) << iw5000rhoC << "   " << std:: scientific  << mw5000rhoC << "   "
	    << std:: fixed << std:: setprecision(6)  << rw5000rhoC << "  " << dr5000rhoC << "   "
	    << std:: setw(7) << iw2000rhoC << "   " << std:: scientific  << mw2000rhoC << "   "
	    << std:: fixed << std:: setprecision(6)  << rw2000rhoC << "  " << dr2000rhoC << "   "
	    << std:: setw(7) << iw1000rhoC << "   " << std:: scientific  << mw1000rhoC << "   "
	    << std:: fixed << std:: setprecision(6) << rw1000rhoC << "  " << dr1000rhoC << "   "      
	    << std:: setw(7) << iw500rhoC << "   " << std:: scientific << mw500rhoC << "   "
	    << std:: fixed << std:: setprecision(6) << rw500rhoC << "  " << dr500rhoC << "   "      
	    << std:: setw(7) << iw200rhoC << "   " << std:: scientific << mw200rhoC << "   "
	    << std:: fixed << std:: setprecision(6) << rw200rhoC << "  " << dr200rhoC << "   "      
	    << std:: setw(7) << iw200rhoB << "   " << std:: scientific  << mw200rhoB << "   "
	    << std:: fixed << std:: setprecision(6) << rw200rhoB << "  " << dr200rhoB      
	    << std:: endl;
  }
  ofilout.close();
}
