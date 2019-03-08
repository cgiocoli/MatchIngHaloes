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
   vesion MiHaloes-1.6
   modidied 11 - 09 - 2017 (cgiocoli@gmail.com) read 5 overdensity cats and match according to the virial one
                                                in this way can be generalized
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
  std:: string pathvir,path1,path2,path3,path4;
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
    filinput >> path1;
    filinput >> sbut;
    filinput >> path2;
    filinput >> sbut;
    filinput >> path3;
    filinput >> sbut;
    filinput >> path4;
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
  std:: string filvir,fil1,fil2,fil3,fil4;
  std:: string ssnap;

  if(snap<10) ssnap = "00"+conv(snap,fINT);
  else if(snap>=10 && snap<100) ssnap = "0"+conv(snap,fINT);
  else ssnap = conv(snap,fINT);

  // set the filenames
  filvir = pathvir + "/SO." + ssnap + ".gv";
  fil1 = path1 + "/SO." + ssnap + ".gv";  
  fil2 = path2 + "/SO." + ssnap + ".gv";  
  fil3 = path3 + "/SO." + ssnap + ".gv";  
  fil4 = path4 + "/SO." + ssnap + ".gv";  
  
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

  f = fexists(fil1);
  if(f){
    std:: cout << fil1 << std:: endl;
  }else{
    std:: cout << fil1 << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }
  
  f = fexists(fil2);
  if(f){
    std:: cout << fil2 << std:: endl;
  }else{
    std:: cout << fil2 << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil3);
  if(f){
    std:: cout << fil3 << std:: endl;
  }else{
    std:: cout << fil3 << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }

  f = fexists(fil4);
  if(f){
    std:: cout << fil4 << std:: endl;
  }else{
    std:: cout << fil4 << " does not exist " << std:: endl;
    std:: cout << " I will STOP here !!! " << std:: endl;
    exit(1);
  }
  
  std:: ifstream ifilvir,ifil1,ifil2,ifil3,ifil4;
  
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

  std:: cout << " 1 " << std:: endl;
  
  // open the first catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m1,r1,xcm1,ycm1,zcm1;
  std:: vector<int> idh1;  
  ifil1.open(fil1.c_str());
  while(ifil1 >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh1.push_back(idh);
    m1.push_back(m);
    r1.push_back(r/boxl/1.e+3);
    xcm1.push_back(x/boxl);
    ycm1.push_back(y/boxl);
    zcm1.push_back(z/boxl);
  }
  ifil1.close();  

  int ***ihead1;
  int *llist1;
  // Allocate memory
  llist1 = new int[xcm1.size()];
  ihead1 = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead1[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead1[i][j] = new int[nlist];
  }
  makeLlist(xcm1,ycm1,zcm1,dl,llist1,ihead1);  

  std:: cout << " 2 " << std:: endl;
  
  // open the second catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m2,r2,xcm2,ycm2,zcm2;
  std:: vector<int> idh2;  
  ifil2.open(fil2.c_str());
  while(ifil2 >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh2.push_back(idh);
    m2.push_back(m);
    r2.push_back(r/boxl/1.e+3);
    xcm2.push_back(x/boxl);
    ycm2.push_back(y/boxl);
    zcm2.push_back(z/boxl);
  }
  ifil2.close();  

  int ***ihead2;
  int *llist2;
  // Allocate memory
  llist2 = new int[xcm2.size()];
  ihead2 = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead2[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead2[i][j] = new int[nlist];
  }
  makeLlist(xcm2,ycm2,zcm2,dl,llist2,ihead2);  
  
  std:: cout << " 3 " << std:: endl;
  
  // open the third catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m3,r3,xcm3,ycm3,zcm3;
  std:: vector<int> idh3;  
  ifil3.open(fil3.c_str());
  while(ifil3 >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh3.push_back(idh);
    m3.push_back(m);
    r3.push_back(r/boxl/1.e+3);
    xcm3.push_back(x/boxl);
    ycm3.push_back(y/boxl);
    zcm3.push_back(z/boxl);
  }
  ifil3.close();  

  int ***ihead3;
  int *llist3;
  // Allocate memory
  llist3 = new int[xcm3.size()];
  ihead3 = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead3[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead3[i][j] = new int[nlist];
  }
  makeLlist(xcm3,ycm3,zcm3,dl,llist3,ihead3);

  std:: cout << " 4 " << std:: endl;

  // open the fourth catalogue and read ID, Mass, Radius and CM pos
  std:: vector<double> m4,r4,xcm4,ycm4,zcm4;
  std:: vector<int> idh4;  
  ifil4.open(fil4.c_str());
  while(ifil4 >> idh >> nsb >> n >> msb >> m >> r >> ibut >> rbut >> x >> y >> z >> vx >> vy >> vz >> vdisp >> etot  >> idmb >> spin){
    idh4.push_back(idh);
    m4.push_back(m);
    r4.push_back(r/boxl/1.e+3);
    xcm4.push_back(x/boxl);
    ycm4.push_back(y/boxl);
    zcm4.push_back(z/boxl);
  }
  ifil4.close();  

  int ***ihead4;
  int *llist4;
  // Allocate memory
  llist4 = new int[xcm4.size()];
  ihead4 = new int**[nlist];
  for (int i = 0; i < nlist; ++i) {
    ihead4[i] = new int*[nlist];
    for (int j = 0; j < nlist; ++j)
      ihead4[i][j] = new int[nlist];
  }
  makeLlist(xcm4,ycm4,zcm4,dl,llist4,ihead4);
    
  std:: cout << " all linked lists done " << std:: endl;

  // number of haloes in the catalogues
  int nhvir = idhvir.size();
  int nh1 = idh1.size();
  int nh2 = idh2.size();
  int nh3 = idh3.size();  
  int nh4 = idh4.size();

  std:: vector<int> idMih1(nhvir);
  std:: vector<int> dMih1(nhvir);
  std:: vector<int> idMih2(nhvir);
  std:: vector<int> dMih2(nhvir);
  std:: vector<int> idMih3(nhvir);
  std:: vector<int> dMih3(nhvir);
  std:: vector<int> idMih4(nhvir);
  std:: vector<int> dMih4(nhvir);

  std:: string filout;
  filout = "MatchingHaloes-1.6_." + ssnap + ".gv";
  std:: ofstream ofilout;
  ofilout.open(filout.c_str());

  for(int i=0;i<nhvir;i++){
    // for each halo compute the distance of the haloes in the other catalogue
    std:: vector<double> d1;
    std:: vector<int> hpos1,hid1;
    std:: vector<double> d2;
    std:: vector<int> hpos2,hid2;
    std:: vector<double> d3;
    std:: vector<int> hpos3,hid3;
    std:: vector<double> d4;
    std:: vector<int> hpos4,hid4;
    int selh1=0;
    int selh2=0;
    int selh3=0;    
    int selh4=0;
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

    std:: cout << " searching for neighbours of halo vir " << i << "/" << nhvir << " in other catalogues " << std:: endl;
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

	  // head halo for this cell in 1
	  int mm = ihead1[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm1[mm] - xi;
	    double dy = ycm1[mm] - yi;
	    double dz = zcm1[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d1.push_back(dist);
	      hpos1.push_back(selh1);
	      hid1.push_back(mm);
	      selh1++;	      
	    }
	    //...take next halo and repeat
	    mm = llist1[mm];
	  }

	  // head halo for this cell in 2
	  mm = ihead2[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm2[mm] - xi;
	    double dy = ycm2[mm] - yi;
	    double dz = zcm2[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d2.push_back(dist);
	      hpos2.push_back(selh2);
	      hid2.push_back(mm);
	      selh2++;	      
	    }
	    //...take next halo and repeat
	    mm = llist2[mm];
	  }

	  // head halo for this cell in 3
	  mm = ihead3[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm3[mm] - xi;
	    double dy = ycm3[mm] - yi;
	    double dz = zcm3[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d3.push_back(dist);
	      hpos3.push_back(selh3);
	      hid3.push_back(mm);
	      selh3++;	      
	    }
	    //...take next halo and repeat
	    mm = llist3[mm];
	  }
	  
	  // head halo for this cell in 4
	  mm = ihead4[ji][ki][li];
	  while(mm>-1){
	    double dx = xcm4[mm] - xi;
	    double dy = ycm4[mm] - yi;
	    double dz = zcm4[mm] - zi;
	    double dist = sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist<=1.5*rvir[i]){
	      d4.push_back(dist);
	      hpos4.push_back(selh4);
	      hid4.push_back(mm);
	      selh4++;	      
	    }
	    //...take next halo and repeat
	    mm = llist4[mm];
	  }
	}
      }
    }
  
    int iw1;
    double mw1,rw1,dr1;
    if(selh1>0){
      std::sort(hpos1.begin(),hpos1.end(), CreateCmpPairs(d1));
      idMih1[i] = hid1[hpos1[0]];
      dMih1[i] = d1[hpos1[0]];
      iw1=idMih1[i]+1;
      mw1=m1[idMih1[i]];
      rw1=r1[idMih1[i]]*boxl;
      dr1=d1[hpos1[0]]*boxl;
    }
    else{
      iw1=0;
      mw1=0.;
      rw1=0.;
      dr1=0.;
    }

    int iw2;
    double mw2,rw2,dr2;
    if(selh2>0){
      std::sort(hpos2.begin(),hpos2.end(), CreateCmpPairs(d2));
      idMih2[i] = hid2[hpos2[0]];
      dMih2[i] = d2[hpos2[0]];
      iw2=idMih2[i]+1;
      mw2=m2[idMih2[i]];
      rw2=r2[idMih2[i]]*boxl;
      dr2=d2[hpos2[0]]*boxl;
    }
    else{
      iw2=0;
      mw2=0.;
      rw2=0.;
      dr2=0.;
    }

    int iw3;
    double mw3,rw3,dr3;
    if(selh3>0){
      std::sort(hpos3.begin(),hpos3.end(), CreateCmpPairs(d3));
      idMih3[i] = hid3[hpos3[0]];
      dMih3[i] = d3[hpos3[0]];
      iw3=idMih3[i]+1;
      mw3=m3[idMih3[i]];
      rw3=r3[idMih3[i]]*boxl;
      dr3=d3[hpos3[0]]*boxl;
    }
    else{
      iw3=0;
      mw3=0.;
      rw3=0.;
      dr3=0.;
    }
    
    int iw4;
    double mw4,rw4,dr4;
    if(selh4>0){
      std::sort(hpos4.begin(),hpos4.end(), CreateCmpPairs(d4));
      idMih4[i] = hid4[hpos4[0]];
      dMih4[i] = d4[hpos4[0]];
      iw4=idMih4[i]+1;
      mw4=m4[idMih4[i]];
      rw4=r4[idMih4[i]]*boxl;
      dr4=d4[hpos4[0]]*boxl;
    }     
    else{
      iw4=0;
      mw4=0.;
      rw4=0.;
      dr4=0.;
    }
            
    ofilout << std:: setw(7) << i+1 << "   " << std:: scientific << mvir[i] << "  " << std:: fixed
	    << std:: setprecision(6) << rvir[i]*boxl << "   "
      	    << std:: setw(7) << iw1 << "   " << std:: scientific  << mw1 << "   "
	    << std:: fixed << std:: setprecision(6)  << rw1 << "  " << dr1 << "   "
	    << std:: setw(7) << iw2 << "   " << std:: scientific  << mw2 << "   "
	    << std:: fixed << std:: setprecision(6)  << rw2 << "  " << dr2 << "   "
	    << std:: setw(7) << iw3 << "   " << std:: scientific  << mw3 << "   "
	    << std:: fixed << std:: setprecision(6)  << rw3 << "  " << dr3 << "   "
	    << std:: setw(7) << iw4 << "   " << std:: scientific  << mw4 << "   "
	    << std:: fixed << std:: setprecision(6) << rw4 << "  " << dr4 
	    << std:: endl;
  }
  ofilout.close();
}
