//code for testing ntt, intt, pmul in FasterNTT
#include <iostream>
using namespace fntt;
using namespace std;

extern void testPoly_mod(){
  cout <<"Task: test for polynomial modulus computation."<<endl;
  int64_t mod = 576460752213245953;
  polynomial a(2049,0);
  a[0]=2047;
  for(int i=1;i<=2048;++i){
    a[i]=mod+i;
  }
  cout <<"After modulus computation: "<<endl;
  modq_poly(a,mod);
  for(int i=1;i<=2048;++i){
    cout <<a[i]<<" ";
  }
  cout <<endl;
}

extern void testPoly_add(){
  cout <<"Task: test for polynomial addition."<<endl;
  int64_t mod = 576460752213245953;
  polynomial a(2049,0);
  polynomial b(2049,0);
  a[0]=2047;
  b[0]=2047;
  for(int i=1;i<=2048;++i){
    a[i]=0;
    b[i]=i;
  }

  cout <<"After polynomial addition: "<<endl;
  add_poly(a,b,mod);
  for(int i=1;i<=2048;++i){
    cout <<a[i]<<" ";
  }
  cout <<endl;
}

extern void testPoly_mul(){
  cout <<"Task: test for polynomial multiplication."<<endl;
  int64_t mod = 576460752213245953;
  polynomial a(2049,0);
  polynomial b(2049,0);
  a[0]=2047;
  b[0]=2047;
  for(int i=1;i<=2048;++i){
    a[i]=0;
    b[i]=1;
  }
  a[1]=1;

  //NTT test
  cout <<"After polymomial multiplication: "<<endl;
  polynomial c=multi_poly_2048(a,b,2048,mod);
  for(int i=1;i<=2048;++i){
    cout <<c[i]<<" ";
  }
  cout <<endl;

}