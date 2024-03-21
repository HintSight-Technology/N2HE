// ================================================
// Clean up 2023 11 10
// Test FasterNTT LUT time
// Single Thread
// Google Benchmark: NO
// ================================================

#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

const int64_t mod = 576460752213245953;     //2^59
const int64_t delta = 1099511627776;  //2^40

int main(){

  cout <<"Task: Test the correctness and time cost for homomorphic look-up table algorithm. "<<endl;
  cout << "----------" << endl;

  //set n,q,t for LWE scheme
  cout << "Parameters for LWE scheme:" << endl;
  const int n1 = 512;                        //vector length
  const int q1 = 1<<13;                      //ciphertext modulus
  const int t1 = 1<<11;                      //plaintext modulus
  const int var1 = 3;                        //-log() of error variance 
  

  cout << "vector length = " << n1 <<", ciphertext modulus = " << q1 <<", variance of error = 2^ -"<<var1 << endl;
  
  //set n,q,k for RLWE scheme
  cout << "Parameters for RLWE scheme:" << endl;
  const int n2 = 2048; //vector length=exp of polynomial=2048
  const int n3 = n2;
  const int k2 = 60;  //k=log(q)
  const int var2 = 3; //-log() of error variance 
  const int64_t b = 1<<30;
  const int logb = 30; //logb=log2(b)-1
  const int64_t q2 = mod; //ciphertext modulus
  cout << "polynomial degree = " << n2 << ", ciphertext modulus: " << q2 <<", variance of error = 2^ -"<<var2  <<", external product based on " << logb << " bits extension" << endl;
  cout << "----------" << endl;

  using R = Rq<n2, 1, FNTT, AVX2>;

  int *x = new int [n1];
  x=LWE32_KeyGen(n1);

  // ======================== Evaluation Key Generation ======================================

  //encrypt the RGSW key
  polynomial s = RLWE64_KeyGen(n2);
  int *lwe_s = new int [n2];
  for (int i = (int)s[0] + 1; i >= 1; --i) {
    lwe_s[(int)s[0] + 1 - i]=s[i];
  }

  //ntt ek0,ek1,ek2
   vector<vector<R>> _ek0,_ek1,_ek2;

  for (int i = 0; i < n1/2; ++i) {
    vector<R> tempek;
    vector<vector<polynomial>> ek0i= key_encrypt_1(n2, q2, k2, var2, x[2*i],x[2*i+1], s, b, logb);
    for(int j=0;j<2*(k2/logb);++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=n2;++h){
          if(ek0i[j][k][h]<0){
            ek0i[j][k][h]+=q2;
            tempp(0,h-1)=(uint64_t)ek0i[j][k][h];
          }
          else{
            tempp(0,h-1)=(uint64_t)ek0i[j][k][h];
          }
        }
        tempp.ntt();
        tempek.push_back(tempp);
      }
    }
    _ek0.push_back(tempek);
  }

  for (int i = 0; i < n1/2; ++i) {
    vector<R> tempek;
    vector<vector<polynomial>> ek1i = key_encrypt_2(n2, q2, k2,var2, x[2*i],x[2*i+1], s, b, logb);
    for(int j=0;j<2*(k2/logb);++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=n2;++h){
          if(ek1i[j][k][h]<0){
            ek1i[j][k][h] += q2;
            tempp(0,h-1)=(uint64_t)ek1i[j][k][h];
          }
          else{
            tempp(0,h-1)=(uint64_t)ek1i[j][k][h];
          }
        }
        tempp.ntt();
        tempek.push_back(tempp);
      }
    }
    _ek1.push_back(tempek);
  }

  for (int i = 0; i < n1/2; ++i) {
    vector<R> tempek;
    vector<vector<polynomial>> ek2i = key_encrypt_3(n2, q2, k2,var2, x[2*i], x[2*i+1], s, b, logb);
    for(int j=0;j<2*(k2/logb);++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=n2;++h){
          if(ek2i[j][k][h]<0){
            ek2i[j][k][h]+=q2;
            tempp(0,h-1)=(uint64_t)ek2i[j][k][h];
          }
          else{
            tempp(0,h-1)=(uint64_t)ek2i[j][k][h];
          }
        }
        tempp.ntt();
        tempek.push_back(tempp);
      }
    }
    _ek2.push_back(tempek);
  }
  
  cout<<"Generate and store evaluation keys in RNS form"<<endl;
  cout << "----------" << endl;

// ======================== Start Testing ======================================

  int pt[10] = {-200, -100, 0, 100, 200, 300, 400, 500, 600, 700};
  cout <<"Evaluation function: ReLU, i.e., f(x) = x if x > 0, f(x) = 0 if x <= 0."<<endl;
  polynomial f = p_f0(n2, delta, q1);
  
  struct timeval tstart1, tend1;
  gettimeofday(&tstart1,NULL);

  for (int i = 0; i < 10; ++i) {
      cout <<"Round: "<<i+1<<", Plaintext = "<<pt[i]<<endl;
      vector<int> ct = LWE32_Enc(q1, n1,var1, pt[i], x);
      cout <<"Plaintext is encrypted."<<endl;

      vector<polynomial> relu=LUT_2048<2048,1,FNTT,AVX2>(i,q1,q2, n2, delta, k2, n1,ct, _ek0, _ek1,_ek2, f,b,logb);
      relu.pop_back();

      vector<int64_t> ans2 = Extract0(q2,delta, n2, relu);

      cout <<"LUT result: "<<LWE64_Dec(q2,n2,ans2,lwe_s)/delta<<endl;
  }
  gettimeofday(&tend1,NULL);

  double time_mod1 = (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);
  time_mod1 = time_mod1/10.0;     
  cout <<"Time for one ReLU function evaluation: "<<time_mod1<<"us. "<<endl;

  return 0;
}