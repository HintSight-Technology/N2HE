
#include <iostream>
using namespace std;

typedef vector<int64_t> polynomial;

extern void testRLWEKeyGen(){
  cout <<"Task: Test for RLWE64 key generation with length 2048. "<<endl;
  polynomial x=RLWE64_KeyGen(2048);
  int num_zero = 0;
  int num_one = 0;
  int num_minusone=0;
  for(int i=1;i<=2048;++i){
    if(x[i] == 0) num_zero++;
    if(x[i] == 1) num_one++;
    if(x[i] == -1) num_minusone++;
  }
  cout<<"Number of 0: "<<num_zero<<", number of 1: "<<num_one<<", number of -1: "<<num_minusone<<", total: "<<num_minusone+num_one+num_zero <<endl;
}

extern void testRLWE64EncDec(){
  cout <<"Task: Test for RLWE encryption and decryption with 64-bit integer. "<<endl;
  polynomial x=RLWE64_KeyGen(2048);
  int64_t q = 576460752213245953;
  polynomial pt(2049,1);
  pt[0]=2047;
  int k=3;
  cout <<"Modulus q = "<<q<<", length n = 2048, variance k = 2^-"<<k  <<endl;
  vector<polynomial> ct = RLWE64_Enc_2048(2048, q,k,pt,x);
  //cout <<"Encrypted. "<<endl;
  polynomial m=RLWE64_Dec_2048(2048,q,x,ct);
  cout <<"Plaintext is a 1+x+x^2+...+x^2047"<<endl;
  cout << "Coefficients of decryption result: "<<endl;
  for(int i=1;i<=m[0];++i){
    cout<<m[i]<<" ";
  }
  cout <<endl;
}


extern void testRLWEKeySwitch(){
  cout <<"Task: Test for RLWE64 key switch. "<<endl;
  int *x=LWE64_KeyGen(2048);
  int64_t q = 576460752213245953;
  int k=3;
  int64_t pt = 1000;
  vector<int64_t> ct = LWE64_Enc(q,2048,k,pt,x);
  int64_t m=LWE64_Dec(q,2048,ct,x);
  cout << "Plaintext = "<<pt <<", decryption result = "<<m <<", key length = 2048. "<<endl;

  polynomial s=RLWE64_KeyGen(512);
  int *lwe_s = new int [512];
    for (int i = 512; i >= 1; --i) {
      lwe_s[512 - i]=(int)s[i];
    }

  vector<vector<vector<polynomial>>> key=LWE_ks_key(2048,512,q,60,k,x,s,2 ,1);
  vector<int64_t> ct_lwe_s=LWE_ks(2048,512,q,60,ct,key,2,1);
  int64_t ans4 = LWE64_Dec(q, 512, ct_lwe_s, lwe_s);

  cout <<"Decryption result after key switch (2048 -> 512) = "<<ans4<<endl;

}
