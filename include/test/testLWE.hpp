//code for testing ntt, intt, pmul in FasterNTT
#include <iostream>
using namespace std;

extern void testLWEKeyGen(){
  cout <<"Task: Test for LWE32 key generation with length 512. "<<endl;
  int *x=LWE32_KeyGen(512);
  int num_zero = 0;
  int num_one = 0;
  for(int i=0;i<512;++i){
    if(x[i] == 0) num_zero++;
    if(x[i] == 1) num_one++;
  }
  cout<<"Number of 0: "<<num_zero<<", number of one: "<<num_one<<", total: "<<num_one+num_zero <<endl;
}

extern void testLWE32EncDec(){
  cout <<"Task: Test for LWE encryption and decryption with 32-bit integer. "<<endl;
  int *x=LWE32_KeyGen(512);
  int q = 8192;
  int pt=200;
  int k=3;
  cout <<"Modulus q = "<<q <<", length n = 512, variance k = 2^-"<<k <<endl;
  vector<int> ct = LWE32_Enc(q,512,k,pt,x);
  int m=LWE32_Dec(q,512,ct,x);
  cout << "Plaintext = "<<pt <<", decryption result (supposed be plaintext + noise) = "<<m<<endl;
}

extern void testLWE64EncDec(){
  cout <<"Task: Test for LWE encryption and decryption with 64-bit integer. "<<endl;
  int *x=LWE64_KeyGen(512);
  int64_t q = 1099511627776+10;
  int64_t pt= 4294967296+3;
  int k=3;
  cout <<"Modulus q = "<<q<<", length n = 512, variance k = 2^-"<<k  <<endl;
  vector<int64_t> ct = LWE64_Enc(q,512,k,pt,x);
  int64_t m=LWE64_Dec(q,512,ct,x);
  cout << "Plaintext = "<<pt <<", decryption result (supposed be plaintext + noise) = "<<m<<endl;
}

extern void testLWE32addct(){
  cout <<"Task: Test for homomorphic addition between two ciphertexts."<<endl;
  int *x=LWE32_KeyGen(512);
  int q = 8192;
  int pt1=200;
  int pt2 = 400;
  int k=3;
  cout <<"Modulus q = "<<q <<", length n = 512, variance k = 2^-"<<k <<endl;
  vector<int> ct1 = LWE32_Enc(q,512,k,pt1,x);
  vector<int> ct2 = LWE32_Enc(q,512,k,pt2,x); 
  vector<int> ctr = LWE32_Add_ct(q,512,ct1,ct2);
   int m=LWE32_Dec(q,512,ctr,x);
  cout << "Plaintext1 = "<<pt1 <<", Plaintext2 = "<<pt2 <<", decryption of homomorphic addition (supposed to be plaintext1 + plaintext2 + noise) = "<<m<<endl;
}

extern void testLWE32ptmulct(){
  cout <<"Task: Test for homomorphic multiplication between one plaintext and one ciphertext."<<endl;
  int *x=LWE32_KeyGen(512);
  int q = 8192;
  int pt1 = 200;
  int pt2 = 5;
  int k=1;
  cout <<"Modulus q = "<<q <<", length n = 512, variance k = 2^-"<<k <<endl;
  vector<int> ct1 = LWE32_Enc(q,512,k,pt1,x); 
  vector<int> ctr = LWE32_Plain_Multi_ct(q,512,ct1,pt2);
   int m=LWE32_Dec(q,512,ctr,x);
  cout << "Plaintext1 = "<<pt1 <<", Plaintext2 = "<<pt2 <<", decryption of homomorphic multiplication between encrypted plaintext1 and unencrypted plaintext2 (supposed to be (plaintext1 * plaintext2) + noise) = "<<m<<endl;
}