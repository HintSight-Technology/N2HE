#include <iostream>
using namespace std;

inline int modq_32(int number, int q) {
  while (number < 0) {
    number += q;
  }
  while (number >= q) {
    number -= q;
  }
  if (number >= (q/2)) {
    number -= q;
  }
  return number;
}

// key generation algorithm
// INPUT: dimension n
// OUTPUT: random binary vector x
extern int* LWE32_KeyGen(int n){
  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int *x = new int [n];
  random_bytes(seed,SEED_LEN);
  int r = gen_bernoulli(x,len_out,seed);
  if (r == 1){
    cout <<"Error in LWE Key Generation: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in LWE Key Generation: Length" <<endl;
  }
  return x;
}

//encryption algorithm
// INPUT: modulus q, dimension n, variance 2^(-k), k >= 1 plaintext m, LWE key x
// OUTPUT：LWE ciphertext (vec(a),b)
extern vector<int> LWE32_Enc(int q, int n, int k, int m, const int* x){
  vector<int> ct;

  //generate random array a
  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int a[len_out];
  random_bytes(seed,SEED_LEN);
  int r = gen_uniform(a, len_out, q, seed);
  if (r == 1){
    cout <<"Error in generation random array a of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation random array a of LWE encryption: modulus" <<endl;
  }

  //generate error array e
  int e[8];
  random_bytes(seed,SEED_LEN);
  r = gen_ternary_var(e, 8, k, seed);
  if (r == 1){
    cout <<"Error in generation a random error of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation a random error of LWE Encryption: modulus" <<endl;
  }

  // b = ( - <a,x> + m + e) mod q
  int b=0;
  for(int i=0;i<n;++i){
    int ax = (a[i]-q/2)*x[i];
    b += ax;
    b = modq_32(b,q);
    ct.push_back(a[i]-q/2);
  }
  b *= (-1); 
  b = modq_32(b,q);
  b += (m+e[0]);
  //b += m;
  b = modq_32(b,q);
  ct.push_back(b);
  //return ciphertext
  return ct;

}


// LWE Decryption
// INPUT: modulus q, dimension n, LWE ciphertext c = (vec(a),b), LWE key x
// OUTPUT: Decryption (b + <a,x>) mod q
extern int LWE32_Dec(int q, int n, const vector<int>& c, const int* x) {
  int ip2 = c[n];
  for (int i = 0; i < n; ++i) {
    if (x[i] == -1) {
      ip2 -= c[i];
    }
    else if (x[i] == 1) {
      ip2 += c[i];
    }
    else;
    ip2 = modq_32(ip2, q);
  }
  int m = modq_32(ip2 , q);
  return m;
}


// Add two LWE ciphertexts  
// INPUT: modulus q, dimension n, LWE ciphertext ct1 and ct2
// OUTPUT: LWE ciphertext (ct1 + ct2) mod q
extern vector<int> LWE32_Add_ct(int q, int n, const vector<int>& ct1, const vector<int>& ct2) {
  vector<int> ct(n+1,0);
  for (register int i = 0; i <= n; ++i) {
    ct[i]=modq_32(ct1[i] + ct2[i],q); 
  }
  return ct;
}

// Add a LWE ciphertext and a plaintext number k  
// INPUT: modulus q, dimension n, LWE ciphertext ct1, plaintext number k
// OUTPUT: LWE ciphertext of dec(ct1) + k
extern vector<int> LWE32_Plain_Add_ct(int q, int n, const vector<int>& ct1, int k) {
  vector<int> ct=ct1;
  ct[n] = modq_32(ct[n]+k, q);
  return ct;
}


// Multiply a ciphertext and a plaintext number k
// INPUT: modulus q, dimension n, ciphertext ct1, plaintext number k
// OUTPUT: LWE ciphertext of k*dec(ct1)
extern vector<int> LWE32_Plain_Multi_ct(int q, int n, const vector<int>& ct1, int k) {
  vector<int> ct;
  int size = ct1.size();
  for (register int i = 0; i <size; ++i) {
    ct.emplace_back(modq_32(ct1[i] * k,q));
  }

  return ct;
}


// modulus Switching, mod q1 -> mod q2  
// INPUT: first modulus q1, second modulus q2, dimension n, LWE ciphertext ct with modulus q1
// OUTPUT: LWE ciphertext with modulus q2
vector<int> LWE32_Rounding(int q1, int q2, int n, const vector<int>& ct) {
  vector<int> ans=ct;
  for (register int i = 0; i <n+1; ++i) {
    ans[i] = q2 * ans[i] / q1;
  }
  return ans;
}


