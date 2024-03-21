// ================================================
// Clean up 2023 11 10

// ================================================

#include <iostream>
using namespace fntt;
using namespace std;

typedef vector<int64_t> polynomial;

// RLWE key generation algorithm
// INPUT: dimension n
// OUTPUT: a degree-(n-1) polynomial s, with coeffients from Uniform Random Distribution on {-1,0,1}
extern polynomial RLWE64_KeyGen(int n) {
  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int *x = new int [n];
  random_bytes(seed,SEED_LEN);
  int r = gen_ternary(x,len_out,seed);
  if (r == 1){
    cout <<"Error in RLWE Key Generation: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in RLWE Key Generation: Length" <<endl;
  }

  polynomial s(n+1,0);
  s[0] = (int64_t)n - 1;
  for (int i = 1; i <= n; ++i) {
    s[i] = x[i-1];
  }
  return s;
}

// RLWE encryption algorithm
// INPUT: dimension n, modulus q, variance 2^(-k), k >= 1, plaintext (polynomial) m, RLWE key s
// OUTPUT: RLWE ciphertext ct = (a, b) = (a, m + e - a * s)
extern vector<polynomial> RLWE64_Enc_2048(int n, int64_t q, int k, const polynomial & m, const polynomial & s) {

  vector<polynomial> ct;

  //generate random a
  double logq=log(q)/log(2.0);
  int int_logq = (int) logq;
  if(logq > (double)int_logq){
    int_logq++;
  }

  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int64_t *array_a = new int64_t [len_out];
  random_bytes(seed,SEED_LEN);
  int r = gen_uniform_int64(array_a, len_out, q, int_logq, seed);
  if (r == 1){
    cout <<"Error in generation random array a of RLWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation random array a of RLWE encryption: modulus" <<endl;
  }

  polynomial a(n+1,0);
  a[0] = (int64_t)n - 1;
  for (int i = 1; i <= n; ++i) {
    a[i]=array_a[i-1]-q/2;
   // a[i]=1;
  }
  ct.push_back(a);

  //compute - a * s
  polynomial as = multi_poly_2048(s,a,n,q);
  multi_scale_poly(-1, as,q);

  //generate  error e
  //generate error array e
  int *array_e = new int [len_out];
  random_bytes(seed,SEED_LEN);
  r = gen_ternary_var(array_e, len_out, k, seed);
  if (r == 1){
    cout <<"Error in generation a random error of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation a random error of LWE Encryption: modulus" <<endl;
  }
  for (int i = 1; i <= m[0]+1; ++i) {
    //as[i] += m[i];
    as[i] += (m[i]+array_e[i-1]);
  }
  modq_poly(as, q);
  ct.push_back(as);

  return ct;
}

extern vector<polynomial> RLWE64_Enc_512(int n, int64_t q, int k, const polynomial & m, const polynomial & s) {

  vector<polynomial> ct;

  //generate random a
  double logq=log(q)/log(2.0);
  int int_logq = (int) logq;
  if(logq > (double)int_logq){
    int_logq++;
  }

  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int64_t *array_a = new int64_t [len_out];
  random_bytes(seed,SEED_LEN);
  int r = gen_uniform_int64(array_a, len_out, q, int_logq, seed);
  if (r == 1){
    cout <<"Error in generation random array a of RLWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation random array a of RLWE encryption: modulus" <<endl;
  }

  polynomial a(n+1,0);
  a[0] = (int64_t)n - 1;
  for (int i = 1; i <= n; ++i) {
    a[i]=array_a[i-1]-q/2;
   // a[i]=1;
  }
  ct.push_back(a);

  //compute - a * s
  polynomial as = multi_poly_512(s,a,n,q);
  multi_scale_poly(-1, as,q);

  //generate  error e
  //generate error array e
  int *array_e = new int [len_out];
  random_bytes(seed,SEED_LEN);
  r = gen_ternary_var(array_e, len_out, k, seed);
  if (r == 1){
    cout <<"Error in generation a random error of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation a random error of LWE Encryption: modulus" <<endl;
  }
  for (int i = 1; i <= m[0]+1; ++i) {
    as[i] += (m[i]+array_e[i-1]);
  }
  modq_poly(as, q);
  ct.push_back(as);

  return ct;
}

// RLWE Decryption algorithm  
// INPUT: dimension n, modulus q, RLWE key s, RLWE ciphertext ct (a, b)
// OUTPUT: polynomial b + a * s
extern polynomial RLWE64_Dec_2048(int n, int64_t q,const polynomial & s, const vector<polynomial> & ct) {
  //compute as
  polynomial as = multi_poly_2048(s, ct[0],n,q);
  //compute b+as
  add_poly(as, ct[1],q);

  return as;
}

extern polynomial RLWE64_Dec_512(int n, int64_t q,const polynomial & s, const vector<polynomial> & ct) {
  //compute as
  polynomial as = multi_poly_512(s, ct[0],n,q);
  //compute b+as
  add_poly(as, ct[1],q);

  return as;
}



//sample RGSW ciphertext from two extended RLWE encryption
// INPUT: dimension n, modulus q, k = log(q), variance 2^(-var), var >= 1, plaintext (polynomial) m, RLWE key s, decomposition base b, log(b)
// OUTPUT: RGSW ciphertext = (~RLWE(m),~RLWE(ms)). They are in the same vector.
vector<vector<polynomial>> RGSWct_2048(int n, int64_t q, int k, int var, const polynomial & m, const polynomial & s, int64_t b, int logb) {

  vector<vector<polynomial>> rgswct;

  //~RLWE(m,2^i)
  int64_t temp2k = b;
  polynomial mi = copy(m);
  for (int i = 0; i < (k/logb); ++i) {
    if (i > 0) {
      multi_scale_poly(temp2k, mi, q);
    }
    vector<polynomial> c1i = RLWE64_Enc_2048(n, q, var, mi, s);
    rgswct.push_back(c1i);
  }

  polynomial ms = multi_poly_2048(s,m,n,q);

  //~RLWE(ms,2^i)
  for (int i = 0; i < (k/logb); ++i) {
    if (i > 0) {
      multi_scale_poly(temp2k, ms, q);
    }
    vector<polynomial> c2i = RLWE64_Enc_2048(n, q, var, ms, s);
    rgswct.push_back(c2i);
  }

  return rgswct;
}


//Operator <> (known polynomial r multiply RGSW ciphertext)
// INPUT: dimension n, modulus q, k = log(q), polynomial r, ~RLWE ciphertext ct, decomposition base b, log(b)
// OUTPUT: r <> ct := (b-decomption of r) slot-wise-multiply (ct)
vector<polynomial> bit_then_multiply_2048(int n, int64_t q, int k, const polynomial & r, const vector<vector<polynomial>> & ct, int b, int logb) {
  //b-decomption of r
  vector<polynomial> rbit = bit_poly(k, r,q,b,logb);

  polynomial ip1 = multi_poly_2048(ct[0][0], rbit[0],n,q);

  polynomial ip2 = multi_poly_2048(ct[0][1], rbit[0],n,q);

  //compute ri*cti
  for (int i = 1; i < (k/logb); ++i) {
    polynomial temp1 = multi_poly_2048(ct[i][0], rbit[i],n,q);
    add_poly(ip1, temp1,q);

    polynomial temp2 = multi_poly_2048(ct[i][1], rbit[i],n,q);
    add_poly(ip2, temp2,q);
  }
  //return (ip1,ip2)
  vector<polynomial> ans;
  ans.push_back(ip1);
  ans.push_back(ip2);

  return ans;
}

// For key switch algorithm
// INPUT: dimension n, modulus q, k = log(q), polynomial r, ~RLWE ciphertext ct, decomposition base b, log(b)
// OUTPUT: r <> ct := (b-decomption of r) slot-wise-multiply (ct)
vector<polynomial> bit_then_multiply_512(int n, int64_t q, int k, const polynomial & r, const vector<vector<polynomial>> & ct, int b, int logb) {


  //b-decomption of r
  vector<polynomial> rbit = bit_poly(k, r,q,b,logb);

  polynomial ip1 = multi_poly_512(ct[0][0], rbit[0],n,q);

  polynomial ip2 = multi_poly_512(ct[0][1], rbit[0],n,q);

  //compute ri*cti
  for (int i = 1; i < (k/logb); ++i) {
    polynomial temp1 = multi_poly_512(ct[i][0], rbit[i],n,q);
    add_poly(ip1, temp1,q);

    polynomial temp2 = multi_poly_512(ct[i][1], rbit[i],n,q);
    add_poly(ip2, temp2,q);
  }
  //return (ip1,ip2)
  vector<polynomial> ans;
  ans.push_back(ip1);
  ans.push_back(ip2);

  return ans;
}


// External product, RLWE * RGSW -> RLWE
// INPUT: dimension n, modulus q, k = log(q), 
// INPUT: RLWE ciphertext rlwe_ct = (a, b), RGSW ciphertext rgsw_ct = (beta, alpha), RLWE key s,
// INPUT: decomposition base b, log(b)
// OUTPUT: an RLWE ciphertext := (a <> alpha + b <> beta)
vector<polynomial> rlwe_multi_rgsw(int n, int64_t q, int k, const vector<polynomial> & rlwe_ct, const vector<vector<polynomial>> & rgsw_ct, int64_t b, int logb) {

  //create beta,alpha from rgsw_ct
  vector<vector<polynomial>> beta(k/logb);
  for (int i = 0; i < (k/logb); ++i) {
    beta[i]=rgsw_ct[i];
  }
  vector<vector<polynomial>> alpha(k/logb);
  for (int i = (k/logb); i < 2*(k/logb); ++i) {
    alpha[i-k/logb]=rgsw_ct[i];
  }


  //compute b <> beta and a <> alpha
  vector<polynomial> b_beta;
  vector<polynomial> a_alpha;
  if(n == 2048){
    b_beta = bit_then_multiply_2048(n, q, k, rlwe_ct[1], beta,b,logb);
    a_alpha = bit_then_multiply_2048(n, q, k, rlwe_ct[0], alpha,b,logb);
  }
  else if (n == 512){
    b_beta = bit_then_multiply_512(n, q, k, rlwe_ct[1], beta,b,logb);
    a_alpha = bit_then_multiply_512(n, q, k, rlwe_ct[0], alpha,b,logb);
  }
  else{
    cout <<"undefined input n in rlwe_multi_rgsw. "<<endl;
  }
  //compute b_beta + a_alpha
  vector<polynomial> ans(2);
  add_poly(b_beta[0], a_alpha[0],q);
  ans[0]=b_beta[0];

  add_poly(b_beta[1], a_alpha[1],q);
  ans[1]=b_beta[1];

  modq_poly(ans[0],q);
  modq_poly(ans[1],q);
  return ans;
}


// Extract the const term of an RLWE encrypted polynomial, return LWE ciphertext 
// INPUT: modulus q, scaler delta, dimension n, RLWE ciphertext RLWE_ct of polynomial k
// OUTPUT: LWE ciphertext of the const term of k
extern vector<int64_t> Extract0(int64_t q, int64_t delta, int n, const vector<polynomial> & RLWE_ct) {

  if(RLWE_ct[0][0] == 0){
    
    vector<int64_t> b(n,0);
    b.push_back(RLWE_ct[1][1]);
    return b;
  }
  else{
    vector<int64_t> a(RLWE_ct[0].begin()+2,RLWE_ct[0].end());
    int size = a.size();
    for (int i = 0; i <size; ++i) {
      a[i] *= -1;
    }
    a.push_back(RLWE_ct[0][1]);

    a.push_back(RLWE_ct[1][1]);

    return a;
  }
}

// ====================================================== Key Switching ===================================================

// KeySwitching based on ~RLWE

extern vector<vector<vector<polynomial>>> LWE_ks_key(int n1, int n2, int64_t q, int k, int var, const int* x, const polynomial& s, int64_t b, int logb) {
  // INPUT: dimension n1, n2 s.t. KeySwitch n1->n2, modulus q, k = log(q)
  // INPUT: LWE key (vec) x, ~RLWE key (polynomial) s 
  // INPUT: decomposition base b, log(b) 

  // OUTPUT: KeySwitching key, a vector of ~RLWE ciphertexts

  vector<vector<vector<polynomial>>> ans;
  int N = n1 / n2;
  if (N * n2 < n1) N++;


    //~RLWE
    

  for (int j=0;j<N;++j){
    vector<vector<polynomial>> temp_ks_key;
    int64_t temp2k = b;

    polynomial mi(n2 + 1,0);
    mi[0] = n2-1;
    
    for (int l=0;l<n2;l++){
      if (j*n2+l < n1) {
        mi[l+1] = x[j*n2+l];
      }
    }

    for (int i = 0; i < (k/logb); ++i) {
      if (i > 0) {
        multi_scale_poly(temp2k, mi, q);
      }

      vector<polynomial> c1i;
      if(n2 == 512){
        c1i = RLWE64_Enc_512(n2, q, var, mi, s);
      }
      else{
        cout <<"undefined input n2 in LWE_ks_key. "<<endl;
      }

      temp_ks_key.push_back(c1i);

    } 

    // end of ~RLWE
    ans.push_back(temp_ks_key);

  }
  return ans;

}

extern vector<int64_t> LWE_ks(int n1, int n2, int64_t q, int k, const vector<int64_t>& lwe_ct, const vector<vector<vector<polynomial>>>& ks_key, int64_t b, int logb){
  // INPUT: dimension n1, n2 s.t. KeySwitch n1->n2, modulus q, k = log(q)
  // INPUT: LWE ciphertext (vec(a),b) with n1, key x, modulus q
  // INPUT: KeySwitching key ks_key, LWE key (vec) x, ~RLWE key (polynomial) s 
  // INPUT: decomposition base b, log(b) 

  // OUTPUT: LWE ciphertext (vec(a)',b') with n2, key s, modulus q 


  vector<polynomial> temp_set_a;
  int N = n1 / n2;
  if (N * n2 < n1) N++;

  for (int j=0;j<N;j++){
    
    polynomial temp_a(n2+1,0);
    temp_a[0] = n2-1;
    temp_a[1] = lwe_ct[j*n2];
    for (int l = 1; l < n2; l++){
      temp_a[n2 - l + 1] -= lwe_ct[j * n2 + l];
    }
    temp_set_a.push_back(temp_a);

  }

  vector<polynomial> rlwe_ct;
  if(n2 == 512){
    rlwe_ct = bit_then_multiply_512( n2,  q,  k,  temp_set_a[0],  ks_key[0],  b,  logb);
    for (int j=1;j<N;j++){
      vector<polynomial> rlwe_ct_j = bit_then_multiply_512( n2,  q,  k, temp_set_a[j], ks_key[j],  b,  logb);
      add_poly(rlwe_ct[0],rlwe_ct_j[0],q);
      add_poly(rlwe_ct[1],rlwe_ct_j[1],q);
    }
  }
  
  else{
    cout <<"undefined input n2 in LWE_ks. "<<endl;
  }

  vector<int64_t> ans = Extract0( q, 0, n2, rlwe_ct);

   ans[n2] += lwe_ct[n1];
   ans[n2] = modq_64(ans[n2],q);

  return ans;

}



