// ================================================
// Clean up 2023 11 10

// ================================================
#include <iostream>
using namespace fntt;
using namespace std;

typedef vector<int64_t> polynomial;

// ====================================================== Evaluation key based on RGSW ====================================
// LWE secret key is encrypted under rlwe secret key


// RGSW encryption of x1 * x2
extern vector<vector<polynomial>> key_encrypt_1(int n, int64_t q, int k,int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2,0);

  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;

  //1x^0 polynomial
  polynomial one(2,0);

  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*x2 == 0) {
     ans = RGSWct_2048(n, q, k, var, zero, s,b, logb);
    }
   else {
     ans = RGSWct_2048(n, q, k, var, one, s,b, logb);
   }
  }
  else{
    cout <<"undefined input n in key_encrypt_1. "<<endl;
  }
  
  return ans;
}

// RGSW encryption of x1 * (1 - x2)
extern vector<vector<polynomial>> key_encrypt_2(int n, int64_t q, int k, int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2,0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2,0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*(1-x2) == 1) {
     ans = RGSWct_2048(n, q, k, var, one, s,b,logb);
    }
    else {
      ans = RGSWct_2048(n, q, k, var, zero, s,b,logb);
    }
  }
    else{
    cout <<"undefined input n in key_encrypt_2. "<<endl;
  }
  return ans;
}

// RGSW encryption of (1 - x1) * x2
extern vector<vector<polynomial>> key_encrypt_3(int n, int64_t q, int k, int var, int x1, int x2, const polynomial& s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2, 0);

  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;

  //1x^0 polynomial
  polynomial one(2, 0);

  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n==2048){
    if (x2 * (1 - x1) == 1) {
      ans = RGSWct_2048(n, q, k, var, one, s, b, logb);
    }
    else {
      ans = RGSWct_2048(n, q, k, var, zero, s, b, logb);
    }
  }
  else{
    cout <<"undefined input n in key_encrypt_3. "<<endl;
  }
  
  return ans;
}



// ======================================================  Look-up-Table  ==========================================================================

// Activation Function 
// ReLU
extern int64_t Evaluate_T(int64_t x) {
  return (x > 0 ? x : 0);
}

//construct f based on Evaluate_T
// INPUT: dimension n, scaler delta, modulus q
// OUTPUT: Blindrotate polynomial f
extern polynomial p_f(int n, int64_t delta, int64_t q) {

  // construct LUT polynomial f in R_q
  // f.exp = n - 1;
  polynomial f(n + 1,0);
  
  f[0] = (int64_t)n - 1;
  int64_t temp_coef;

  //f: divide q to n parts. 
  for (int j = 0; j < n; j++) {
    if (j == 0) {
      temp_coef = Evaluate_T(0);
    }
    else if (j < (n / 2 )) {
      temp_coef = delta*Evaluate_T(-1 * (j * q / (4 * (int64_t)n)));

    }
    else {
      temp_coef = (-1)* delta * Evaluate_T(((int64_t)n - j) * q / (4 * (int64_t)n));
    }
    f[j + 1] = temp_coef;
  }

  return f;
}

extern polynomial p_f0(int n, int64_t delta, int64_t q) {

  // construct LUT polynomial f in R_q
  // f.exp = n - 1;
  polynomial f(n + 1,0);
  
  f[0] = (int64_t)n - 1;
  int64_t temp_coef;

  //f: divide q to n parts. 
  for (int j = 0; j < n; j++) {
    if (j == 0) {
      temp_coef = Evaluate_T(0);
    }
    else if (j < (n / 2 )) {
      temp_coef = delta*Evaluate_T(-1 * (j * q / (2 * (int64_t)n)));

    }
    else {
      temp_coef = (-1)* delta * Evaluate_T(((int64_t)n - j) * q / (2 * (int64_t)n));
    }
    f[j + 1] = temp_coef;
  }

  return f;
}


// LUT
// Input LWE ciphertext ct_LWE of m
// OUTPUT RLWE ciphertext of X^m * f
template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
extern vector<polynomial> LUT_2048(int index, int q1, int64_t q, int n, int64_t delta, int k, int N, vector<int> const &ct_LWE, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek0, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek1, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek2, polynomial const& f, int64_t b, int logb) {
  // INPUT: id of the current node index
  // INPUT: modulus q1 and q, dimension n and N, scaler delta, k = log(q)
  // INPUT: Look-up LWE ciphertext ct_LWE
  // INPUT: RGSW Evaluation Key ek0,ek1,ek2, RLWE key s
  // INPUT: Blindrotate polynomial f
  // INPUT: decompostion base b, log(b)
  // OUTPUT: RLWE ciphertext X^(dec(ct_LWE)) * f

  //insert index
  using R=Rq<NN, qNum, NTT, Impl>;
  polynomial p_index(2, 0);
  p_index[1] = index;

  //change the ciphertext modulus of LWE ct from q to 2n
  
  vector<int64_t> a_n(N);
  
  int64_t b_n = 2 * (int64_t)n * ct_LWE[N] / q1;
  int64_t b_n2 = 2 * (int64_t)n * ct_LWE[N] / (q1/2);
  if (b_n2 - 2*b_n >= 1 && b_n2>0) {
    b_n++;
  }
  else if (b_n2 - 2*b_n <= -1 && b_n2 < 0) {
    b_n--;
  }
  
  for (int i = 0; i < N; ++i) {
    a_n[i]= 2 * (int64_t)n * ct_LWE[i] / (q1);
    int64_t tempani= 2 * (int64_t)n * ct_LWE[i] / (q1/2);
    if (tempani - 2*a_n[i] >= 1 && tempani>0) {
      a_n[i]++;
    }
    else if (tempani - 2*a_n[i] <= -1 && tempani < 0) {
      a_n[i]--;
    }
    
  }

  //initialize AC0=(f*X^b,0)
  vector<polynomial> AC0(2);
  //compute X^{b_n mod n}
  
  int64_t temp = -1;
  if (b_n >= 0) {
    temp = b_n;
  }
  else {
    while (temp < 0) {
      temp = b_n + n;
    }
  }
  polynomial xb((int)temp+2,0);
  xb[0] = temp;
  if (b_n >= 0) {
    //xb.exp = b_n;
    xb[(int)temp+1] = 1;
  }
  else if (b_n < 0) {
    //xb.exp = n + b_n;
    xb[(int)temp+1] = -1;
  }
  else;

  //compute AC[0]=0X^0, AC0[1]=f*X^b
  polynomial zero_p(n+1,0);
  AC0[0]=zero_p;
  AC0[1]=(multi_poly_2048(f, xb, n, q));

  for(int i = 0; i < N/2 ; ++i){

    //generate ntt(x^a-1)
    R xa1,xa2,xa3;
    XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i]+a_n[2*i+1],xa1);
    XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i],xa2);
    XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i+1],xa3);

    xa1.ntt();
    xa2.ntt();
    xa3.ntt(); 

    //compute ntt(X^a-1)*ntt(ek), store in Xa_times_ek
    vector<vector<R>> Xa_times_ek(2*(k/logb),vector<R>(2));

    for(int j1 = 0; j1 < 2*(k/logb) ; ++j1){
      for(int j2 = 0; j2 < 2 ; ++j2){
        R tempxaek1,tempxaek2,tempxaek3;
        R tempxaeksum1;
        R::pmul(tempxaek1, xa1, _ek0[i][j1*2+j2]);
        R::pmul(tempxaek2, xa2, _ek1[i][j1*2+j2]);

        R::add (tempxaeksum1,tempxaek1,tempxaek2);

        R::pmul(tempxaek3, xa3, _ek2[i][j1*2+j2]);
        R::add(Xa_times_ek[j1][j2],tempxaeksum1,tempxaek3);

      }
    }

    //decompose AC0 to 2*k/logb
    vector<polynomial> de_AC00 = bit_poly(k,AC0[0],q,b,logb);
    vector<polynomial> de_AC01 = bit_poly(k,AC0[1],q,b,logb);

    //for test, output the size after decomposition

    //ntt(de_AC00, de_AC01)
    vector<R> ntt_de_AC00(k/logb);
    vector<R> ntt_de_AC01(k/logb);
    for(int j1 = 0; j1 < k/logb ; ++j1){
      R tempac00,tempac01;
      for(int j2 = 0; j2 < n ; ++j2){
        if(j2 <= de_AC00[j1][0]){
          if(de_AC00[j1][j2+1] >= 0){
            tempac00(0,j2) = (uint64_t)de_AC00[j1][j2+1];
          }
          else{
            tempac00(0,j2) = (uint64_t)(de_AC00[j1][j2+1]+q);
          }
        }
        else{
          tempac00(0,j2) = 0;
        }
        if(j2 <= de_AC01[j1][0]){
          if(de_AC01[j1][j2+1] >= 0){
            tempac01(0,j2) = (uint64_t)de_AC01[j1][j2+1];
          }
          else{
            tempac01(0,j2) = (uint64_t)(de_AC01[j1][j2+1]+q);
          }
        }
        else{
          tempac01(0,j2) = 0;
        }
      }
      tempac00.ntt();
      tempac01.ntt();
      ntt_de_AC00[j1]=tempac00;
      ntt_de_AC01[j1]=tempac01;
    }

    //b <> beta = ntt_de_AC01 <> Xa_times_ek[first k/logb]
    vector<R> b_beta;
    R tempac01ek0,tempac01ek1;
    R::pmul(tempac01ek0,ntt_de_AC01[0],Xa_times_ek[0][0]);
    R::pmul(tempac01ek1,ntt_de_AC01[0],Xa_times_ek[0][1]);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){
      R ac01ek0,sum0;
      R::pmul(ac01ek0,ntt_de_AC01[j1],Xa_times_ek[j1][0]);
      R::add(sum0,tempac01ek0,ac01ek0);
      tempac01ek0 = sum0;

      R ac01ek1,sum1;
      R::pmul(ac01ek1,ntt_de_AC01[j1],Xa_times_ek[j1][1]);
      R::add(sum1,tempac01ek1,ac01ek1);
      tempac01ek1 = sum1;
    }

    b_beta.push_back(tempac01ek0);
    b_beta.push_back(tempac01ek1);

    //a <> alpha = ntt_de_AC00 <> Xa_times_ek[second k/logb]
    vector<R> a_alpha;
    R tempac00ek0, tempac00ek1;
    R::pmul(tempac00ek0,ntt_de_AC00[0],Xa_times_ek[k/logb][0]);
    R::pmul(tempac00ek1,ntt_de_AC00[0],Xa_times_ek[k/logb][1]);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){
      R ac00ek0,sum0;
      R::pmul(ac00ek0,ntt_de_AC00[j1],Xa_times_ek[j1+k/logb][0]);
      R::add(sum0,tempac00ek0,ac00ek0);
      tempac00ek0 = sum0;

      R ac00ek1,sum1;
      R::pmul(ac00ek1,ntt_de_AC00[j1],Xa_times_ek[j1+k/logb][1]);
      R::add(sum1,tempac00ek1,ac00ek1);
      tempac00ek1 = sum1;
    }

    a_alpha.push_back(tempac00ek0);
    a_alpha.push_back(tempac00ek1);


    //a_alpha + b_beta
    vector<R> new_AC0(2);
    R::add(new_AC0[0],a_alpha[0],b_beta[0]);
    R::add(new_AC0[1],a_alpha[1],b_beta[1]);

    //intt(new_AC0)
    new_AC0[0].intt();
    new_AC0[1].intt();


    //update AC0
    int AC00_index = 0;
    int AC01_index = 0;
    for(int j1 = 0; j1 < n ; ++j1){
      AC0[0][j1+1] += (int64_t)new_AC0[0](0,j1);
      while(AC0[0][j1+1] >= q/2){
        AC0[0][j1+1] -= q;
      }
      if(AC0[0][j1+1] != 0){
        AC00_index = j1;
      }

      AC0[1][j1+1] += (int64_t)new_AC0[1](0,j1);
      while(AC0[1][j1+1] >= q/2){
        AC0[1][j1+1] -= q;
      }
      if(AC0[1][j1+1] != 0){
        AC01_index = j1;
      }

    }
    AC0[0][0] = AC00_index;
    AC0[1][0] = AC01_index;
    modq_poly(AC0[0],q);
    modq_poly(AC0[1],q);
  }

  //return AC0
  AC0.push_back(p_index);
  return AC0;

}


// LUT
// Input LWE ciphertext ct_LWE of m
// OUTPUT RLWE ciphertext of X^m * f
template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> LUT_2048_dec_to_2(int index, int q1, int64_t q, int n, int64_t delta, int k, int N, vector<int> const &ct_LWE, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek0, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek1, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek2, polynomial const& f, int64_t b, int logb) {
  // INPUT: id of the current node index
  // INPUT: modulus q1 and q, dimension n and N, scaler delta, k = log(q)
  // INPUT: Look-up LWE ciphertext ct_LWE
  // INPUT: RGSW Evaluation Key ek0,ek1,ek2, RLWE key s
  // INPUT: Blindrotate polynomial f
  // INPUT: decompostion base b, log(b)
  // OUTPUT: RLWE ciphertext X^(dec(ct_LWE)) * f

  using R=Rq<NN, qNum, NTT, Impl>;

  // insert index
  // index: to mark the id of the current node so that multithreading can be in order
  polynomial p_index(2, 0);
  p_index[1] = index;

  // change the ciphertext modulus of LWE ct from q1 to 2n
  // (vec(a), b) -> (vec(a_n), b_n)
  
  vector<int64_t> a_n(N);
  
  int64_t b_n = 2 * (int64_t)n * ct_LWE[N] / (1<<13);
  int64_t b_n2 = 2 * (int64_t)n * ct_LWE[N] / (1 << 12);
  if (b_n2 - 2*b_n >= 1 && b_n2>0) {
    b_n++;
  }
  else if (b_n2 - 2*b_n <= -1 && b_n2 < 0) {
    b_n--;
  }
  
 for (int i = 0; i < N; ++i) {
    a_n[i]= 2 * (int64_t)n * ct_LWE[i] / (1<<13);
    int64_t tempani= 2 * (int64_t)n * ct_LWE[i] / (1 << 12);
    if (tempani - 2*a_n[i] >= 1 && tempani>0) {
      a_n[i]++;
    }
    else if (tempani - 2*a_n[i] <= -1 && tempani < 0) {
      a_n[i]--;
    }
  }

  //initialize AC0=(f*X^b,0)
  // polynomial w(2049,0);
  vector<polynomial> AC0(2);
  //compute X^{b_n}
  
  int64_t temp = -1;

  if (b_n >= 0) {

    temp = b_n;
  }
  else {

    while (temp < 0) {
      temp = b_n + n;
    }

  }
  polynomial xb((int)temp+2,0);
  xb[0] = temp;
  if (b_n >= 0) {
    //xb.exp = b_n;
    xb[(int)temp+1] = 1;
  }
  else if (b_n < 0) {
    //xb.exp = n + b_n;
    xb[(int)temp+1] = -1;
  }
  else;

  //compute AC[0]=0X^0, AC0[1]=f*X^b
  polynomial zero_p(2049,0);
  AC0[0]=zero_p;


  AC0[1]=multi_poly_2048(f, xb, n, q);

  
  
  // Trivial case: if ciphertext is 0 0 0 0 0 0 0 0 return 0
  int iszero=1;
  for(int i=0;i<N;++i){
    if(a_n[i] != 0){
      iszero=0;
      break;
    }
  }

  if(iszero==1){
    cout <<"all zero"<<endl;
    AC0.push_back(p_index);
    return AC0;
  }


  else{
    // Normal case

    R xa1,xa2,xa3;
    vector<R> ans(8);
    for (int i = 0; i < N/2; ++i) {

      //generate xa1,xa2,xa3
      
      XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i]+a_n[2*i+1],xa1);
      XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i],xa2);
      XA_minus_1_2048<2048,1,FNTT,AVX2>(n,q,a_n[2*i+1],xa3);

      xa1.ntt();
      xa2.ntt();
      xa3.ntt();
      
      //xa1,xa2,xa3 \in [0,4q)
      //_ek0,_ek1_ek2 \in [0,4q)
      //pmul_lazy: tempxaek1,tempxaek2,tempxaek3 \in [0,2q)
      //pmul: tempxaek1,tempxaek2,tempxaek3 \in [0,q)
      for(int j=0;j<8;++j){ 
        R tempxaek1,tempxaek2,tempxaek3,tempxaeksum1,tempxaeksum2;
        R::pmul(tempxaek1, xa1, _ek0[i][j]);
        R::pmul(tempxaek2, xa2, _ek1[i][j]);
        R::pmul(tempxaek3, xa3, _ek2[i][j]);


        R::add(tempxaeksum1,tempxaek1,tempxaek2);
        R::add(tempxaeksum2,tempxaeksum1,tempxaek3);

        ans[j]=tempxaeksum2;
      }

      // external product
      //decomposition AC0[0]->(_ac00,_ac01); AC0[1]->(_ac10,_ac11);

      R _ac00,_ac01,_ac10,_ac11;
      for(int ii=0;ii<2048;++ii){
        int64_t tempa;
        if(ii > AC0[0][0]){
          tempa = 0;
        }
        else{
          tempa=AC0[0][ii+1];
        }
        while (tempa<0){
          tempa += q;
        }
        int64_t tempmod=tempa%b;
        if(tempmod < 0) cout <<"dec ac wrong"<<endl;
        int64_t tempq=(tempa-tempmod) >> logb;
        _ac00(0,ii)=(uint64_t)tempq;
        _ac01(0,ii)=(uint64_t)tempmod;
      }

      for(int ii=0;ii<2048;++ii){
        int64_t tempa;
        if(ii > AC0[1][0]){
          tempa = 0;
        }
        else{
          tempa=AC0[1][ii+1];
        }
       // tempa=AC0[1][i+1];
        while (tempa<0){
          tempa += q;
        }
        int64_t tempmod=tempa%b;
        if(tempmod < 0) cout <<"dec ac wrong"<<endl;
        int64_t tempq=(tempa-tempmod)>>logb;
        _ac10(0,ii)=(uint64_t)tempq;
        _ac11(0,ii)=(uint64_t)tempmod;
      }

      //_ac00,_ac01,_ac11,_ac10 -> ntt
      _ac00.ntt();
      _ac01.ntt();
      _ac10.ntt();
      _ac11.ntt();

      // a <> alpha
      //_ac01*(ans[4],ans[5])+_ac00(ans[6],ans[7])
      R temp1,temp2;
      R aalpha1,aalpha2;

      R::pmul(temp1, _ac01, ans[4]);
      R::pmul(temp2, _ac00, ans[6]);

      //aalpha1 \in [0,4q)
      R::add(aalpha1,temp1,temp2);

      R temp3,temp4;
      R::pmul(temp3, _ac01, ans[5]);
      R::pmul(temp4, _ac00, ans[7]);

      //aalpha2 \in [0,4q)
      R::add(aalpha2,temp3,temp4);

      // b <> beta
      //_ac11*(ans[0],ans[1])+_ac10*(ans[2],ans[3])
      R bbeta1,bbeta2;

      R temp5,temp6;
      R::pmul(temp5, _ac11, ans[0]);
      R::pmul(temp6, _ac10, ans[2]);

      //bbeta1 \in [0,4q)
      R::add(bbeta1,temp5,temp6);

      R temp7,temp8;
      R::pmul(temp7, _ac11, ans[1]);
      R::pmul(temp8, _ac10, ans[3]);

      //bbeta2 \in [0,4q)
      R::add(bbeta2,temp7,temp8);

      //a<>alpha+b<>beta
      R plus1,plus2;
      //plus1,plus2 \in [0,8q)
      R::add(plus1,aalpha1,bbeta1);
      R::add(plus2,aalpha2,bbeta2);

      //plus1.intt, plus2.intt
      plus1.intt();
      plus2.intt();

      //AC0[0] += plus1, AC0[1] += plus2
     
      int AC00_index = 0;
      int AC01_index = 0;
      for(int ii=1;ii<2049;++ii){
        if(plus1(0,ii-1)>q || plus2(0,ii-1)>q){
          cout <<"intt > q"<<endl;
        }
        AC0[0][ii]+= (int64_t)plus1(0,ii-1);
        AC0[1][ii]+= (int64_t)plus2(0,ii-1);
        
        while(AC0[0][ii] >= q){
          AC0[0][ii] -= q;
        }
        while(AC0[1][ii] >= q){
          AC0[1][ii] -= q;
        }
        
        if(AC0[0][ii] != 0){
          AC00_index = ii-1;
        }
        if(AC0[1][ii] != 0){
          AC01_index = ii-1;
        }
      }
      AC0[0][0]=AC00_index;
      AC0[1][0]=AC01_index;

      modq_poly(AC0[0],q);
      modq_poly(AC0[1],q);     
    }

    //insert index
    AC0.push_back(p_index);

    return AC0;
  }
}

template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> LUT_2048_slow(int index,int q1, int64_t q, int n, int64_t delta, int k, int N, vector<int> const &ct_LWE, vector<vector<vector<polynomial>>> const &ek0, vector<vector<vector<polynomial>>> const &ek1, vector<vector<vector<polynomial>>> const& ek2, polynomial const& f, int b, int logb) {

  //insert index
  using R=Rq<NN, qNum, NTT, Impl>;
  polynomial p_index(2, 0);
  p_index[1] = index;
  

  //change the ciphertext modulus of LWE ct from q to n
  
  vector<int64_t> a_n(N);
  
  int64_t b_n = 2 * (int64_t)n * ct_LWE[N] / (1<<13);
  int64_t b_n2 = 2 * (int64_t)n * ct_LWE[N] / (1 << 12);
  if (b_n2 - 2*b_n >= 1 && b_n2>0) {
    b_n++;
  }
  else if (b_n2 - 2*b_n <= -1 && b_n2 < 0) {
    b_n--;
  }


  for (int i = 0; i < N; ++i) {    
    a_n[i]= 2 * (int64_t)n * ct_LWE[i] / (1<<13);
    int64_t tempani= 2 * (int64_t)n * ct_LWE[i] / (1 << 12);
    if (tempani - 2*a_n[i] >= 1 && tempani>0) {
      a_n[i]++;
    }
    else if (tempani - 2*a_n[i] <= -1 && tempani < 0) {
      a_n[i]--;
    }
  }

  
  //initialize AC0=(f*X^b,0)
  vector<polynomial> AC0(2);
  //compute X^{b_n mod n}
  
  int64_t temp = -1;
  if (b_n >= 0) {
    temp = b_n;
  }
  else {
    while (temp < 0) {
      temp = b_n + n;
    }
  }
  polynomial xb((int)temp+2,0);
  xb[0] = temp;
  if (b_n >= 0) {
    //xb.exp = b_n;
    xb[(int)temp+1] = 1;
  }
  else if (b_n < 0) {
    //xb.exp = n + b_n;
    xb[(int)temp+1] = -1;
  }
  else;
  
  //compute AC[0]=0X^0, AC0[1]=f*X^b
  polynomial zero_p(2,0);
  AC0[0]=zero_p;


  AC0[1]=(multi_poly_2048(f, xb, n, q));
  
  for (int i = 0; i < N/2; ++i) {
    vector<vector<polynomial>> ans(2*(k/logb),vector<polynomial>(2,zero_p));

    for (int j1 = 0; j1 < 2 * (k / logb); ++j1) {
      for (int j2 = 0; j2 < 2; ++j2) {
        polynomial gsw= quick_mul_poly(n, q, (a_n[2*i]+a_n[2*i+1]), ek0[i][j1][j2]);
        add_poly(gsw, quick_mul_poly(n, q, (a_n[2 * i]), ek1[i][j1][j2]),q);
        add_poly(gsw, quick_mul_poly(n, q, (a_n[2 * i+1]), ek2[i][j1][j2]), q);
        add_poly(ans[j1][j2], gsw, q);
      }
    }
    vector<polynomial> ansAC0 = rlwe_multi_rgsw(n, q, k, AC0, ans, b, logb);
    for (int j = 0; j < 2; ++j) {
      add_poly(AC0[j], ansAC0[j], q);
    }
  }

  //insert index
  AC0.push_back(p_index);

  return AC0;

}







