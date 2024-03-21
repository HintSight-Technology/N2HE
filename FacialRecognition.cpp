// ================================================
// Clean up 2023 11 10
// Facial Recognition for 30 people
// ================================================

#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

// Switches
const bool ENABLE_MULTITHREADING = true;


// ================================================

// parameters for Neural Network
                     
const int nh2 = 30;                     // #Nodes in Hidden Layer
const int classes = 30;

double w1[512][nh2];                     // (DOUBLE) Weights between Input Layer and Hidden Layer

double w2[nh2][classes];                      // (DOUBLE) Weights between Output Layer and Hidden Layer

double bias1[nh2];                       // (DOUBLE) Bias in Hidden Layer

double bias2[classes];                       // (DOUBLE) Bias in Output Layer

int scaler0 = 16;

int scaler1 = 16;                       // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int scaler2 = 4;                        // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int ww1[512][nh2];                        // (int) Weights between Input Layer and Hidden Layer

int ww2[nh2][classes];                        // (int) Weights between Output Layer and Hidden Layer

int ibias1[nh2];                          // (int) Bias in Hidden Layer

int ibias2[classes];                          // (int) Bias in Output Layer

int64_t hat = -3000;                      // threshold for unknown people (max < hat)
int64_t threshold = 4000;               // threshold for unknown people (max - 2nd max)


// Test set 

vector< vector<int> > test_vectors;

vector<int> test_labels;

//parameter for NTT

const int64_t mod = 576460752213245953;

// ================================================

// Define Type polynomial: 
// polynomial p is an array of size (deg(p) + 2): [deg(p), coef_0(p), coef_1(p), coef_2(p),..., coef_{deg(p)}(p)]
typedef vector<int64_t> polynomial;

// ================================================
template<size_t N, size_t qNum, enum NTTSelector NTT, enum ImplSelector Impl>
struct Params {
  using R= Rq<N,qNum,NTT,Impl>;
  int index;
  int n;
  int k;
  int NN;
  int q1;
  int64_t b;
  int logb;
  int64_t delta;
  int64_t q2;
  vector<vector<int>> const &ct_LWE;
  vector<vector<R>> const &ek0;
  vector<vector<R>> const &ek1;
  vector<vector<R>> const &ek2;
  //polynomial const &s;
  polynomial const &f;
  
};

// Define to sort
bool cmp(vector<polynomial> a, vector<polynomial> b) {
   return a[2][1] < b[2][1];
 }


// ======================================================  Neural Network  ==========================================================================

//read test vectors and labels
void read_vector() {
  ifstream infile;
  infile.open("../FacicalRecognition_data/testImgAll.txt");
  if(!infile.is_open()){
    cout <<"Cannot open file testImgAll.txt"<<endl;
  }
   int num = 0;
   infile >> num;
  //int num = 1;
  for (int i = 0; i < num; ++i) {
    vector<int> ivector;
    for (int j = 0; j < 512; ++j) {
      double t;
      infile >> t;
      t = t * scaler0;
      int tt=(int)t;
      if (t > 0 && t - tt > 0.5) {
        ivector.push_back( tt + 1);
      }
      else if (t < 0 && tt - t>0.5) {
        ivector.push_back( tt - 1);
      }
      else {
        ivector.push_back( tt);
      }
    }
    test_vectors.push_back(ivector);
  }
  infile.close();
}


// Import trained Neural Network
// Turn it into int Neural Network

void read_w() {

  ifstream infile;
  infile.open("../FacicalRecognition_data/cnn30.txt");
  if(!infile.is_open()){
    cout <<"Cannot open file cnn30.txt"<<endl;
  }
  double round_w = 1.0;
  for (int i = 0; i < nh2; i++)
  {
    for (int j = 0; j < 512; j++)
    {
      infile >> w1[j][i];
      double t = w1[j][i] * scaler1;
      ww1[j][i] = ((int)(w1[j][i] * scaler1));
      if (t - ww1[j][i] > 0.5 && t > 0) {
        ww1[j][i]++;
      }
      else if (t < 0 && w1[j][i] - t>0.5) {
        ww1[j][i]--;
      }
    }
  }
  for (int j = 0; j < classes; j++)
  {
    for (int i = 0; i < nh2; i++)
    {
      infile >> w2[i][j];
      double t = w2[i][j] * scaler2;
      int tt = (int)(w2[i][j] * scaler2);
      if (t > 0 && t - tt > 0.5) {
        ww2[i][j] = tt + 1;
      }
      else if (t < 0 && tt - t>0.5) {
        ww2[i][j] = tt - 1;
      }
      else {
        ww2[i][j] = tt;
      }
    }
  }
  for (int i = 0; i < nh2; i++)
  {
    infile >> bias1[i];
    double t = bias1[i] * scaler1;
    int tt = (int)(bias1[i] * scaler1);
    if (t > 0 && t - tt > 0.5) {
      ibias1[i] = tt + 1;
    }
    else if (t < 0 && tt - t>0.5) {
      ibias1[i] = tt - 1;
    }
    else {
      ibias1[i] = tt;
    }
  }
  for (int i = 0; i < classes; i++)
  {
    infile >> bias2[i];
    double t = bias2[i] * scaler2;
    int tt = (int)(bias2[i] * scaler2);
    if (t > 0 && t - tt > 0.5) {
      ibias2[i] = tt + 1;
    }
    else if (t < 0 && tt - t>0.5) {
      ibias2[i] = tt - 1;
    }
    else {
      ibias2[i] = tt;
    }
  }
  infile.close();
}

// label index
string print_name(int index) {
  return "name"+to_string(index);
}

// END OF Import trained Neural Network

// Look up table function
template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> IP1_LUT(int index, int q1, int64_t q2, int n, int64_t delta, int k, int N, vector<vector<int>> const &ct_LWE, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek0, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek1, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek2, polynomial const& f, int64_t b, int logb) {
  using R=Rq<NN, qNum, NTT, Impl>;
  vector<int> ct_ip1;
    for (register int i = 0; i < 512; ++i) {
        vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, 512, ct_LWE[i], ww1[i][index]);
        if (i == 0) {
          ct_ip1=ct_iip0;
        }
        else {
          ct_ip1 = LWE32_Add_ct(q1, 512, ct_iip0, ct_ip1);
        }
      }

      ct_ip1 = LWE32_Plain_Add_ct(q1, 512, ct_ip1, 2*scaler0*ibias1[index]);
      return LUT_2048<2048,1,FNTT,AVX2>(index,q1, q2, n, delta, k, N,ct_ip1, _ek0, _ek1,_ek2, f,b,logb);

}


int main(){

  cout <<"Task: Privacy-preserving facial recognition. "<<endl;
  int start_index = 0;
  int test_count = 37;
  
  // ======================== Initialization and Information Output ==========================
  read_vector();

  read_w();


  srand(time(0));
  clock_t start, end;

  //set n,q,t for LWE scheme
  cout << "Parameters for LWE scheme:" << endl;
  const int n1 = 512;                        //vector length
  const int q1 = 1<<13;                      //ciphertext modulus
  const int var1 = 3;                          //variance of error
  
  cout << "vector length = " << n1 <<", ciphertext modulus = " << q1 <<", variance of error = 2^ -"<<var1 << endl;
  cout << "----------" << endl;

  
  //set n,q,k for RLWE scheme
  cout << "Parameters for RLWE scheme:" << endl;
  const int n2 = 1<<11; //vector length=exp of polynomial=2048
  const int n3 = n2;
  const int var2=3;                      //variance of error
  const int k2 = 60; //k=log(q)
  const int64_t b = 1<<30;
  const int logb = 30; //logb=log2(b)-1
  const int64_t q2 = mod; //ciphertext modulus
  const int64_t delta = 1099511627776;       //delta * ip1 < q/4 && delta*ip2 < q/2
  cout << "polynomial degree = " << n2 << ", ciphertext modulus: " << q2 <<", variance of error = 2^ -"<<var2  <<", external product based on " << logb << " bits extension" << endl;
  cout << "----------" << endl;

  using R = Rq<n2, 1, FNTT, AVX2>;
  int *x = new int [n1];
  x=LWE32_KeyGen(n1);
  //generate LWE sk and read test images



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
    for(int j=0;j<4;++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=2048;++h){
          if(ek0i[j][k][h]<0){
            tempp(0,h-1)=ek0i[j][k][h]+q2;
          }
          else{
            tempp(0,h-1)=ek0i[j][k][h];
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
    for(int j=0;j<4;++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=2048;++h){
          if(ek1i[j][k][h]<0){
            tempp(0,h-1)=ek1i[j][k][h]+q2;
          }
          else{
            tempp(0,h-1)=ek1i[j][k][h];
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
    vector<vector<polynomial>> ek2i = key_encrypt_3(n2, q2, k2,var2, x[2 * i], x[2 * i + 1], s, b, logb);
    for(int j=0;j<4;++j){
      for(int k=0;k<2;++k){
        R tempp;
        for(int h=1;h<=2048;++h){
          if(ek2i[j][k][h]<0){
            tempp(0,h-1)=ek2i[j][k][h]+q2;
          }
          else{
            tempp(0,h-1)=ek2i[j][k][h];
          }
        }
        tempp.ntt();
        tempek.push_back(tempp);
      }
    }
    _ek2.push_back(tempek);
  }
  
 // cout<<"Evaluation key generation done. Store in RNS form."<<endl;
 // cout << "----------" << endl;
  double ave_time = 0;
  struct timeval tstart, tend;

  // =========================== Start Testing ===============================================

  cout <<"Test start from photo No. "<<start_index <<", to No. "<< start_index+test_count<<endl;
  for (int tt = start_index; tt < start_index + test_count; ++tt) {
    cout <<"Test no: "<<tt <<": ";
    //encrypt image i, scale from {0,1} to {0,2}

    vector<vector<int>> ct_i;
    for (int j = 0; j < 512; ++j) {
      vector<int>ctj = LWE32_Enc(q1, n1,var1, 2*test_vectors[tt][j], x);
      ct_i.push_back(ctj);
    }
    int true_label = 0;

    //compute the first inner-product
    //LUT, construct polynomial f
    polynomial f = p_f(n2, delta, q1);
    
    vector<vector<int64_t>>ct_lut;
    double time2 = 0;


    if (ENABLE_MULTITHREADING == 1){
      //multithreading

      //construct inputs
      vector<Params<n2, 1, FNTT, AVX2>> inputs;
      for (int i = 0; i < nh2; ++i) {
        Params<n2, 1, FNTT, AVX2> a = { i,  n2,  k2, n1, q1, b,logb, delta, q2, ct_i, _ek0, _ek1, _ek2, f };
        inputs.push_back(a);
      }

      condition_variable cv;
      mutex mu;
      vector<vector<polynomial>> results;
      vector<thread> w;

      gettimeofday(&tstart,NULL);
      for (auto& input : inputs) {
        w.emplace_back([&]() {
          auto r = IP1_LUT<2048,1,FNTT,AVX2>(input.index, input.q1, input.q2, input.n, input.delta, input.k, input.NN, input.ct_LWE, input.ek0, input.ek1,input.ek2, input.f, input.b, input.logb);
          mu.lock();
          results.push_back(std::move(r));
          mu.unlock();
          cv.notify_one();
          });
      }

      {
        unique_lock<mutex> lk(mu);
        cv.wait(lk, [&] { return results.size() >= nh2; });
      }

      for (auto& th : w) {
        th.join();
      }
      gettimeofday(&tend,NULL);

      // sort hidden layer output in order
      sort(results.begin(), results.end(), cmp);
  
      for (int i = 0; i < nh2; ++i) {
        results[i].pop_back();
        vector<int64_t> ans2 = Extract0(q2, delta, n2, results[i]);
        ct_lut.push_back(ans2);
      }
      
      time2 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;      // time cost of LUT (multithreading)
    }

    else {
      // single 
      vector<vector<int>> ct_ip1;
      gettimeofday(&tstart,NULL);

      // First inner product
      for (register int i = 0; i < 512; ++i) {
        for (register int j = 0; j < nh2; ++j) {
          vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, n1, ct_i[i], ww1[i][j]);
          if (i == 0) {
            ct_ip1.push_back(ct_iip0);
          }
          else {
            ct_ip1[j] = LWE32_Add_ct(q1, n1, ct_iip0, ct_ip1[j]);
          }
        }
      }

      // First fc layer, add bias
      for (int i = 0; i < nh2; ++i) {
        ct_ip1[i] = LWE32_Plain_Add_ct(q1, n1, ct_ip1[i], 2*scaler0*ibias1[i]);
      }

      gettimeofday(&tend,NULL);
      double time1 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;     // time cost of the first inner-product and bias
              
      gettimeofday(&tstart,NULL);
      // LUT in hidden layer
      for (int i = 0; i < nh2; ++i) {
        vector<polynomial> relu=LUT_2048<2048,1,FNTT,AVX2>(i,q1, q2, n2, delta, k2, n1,ct_ip1[i], _ek0, _ek1,_ek2,f,b,logb);
        relu.pop_back();
        vector<int64_t> ans2 = Extract0(q2,delta, n2, relu);
        ct_lut.push_back(ans2);
      }

      gettimeofday(&tend,NULL);

      time2 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;      // time cost of LUT (single thread)
      time2 += time1;      
    }

    //compute the second inner-product
    vector<vector<int64_t>> ct_ip2;

    gettimeofday(&tstart,NULL);
    for (int i = 0; i < classes; ++i) {
      vector<int64_t> ct_iip0 = LWE64_Plain_Multi_ct(q2, n3, ct_lut[0], ww2[0][i]);
      for (int j = 1; j < nh2; ++j) {
        vector<int64_t> ct_iipj = LWE64_Plain_Multi_ct(q2, n3, ct_lut[j], ww2[j][i]);
        ct_iip0 = LWE64_Add_ct(q2, n3, ct_iip0, ct_iipj);
      }
      ct_iip0 = LWE64_Plain_Add_ct(q2, n3, ct_iip0, ibias2[i] * (2*delta*scaler0)); // add bias
      ct_ip2.push_back(ct_iip0); 
    }

    

    //decrypt
    int64_t max = -1*mod;
    int64_t max2 = max;
    int cindex = 0;
    int cindex2 = 0;
    double likelihood[classes];

    //check vector
    for (int i = 0; i < classes; ++i) {
      int64_t tempr = LWE64_Dec(q2, n3, ct_ip2[i], lwe_s);
      likelihood[i] = 1.0 * tempr / delta;
      if (likelihood[i] > max) {
        max = (int64_t) likelihood[i];
        cindex = i;
      }
    }
    
    for (int i = 0; i < classes; ++i) {
      if (likelihood[i] > max2 && i != cindex) {
        max2 = (int64_t) likelihood[i];
        cindex2 = i;
      }
    }

    gettimeofday(&tend,NULL);
    double time3 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;

    double total_time = time2+time3;
    ave_time += total_time;

    // output 
     if (max < hat && max-max2<threshold) {
      cout <<  "Prediction: This people is not in dataset, time: " << total_time << "s. " << endl;
     } 
     else {
      cout <<  "Prediction: " << print_name(cindex) << ", time: " << total_time << "s. " << endl;
     }
  }
  
  // =========================== End of Testing ==============================================
  cout << "----------" << endl;
  cout <<"Average time: "<<ave_time/test_count<<"s. "<<endl;
  return 0;

}



