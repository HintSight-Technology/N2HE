// ================================================
// Clean up 2023 11 10
// MNIST, BP network
// ================================================

#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

// Switches
const bool ENABLE_MULTITHREADING = true;


// ================================================

// parameters for Neural Network

const int nh = 30;                      // #Nodes in Hidden Layer
const int classes = 10;

double w1[784][nh];                     // (DOUBLE) Weights between Input Layer and Hidden Layer

double w2[nh][classes];                      // (DOUBLE) Weights between Output Layer and Hidden Layer

double bias1[nh];                       // (DOUBLE) Bias in Hidden Layer

double bias2[classes];                       // (DOUBLE) Bias in Output Layer

int scaler1 = 16;                       // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int scaler2 = 8;                        // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int ww1[784][nh];                       // (int) Weights between Input Layer and Hidden Layer

int64_t ww2[nh][classes];                        // (int) Weights between Output Layer and Hidden Layer

int ibias1[nh];                         // (int) Bias in Hidden Layer

int64_t ibias2[classes];                         // (int) Bias in Output Layer


// Test set 

vector<int> label;     

vector< vector<int> > test_images;

vector<int> test_labels;

int testNumber = 10000;

//test image index to index+count

int start_index = 0;

int test_count = 100;

//parameter for NTT

const int64_t mod = 576460752213245953;     //2^59
const int64_t delta = 1099511627776;  //2^40

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
  polynomial const &f;
};

// Define to sort
bool cmp(vector<polynomial> a, vector<polynomial> b) {
   return a[2][1] < b[2][1];
 }


// ======================================================  Neural Network  ==========================================================================

// Import MNIST dataset 

int reverse_int(int i) {
  unsigned char ch1, ch2, ch3, ch4;
  ch1 = i & 255;
  ch2 = (i >> 8) & 255;
  ch3 = (i >> 16) & 255;
  ch4 = (i >> 24) & 255;
  return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}

vector<vector<int>> read_test_images() {
  ifstream file("../MNIST_data/t10k-images-idx3-ubyte", ios::binary);
  if (file.is_open()) {
    int magic_number = 0;
    int number_of_images = 0;
    int row = 0;
    int col = 0;

    file.read((char*)&magic_number, sizeof(magic_number));
    file.read((char*)&number_of_images, sizeof(number_of_images));
    file.read((char*)&row, sizeof(row));
    file.read((char*)&col, sizeof(col));

    magic_number = reverse_int(magic_number);
    number_of_images = reverse_int(number_of_images);
    row = reverse_int(row);
    col = reverse_int(col);
    number_of_images = testNumber;
    for (int i = 0; i < number_of_images; i++) {
      vector<int> this_image;
      for (int r = 0; r < row; r++) {
        for (int c = 0; c < col; c++) {
          unsigned char pixel = 0;
          file.read((char*)&pixel, sizeof(pixel));
          this_image.push_back(pixel != 0 ? 1 : 0);  // Turn the image to binary
        }
      }
      test_images.push_back(this_image);
    }
  }
  else{
    cout <<"cannot open file t10k-images-idx3-ubyte. "<<endl;
  }
  return test_images;
}

vector<int> read_test_labels() {
  ifstream file;
  file.open("../MNIST_data/t10k-labels-idx1-ubyte", ios::binary);
  if (file.is_open()) {
    int magic_number = 0;
    int number_of_images = 0;

    file.read((char*)&magic_number, sizeof(magic_number));
    file.read((char*)&number_of_images, sizeof(number_of_images));

    magic_number = reverse_int(magic_number);
    number_of_images = reverse_int(number_of_images);
    number_of_images = testNumber;
    for (int i = 0; i < number_of_images; i++) {
      unsigned char label = 0;
      file.read((char*)&label, sizeof(label));
      test_labels.push_back((int)label);
    }
  }
  else{
    cout <<"cannot open file t10k-labels-idx1-ubyte. "<<endl;
  }
  return test_labels;
}

// END OF Import MNIST dataset 

// ================================================

// Import trained Neural Network
// Turn it into int Neural Network

void read_w() {

  ifstream infile;
  infile.open("../MNIST_data/bp.txt");
  if(!infile.is_open()){
    cout <<"cannot open file bp.txt. "<<endl;
  }
  int nh;
  double round_w = 1.0;
  infile >> nh;
  for (int i = 0; i < 784; i++)
  {
    for (int j = 0; j < nh; j++)
    {
      infile >> w1[i][j];
      double t = w1[i][j] * scaler1;
      ww1[i][j] = ((int)(w1[i][j] * scaler1));
      if (t - ww1[i][j] > 0.5 && t > 0) {
        w1[i][j]++;
      }
      else if (t < 0 && w1[i][j] - t>0.5) {
        w1[i][j]--;
      }
    }
  }for (int i = 0; i < nh; i++)
  {
    for (int j = 0; j < 10; j++)
    {
      infile >> w2[i][j];
      double t = w2[i][j] * scaler2;
      int64_t tt = (int64_t)(w2[i][j] * scaler2);
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
  for (int i = 0; i < nh; i++)
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
  for (int i = 0; i < 10; i++)
  {
    infile >> bias2[i];
    double t = bias2[i] * scaler2;
    int64_t tt = (int64_t)(bias2[i] * scaler2);
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
// END OF Import trained Neural Network

// plaintext test accuracy
void test_pt(){
  int cnt = 0;
  int largeip1 = 0;
  int largeip2 = 0;

  for(int i=start_index;i<start_index+test_count;++i){
    int true_label = test_labels[i];
    int64_t relu[nh];
    for(int j = 0;j<nh;++j){
      int ip = 0;
      for(int k = 0;k<784;++k){
        ip += (2*ww1[k][j]*test_images[i][k]);
        ip = modq_32(ip,8192);
      }
      ip += (2*ibias1[j]);
      if(ip > 2048 || ip < -2048){
        largeip1 ++;
      }
      ip = modq_32(ip,8192);
      if(ip < 0){
        ip = 0;
      }
      
      relu[j]=(int64_t)ip*delta;
    }
    int64_t output[10];
    for(int j = 0; j<classes;++j){
      int64_t ip2 = 0;
      for(int k = 0 ; k < nh ; ++k){
        ip2 += (2*ww2[k][j]*relu[k]);
        ip2 = modq_64(ip2,mod);
      }
      ip2 += (2*scaler1*delta*ibias2[j]);
      if(ip2 > mod/2 || ip2 < -mod/2){
        largeip2++;
      }
      ip2 = modq_64(ip2,mod);
      
      output[j]=ip2;
    }

    int64_t max = -1*mod;
    int cindex = 0;

    for (int i = 0; i < classes; ++i) {
      if (output[i] > max) {
        max = output[i];
        cindex = i;
      }
    }
    
  if (cindex == true_label) {
      // when correct
      cnt++;
    }
  }
  cout <<"Test in plaintext, Correctness = "<<((double)cnt / test_count) * 100 << "%" << endl;
}


// LUT functions
template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> IP1_LUT(int index, int q1, int64_t q2, int n, int64_t delta, int k, int N, vector<vector<int>> const &ct_LWE, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek0, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek1, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek2,  polynomial const& f, int64_t b, int logb) {
  using R=Rq<NN, qNum, NTT, Impl>;
  vector<int> ct_ip1;
    for (register int i = 0; i < 784; ++i) {
        vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i], ww1[i][index]);
        if (i == 0) {
          ct_ip1=ct_iip0;
        }
        else {
          ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
        }
      }
      ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*ibias1[index]);
      return LUT_2048<2048,1,FNTT,AVX2>(index,q1,q2, n, delta, k, N,ct_ip1, _ek0, _ek1,_ek2, f,b,logb);
}

// LUT functions
template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> IP1_LUT_slow(int index, int q1, int64_t q2, int n, int64_t delta, int k, int N, vector<vector<int>> const &ct_LWE, const vector<vector<vector<polynomial>>> &_ek0, const vector<vector<vector<polynomial>>> &_ek1, const vector<vector<vector<polynomial>>> &_ek2, polynomial const& f, int64_t b, int logb) {
  using R=Rq<NN, qNum, NTT, Impl>;
  vector<int> ct_ip1;
    for (register int i = 0; i < 784; ++i) {
        vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i], ww1[i][index]);
        if (i == 0) {
          ct_ip1=ct_iip0;
        }
        else {
          ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
        }
      }
      ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*ibias1[index]);
      return LUT_2048_slow<2048,1,FNTT,AVX2>(index,q1,q2, n, delta, k, N,ct_ip1, _ek0, _ek1,_ek2, f,b,logb);
}


int main(){
  cout <<"Task: Privacy-preserving neural network inference on MNIST dataset."<<endl;

 // ======================== Initialization and Information Output ==========================
  label=read_test_labels();   
  vector<vector<int>> image = read_test_images();             
  read_w();

  for(int test_time = 0; test_time < 1 ; ++test_time){
    cout << "Parameters for image: " << endl;
    cout << "scaler1 = " << scaler1 << " , scaler2 = " << scaler2 << ", image size: {0,2}^784" << endl;
    cout << "----------" << endl;
    srand(time(0));
    clock_t start, end;

    test_pt(); //test plaintext accuracy
    cout << "----------" << endl;

    //set n,q,t for LWE scheme
    cout << "Parameters for LWE scheme:" << endl;
    const int n1 = 512;                        //vector length
    const int q1 = 1<<13;                      //ciphertext modulus
    const int t1 = 1<<11;                      //plaintext modulus
    const int var1 = 3;

    cout << "vector length = " << n1 <<", ciphertext modulus = " << q1 <<", variance of error = 2^ -"<<var1 << endl;
    cout << "----------" << endl;

    

    //set n,q,k for RLWE scheme
    cout << "Parameters for RLWE scheme:" << endl;
    const int n2 = 2048; //vector length=exp of polynomial=2048
    const int n3 = n2;
    const int k2 = 60; //k=log(q)
    const int var2 = 3;
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
      for(int j=0;j<4;++j){
        for(int k=0;k<2;++k){
          R tempp;
          for(int h=1;h<=2048;++h){
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
      for(int j=0;j<4;++j){
        for(int k=0;k<2;++k){
          R tempp;
          for(int h=1;h<=2048;++h){
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
      vector<vector<polynomial>> ek2i = key_encrypt_3(n2, q2, k2,var2, x[2 * i], x[2 * i + 1], s, b, logb);
      for(int j=0;j<4;++j){
        for(int k=0;k<2;++k){
          R tempp;
          for(int h=1;h<=2048;++h){
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
    

  // =========================== Start Testing ===============================================
    
    double ave_time = 0;
    double ave_ip1_lut_time = 0;
    double ave_ip2_time = 0;
    struct timeval tstart, tend;
    int correct_num=0;
    int correct_num_test = 0;

    cout <<"Test start from image No. "<<start_index <<", to No. "<< start_index+test_count-1<<endl;

    for (int tt = start_index; tt < start_index + test_count; ++tt) {
      //encrypt image i, scale from {0,1} to {0,2}

      vector<vector<int>> ct_i;
      for (int j = 0; j < 784; ++j) {
        vector<int>ctj = LWE32_Enc(q1, n1,var1, 2*image[tt][j], x);
        ct_i.push_back(ctj);
      }
      int true_label = label[tt];

      //compute the first inner-product
      //lut, construct polynomial f
      polynomial f = p_f(n2, delta, q1);
      
      int dec_ip1[nh];
      vector<vector<int64_t>>ct_lut;
      double time2 = 0;


      if (ENABLE_MULTITHREADING == 1){
        //multithreading
        gettimeofday(&tstart,NULL);

        //construct inputs
        vector<Params<n2, 1, FNTT, AVX2>> inputs;
        for (int i = 0; i < nh; ++i) {
          Params<n2, 1, FNTT, AVX2> a = { i,  n2,  k2, n1, q1, b,logb, delta, q2, ct_i, _ek0, _ek1, _ek2, f };
          inputs.push_back(a);
        }

        condition_variable cv;
        mutex mu;
        vector<vector<polynomial>> results;
        vector<thread> w;

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
          cv.wait(lk, [&] { return results.size() >= nh; });
        }

        for (auto& th : w) {
          th.join();
        }
        
        sort(results.begin(), results.end(), cmp);
    
        for (int i = 0; i < nh; ++i) {
          results[i].pop_back();
          vector<int64_t> ans2 = Extract0(q2, delta, n2, results[i]);
          ct_lut.push_back(ans2);
        }

        gettimeofday(&tend,NULL);
        time2 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;      // time cost of LUT (multithreading)
      }

      else {
        // single 
        vector<vector<int>> ct_ip1;
        gettimeofday(&tstart,NULL);

        for (register int i = 0; i < 784; ++i) {
          for (register int j = 0; j < nh; ++j) {
            vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, n1, ct_i[i], ww1[i][j]);
            if (i == 0) {
              ct_ip1.push_back(ct_iip0);
            }
            else {
              ct_ip1[j] = LWE32_Add_ct(q1, n1, ct_iip0, ct_ip1[j]);
            }
          }
        }
      
        for (int i = 0; i < nh; ++i) {
          ct_ip1[i] = LWE32_Plain_Add_ct(q1, n1, ct_ip1[i], 2*ibias1[i]);
        }

        gettimeofday(&tend,NULL);
        double time1 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;     // time cost of the first inner-product
            
        gettimeofday(&tstart,NULL);
        for (int i = 0; i < nh; ++i) {
          vector<polynomial> relu=LUT_2048<2048,1,FNTT,AVX2>(i,q1,q2, n2, delta, k2, n1,ct_ip1[i], _ek0, _ek1,_ek2, f,b,logb);
          relu.pop_back();
          vector<int64_t> ans2 = Extract0(q2,delta, n2, relu);
          ct_lut.push_back(ans2);
        }
        gettimeofday(&tend,NULL);
        time2 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;      // time cost of LUT (singlethreading)
        time2 += time1;      
      }

      //compute the second inner-product
      vector<vector<int64_t>> ct_ip2;

      gettimeofday(&tstart,NULL);
      for (int i = 0; i < classes; ++i) {
        vector<int64_t> ct_iip0 = LWE64_Plain_Multi_ct(q2, n2, ct_lut[0], ww2[0][i]);
        for (int j = 1; j < nh; ++j) {
          vector<int64_t> ct_iipj = LWE64_Plain_Multi_ct(q2, n2, ct_lut[j], ww2[j][i]);
          ct_iip0 = LWE64_Add_ct(q2, n2, ct_iip0, ct_iipj);
        }
        ct_iip0 = LWE64_Plain_Add_ct(q2, n2, ct_iip0, ibias2[i]*delta );
        ct_ip2.push_back(ct_iip0);
      }

      //decrypt
      int64_t max = -1*mod;
      int cindex = 0;

      //check vector
      for (int i = 0; i < classes; ++i) {
        int64_t tempr = LWE64_Dec(q2, n2, ct_ip2[i], lwe_s);
        if (tempr > max) {
          max = tempr;
          cindex = i;
        }
      }
      
      if (cindex == true_label) {
        // when correct
        correct_num++;
      }
      gettimeofday(&tend,NULL);
      double time3 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;

      // sum time
      double total_time = time2+time3;
      ave_time += total_time;
      ave_ip1_lut_time += time2;
      ave_ip2_time += time3;
    }
    
    // =========================== End of Testing ==============================================
    cout << "Correctness is: " << ((double)correct_num / test_count) * 100 << "%" << endl;
    cout <<"Average time: "<<ave_time/test_count<<"s. "<<endl;
    cout <<"Average time for ip1+lut: "<<ave_ip1_lut_time/test_count<<"s. "<<endl;
    cout <<"Average time for ip2: "<<ave_ip2_time/test_count<<"s. "<<endl;
  }
  return 0;
}



