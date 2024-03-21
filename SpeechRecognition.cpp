// ================================================
// Clean up 2023 11 10
// speaker verification based on MFCC, BP network
// ================================================

#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

// Switches
const bool ENABLE_MULTITHREADING = true;


// ================================================

// parameters for Neural Network
                    
const int nh2 = 30;											// #Nodes in Hidden Layer
const int classes = 34;

double w1[40][nh2];                     // (DOUBLE) Weights between Input Layer and Hidden Layer

double w2[nh2][classes];                      // (DOUBLE) Weights between Output Layer and Hidden Layer

double bias1[nh2];                       // (DOUBLE) Bias in Hidden Layer

double bias2[classes];                       // (DOUBLE) Bias in Output Layer

int scaler0 = 8;

int scaler1 = 8;                       // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int scaler2 = 4;                        // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int ww1[40][nh2];												// (int) Weights between Input Layer and Hidden Layer

int ww2[nh2][classes];												// (int) Weights between Output Layer and Hidden Layer

int ibias1[nh2];											  	// (int) Bias in Hidden Layer

int ibias2[classes];									  			// (int) Bias in Output Layer

// ================================================

// Test set 
   

vector< vector<int> > test_vectors;

vector<int> test_labels;


// ================================================
//parameter for NTT

const int64_t mod = 576460752213245953;
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

//read test vectors and labels
void read_vector() {
	ifstream infile;
	infile.open("../SpeechRecognition_data/speech-all.txt");
	if(!infile.is_open()){
    cout <<"Cannot open file speech-all.txt"<<endl;
  }
	 int num = 1700;
	for (int i = 0; i < num; ++i) {
		vector<int> ivector;
		for (int j = 0; j < 40; ++j) {
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
			// cout <<endl;
		}
		test_vectors.push_back(ivector);
	}
	infile.close();
}


// Import trained Neural Network
// Turn it into int Neural Network

void read_w() {

	ifstream infile;
	infile.open("../SpeechRecognition_data/speech-30.txt");
	if(!infile.is_open()){
    cout <<"Cannot open file speech-30.txt"<<endl;
  }
	double round_w = 1.0;
	for (int i = 0; i < 40; i++)
	{
		for (int j = 0; j < nh2; j++)
		{
			infile >> w1[i][j];
			double t = w1[i][j] * scaler1;
			ww1[i][j] = ((int)(w1[i][j] * scaler1));
			if (t - ww1[i][j] > 0.5 && t > 0) {
				ww1[i][j]++;
			}
			else if (t < 0 && w1[i][j] - t>0.5) {
				ww1[i][j]--;
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
	for (int j = 0; j < nh2; j++)
	{
		for (int i = 0; i < classes; i++)
		{
			infile >> w2[j][i];
			double t = w2[j][i] * scaler2;
			int tt = (int)(w2[j][i] * scaler2);
			if (t > 0 && t - tt > 0.5) {
				ww2[j][i] = tt + 1;
			}
			else if (t < 0 && tt - t>0.5) {
				ww2[j][i] = tt - 1;
			}
			else {
				ww2[j][i] = tt;
			}
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


// END OF Import trained Neural Network

template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
vector<polynomial> IP1_LUT(int index, int q1, int64_t q2, int n, int64_t delta, int k, int N, vector<vector<int>> const &ct_LWE, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek0, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek1, const vector<vector<Rq<NN,qNum,NTT,Impl>>> &_ek2, polynomial const& f, int64_t b, int logb) {
  using R=Rq<NN, qNum, NTT, Impl>;
  vector<int> ct_ip1;
    for (register int i = 0; i < 40; ++i) {
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

	cout <<"Task: Privacy-preserving speech identity recognition."<<endl;

	int start_index = 0;
	int test_count = 1700;
	
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
	
	// ======================== Initialization and Information Output ==========================



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
		vector<vector<polynomial>> ek0i= key_encrypt_1(n2, q2, k2,var2, x[2*i],x[2*i+1], s, b, logb);
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
	

	double ave_time = 0;
	struct timeval tstart, tend;

	// =========================== Start Testing ===============================================
	cout <<"Test start from audio No. "<<0 <<", to No. "<< 33<<endl;
	for(int ttt=0;ttt<34;++ttt){
		cout <<"Test for No." <<ttt << endl;
		int result_count[34]={0};
		int result_count2[5]={0};
		ave_time=0;
		for(int tt1=0;tt1<5;++tt1){
			for (int tt = 0; tt < 10; ++tt) {

				//encrypt image i, scale from {0,1} to {0,2}

				vector<vector<int>> ct_i;
				for (int j = 0; j < 40; ++j) {
					vector<int>ctj = LWE32_Enc(q1, n1, var1, 2*test_vectors[ttt*50+tt1*10+tt][j], x);
					ct_i.push_back(ctj);
				}
				int true_label = 0;

				//compute the first inner-product
				//lut, construct polynomial f
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

					for (register int i = 0; i < 40; ++i) {
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
					for (int i = 0; i < nh2; ++i) {
						ct_ip1[i] = LWE32_Plain_Add_ct(q1, n1, ct_ip1[i], 2*scaler0*ibias1[i]);
					}

					gettimeofday(&tend,NULL);
					double time1 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;     // time cost of the first inner-product
					
					gettimeofday(&tstart,NULL);
					for (int i = 0; i < nh2; ++i) {
						vector<polynomial> relu=LUT_2048<2048,1,FNTT,AVX2>(i,q1, q2, n2, delta, k2, n1,ct_ip1[i], _ek0, _ek1,_ek2,f,b,logb);
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
					vector<int64_t> ct_iip0 = LWE64_Plain_Multi_ct(q2, n3, ct_lut[0], ww2[0][i]);
					for (int j = 1; j < nh2; ++j) {
						vector<int64_t> ct_iipj = LWE64_Plain_Multi_ct(q2, n3, ct_lut[j], ww2[j][i]);
						ct_iip0 = LWE64_Add_ct(q2, n3, ct_iip0, ct_iipj);
					}
					ct_iip0 = LWE64_Plain_Add_ct(q2, n3, ct_iip0, ibias2[i] * (2*delta*scaler0));
					ct_ip2.push_back(ct_iip0);
				}


				//decrypt
				int64_t max = -1*mod;
				int cindex = 0;
				for (int i = 0; i < classes; ++i) {
					int64_t tempr = LWE64_Dec(q2, n3, ct_ip2[i], lwe_s);
					if (tempr > max) {
						max = tempr;
						cindex = i;
					}
				}
				result_count[cindex]++;
				gettimeofday(&tend,NULL);
				double time3 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;
				double total_time = time2+time3;
				ave_time += total_time;
			}
			int max_result=-1;
			int max_index=-1;
			for(int i=0;i<34;++i){
				if(result_count[i]>max_result){
					max_result=result_count[i];
					max_index=i;
				}
			}
			result_count2[tt1]=max_index;
		}

		cout<<"Prediction result over 5*10 tests: "<<result_count2[0]<<" "<<result_count2[1]<<" "<<result_count2[2]<<" "<<result_count2[3]<<" "<<result_count2[4]<<endl;
		//cout <<"Final result: "<<final_index<<endl;
		cout <<"Total time: "<<ave_time<<"s. "<<endl;
	}
	// =========================== End of Testing ==============================================
	return 0;
}


