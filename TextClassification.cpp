// ================================================
// Clean up 2023 11 10
// Test binary classification, TextCNN + linear
// ================================================

#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

// ================================================================================================================================================


// Switches
const bool ENABLE_MULTITHREADING = false;


// ================================================

// parameters for Neural Network
                     
const int nh2 = 2;											// #Nodes in Hidden Layer

const int feature_dim = 300;

double w1[feature_dim][nh2];                     // (DOUBLE) Weights between Input Layer and Hidden Layer

double bias1[nh2];                       // (DOUBLE) Bias in Hidden Layer

int scaler0 = 30;

int scaler1 = 6;                       // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int ww1[feature_dim][nh2];												// (int) Weights between Input Layer and Hidden Layer


int ibias1[nh2];											  	// (int) Bias in Hidden Layer

// Test set    

vector< vector<int> > test_vectors;

vector<int> test_labels;


// ======================================================  Neural Network  ==========================================================================

//read test vectors and labels
void read_vector() {
	ifstream infile;
	infile.open("../TextClassification_data/test.txt");
	if(!infile.is_open()){
    	cout <<"Cannot open file test.txt"<<endl;
  	}
	int num = 100;
	for (int i = 0; i < num; ++i) {
		vector<int> ivector;
		for (int j = 0; j < feature_dim; ++j) {
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
	infile.open("../TextClassification_data/cnn-text.txt"); // bias (1), weight (feature_dim)
	if(!infile.is_open()){
    	cout <<"Cannot open file cnn-text.txt"<<endl;
  	}
	double round_w = 1.0;
	for (int i = 0; i < nh2; i++)
	{
		for (int j = 0; j < feature_dim; j++)
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

	infile.close();
}


// END OF Import trained Neural Network

int main(){

	cout <<"Task: Privacy-preserving text classification."<<endl;
	int start_index = 0;
	int test_count = 100;

	// ======================== Initialization and Information Output ==========================
	read_vector();

	read_w();


	srand(time(0));
	clock_t start, end;

	//set n,q,t for LWE scheme
	 cout << "Parameter setting for LWE scheme:" << endl;
	const int n1 = 512;                        //vector length
	const int q1 = 1<<13;                      //ciphertext modulus
	const int var1 = 3;                          //variance of error
  
 	cout << "vector length = " << n1 <<", ciphertext modulus = " << q1 <<", variance of error = 2^ -"<<var1 << endl;
 	cout << "----------" << endl;

	int *x = new int [n1];
  	x=LWE32_KeyGen(n1);

	double total_time = 0;
	struct timeval tstart, tend;


	// =========================== Start Testing ===============================================
	cout <<"Test start from text No. "<<start_index <<", to No. "<< start_index+test_count<<endl;
	int err = 0;
	for (int tt = start_index; tt < start_index + test_count; ++tt) {

		//encrypt image i, scale from {0, 1} to {0, 2}

		vector<vector<int>> ct_i;
		for (int j = 0; j < feature_dim; ++j) {
			vector<int>ctj = LWE32_Enc(q1, n1, var1, 2*test_vectors[tt][j], x);
			ct_i.push_back(ctj);
		}
		int true_label = 0;

		//compute the first inner-product
		vector<vector<int>> ct_ip1;
		gettimeofday(&tstart,NULL);

		for ( int i = 0; i < feature_dim; ++i) {
			for ( int j = 0; j < nh2; ++j) {
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

		//decrypt 1-length output vector
		int dec1 =  LWE32_Dec(q1, n1, ct_ip1[1], x);
		int dec0 =  LWE32_Dec(q1, n1, ct_ip1[0], x);

		// check threshold
		int ans;
		int alg;
		if (tt<50) {
			ans = 0;
		}
		else{
			ans = 1;
		}
		if (dec0 > dec1) {
			alg = 0;
		}
		else{
			alg = 1;
		}
		err += (ans!=alg);
		gettimeofday(&tend,NULL);
		double time1 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;     // time cost of the first inner-product
		total_time += time1;
	}
	
	// =========================== End of Testing ==============================================
	cout << "Correctness: " << ((double)(test_count-err)/(double)test_count)*100 <<"%" << endl;
	cout << "Average Time: " << total_time/test_count << endl;
	return 0;
}


