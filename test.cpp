
#include <iostream>
#include "include.hpp"

using namespace fntt;
using namespace std;

int main(){
  //random number generator test
  cout << "====================================================="<<endl;
  RNG();
  cout <<endl;


  //LWE test
  cout << "====================================================="<<endl;
  cout <<"Tests of LWE encryption scheme: "<<endl;
  cout <<endl;
  testLWEKeyGen();
  cout <<endl;
  testLWE32EncDec();
  cout <<endl;
  testLWE64EncDec();
  cout <<endl;
  testLWE32addct();
  cout <<endl;
  testLWE32ptmulct();
  cout <<endl;


  //RLWE test
  cout << "====================================================="<<endl;
  cout <<"Tests of RLWE encryption scheme: "<<endl;
  cout<<endl;
  testRLWEKeyGen();
  cout <<endl;
  testRLWE64EncDec();
  cout <<endl;
  testRLWEKeySwitch();
  cout <<endl;
  return 0;
}