//code for testing ntt, intt, pmul in FasterNTT
#include <iostream>
using namespace std;

extern void RNG(){
  cout <<"Task: Test for generating random numbers from uniform distribution: "<<endl;
  unsigned int len_out = 10;
  unsigned char seed[SEED_LEN];
  int64_t out[len_out];
  random_bytes(seed,SEED_LEN);
  int64_t p = 1099511627776 + 1;
  cout <<"Generate "<<len_out<<" random numbers in Zp, p = "<<p<<endl;
  int a = gen_uniform_int64(out, len_out, p, 41, seed);
   for (int i=0;i<len_out;++i){
    printf("%lld \n",out[i]);
   }
}