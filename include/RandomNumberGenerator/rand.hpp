// ================================================
// Clean up 2023 11 10
// lac_param.h and rng.hpp are external
// ================================================

#include <stdint.h>
#include "rng.hpp"
#include "lac_param.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <openssl/rand.h>
#include <openssl/aes.h>
#include <openssl/sha.h>
#include <openssl/crypto.h>
#include <openssl/evp.h>
#include <iostream>
using namespace std;

//random bytes
int random_bytes(unsigned char *r, unsigned int len)
{
  //check parameter
  if(r==NULL)
  {
    return 1;
  }
  // call the random function 
  RAND_bytes(r,len);
  return 0;
}


//pseudo-random bytes
int pseudo_random_bytes(unsigned char *r, unsigned int len, const unsigned char *seed)
{
  
  int c_len;
  unsigned char data[AES_BLOCK_SIZE],c[AES_BLOCK_SIZE];
//  unsigned int *p=(unsigned int *)data;
  int i,loop=len/AES_BLOCK_SIZE;
  //check  parameter
  if(r==NULL || seed==NULL)
  {
    return 1;
  }
  memset(r,0,len);
  EVP_CIPHER_CTX *ctx;
  ctx = EVP_CIPHER_CTX_new();
  EVP_EncryptInit_ex(ctx, EVP_aes_256_ctr(), NULL, seed, NULL);
  memset(data,0,AES_BLOCK_SIZE);
  for(i=0;i<loop;i++)
  {
    EVP_EncryptUpdate(ctx, r+i*AES_BLOCK_SIZE, &c_len, data, AES_BLOCK_SIZE);
  }
  //check tail
  if(len%AES_BLOCK_SIZE>0)
  {
    EVP_EncryptUpdate(ctx, c, &c_len, data, AES_BLOCK_SIZE);
  }
  memcpy(r+loop*AES_BLOCK_SIZE,c,len%AES_BLOCK_SIZE);
  EVP_CIPHER_CTX_free(ctx);
  
   
//  aes256ctr_prf(r,len,seed,0);
  return 0;
}



//hash
int hash_h(const unsigned char *in, unsigned int len_in, unsigned char * out)
{
  //check  parameter
  if(in==NULL || out==NULL)
  {
    return 1;
  } 
  
  SHA256(in,len_in,out); // MUST BE len_in <= MESSAGE_LEN
  
  return 0;
}


//hash
int hash_to_k(const unsigned char *in, unsigned int len_in, unsigned char * out)
{
  //check  parameter
  if(in==NULL || out==NULL)
  {
    return 1;
  } 
  unsigned char tmp_out[32];
  

  SHA256(in,len_in,tmp_out);
  memcpy(out,tmp_out,MESSAGE_LEN);
  
  return 0;
}


//generate seed
int gen_seed(unsigned char *in, unsigned int len_in, unsigned char * out)
{
  //check  parameter
  if(in==NULL || out==NULL)
  {
    return 1;
  }
    
  SHA256(in,len_in,out);
  return 0;
}

// =========================================
// wangxn Added
// generate Bernoulli, for lwe key
// for simplicity only work when 8 | len_out
extern int gen_bernoulli(int *out, unsigned int len_out, unsigned char *seed){
  if(out == NULL){
    return 1;
  }
  if (len_out % 8 != 0){
    // for simplicity only work when 8 | len_out
    return 2;
  }

  // generate random array, each element is 8 bits
  unsigned int len = len_out >> 3; 
  unsigned char buf[len];
  pseudo_random_bytes(buf,len,seed);
  
  // generate output
  unsigned char tmp;
  for (int i=0;i<len;++i){
    tmp = buf[i];
    for (int j=0;j<8;++j){
      out[i*8+j] = tmp & (0x01); // extract each bit in buf
      tmp = tmp >> 1;
    }
  }
  return 0;
}

// generate ternary, for rlwe key
// for simplicity only work when 8 | len_out
extern int gen_ternary(int *out, unsigned int len_out, unsigned char *seed){
  int out1[len_out], out2[len_out];
  
  if(out == NULL){
    return 1;
  }
  if (len_out % 8 != 0){
    // for simplicity only work when 8 | len_out
    return 2;
  }

  unsigned char seed1[SEED_LEN];
  unsigned char seed2[SEED_LEN];
  pseudo_random_bytes(seed1,SEED_LEN,seed);
  pseudo_random_bytes(seed2,SEED_LEN,seed1);
  
  int a1 = gen_bernoulli(out1, len_out, seed1);
  int a2 = gen_bernoulli(out2, len_out, seed2);
  for (int i=0;i<len_out;++i){
    out[i] = out1[i] - out2[i]; 
  }
  return 0;
}

// generate ternary distribution with variance 2^(-k), k >= 1
// for simplicity only work when 8 | len_out
int gen_ternary_var(int *out, unsigned int len_out, unsigned int k, unsigned char* seed){
  if(out == NULL){
    return 1;
  }
  if (len_out % 8 != 0){
    // for simplicity only work when 8 | len_out
    return 2;
  }


  // generate random array, each element is 8 bits
  unsigned int len = len_out * (k + 1); // each slot needs (k+1) Bernoulli variables.
  int bernoulli[len];
  int a = gen_bernoulli(bernoulli, len, seed);
  

  // generate output
  
  int tmp;
  for (int i=0;i<len_out;++i){
    tmp = 1;
    for (int j=0;j<k;++j){
      tmp = tmp * bernoulli[i*(k+1)+j]; 
    }
    if (tmp == 0) {
      out[i] = 0;
    }
    else if (bernoulli[i*(k+1)+k] == 0) {
      out[i] = -1;
    }
    else {
      out[i] = 1;
    }
  }
  return 0;
}

// generate uniform distribution on Zp, with int p
    // for simplicity only work when p = 2^k
    // one can simply extend to general p
extern int gen_uniform(int *out, unsigned int len_out, int p, unsigned char *seed){
  if(out == NULL){
    return 1;
  }
  
  if ((p & (p - 1)) != 0){
    // for simplicity only work when p = 2^k
    // one can simply extend to general p
    return 2;
  }
  

  unsigned int bit_length = 0; 
  p = p-1;
  while (p){
    ++bit_length;
    p = p >> 1; 
  }
  // number of random bytes for each element in Zp
  unsigned int number_of_bytes = (bit_length >> 3) + 1; 
  // number of unused bits
  unsigned int remaining_bits = (number_of_bytes << 3) - bit_length; 


  unsigned int total_len = number_of_bytes * len_out; // total number of random bytes
  unsigned char buf[total_len];
  pseudo_random_bytes(buf,total_len,seed);

  // generate output
  int tmp;
  for (int i=0;i<len_out;++i){
    tmp = 0;
    for (int j=0;j<number_of_bytes;++j){
      // combile number_of_bytes
      tmp = (tmp << 8) + buf[i*number_of_bytes+j];
    }
    out[i] = (tmp >> remaining_bits) - ((p + 1) >> 1);
  }
  return 0;
}

// generate uniform distribution on Zp, with int64 p, and k=ceil(log(p))
extern int gen_uniform_int64(int64_t *out, unsigned int len_out, int64_t p, unsigned int k, unsigned char *seed){
  if(out == NULL || seed == NULL){
    return 1;
  }

  unsigned int bit_length = k;

  // number of random bytes for each element in Zp
  unsigned int number_of_bytes = (bit_length >> 3) + 1; 
  // number of unused bits
  unsigned int remaining_bits = (number_of_bytes << 3) - bit_length; 

  unsigned int total_len = number_of_bytes * len_out; // total number of random bytes
  unsigned char buf0[total_len];
  unsigned char buf[total_len];
  pseudo_random_bytes(buf0,total_len,seed);
  hash_h(seed,MESSAGE_LEN,buf); // back up

  // generate output
  uint64_t tmp;
  for (int i=0;i<len_out;++i){
    tmp = 0;
    for (int j=0;j<number_of_bytes;++j){
      // combile number_of_bytes
      tmp = (tmp << 8) + buf0[i*number_of_bytes+j];
    }
    out[i] = (int64_t)(tmp >> remaining_bits);
  }

  // check each out[i] in Zp
  int hash_index = 0; 
  for (int i=0;i<len_out;++i){
    while (out[i] >= p) {
      uint64_t tmp = 0;
      for (int j=0;j<number_of_bytes;++j){
        tmp = (tmp << 8) + buf[hash_index+j];
      }
      hash_index += number_of_bytes;
      out[i] = (int64_t)(tmp >> remaining_bits);
      if (hash_index+number_of_bytes >= MESSAGE_LEN) {
        hash_h(buf,MESSAGE_LEN,buf);
        hash_index = 0;
      }
    }
  } 
  return 0;
}