# P2N3 (Privacy-preserving neural network of NTU)
## Introduction
P2N3 (Privacy-preserving neural network of NTU) is a C++ open-source library which implements an optimized fully homomorphic encryption (FHE) scheme for privacy-preserving neural networks.  
Our optimized FHE scheme enables us the ability to perform weighted sums and convolutions on the approximate LWE-based additive homomorphic encryption scheme, and to evaluate non-polynomial activations on FHEW ciphertexts. This FHE scheme has the following properties: 
- It can be applied to neural networks of arbitrary depth.
- It supports many kinds of widely used activations, such as the most popular ReLU. 
- When applied to inference of privacy-preserving neural networks, it is fast and accurate.   

## Prerequisites
- [OpenSSL](https://www.openssl.org/)

## Installation
```
mkdir build && cd build
cmake ..
make
```


## Examples
We offer some applications of our FHE scheme. It includes inference on MNIST dataset, facial recognition, speaker verification, text classification and object classification.   
REMARK: All the image/audio data are stored in the form of extracted feature vector (in plaintext). Please refer to the detailed experiments in our publication [[1]](https://ieeexplore.ieee.org/document/10398424).  

## License
This software is distributed under the BSD-3-Clause-Clear license. 
