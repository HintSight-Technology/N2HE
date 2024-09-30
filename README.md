# N2HE (Neural Network Homomorphic Encryption)
## Introduction
N2HE (**N**eural **N**etwork **H**omomorphic **E**ncryption) is a C++ open-source library which implements an optimized fully homomorphic encryption (FHE) scheme for privacy-preserving neural networks.  
Our optimized FHE scheme enables us the ability to perform weighted sums and convolutions on the approximate LWE-based additive homomorphic encryption scheme, and to evaluate non-polynomial activations on FHEW ciphertexts. This FHE scheme has the following properties: 
- It can be applied to neural networks of arbitrary depth.
- It supports many kinds of widely used activations, such as the most popular ReLU. 
- When applied to inference of privacy-preserving neural networks, it is fast and accurate.

It is the implementation of publication [[1]](https://ieeexplore.ieee.org/document/10398424). For its application in other fields such as pandemic modelling, please refer to [[2]](https://ieeexplore.ieee.org/document/10376224). 

[[1]](https://ieeexplore.ieee.org/document/10398424) K.Y. Lam, X. Lu, L. Zhang, X. Wang, H. Wang and S.Q. Goh. “Efficient FHE-based Privacy-Enhanced Neural Network for Trustworthy AI-as-a-Service”, IEEE Transactions on Dependable and Secure Computing, IEEE.   
[[2]](https://ieeexplore.ieee.org/document/10376224) L. Zhang, X. Wang, J. Wang, R. Pung, H. Wang and K.Y. Lam. “An Efficient FHE-Enabled Secure Cloud-Edge Computing Architecture for IoMTs Data Protection With Its Application to Pandemic Modelling”, IEEE Internet of Things Journal, IEEE.

[N2HE-hexl](https://github.com/HintSight-Technology/P2N3-hexl) is another version of our scheme, using Intel Homomorphic Encryption (HE) Acceleration Library [hexl](https://github.com/intel/hexl). 

## Prerequisites
- [OpenSSL](https://www.openssl.org/)  3.2.1 or later

## Installation
Supported platforms: Linux and MacOS.  

```
mkdir build && cd build
cmake ..
make
```
This library uses third-party codes: 
- [FasterNTT](https://github.com/nict-sfl/FasterNTT), a very fast and advanced library for polynomial multiplication based on Number-theoretic Transform (NTT). 

## Document
Please refer to [N2HE_Programmer_Guide.pdf](https://github.com/HintSight-Technology/N2HE/blob/main/N2HE_Programmer_Guide.pdf) in this repo. 

## Examples
We offer some applications of our FHE scheme. It includes inference on MNIST dataset, facial recognition, speaker verification, text classification and object classification.   
REMARK: All the image/audio data are stored in the form of extracted feature vector (in plaintext). Please refer to the detailed experiments in our publication [[1]](https://ieeexplore.ieee.org/document/10398424).  

[[1]](https://ieeexplore.ieee.org/document/10398424) K.Y. Lam, X. Lu, L. Zhang, X. Wang, H. Wang and S.Q. Goh. “Efficient FHE-based Privacy-Enhanced Neural Network for Trustworthy AI-as-a-Service”, IEEE Transactions on Dependable and Secure Computing, IEEE.   

## License
This software is distributed under the BSD-3-Clause-Clear license. 
