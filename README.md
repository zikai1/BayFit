# BayFit
This repository contains our official implementation of _A Bayesian Approach Toward Robust and Multidimensional Ellipsoid Fitting_ 
## 1. Motivation
Ellipsoid fitting is an important foundamental problem, which has various applications in computer vision, computer graphics, and biomedical data analysis. Most existing approaches adopt the least-squares principle to find the target parameter. This manner generates satisfactory estimates under simple or clean scenes, but suffers from outliers and spatial dimensions. Unlike predecessor algorithms, we adopt the Bayesian parameter estimate process to solve this problem. Our method is ellipsoid-specific, robust against outliers, and can be generalized to high dimensional spaces. 


## 2. Usage
Step 1:

`git clone https://github.com/zikai1/BayFit` or directly `download the source files`;


Step 2:

Compile the C++ files by the mex operation in Matlab command line as follows. To this end, you are recommended to install the `MinGW64 Compiler (C)` or `Microsoft Visual C++ 2019 (C)`. Once either one is successfully installed, then perform the subsequent steps:

(1) Setup the compile by `mex -setup`

(2) Then choose the compile language designed for C++  `mex -setup C++`

(3) `mex knn_cpp.cpp`

(4) `mex precompute.cpp`


Step 3:

`Run "demo.m" to see demo examples.`




## 3. Contact
If you have any question, please [submit an issue](https://github.com/zikai1/BayFit/issue) or [contact me](myzhao@baai.ac.cn).


