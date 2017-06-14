Borri_NCS
=========
This project is an implementaion of the NCS abstractions as they appeared in the paper:

Alessandro Borri, Giordano Pola, Maria Domenica Di Benedetto,  "A Symbolic Approach to the Design of Nonlinear Networked Control Systems"
Available on:: https://arxiv.org/abs/1203.1069
(Submitted on 5 Mar 2012 (v1), last revised 10 Mar 2012 (this version, v2))

we provide to examples:
1- The double integrator dynamics
2- The vehicle dynamics

Supported Systems:
------------------
- It is tested and can run on Linux.
- It can run on Mac (not tested)
- It can run on Windows 10 with Ubuntu bash installed (not tested)
- It can run on Windows (any version) with MSYS-2

HOW TO RUN:
------------
1- You need to have active installation of CUDD library version 3.0.0
2- You need to have a gcc/g++ compiler
3- Edit the Makefile to point out to your CUDD library includes and libs
4- Edit the system.hh file to select between the Double-Integrator and the Vehicle dynamics by manipulating the defines.
5- Build the binaries
6- Run it


