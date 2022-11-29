# Ellipsoid based Anisotropic Point Scatterer Model
Copyright (c) 2021 Eric Huang  
School of Electrical and Computer Engineering  
Georgia Institute of Technology  

The Matlab code associates with the paper "Anisotropic Scatterer Models for Representing RCS of Complex Objects".

The need to simulate complex electromagnetic (EM) wave interactions by multiple radar targets, transmitters, and receivers to better study the performance of radar systems, antenna designs, and/or stealth technologies has grown over time. High performance computing (HPC) based emulators can be used to model the scattering from multiple stationary and moving targets for radar applications. These emulators rely on the Radar Cross Section (RCS) of the targets being available in complex scenarios. Representing the RCS using tables generated from EM simulations is often times cumbersome leading to large storage requirement. An alternative approach is to represent the targets as a collection of isotropic or anisotropic scatterers. In this paper we present a method to represent the RCS of complex targets using a 3D anisotropic scatterer model, where we use the analytical RCS representation of a large ellipsoid as the basis function to determine the angular dependency of the RCS from each scatterer. The scatterer model that best represents the RCS data is obtained by solving a least square inverse problem. To improve the correlation with EM solvers, we further break down the optimization problem by considering shadowing effect and use multiple models, each representing a subset of the RCS data. The results show that the scatterer model can effectively represent the RCS data of complex targets.

For questions, queries and bug reports, please feel free to contact: huangeric@gatech.edu

## Note:
The code is tested using Matlab R2019b. The aircraft geometry is obtained from R. Okada, “B787-8 dreamliner.” Online. https://grabcad.com/library/ b787-8-dreamliner-1 Accessed October 27, 2020.
