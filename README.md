# Prototype Installation of CLRT (Continuous Time Lagrangian Reachability)
Authors: Md Ariful Islam and Jacek Cyranka

Contact: ariful.islam@ttu.edu

Installation
============

Required Packages
-----------------

 - [GiNac]
 - [CAPD]
 - [Boost] (Link with libc++)
 - [dReal]
 - [NLOPT]
 - [Eigen Package]

[GiNac]: https://ginac.de/
[CAPD]: http://capd.sourceforge.net/capdDynSys/
[Boost]: http://www.boost.org/
[dReal]: http://dreal.cs.cmu.edu
[NLOPT]: https://nlopt.readthedocs.io/
[Eigen Package]: http://eigen.tuxfamily.org


Compile
-------

    cd <src_dir>
    make

Usage
=====

The command line is as follows:
 
    mainCLRT <model_name> <Time_horizon> <ist_th> <norm_th> <delta>
 
 where:
 
- <model_name>: Name of Model, e.g., bruss, fvdp, robot etc
- TIME_horizon > 0 : Total computation duration
- <ist_th> (Optional): Threshold for IST in (0,1) with default value 0.01
- <norm_th> (Optional): Threshold for norm change (0, 0.1) with default value 0.001
- Example run: ./mainCLRT bruss 5 0.02 0.001 1e-5

The output will be stored in following two files in a directory named outputs:
- dset_<model>_<ist_th>_<norm_th>.txt: contains the discrete reachset (time, center, radius, norm)
- cset_<model>_<ist_th>_<norm_th>.txt: contains the continuous reachset (time-interval, center, radius, norm)
       

