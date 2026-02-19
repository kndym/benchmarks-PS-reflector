<img src="images/romsoclogo-logo.png" alt="ROMSOC logo"  width="150"/>

# Benchmarks for Point Source to Far Field reflector problem
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5171811.svg)](https://doi.org/10.5281/zenodo.5171811) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ROMSOC/benchmarks-ps-reflector/HEAD?labpath=Benchmark.ipynb)

This code is developed to simulate the benchmark cases of model hierarchies for the Point Source to Far Field reflector problem. In order to use the code, one should choose desired benchmark case in the main.cpp file, and then run the makefile.

The code is arranged in the following folders:

"Benchmarks" - Contains two different testcase header files representing different problems. The desired one should be included in the main.cpp file by uncommenting the corresponding command.

"PushForward" - Contains a header file "Pushforward_of_RefRegular.h" which contains the methods for computing the reflection from the obtained reflector and another header file "PushForward_Cloud_128.h" containing a pointcloud on the portion of the sphere, with approximately 128*128 points, sampled using quasi monte-carlo methods.


"QuasiMonteCarlo" - Contains a header file  which defines necessary data containers and functions to work with quasi monte-carlo grids and another header file containing such grid.

"SmallGrid" - Contains header file with a coarser quasi monte-carlo grid and another header file with complete code for computing benchmark test case on the coarse grid, in order to obtain better initialization for the finer grid.

main.cpp file contains the rest of the code, including the different versions of sinkhorn algorithm and interchanges between structured and quasi monte-carlo grids.



### RUN ON LOCAL HOSTS
This code can be adjusted to any OS, but in its current form it can be run on Linux only, since it uses file and folder creation commands. 
To run this code, Intel's MKL is necessary. 

Currently it is available at https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html
In order to run the code, edit the Makefile by adjusting all mkl paths according to your installation.

### RUN JUPYTER NOTEBOOK
The entire benchmark repository can be executed in a provided Docker container where a full installation of Intel OneAPI is available. Once you have clone or downloaded this repository, to build the container just type
```bash
docker build -t benchmarks-ps-reflector . 
```
and for running it locally:
```bash
docker run -it --rm -p 8888:8888 benchmarks-ps-reflector jupyter-lab --ip=0.0.0.0 --port=8888 --allow-root
```

Alternatively, user-friendly Jupyter Notebooks could be used to run different benchmarks on the cloud. For instance, the benchmark is available at:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ROMSOC/benchmarks-ps-reflector/HEAD?labpath=Benchmark.ipynb). Please, notice that mybinder cloud computations are limited to 2GB of RAM memory.

## Disclaimer
In downloading this SOFTWARE you are deemed to have read and agreed to the following terms:
This SOFTWARE has been designed with an exclusive focus on civil applications. It is not to be used
for any illegal, deceptive, misleading or unethical purpose or in any military applications. This includes ANY APPLICATION WHERE THE USE OF THE SOFTWARE MAY RESULT IN DEATH,
PERSONAL INJURY OR SEVERE PHYSICAL OR ENVIRONMENTAL DAMAGE. Any redistribution of the software must retain this disclaimer. BY INSTALLING, COPYING, OR OTHERWISE
USING THE SOFTWARE, YOU AGREE TO THE TERMS ABOVE. IF YOU DO NOT AGREE TO
THESE TERMS, DO NOT INSTALL OR USE THE SOFTWARE

## Acknowledgments
<img src="images/EU_Flag.png" alt="EU Flag"  width="150" height="100" />

The ROMSOC project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Grant Agreement No. 765374.
This repository reflects the views of the author(s) and does not necessarily reflect the views or policy of the European Commission. The REA cannot be held responsible for any use that may be made of the information this repository contains.
