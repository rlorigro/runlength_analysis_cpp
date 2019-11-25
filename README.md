# runlength_analysis_cpp
This library was built to replace rlorigro/runlength_analysis. Given raw sequences, or an assembly, this library will generate a confusion matrix of homopolymer length predictions. Some evaluation tools are included to aid in the development of runlength related methods.

## Installation

#### Download repo:
```
git clone https://github.com/rlorigro/runlength_analysis_cpp.git
cd runlength_analysis_cpp/
sudo sh install_prerequisites_ubuntu_18.sh 
```

#### If minimap2 not installed and added to path:
```
cd your/preferred/installation/directory
sudo sh path/to/runlength_analysis_cpp/install_minimap2.sh
```

#### If cmake is not installed:
```
Figure that out
```

#### Finally:
```
sudo apt-get install cmake 
mkdir build
cd build
cmake ..
make -j 8
```
