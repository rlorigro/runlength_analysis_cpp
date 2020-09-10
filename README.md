# runlength_analysis_cpp
This library was built to replace rlorigro/runlength_analysis. Given raw sequences, or an assembly, this library will generate a confusion matrix of homopolymer length predictions. Some evaluation tools are included to aid in the development of runlength related methods.

## Installation

#### Download repo and run library installation script:
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
wget https://github.com/Kitware/CMake/releases/download/v3.14.4/cmake-3.14.4-Linux-x86_64.sh
sudo mkdir /opt/cmake 
sudo sh cmake-3.14.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license 
sudo ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
cmake --version
```

#### Finally:
```
mkdir build
cd build
cmake ..
make -j [n_threads]
```

## Usage

#### Generating a MarginPolish RLE Matrix

To generate a MarginPolish RLE Matrix you need to use the `measure_runlength_distribution_from_fasta` exectuable and the `convert_frequency_matrix_to_marginpolish.py` script from this repository, and a reference and .fasta (not .fastq) file with your reads.

From the `build/` directory:
```
./measure_runlength_distribution_from_fasta --ref reference.fasta --sequences reads.fasta
python ../scripts/convert_frequency_matrix_to_marginpolish.py -i output/length_frequency_matrix_directional.csv -o output/marginPolish
```

This will create a file in `output/marginPolish`.  The contents of this file can be placed into the "repeatCountSubstitutionMatrix" node of the polish parameters.
