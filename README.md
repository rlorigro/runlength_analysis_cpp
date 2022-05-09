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

#### Generating a Shasta config file 

Acquire a FASTA file containing read sequences.

Run `measure_runlength_distribution_from_fasta` using the executable in your `build` directory:

```
/home/ubuntu/software/runlength_analysis_cpp/build/measure_runlength_distribution_from_fasta \
--minimap_preset map-ont \
--max_threads 92 \
--minimum_match_length 11 \
--output_dir /home/ubuntu/data/human/hg005/reads/rle_run1/ \
--ref /home/ubuntu/data/human/reference/hg005/HG005.paternal.f1_assembly_v2_genbank.fasta \
--sequences /home/ubuntu/data/human/hg005/reads/05_25_21_R941_GM24631_9X_2_Guppy_5.0.7_prom_sup.fasta
```

To avoid regions of the alignment that contain misalignments as a result of short repeats (dinucleotide or others), it is highly recommended that the `minimum_match_length` parameter is used. See below in the "x_repeat" channel for an example of misalignments in RLE-space:

![image](https://user-images.githubusercontent.com/28764332/167509328-f9106e88-1379-4494-8da0-5a40fe7a4ee5.png)


For ONT data, be sure to use the `map-ont` preset, or for assemblies the `asm20` preset, etc.

Once this finishes running, a CSV will be generated containing confusion matrices for each nucleotide's runlength from 1-50bp. This is used to generate a table of normalized likelihoods (shasta's Bayesian ConsensusCaller config file). For example:

```
../scripts/convert_frequency_matrix_to_shasta_config.py \
--prior human \
--name HG005_wg_guppy_5.0.7_with_10_pseudocount_4-6-2022 \
--input /home/ryan/data/human/hg005/rle/guppy5/run2/length_frequency_matrix_nondirectional.csv \
--output /home/ryan/data/human/hg005/rle/guppy5/run2/shasta_bayesian_config_nondirectional_zero-fix/ \
--pseudocount 10
```
This generates an output directory and in it is a CSV which can be used directly by Shasta.

You can also edit the plot_matrix script to generate comparisons of frequency matrices for debugging/QC purposes:
![image](https://user-images.githubusercontent.com/28764332/167510644-abc86926-15f9-46c8-b86b-bb26c7ed369b.png)



#### Generating a MarginPolish RLE Matrix

To generate a MarginPolish RLE Matrix you need to use the `measure_runlength_distribution_from_fasta` exectuable and the `convert_frequency_matrix_to_marginpolish.py` script from this repository, and a reference and .fasta (not .fastq) file with your reads.

From the `build/` directory:
```
./measure_runlength_distribution_from_fasta --ref reference.fasta --sequences reads.fasta
python ../scripts/convert_frequency_matrix_to_marginpolish.py -i output/length_frequency_matrix_directional.csv -o output/marginPolish
```

This will create a file in `output/marginPolish`.  The contents of this file can be placed into the "repeatCountSubstitutionMatrix" node of the polish parameters.
