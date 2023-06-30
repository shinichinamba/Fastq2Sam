# fastq2sam

fastq2sam - Converting paired-end fastq files with illumina-format quality scores into an unmapped sam/bam file while properly handling the RG tags.

## installation

The SeqAn3 library requires g++>=11 and the C++20 support.

### in centos

To use cmake in WSL, see https://stackoverflow.com/questions/62879479/every-call-to-configure-file-fails-on-wsl-configure-file-problem-configuring-fi

Then, in my local env, I ran the following commands:
```
conda create --name gxx13 python=3.7 
conda activate gxx13
conda install -c conda-forge -c anaconda gxx_linux-64 gcc_linux-64 sysroot_linux-64 #  bzip2 zlib
mkdir build && cd $_
# export CMAKE_PREFIX_PATH="/usr/lib64/:$CMAKE_PREFIX_PATH" # I wrote this in ~/.bashrc, but it might be unnecessary
cmake -DCMAKE_CXX_COMPILER=$CONDA_PREFIX/bin/g++ .. # I specified /usr/local/bin/g++-13 in my Mac.
make
make test
```

Although the binary was not static, we can also run it in Shirokane.

In Zeroone, We can run the same binary after running `export LD_LIBRARY_PATH="/work1/home/snamba/conda_pack/gxx13/lib:$LD_LIBRARY_PATH"`
