# fastq2sam

fastq2sam - Converting paired-end fastq files with illumina-format quality scores into an unmapped sam/bam file while properly handling the RG tags.

## installation

This software depends on the SeqAn3 and Sharg-parser libraries, and the SeqAn3 library requires g++>=11 and the C++20 support.

### set up in centos

To use cmake in WSL, see https://stackoverflow.com/questions/62879479/every-call-to-configure-file-fails-on-wsl-configure-file-problem-configuring-fi

Then, in my local env, I ran the following commands:
```
conda create --name gxx13 python=3.7 
conda activate gxx13
conda install -c conda-forge -c anaconda gxx_linux-64 gcc_linux-64 sysroot_linux-64 bzip2 zlib
# export CMAKE_PREFIX_PATH="/usr/lib64/:$CMAKE_PREFIX_PATH" # I wrote this in ~/.bashrc, but it might be unnecessary
```

### Release or Debug mode

If you prefer static build, edit the line 4 of the CMakeLists.txt so that
```
set (CMAKE_BUILD_TYPE Release)
```

We note that in Mac, there is no static file for some dependent libraries, so please set `CMAKE_BUILD_TYPE` as `Debug`.

The software does not use libbz2, but the SeqAn3 and Sharg-parser libraries link libbz2 by default.
You can tell these two libraries to skip linking libbz2 so that you need not to prepare a static library for libbz2.
To do so, edit the following config files:

1. seqan3/seqan3-config.cmake
```
option (SEQAN3_NO_BZIP2 "Don't use BZip2, even if present." ON)
```

1. sharg-parser/sharg-config.cmake
```
option (SHARG_NO_BZIP2 "Don't use BZip2, even if present." ON)
```


### build

```
mkdir build && cd $_
cmake .. 
# If your default C++ compiler is not g++>=11, please specify the compiler by e.g., -DCMAKE_CXX_COMPILER=$CONDA_PREFIX/bin/g++.
# (In my Mac OS, I specified /usr/local/bin/g++-13)
make
make test
```