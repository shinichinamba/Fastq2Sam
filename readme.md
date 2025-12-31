# fastq2sam / split_fastq

fastq2sam - Converting paired-end fastq files with illumina/DNBSEQ-format quality scores into an unmapped sam/bam file while properly handling the RG tags.

split_fastq - Splitting paired-end illumina/DNBSEQ fastq files by read groups.

## installation

### static binary

Linux users can download static binary files from the release page.

### from source

```
git clone --recurse-submodules https://github.com/shinichinamba/Fastq2Sam.git
# git submodule update --init --recursive
```
This software depends on the SeqAn3 and Sharg-parser libraries, and the SeqAn3 library requires g++>=11 and the C++20 support.

#### build

```
mkdir build && cd $_
cmake .. 
# If your default C++ compiler is not g++>=11, please specify the compiler by e.g., -DCMAKE_CXX_COMPILER=$CONDA_PREFIX/bin/g++.
# (In my Mac OS, I specified /usr/local/bin/g++-13 or /opt/homebrew/bin/g++)
make
make test
```

#### Release or Debug mode

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
option (SEQAN3_NO_CEREAL "Don't use cereal, even if present." ON)
```

1. sharg-parser/sharg-config.cmake
```
option (SHARG_NO_BZIP2 "Don't use BZip2, even if present." ON)
```

#### set up in centos/WSL2

To use cmake in WSL2, see https://stackoverflow.com/questions/62879479/every-call-to-configure-file-fails-on-wsl-configure-file-problem-configuring-fi

Then, in my local environment for example, I ran the following commands:
```
conda create --name gxx13 python=3.7 
conda activate gxx13
conda install -c conda-forge -c anaconda gxx_linux-64 gcc_linux-64 sysroot_linux-64 bzip2 zlib
# export CMAKE_PREFIX_PATH="/usr/lib64/:$CMAKE_PREFIX_PATH" # I wrote this in ~/.bashrc, but it might be unnecessary
```

## history

### 2025/12/28 v0.1.1
* improved read group parsing

### 2025/12/28 v0.1.0
* added the executable "split_fastq"

### 2025/2/25 v0.0.4
* supported the files using irregular characters for zero quality

### 2024/11/6 v0.0.3
* added the option "--phred"

### 2024/2/9 v0.0.2
* fixed an error which was raised for too short (< 2 characters in default) read names
* supported the DNBSEQ-format read names

### 2023/7/26 v0.0.1
* initial release