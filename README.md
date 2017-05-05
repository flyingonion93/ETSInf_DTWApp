# DTWApp

## Configuration
The package is delivered as a standard autoconf program, so in order to configure it to the target device there are some prerequisites to be met. Autoconf is required for all versions, but there are some platform dependant libraries:
CPU
* Basic Linear Algebra libraries (OpenBLAS, Intel MKL)
* FFT solver (FFTW)

GPU
* CUDA runtime
* Linear Algebra libraries (CuBLAS)
* FFT solver (CuFFT)

## Installation
The procedure to install the package is the same as with every autoconf package, but this software is delivered at a minimum state. That means the first thing to do is to generate the installer package by using autotools. The procedure consist on executing 
```
aclocal
autoconf
automake --add-mising
``` 

Once the package has been generated the procedure is the autoconf standard
```
./configure
make
sudo make install
```

When calling the configuration script the user can see the environment variables for the package. Besides the common ones, such as the compiler used or its flags, there are four defined to set the target platform.

### IMP
Determines the implementation that is going to compile. The options are **CPU**, detecting automatically if its x86_64 or ARM, or **GPU**. The default value is CPU

### FFT
Used to set the FFT solver that is going to be used. This variable only affects the CPU version, because the GPU one is based exclusively on CuFFT. The package currently supports FFTW3 library, so the only option is FFTW, and is also the default one.

### ALG
The linear algebra library that the software uses. It also affects exclusively the CPU version. The package supports the OpenBLAS library natively, but it can also work with the Intel MKL by modifying some flags before calling the configure script. The default value is BLAS.

### PRECISSION
Sets if the floating point values are simple or double precission. This also determines if it must compiles some libraries in a specific manner. The options are **SIMPLE** or **DOUBLE** and this first one is the default value.

## About
This software was developed as a final degree project for the Compute Science Degree at the Universitat Politecnica de Valencia. It can be forked and modified at will, as long as it's properly mentioned as a source.
There is also available a [report](https://riunet.upv.es/bitstream/handle/10251/76389/G%C3%93MEZ%20-%20Un%20aplicaci%C3%B3n%20para%20el%20alineamiento%20de%20partituras%20en%20tiempo%20real..pdf?sequence=1&isAllowed=y) (in spanish) describing the software architecture, the testing methodology and how it performs on systems such as Desktop PC or ARM based devices for CPU execution or CUDA based devices for GPU execution.
