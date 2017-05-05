AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([./macros])
AC_INIT([dtw], [1.00])
AM_INIT_AUTOMAKE([subdir-objects no-define])

# Check user values (when calling configure)
AC_ARG_VAR([IMP], [The implementation of the application that must be compiled.\nThe options are CPU and GPU. The default value is CPU]
)

AC_ARG_VAR([FFT], [The library that will be used to process the Fast Fourier Transform by the application. The default value is FFTW]
)

AC_ARG_VAR([ALG], [The algebra library that will be used by the application. The default value is BLAS]
)

AC_ARG_VAR([PRECISSION], [The precission of the floating point values. If simple is declared the macro SIMPLE_PRECISSION will be declared when compiled. The default value is double]
)

AC_ARG_ENABLE( [docs],
	AS_HELP_STRING([--disable-docs], [Disables the generation of the project documentation]) 
)

AS_IF(	[test "x$IMP" = x], [IMP="CPU"],
		[test "x$IMP" = xCPU], [],
		[test "x$IMP" = xGPU], [],
		[echo Parameter error. IMPLEMENTATION not recognized]
)
AS_IF(	[test "x$FFT" = x], [FFT="FFTW"],
		[test "x$FFT" = xFFTW], [],
		[echo Parameter error. FFT library not recognized]
)
AS_IF(	[test "x$ALG" = x], [ALG="BLAS"],
		[test "x$ALG" = xBLAS], [],
		[echo Parameter error. Algebra library not recognized]
)
AS_IF(	[test "x$PRECISSION" = x], [PRECISSION="double"],
		[test "x$PRECISSION" = xsimple], [],
		[test "x$PRECISSION" = xdouble], [],
		[echo Parameter error. Precission must be simple or double]
)
[echo Implementation = $IMP]
[echo Fft = $FFT]
[echo Alg = $ALG]
[echo Prec = $PRECISSION]

# Some platform names
desktop="x86_64"
ARM64="aarch64"

LIBS=""
INCLUDE=""
CFLAGS=""

# Check for non-architecture-dependent programs.
AS_IF([test "x$enable_docs" != "xno"],[
	AC_CHECK_PROGS([DOXYGEN], [doxygen] )
	if test -z "$DOXYGEN";
	then
		 AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
	fi
	AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"]) 
	AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([Doxyfile])])
])

if test "x$IMP" = "xCPU";
then
	AC_DEFINE([CPU])
	platform=$(uname -m)
	AS_IF([test "$platform" = "$desktop"], [CFLAGS="-m64 $CFLAGS"],
		[test "x$platform" = "x$ARM64"], [CFLAGS="-march=armv8-a, -mfpu=neon-fp-armv8 -mfloat -mfloat-abi=hard -DARM32 -Wno-unused-result -Wno-unused-but-set-variable $CFLAGS"],
		[CFLAGS="-mcpu=cortex-a15 -DARM32 -Wno-unused-result -Wno-unused-but-set-variable $CFLAGS"])

	CFLAGS="-std=c99 -Wall -O3 -fopenmp $CFLAGS"
	if test "x$PRECISSION" = "xsimple"		
	then
		AC_DEFINE([SIMPLE_PRECISSION])
	fi

	# Check for programs.
	AC_PROG_CC

	# Check for libraries
														# OMP library
	AC_OPENMP([echo OMP library found])
														# FFT library
	if test "x$FFT" = "xFFTW"
	then
		AC_DEFINE([FFTW])
		# Depending on the floating point type we must check for FFTW or FFTWF
		if test "x$PRECISSION" = "xsimple"
		then
			LIBS="-L/opt/fftwf/lib -lfftw3_omp -lfftw3f $LIBS"
			INCLUDE="-I./ -I/opt/fftwf/include $INCLUDE"
		else
			LIBS="-L/opt/fftw/lib -lfftw3_omp -lfftw3 $LIBS"
			INCLUDE="-I./ -I/opt/fftwf/include $INCLUDE"
		fi
		AC_CHECK_HEADERS([fftw3.h])
	fi
														# Algebra library
	if test "x$ALG" = "xBLAS"
	then
		AC_DEFINE([BLAS])
		LIBS="-L/opt/openblas/lib -lopenblas $LIBS"
		INCLUDE="-I/opt/openblas/include $INCLUDE"
		AC_CHECK_HEADERS([cblas.h])
	fi
														# Math library
	AC_CHECK_LIB([m], [pow], [echo Math library found], [echo Math library not found])
	[LIBS="-lm $LIBS"]

	# Check for header files
	AC_CHECK_HEADERS([stdio.h stdarg.h omp.h stdlib.h math.h string.h fftw3.h cblas.h])

	# Check for typedefs, structures and compiler characteristics
	AC_TYPE_SIZE_T
fi
#if test "x$IMP" = "xGPU";
[echo -------------------------------------]
[echo $LIBS]
[echo $INCLUDE]
[echo $CFLAGS]
[echo -------------------------------------]
doxygen Doxyfile.in
AC_OUTPUT([Makefile] [src/Makefile])
