#if CPU
	#include "../launcher/DTWExec_CPU.h"
#elif GPU
	#include "../launcher/DTWExec_GPU.h"
#endif

int main( int argc, char *argv[] )
{
	DTWExec( argc, argv );
	return 0;
}
