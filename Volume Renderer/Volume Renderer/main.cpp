#include "colour.h"
#include "Image.h"

#include <stdio.h>
int main( int argc, const char* argv[] )
{
	// Prints each argument on the command line.
	printf( "Volume Renderer\n");
	for( int i = 0; i < argc; i++ )
	{
		printf( "arg %d: %s\n", i, argv[i] );
	}
}