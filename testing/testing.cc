// A testing file to make sure that libgenetio functionality is built correctly....
#include <stdio.h>
#include <string>

#include "../marker.h"

int main(int argc, char **argv){
	if (argc != 2){
		fprintf(stderr, "Usage: %s [bim/map file]\n", argv[0]);
		return 1;
	}

	// Call some of the genetio methods in here....


	//1. Read in a .bim file (> 2 million SNPs)
	char * chr;
	strcpy(chr, "1");
	Marker::readBIMFile(argv[1], chr, 0, 2000000);

















}