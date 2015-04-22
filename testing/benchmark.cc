#include <stdio.h>
#include <string>
#include "../dynarray.h"

// TODO: decide whether option 3 || 4 is faster


int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s [file] [type]\n", argv[0]);
    return 1;
  }

  int type = atoi(argv[2]);

  FILE *in = fopen(argv[1], "r");
  if (!in) {
    fprintf(stderr, "Couldn't open %s\n", argv[1]);
    perror(argv[1]);
    return 1;
  }


  if (type == 1) {
    dynarray<char *> lines;

    int lineNum = 0;
    while (1) {
      char *line = new char[1024];
      char *ret = fgets(line, 1024, in);

      if (ret == NULL) {
	break;
      }

      lines.append(line);
      lineNum++;

      if (lineNum % 1024 == 0)
	printf("On line %d\r", lineNum);
    }

    printf("Read a total of %d lines\n", lineNum);
  }
  else if (type == 2) {
    dynarray<std::string *> lines;

    char c;
    int lineNum = 0;
    do {
      std::string *line = new std::string;

      while ((c = fgetc(in)) && c != '\n' && c != EOF)
	line += c;

      lines.append(line);
      lineNum++;

      if (lineNum % 1024 == 0)
	printf("On line %d\r", lineNum);
    } while (c != EOF);

    printf("Read a total of %d lines\n", lineNum);
  }
  else if (type == 3) {
    char *data = new char[1024*1024*400]; // 400 MiB

    int kbread = 0;
    size_t last_count;
    for(int i = 0; i < 1024 * 400; i++) {
      last_count = fread(&data[i*1024], 1, 1024, in);
      if (last_count < 1024)
	break; // read it all
      kbread++;

      if (kbread % 16 == 0) {
	printf("Read %d kb\r", kbread);
      }
    }

    printf("Read a total of %d kb\n", kbread);

    int lineNum = 0;
    for(int k = 0; k < kbread; k++) {
      for(int c = 0; c < 1024; c++) {
	if (data[k*1024 + c] == '\n')
	  lineNum++;

      if (lineNum % 1024 == 0)
	printf("On line %d\r", lineNum);
      }
    }

    printf("Read a total of %d lines\n", lineNum);
  }
  else {
    const int SIZE = 1024;
    char data1[SIZE];
    char data2[SIZE];

    char *cur = data1;
    char *prev = data2;

    size_t last_count = fread(prev, 1, SIZE, in);

    int lineNum = 0;
    while (last_count == SIZE) {
      last_count = fread(cur, 1, SIZE, in);

      for(int c = 0; c < 1024; c++) {
	if (prev[c] == '\n')
	  lineNum++;

	if (lineNum % SIZE == 0)
	  printf("On line %d\r", lineNum);
      }

//      fwrite(prev, 1, SIZE, stderr);
//      fprintf(stderr, "\n:\n");

      char *tmp = prev;
      prev = cur;
      cur = tmp;
    }

    for(unsigned int c = 0; c < last_count; c++)
      if (prev[c] == '\n')
	lineNum++;

//    fwrite(prev, 1, SIZE, stderr);
//    fprintf(stderr, "\n:\n");

    printf("Read a total of %d lines\n", lineNum);
  }

  return 0;
}
