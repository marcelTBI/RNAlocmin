#ifndef __RNALOCMIN_H
#define __RNALOCMIN_H

#include <string>
#include <set>

#include "globals.h"

using namespace std;


// reads a line no matter how long
char* my_getline(FILE *fp);

// pt to str
string pt_to_str(short *pt);

// if the structure has lone pairs
int find_lone_pair(string &str);

// if the structure has lone pairs
int find_lone_pair(short* str);

// print rates/saddles to a file
void print_rates(char *filename, double temp, int num, float *energy_barr, vector<int> &output_en, bool only_saddles = false);

// just encapsulation
int move_set(struct_en &input);


#endif
