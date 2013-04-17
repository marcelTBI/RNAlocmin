#ifndef __MOVE_SET_H
#define __MOVE_SET_H

#include <string>
#include <set>

#include "globals.h"

using namespace std;


//#################################################### SMALL (HELPER) FUNCTIONS

// count how many times move_set has been called
int count_move();

// abs(float)
inline float abs(float a) { return (a<0?-a:a); }

// reads a line no matter how long
char* my_getline(FILE *fp);

// pt to str
string pt_to_str(short *pt);

// if the structure has lone pairs
int find_lone_pair(string &str);

// if the structure has lone pairs
int find_lone_pair(short* str);


// ################################################### BIG (IMPORTANT) FUNCTIONS

// move to deepest(first) neighbour (returns minimum in str) returns 0 when minimum has been found
int move_set(hash_entry &str);

// move to random neighbour (returns minimum in str) returns 0 when minimum has been found
int move_rand(hash_entry &str);

// browse neighbours and apply funct (from Opt) on all of them (returns minimum) return saddle energy in saddle_en (in flooding) return if funct(str) returns true
hash_entry *browse_neighs(hash_entry &str, int &saddle_en);

// print rates/saddles to a file
void print_rates(char *filename, double temp, int num, float *energy_barr, vector<int> &output_en, bool only_saddles = false);


#endif
