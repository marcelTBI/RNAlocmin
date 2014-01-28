#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <map>
#include <vector>
#include <algorithm>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include "RNAlocmin.h"

// reads a line no matter how long
char* my_getline(FILE *fp)
{
  char s[512], *line, *cp;
  line = NULL;
  do {
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0';
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char));
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1);
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}

// find the structure's lone pairs
int find_lone_pair(string &str)
{
  for(unsigned int i=0; i<str.length(); i++) {
    if (str[i]=='(') {
      if (i+1==str.length() || str[i+1]!='(') {
        return i;
      } else while (i+1!=str.length() && str[i+1]=='(') i++;
    }

    if (str[i]==')') {
      if (i+1==str.length() || str[i+1]!=')') {
        return i;
      } else while (i+1!=str.length() && str[i+1]==')') i++;
    }
  }

  return -1;
}

// if the structure has lone pairs
int find_lone_pair(short* str)
{
  for(int i=1; i<str[0]; i++) {
    if (str[i]==0) continue; // '.'

    if (str[i]>str[str[i]]) {  // '('
      if (i+1==str[0] || str[i+1]==0 || str[i+1]<str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]>str[str[i+1]]) i++;
    }

    if (str[i]<str[str[i]]) {  // ')'
      if (i+1==str[0] || str[i+1]==0 || str[i+1]>str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]<str[str[i+1]]) i++;
    }
  }

  return -1;
}

// print rates to a file
void print_rates(char *filename, double temp, int num, float *energy_barr, vector<int> &output_en, bool only_saddles)
{
  FILE *rates;
  rates = fopen(filename, "w");
  if (rates==NULL) {
    fprintf(stderr, "ERROR: couldn't open file \"%s\" for rates! (using stderr instead)\n", filename);
    rates = stderr;
  }
  double _kT = 0.00198717*(273.15 + temp);
  for (int i=0; i<num; i++) {
    for (int j=0; j<num; j++) {
      float res = 0.0;
      if (i!=j) {
        // Arhenius kinetics (as A method in treekin)
        res = 1.0*exp(-(energy_barr[i*num+j]-(output_en[i]/100.0))/_kT);
      }
      if (only_saddles) fprintf(rates, "%6.2f ", i==j?output_en[i]/100.0:energy_barr[i*num+j]);
      else              fprintf(rates, "%10.4g ", res);
    }
    fprintf(rates, "\n");
  }
  fclose(rates);
}

// pt to str
string pt_to_str(short *pt)
{
  string str;
  str.resize(pt[0]);
  //fprintf(stderr, "pt_to_str %d\n", pt[0]);
  for (int i=1; i<=pt[0]; i++) {
    if (pt[i]==0) str[i-1]='.';
    else if (pt[i]<i) str[i-1]=')';
    else str[i-1]='(';
  }
  return str;
}

// encapsulation
int move_set(struct_en &input)
{
  // call the coresponding method
  int verbose = (Opt.verbose_lvl-2<0?0:Opt.verbose_lvl-2);
  if (Opt.rand) input.energy = move_adaptive(Enc.seq, input.structure, Enc.s0, Enc.s1, verbose);
  else {
    if (Opt.first) input.energy = move_first(Enc.seq, input.structure, Enc.s0, Enc.s1, verbose, Opt.shift, Opt.noLP);
    else input.energy = move_gradient(Enc.seq, input.structure, Enc.s0, Enc.s1, verbose, Opt.shift, Opt.noLP);
  }
  return input.energy;
}