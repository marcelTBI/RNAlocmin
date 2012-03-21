#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "hash_util.h"

void copy_arr(short *dest, short *src)
{
  if (!src || !dest) {
    fprintf(stderr, "Empty pointer in copying\n");
    return;
  }
  memcpy(dest, src, sizeof(short)*(src[0]+1));
}

short *allocopy(short *src)
{
  short *res = (short*) space(sizeof(short)*(src[0]+1));
  copy_arr(res, src);
  return res;
}

hash_entry *copy_entry(const hash_entry *he)
{
  hash_entry *he_n = (hash_entry*) space(sizeof(hash_entry));
  he_n->structure = allocopy(he->structure);
  he_n->energy = he->energy;

  return he_n;
}

void free_entry(hash_entry *he)
{
  if (he->structure) free(he->structure);
  free(he);
}

void print_stats(unordered_map<hash_entry, int, hash_fncts> structs)
{
  double mean = 0.0;
  int count = 0;
  double entropy = 0.0;
  unordered_map<hash_entry, int, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    count += it->second;
    mean += (it->first.energy)*(it->second);
    entropy += it->second*log(it->second);
  }

  mean /= (double)count*100.0;
  entropy = entropy/(double)count - log(count);

  fprintf(stderr, "Mean  : %.3f\n"
                  "Entrpy: %.3f\n", mean, entropy);
}

// free hash
void free_hash(unordered_map<hash_entry, int, hash_fncts> &structs)
{
  unordered_map<hash_entry, int, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    free(it->first.structure);
  }
  structs.clear();
}

// free hash
void free_hash(unordered_set<hash_entry*, hash_fncts2, hash_eq> &structs)
{
  unordered_set<hash_entry*, hash_fncts2, hash_eq>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    free((*it)->structure);
    free(*it);
  }
  structs.clear();
}
