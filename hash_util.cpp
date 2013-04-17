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

// does not depend on number - just on dot-bracket notation:
bool compf_short (const short *lhs, const short *rhs) {
  int i=1;
  char l,r;
  while (i<=lhs[0]) {
    l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?')':'('));
    r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?')':'('));
    if (l != r) {
      //fprintf (stderr, "%c %c %d %d\n", l, r, l, r);
      break;
    }
    i++;
  }
  return (i<=lhs[0] && l<r);
}

bool compf_entries (const hash_entry *lhs, const hash_entry *rhs) {
  if (lhs->energy!=rhs->energy) return lhs->energy<rhs->energy;
  return compf_short(lhs->structure, rhs->structure);
}
bool compf_entries2 (const hash_entry &lhs, const hash_entry &rhs)
{
  if (lhs.energy!=rhs.energy) return lhs.energy<rhs.energy;
  return compf_short(lhs.structure, rhs.structure);
}

bool compf_short_rev (const short *lhs, const short *rhs) {
  int i=1;
  char l,r;
  while (i<=lhs[0]) {
    l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?')':'('));
    r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?')':'('));
    if (l != r) break;
    i++;
  }
  return (i<=lhs[0] && l>r);
}

bool compf_entries_rev (const hash_entry *lhs, const hash_entry *rhs) {
  if (lhs->energy!=rhs->energy) return lhs->energy>rhs->energy;
  return compf_short_rev(lhs->structure, rhs->structure);
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

void print_stats(unordered_map<hash_entry, gw_struct, hash_fncts> &structs)
{
  double mean = 0.0;
  int count = 0;
  double entropy = 0.0;
  unordered_map<hash_entry, gw_struct, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    count += it->second.count;
    mean += (it->first.energy)*(it->second.count);
    entropy += it->second.count*log(it->second.count);
  }

  mean /= (double)count*100.0;
  entropy = entropy/(double)count - log(count);

  fprintf(stderr, "Mean  : %.3f (Entrpy: %.3f)\n", mean, entropy);
}
//#include "move_set.h"

void add_stats(unordered_map<hash_entry, gw_struct, hash_fncts> &structs, map<hash_entry, int, comps_entries> &output)
{
  unordered_map<hash_entry, gw_struct, hash_fncts>::iterator it;
  for (it=structs.begin(); it!=structs.end(); it++) {
    // add stats:
    //fprintf(stderr, "struct: %s %6.2f %d\n", pt_to_str(it->second.he.structure).c_str(), it->second.he.energy/100.0, it->second.count);

    if (output.count(it->second.he) == 0) {
      fprintf(stderr, "ERROR: output does not contain structure it should!!!\n");
      exit(EXIT_FAILURE);
    }
    output[it->second.he] += it->second.count-1;
  }
}

// free hash
void free_hash(unordered_map<hash_entry, gw_struct, hash_fncts> &structs)
{
  unordered_map<hash_entry, gw_struct, hash_fncts>::iterator it;
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
