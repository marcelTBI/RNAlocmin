#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <queue>

extern "C" {
  #include "hash_util.h"
  #include "utils.h"
}

#include "flood.h"

extern short* allocopy(short *src);
extern void copy_arr(short *dest, short *src);

// global priority queue for stuff in flooding (does not hold memory, thus should not be freed (thats done in hash))
priority_queue<hash_entry*, vector<hash_entry*>, compare_queue> neighs;
int energy_lvl;
bool debugg;

bool compare_vect (const hash_entry &lhs, const hash_entry &rhs){
  // first energies
  if (lhs.energy != rhs.energy) {
    return lhs.energy<rhs.energy;
  }
  // then structures (here we have structures as numbers, but we want to compare them as chars in bractet dot notation: "()." )
  int i=1;
  char l=0,r=0;
  while (i<=lhs.structure[0]) {
    l = (lhs.structure[i]==0?'.':(lhs.structure[i]<lhs.structure[lhs.structure[i]]?'(':')'));
    r = (rhs.structure[i]==0?'.':(rhs.structure[i]<rhs.structure[rhs.structure[i]]?'(':')'));
    if (l != r) break;
    i++;
  }
  return (i<=lhs.structure[0] && l<r);
}

// function to do on all the items...
bool flood_func(hash_entry &input)
{
  // have we seen him?
  hash_entry *seen =  (hash_entry*)lookup_hash(&input);
  if (seen != NULL) {
    // nothing to do with already processed structure
    if (debugg) fprintf(stderr, "       already seen: %s %.2f\n", pt_to_str(struc).c_str(), energy/100.0);
    return false;
  } else {
    // found escape? (its energy is lower than our energy lvl and we havent seen it)
    if (input->energy < energy_lvl) {
      // ends flood and return it as a structure to walk down
      if (debugg) fprintf(stderr, "       escape  : %s %.2f\n", pt_to_str(struc).c_str(), energy/100.0);
      return true;
    } else {
      if (debugg) fprintf(stderr, "       adding  : %s %.2f\n", pt_to_str(struc).c_str(), energy/100.0);
      // just add it to the queue... and to hash
      hash_entry *hee = copy_entry(input);
      neighs.push(hee);
      write_hash(hee); // create hash_entry to write (hash stores pointers -> should be emptied)
      return false;
    }
  }
}

hash_entry* flood(hash_entry &he, int &saddle_en)
{
  int count = 0;
  debugg = Opt.verbose_lvl>3;

  // init priority queue
  while (!neighs.empty()) {
    //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
    //free(neighs.top().structure);
    neighs.pop();
  }

  // init hash
  initialize_hash();

  // add the first structure
  hash_entry *bottom = copy_entry(he);
  neighs.push(bottom);
    // also to hash
  write_hash(bottom);

  // save deg.first + create deg
  bool first = Opt.first;
  Opt.first = true;
  Opt.f_point = flood_func;

  // returning hash_entry
  hash_entry *escape = NULL;

  // FLOOOD!
  while ((int)hash_size() < Opt.floodMax) {
    // should not be empty
    if (neighs.empty()) break;

    // get structure (pop it in the end)
    hash_entry *he_tmp = neighs.top();
    neighs.pop();
    energy_lvl = he_tmp->energy;

    if (Opt.verbose_lvl>2) fprintf(stderr, "  neighbours of: %s %.2f\n", pt_to_str(Enc.pt).c_str(), he_tmp->energy/100.0);

    escape = browse_neighs(*he_tmp, saddle_en);

    if (Opt.verbose_lvl>2) fprintf(stderr, "sad= %6.2f (%c): %s %.2f\n", saddle_en/100.0, (escape?'t':'f'), pt_to_str(Enc.pt).c_str(), Enc.energy/100.0);

    // did we find exit from basin?
    if (escape) {
      break;
    }

    count++;
  }

  // restore deg options
  Opt.first = first;
  Opt.f_point = NULL;

  // destroy queue
  while (!neighs.empty()) {
    neighs.pop();
  }

  // destroy hash
  kill_hash();

  return escape;
}


