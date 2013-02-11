#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <queue>

extern "C" {
  #include "utils.h"
}

#include "flood.h"

// global priority queue for stuff in flooding (does not hold memory - memory is in hash)
priority_queue<hash_entry*, vector<hash_entry*>, comps_entries_rev> neighs;
int energy_lvl;
bool debugg;
int top_lvl;
// hash for the flooding
unordered_set<hash_entry*, hash_fncts2, hash_eq> hash_flood (HASHSIZE);
unordered_set<hash_entry*, hash_fncts2, hash_eq>::iterator it_hash;


// function to do on all the items...
bool flood_func(hash_entry &input)
{
  // have we seen him?
  it_hash = hash_flood.find(&input);
  if (it_hash != hash_flood.end()) {
    // nothing to do with already processed structure
    if (debugg) fprintf(stderr,     "   already seen: %s %.2f\n", pt_to_str(input.structure).c_str(), input.energy/100.0);
    return false;
  } else {
    // found escape? (its energy is lower than our energy lvl and we havent seen it)
    if (input.energy < energy_lvl) {
      // ends flood and return it as a structure to walk down
      if (debugg) fprintf(stderr,   "       escape  : %s %.2f\n", pt_to_str(input.structure).c_str(), input.energy/100.0);
      return true;
    } else {
      if (input.energy > top_lvl) {
        if (debugg) fprintf(stderr, "energy too high: %s %.2f\n", pt_to_str(input.structure).c_str(), input.energy/100.0);
        return false;
      } else {
        if (debugg) fprintf(stderr, "       adding  : %s %.2f\n", pt_to_str(input.structure).c_str(), input.energy/100.0);
        // just add it to the queue... and to hash
        hash_entry *he_tmp = (hash_entry*)space(sizeof(hash_entry));
        he_tmp->structure = allocopy(input.structure);
        he_tmp->energy = input.energy;
        neighs.push(he_tmp);
        hash_flood.insert(he_tmp);
        return false;
      }
    }
  }
}

hash_entry* flood(const hash_entry &he, int &saddle_en, int maxh)
{
  int count = 0;
  debugg = Opt.verbose_lvl>2;

  // init priority queue
  while (!neighs.empty()) {
    //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
    neighs.pop();
  }

  // init hash
  free_hash(hash_flood);

  // add the first structure to hash, get its adress and add it to priority queue
  {
    hash_entry *he_tmp = (hash_entry*)space(sizeof(hash_entry));
    he_tmp->structure = allocopy(he.structure);
    he_tmp->energy = he.energy;
    neighs.push(he_tmp);
    hash_flood.insert(he_tmp);
  }

  // save deg.first + create deg
  bool first = Opt.first;
  Opt.first = true;
  Opt.f_point = flood_func;

  // returning hash_entry
  hash_entry *escape = NULL;

  // if minh specified, assign top_lvl
  if (maxh>0) {
    top_lvl = he.energy + maxh;
  } else {
    top_lvl = 1e9;
  }

  // FLOOOD!
  while ((int)hash_flood.size() < Opt.floodMax) {
    // should not be empty (only when maxh specified)
    if (neighs.empty()) break;

    // get structure
    hash_entry *he_top = neighs.top();
    neighs.pop();
    energy_lvl = he_top->energy;

    if (Opt.verbose_lvl>2) fprintf(stderr, "  neighbours of: %s %.2f\n", pt_to_str(he_top->structure).c_str(), he_top->energy/100.0);

    hash_entry he_browse = *he_top;   // maybe not necessary...
    he_browse.structure = allocopy(he_top->structure);
    escape = browse_neighs(he_browse, saddle_en);
    free(he_browse.structure);

    if (escape && Opt.verbose_lvl>2) fprintf(stderr, "sad= %6.2f    : %s %.2f\n", saddle_en/100.0, pt_to_str(escape->structure).c_str(), escape->energy/100.0);

    // did we find exit from basin?
    if (escape) {
      break;
    }

    count++;
  }

  // restore deg options
  Opt.first = first;
  Opt.f_point = NULL;

  // return status in saddle_en :/
  if (!escape) {
    saddle_en = (neighs.empty() ? 1 : 0);
  }

  // destroy queue
  while (!neighs.empty()) {
    //fprintf(stderr, "-neighs size: %d\n", (int)neighs.size());
    neighs.pop();
  }

  // destroy hash
  free_hash(hash_flood);

  return escape;
}


