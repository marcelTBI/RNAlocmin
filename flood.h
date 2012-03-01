#ifndef __FLOOD_H
#define __FLOOD_H

#include "hash_util.h"
#include "move_set.h"

// comparators for hash_lists (structure) (map and queue use different :-/ )
struct compare_map {
  bool operator() (const hash_entry &lhs, const hash_entry &rhs) const {
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
    return (i<=lhs.structure[0] && l<r); /// isnt this wrong??? i>lhs
  }
};

struct compare_queue {
  bool operator() (const hash_entry *lhs, const hash_entry *rhs) const {
    // first energies
    if (lhs->energy != rhs->energy) {
      return lhs->energy>rhs->energy;
    }
    // then structures (here we have structures as numbers, but we want to compare them as chars in bractet dot notation: "()." )
    int i=1;
    char l=0,r=0;
    while (i<=lhs->structure[0]) {
      l = (lhs->structure[i]==0?'.':(lhs->structure[i]<lhs->structure[lhs->structure[i]]?'(':')'));
      r = (rhs->structure[i]==0?'.':(rhs->structure[i]<rhs->structure[rhs->structure[i]]?'(':')'));
      if (l != r) break;
      i++;
    }
    return (i<=lhs->structure[0] && l>r);
  }
};

// hash_entry comparator function
bool compare_vect (const hash_entry &lhs, const hash_entry &rhs);

// flood the structure - return one below saddle structure (should be freed then) energy of saddle is in "saddle_en"
hash_entry* flood(encoded &enc, int energy, options &opt, int &saddle_en);

#endif
