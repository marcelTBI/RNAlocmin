#ifndef __FLOOD_H
#define __FLOOD_H

#include "hash_util.h"
#include "move_set.h"


// hash_entry comparator function
bool compare_vect (const hash_entry &lhs, const hash_entry &rhs);

// flood the structure - return one below saddle structure (should be freed then) energy of saddle is in "saddle_en"
  // maxh = height of flood (0 = infinity)
  // if returns NULL - in saddle_en is fail status - 1 for maxh reached, 0 otherwise
hash_entry* flood(const hash_entry &str, int &saddle_en, int maxh = 0);

#endif
