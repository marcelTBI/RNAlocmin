/* stolen from Barriers 1.5.2 */
/* hash_util.h */

#ifndef _hash_util_h
#define _hash_util_h

extern void * lookup_hash (void *x);
extern int write_hash (void *x);
extern void delete_hash (void *x);
extern void kill_hash();
extern void initialize_hash();
extern unsigned int hash_size();

// print stats about current hash
extern void print_stats();

typedef struct _hash_entry {
  short *structure;    /* my structure */
  int energy;       /* my energy */

  int count; // for counting in descend (not used in flooding)
} hash_entry;

#endif

/* End of file */
