#ifndef _hash_util_h
#define _hash_util_h

#include <unordered_map>
#include <unordered_set>

extern "C" {
  #include "utils.h"
}

using namespace std;

// copying of arrays
short* allocopy(short *src);
void copy_arr(short *desc, short *src);


typedef struct _hash_entry {
  short *structure;    /* my structure */
  int energy;       /* my energy */

  bool operator==(const _hash_entry &second) const{
    int i=1;
    while (i<=structure[0] && structure[i]==second.structure[i]) {
      i++;
    }
    if (i>structure[0]) return true;
    else return false;
  }

} hash_entry;

#ifndef HASHBITS
#define HASHBITS 24
#endif
#define HASHSIZE (((unsigned) 1<<HASHBITS)-1)

#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

struct hash_eq {
  bool operator()(const hash_entry *lhs, const hash_entry *rhs) const{
    int i=1;
    while (i<=lhs->structure[0] && lhs->structure[i]==rhs->structure[i]) {
      i++;
    }
    if (i>lhs->structure[0]) return true;
    else return false;
  }
};

struct hash_fncts{
  size_t operator()(const hash_entry &x) const {

  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x.structure;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }
};

struct hash_fncts2{
  size_t operator()(const hash_entry *x) const {

  register short *k;        /* the key */
  register unsigned  length;   /* the length of the key */
  register unsigned  initval=0;  /* the previous hash, or an arbitrary value */
  register unsigned a,b,c,len;

  /* Set up the internal state */
  k = x->structure;
  len = length = (unsigned) k[0];
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] +((unsigned)k[1]<<8) +((unsigned)k[2]<<16) +((unsigned)k[3]<<24));
    b += (k[4] +((unsigned)k[5]<<8) +((unsigned)k[6]<<16) +((unsigned)k[7]<<24));
    c += (k[8] +((unsigned)k[9]<<8) +((unsigned)k[10]<<16)+((unsigned)k[11]<<24));
    mix(a,b,c);
    k += 12; len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch(len) {             /* all the case statements fall through */
    case 11: c+=((unsigned)k[10]<<24);
    case 10: c+=((unsigned)k[9]<<16);
    case 9 : c+=((unsigned)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b+=((unsigned)k[7]<<24);
    case 7 : b+=((unsigned)k[6]<<16);
    case 6 : b+=((unsigned)k[5]<<8);
    case 5 : b+=k[4];
    case 4 : a+=((unsigned)k[3]<<24);
    case 3 : a+=((unsigned)k[2]<<16);
    case 2 : a+=((unsigned)k[1]<<8);
    case 1 : a+=k[0];
     /* case 0: nothing left to add */
  }
  mix(a,b,c);
   /*-------------------------------------------- report the result */
  return (c & HASHSIZE);
  }
};

// print stats about hash
void print_stats(unordered_map<hash_entry, int, hash_fncts> structs);

// free hash
void free_hash(unordered_map<hash_entry, int, hash_fncts> &structs);
void free_hash(unordered_set<hash_entry*, hash_fncts2, hash_eq> &structs);

// entry handling
hash_entry *copy_entry(const hash_entry *he);
void free_entry(hash_entry *he);
#endif
