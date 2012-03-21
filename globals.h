#ifndef __GLOBALS_H
#define __GLOBALS_H

#include <set>

using namespace std;

extern "C" {
  #include "RNAlocmin_cmdline.h"
}
#include "hash_util.h"

// minimal gap for loop
#define MINGAP 3

// for energy_of_move
class Encoded {

public:
  //short *pt;    // structure
  short *s0;
  short *s1;

  char  *seq;

  // moves
  int   bp_left;
  int   bp_right;
  int   bp_left2;   // if noLP is enabled (and for shift moves)
  int   bp_right2;

  // last energy
  int last_en;

public:
  Encoded();
  ~Encoded();

  void Init(const char* seq);
  short *Struct(const char *str);  // should be freed!!!

  void Forget();
  inline void Move(hash_entry &he, bool first = true, bool second = true);
  void UndoMove(hash_entry &he, bool first = true, bool second = true);

  // energy calculations on structures
  int Energy(hash_entry &he);
  int EnergyOfMove(hash_entry &he);
};

// cute nice options singleton class
class Options {
  // options
public:
  float minh;   // leave out shallow minima (should be relativelly small)
  bool noLP;    // no lone pairs
  bool EOM;     // use energy_of_move
  bool first;   // use first descent, not deepest
  bool shift;   // use shifts?
  int verbose_lvl; // level of verbosity
  int floodMax; // cap for flooding

  // pointer to function used on every neighbour (in update deepest)
  bool (*f_point) (hash_entry &);  // (I dont like it either, but it was easiest way to program it...:/ )

public:
  Options();

  // return 0 if success
  int Init(gengetopt_args_info &args_info);
};

// comparator for structures (in notation as used in energy_of_move)
struct setcomp {
  bool operator() (const short *lhs, const short *rhs) const {
    // here we have structures as numbers, but we want to compare them as chars in bractet dot notation: "()."
    int i=1;
    char l=0,r=0;
    while (i<=lhs[0]) {
      l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?'(':')'));
      r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?'(':')'));
      if (l != r) break;
      i++;
    }
    return (i<=lhs[0] && l<r);

    /*
    int i=1;
    while (i<=lhs[0] && lhs[i]==rhs[i]) {
      i++;
    }
    return (i<=lhs[0] && lhs[i]>rhs[i]);*/
  }
};

// structure for degeneracy
class Degen {
  // for degeneracy (structures with equal energies)
public:
  int current;    // all structures here have this energy
  set<short*, setcomp> processed;
  set<short*, setcomp> unprocessed;

public:
  Degen();
  ~Degen();
  void Clear();
};

// some singleton objects
extern Degen Deg;
extern Options Opt;
extern Encoded Enc;

#endif
