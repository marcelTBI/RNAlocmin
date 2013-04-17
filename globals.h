#ifndef __GLOBALS_H
#define __GLOBALS_H

#include <set>
#include <vector>

using namespace std;

extern "C" {
  #include "RNAlocmin_cmdline.h"
}
#include "hash_util.h"

// minimal gap for loop
#define MINGAP 3

// some global counters
static int num_moves = 0;
static int seq_len;

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

  // all possible moves
  vector<int> moves_from;
  vector<int> moves_to;

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

  // permute possible moves
  void Permute();

  void PossMoves(hash_entry &str);
};

// cute options singleton class
class Options {
  // options
public:
  float minh;   // leave out shallow minima (should be relativelly small)
  bool noLP;    // no lone pairs
  bool EOM;     // use energy_of_move
  bool first;   // use first descent, not deepest
  bool rand;    // use random walk, not deepest
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

// structure for degeneracy
class Degen {
  // for degeneracy (structures with equal energies)
public:
  int current;    // all structures here have this energy
  set<short*, comps_short> processed;
  set<short*, comps_short> unprocessed;

public:
  Degen();
  ~Degen();
  void Clear();
};

// some good functions

  // compatible base pair?
inline bool compat(char a, char b) {
  if (a=='A' && b=='U') return true;
  if (a=='C' && b=='G') return true;
  if (a=='G' && b=='U') return true;
  if (a=='U' && b=='A') return true;
  if (a=='G' && b=='C') return true;
  if (a=='U' && b=='G') return true;
  // and with T's
  if (a=='A' && b=='T') return true;
  if (a=='T' && b=='A') return true;
  if (a=='G' && b=='T') return true;
  if (a=='T' && b=='G') return true;
  return false;
}

// try insert base pair (i,j)
inline bool try_insert(const short *pt, const char *seq, int i, int j)
{
  if (i<=0 || j<=0 || i>pt[0] || j>pt[0]) return false;
  return (j-i>MINGAP && pt[j]==0 && pt[i]==0 && compat(seq[i-1], seq[j-1]));
}

// try insert base pair (i,j)
inline bool try_insert(const char *seq, int i, int j)
{
  if (i<=0 || j<=0) return false;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]));
}

// some singleton objects
extern Degen Deg;
extern Options Opt;
extern Encoded Enc;

#endif
