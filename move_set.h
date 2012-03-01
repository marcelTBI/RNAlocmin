#ifndef __MOVE_SET_H
#define __MOVE_SET_H

#include <string>
#include <set>

using namespace std;

// minimal gap for loop
#define MINGAP 3

typedef struct _encoded {
  // for energy_of_move
  short *pt;    // structure
  short *s0;
  short *s1;

  const char  *seq;

  int   bp_left;
  int   bp_right;
  int   bp_left2;   // if noLP is enabled (and for shift moves)
  int   bp_right2;
} encoded;

// comparator for structures (in notation as used in energy_of_move)
struct setcomp {
  bool operator() (const short *lhs, const short *rhs) const {
    // here we have structures as numbers, but we want to compare them as chars in bractet dot notation: "()."
    int i=1;
    char l,r;
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

// cute nice options structure
typedef struct _options {
  // options
  float minh;
  bool noLP;    // no lone pairs
  bool EOM;     // use energy_of_move
  bool first;   // use first descent, not deepest
  bool shift;   // use shifts?
  int verbose_lvl; // level of verbosity

  // pointer to function used on every neighbour (in update deepest)
  bool (*f_point) (short *, int);  // (I dont like it either, but it was easiest way to program it...:/ )

} options;

// structure for degeneracy + options management (almost everything to avoid using global variables - maybe bad idea... but i really hate them)
typedef struct _degen {
  // for degeneracy (structures with equal energies)
  float current;
  set<short*, setcomp> processed;
  set<short*, setcomp> unprocessed;

  options *opt;

} degen;

//#################################################### SMALL (HELPER) FUNCTIONS

// count how many times move_set has been called
int count_move();

// create encoded structure (user must free the structure by calling free_encode() once after each encode_seq)
encoded *encode_seq(const char* seq);
void encode_str(encoded *enc, const char *str);
void free_encode(encoded *enc);

// abs(float)
inline float abs(float a) { return (a<0?-a:a); }

// reads a line no matter how long
char* my_getline(FILE *fp);

// pt to str
string pt_to_str(short *pt);

// erase set (free memory)
void erase_set(set<short*, setcomp> &_set);

// if the structure has lone pairs
int find_lone_pair(string &str);

// if the structure has lone pairs
int find_lone_pair(short* str);


// ################################################### BIG (IMPORTANT) FUNCTIONS

// move to deepest neighbour (returns minimum in enc, deepest) returns 0 when minimum has been found
int move_set(encoded &enc, int &deepest, degen &deg);

// browse neighbours and apply funct on all of them (returns minimum in enc, energy) return saddle energy in saddle_en (in flooding) return if funct(str) returns true
bool browse_neighs(encoded &enc, int &energy, degen &deg, int &saddle_en);

// print rates to a file
void print_rates(char *filename, double temp, int num, float *energy_barr, vector<int> &output_en);


#endif
