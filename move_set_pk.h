#ifndef __MOVE_SET_H
#define __MOVE_SET_H

/* used data structure*/
typedef struct _struct_en{
  int energy;        /* energy in 10kcal/mol*/
  short *structure;  /* structure in energy_of_move format*/
} struct_en;

enum PK_TYPE = {S, H, K, L ,M}; // types of pseudoknot - S is PKfree structure

class Bpair {
  int start;  // left parenthesis
  int parent; // enclosed in multiloop starting in?
  int right_most; // rightmost nested sibling
  vector<int> left_cross; // crossings from left side
  vector<int> right_cross; // crossings from right side

  Bpair(Pseudoknot &pknot, int left, int right);
}

class Pseudoknot {

  // data
  PK_TYPE type;
  int energy_penalty;

  vector<Bpair> bpairs;

  // helper
  map<int, int> points;  // points to starts and ends of bps.

  Pseudoknot();
  int AddBpair(int left, int right);
  int RemoveBpair(int left);
}




/* free_se */
void free_se(struct_en *to_free);

/* prints structure*/
void print_stren(FILE *out, struct_en *str);
void print_str(FILE *out, short *str);

/* copying functions*/
void copy_arr(short *dest, const short *src); /*just copy*/
short *allocopy(const short *src);            /*copy and make space*/
void copy_se(struct_en *dest, const struct_en *src); /*just copy*/
struct_en *allocopy_se(const struct_en *input); /*copy and make space*/

enum MOVE_TYPE {GRADIENT, FIRST, ADAPTIVE};

/* walking methods (verbose_lvl 0-2, shifts = use shift moves? noLP = no lone pairs? (not compatible with shifts))
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
    methods:  deepest - lowest energy structure is used
              first - first found lower energy structure is used
              rand - random lower energy structure is used
    returns local minima structure in ptable and its energy in 10kcal/mol as output */

int move_gradient( char *seq,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level,
                  int shifts,
                  int noLP);
int move_first( char *seq,
                short *ptable,
                short *s,
                short *s1,
                int verbosity_level,
                int shifts,
                int noLP);
int move_adaptive(  char *seq,
                short *ptable,
                short *s,
                short *s1,
                int verbosity_level);

/* standardized method that encapsulates above "_pt" methods
  input:  seq - sequence
          struc - structure in dot-bracket notation
          type - type of move selection according to MOVE_TYPE enum
  return: energy of LM
          structure of LM in struc in bracket-dot notation
*/
int move_standard(char *seq,
                  char *struc,
                  enum MOVE_TYPE type,
                  int verbosity_level,
                  int shifts,
                  int noLP);


/* browse_neighbours and perform funct function on each of them (used mainly for user specified flooding)
    input:    seq - sequence
              ptable - structure encoded with make_pair_table() from pair_mat.h (in output the structure set by "funct")
              s, s1 - sequence encoded with encode_sequence from pair_mat.h
              funct - function (structure from the neighbourhood, structure from input) to perform on every structure in neigbourhood (if the function returns non-zero, the iteration through neighbourhood stops.)
                      function takes 2 arguments: first is input - structre in iteration, second is input/output - at first it is the structure from input "ptable" with energy, then it is as structure set with funct on subsequent runs
                        this goes as output of whole function (energy as return value, structure as ptable)
    returns energy of the structure funct sets as the second argument*/
int browse_neighs_pt( char *seq,
                   short *ptable,
                   short *s,
                   short *s1,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*));

int browse_neighs( char *seq,
                   char *struc,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*));

#endif



