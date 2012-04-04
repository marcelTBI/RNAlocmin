#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include "globals.h"

// some singleton objects
Degen Deg;
Options Opt;
Encoded Enc;

// create encoded structure
Encoded::Encoded()
{
  // inicialize encode_sequence
  make_pair_matrix();

  bp_left = 0;
  bp_right = 0;
  bp_left2 = 0;
  bp_right2 = 0;

  seq = NULL;
}

Encoded::~Encoded()
{
  if (s0) free(s0);
  if (s1) free(s1);
  if (seq) free(seq);
}

void Encoded::Init(const char *seq)
{
  if (this->seq) free(this->seq);
  this->seq = (char*) space((strlen(seq)+1)*sizeof(char));
  strcpy(this->seq, seq);

  srand(time(NULL));

  s0   = encode_sequence(seq, 0);
  s1   = encode_sequence(seq, 1);
}

void Encoded::PossMoves(hash_entry &str)
{
  moves_to.clear();
  moves_from.clear();
  moves_to.reserve(2*str.structure[0]);
  moves_from.reserve(2*str.structure[0]);

  // generate all possible moves (less than n^2)
  for (int i=1; i<=str.structure[0]; i++) {
    if (str.structure[i]!=0) {
      if (str.structure[i]<i) continue;
      moves_from.push_back(-i);
      moves_to.push_back(-str.structure[i]);
      //fprintf(stderr, "add  d(%d, %d)\n", i, str.structure[i]);
    } else {
      for (int j=i+1; j<=str.structure[0]; j++) {
        //fprintf(stderr, "check (%d, %d)\n", i, j);
        if (str.structure[j]==0) {
          if (try_insert(seq,i,j)) {
            moves_from.push_back(i);
            moves_to.push_back(j);
            //fprintf(stderr, "add  i(%d, %d)\n", i, j);
            continue;
          }
        } else if (str.structure[j]>j) { // '('
          j = str.structure[j];
        } else break;
      }
    }
  }
}

void Encoded::Permute()
{
    if (Opt.verbose_lvl>2) {
    for (unsigned int i=0; i<moves_from.size(); i++) {
      fprintf(stderr, "(%d %d) ", moves_from[i], moves_to[i]);
    }
    fprintf(stderr, "\n");
  }
  for (unsigned int i=0; i<moves_from.size()-1; i++) {
    int rnd = rand();
    rnd = rnd % (moves_from.size()-i) + i;
    swap(moves_from[i], moves_from[rnd]);
    swap(moves_to[i], moves_to[rnd]);
  }


}

short *Encoded::Struct(const char *str)
{
  return make_pair_table(str);
}

int Encoded::Energy(hash_entry &he)
{
  return energy_of_structure_pt(seq, he.structure, s0, s1, 0);
}

int Encoded::EnergyOfMove(hash_entry &he)
{
  int tmp_en;
  if (Opt.EOM) {
    tmp_en = he.energy + energy_of_move_pt(he.structure, s0, s1, bp_left, bp_right);
    Move(he, true, false);
    if (bp_left2 != 0) {
      tmp_en += energy_of_move_pt(he.structure, s0, s1, bp_left2, bp_right2);
      Move(he, false, true);
    }
  } else {
    Move(he);
    tmp_en = energy_of_structure_pt(seq, he.structure, s0, s1, 0);
  }
  last_en = he.energy;
  he.energy = tmp_en;
  return tmp_en;
}

inline void do_move(short *pt, int bp_left, int bp_right)
{
  // delete
  if (bp_left<0) {
    pt[-bp_left]=0;
    pt[-bp_right]=0;
  } else { // insert
    pt[bp_left]=bp_right;
    pt[bp_right]=bp_left;
  }
}

inline void Encoded::Move(hash_entry &he, bool first, bool second)
{
  if (first && bp_left != 0) {
    do_move(he.structure, bp_left, bp_right);
  }
  if (second && bp_left2 != 0) {
    do_move(he.structure, bp_left2, bp_right2);
  }
}

void Encoded::UndoMove(hash_entry &he, bool first, bool second)
{
  if (second && bp_left2 != 0) {
    do_move(he.structure, -bp_left2, -bp_right2);
  }
  if (first && bp_left != 0) {
    do_move(he.structure, -bp_left, -bp_right);
  }

  he.energy = last_en;
  last_en = 1e5;

  Forget();
}

void Encoded::Forget()
{
  bp_left=0;
  bp_right=0;
  bp_left2=0;
  bp_right2=0;
}

Options::Options()
{
}

int Options::Init(gengetopt_args_info &args_info)
{
  int ret = 0;

  if (args_info.min_num_arg<0) {
    fprintf(stderr, "Number of local minima should be non-negative integer (min-num)\n");
    ret = -1;
  }

  if (args_info.find_num_given && args_info.find_num_arg<=0) {
    fprintf(stderr, "Number of local minima should be positive integer (find-num)\n");
    ret = -1;
  }

  if (args_info.verbose_lvl_arg<0 || args_info.verbose_lvl_arg>4) {
    if (args_info.verbose_lvl_arg<0) args_info.verbose_lvl_arg = 0;
    else args_info.verbose_lvl_arg = 4;
    fprintf(stderr, "WARNING: level of verbosity is not in range (0-4), setting it to %d\n", args_info.verbose_lvl_arg);
  }

  if (args_info.temp_arg<-273.15) {
    fprintf(stderr, "Temperature cannot be below absolute zero\n");
    ret = -1;
  }

  if (args_info.floodMax_arg<0) {
    fprintf(stderr, "Flood cap must be non-negative\n");
    ret = -1;
  }

  if (args_info.floodPortion_arg<0.0 || args_info.floodPortion_arg>1.0) {
    args_info.floodPortion_arg = (args_info.floodPortion_arg<0.0 ? 0.0 : 1.0);
    fprintf(stderr, "WARNING: floodPortion is not in range (0.0-1.0), setting it to %.1f\n", args_info.floodPortion_arg);
  }

  if (args_info.depth_arg<=0) {
    fprintf(stderr, "Depth of findpath search should be positive integer\n");
    ret = -1;
  }

  if (args_info.minh_arg<0.0) {
    fprintf(stderr, "Depth of findpath search should be non-negative number\n");
    ret = -1;
  }

  if (args_info.numIntervals_arg<0) {
    fprintf(stderr, "Number of intervals should be non-negative number\n");
    ret = -1;
  }

  if (ret ==-1) return -1;

  // adjust options
  minh = (int)(args_info.minh_arg*100);
  minhall = args_info.minhall_flag;
  noLP = args_info.noLP_flag;
  EOM = !args_info.useEOS_flag;
  first = args_info.walk_arg[0]=='F';
  rand = args_info.walk_arg[0]=='R';
  f_point = NULL;
  shift = args_info.move_arg[0]=='S';
  verbose_lvl = args_info.verbose_lvl_arg;
  floodMax = args_info.floodMax_arg;
  rna2d = args_info.rna2D_flag;

  return ret;
}

Degen::Degen()
{
  current = 0;
}

Degen::~Degen()
{
  Clear();
}

void Degen::Clear()
{
  // unprocessed
  for (set<short*, setcomp>::iterator it=unprocessed.begin(); it!=unprocessed.end(); it++) {
    if (*it) free(*it);
  }
  unprocessed.clear();

  // processed
  for (set<short*, setcomp>::iterator it=processed.begin(); it!=processed.end(); it++) {
    if (*it) free(*it);
  }
  processed.clear();
}
