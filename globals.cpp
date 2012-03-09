#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include "globals.h"

// get global degen structure
static Degen *get_degen()
{
  return &Deg;
}

// get global options
static Options *get_opt()
{
  return &Opt;
}

// get current encoded
static Encoded *get_enc()
{
  return &Enc;
}

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
  this->seq = (char*) space(strlen(seq)*sizeof(char));
  strcpy(this->seq, seq);

  s0   = encode_sequence(seq, 0);
  s1   = encode_sequence(seq, 1);
}

void Encoded::Struct(const char *str)
{
  if (pt) free(pt);
  pt = make_pair_table(str);
}

int Encoded::Energy(hash_entry &he)
{
  return energy_of_structure_pt(seq, he.pt, s0, s1, 0);
}

int Encoded::EnergyOfMove(hash_entry &he)
{
  int tmp_en;
  if (Opt.EOM) {
    tmp_en = he.energy + energy_of_move_pt(he.pt, s0, s1, bp_left, bp_right);
    Move(he, true, false);
    if (bp_left2 != 0) {
      tmp_en += energy_of_move_pt(he.pt, s0, s1, bp_left2, bp_right2);
      Move(he, false, true);
    }
  } else {
    Move(he);
    tmp_en = energy_of_structure_pt(seq, he.pt, s0, s1, 0);
  }

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

inline void Encoded::Move(hash_entry &he, bool first = true, bool second = true)
{
  if (first && bp_left != 0) {
    do_move(he.pt, bp_left, bp_right);
  }
  if (second && bp_left2 != 0) {
    do_move(he.pt, bp_left2, bp_right2);
  }
}

void Encoded::UndoMove(hash_entry &he, bool first = true, bool second = true)
{
  if (second && bp_left2 != 0) {
    do_move(he.pt, -bp_left2, -bp_right2);
  }
  if (first && bp_left != 0) {
    do_move(he.pt, -bp_left, -bp_right);
  }

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
/*option "move"               m "Move set:\nI ==> insertion & deletion of base pair\nS ==> I&D& switch base pair" values="I","S" default="I" no
option "min-num"            n "Maximal number of local minima returned" int default="100" no
option "find-num"           - "Maximal number of local minima found \n  (default = unlimited - crawl through whole file)" int no
option "seq"                s "Sequence file in FASTA format" string default="seq.txt"
option "verbose-lvl"        v "Level of verbosity (0 = nothing, 3 = full)" int default="0" no
option "rates"              r "Create rates for treekin" flag off
option "rates-file"         f "File where to write rates" string default="rates.out" no
option "temp"               T "Temperature in Celsius (only for rates)" double default="37.0" no
option "depth"              d "Depth of findpath search (higher values increase running time)" int default="10" no
option "minh"               - "Print only minima with energy barrier greater than this" double default="0.0" no
option "noLP"               - "Work with canonical RNA structures (w/o isolated base pairs)" flag off
option "bartree"            b "Generate possible barrier tree" flag off
option "useEOS"             e "Use energy_of_structure_pt calculation instead of energy_of_move (slower, it should not affect results)" flag off
option "useFirst"           - "Use first found lower energy structure instead of deepest" flag off
option "floodPortion"       - "Fraction of minima to flood\n(0.0 -> no flood; 1.0 -> try to flood all of them)" double default="0.95" no
option "floodMax"           - "Flood cap - how many structures to flood in one basin" int default="1000" no
*/

  int ret = 0;

  if (args_info.min_num_arg<=0) {
    fprintf(stderr, "Number of local minima should be positive integer (min-num)\\n");
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

  if (ret ==-1) return -1;

  // adjust options
  minh = args_info.minh_arg;
  noLP = args_info.noLP_flag;
  EOM = !args_info.useEOS_flag;
  first = args_info.useFirst_flag;
  f_point = NULL;
  shift = args_info.move_arg[0]=='S';
  verbose_lvl = args_info.verbose_lvl_arg;
  floodMax = args_info.floodMax_arg;

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

void copy_arr(short *dest, short *src)
{
  if (!src || !dest) {
    fprintf(stderr, "Empty pointer in copying\n");
    return;
  }
  memcpy(dest, src, sizeof(short)*(src[0]+1));
}

short *allocopy(short *src)
{
  short *res = (short*) space(sizeof(short)*(src[0]+1));
  copy_arr(res, src);
  return res;
}
