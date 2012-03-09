#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <map>
#include <vector>
#include <algorithm>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include "move_set.h"

// ############################## DECLARATION #####################################
// private functions & declarations

static int cnt_move = 0;
int count_move() {return cnt_move;}

// compatible base pair?
inline bool compat(char a, char b);

// done with all structures along the way to deepest
int update_deepest(Encoded &Enc, float &deepest, short *pt, degen &deg, bool verbose);

// if the base is lone
inline bool lone_base(short *pt, int i);

// try insert base pair (i,j)
inline bool try_insert(short *pt, const char *seq, int i, int j);

// ############################## IMPLEMENTATION #####################################

// reads a line no matter how long
char* my_getline(FILE *fp)
{
  char s[512], *line, *cp;
  line = NULL;
  do {
    if(fgets(s, 512, fp) == NULL) break;
    cp = strchr(s, '\n');
    if(cp != NULL) *cp = '\0';
    if(line == NULL) line = (char *) calloc(strlen(s) + 1, sizeof(char));
    else line = (char *) realloc(line, strlen(s) + strlen(line) + 1);
    strcat (line, s);
  } while (cp == NULL);
  return (line);
}

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

// done with all structures along the way to deepest
int update_deepest(int &deepest, short *min_pt)
{
  // debug
  /*if (deg.opt->f_point) {
    fprintf(stderr, "UD: %s %.2f (%d, %d) (%d, %d)\n", pt_to_str(Enc.pt).c_str(), deepest/100.0, Enc.bp_left, Enc.bp_right, Enc.bp_left2, Enc.bp_right2);
  }*/

  bool verbose = Opt.verbose_lvl>3;

  int tmp_en;
  tmp_en = Enc.EnergyOfMove();

  //use f_point if we have it; (flooding)
  if (deg.opt->f_point) {
    bool end = deg.opt->f_point(Enc.pt, tmp_en);

    // if we continue - revert changes
    if (!end) {
      undo_moves(Enc);
      return 0;
    } else {
      deepest = tmp_en;
      forget(Enc);
      return 1;
    }
  }

  if (verbose) fprintf(stderr, "  %s %d\n", pt_to_str(Enc.pt).c_str(), tmp_en);
  // better deepest
  if (tmp_en < deepest) {
    deepest = tmp_en;
    copy_arr(min_pt, Enc.pt);

    // delete degeneracy
    Deg.Clear();

    Enc.UndoMove();

    return 1;
  }

  // degeneracy
  if ((abs(deepest - deg.current) <= 0) && (abs(tmp_en - deg.current) <= 0)) {
    if (deg.processed.count(Enc.pt)==0 && deg.unprocessed.count(Enc.pt)==0) {
      deg.unprocessed.insert(allocopy(Enc.pt));
    }
  }
  Enc.UndoMove();
  return 0;
}


// deletions move set
int deletions(Encoded &Enc, int &deepest, short *min_pt, degen &deg, bool verbose)
{
  int cnt = 0;
  short *pt = Enc.pt;
  int len = pt[0];

  for (int i=1; i<=len; i++) {
    if (pt[i]>pt[pt[i]]) {  // '('
      Enc.bp_left=-i;
      Enc.bp_right=-pt[i];

      //if nolp enabled, make (maybe) 2nd delete
      if (deg.opt->noLP) {
        int lone = -1;
        if (lone_base(pt, i-1)) lone=i-1;
        else if (lone_base(pt, i+1)) lone=i+1;
        else if (lone_base(pt, pt[i]-1)) lone=pt[i]-1;
        else if (lone_base(pt, pt[i]+1)) lone=pt[i]+1;

        // check
        if (lone != -1 && (pt[lone]==0 || pt[pt[lone]]==0)) {
          fprintf(stderr, "WARNING: pt[%d(or %d)]!=\'.\'", lone, pt[lone]);
        }

        if (lone != -1) {
          Enc.bp_left2=-lone-1;
          Enc.bp_right2=-pt[lone]-1;
        }
        if (!lone_base(pt, pt[lone]-1) && !lone_base(pt, pt[lone]+1)) {
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }
      } else {  // nolp not enabled
        cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
        // in case useFirst is on and structure is found, end
        if (deg.opt->first && cnt > 0) return cnt;
      }
    }
  }
  return cnt;
}

// insertions move set
int insertions(Encoded &Enc, int &deepest, short *min_pt, degen &deg, bool verbose)
{
  int cnt = 0;
  short *pt = Enc.pt;
  int len = pt[0];

  for (int i=1; i<=len; i++) {
    if (pt[i]==0) {
      for (int j=i+1; j<=len; j++) {
        // end if found closing bracket
        if (pt[j]!=0 && pt[j]<j) break;  //')'
        if (pt[j]!=0 && pt[j]>j) {       //'('
          j = pt[j];
          continue;
        }
        // if conditions are met, do insert
        if (try_insert(pt, Enc.seq, i, j)) {
          Enc.bp_left=i;
          Enc.bp_right=j;

          if (deg.opt->noLP) {
            // if lone bases occur, try inserting one another base
            if (lone_base(pt, i) || lone_base(pt, j)) {
              // inside
              if (try_insert(pt, Enc.seq, i+1, j-1)) {
                  Enc.bp_left2=i+1;
                  Enc.bp_right2=j-1;
                cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
                // in case useFirst is on and structure is found, end
                if (deg.opt->first && cnt > 0) return cnt;
              } else  //outside
              if (try_insert(pt, Enc.seq, i-1, j+1)) {
                Enc.bp_left2=i-1;
                Enc.bp_right2=j+1;
                cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
                // in case useFirst is on and structure is found, end
                if (deg.opt->first && cnt > 0) return cnt;
              }
            } else {
              cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
              // in case useFirst is on and structure is found, end
              if (deg.opt->first && cnt > 0) return cnt;
            }
          } else {
            cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
            // in case useFirst is on and structure is found, end
            if (deg.opt->first && cnt > 0) return cnt;
          }
        }
      }
    }
  }
  return cnt;
}

//shift move set
int shifts(int &deepest, short *min_pt)
{
  int cnt = 0;
  int brack_num = 0;
  short *pt = Enc.pt;
  int len = pt[0];

  bool verbose = Opt.verbose_lvl>3;

  for (int i=1; i<=len; i++) {
    if (pt[i]!=0 && pt[i]>i) {  //'('
      int j=pt[i];

      // outer switch left
      if (verbose) fprintf(stderr, "%2d bracket %2d position, outer switch left\n", brack_num+1, i);
      for (int k=i-1; k>0; k--) {
        if (pt[k]!=0 && pt[k]>k/*'('*/) break;
        if (pt[k]!=0 && pt[k]<k/*')'*/) {
          k = pt[k];
          continue;
        }
        // checks
        if (pt[k]!=0) {
          fprintf(stderr, "WARNING: \'%c\'should be \'.\' at pos %d!\n", pt[k], k);
        }

        // switch (i,j) to (k,j)
        if (j-k>MINGAP && compat(Enc.seq[k-1], Enc.seq[j-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=k;
          Enc.bp_right2=j;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }

        // switch (i,j) to (k,i)
        if (i-k>MINGAP && compat(Enc.seq[i-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=k;
          Enc.bp_right2=i;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;

        }
      }

      // outer switch right
      if (verbose) fprintf(stderr, "%2d bracket %2d position, outer switch right\n", brack_num+1, i);
      for (int k=j+1; k<=len; k++) {
        if (pt[k]!=0 && pt[k]<k/*')'*/) break;
        if (pt[k]!=0 && pt[k]>k/*'('*/) {
          k = pt[k];
          continue;
        }

        // check
        if (pt[k]!=0) {
          fprintf(stderr, "WARNING: \'%c\'should be \'.\' at pos %d!\n", pt[k], k);
        }
        // switch (i,j) to (i,k)
        if (k-i>MINGAP && compat(Enc.seq[i-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=i;
          Enc.bp_right2=k;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }
        // switch (i,j) to (j,k)
        if (k-j>MINGAP && compat(Enc.seq[j-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=j;
          Enc.bp_right2=k;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }
      }

      if (verbose) fprintf(stderr, "%2d bracket %2d position, inner switch\n", brack_num+1, i);
      // inner switch
      for (int k=i+1; k<j; k++) {
        // jump to end of the sub-bracketing
        if (pt[k]!=0 && pt[k]>k/*'('*/) {
            k=pt[k];
            continue;
        }

        // left switch (i,j) to (k,j)
        if (j-k>MINGAP && compat(Enc.seq[k-1], Enc.seq[j-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=k;
          Enc.bp_right2=j;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }

        // right switch (i,j) to (i,k)
        if (k-i>MINGAP && compat(Enc.seq[i-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=i;
          Enc.bp_right2=k;
          cnt += update_deepest(Enc, deepest, min_pt, deg, verbose);
          // in case useFirst is on and structure is found, end
          if (deg.opt->first && cnt > 0) return cnt;
        }
      } // end inner switch for
      brack_num++;
    } // end if (pt[i]=='(')
  } // end for in switches
  return cnt;
}

// move to deepest (or first) neighbour
int move_set(hash_entry &str)
{
  // count how many times called
  cnt_move++;

  // count better neighbours
  int cnt = 0;

  // deepest descent
  hash_entry min;
  min.structure = allocopy(str.structure);
  min.energy = str.energy;

  if (deg.opt->verbose_lvl>3) fprintf(stderr, "\n  start of MS:\n  %s %d\n\n", pt_to_str(min.structure).c_str(), deepest);

  // if using first dont do all of them
  bool end = false;
  // insertions
  if (!end) cnt += insertions(str, min);
  if (Opt.first && cnt>0) end = true;
  if (Opt.verbose_lvl>3) fprintf(stderr, "\n");

  // deletions
  if (!end) cnt += deletions(str, min);
  if (Opt.first && cnt>0) end = true;

  // shifts (only if enabled + noLP disabled)
  if (!end && Opt.shift && !Opt.noLP) {
    cnt += shifts(str, min);
    if (Opt.first && cnt>0) end = true;
  }

  if (Opt.verbose_lvl>3) fprintf(stderr, "\n  %s\n  %s\n", Enc.seq, pt_to_str(min.structure).c_str());

  // if degeneracy occurs, solve it!
  if (!end && Deg.unprocessed.size()>0) {
    Deg.current = str.energy;  // should be true nevertheless
    Deg.processed.insert(allocopy(str.structure));
    // take first to structure
    if (str.structure) free(str.structure);
    str.structure = (*Deg.unprocessed.begin());
    Deg.unprocessed.erase(Deg.unprocessed.begin());
    cnt += move_set(str);
  } else {
    // write output to str
    copy_arr(str.structure, min.structure);
    str.energy = deepest;
  }
  // release minimal
  free(min.structure);

  // resolve degeneracy in local minima
  if (Deg.processed.size()>0) {
    Deg.processed.insert(str.stucture);
    str.structure = *Deg.processed.begin();
    str.energy = Deg.current; // should be true nevertheless
    Deg.processed.erase(Deg.processed.begin());

    erase_set(Deg.processed);
  }

  return cnt;
}

// browse neighbours without degeneracy (almost the same as move_set) assume Deg.first = true
hash_entry *browse_neighs(hash_entry &str, int &saddle_en)
{
  // count how many times called
  cnt_move++;

  // unused, just for necessity of functions
  short *min_pt = Enc.pt;

  // used
  Deg.current = str.energy;

  // did we just find escape from basin?
  bool escape = false;

  //if (deg.opt->verbose_lvl>3) fprintf(stderr, "\n  browse neighs:\n  %s %d\n\n", pt_to_str(Enc.pt).c_str(), energy);

  escape = insertions(energy, min_pt);
  //if (deg.opt->verbose_lvl>3) fprintf(stderr, "\n");

  // deletions
  if (!escape) escape = deletions(energy, min_pt);

  // shifts (only if enabled + noLP disabled)
  if (Opt.shift && !Opt.noLP) {
    if (!escape) escape = shifts(energy, min_pt);
  }

  if (escape) saddle_en = Deg.current;

  return escape;
}


//check if base is lone
bool lone_base(short *pt, int i)
{
  if (i<=0 || i>pt[0]) return false;
  // is not a base pair
  if (pt[i]==0) return false;

  // base is lone:
  if (i-1>0) {
    // is base pair and is the same bracket
    if (pt[i-1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return false;
  }

  if (i+1<=pt[0]) {
    if (pt[i+1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return false;
  }

  return true;
}

// find the structure's lone pairs
int find_lone_pair(string &str)
{
  for(unsigned int i=0; i<str.length(); i++) {
    if (str[i]=='(') {
      if (i+1==str.length() || str[i+1]!='(') {
        return i;
      } else while (i+1!=str.length() && str[i+1]=='(') i++;
    }

    if (str[i]==')') {
      if (i+1==str.length() || str[i+1]!=')') {
        return i;
      } else while (i+1!=str.length() && str[i+1]==')') i++;
    }
  }

  return -1;
}

// if the structure has lone pairs
int find_lone_pair(short* str)
{
  for(int i=1; i<str[0]; i++) {
    if (str[i]==0) continue; // '.'

    if (str[i]>str[str[i]]) {  // '('
      if (i+1==str[0] || str[i+1]==0 || str[i+1]<str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]>str[str[i+1]]) i++;
    }

    if (str[i]<str[str[i]]) {  // ')'
      if (i+1==str[0] || str[i+1]==0 || str[i+1]>str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]<str[str[i+1]]) i++;
    }
  }

  return -1;
}

// try insert base pair (i,j)
bool try_insert(short *pt, const char *seq, int i, int j)
{
  if (i<=0 || j<=0 || i>pt[0] || j>pt[0]) return false;
  return (j-i>MINGAP && pt[j]==0 && pt[i]==0 && compat(seq[i-1], seq[j-1]));
}


// print rates to a file
void print_rates(char *filename, double temp, int num, float *energy_barr, vector<int> &output_en)
{
  FILE *rates;
  rates = fopen(filename, "w");
  if (rates==NULL) {
    fprintf(stderr, "ERROR: couldn't open file \"%s\" for rates! (using stderr instead)\n", filename);
    rates = stderr;
  }
  double _kT = 0.00198717*(273.15 + temp);
  for (int i=0; i<num; i++) {
    for (int j=0; j<num; j++) {
      float res = 0.0;
      if (i!=j) {
        // Arhenius kinetics (as A method in treekin)
        res = 1.0*exp(-(energy_barr[i*num+j]-(output_en[i]/100.0))/_kT);
      }
      fprintf(rates, "%10.4g ", res);
    }
    fprintf(rates, "\n");
  }
  fclose(rates);
}

// pt to str
string pt_to_str(short *pt)
{
  string str;
  str.resize(pt[0]);
  for (int i=1; i<=pt[0]; i++) {
    if (pt[i]==0) str[i-1]='.';
    else if (pt[i]<i) str[i-1]=')';
    else str[i-1]='(';
  }
  return str;
}
