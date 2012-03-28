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

// done with all structures along the way to deepest
int update_deepest(hash_entry &str, hash_entry &min);

// if the base is lone
inline bool lone_base(short *pt, int i);

// can move be done on this structure? (move is in the Enc)
bool check_insert(hash_entry &str, int i, int j);

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

// done with all structures along the way to deepest
int update_deepest(hash_entry &str, hash_entry &min)
{
  // debug
  /*if (deg.opt->f_point) {
    fprintf(stderr, "UD: %s %.2f (%d, %d) (%d, %d)\n", pt_to_str(Enc.pt).c_str(), deepest/100.0, Enc.bp_left, Enc.bp_right, Enc.bp_left2, Enc.bp_right2);
  }*/

  bool verbose = Opt.verbose_lvl>3;

  int tmp_en;
  tmp_en = Enc.EnergyOfMove(str);

  //use f_point if we have it; (flooding)
  if (Opt.f_point) {
    bool end = Opt.f_point(str);

    // if we end - remember one below in min
    if (end) {
      min.energy = tmp_en;
      copy_arr(min.structure, str.structure);
    }
    Enc.UndoMove(str);
    return (end?1:0);
  }

  if (verbose) fprintf(stderr, "  %s %d\n", pt_to_str(str.structure).c_str(), tmp_en);
  // better deepest
  if (tmp_en < min.energy) {
    min.energy = tmp_en;
    copy_arr(min.structure, str.structure);

    // delete degeneracy
    Deg.Clear();

    Enc.UndoMove(str);
    return 1;
  }

  // degeneracy
  if ((str.energy == min.energy) && (Deg.current == min.energy)) {
    if (Deg.processed.count(str.structure)==0 && Deg.unprocessed.count(str.structure)==0) {
      Deg.unprocessed.insert(allocopy(str.structure));
    }
  }
  Enc.UndoMove(str);
  return 0;
}


// deletions move set
int deletions(hash_entry &str, hash_entry &minim)
{
  int cnt = 0;
  short *pt = str.structure;
  int len = pt[0];

  for (int i=1; i<=len; i++) {
    if (pt[i]>pt[pt[i]]) {  // '('
      Enc.bp_left=-i;
      Enc.bp_right=-pt[i];

      //if nolp enabled, make (maybe) 2nd delete
      if (Opt.noLP) {
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
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
        }
      } else {  // nolp not enabled
        cnt += update_deepest(str, minim);
        // in case useFirst is on and structure is found, end
        if (Opt.first && cnt > 0) return cnt;
      }
    }
  }
  return cnt;
}

// insertions move set
int insertions(hash_entry &str, hash_entry &minim)
{
  int cnt = 0;
  short *pt = str.structure;
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

          if (Opt.noLP) {
            // if lone bases occur, try inserting one another base
            if (lone_base(pt, i) || lone_base(pt, j)) {
              // inside
              if (try_insert(pt, Enc.seq, i+1, j-1)) {
                Enc.bp_left2=i+1;
                Enc.bp_right2=j-1;
                cnt += update_deepest(str, minim);
                // in case useFirst is on and structure is found, end
                if (Opt.first && cnt > 0) return cnt;
              } else  //outside
              if (try_insert(pt, Enc.seq, i-1, j+1)) {
                Enc.bp_left2=i-1;
                Enc.bp_right2=j+1;
                cnt += update_deepest(str, minim);
                // in case useFirst is on and structure is found, end
                if (Opt.first && cnt > 0) return cnt;
              }
            } else {
              cnt += update_deepest(str, minim);
              // in case useFirst is on and structure is found, end
              if (Opt.first && cnt > 0) return cnt;
            }
          } else {
            cnt += update_deepest(str, minim);
            // in case useFirst is on and structure is found, end
            if (Opt.first && cnt > 0) return cnt;
          }
        }
      }
    }
  }
  return cnt;
}

//shift move set
int shifts(hash_entry &str, hash_entry &minim)
{
  int cnt = 0;
  int brack_num = 0;
  short *pt = str.structure;
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
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
        }

        // switch (i,j) to (k,i)
        if (i-k>MINGAP && compat(Enc.seq[i-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=k;
          Enc.bp_right2=i;
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;

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
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
        }
        // switch (i,j) to (j,k)
        if (k-j>MINGAP && compat(Enc.seq[j-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=j;
          Enc.bp_right2=k;
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
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
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
        }

        // right switch (i,j) to (i,k)
        if (k-i>MINGAP && compat(Enc.seq[i-1], Enc.seq[k-1])) {
          Enc.bp_left=-i;
          Enc.bp_right=-j;
          Enc.bp_left2=i;
          Enc.bp_right2=k;
          cnt += update_deepest(str, minim);
          // in case useFirst is on and structure is found, end
          if (Opt.first && cnt > 0) return cnt;
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
  Deg.current = str.energy;

  if (Opt.verbose_lvl>3) fprintf(stderr, "\n  start of MS:\n  %s %d\n\n", pt_to_str(min.structure).c_str(), min.energy);

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
    if(Deg.current!=str.energy) throw;
    Deg.processed.insert(allocopy(str.structure));
    // take first to structure
    if (str.structure) free(str.structure);
    str.structure = (*Deg.unprocessed.begin());
    Deg.unprocessed.erase(Deg.unprocessed.begin());
    cnt += move_set(str);
  } else {
    // write output to str
    copy_arr(str.structure, min.structure);
    str.energy = min.energy;
  }
  // release minimal
  free(min.structure);

  // resolve degeneracy in local minima
  if (Deg.processed.size()>0) {
    Deg.processed.insert(str.structure);
    str.structure = *Deg.processed.begin();
    if(Deg.current!=str.energy) throw;
    Deg.processed.erase(Deg.processed.begin());

    if (!Deg.unprocessed.empty()) fprintf(stderr, "WARNING: Deg.unprocessed not empty!!!\n");
    Deg.Clear();
  }

  return cnt;
}

int move_rand(hash_entry &str){
  // count how many times called
  cnt_move++;

  // count better neighbours
  int cnt = 0;

  // deepest descent
  hash_entry min;
  min.structure = allocopy(str.structure);
  min.energy = str.energy;
  Deg.current = str.energy;

  if (Opt.verbose_lvl>3) fprintf(stderr, "\n  start of MR:\n  %s %d\n\n", pt_to_str(str.structure).c_str(), str.energy);

  // construct and permute possible moves
  Enc.PossMoves(str);
  Enc.Permute();

  // crawl through possible ones
 /* for (unsigned int i=0; i<Enc.moves_from.size(); i++) {
    if (str.structure[Enc.moves_from[i]]==0 && str.structure[Enc.moves_to[i]]==0) { // insert
      if (!check_insert(str, Enc.moves_from[i], Enc.moves_to[i])) continue;
      Enc.bp_left = Enc.moves_from[i];
      Enc.bp_right = Enc.moves_to[i];
    } else if (str.structure[Enc.moves_from[i]]==Enc.moves_to[i]) {  // delete
      Enc.bp_left = -Enc.moves_from[i];
      Enc.bp_right = -Enc.moves_to[i];
    } else continue;*/

  for (unsigned int i=0; i<Enc.moves_from.size(); i++) {
    Enc.bp_left = Enc.moves_from[i];
    Enc.bp_right = Enc.moves_to[i];
    //if (Opt.verbose_lvl>3) fprintf(stderr, "\n  start of MR:\n  %s %d\n\n", pt_to_str(str.structure).c_str(), str.energy);
    cnt = update_deepest(str, min);
    if (cnt) break;
  }

  // if degeneracy occurs, solve it!
  if (!cnt && Deg.unprocessed.size()>0) {
    if(Deg.current!=str.energy) throw;
    Deg.processed.insert(allocopy(str.structure));
    // take first to structure
    if (str.structure) free(str.structure);
    str.structure = (*Deg.unprocessed.begin());
    Deg.unprocessed.erase(Deg.unprocessed.begin());
    cnt += move_rand(str);
  } else {
    // write output to str
    copy_arr(str.structure, min.structure);
    str.energy = min.energy;
  }
  // release minimal
  free(min.structure);

  // resolve degeneracy in local minima
  if (Deg.processed.size()>0) {
    Deg.processed.insert(str.structure);
    str.structure = *Deg.processed.begin();
    if(Deg.current!=str.energy) throw;
    Deg.processed.erase(Deg.processed.begin());

    if (!Deg.unprocessed.empty()) fprintf(stderr, "WARNING: Deg.unprocessed not empty!!!\n");
    Deg.Clear();
  }

  return cnt;
}

// browse neighbours without degeneracy (almost the same as move_set) assume Deg.first = true
hash_entry *browse_neighs(hash_entry &str, int &saddle_en)
{
  // count how many times called
  cnt_move++;

  // returned
  hash_entry *min = (hash_entry*)space(sizeof(hash_entry));
  min->energy = str.energy;
  min->structure = allocopy(str.structure);


  // did we just find escape from basin?
  bool escape = false;

  //if (deg.opt->verbose_lvl>3) fprintf(stderr, "\n  browse neighs:\n  %s %d\n\n", pt_to_str(Enc.pt).c_str(), energy);

  escape = insertions(str, *min);
  //if (deg.opt->verbose_lvl>3) fprintf(stderr, "\n");

  // deletions
  if (!escape) escape = deletions(str, *min);

  // shifts (only if enabled + noLP disabled)
  if (Opt.shift && !Opt.noLP) {
    if (!escape) escape = shifts(str, *min);
  }

  if (escape) saddle_en = str.energy;

  if (!escape) {
    free_entry(min);
    min = NULL;
  }

  return min;
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
  //fprintf(stderr, "pt_to_str %d\n", pt[0]);
  for (int i=1; i<=pt[0]; i++) {
    if (pt[i]==0) str[i-1]='.';
    else if (pt[i]<i) str[i-1]=')';
    else str[i-1]='(';
  }
  return str;
}

// can insert be done on this structure?
bool check_insert(hash_entry &str, int i, int j)
{
  int cnt = 0;
  for (i++; i<j; i++) {
    if (str.structure[i]==0) continue;
    if (str.structure[i]>str.structure[str.structure[i]]) { //'('
      cnt++;
    } else { // ')'
      cnt--;
      if (cnt<0) return false;
    }
  }
  if (cnt!=0) return false;
  return true;
}
