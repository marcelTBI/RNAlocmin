#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "move_set.h"

#include "pair_mat.h"
#include "fold.h"

/* maximum degeneracy value - if degeneracy is greater than this, program segfaults*/
#define MAX_DEGEN 100
#define MINGAP 3


#define space(a) malloc(a)


class Bpair {
  int start;  // left parenthesis
  int next_left; // enclosed in multiloop starting in? (next left nested sibling)
  int next_right; // next right nested sibling
  //int nestings; // number of nestings that I encapsule (0 means im topmost nesting - left_most = this) ??? we don't need this do we?
  set<int> left_cross; // crossings from left side
  set<int> right_cross; // crossings from right side

  const bool operator<(const Bpair &second) const {
    return start < second.start;
  }

Bpair::Bpair(int left)
{
  // init
  start = left;
  next_left = left;
  next_right = left;
}

class Pseudoknot {

  // data
  PK_TYPE type;
  int energy_penalty;

  map<int, Bpair> bpairs;  // map starts of bpairs to the data

  // helper
  map<int, int> points;  // points to starts and ends of bps. -> points to starts always points[8]=5 means bp 5-8 is there and there should be also points[5]=5

  Pseudoknot::Pseudoknot()
  {
    type = S;
  }

  int Pseudoknot::AddBpair(short *str, int left, int right)
  {
    // create new:
    Bpair new_bp(left);

    // some helpers
    int left_reach = 0; // left reach of the nesting - where can I furthest go with nestings
    int right_reach = str[0]+1; // right --||--

    // find crossings and nestings
    map<int, int>::iterator it = lower_bound(points.begin(), points.end(), left); // get first point inside the basepair
    for (; it!=points.end() && it->first<=right; it++) {
      int l = it->second;
      int r = str[it->second];

      // 3 options:
      if (l<left && r<right) { //left cross
        new_bp.left_cross.insert(l);
        left_reach = max(left_reach, l);
      } else if (l>left && r>right) { // right cross
        new_bp.right_cross.insert(l);
        right_reach = min(right_reach, r);
      } else if (l>left && r<right) { // nested, need to resolve left sibling later
        // we will not keep info other than next right next left
        if (new_bp.next_right == left) new_bp.next_right = l;
      }
    }

    // left_most missing only: (we are searching in (left_reach, left) region backwards)
    map<int, int>::reverse_iterator it (lower_bound(points.begin(), points.end(), left);
    for (; it!=points.rend() && it->first>left_reach; it++) {
      int l = it->second;
      int r = str[it->second];

      // check if it is nested.
      if (l<left && l>left_reach && r>right && r<right_reach) {
        new_bp.next_left = l;
        break;
      }
    }

    // now insert it into the structure
    bpairs.insert(make_pair(left, new_bp));

    // update points
    points[left] = left;
    points[right] = left;

    // update nested list:
    if (new_bp.next_left != left) {
      bpairs[new_bp.next_left].next_right = left;
    }
    if (new_bp.next_right != left) {
      bpairs[new_bp.next_right].next_left = left;
    }

    // update crossings lists:
    for (set<int>::iterator it = new_bp.left_cross.begin(); it!=new_bp.left_cross.end(); it++) {
      bpairs[it->start].right_cross.insert(start);
    }
    for (set<int>::iterator it = new_bp.right_cross.begin(); it!=new_bp.right_cross.end(); it++) {
      bpairs[it->start].left_cross.insert(start);
    }

    // do we change type? and energy penalty
      //-- still has to be done
    // non nested:
    if (new_bp.next_left != left || new_bp_next_right != left) {
      // maybe the type of PK has changed...
    }

  }
  int Pseudoknot::RemoveBpair(int left)
  {
    map<int, Bpair>:iterator to_remove;
    // check existence
    if ((to_remove = bpairs.find(left))!=bpairs.end()) {
      bpairs.erase(to_remove);
      // TODO
    }
  }

  // helpers
  int Pseudoknot::Start() {
    if (points.size()==0) return 0;
    return points.begin()->first;
  }
  int Pseudoknot::End() {
    if (points.size()==0) return 0;
    return points.rbegin()->first;
  }
}


int compare(short *lhs, short *rhs)
{

  //printf("%d ", (int)lhs[0]);

  int i=1;
  char l=0,r=0;
  while (i<=lhs[0]) {
    l = (lhs[i]==0?'.':(lhs[i]<lhs[lhs[i]]?')':'('));
    r = (rhs[i]==0?'.':(rhs[i]<rhs[rhs[i]]?')':'('));
    if (l != r) break;
    i++;
  }

  return (i<=lhs[0] && l<r);
}

int find_min(short *arr[MAX_DEGEN], int begin, int end) {
  short *min = arr[begin];
  short min_num = begin;
  int i;

  for (i=begin+1; i<end; i++) {
    if (compare(arr[i], min)) {
      min = arr[i];
      min_num = i;
    }
  }
  return min_num;
}

int equals(const short *first, const short *second)
{
  int i=1;
  while (i<=first[0] && first[i]==second[i]) {
    i++;
  }
  if (i>first[0]) return 1;
  else return 0;
}

void copy_arr(short *dest, const short *src)
{
  if (!src || !dest) {
    fprintf(stderr, "Empty pointer in copying\n");
    return;
  }
  memcpy(dest, src, sizeof(short)*(src[0]+1));
}

short *allocopy(const short *src)
{
  short *res = (short*) space(sizeof(short)*(src[0]+1));
  copy_arr(res, src);
  return res;
}

void copy_se(struct_en *dest, const struct_en *src) {
	if (!dest || !src) return ;
	dest->energy = src->energy;
	copy_arr(dest->structure, src->structure);
}

struct_en *allocopy_se(const struct_en *input) {
	struct_en *output = (struct_en*) space(sizeof(struct_en));
	output->energy = input->energy;
	output->structure = allocopy(input->structure);
	return output;
}

void free_se(struct_en *to_free) {
	if (!to_free) return;
	free(to_free->structure);
	free(to_free);
	to_free = NULL;
}

/* ############################## DECLARATION #####################################*/
/* private functions & declarations*/

static int cnt_move = 0;
int count_move() {return cnt_move;}

/*check if base is lone*/
int lone_base(short *pt, int i);

/* can move be done on this structure? (move is in the Enc)*/
int check_insert(struct_en *str, int i, int j);

/* internal struct with moves, sequence, degeneracy and options*/
typedef struct _Encoded {
  /* sequence*/
  short *s0;
  short *s1;

  const char  *seq;

  /* moves*/
  int   bp_left;
  int   bp_right;
  int   bp_left2;   /* if noLP is enabled (and for shift moves)*/
  int   bp_right2;

  /* options*/
  int noLP;
  int verbose_lvl;
  int first;
  int shift;

  /* degeneracy*/
  int begin_unpr;
  int begin_pr;
  int end_unpr;
  int end_pr;
  short *processed[MAX_DEGEN];
  short *unprocessed[MAX_DEGEN];
  int current_en;

  /* moves in random (needs to be freed afterwards)*/
  int *moves_from;
  int *moves_to;
  int num_moves;

  /* function for flooding */
  int (*funct) (struct_en*, struct_en*);


} Encoded;

/* frees all things allocated by degeneracy...*/
void free_degen(Encoded *Enc)
{
  int i;
  for (i=Enc->begin_unpr; i<Enc->end_unpr; i++) {
    if (Enc->unprocessed[i]) {
      free(Enc->unprocessed[i]);
      Enc->unprocessed[i]=NULL;
    }
  }
  for (i=Enc->begin_pr; i<Enc->end_pr; i++) {
    if (Enc->processed[i]) {
      free(Enc->processed[i]);
      Enc->processed[i]=NULL;
    }
  }
  Enc->begin_pr=0;
  Enc->begin_unpr=0;
  Enc->end_pr=0;
  Enc->end_unpr=0;
}

/* ############################## IMPLEMENTATION #####################################*/

/* reads a line no matter how long*/
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

inline void do_move(short *pt, int bp_left, int bp_right)
{
  /* delete*/
  if (bp_left<0) {
    pt[-bp_left]=0;
    pt[-bp_right]=0;
  } else { /* insert*/
    pt[bp_left]=bp_right;
    pt[bp_right]=bp_left;
  }
}

/* done with all structures along the way to deepest*/
int update_deepest(Encoded *Enc, struct_en *str, struct_en *min)
{
  /* apply move + get its energy*/
  int tmp_en;
  tmp_en = str->energy + energy_of_move_pt(str->structure, Enc->s0, Enc->s1, Enc->bp_left, Enc->bp_right);
  do_move(str->structure, Enc->bp_left, Enc->bp_right);
  if (Enc->bp_left2 != 0) {
    tmp_en += energy_of_move_pt(str->structure, Enc->s0, Enc->s1, Enc->bp_left2, Enc->bp_right2);
    do_move(str->structure, Enc->bp_left2, Enc->bp_right2);
  }
  int last_en = str->energy;
  str->energy = tmp_en;


  /* use f_point if we have it */
  if (Enc->funct) {
    int end = Enc->funct(str, min);

    // undo moves
    if (Enc->bp_left2!=0) do_move(str->structure, -Enc->bp_left2, -Enc->bp_right2);
    do_move(str->structure, -Enc->bp_left, -Enc->bp_right);
    str->energy = last_en;
    Enc->bp_left=0;
    Enc->bp_right=0;
    Enc->bp_left2=0;
    Enc->bp_right2=0;

    return (end?1:0);
  }

  if (Enc->verbose_lvl>1) { fprintf(stderr, "  "); print_str(stderr, str->structure); fprintf(stderr, " %d\n", tmp_en); }

  /* better deepest*/
  if (tmp_en < min->energy) {
    min->energy = tmp_en;
    copy_arr(min->structure, str->structure);

    /* delete degeneracy*/
    free_degen(Enc);

    /* undo moves*/
    if (Enc->bp_left2!=0) do_move(str->structure, -Enc->bp_left2, -Enc->bp_right2);
    do_move(str->structure, -Enc->bp_left, -Enc->bp_right);
    str->energy = last_en;
    Enc->bp_left=0;
    Enc->bp_right=0;
    Enc->bp_left2=0;
    Enc->bp_right2=0;
    return 1;
  }

  /* degeneracy*/
  if ((str->energy == min->energy) && (Enc->current_en == min->energy)) {
    int found = 0;
    int i;
    for (i=Enc->begin_pr; i<Enc->end_pr; i++) {
      if (equals(Enc->processed[i], str->structure)) {
        found = 1;
        break;
      }
    }
    for (i=Enc->begin_unpr; !found && i<Enc->end_unpr; i++) {
      if (equals(Enc->unprocessed[i], str->structure)) {
        found = 1;
        break;
      }
    }

    if (!found) {
      //print_stren(stderr, str); // fprintf(stderr, " %6.2f\n", str->energy);
      Enc->unprocessed[Enc->end_unpr]=allocopy(str->structure);
      Enc->end_unpr++;
    }
  }

  /* undo moves*/
  if (Enc->bp_left2!=0) do_move(str->structure, -Enc->bp_left2, -Enc->bp_right2);
  do_move(str->structure, -Enc->bp_left, -Enc->bp_right);
  str->energy = last_en;
  Enc->bp_left=0;
  Enc->bp_right=0;
  Enc->bp_left2=0;
  Enc->bp_right2=0;
  return 0;
}


/* deletions move set*/
int deletions(Encoded *Enc, struct_en *str, struct_en *minim)
{
  int cnt = 0;
  short *pt = str->structure;
  int len = pt[0];
  int i;

  for (i=1; i<=len; i++) {
    if (pt[i]>pt[pt[i]]) {  /* '('*/
      Enc->bp_left=-i;
      Enc->bp_right=-pt[i];

      /*if nolp enabled, make (maybe) 2nd delete*/
      if (Enc->noLP) {
        int lone = -1;
        if (lone_base(pt, i-1)) lone=i-1;
        else if (lone_base(pt, i+1)) lone=i+1;
        else if (lone_base(pt, pt[i]-1)) lone=pt[i]-1;
        else if (lone_base(pt, pt[i]+1)) lone=pt[i]+1;

        /* check*/
        if (lone != -1 && (pt[lone]==0 || pt[pt[lone]]==0)) {
          fprintf(stderr, "WARNING: pt[%d(or %d)]!=\'.\'", lone, pt[lone]);
        }

        if (lone != -1) {
          Enc->bp_left2=-lone-1;
          Enc->bp_right2=-pt[lone]-1;
        }
        if (!lone_base(pt, pt[lone]-1) && !lone_base(pt, pt[lone]+1)) {
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }
      } else {  /* nolp not enabled*/
        cnt += update_deepest(Enc, str, minim);
        /* in case useFirst is on and structure is found, end*/
        if (Enc->first && cnt > 0) return cnt;
      }
    }
  }
  return cnt;
}

  /* compatible base pair?*/
inline bool compat(char a, char b) {
  if (a=='A' && b=='U') return true;
  if (a=='C' && b=='G') return true;
  if (a=='G' && b=='U') return true;
  if (a=='U' && b=='A') return true;
  if (a=='G' && b=='C') return true;
  if (a=='U' && b=='G') return true;
  /* and with T's*/
  if (a=='A' && b=='T') return true;
  if (a=='T' && b=='A') return true;
  if (a=='G' && b=='T') return true;
  if (a=='T' && b=='G') return true;
  return false;
}

/* try insert base pair (i,j)*/
inline bool try_insert(const short *pt, const char *seq, int i, int j)
{
  if (i<=0 || j<=0 || i>pt[0] || j>pt[0]) return false;
  return (j-i>MINGAP && pt[j]==0 && pt[i]==0 && compat(seq[i-1], seq[j-1]));
}

// try insert base pair (i,j)
inline bool try_insert_seq(const char *seq, int i, int j)
{
  if (i<=0 || j<=0) return false;
  return (j-i>MINGAP && compat(seq[i-1], seq[j-1]));
}

/* insertions move set*/
int insertions(Encoded *Enc, struct_en *str, struct_en *minim)
{
  int cnt = 0;
  short *pt = str->structure;
  int len = pt[0];
  int i,j;

  for (i=1; i<=len; i++) {
    if (pt[i]==0) {
      for (j=i+1; j<=len; j++) {
        /* end if found closing bracket*/
        if (pt[j]!=0 && pt[j]<j) break;  /*')'*/
        if (pt[j]!=0 && pt[j]>j) {       /*'('*/
          j = pt[j];
          continue;
        }
        /* if conditions are met, do insert*/
        if (try_insert(pt, Enc->seq, i, j)) {
          Enc->bp_left=i;
          Enc->bp_right=j;

          if (Enc->noLP) {
            /* if lone bases occur, try inserting one another base*/
            if (lone_base(pt, i) || lone_base(pt, j)) {
              /* inside*/
              if (try_insert(pt, Enc->seq, i+1, j-1)) {
                Enc->bp_left2=i+1;
                Enc->bp_right2=j-1;
                cnt += update_deepest(Enc, str, minim);
                /* in case useFirst is on and structure is found, end*/
                if (Enc->first && cnt > 0) return cnt;
              } else  /*outside*/
              if (try_insert(pt, Enc->seq, i-1, j+1)) {
                Enc->bp_left2=i-1;
                Enc->bp_right2=j+1;
                cnt += update_deepest(Enc, str, minim);
                /* in case useFirst is on and structure is found, end*/
                if (Enc->first && cnt > 0) return cnt;
              }
            } else {
              cnt += update_deepest(Enc, str, minim);
              /* in case useFirst is on and structure is found, end*/
              if (Enc->first && cnt > 0) return cnt;
            }
          } else {
            cnt += update_deepest(Enc, str, minim);
            /* in case useFirst is on and structure is found, end*/
            if (Enc->first && cnt > 0) return cnt;
          }
        }
      }
    }
  }
  return cnt;
}

/*shift move set*/
int shifts(Encoded *Enc, struct_en *str, struct_en *minim)
{
  int cnt = 0;
  int brack_num = 0;
  short *pt = str->structure;
  int len = pt[0];
  int i, k;

  for (i=1; i<=len; i++) {
    if (pt[i]!=0 && pt[i]>i) {  /*'('*/
      int j=pt[i];

      /* outer switch left*/
      if (Enc->verbose_lvl>1) fprintf(stderr, "%2d bracket %2d position, outer switch left\n", brack_num+1, i);
      for (k=i-1; k>0; k--) {
        if (pt[k]!=0 && pt[k]>k/*'('*/) break;
        if (pt[k]!=0 && pt[k]<k/*')'*/) {
          k = pt[k];
          continue;
        }
        /* checks*/
        if (pt[k]!=0) {
          fprintf(stderr, "WARNING: \'%c\'should be \'.\' at pos %d!\n", pt[k], k);
        }

        /* switch (i,j) to (k,j)*/
        if (j-k>MINGAP && compat(Enc->seq[k-1], Enc->seq[j-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=k;
          Enc->bp_right2=j;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }

        /* switch (i,j) to (k,i)*/
        if (i-k>MINGAP && compat(Enc->seq[i-1], Enc->seq[k-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=k;
          Enc->bp_right2=i;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;

        }
      }

      /* outer switch right*/
      if (Enc->verbose_lvl>1) fprintf(stderr, "%2d bracket %2d position, outer switch right\n", brack_num+1, i);
      for (k=j+1; k<=len; k++) {
        if (pt[k]!=0 && pt[k]<k/*')'*/) break;
        if (pt[k]!=0 && pt[k]>k/*'('*/) {
          k = pt[k];
          continue;
        }

        /* check*/
        if (pt[k]!=0) {
          fprintf(stderr, "WARNING: \'%c\'should be \'.\' at pos %d!\n", pt[k], k);
        }
        /* switch (i,j) to (i,k)*/
        if (k-i>MINGAP && compat(Enc->seq[i-1], Enc->seq[k-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=i;
          Enc->bp_right2=k;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }
        /* switch (i,j) to (j,k)*/
        if (k-j>MINGAP && compat(Enc->seq[j-1], Enc->seq[k-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=j;
          Enc->bp_right2=k;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }
      }

      if (Enc->verbose_lvl>1) fprintf(stderr, "%2d bracket %2d position, inner switch\n", brack_num+1, i);
      /* inner switch*/
      for (k=i+1; k<j; k++) {
        /* jump to end of the sub-bracketing*/
        if (pt[k]!=0 && pt[k]>k/*'('*/) {
            k=pt[k];
            continue;
        }

        /* left switch (i,j) to (k,j)*/
        if (j-k>MINGAP && compat(Enc->seq[k-1], Enc->seq[j-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=k;
          Enc->bp_right2=j;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }

        /* right switch (i,j) to (i,k)*/
        if (k-i>MINGAP && compat(Enc->seq[i-1], Enc->seq[k-1])) {
          Enc->bp_left=-i;
          Enc->bp_right=-j;
          Enc->bp_left2=i;
          Enc->bp_right2=k;
          cnt += update_deepest(Enc, str, minim);
          /* in case useFirst is on and structure is found, end*/
          if (Enc->first && cnt > 0) return cnt;
        }
      } /* end inner switch for*/
      brack_num++;
    } /* end if (pt[i]=='(')*/
  } /* end for in switches*/
  return cnt;
}

/* move to deepest (or first) neighbour*/
int move_set(Encoded *Enc, struct_en *str)
{
  /* count how many times called*/
  cnt_move++;

  /* count better neighbours*/
  int cnt = 0;

  /* deepest descent*/
  struct_en min;
  min.structure = allocopy(str->structure);
  min.energy = str->energy;
  Enc->current_en = str->energy;

  if (Enc->verbose_lvl>0) { fprintf(stderr, "  start of MS:\n  "); print_str(stderr, str->structure); fprintf(stderr, " %d\n\n", str->energy); }

  /* if using first dont do all of them*/
  bool end = false;
  /* insertions*/
  if (!end) cnt += insertions(Enc, str, &min);
  if (Enc->first && cnt>0) end = true;
  if (Enc->verbose_lvl>1) fprintf(stderr, "\n");

  /* deletions*/
  if (!end) cnt += deletions(Enc, str, &min);
  if (Enc->first && cnt>0) end = true;

  /* shifts (only if enabled + noLP disabled)*/
  if (!end && Enc->shift && !Enc->noLP) {
    cnt += shifts(Enc, str, &min);
    if (Enc->first && cnt>0) end = true;
  }

  /* if degeneracy occurs, solve it!*/
  if (!end && (Enc->end_unpr - Enc->begin_unpr)>0) {
    Enc->processed[Enc->end_pr] = str->structure;
    Enc->end_pr++;
    str->structure = Enc->unprocessed[Enc->begin_unpr];
    Enc->unprocessed[Enc->begin_unpr]=NULL;
    Enc->begin_unpr++;
    cnt += move_set(Enc, str);
  } else {
    /* write output to str*/
    copy_arr(str->structure, min.structure);
    str->energy = min.energy;
  }
  /* release minimal*/
  free(min.structure);

  /* resolve degeneracy in local minima*/
  if ((Enc->end_pr - Enc->begin_pr)>0) {
    Enc->processed[Enc->end_pr]=str->structure;
    Enc->end_pr++;

    int min = find_min(Enc->processed, Enc->begin_pr, Enc->end_pr);
    short *tmp = Enc->processed[min];
    Enc->processed[min] = Enc->processed[Enc->begin_pr];
    Enc->processed[Enc->begin_pr] = tmp;
    str->structure = Enc->processed[Enc->begin_pr];
    Enc->begin_pr++;
    free_degen(Enc);
  }

  if (Enc->verbose_lvl>1 && !(Enc->first)) { fprintf(stderr, "\n  end of MS:\n  "); print_str(stderr, str->structure); fprintf(stderr, " %d\n\n", str->energy); }

  return cnt;
}

void construct_moves(Encoded *Enc, short *structure)
{
  /* generate all possible moves (less than n^2)*/
  Enc->num_moves = 0;
  int i;
  for (i=1; i<=structure[0]; i++) {
    if (structure[i]!=0) {
      if (structure[i]<i) continue;
      Enc->moves_from[Enc->num_moves]=-i;
      Enc->moves_to[Enc->num_moves]=-structure[i];
      Enc->num_moves++;
      //fprintf(stderr, "add  d(%d, %d)\n", i, str.structure[i]);
    } else {
      int j;
      for (j=i+1; j<=structure[0]; j++) {
        //fprintf(stderr, "check (%d, %d)\n", i, j);
        if (structure[j]==0) {
          if (try_insert_seq(Enc->seq,i,j)) {
            Enc->moves_from[Enc->num_moves]=i;
            Enc->moves_to[Enc->num_moves]=j;
            Enc->num_moves++;
            //fprintf(stderr, "add  i(%d, %d)\n", i, j);
            continue;
          }
        } else if (structure[j]>j) { // '('
          j = structure[j];
        } else break;
      }
    }
  }

  /* permute them */
  for (i=0; i<Enc->num_moves-1; i++) {
    int rnd = rand();
    rnd = rnd % (Enc->num_moves-i) + i;
    int swp;
    swp = Enc->moves_from[i];
    Enc->moves_from[i]=Enc->moves_from[rnd];
    Enc->moves_from[rnd]=swp;
    swp = Enc->moves_to[i];
    Enc->moves_to[i]=Enc->moves_to[rnd];
    Enc->moves_to[rnd]=swp;
  }
}

int move_rset(Encoded *Enc, struct_en *str)
{
  /* count how many times called*/
  cnt_move++;

  /* count better neighbours*/
  int cnt = 0;

  /* deepest descent*/
  struct_en min;
  min.structure = allocopy(str->structure);
  min.energy = str->energy;
  Enc->current_en = str->energy;

  if (Enc->verbose_lvl>0) { fprintf(stderr, "  start of MR:\n  "); print_str(stderr, str->structure); fprintf(stderr, " %d\n\n", str->energy); }

  // construct and permute possible moves
  construct_moves(Enc, str->structure);

  /* find first lower one*/
  int i;
  for (i=0; i<Enc->num_moves; i++) {
    Enc->bp_left = Enc->moves_from[i];
    Enc->bp_right = Enc->moves_to[i];
    cnt = update_deepest(Enc, str, &min);
    if (cnt) break;
  }

  /* if degeneracy occurs, solve it!*/
  if (!cnt && (Enc->end_unpr - Enc->begin_unpr)>0) {
    Enc->processed[Enc->end_pr] = str->structure;
    Enc->end_pr++;
    str->structure = Enc->unprocessed[Enc->begin_unpr];
    Enc->unprocessed[Enc->begin_unpr]=NULL;
    Enc->begin_unpr++;
    cnt += move_rset(Enc, str);
  } else {
    /* write output to str*/
    copy_arr(str->structure, min.structure);
    str->energy = min.energy;
  }
  /* release minimal*/
  free(min.structure);

  /* resolve degeneracy in local minima*/
  if ((Enc->end_pr - Enc->begin_pr)>0) {
    Enc->processed[Enc->end_pr]=str->structure;
    Enc->end_pr++;

    int min = find_min(Enc->processed, Enc->begin_pr, Enc->end_pr);
    short *tmp = Enc->processed[min];
    Enc->processed[min] = Enc->processed[Enc->begin_pr];
    Enc->processed[Enc->begin_pr] = tmp;
    str->structure = Enc->processed[Enc->begin_pr];
    Enc->begin_pr++;
    free_degen(Enc);
  }

  return cnt;
}

/*check if base is lone*/
int lone_base(short *pt, int i)
{
  if (i<=0 || i>pt[0]) return 0;
  /* is not a base pair*/
  if (pt[i]==0) return 0;

  /* base is lone:*/
  if (i-1>0) {
    /* is base pair and is the same bracket*/
    if (pt[i-1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return 0;
  }

  if (i+1<=pt[0]) {
    if (pt[i+1]!=0 && ((pt[i-1]<pt[pt[i-1]]) == (pt[i]<pt[pt[i]]))) return 0;
  }

  return 1;
}

/* if the structure has lone pairs*/
int find_lone_pair(short* str)
{
  int i;
  for(i=1; i<str[0]; i++) {
    if (str[i]==0) continue; /* '.'*/

    if (str[i]>str[str[i]]) {  /* '('*/
      if (i+1==str[0] || str[i+1]==0 || str[i+1]<str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]>str[str[i+1]]) i++;
    }

    if (str[i]<str[str[i]]) {  /* ')'*/
      if (i+1==str[0] || str[i+1]==0 || str[i+1]>str[str[i+1]]) {
        return i;
      } else while (i+1!=str[0] && str[i+1]!=0 && str[i+1]<str[str[i+1]]) i++;
    }
  }

  return -1;
}

int move_standard(char *seq,
                  char *struc,
                  enum MOVE_TYPE type,
                  int verbosity_level,
                  int shifts,
                  int noLP)
{
  make_pair_matrix();

  short int *s0 = encode_sequence(seq, 0);
  short int *s1 = encode_sequence(seq, 1);
  short int *str = make_pair_table(struc);

  int energy = 0;
  switch (type){
  case GRADIENT: energy = move_gradient(seq, str, s0, s1, verbosity_level, shifts, noLP); break;
  case FIRST: energy = move_first(seq, str, s0, s1, verbosity_level, shifts, noLP); break;
  case ADAPTIVE: energy = move_adaptive(seq, str, s0, s1, verbosity_level); break;
  }

  int i=1;
  for (; i<=str[0]; i++) {
    if (str[i]==0) struc[i-1]='.';
    else if (str[i]>str[str[i]]) struc[i-1]='(';
      else struc[i-1]=')';
  }

  free(s0);
  free(s1);
  free(str);

  return energy;
}

int move_gradient( char *string,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level,
                  int shifts,
                  int noLP)
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = string;
  enc.s0 = s;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;
  enc.bp_left2=0;
  enc.bp_right2=0;

  /* options*/
  enc.noLP=noLP;
  enc.verbose_lvl=verbosity_level;
  enc.first=0;
  enc.shift=shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  struct_en str;
  str.structure = allocopy(ptable);
  str.energy = energy_of_structure_pt(enc.seq, str.structure, enc.s0, enc.s1, 0);

  while (move_set(&enc, &str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  copy_arr(ptable, str.structure);
  free(str.structure);

  return str.energy;
}

int move_first(   char *string,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level,
                  int shifts,
                  int noLP)
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = string;
  enc.s0 = s;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;
  enc.bp_left2=0;
  enc.bp_right2=0;

  /* options*/
  enc.noLP=noLP;
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.shift=shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  struct_en str;
  str.structure = allocopy(ptable);
  str.energy = energy_of_structure_pt(enc.seq, str.structure, enc.s0, enc.s1, 0);

  while (move_set(&enc, &str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  copy_arr(ptable, str.structure);
  free(str.structure);

  return str.energy;
}

int move_adaptive(char *string,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level)
{
  srand(time(NULL));

  cnt_move = 0;

  Encoded enc;
  enc.seq = string;
  enc.s0 = s;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;
  enc.bp_left2=0;
  enc.bp_right2=0;

  /* options*/
  enc.noLP=0;
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.shift=0;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=NULL;

  // allocate memory for moves
  enc.moves_from = (int*) space(ptable[0]*ptable[0]*sizeof(int));
  enc.moves_to = (int*) space(ptable[0]*ptable[0]*sizeof(int));

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  struct_en str;
  str.structure = allocopy(ptable);
  str.energy = energy_of_structure_pt(enc.seq, str.structure, enc.s0, enc.s1, 0);

  while (move_rset(&enc, &str)!=0) {
    free_degen(&enc);
  }
  free_degen(&enc);

  copy_arr(ptable, str.structure);
  free(str.structure);
  free(enc.moves_from);
  free(enc.moves_to);

  return str.energy;
}

int browse_neighs( char *seq,
                   char *struc,
                   int verbosity_level,
                   int shifts,
                   int noLP,
                   int (*funct) (struct_en*, struct_en*))

{
  make_pair_matrix();

  short int *s0 = encode_sequence(seq, 0);
  short int *s1 = encode_sequence(seq, 1);
  short int *str = make_pair_table(struc);

  int res = browse_neighs_pt(seq, str, s0, s1, verbosity_level, shifts, noLP, funct);

  int i=1;
  for (; i<=str[0]; i++) {
    if (str[i] == 0) {
      struc[i-1] = '.';
    } else if (str[i]>str[str[i]]) {
      struc[i-1] = '(';
    } else {
      struc[i-1] = ')';
    }
  }

  free(s0);
  free(s1);
  free(str);

  return res;
}

int browse_neighs_pt(   char *string,
                  short *ptable,
                  short *s,
                  short *s1,
                  int verbosity_level,
                  int shifts,
                  int noLP,
                  int (*funct) (struct_en*, struct_en*))
{
  cnt_move = 0;

  Encoded enc;
  enc.seq = string;
  enc.s0 = s;
  enc.s1 = s1;

  /* moves*/
  enc.bp_left=0;
  enc.bp_right=0;
  enc.bp_left2=0;
  enc.bp_right2=0;

  /* options*/
  enc.noLP=noLP;
  enc.verbose_lvl=verbosity_level;
  enc.first=1;
  enc.shift=shifts;

  /* degeneracy*/
  enc.begin_unpr=0;
  enc.begin_pr=0;
  enc.end_unpr=0;
  enc.end_pr=0;
  enc.current_en=0;

  // function
  enc.funct=funct;

  int i;
  for (i=0; i<MAX_DEGEN; i++) enc.processed[i]=enc.unprocessed[i]=NULL;

  struct_en str;
  str.structure = allocopy(ptable);
  str.energy = energy_of_structure_pt(enc.seq, str.structure, enc.s0, enc.s1, 0);

  move_set(&enc, &str);
  free_degen(&enc);

  copy_arr(ptable, str.structure);
  free(str.structure);

  return str.energy;
}

/* printf*/
void print_stren(FILE *out, struct_en *str) {
  print_str(out, str->structure);
  fprintf(out, " %6.2f\n", str->energy/100.0);
}

void print_str(FILE *out, short *str) {
  int i;
  for (i=1; i<=str[0]; i++) {
    if (str[i]==0) fprintf(out, ".");
    else if (str[i]<i) fprintf(out, ")");
    else fprintf(out, "(");
  }
}

/*
// sample usage:
int main() {
  char seq[20] = "ACCCCCCTCTGTAGGGGGA";
  char str[20] = ".((.(.........).)).";

  // move to the local minimum and display it
  int energy = move_standard(seq, str, GRADIENT, 0, 0, 0);
  fprintf(stdout, "%s %6.2f\n\n", str, energy/100.0);

  //now create an array of every structure in neighbourhood of str structure
  struct_en *list = NULL;
  int list_length = 0;

  int get_list(struct_en *new_one, struct_en *old_one)
  {
    // enlarge the list
    list_length++;
    list = (struct_en*) realloc(list, list_length*sizeof(struct_en));

    // copy the structure
    list[list_length-1].energy = new_one->energy;
    list[list_length-1].structure = allocopy(new_one->structure);

    // we want to continue -> return 0
    return 0;
  }
  browse_neighs(seq, str, 0, 0, 0, get_list);

  // print them and free the memory:
  int i;
  for (i=0; i<list_length; i++) {
    print_stren(stdout, &list[i]);
    free(list[i].structure);
  }
  free(list);

  return 0;
}*/


