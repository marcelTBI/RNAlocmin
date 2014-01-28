#ifndef __MOVE_SET_PK_H
#define __MOVE_SET_PK_H

#include <set>
#include <map>
#include <string>

using namespace std;

enum  PK_TYPE         {NPK, PK_H, PK_K, PK_L, PK_M}; // types of pseudoknot - NPK is Kfree structure
const int beta1_pen[] = {0,  960, 1260, 1460, 1760};
const int beta2_pen[] = {0,   10,   41,   10,   10};
const int beta3_pen[] = {0,   10,   10,   10,   10};

class Bpair {
public:
  int start;  // left parenthesis
  int end;    // right parenthesis
  int next_left; // enclosed in multiloop starting in? (next left nested sibling)
  int next_right; // next right nested sibling
  //int nestings; // number of nestings that I encapsule (0 means im topmost nesting - left_most = this) ??? we don't need this do we?
  set<int> left_cross; // crossings from left side
  set<int> right_cross; // crossings from right side

  const bool operator<(const Bpair &second) const {
    if (start == second.start) return end < second.end;
    return start < second.start;
  }

public:
  Bpair(); // forbidden to use!!
  Bpair(int left, int right);
};

struct pk_info {
  int num_bp;    // bp total
  int num_nn_bp; // bp non-nested
  int start;
  int end;

  pk_info() {
    num_bp = 0;
    num_nn_bp = 0;
    start = 0;
    end = 0;
  }
};

class Pseudoknot {

  // data
  PK_TYPE type;
  int energy_penalty;

  map<int, Bpair> bpairs;

  int num_cross;    // total number of crossings in the PK
  //int num_ind_cross; // num of the independent crossings - non-nested ones

  // helper
  map<int, int> points;  // points to starts and ends of bps.

  // pk info:
  pk_info pki;

public:
  Pseudoknot();
  int AddBpair(int left, int right);
  int RemoveBpair(int left);

  // helpers
  int Start();
  int End();
  int Clear(); // returns the energy penalty freed

  // clears/add back all crossing pairs with base pair of choice from structure
  int ClearNeighsOfBP(short *str, int left);
  int AddNeighsOfBP(short *str, int left);

  pk_info FindPKrange(int point);
};

string pt_to_str_pk(short *str);
int try_pk();
#endif



