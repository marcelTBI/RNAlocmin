#ifndef __MOVE_SET_PK_H
#define __MOVE_SET_PK_H

#include <set>
#include <map>

using namespace std;

enum  PK_TYPE {S, H, K, L ,M}; // types of pseudoknot - S is PKfree structure
float penalties[] = {0.0, 8.0, 12.0, 12.0, 12.0};

class Bpair {
public:
  int start;  // left parenthesis
  int next_left; // enclosed in multiloop starting in? (next left nested sibling)
  int next_right; // next right nested sibling
  //int nestings; // number of nestings that I encapsule (0 means im topmost nesting - left_most = this) ??? we don't need this do we?
  set<int> left_cross; // crossings from left side
  set<int> right_cross; // crossings from right side

  const bool operator<(const Bpair &second) const {
    return start < second.start;
  }

public:
  Bpair(); // forbidden to use!!
  Bpair(int left);
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

public:
  Pseudoknot();
  int AddBpair(short *str, int left, int right);
  int RemoveBpair(short *str, int left);

  // helpers
  int Start();
  int End();
  int Clear(); // returns the energy penalty freed
};

#endif



