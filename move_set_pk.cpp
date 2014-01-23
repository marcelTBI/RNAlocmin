#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include <iostream>

extern "C" {
  #include "pair_mat.h"
  #include "fold.h"
}

#include <vector>
#include <string>
#include <stack>
#include <queue>

#include "move_set_pk.h"

Bpair::Bpair(int left, int right)
{
  // init
  start = left;
  end = right;
  next_left = left;
  next_right = left;
}

Bpair::Bpair()
{
  start = -1;
  end = -1;
  fprintf(stderr, "wrong initialisation ERROR\n");
  throw;
}

Pseudoknot::Pseudoknot()
{
  type = NPK;
  num_cross = 0;
}

int Pseudoknot::AddBpair(int left, int right)
{
  // create new:
  Bpair new_bp(left, right);

  // some helpers
  int left_reach = 0; // left reach of the nesting - where can I furthest go with nestings
  int right_reach = INT_MAX; // right --||--

  // find crossings and nestings
  std::pair<const int, int> pr(left, 0);
  map<int, int>::iterator it = lower_bound(points.begin(), points.end(), pr); // get first point inside the basepair
  for (; it!=points.end() && it->first<=right; it++) {
    int l = it->second;
    int r = bpairs[it->second].end;

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
  std::pair<const int, int> pr2(left, 0);
  map<int, int>::reverse_iterator rit (lower_bound(points.begin(), points.end(), pr2));
  for (; rit!=points.rend() && rit->first>left_reach; rit++) {
    int l = rit->second;
    int r = bpairs[rit->second].end;

    // check if it is nested.
    if (l<left && l>left_reach && r>right && r<right_reach) {
      new_bp.next_left = l;
      break;
    }
  }

  // now insert it into the structure
  bpairs.insert(make_pair(left, new_bp));
  num_cross += new_bp.left_cross.size() + new_bp.right_cross.size();

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
    bpairs[*it].right_cross.insert(left);
  }
  for (set<int>::iterator it = new_bp.right_cross.begin(); it!=new_bp.right_cross.end(); it++) {
    bpairs[*it].left_cross.insert(left);
  }

  // do we change type? and energy penalty
    //-- still has to be done
  // non nested:
  if (new_bp.next_left != left || new_bp.next_right != left) {
    // maybe the type of PK has changed...
    // now just S->H and H -> nonH (disabled)
    bool nonH = false;
    if (new_bp.left_cross.size()>0 && new_bp.right_cross.size()>0) {
      nonH = true;
    }
    for (set<int>::iterator it = new_bp.left_cross.begin(); it!=new_bp.left_cross.end(); it++) {
      if (bpairs[*it].left_cross.size()>0) nonH = true;
    }
    for (set<int>::iterator it = new_bp.right_cross.begin(); it!=new_bp.right_cross.end(); it++) {
      if (bpairs[*it].right_cross.size()>0) nonH = true;
    }
    // non-allowed
    if (nonH) {
      RemoveBpair(left);
      return 0;
    }
    // S->H
    if (type == NPK && num_cross>0) {
      type = PK_H;
      int tmp = -energy_penalty + penalties[type];
      energy_penalty = penalties[type];
      return tmp;
    }
  }
  return 0;
}
int Pseudoknot::RemoveBpair(int left)
{
  map<int, Bpair>::iterator to_remove;
  // check existence
  if ((to_remove = bpairs.find(left))!=bpairs.end()) {
    Bpair &rem_bp = to_remove->second;
    num_cross -= rem_bp.left_cross.size(), rem_bp.right_cross.size();

    // remove the crosings:
    for (set<int>::iterator it = rem_bp.left_cross.begin(); it!=rem_bp.left_cross.end(); it++) {
      bpairs[*it].right_cross.erase(left);
    }
    for (set<int>::iterator it = rem_bp.right_cross.begin(); it!=rem_bp.right_cross.end(); it++) {
      bpairs[*it].left_cross.erase(left);
    }

    // update nested list:
    if (rem_bp.next_left != left) {
      Bpair &update = bpairs[rem_bp.next_left];
      if (rem_bp.next_right == rem_bp.start) {
        update.next_right = update.start;
      } else update.next_right = rem_bp.next_right;
    }
    if (rem_bp.next_right != left) {
      Bpair &update = bpairs[rem_bp.next_right];
      if (rem_bp.next_left == rem_bp.start) {
        update.next_left = update.start;
      } else update.next_left = rem_bp.next_left;
    }

    // update points:
    points.erase(left);
    points.erase(rem_bp.end);

    // erase the base pair
    bpairs.erase(to_remove);

    // now reset type:
     // still TODO!!!
      //now just H -> S:
    if (num_cross == 0) {
      type = NPK;
      int tmp = -energy_penalty + penalties[type];
      energy_penalty = penalties[type];
      return tmp;
    }
  }
  return 0;
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

int Pseudoknot::Clear() {
  points.clear();
  bpairs.clear();
  type = NPK;
  int tmp = -energy_penalty + penalties[type];
  energy_penalty = penalties[type];
  return tmp;
}

int Pseudoknot::ClearNeighsOfBP(short *str, int left)
{
  int res = 0;
  auto it = bpairs.find(abs(left));
  if (it==bpairs.end()) {
    return -1;
  }

  Bpair &bp = it->second;
  for (auto sit=bp.left_cross.begin(); sit!=bp.left_cross.end(); sit++) {
    str[*sit] = 0;  // start
    str[bpairs[*sit].end] = 0; // end
    res++;
  }
  for (auto sit=bp.right_cross.begin(); sit!=bp.right_cross.end(); sit++) {
    str[*sit] = 0;  // start
    str[bpairs[*sit].end] = 0; // end
    res++;
  }

  return res;
}

int Pseudoknot::AddNeighsOfBP(short *str, int left)
{
  int res = 0;
  auto it = bpairs.find(abs(left));
  if (it==bpairs.end()) return -1;

  Bpair &bp = it->second;
  for (auto sit=bp.left_cross.begin(); sit!=bp.left_cross.end(); sit++) {
    int end = bpairs[*sit].end;
    str[*sit] = end;  // start
    str[end] = *sit; // end
    res++;
  }
  for (auto sit=bp.right_cross.begin(); sit!=bp.right_cross.end(); sit++) {
    int end = bpairs[*sit].end;
    str[*sit] = end;  // start
    str[end] = *sit; // end
    res++;
  }

  return res;
}

int move_PK(Pseudoknot &PKstruct, short *str, char *seq, short *s0, short *s1, int left, int right)
{
  //init
  bool deletion = (left<0);
  int energy = 0;

  // deletion
  if (deletion) {
    // first remove it from struct
    str[-left] = 0;
    str[-right] = 0;

    // now recompute the energy
    // energy_of_move energy: (here we ask only for insertion energy: if deletion, then we perform deletion and ask for energy of reinsertion)
    // first we have to remove all conflicting edges:
    int removed = PKstruct.ClearNeighsOfBP(str, -left);
    fprintf(stderr, "%s cleared=%d\n", pt_to_str_pk(str).c_str(), removed);
    energy -= energy_of_move_pt(str, s0, s1, -left, -right);
    int added = PKstruct.AddNeighsOfBP(str, -left);
    fprintf(stderr, "%s added=%d (%8.2F)\n", pt_to_str_pk(str).c_str(), added, energy/100.0);

    //at last update pseudoknot
    energy += PKstruct.RemoveBpair(-left);
  } else { // insertion

    // first update the Pseudoknot
    energy += PKstruct.AddBpair(left, right);

    // energy_of_move energy: (here we ask only for insertion energy: if deletion, then we perform deletion and ask for energy of reinsertion)
      // first we have to remove all conflicting edges:
    int removed = PKstruct.ClearNeighsOfBP(str, left);
    fprintf(stderr, "%s cleared=%d\n", pt_to_str_pk(str).c_str(), removed);
    energy += energy_of_move_pt(str, s0, s1, left, right);
    int added = PKstruct.AddNeighsOfBP(str, left);
    fprintf(stderr, "%s added=%d (%8.2F)\n", pt_to_str_pk(str).c_str(), added, energy/100.0);


    // lastly, update structure
    str[left] = right;
    str[right] = left;
  }

  return energy;
}

/*class Stack
{
  vector<int> stck;

  int top() {return stck[stck.size()-1];}
  void push(int smth) {stck.push_back(smth);}
  bool empty() {return stck.size()==0;}
  void pop() {stck.pop_back();}
};*/

string pt_to_str_pk(short *str)
{
  // we can do with 4 types now:
  const int maxp = 4;
  const char ptl[] = {"([{<"};
  const char ptr[] = {")]}>"};

  // result:
  string res;
  res.reserve(str[0]);

  // stack of end points of each parenthesis
  vector<stack<int> > pars(maxp);

  // temprary structure
  vector<int> tmp(str[0]+1, 0);

  // precomputation
  for (int i=1; i<=str[0]; i++) {
    if (str[i]>i) {
      int j=0;
      //fprintf(stderr, "%d %d \n", (int)pars[j].size(), pars[j].empty());
      while (j<maxp && !pars[j].empty() && pars[j].top()<str[i]) {
        j++;
      }
      if (j==maxp) {
        fprintf(stderr, "Cannot print it with %d types of parentheses!!!\n", maxp);
        throw;
        return res;
      } else {
        // jth parenthesis:
        pars[j].push(str[i]);
        tmp[i] = j;
        tmp[str[i]] = j;
      }
    } else if (str[i]>0 && str[i]<i) {
      pars[tmp[i]].pop();
    }
  }

  // filling the result:
  for (int i=1; i<=str[0]; i++) {
    if (str[i]==0) res += '.';
    else if (str[i]>i) res += ptl[tmp[i]];
    else res += ptr[tmp[i]];
  }
  res+='\0';

  return res;
}

int try_pk()
{
  char *seq = "CAAUUUUCUGAAAAUUUUCAC";
  char *str = ".....................";

//CAAUUUUCUGAAAAUUUUCAC
//123456789012345678901
//:(((((::[[[))))):]]]:

  vector<std::pair<int, int> > moves;
  moves.push_back(make_pair(2,16));
  moves.push_back(make_pair(3,15));
  moves.push_back(make_pair(9,20));
  moves.push_back(make_pair(-2,-16));
  moves.push_back(make_pair(-3,-15));
  moves.push_back(make_pair(10,19));
  moves.push_back(make_pair(11,17));
  moves.push_back(make_pair(4,14));

  make_pair_matrix();
  short *pt = make_pair_table(str);
  short *s0 = encode_sequence(seq, 0);
  short *s1 = encode_sequence(seq, 1);

  int energy = energy_of_structure_pt(seq, pt, s0, s1, 0);

  Pseudoknot pk;
  fprintf(stderr, "%s %8.3f\n", pt_to_str_pk(pt).c_str(), energy/100.0);
  for (unsigned int i=0; i<moves.size(); i++) {
    fprintf(stderr, "MOVE: %d %d\n", moves[i].first, moves[i].second);
    energy += move_PK(pk, pt, seq, s0, s1, moves[i].first, moves[i].second);
    fprintf(stderr, "%s %8.3f\n", pt_to_str_pk(pt).c_str(), energy/100.0);
  }

  free(pt);
  free(s0);
  free(s1);

  return 0;
}

