#include <stdio.h>

#include <queue>

#include "barrier_tree.h"

using namespace std;

typedef struct {
  float barrier;
  int i;
  int j;
} energy_pair;

struct comparator {
    bool operator()(const energy_pair& __x, const energy_pair& __y) const
      { return __x.barrier > __y.barrier; }
};


// find father of node i (uncomputed -> return i)
int findfather(nodeT *nodes, int i) {
  if (nodes[i].father == -1) return i;
  else return findfather(nodes, nodes[i].father);
}

// make barrier tree
int make_tree(int n, float *energy_barr, nodeT *nodes)
{
  priority_queue<energy_pair, vector<energy_pair>, comparator> saddles;
  for (int i=0; i<n; i++) {
    if (nodes[i].father!=-1) continue;
    for (int j=i+1; j<n; j++) {
      if (nodes[j].father!=-1) continue;
      if (i!=j) {
        energy_pair ep;
        ep.barrier = energy_barr[i*n+j];
        ep.i=i;
        ep.j=j;
        saddles.push(ep);
      }
    }
  }

  // max_height
  float max_height = -1e10;

  // compute all except one nodes
  while (!saddles.empty()) {
    energy_pair ep = saddles.top();
    saddles.pop();

    // if not already computed
    int fatheri = findfather(nodes, ep.i);
    int fatherj = findfather(nodes, ep.j);
    if (fatheri != fatherj) {
      // merge i and j by their fathers
      if (fatheri>fatherj) swap(fatheri, fatherj);
      nodes[fatherj].saddle_height = ep.barrier;
      add_father(nodes, fatherj, fatheri);
      //nodes[fatherj].father = fatheri;
      if (ep.barrier>max_height) max_height = ep.barrier;
    }
  }

  // finish the last one
  nodes[0].saddle_height = max_height + 0.1;

  return 0;
}


void add_father(nodeT *nodes, int child, int father)
{
  //search for old father
  int old_father = nodes[child].father;
  set<int>::iterator it;
  if (old_father !=-1 && (it = nodes[old_father].children.find(child))!=nodes[old_father].children.end()) {
    nodes[old_father].children.erase(it);
  }

  // add new one
  nodes[child].father = father;
  nodes[father].children.insert(child);

  // recompute others
  set<int> tmp = nodes[child].children;
  for (it=tmp.begin(); it!=tmp.end(); it++) {
    if (nodes[*it].saddle_height > nodes[child].saddle_height) {

      //nodes[*it].father = father;
      add_father(nodes, *it, father);
    }
  }
}
