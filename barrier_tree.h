#include "treeplot.h"

// make barrier tree
int make_tree(int n, float *energy_bar, nodeT *nodes);

// recompute single father change
void add_father(nodeT *nodes, int child, int father, double color);
