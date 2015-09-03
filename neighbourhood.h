#ifndef __NEIGHBOURHOOD_H
#define __NEIGHBOURHOOD_H
#include <vector>

struct Neigh
{
  int i;
  int j;

  int energy_change; // = INTMAX unless evaluated

  Neigh(int i, int j, int energy = INT_MAX);
};

struct Loop
{
  // enclosed by:
  int left;
  int right;

  int energy; // = INTMAX unless evaluated

  /*Neigh *neighs;
  int num_neighs;*/
  std::vector<Neigh> neighs;

  Loop(int i, int j);
  //Loop(Loop &second);

  // for degeneracy:
  std::vector<Neigh> lowNs;

  int GenNeighs(char *seq, short *pt);
  int EvalLoop(short *pt, short *s0, short *s1, bool inside);
};

class Neighborhood
{
  short *pt;
  static char *seq;
  static short *s0;
  static short *s1;

  std::vector<Loop*> loops;

  int energy; // = INTMAX until not evaluated;

  // for enumeration:
  int loopnum;
  int neighnum;

  // for degeneracy:
  static std::vector<Neighborhood> degen_todo;
  static std::vector<Neighborhood> degen_done;

public:
  Neighborhood(char *seq, short *pt);
  Neighborhood(Neighborhood &second);
  ~Neighborhood();

  bool const operator==(const Neighborhood &second) const {
    int i=0;
    while (i<pt[0] && pt[i]==second.pt[i]) i++;
    return i==pt[0];
  }

  // move the neighbourhood:
  int AddBase(int i, int j, bool reeval = true);
  int RemBase(int i, int j, bool reeval = true);

  int RemEnergy(short *pt, int loop, int last_loop = -1);

  // printing:
  int PrintNeighs();
  int PrintEnum();
  void PrintStr();

  // eval:
  int EvalNeighs(bool full); // evaluate the neighbourhood energies and store it efficiently

  // gradient descent:
  int MoveLowest(bool reeval = true);  // move to lowest possible bpair (gradient walk)

  // enumerating neighbors:
  void StartEnumerating();
  bool NextNeighbor(Neigh &res, bool with_energy);

  // degeneracy:
  bool AddDegen(Neigh &neigh, bool lower);
};

void test();

#endif



