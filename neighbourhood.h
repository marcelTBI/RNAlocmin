#ifndef __NEIGHBOURHOOD_H
#define __NEIGHBOURHOOD_H
#include <vector>
#include <string>

struct Neigh
{
  int i;
  int j;

  int energy_change; // = INTMAX unless evaluated

  Neigh(int i, int j, int energy = INT_MAX);
  Neigh();
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

  int GenNeighs(char *seq, short *pt);  // return next loop inside, -1 if not found
  int EvalLoop(short *pt, short *s0, short *s1, bool inside); // return energy of loop (as from loop_energy() )
};

class Neighborhood
{
private:
  static char *seq;
  static short *s0;
  static short *s1;

  std::vector<Loop*> loops;

  // for enumeration:
  int loopnum;
  int neighnum;
  std::vector<int> top_loop;
  bool deletes;

  // for degeneracy:
  static int energy_deg;
  static std::vector<Neighborhood*> degen_todo;
  static std::vector<Neighborhood*> degen_done;

public:
  short *pt;
  int energy; // = INTMAX until not evaluated;
  static int debug;

public:
  Neighborhood(char *seq, short *s0, short *s1, short *pt, bool eval = true);
  Neighborhood(const Neighborhood &second);
  ~Neighborhood();

  void Free();
  void SoftCopy(const Neighborhood &second, bool free_it = true);
  void HardCopy(const Neighborhood &second);

  bool const operator==(const Neighborhood &second) const {
    int i=0;
    while (i<pt[0] && pt[i]==second.pt[i]) i++;
    return i==pt[0];
  }

  bool const operator<(const Neighborhood &second) const;

  // move the neighbourhood:
  int AddBase(int i, int j, bool reeval = true);  // return change in the number of neighbors
  int RemBase(int i, int j, bool reeval = true);  // return change in the number of neighbors
  int ApplyNeigh(Neigh &neigh, bool reeval = true);  // return change in the number of neighbors

  int RemEnergy(short *pt, int loop, int last_loop = -1); // return the energy of a loop removal

  // printing:
  int PrintNeighs(); // return count neighbors
  int PrintEnum(bool inserts_first = true); // return count neighbors
  void PrintStr();

  // eval:
  int EvalNeighs(bool full); // evaluate the neighbourhood energies and store it efficiently, return energy of us

  // gradient descent:
  int MoveLowest(bool first = false, bool reeval = true);  // move to lowest possible bpair (gradient walk or adaptive first found walk if first = true), return CHANGE in energy
  int MoveRandom(bool reeval = true);   // move to random lowest bpair return CHANGE in energy

  // enumerating neighbors:
  void StartEnumerating(bool inserts_first = true);
  void IncreaseCount(bool inserts_first = true);
  inline void MoveLoop();
  bool NextNeighbor(Neigh &res, bool inserts_first = true, bool with_energy = true);  // return True if success, False if end of enumerating

  // degeneracy:
  bool AddDegen(Neigh &neigh);  // return True if added, False if already found.
  static void ClearDegen();
  int SolveDegen(bool random, bool reeval, int lowest = 0, bool first = false);

  static void ClearStatic();

  // debug
  std::string GetPT(Neigh &next);
};

void test();

#endif



