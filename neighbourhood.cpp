#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "neighbourhood.h"
#include "RNAlocmin.h"

extern "C" {
  #include "move_set_inside.h"
  #include "fold.h"
  #include "pair_mat.h"
}

#define MINGAP 3

void error_message(char *str, int i = -1, int j = -1, int k = -1, int l = -1)
{
  fprintf(stderr, str, i, j, k, l);
  exit(EXIT_FAILURE);
}

Neigh::Neigh(int i, int j, int energy)
{
  this->i = i;
  this->j = j;

  energy_change = energy;
}

Neigh::Neigh():Neigh(-1,-1)
{
}

Loop::Loop(int i, int j)
{
  this->left = i;
  this->right = j;

  energy = INT_MAX;
}

/*Loop::Loop(Loop &second)
{
  *this = second;
}*/

inline bool
compat(char a, char b){

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

int Loop::GenNeighs(char *seq, short *pt)
{
  int res = -1;
  neighs.clear();

  for (int i=left+1; i<right; i++) {
    if (pt[i]>i) {   // '('
      res = i;
      i = pt[i];
      continue;
    }
    for (int j=i+1; j<right; j++) {
      if (pt[j]>j) {
        j = pt[j];
        continue;
      }
      if (j-i>MINGAP && pt[j] == 0 && compat(seq[i-1], seq[j-1])) {
        neighs.push_back(Neigh(i,j));
      }
    }
  }

  return res;
}

int Loop::EvalLoop(short *pt, short *s0, short *s1, bool inside)
{
  energy = loop_energy(pt, s0, s1, left);

  if (inside) {
    for (int i=0; i<(int)neighs.size(); i++) {
      pt[neighs[i].i] = neighs[i].j;
      pt[neighs[i].j] = neighs[i].i;
      neighs[i].energy_change = -energy + loop_energy(pt, s0, s1, neighs[i].i) + loop_energy(pt, s0, s1, left);
      pt[neighs[i].i] = 0;
      pt[neighs[i].j] = 0;
    }
  }

  return energy;
}

Neighborhood::Neighborhood(char *seq, short *pt)
{
  this->pt = allocopy(pt);
  this->seq = (char*)malloc((strlen(seq)+1)*sizeof(char));
  strcpy(this->seq, seq);

  make_pair_matrix();
  update_fold_params();

  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);

  energy = INT_MAX;

   // create array of loops:
  loops.resize(pt[0]+1);
  for (int i=0; i<(int)loops.size(); i++) loops[i] = NULL;

  // generate the external loop
  Loop *newone = new Loop(0, pt[0]);
  loops[0] = newone;
  int i = newone->GenNeighs(seq, pt);
  if (i==-1) return;

  // generate the neighbourhood (inserts)
  for (; i<pt[0]; i++) {
    if (pt[i]>i) {
      Loop *newone = new Loop(i, pt[i]);
      loops[i] = newone;
      int k = newone->GenNeighs(seq, pt);
      if (k!=-1) i = k-1;
      else i = pt[i];
    }
  }
}

Neighborhood::Neighborhood(Neighborhood &second)
{
  this->pt = allocopy(second.pt);

  this->energy = second.energy;
  this->loopnum = second.loopnum;
  this->neighnum = second.neighnum;

  loops.resize(second.loops.size(), NULL);
  for (int i=0; i<(int)second.loops.size(); i++) {
    if (second.loops[i]) loops[i] = new Loop(*second.loops[i]);
  }

}

Neighborhood::~Neighborhood()
{
  free(pt);
  free(seq);
  if (s0) free(s0);
  if (s1) free(s1);
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) {
      delete loops[i];
    }
  }
}

inline int find_enclosing(short *pt, int i)
{
  int beg = i-1;
  for (; beg>0; beg--) {
    if (pt[beg]==0) continue;
    if (pt[beg]>beg) break;
    if (pt[beg]<beg) beg = pt[beg];
  }
  return beg;
}

int Neighborhood::AddBase(int i, int j, bool reeval)
{
  // find enclosing loop (can be better)
  int beg = find_enclosing(pt, i);

  int size = -loops[beg]->neighs.size();

  // insert it + generate new neighbors
  Loop* newloop = new Loop(i,j);
  if (loops[i]) error_message("Loop %3d already set!!!", i);
  loops[i] = newloop;

  newloop->GenNeighs(seq, pt);
  pt[i] = j;
  pt[j] = i;
  if (reeval) newloop->EvalLoop(pt, s0, s1, true);

  // delete the neighbors that are wrong now (can be better)
  loops[beg]->GenNeighs(seq, pt);
  if (reeval) loops[beg]->EvalLoop(pt, s0, s1, true);

  size += loops[i]->neighs.size() + loops[beg]->neighs.size();

  return size;
}

int Neighborhood::RemBase(int i, int j, bool reeval)
{
  // loop exists?
  if (loops[i] == NULL) error_message("There is no loop at point %d!!!\n", i);
  if (loops[i]->left!=i || loops[i]->right!=j) error_message("Different end: removing (%d, %d); exists (%d, %d)\n", i, j, loops[i]->left, loops[j]->right);

  // find the upper loop:
  int upper = find_enclosing(pt, i);
  int size = -loops[upper]->neighs.size() - loops[i]->neighs.size();

  // delete this one
  delete loops[i];
  loops[i] = NULL;
  pt[i] = 0;
  pt[j] = 0;

  // recompute the upper one:
  loops[upper]->GenNeighs(seq, pt);
  if (reeval) loops[upper]->EvalLoop(pt, s0, s1, true);
  size += loops[upper]->neighs.size();

  return size;
}

int Neighborhood::ApplyNeigh(Neigh &neigh, bool reeval)
{
  if (neigh.i>0) return AddBase(neigh.i, neigh.j, reeval);
  else return RemBase(-neigh.i, -neigh.j, reeval);
}

int Neighborhood::PrintNeighs()
{
  int res = 0;
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) {
      if (i!=0) res++; // one delete move per loop
      fprintf(stdout, "Loop %3d %3d - %5d (%d neighbors):\n", loops[i]->left, loops[i]->right, loops[i]->energy, (int)loops[i]->neighs.size());
      for (int j=0; j<(int)loops[i]->neighs.size(); j++) {
        fprintf(stdout, "  %3d %3d %5d\n", loops[i]->neighs[j].i, loops[i]->neighs[j].j, loops[i]->neighs[j].energy_change);
        res++;
      }
    }
  }
  return res;
}

int Neighborhood::PrintEnum()
{
  int res = 0;
  Neigh tmp(0,0,0);
  StartEnumerating();
  while (NextNeighbor(tmp, true)) {
    fprintf(stdout, "  %3d %3d %5d\n", tmp.i, tmp.j, tmp.energy_change);
  }
  return res;
}

void Neighborhood::PrintStr()
{
  fprintf(stdout, "%s %6.2f\n", pt_to_str(pt).c_str(), energy/100.0);
}

int Neighborhood::EvalNeighs(bool full)
{
  energy = 0;
  for (int i=0; i<(int)loops.size(); i++) {
    if (loops[i]) energy += loops[i]->EvalLoop(pt, s0, s1, full);
  }

  return energy;
}

void addToLowest(std::vector<Neigh> &low, Neigh &lowest, bool clear)
{
  if (clear) low.clear();
  low.push_back(lowest);
}

int Neighborhood::RemEnergy(short *pt, int loop, int last_loop)
{
  // find last loop if not provided
  if (last_loop == -1) {
    for (last_loop = loop-1; last_loop>0; last_loop--) {
      if (pt[last_loop] == 0) continue;
      if (pt[last_loop] > last_loop) break;
    }
  }

  // resolve energy:
  pt[loops[loop]->left] = 0;
  pt[loops[loop]->right] = 0;
  int change = -loops[loop]->energy - loops[last_loop]->energy + loop_energy(pt, s0, s1, loops[last_loop]->left);
  pt[loops[loop]->left] = loops[loop]->right;
  pt[loops[loop]->right] = loops[loop]->left;

  return change;
}

int Neighborhood::MoveLowest(bool reeval)
{
  int lowest = 0;
  StartEnumerating();
  Neigh next;
  Neigh lowest_n;
  while (NextNeighbor(next, true)) {
    if (next.energy_change == lowest) {
      AddDegen(next);
    }
    if (next.energy_change < lowest) {
      ClearDegen();
      lowest = next.energy_change;
      lowest_n = next;
    }
  }

  // resolve degeneracy
  if (degen_todo.size() > 0) {
    degen_done.push_back(*this);
    for (int i=0; i<(int)degen_todo.size(); i++) {
      int degen_en = degen_todo[i].MoveLowest(reeval);
      if (degen_en < lowest) {
        *this = degen_todo[i];
        ClearDegen();
        return degen_en;
      }
    }

  if (degen_done.size() > 0) {
    // chose the lowest one lexicographically:
    Neighborhood &res = degen_done[i];
    for (int i=1; i<(int)degen_done.size(); i++) {
      if (degen_done[i] < res) res = degen_done[i];
    }
    *this = res;
    ClearDegen();
    return degen_en;
  }

  // apply it:
  if (lowest < 0) {
    ApplyNeigh(lowest_n);
  }

  return lowest;
}


void Neighborhood::StartEnumerating()
{
  loopnum = 0;
  neighnum = 0;
}

bool Neighborhood::NextNeighbor(Neigh &res, bool with_energy)
{
  // first get end:
  if ((int)loops.size() <= loopnum) return false;

  // store last loop:
  int last_loop = -1;

  // deletes and inserts:
  if (neighnum == -1) {
    // deletes
    res = Neigh(-loops[loopnum]->left, -loops[loopnum]->right, with_energy?RemEnergy(pt, loopnum, last_loop):INT_MAX);
  } else {
    // inserts
    res = loops[loopnum]->neighs[neighnum];
  }

  // increase the count
  neighnum++;
  if (neighnum >= (int)loops[loopnum]->neighs.size()) {
    neighnum = -1;
    last_loop = loopnum;
    loopnum++;
    while (loopnum<(int)loops.size() && loops[loopnum]==NULL) loopnum++;
  }

  //if (loop.size() <= loopnum) return false;
  return true;
}

void Neighborhood::ClearDegen()
{
  degen_done.clear();
  degen_todo.clear();
}

bool Neighborhood::AddDegen(Neigh &neigh)
{
  int res = false;

  // check if already there:
  ApplyNeigh(neigh);
  for (int i=0; i<(int)degen_todo.size(); i++) {
    if (degen_todo[i] == *this) {
      res = true;
      break;
    }
  }
  if (!res) {
    for (int i=0; i<(int)degen_done.size(); i++) {
      if (degen_done[i] == *this) {
        res = true;
        break;
      }
    }
  }
  if (!res)  {
    degen_todo.push_back(*this);
  }
  ApplyNeigh(Neigh(-neigh.i; -neigh.j));

  return res;
}

extern "C" {
  #include "utils.h"
}

void test()
{
//char num[] = "123456789012345678901234567890123456";
  char seq[] = "CCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGG";
  char str[] = "....................................";
  short *pt = make_pair_table(str);

  Neighborhood nh(seq, pt);
  nh.EvalNeighs(true);
  nh.PrintNeighs();

  int size = nh.AddBase(15,25,true);
  fprintf(stderr, "adding %d %d added %d neighbours\n", 15,25,size);
  nh.PrintNeighs();
  nh.PrintEnum();
  nh.PrintStr();
  nh.MoveLowest();
  nh.PrintStr();

  free(pt);
  free_arrays();
}
