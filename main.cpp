#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>

extern "C" {
  #include "fold.h"
  #include "findpath.h"
  #include "RNAlocmin_cmdline.h"
  #include "utils.h"
  #include "read_epars.h"
  #include "fold_vars.h"
}

#include "hash_util.h"
#include "globals.h"
#include "move_set.h"
#include "flood.h"

#include "barrier_tree.h"

using namespace std;

inline bool isStruct(char *p)
{
  // check first two chars - should be enough
  if ((p[0]=='.' || p[0]=='(' || p[0]==')') && (p[1]=='.' || p[1]=='(' || p[1]==')')) return true;
  else return false;
}

inline bool isEnergy(char *p, float &energy)
{
  if (sscanf(p, "%f", &energy) == 1 && energy<100.0 && energy >-200.0) return true;
  else return false;
}

int move(unordered_map<hash_entry, int, hash_fncts> &structs, map<hash_entry, int, compare_map> &output, set<hash_entry, compare_map> &output_shallow)
{
  // count moves
  num_moves++;

  // read a line
  char *line = my_getline(stdin);
  if (line == NULL) return -1;
  if (line[0]=='>') {
    free(line);
    return 0;
  }

  float energy=1e10;

  // process lines
  char *p = line;
  char *sep = " \t\n";
  char *temp;

  bool struct_found = false;
  bool energy_found = false;

  p = strtok(line, sep);
  while(p!=NULL && !(struct_found && energy_found)) {
    //fprintf(stderr, "%s\n", p);
    if (isStruct(p)) {
      if (struct_found) fprintf(stderr, "On line \"%s\" two structure-like sequences found!\n", line);
      else {
        temp = p;
      }
      struct_found = true;
    } else {
      if (isEnergy(p, energy)) {
        energy_found = true;
      }
    }

    p = strtok(NULL, sep);
  }
  p = temp;
  if (!struct_found) {
    free(line);
    return 0;
  }

  // find length of structure
  int len=0;
  while (p[len]!='\0' && p[len]!=' ') len++;

  if (len!=seq_len) {
    fprintf(stderr, "Unequal lengths:\n(structure) %s\n (sequence) %s\n", p, Enc.seq);
    free(line);
    return -1;
  }

  // was it before?
  hash_entry str;
  str.structure = Enc.Struct(p);
  free(line);
  unordered_map<hash_entry, int, hash_fncts>::iterator it_s = structs.find(str);
  //hash_entry *tmp_h = (hash_entry*)lookup_hash(str);
  // if it was - release memory + get another
  if (it_s != structs.end()) {
    it_s->second++;
    free(str.structure);
    return 0;
  } else {
    // find energy only if not in input (not working - does energy_of_move require energy_of_struct run first???)
    str.energy = Enc.Energy(str);
    /*if (1 || !energy_found) str.energy = Enc.Energy(str);
    else str.energy = (int)(energy*100.0+(energy<0.0 ? -0.5 : 0.5));*/

    //if (en != str.energy) fprintf(stderr, "%d %d\n", en, str.energy);

    // insert into hash
    structs[str] = 1;
    str.structure = allocopy(str.structure);
  }

  //is it canonical (noLP)
  if (Opt.noLP && find_lone_pair(str.structure)!=-1) {
    fprintf(stderr, "WARNING: structure \"%s\" has lone pairs, skipping...\n", pt_to_str(str.structure).c_str());
    return -2;
  }

  //debugging
  if (Opt.verbose_lvl>2) fprintf(stderr, "processing: %d %s\n", num_moves, pt_to_str(str.structure).c_str());

  // descend
  int i;
  while ((i = (Opt.rand? move_rand(str) : move_set(str)))!=0) {
    Deg.Clear();
  }
  Deg.Clear();

  // discard shallow ones (check if we have them already)
  if (Opt.minhall && Opt.minh>0 && (output.find(str) == output.end()) && (output_shallow.find(str) == output_shallow.end())) {

    // flood it a little
    int status;
    hash_entry *escape = flood(str, status, Opt.minh);

    // while not found non-shallow one
    while (escape) {

      //add to shallow set
      if (output_shallow.find(str)!=output_shallow.end()) break;
      output_shallow.insert(str);
      str.structure = allocopy(str.structure);
      if (Opt.verbose_lvl>1) fprintf(stderr, "shallow:  %s %d\n", pt_to_str(str.structure).c_str(), str.energy);

      // find another minima
      while (move_set(*escape)!=0) {
        Deg.Clear();
      }
      Deg.Clear();

      // try to flood again
      hash_entry *he_tmp = flood(*escape, status, Opt.minh);
      free_entry(escape);
      escape = he_tmp;
    }
  }

  if (Opt.verbose_lvl>2) fprintf(stderr, "\n  %s %d\n", pt_to_str(str.structure).c_str(), str.energy);

  // save for output
  map<hash_entry, int, compare_map>::iterator it;
  if ((it = output.find(str)) != output.end()) {
    it->second++;
    free(str.structure);
  } else {
    str.num = output.size();
    output.insert(make_pair(str, 1));
  }

  return 1;
}

int main(int argc, char **argv)
{
  clock_t clck1 = clock();

  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "ERROR: argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  //check for bad arguments
  if (Opt.Init(args_info) !=0 ) {
    fprintf(stderr, "ERROR: one or more bad arguments, exiting...\n");
    exit(EXIT_FAILURE);
  }

  // read parameter file
  if (args_info.paramFile_given) {
    read_parameter_file(args_info.paramFile_arg);
  }

  // dangle setup
  if (args_info.dangles_given) {
    dangles = args_info.dangles_arg;
  }

  // read sequence
  FILE *fseq;
  fseq = fopen(args_info.seq_arg, "r");
  if (fseq == NULL) {
    fprintf(stderr, "Cannot open file \"%s\".\n", args_info.seq_arg);
    exit(EXIT_FAILURE);
  }
  char *seq;
  char *name =NULL;
  name = my_getline(fseq);
  if (name == NULL) {
    fprintf(stderr, "File \"%s\" empty.\n", args_info.seq_arg);
    exit(EXIT_FAILURE);
  }
  seq = my_getline(fseq);
  if (seq == NULL || name[0]!='>') {
    //fprintf(stderr, "WARNING: File \"%s\" not in FASTA format. Using it as a sequence.\n", args_info.seq_arg);
    if (seq!=NULL) free(seq);
    seq = name;
    name = NULL;
  }
  fclose(fseq);
  if (args_info.verbose_lvl_arg>1) fprintf(stderr, "%s\n", seq);

  seq_len = strlen(seq);

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Time to initialize: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // ########################## main loop - reads structures from RNAsubopt and process them
  int count = 0;  //num of local minima
  Enc.Init(seq);

  // keep track of structures & statistics
  map<hash_entry, int, compare_map> output; // structures plus energies to output (+ how many hits has this minima)
  set<hash_entry, compare_map> output_shallow; // shallow structures (if minh specified)
  int not_canonical = 0;
  // hash
  unordered_map<hash_entry, int, hash_fncts> structs (HASHSIZE);
  while (!args_info.find_num_given || count != args_info.find_num_arg) {
    int res = move(structs, output, output_shallow);
    if (res==0)   continue;
    if (res==-1)  break;
    if (res==-2)  not_canonical++;
    if (res==1)   count=output.size();
  }

  // free hash
  //int num_of_structures = hash_size();
  if (args_info.verbose_lvl_arg>0) print_stats(structs);
  free_hash(structs);

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Main loop (deepest descent from RNAsubopt): %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // vectors for rates computation
  int num = ((count > args_info.min_num_arg && args_info.min_num_arg!=0) ? args_info.min_num_arg : count);
  vector<string> output_str;
  output_str.resize(num);
  vector<hash_entry> output_he;
  hash_entry h;
  h.structure = NULL;
  output_he.resize(num, h);
  vector<int> output_en;
  output_en.resize(num);
  vector<int> output_num;
  output_num.resize(num);

  // threshold for flooding
  int threshold;

  // insert structures into vectors (insert only <num> minima)
  {
    int i=0;
    for (map<hash_entry, int, compare_map>::iterator it=output.begin(); it!=output.end(); it++) {
      // if not enough minima
      if (i<num) {
        // first check if the output is not shallow
        if (!Opt.minhall && Opt.minh>0 && i<num) {
          int saddle;
          hash_entry *escape = flood(it->first, saddle, Opt.minh);

          // shallow
          if (escape) {
            if (args_info.verbose_lvl_arg>0) {
              fprintf(stderr, "shallow: %s %6.2f\n", pt_to_str(it->first.structure).c_str(), it->first.energy/100.0);
            }
            free_entry(escape);
            free(it->first.structure);
            continue;
          }
        }
        // then add it to outputs.
        if (args_info.noSort_flag) {
          if (it->first.num<num) {
            output_str[it->first.num]=pt_to_str(it->first.structure);
            output_num[it->first.num]=it->second;
            output_he[it->first.num]=it->first;
            output_en[it->first.num]=it->first.energy;
            //printf("%d ", it->first.num);
            i++;
          }
        } else {
          output_str[i]=pt_to_str(it->first.structure);
          output_num[i]=it->second;
          output_he[i]=it->first;
          output_en[i]=it->first.energy;
          i++;
        }
      } else { // we have enough minima
        free(it->first.structure);
      }
    }
    output.clear();

    // erase possible NULL elements...
    if (args_info.noSort_flag) {
      for (int j=output_he.size()-1; j>=0; j--) {
        if (output_he[j].structure==NULL) {
          output_he.erase(output_he.begin()+j);
          output_str.erase(output_str.begin()+j);
          output_en.erase(output_en.begin()+j);
          output_num.erase(output_num.begin()+j);
        }
      }
    }

    // threshold for flooding
    vector<int> tmp = output_num;
    sort(tmp.begin(), tmp.end());
    int thr = num*args_info.floodPortion_arg;
    thr--;
    threshold = (thr<0 ? 0 : tmp[thr]);
  }

    // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Discarding shallow minima: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // array of energy barriers
  float *energy_barr = NULL;

  // find saddles - fill energy barriers
  nodeT nodes[num];
  if (args_info.rates_flag || args_info.bartree_flag) {
    energy_barr = (float*) malloc(num*num*sizeof(float));
    for (int i=0; i<num*num; i++) energy_barr[i]=1e10;

    // fill nodes
    for (int i=0; i<num; i++) {
      nodes[i].father = -1;
      nodes[i].height = output_en[i]/100.0;
      nodes[i].label = NULL;
      nodes[i].color = 0.0;
      nodes[i].saddle_height = 1e10;
    }

    int flooded = 0;
    // first try to flood the highest bins
    for (int i=num-1; i>=0; i--) {
      // flood only if low number of walks ended there
      if (output_num[i]<=threshold) {
        //copy_arr(Enc.pt, output_he[i].structure);
        if (args_info.verbose_lvl_arg>2) fprintf(stderr,   "flooding  (%3d): %s %.2f\n", i, output_str[i].c_str(), output_he[i].energy/100.0);

        int saddle;
        hash_entry *he = flood(output_he[i], saddle);

        // print info
        if (args_info.verbose_lvl_arg>1) {
          if (he) {
            fprintf(stderr, "below     (%3d): %s %.2f\n"
                            "en: %7.2f  is: %s %.2f\n", i,
                    output_str[i].c_str(), output_he[i].energy/100.0, saddle/100.0,
                    pt_to_str(he->structure).c_str(), he->energy/100.0);
          } else {
            fprintf(stderr, "unsucesful(%3d): %s %.2f\n", i,
                    output_str[i].c_str(), output_he[i].energy/100.0);
          }
        }

        // if flood succesfull - walk down to find father minima
        if (he) {
          // walk down
          while (move_set(*he)!=0) {
            Deg.Clear();
          };
          Deg.Clear();

          // now check if we have the minimum already (hopefuly yes ;-) )
          vector<hash_entry>::iterator it;
          it = lower_bound(output_he.begin(), output_he.end(), *he, compare_vect);

          if (args_info.verbose_lvl_arg>1) fprintf(stderr, "minimum: %s %.2f\n", pt_to_str(he->structure).c_str(), he->energy/100.0);
          // we dont need he again
          free_entry(he);

          if (it!=output_he.end()) {
            int pos = (int)(it-output_he.begin());
            if (args_info.verbose_lvl_arg>1) fprintf(stderr, "found father at pos: %d\n", pos);

            // which one is father?
            bool res_higher = compare_vect(output_he[i], output_he[pos]); // is found minima higher in energy than flooded min?
            int father = (res_higher ? i:pos);
            int child = (res_higher ? pos:i);

            nodes[i].saddle_height = saddle/100.0; // this is always true
            nodes[pos].saddle_height = (nodes[pos].saddle_height > saddle/100.0 ? saddle/100 : nodes[pos].saddle_height);
            //nodes[i].color = 0.5;

            // some father issues
            int chng_child, chng_father;
            if (nodes[child].father == -1 || nodes[child].father == father) {
              chng_child = child;
              chng_father = father;
            } else { // going to change a father (must be sure we are doing it right)
              int old_father = nodes[child].father;
              if (nodes[old_father].height > nodes[father].height) {
                chng_child = old_father;
                chng_father = father;
              } else {
                chng_child = father;
                chng_father = old_father;
              }
            }

            // add new one and recompute existing fathers
            if (nodes[chng_child].father != father) {
              add_father(nodes, chng_child, chng_father, 0.0);
            }

            flooded++;
            energy_barr[i*num+pos] = energy_barr[pos*num+i] = saddle/100.0;
          }
        }
      }
    }

    // time?
    if (args_info.verbose_lvl_arg>0) {
      fprintf(stderr, "Flood(%d(%d)/%d): %.2f secs.\n", flooded, (int)(num*args_info.floodPortion_arg), num, (clock() - clck1)/(double)CLOCKS_PER_SEC);
      if (args_info.verbose_lvl_arg>1) {
        fprintf(stderr, "Minima left to findpath (their father = -1): ");
        for (int i=0; i<num; i++) {
          if (nodes[i].father == -1) fprintf(stderr, "%d ", i);
        }
        fprintf(stderr, "\n");
      }
      clck1 = clock();
    }

    // for others, just do findpath.
    int findpath = 0;
    for (int i=0; i<num; i++) {
      if (nodes[i].father != -1) continue;

      for (int j=0; j<num; j++) {

        // check if we already
        if (energy_barr[i*num+j]<1e9) continue;
        if (nodes[j].father!=-1) continue;

        if (i>j) {
          energy_barr[i*num+j] = energy_barr[j*num+i]; // assume symetry
        } else if (i==j) {
          energy_barr[i*num+j] = 0.0;
        } else {
          // find (maybe suboptimal) saddle energy
          energy_barr[i*num+j] = find_saddle(seq, output_str[i].c_str(), output_str[j].c_str(), args_info.depth_arg)/100.0;
          findpath++;
        }
      }
    }
    // debug output
    if (args_info.verbose_lvl_arg>2) {
      fprintf(stderr, "Energy barriers:\n");
      //bool symmetric = true;
      for (int i=0; i<num; i++) {
        for (int j=0; j<num; j++) {
          fprintf(stderr, "%6.2f ", energy_barr[i*num+j]);
          //if (energy_barr[i*num+j] != energy_barr[j*num+i]) symmetric = false;
        }
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
      //fprintf(stderr, "%s", (symmetric? "":"non-symmetric energy barriers!!\n"));
    }

    // time?
    if (args_info.verbose_lvl_arg>0) {
      fprintf(stderr, "Findpath(%d/%d): %.2f secs.\n", findpath, num*(num-1)/2, (clock() - clck1)/(double)CLOCKS_PER_SEC);
      clck1 = clock();
    }
  }

  // create rates for treekin
  if (args_info.rates_flag) {
    print_rates(args_info.rates_file_arg, args_info.temp_arg, num, energy_barr, output_en);
  }

  // generate barrier tree?
  if (args_info.bartree_flag) {

    //PS_tree_plot(nodes, num, "tst.ps");

    // make tree (fill missing nodes)
    make_tree(num, energy_barr, nodes);

    // plot it!
    PS_tree_plot(nodes, num, args_info.barr_name_arg);
  }

    // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Rates + barrier tree generation: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // printf output
  printf("     %s\n", seq);
  for (int i=0; i<num; i++) {
    printf("%4d %s %6.2f %6d", i+1, output_str[i].c_str(), output_en[i]/100.0, output_num[i]);
    if (args_info.bartree_flag) {
      printf(" %4d %6.2f\n", nodes[i].father+1, nodes[i].saddle_height-nodes[i].height);
    } else printf("\n");
  }

  // Jing's visualisation
  if (args_info.numIntervals_arg>0) {
    // whole range of our minima
    int range = output_en[num-1] - output_en[0];
    float dE = range/(float)args_info.numIntervals_arg;

    int max_samples = 0;
    for (int i=0; i<num; i++) {
      if (output_num[i]>max_samples) max_samples = output_num[i];
    }

    // function f(i,j)
    int visual[args_info.numIntervals_arg][max_samples];
    for (int i=0; i<max_samples; i++) {
      for (int j=0; j<args_info.numIntervals_arg; j++) {
        visual[j][i]=0;
      }
    }
    // fill the array
    int curr_en = 0;
    for (int i=0; i<num; i++) {
      float curr_range = (float)(output_en[i]-output_en[0]);
      while (curr_range>(curr_en+1)*dE) curr_en++;
      visual[curr_en][output_num[i]-1]++;
    }

    FILE *file;
    char filename[] = "visuals.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
      fprintf(stderr, "Cannot open file \"%s\"\n", filename);
    } else {

      // header
      fprintf(file, " #samp " );
      for (int j=0; j<args_info.numIntervals_arg; j++) fprintf(file, "%6.1f", (output_en[0]+(j+1)*dE)/100.0);
      fprintf(file, "\n\n");

      // data
      for (int i=0; i<max_samples; i++) {
        fprintf(file, "%6d", i+1);
        for (int j=0; j<args_info.numIntervals_arg; j++) {
          fprintf(file, "%6d", visual[j][i]);
        }
        fprintf(file, "\n");
      }
      fclose(file);
    }
  }


  // release resources
  if (energy_barr!=NULL) free(energy_barr);
  if (seq!=NULL) free(seq);
  if (name!=NULL) free(name);
  cmdline_parser_free(&args_info);
  for(unsigned int i=0; i<output_he.size(); i++) {
    free(output_he[i].structure);
  }

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Printing results + freeing args: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  return 0;
}
