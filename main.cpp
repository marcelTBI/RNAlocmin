#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
  #include "hash_util.h"
}

#include "move_set.h"
#include "flood.h"

#include "barrier_tree.h"

extern short* allocopy(short *src);
extern void copy_arr(short *desc, short *src);

using namespace std;

// are arguments bad?
int bad_arguments(gengetopt_args_info &args_info)
{
/*option "move"               m "Move set:\nI ==> insertion & deletion of base pair\nS ==> I&D& switch base pair" values="I","S" default="I" no
option "min-num"            n "Maximal number of local minima returned" int default="100" no
option "find-num"           - "Maximal number of local minima found \n  (default = unlimited - crawl through whole file)" int no
option "seq"                s "Sequence file in FASTA format" string default="seq.txt"
option "verbose-lvl"        v "Level of verbosity (0 = nothing, 3 = full)" int default="0" no
option "rates"              r "Create rates for treekin" flag off
option "rates-file"         f "File where to write rates" string default="rates.out" no
option "temp"               T "Temperature in Celsius (only for rates)" double default="37.0" no
option "depth"              d "Depth of findpath search (higher values increase running time)" int default="10" no
option "minh"               - "Print only minima with energy barrier greater than this" double default="0.0" no
option "noLP"               - "Work with canonical RNA structures (w/o isolated base pairs)" flag off
option "bartree"            b "Generate possible barrier tree" flag off
option "useEOS"             e "Use energy_of_structure_pt calculation instead of energy_of_move (slower, it should not affect results)" flag off
option "useFirst"           - "Use first found lower energy structure instead of deepest" flag off
option "floodPortion"       - "Fraction of minima to flood\n(0.0 -> no flood; 1.0 -> try to flood all of them)" double default="0.95" no*/

  if (args_info.min_num_arg<=0) {
    fprintf(stderr, "Number of local minima should be positive integer (min-num)\\n");
    return -1;
  }

  if (args_info.find_num_given && args_info.find_num_arg<=0) {
    fprintf(stderr, "Number of local minima should be positive integer (find-num)\n");
    return -1;
  }

  if (args_info.verbose_lvl_arg<0 || args_info.verbose_lvl_arg>4) {
    if (args_info.verbose_lvl_arg<0) args_info.verbose_lvl_arg = 0;
    else args_info.verbose_lvl_arg = 4;
    fprintf(stderr, "WARNING: level of verbosity is not in range (0-4), setting it to %d\n", args_info.verbose_lvl_arg);
  }

  if (args_info.temp_arg<-273.15) {
    fprintf(stderr, "Temperature cannot be below absolute zero\n");
    return -1;
  }

  if (args_info.floodPortion_arg<0.0 || args_info.floodPortion_arg>1.0) {
    args_info.floodPortion_arg = (args_info.floodPortion_arg<0.0 ? 0.0 : 1.0);
    fprintf(stderr, "WARNING: floodPortion is not in range (0.0-1.0), setting it to %.1f\n", args_info.floodPortion_arg);
    return -1;
  }

  if (args_info.depth_arg<=0) {
    fprintf(stderr, "Depth of findpath search should be positive integer\n");
    return -1;
  }

  if (args_info.minh_arg<0.0) {
    fprintf(stderr, "Depth of findpath search should be non-negative number\n");
    return -1;
  }

  // everything ok
  return 0;
}

static int num_moves = 0;
static int seq_len;

int move(encoded &enc, char *seq, map<hash_entry, int, compare_map> &output, options &opt)
{
  // count moves
  num_moves++;

  // read a structure
  char *structure = my_getline(stdin);
  if (structure == NULL) return -1;
  if (structure[0]=='>') {
    free(structure);
    return 0;
  }
  // find length of structure
  int len=0;
  while (structure[len]!='\0' && structure[len]!=' ') len++;

  if (len!=seq_len) {
    fprintf(stderr, "Unequal lengths:\n(structure) %s\n (sequence) %s\n", structure, seq);
    free(structure);
    return -1;
  }

  // convert to pt
  encode_str(&enc, structure);
  free(structure);

  // was it before?
  hash_entry *hee = (hash_entry*)space(sizeof(hash_entry));
  hee->structure = enc.pt;
  hash_entry *tmp_h = (hash_entry*)lookup_hash(hee);
  if (tmp_h) {
    free(hee);
    tmp_h->count++;
    return 0;
  } else {
    hee->structure = allocopy(enc.pt);
    hee->count = 1;
    write_hash(hee);
  }

  //is it canonical (noLP)
  if (opt.noLP && find_lone_pair(enc.pt)!=-1) {
    fprintf(stderr, "WARNING: structure \"%s\" has lone pairs, skipping...\n", pt_to_str(enc.pt).c_str());
    return -2;
  }

  //debugging
  if (opt.verbose_lvl>2) fprintf(stderr, "processing: %d %s\n", num_moves, pt_to_str(enc.pt).c_str());

  // find energy and set options
  degen deg;
  deg.opt = &opt;

  int energy = energy_of_structure_pt(seq, enc.pt, enc.s0, enc.s1, 0);
  hee->energy = energy;

  // deepest descend
  while (move_set(enc, energy, deg)!=0) {
    erase_set(deg.unprocessed);
    erase_set(deg.processed);
  }
  erase_set(deg.unprocessed);
  erase_set(deg.processed);

  if (opt.verbose_lvl>2) fprintf(stderr, "\n  %s %d\n", pt_to_str(enc.pt).c_str(), energy);

  // save for output
  //string stro = pt_to_str(enc.pt);
  map<hash_entry, int, compare_map>::iterator it;
  hash_entry he;
  he.structure = enc.pt;
  he.energy = energy;
  he.count = 0;
  if ((it = output.find(he)) != output.end()) {
    it->second++;
  } else {
    he.structure = allocopy(enc.pt);
    output.insert(make_pair(he, 1));
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
  if (bad_arguments(args_info) !=0 ) {
    fprintf(stderr, "ERROR: one or more bad arguments, exiting...\n");
    exit(EXIT_FAILURE);
  }

  // adjust options
  options opt;
  opt.minh = args_info.minh_arg;
  opt.noLP = args_info.noLP_flag;
  opt.EOM = !args_info.useEOS_flag;
  opt.first = args_info.useFirst_flag;
  opt.f_point = NULL;
  opt.shift = args_info.move_arg[0]=='S';
  opt.verbose_lvl = args_info.verbose_lvl_arg;

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

  // keep track of structures & statistics
  map<hash_entry, int, compare_map> output; // structures plus energies to output (+ how many hits has this minima)
  int not_canonical = 0;

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Time to initialize: %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // ########################## main loop - reads structures from RNAsubopt and process them
  int count = 0;  //num of local minima
  encoded *enc=NULL;
  enc = encode_seq(seq);
  while (!args_info.find_num_given || count != args_info.find_num_arg) {
    int res = move(*enc, seq, output, opt);
    if (res==0)   continue;
    if (res==-1)  break;
    if (res==-2)  not_canonical++;
    if (res==1)   count=output.size();
  }

  // free hash
  //int num_of_structures = hash_size();
  if (args_info.verbose_lvl_arg>0) print_stats();
  kill_hash();

  // time?
  if (args_info.verbose_lvl_arg>0) {
    fprintf(stderr, "Main loop (deepest descent from RNAsubopt): %.2f secs.\n", (clock() - clck1)/(double)CLOCKS_PER_SEC);
    clck1 = clock();
  }

  // vectors for rates computation
  int num = (count > args_info.min_num_arg ? args_info.min_num_arg : count);
  vector<string> output_str;
  output_str.resize(num);
  vector<hash_entry> output_he;
  output_he.resize(num);
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
      if (i<num) {
        output_str[i]=pt_to_str(it->first.structure);
        output_num[i]=it->second;
        output_he[i]=it->first;
        output_en[i]=it->first.energy;
        i++;
      } else { // we have enough minima
        free(it->first.structure);
      }
    }
    output.clear();

    // threshold for flooding
    vector<int> tmp = output_num;
    sort(tmp.begin(), tmp.end());
    int thr = num*args_info.floodPortion_arg;
    thr--;
    threshold = (thr<0 ? 0 : tmp[thr]);
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
    }

    int flooded = 0;
    // first try to flood the highest bins
    for (int i=num-1; i>=0; i--) {
      // flood only if low number of walks ended there
      if (output_num[i]<=threshold) {
        copy_arr(enc->pt, output_he[i].structure);
        fprintf(stderr,   "flooding  (%3d): %s %.2f\n", i, output_str[i].c_str(), output_he[i].energy/100.0);

        int saddle;
        hash_entry *he = flood(*enc, output_he[i].energy, opt, saddle);

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
          // setup move_set
          degen deg;
          deg.opt = &opt;
          copy_arr(enc->pt, he->structure);
          int en = he->energy;
          free(he->structure);
          free(he);
          // walk down
          while (move_set(*enc, en, deg)!=0) {
            erase_set(deg.processed);
            erase_set(deg.unprocessed);
          };
          erase_set(deg.processed);
          erase_set(deg.unprocessed);

          // now check if we have the minimum already (hopefuly yes ;-) )
          hash_entry he_tmp;
          he_tmp.structure = enc->pt;
          he_tmp.energy = en;
          vector<hash_entry>::iterator it;
          it = lower_bound(output_he.begin(), output_he.end(), he_tmp, compare_vect);

          if (args_info.verbose_lvl_arg>1) fprintf(stderr, "minimum: %s %.2f\n", pt_to_str(he_tmp.structure).c_str(), he_tmp.energy/100.0);

          if (it!=output_he.end()) {
            int pos = (int)(it-output_he.begin());
            fprintf(stderr, "found father at pos: %d\n", pos);

            // which one is father?
            bool res_higher = compare_vect(output_he[i], output_he[pos]); // is found minima higher in energy than flooded min?
            int father = (res_higher ? i:pos);
            int child = (res_higher ? pos:i);

            nodes[i].saddle_height = saddle/100.0; // this is always true

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
              add_father(nodes, chng_child, chng_father);
            }

            flooded++;
            energy_barr[i*num+pos] = energy_barr[pos*num+i] = saddle/100.0;
          } else {
            fprintf(stderr, "minimum not found: %s %.2f\n", pt_to_str(he_tmp.structure).c_str(), he_tmp.energy/100.0);
          }
        } else {
          fprintf(stderr, "flood unsuccesful!\n");
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

    // make tree (fill missing nodes)
    make_tree(num, energy_barr, nodes);

    // plot it!
    PS_tree_plot(nodes, num, "treeLocmin.ps");
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

/*
  // print structures + count (sorted)
  if (args_info.verbose_lvl_arg>0) {
    vector<pair<int, string> > out;
    for (map<string, int>::iterator it=structures.begin(); it!=structures.end();it++) {
      out.push_back(make_pair(it->second, it->first));
      //fprintf(stderr, "%s %d\n", it->first.c_str(), it->second);
    }
    if (args_info.noLP_flag) {
      fprintf(stderr, "move_set called %d times on %d canonical structures (%d total structures).\n", count_move(), (int)out.size()-not_canonical, (int)out.size());
    } else {
      fprintf(stderr, "move_set called %d times on %d structures.\n", count_move(), (int)out.size());
    }
    if (args_info.verbose_lvl_arg>2) {
      sort(out.begin(), out.end());
      for(int i=out.size()-1; i>=0; i--) {
        fprintf(stderr, "%s %d\n", out[i].second.c_str(), out[i].first);
      }
      fprintf(stderr, "\n");
    }
  }
*/
/*  // print mean energy + entropy of generation
  if (args_info.verbose_lvl_arg>0) {
    double ME = 0.0;
    double Z = 0.0;
    double E = 0.0;

    double meanE = 0.0;
    double entropy = 0.0;

    int all = 0;
    double _kT = 0.00198717*(273.15 + args_info.temp_arg);

    for (map<string, int>::iterator it=structures.begin(); it!=structures.end(); it++) {
      float energy = energy_of_structure(seq, it->first.c_str(), 0);
      all += it->second;
      Z += exp(-energy/_kT);
    }
    for (map<string, int>::iterator it=structures.begin(); it!=structures.end(); it++) {
      float energy = energy_of_structure(seq, it->first.c_str(), 0);
      double expe = exp(-energy/_kT);

      meanE += energy*(it->second);
      entropy -= it->second/(double)all*(log(it->second/(double) all));

      ME += energy*expe;


      E += expe/Z*log(exp(energy/_kT)/Z);
    }

    meanE /= (double)all;
    ME /= Z;

    fprintf(stderr, "mean energy : %lf\n", meanE);
    fprintf(stderr, "entropy     : %lf\n", entropy);
    fprintf(stderr, "minima found: %d\n", count);
  }*/

  // release resources
  if (energy_barr!=NULL) free(energy_barr);
  if (seq!=NULL) free(seq);
  if (name!=NULL) free(name);
  if (enc) free_encode(enc);
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