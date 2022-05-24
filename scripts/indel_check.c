#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "../klib/kstring.h"
#include "../klib/khash.h"
#include "../klib/kvec.h"

// to compile
// gcc -Wall -g -O2 indel_check.c ../klib/kstring.c  -o indel_check
// int exclamationCheck = strchr(str, '!') != NULL;
struct cigar_info {
    kvec_t(char) vop; // M, S etc
    kvec_t(int) vlen; // length
};

struct cigar_info split_cigar (char *cigar) {
  struct cigar_info r;// = {vop, vlen };
  kv_init(r.vop);
  kv_init(r.vlen);
  char *numbers = "0123456789";
  char tmp[6] = "";
  int n = 0;
  for (int i = 0; i < strlen(cigar); i++) {
    int isNum = strchr(numbers, cigar[i]) != NULL;
    if (isNum) {
        tmp[n] = cigar[i];
        n++;
    } else {
      kv_push(char, r.vop, cigar[i]); // append
      kv_push(int, r.vlen, atoi(tmp));
      // tmp[0] = 0; // wrong way
      memset(tmp, 0, sizeof(tmp)); // reset tmp
      n = 0;
    }
  }

  return r;
}

// fn parse_cigar (cigar: &str, ref_pos: isize, same_strand: bool, read_len: isize)
int * parse_cigar(char *cigar, int ref_pos, int same_strand, size_t read_len){
  struct cigar_info r = split_cigar(cigar);
  int read_pos1 = 0; // left start
  int read_pos2 = -1; // right end
  int ref_pos1 = ref_pos; // left start
  int ref_pos2 = ref_pos - 1; // right end
  int nmatch = 0; // number of M, if match showed up, then no more S or H
  for (int i =0; i < r.vop.n; i++) {
    int num = r.vlen.a[i];
    char op = r.vop.a[i];
    if (op == 'M' || op == '=' || op == 'X'){
      read_pos2 += num;
      ref_pos2 += num;
      nmatch += 1;
    } else if (op == 'S' || op == 'H') {
      if (nmatch == 0) {
        read_pos1 += num;
        read_pos2 += num;
      }
    } else if (op == 'I'){
      read_pos2 += num;
    } else if (op == 'D' || op == 'N') {
      ref_pos2 += num - 1;
    }
  }
  static int rr[4];
  if (same_strand) {
    rr[0] = read_pos1; rr[1]=read_pos2; rr[2]= ref_pos1; rr[3] = ref_pos2;
  } else {
    rr[0] = read_len - read_pos2 - 1; rr[1] = read_len - read_pos1 - 1; rr[2] = ref_pos1; rr[3] = ref_pos2;
  }
  // destroy
  kv_destroy(r.vop);
  kv_destroy(r.vlen);

  return rr;
}

// string slice
void slice(const char *str, char *result, size_t start, size_t end)
{
    strncpy(result, str + start, end - start);
    result[end - start] = 0;
}

KHASH_MAP_INIT_STR(str, int)      // instantiate structs and methods
// fn parse_line (line: &str, map: &mut HashMap<String, isize>, no_small_indels: bool, debug: bool) 
int parse_line(kstring_t *ks, khash_t(str) *h, int debug){
  int *ff, n, i;
  // char *line = ks->s;
  ff = ksplit(ks, '\t', &n);
  char *read_id = ks->s + ff[0];
  int flag = atoi(ks->s + ff[1]);
  char *chrom = ks->s + ff[2];
  int pos  = atoi(ks->s + ff[3]) - 1; // 0 based
  char *cigar  = ks->s + ff[5];
  char *read_seq= ks->s + ff[9];
  size_t read_len = strlen(read_seq);
  char *strand = flag & 0x10 ? "-" : "+";
  // let read_id = ff[0];
  if (strcmp(cigar, "*") == 0 ) { // no mapping
      return 0;
  }
  int has_sa = 0;
  for (i = 11; i < n; ++i){
    if (strstr(ks->s + ff[i], "SA:Z") != NULL) {
      has_sa = 1;
      break;
    }
  }
  if (has_sa && strstr(cigar, "H") == NULL){
    char *sa_info = ks->s + ff[i];
    char *token = strtok(sa_info, ":"); // first string
    token = strtok(NULL, ":"); // 2nd string
    token = strtok(NULL, ":"); // 3rd string
    // printf("token is %s\n", token);
    kstring_t s = { 0, 0, NULL };
    kputs(token, &s); // string to Kstring
    int *ff2 = ksplit(&s, ',', &n);
    char *sa_chrom  = s.s + ff2[0];
    int sa_pos  = atoi(s.s + ff2[1]) - 1; // 0-based this is close to the border on the left, may need to adjust
    char *sa_strand = s.s + ff2[2];
    char *sa_cigar  = s.s + ff2[3];
    // free: should not free
    // free(s.s); free(ff2);free(ff);
    int all_pos1[4];
    int all_pos2 [4];
    if (strcmp(chrom,sa_chrom)==0 && strcmp(strand, sa_strand) ==0) { // potential big deletion, could be insertion too, but update later
      if (sa_pos > pos) { // SA is on the right
        memcpy(all_pos1, parse_cigar(cigar, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos2, parse_cigar(sa_cigar, sa_pos, 1, read_len), 4 * sizeof(int));
      } else {
        memcpy(all_pos2, parse_cigar(cigar, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos1, parse_cigar(sa_cigar, sa_pos, 1, read_len), 4 * sizeof(int));
      }
      if (debug){
        printf("potential big deletion\n%s\t%s\n", read_id, chrom);
        printf("all_pos1: [%d, %d, %d, %d]\n", all_pos1[0], all_pos1[1], all_pos1[2], all_pos1[3]);
        printf("all_pos2: [%d, %d, %d, %d]\n", all_pos2[0], all_pos2[1], all_pos2[2], all_pos2[3]);
      }
      if (all_pos1[0] > all_pos2[1] || all_pos1[3] >= all_pos2[2] || all_pos1[1] >= all_pos2[1]) {return 0;}
      int read_pos1 = all_pos1[1];
      int ref_pos1  = all_pos1[3];
      int read_pos2 = all_pos2[0];
      int ref_pos2  = all_pos2[2];
      int shift = read_pos1 >= read_pos2 ? read_pos1 - read_pos2 + 1 : 0;
      int del_end_pos = ref_pos2 + shift;
      if (ref_pos1 < del_end_pos) {
        if (read_pos2+shift+1 > read_len) {
            return 0;
        }
        char alt_seq[read_pos2+shift+1-read_pos1];
        slice(read_seq, alt_seq, read_pos1, read_pos2+shift+1);
        // printf("alt_seq is %s\n", alt_seq);
        int mut_size = del_end_pos - ref_pos1 - 1;
        kstring_t kk = { 0, 0, NULL };
        ksprintf(&kk, "%s\t%d\t%d\t%s\t%d\tbig_indel", chrom, ref_pos1+1, del_end_pos+1, alt_seq, mut_size);
        // printf("key is %s\n", kk.s);
        khint_t k;
        int absent;
        k = kh_put(str, h, kk.s, &absent);
        if (!absent) {
          kh_value(h, k) += 1; // set the value
        } else {
          kh_value(h, k) = 1;
        }
      }
    }else if (strcmp(chrom,sa_chrom)==0 && strcmp(strand, sa_strand) !=0) {// inversions
      if (sa_pos > pos) { // SA is on the right
        memcpy(all_pos1, parse_cigar(cigar, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos2, parse_cigar(sa_cigar, sa_pos, 0, read_len), 4 * sizeof(int));
      } else {
        memcpy(all_pos2, parse_cigar(cigar, pos, 0, read_len), 4 * sizeof(int));
        memcpy(all_pos1, parse_cigar(sa_cigar, sa_pos, 1, read_len), 4 * sizeof(int));
      }
      if (debug){
        printf("potential inversion\n%s\t%s\n", read_id, chrom);
        printf("all_pos1: [%d, %d, %d, %d]\n", all_pos1[0], all_pos1[1], all_pos1[2], all_pos1[3]);
        printf("all_pos2: [%d, %d, %d, %d]\n", all_pos2[0], all_pos2[1], all_pos2[2], all_pos2[3]);
      }
      int read_pos1 = all_pos1[1];
      int read_pos2 = all_pos2[0];
      int ref_pos1  = all_pos1[3];
      int ref_pos2  = all_pos2[3];
      int shift = read_pos1 >= read_pos2 ? read_pos1 - read_pos2 + 1 : 0;
      int del_end_pos = ref_pos2 - shift;
      if (all_pos1[0] > all_pos2[0]){ // count from right
          read_pos1 = all_pos1[0];
          read_pos2 = all_pos2[1];
          ref_pos1  = all_pos1[2];
          ref_pos2  = all_pos2[2];
          shift = read_pos2 >= read_pos1 ? read_pos2 - read_pos1 + 1 : 0;
          del_end_pos = ref_pos2 + shift;
      }
      if (ref_pos1 <= del_end_pos) {
        int mut_size = del_end_pos - ref_pos1 - 1;
        kstring_t kk = { 0, 0, NULL };
        ksprintf(&kk, "%s\t%d\t%d\tinversion\t%d\tinv", chrom, ref_pos1+1, del_end_pos+1,  mut_size);
        khint_t k;
        int absent;
        k = kh_put(str, h, kk.s, &absent);
        if (!absent) {
          kh_value(h, k) += 1; // set the value
        } else {
          kh_value(h, k) = 1;
        }
      }
    }
  }

  return 0;
}

// main function
int main (int argc, char **argv)
{
  int min_cov = 1;
  int debug = 0;
  int c;

  opterr = 0;

  while ((c = getopt (argc, argv, "dc:")) != -1)
    switch (c)
      {
      case 'd':
        debug = 1;
        break;
      case 'c':
        min_cov = atoi(optarg);
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort ();
      }

  fprintf (stderr, "debug = %d, min_cov = %d\n", debug, min_cov);

  FILE *input;
  if (optind < argc) input = fopen(argv[optind], "r");
  else input = stdin;

  kstring_t ks = { 0, 0, NULL };
  khash_t(str) *h; // hash for mutations
  h = kh_init(str);
  if (input) {
    for (ks.l = 0; kgetline(&ks, (kgets_func *)fgets, input) == 0; ks.l = 0){
      // printf("new line is %s\n",  ks.s);
      parse_line(&ks, h, debug);
    }
    fclose(input);
  }
	free(ks.s);
  //print output
  khint_t k;
  printf("chrom\tref_start\tref_end\talt\tsize\ttype\tmutCov\n");
  for (k = kh_begin(h); k != kh_end(h); ++k) { // traverse
    	if (!kh_exist(h,k)) continue;
    	const char *kk = kh_key(h, k);
    	int vv = kh_val(h, k);
    	if (vv >= min_cov) printf("%s\t%d\n", kk, vv);
    }

  // free memory
  kh_destroy(str, h);
  return 0;
}