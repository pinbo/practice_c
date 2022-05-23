#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "../klib/kstring.h"
#include "../klib/khash.h"
#include "../klib/kvec.h"

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
    //   tmp[0] = 0; // wrong
      memset(tmp, 0, sizeof(tmp));
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


int main (int argc, char **argv)
{
//   char *cigar = "8M1D108M";
  if (argc > 1){
    char *cigar = argv[1];
    printf("CIGAR is %s\n", cigar);
    struct cigar_info r = split_cigar(cigar);
    for (int i =0; i < r.vop.n; i++){
        printf("%c has length %d\n", r.vop.a[i], r.vlen.a[i]);
    }
    // get positions
    int *pp;
    printf("CIGAR is %s\n", cigar);
    pp = parse_cigar(cigar, 1, 1, 150);
    for (int i = 0; i < 4; i++ ) {
        // printf( "*(p + %d) : %d\n", i, *(pp + i));
        printf( "*(p + %d) : %d\n", i, pp[i]);
    }
  }
  return 0;
}