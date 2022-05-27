#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "../klib/kstring.h"
#include "../klib/khash.h"
#include "../klib/kvec.h"
#include "../klib/kseq.h"

// to compile
// gcc -Wall -g -O2 snp_indel_check.c ../klib/kstring.c  -o snp_indel_check -lz


// typedef kvec_t(int) kvecn; // int or char vector
// typedef kvec_t(char *) kvecs; // string vector

// string slice
void slice(const char *str, char *result, size_t start, size_t end)
{
    strncpy(result, str + start, end - start);
    result[end - start] = 0;
}

// void slice(const char * str, char * buffer, size_t start, size_t end)
// {
//     size_t j = 0;
//     for ( size_t i = start; i < end; ++i ) {// not include end
//         buffer[j++] = str[i];
//     }
//     buffer[j] = 0;
// }

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(fasta, char *)      // instantiate structs and methods
KHASH_MAP_INIT_STR(str, int)      // instantiate structs and methods

// fn parse_line (line: &str, map: &mut HashMap<String, isize>, no_small_indels: bool, debug: bool) 
khash_t(fasta) * read_fasta(char *infile){
  gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(infile, "r");
	seq = kseq_init(fp);
  khint_t k;
  int absent;
  khash_t(fasta) *h; // hash for mutations
  h = kh_init(fasta);

	while ((l = kseq_read(seq)) >= 0) {
    // printf("name: %s\n", seq->name.s);
    // if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    // printf("seq: %s\n", seq->seq.s);
    // if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
    k = kh_put(fasta, h, seq->name.s, &absent);
    if (absent) {
      kh_key(h, k) = strdup(seq->name.s); // strdup will malloc, so need to be freed
      kh_val(h, k) = strdup(seq->seq.s);
    }
	}
	// printf("return value: %d\n", l);
	kseq_destroy(seq);
	gzclose(fp);
	return h;
}


// int exclamationCheck = strchr(str, '!') != NULL;
typedef struct {
    kvec_t(char) vop; // M, S etc
    kvec_t(int) vlen; // length
} cigar_info;

cigar_info split_cigar (char *cigar) {
  cigar_info r;// = {vop, vlen };
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

// compare two string to find mutations
// s1 and s2 should have the same length
int putsnps(char *s1, char *s2, khash_t(str) *h, char * chrom, int ref_pos)
{
  size_t dl = strlen(s1); // delim length
  int i;
  for (i=0; i<dl; i++){
    if (s1[i]!=s2[i]) {
      // kv_push(int, *array);
      kstring_t kk = { 0, 0, NULL };
      ksprintf(&kk, "%s\t%d\t%d\t%c\t%c\t0\tsnp", chrom, ref_pos+i+1, ref_pos+i+1, s1[i], s2[i]);
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
  return 0;
}

int putindels(char *ref_seq, char *read_seq, khash_t(str) *h, char *chrom, int ref_pos, int read_pos, int indel_size)
{ 
  kstring_t kk = { 0, 0, NULL };
  if (indel_size > 0){// insertion
    char alt_seq[indel_size+3];
    slice(read_seq, alt_seq, read_pos-1, read_pos+indel_size+1);
    ksprintf(&kk, "%s\t%d\t%d\t%c%c\t%s\t%d\tins", chrom, ref_pos, ref_pos+1, ref_seq[ref_pos-1],ref_seq[ref_pos], alt_seq, indel_size);
  } else {//deletion
    char alt_seq[-indel_size+3];
    slice(ref_seq, alt_seq, ref_pos, ref_pos-indel_size+1);
    ksprintf(&kk, "%s\t%d\t%d\t%s\t%c\t%d\tdel", chrom, ref_pos+1, ref_pos-indel_size+1, alt_seq, ref_seq[ref_pos], indel_size);
  }
  khint_t k;
  int absent;
  k = kh_put(str, h, kk.s, &absent);
  if (!absent) {
    kh_value(h, k) += 1; // set the value
  } else {
    kh_value(h, k) = 1;
  }
  return 0;
}

// fn parse_cigar (cigar: &str, ref_pos: isize, same_strand: bool, read_len: isize)
int * parse_snp(cigar_info *r, char *ref_seq, char *read_seq, khash_t(str) *h, char *chrom, int ref_pos, int debug){
  // cigar_info r = split_cigar(cigar);
  int read_pos1 = 0; // left start
  int read_pos2 = 0; // right end
  int ref_pos1 = ref_pos; // left start
  int ref_pos2 = ref_pos; // right end
  int nmatch = 0; // number of M, if match showed up, then no more S or H
  for (int i =0; i < r->vop.n; i++) {
    int num = r->vlen.a[i];
    char op = r->vop.a[i];
    if (op == 'M' || op == '=' || op == 'X'){
      read_pos2 += num;
      ref_pos2 += num;
      nmatch += 1;
      char ref_slice[num+1];
      char read_slice[num+1];
      slice(ref_seq, ref_slice, ref_pos1, ref_pos2);
      slice(read_seq, read_slice, read_pos1, read_pos2);
      if (debug){
        printf("ref_seq  is %s\n", ref_seq);
        printf("ref_pos1, ref_pos2 are %d and %d\n", ref_pos1, ref_pos2);
        printf("read_pos1, read_pos2 are %d and %d\n", read_pos1, read_pos2);
        printf("ref_slice  is %s\n", ref_slice);
        printf("read_slice is %s\n", read_slice);
      }
      putsnps(ref_slice, read_slice, h, chrom, ref_pos1);
      read_pos1 = read_pos2; // for next match
      ref_pos1 = ref_pos2;
    } else if (op == 'S' || op == 'H') {
      if (nmatch == 0) {
        read_pos1 += num;
        read_pos2 += num;
      }
    } else if (op == 'I'){
      read_pos2 += num;
      putindels(ref_seq, read_seq, h, chrom, ref_pos1, read_pos1, num);
      read_pos1 = read_pos2;
    } else if (op == 'D' || op == 'N') {
      ref_pos2 += num;
      putindels(ref_seq, read_seq, h, chrom, ref_pos1-1, read_pos1-1, -num);
      ref_pos1 = ref_pos2;
    }
  }

  return 0;
}

// fn parse_cigar (cigar: &str, ref_pos: isize, same_strand: bool, read_len: isize)
int * parse_cigar(cigar_info *r, int ref_pos, int same_strand, size_t read_len){
  // struct cigar_info r = split_cigar(cigar);
  int read_pos1 = 0; // left start
  int read_pos2 = -1; // right end
  int ref_pos1 = ref_pos; // left start
  int ref_pos2 = ref_pos - 1; // right end
  int nmatch = 0; // number of M, if match showed up, then no more S or H
  for (int i =0; i < r->vop.n; i++) {
    int num = r->vlen.a[i];
    char op = r->vop.a[i];
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
  // kv_destroy(r.vop);
  // kv_destroy(r.vlen);

  return rr;
}


// fn parse_line (line: &str, map: &mut HashMap<String, isize>, no_small_indels: bool, debug: bool) 
int parse_line(kstring_t *ks, khash_t(str) *h, int debug, khash_t(fasta) *fh){
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
  // printf("cigar is %s\n", cigar);
  // printf("read_seq is %s\n", read_seq);
  // let read_id = ff[0];
  if (strcmp(cigar, "*") == 0 ) { // no mapping
      return 0;
  }
  cigar_info r = split_cigar(cigar);
  // SNPs and small indels
  khint_t k;
  k = kh_get(fasta, fh, chrom);
  char *ref_seq = kh_val(fh, k);
  parse_snp(&r, ref_seq, read_seq, h, chrom, pos, debug);

  // big deletions
  int has_sa = 0;
  for (i = 11; i < n; ++i){
    if (strstr(ks->s + ff[i], "SA:Z") != NULL) {
      has_sa = 1;
      break;
    }
  }
  // SNPs or small indels

  // big deletions or inversions
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
    cigar_info sr = split_cigar(sa_cigar);
    // free: should not free
    // free(s.s); free(ff2);free(ff);
    int all_pos1[4];
    int all_pos2 [4];
    if (strcmp(chrom,sa_chrom)==0 && strcmp(strand, sa_strand) ==0) { // potential big deletion, could be insertion too, but update later
      if (sa_pos > pos) { // SA is on the right
        memcpy(all_pos1, parse_cigar(&r, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos2, parse_cigar(&sr, sa_pos, 1, read_len), 4 * sizeof(int));
      } else {
        memcpy(all_pos2, parse_cigar(&r, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos1, parse_cigar(&sr, sa_pos, 1, read_len), 4 * sizeof(int));
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
        ksprintf(&kk, "%s\t%d\t%d\t\t%s\t%d\tbig_indel", chrom, ref_pos1+1, del_end_pos+1, alt_seq, mut_size);
        // printf("key is %s\n", kk.s);
        // khint_t k;
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
        memcpy(all_pos1, parse_cigar(&r, pos, 1, read_len), 4 * sizeof(int));
        memcpy(all_pos2, parse_cigar(&sr, sa_pos, 0, read_len), 4 * sizeof(int));
      } else {
        memcpy(all_pos2, parse_cigar(&r, pos, 0, read_len), 4 * sizeof(int));
        memcpy(all_pos1, parse_cigar(&sr, sa_pos, 1, read_len), 4 * sizeof(int));
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
        ksprintf(&kk, "%s\t%d\t%d\t\tinversion\t%d\tinv", chrom, ref_pos1+1, del_end_pos+1,  mut_size);
        // khint_t k;
        int absent;
        k = kh_put(str, h, kk.s, &absent);
        if (!absent) {
          kh_value(h, k) += 1; // set the value
        } else {
          kh_value(h, k) = 1;
        }
      }
    }
    // free sr
    kv_destroy(sr.vop);
    kv_destroy(sr.vlen);
  }

  // destroy r
  kv_destroy(r.vop);
  kv_destroy(r.vlen);

  return 0;
}

// main function
int main (int argc, char **argv)
{
  int min_cov = 1;
  int debug = 0;
  // int call_snp = 0; // whether to call snp and small indels
  int c;
  char *fasta_file = NULL;

  opterr = 0;

  while ((c = getopt (argc, argv, "dc:f:")) != -1)
    switch (c)
      {
      case 'd':
        debug = 1;
        break;
      case 'c':
        min_cov = atoi(optarg);
        break;
    case 'f':
        fasta_file = optarg;
        break;
      case '?':
        if (optopt == 'c' || optopt == 'f')
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

  // read fasta
  if (fasta_file == NULL){
    fprintf(stderr, "Please provide a fasta file with template sequences (-f your_sequence.fa)\n");
    fprintf(stderr,
      "Usage: snp_indel_check [options] -f ref.fasta <aln.sam>\n"
      "or:    samtools view aln.bam | snp_indel_check [options] -f ref.fasta\n"
      "Options:\n"
      "  -d debug              print extra information for debugging \n"
      "  -c coverage [int]     minimum coverage for a variant (defualt 1)\n"
      "  -f fasta file name    your reference sequences\n");
    return 1;
  }
  khash_t(fasta) *fh = read_fasta(fasta_file);

  kstring_t ks = { 0, 0, NULL };
  khash_t(str) *h; // hash for mutations
  h = kh_init(str);
  if (input) {
    for (ks.l = 0; kgetline(&ks, (kgets_func *)fgets, input) == 0; ks.l = 0){
      // printf("new line is %s\n",  ks.s);
      parse_line(&ks, h, debug, fh);
    }
    fclose(input);
  }
	free(ks.s);
  //print output
  khint_t k;
  printf("chrom\tref_start\tref_end\tref\talt\tsize\ttype\tmutCov\n");
  for (k = kh_begin(h); k != kh_end(h); ++k) { // traverse
    	if (!kh_exist(h,k)) continue;
    	const char *kk = kh_key(h, k);
    	int vv = kh_val(h, k);
    	if (vv >= min_cov) printf("%s\t%d\n", kk, vv);
      free((char*)kh_key(h, k));
    }

  // free memory
  // for (k = 0; k < kh_end(h); ++k)
  //   if (kh_exist(h, k))
  //     free((char*)kh_key(h, k));
  kh_destroy(str, h);
  
  // free fasta hash
  for (k = 0; k < kh_end(fh); ++k)
    if (kh_exist(fh, k)){
      free((char*)kh_key(fh, k));
      free((char*)kh_val(fh, k));
    }
  kh_destroy(fasta, fh);
  return 0;
}