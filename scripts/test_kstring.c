#include <stdio.h>
#include <ctype.h>
#include "../klib/kstring.h"

int kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list args;
	int l;
	va_copy(args, ap);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(args);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		va_copy(args, ap);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args);
		va_end(args);
	}
	s->l += l;
	return l;
}

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = kvsprintf(s, fmt, ap);
	va_end(ap);
	return l;
}

// s MUST BE a null terminated string; l = strlen(s)
int ksplit_core(char *s, int delimiter, int *_max, int **_offsets)
{
	int i, n, max, last_char, last_start, *offsets, l;
	n = 0; max = *_max; offsets = *_offsets;
	l = strlen(s);
	
#define __ksplit_aux do {						\
		if (_offsets) {						\
			s[i] = 0;					\
			if (n == max) {					\
				int *tmp;				\
				max = max? max<<1 : 2;			\
				if ((tmp = (int*)realloc(offsets, sizeof(int) * max))) {  \
					offsets = tmp;			\
				} else	{				\
					free(offsets);			\
					*_offsets = NULL;		\
					return 0;			\
				}					\
			}						\
			offsets[n++] = last_start;			\
		} else ++n;						\
	} while (0)

	for (i = 0, last_char = last_start = 0; i <= l; ++i) {
		if (delimiter == 0) {
			if (isspace(s[i]) || s[i] == 0) {
				if (isgraph(last_char)) __ksplit_aux; // the end of a field
			} else {
				if (isspace(last_char) || last_char == 0) last_start = i;
			}
		} else {
			if (s[i] == delimiter || s[i] == 0) {
				if (last_char != 0 && last_char != delimiter) __ksplit_aux; // the end of a field
			} else {
				if (last_char == delimiter || last_char == 0) last_start = i;
			}
		}
		last_char = s[i];
	}
	*_max = max; *_offsets = offsets;
	return n;
}

int kgetline(kstring_t *s, kgets_func *fgets_fn, void *fp)
{
	size_t l0 = s->l;

	while (s->l == l0 || s->s[s->l-1] != '\n') {
		if (s->m - s->l < 200) ks_resize(s, s->m + 200);
		if (fgets_fn(s->s + s->l, s->m - s->l, fp) == NULL) break;
		s->l += strlen(s->s + s->l);
	}

	if (s->l == l0) return EOF;

	if (s->l > l0 && s->s[s->l-1] == '\n') {
		s->l--;
		if (s->l > l0 && s->s[s->l-1] == '\r') s->l--;
	}
	s->s[s->l] = '\0';
	return 0;
}

int main(int argc, char **argv)
{
	// kstring_t *s;
	// int *fields, n, i;
	// s = (kstring_t*)calloc(1, sizeof(kstring_t));
	// // test ksprintf()
    kstring_t s = { 0, 0, NULL };
	ksprintf(&s, " abc\tdefg:    %d ", 100);
	printf("'%s'\n", s.s);
    char *p = "opq";
    kputs(p, &s);
    printf("'%s'\n", s.s);

	// // test ksplit()
	// // fields = ksplit(s, 0, &n);
    // fields = ksplit(s, '\t', &n);
	// for (i = 0; i < n; ++i)
	// 	printf("field[%d] = '%s'\n", i, s->s + fields[i]);
    // // free
	// free(s->s); free(s); free(fields);
    //
    kstring_t ks = { 0, 0, NULL };
    int *fields, n;
    if (argc > 1) {
        FILE *f = fopen(argv[1], "r");
        if (f) {
            for (ks.l = 0; kgetline(&ks, (kgets_func *)fgets, f) == 0; ks.l = 0){
                // printf("new line is %s\n",  ks.s);
                fields = ksplit(&ks, '\t', &n);
                char *chrom = ks.s + fields[2];
                printf("field[2] = '%s'\n",  chrom);
            }
            fclose(f);
        }
    }

	free(ks.s);

	return 0;
}
