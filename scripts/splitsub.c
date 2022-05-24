#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../klib/kvec.h"

// please just use splitsubv3 as an example
//v1 use kvec, but always has a warning
// split by a substring, not by char
// to fix the warning, seems because it is a macro of "struct <anonymous> *"
// I have to define my kvec to a type, the kvec_t funtions can still be used since kvec has the same structure

// struct { size_t n, m; type *a; }
// typedef struct {
//   size_t n, m; // n is number of elements, m is the capacity
//   char * *a;
// } kvec;

typedef kvec_t(char *) kvec;

// int splitsub(char *str, const char *delim, kvec_t(char *) *array)
int splitsub(char *str, const char *delim, kvec *array)
{
    size_t dl = strlen(delim); // delim length
    kv_push(char *, *array, str);
    char *p, *tmp;
    p = tmp = str;
    while (p != NULL){
        p = strstr(tmp, delim);
        if (p == NULL) return 0;
        *p = '\0';
        tmp = p + dl;
        printf("tmp is %s \n", tmp);
        kv_push(char *, *array, tmp);
    }
    return 0;
}

int main(void)
{
    char str[] = "ABcdxxefgxxhhhhhxx";
    // kvec_t(char *) sublist;
    kvec sublist;
    kv_init(sublist);
    splitsub(str, "xx", &sublist);
    for (int i = 0; i<sublist.n; i++) printf("%d substring is: %s\n", i+1, sublist.a[i]);
    kv_destroy(sublist);
    return 0;
}