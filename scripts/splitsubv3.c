#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../klib/kstring.h"
// v3: add splitsub to kstring.h and kstring.c
// gcc -Wall -g splitsubv3.c ../klib/kstring.c -o splitsubv3

int main(void)
{
    char str[] = "ABcdxxefgxxhxxiixx";
    printf("raw string is %s\n", str);
    kvec sublist;
    initArray(&sublist,2); // inilize a vector of length 2, you can give any numbers
    printf("sublist has capacity %ld\n", sublist.m);
    splitsub(str, "xx", &sublist);
    printf("finished spliting\n");
    for (int i = 0; i < sublist.n; i++) printf("%d substring is: %s\n", i+1, sublist.a[i]);
    freeArray(&sublist);

    return 0;
}