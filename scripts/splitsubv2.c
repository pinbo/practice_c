#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// v2: self defined struct
typedef struct {
  char **array;
  size_t n; // array length
  size_t size; // capacity
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(char*));
  a->n = 0;
  a->size = initialSize;
}

void insertArray(Array *a, char* element) {
  // a->n is the number of used entries, because a->array[a->n++] updates a->n only *after* the array has been accessed.
  // Therefore a->n can go up to a->size 
  if (a->n == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(char*));
  }
  a->array[a->n++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->n = a->size = 0;
}

// split by a substring, not by char
int splitsub(char *str, const char *delim, Array *array)
{
    size_t dl = strlen(delim); // delim length
    insertArray(array, str);
    char *tmp;
    char *p;
    p = tmp = str;
    while (p != NULL){
        p = strstr(tmp, delim);
        if (p == NULL) return 0;
        *p = '\0';
        tmp = p + dl;
        printf("tmp is %s \n", tmp);
        insertArray(array, tmp);
    }
    return 0;
}

int main(void)
{
    char str[] = "ABcdxxefgxxhxxiixx";
    printf("raw string is %s\n", str);
    Array sublist;
    initArray(&sublist,2); // inilize a vector of length 2, you can give any numbers
    printf("sublist has capacity %ld\n", sublist.size);
    splitsub(str, "xx", &sublist);
    printf("finished spliting\n");
    for (int i = 0; i < sublist.n; i++) printf("%d substring is: %s\n", i+1, sublist.array[i]);
    freeArray(&sublist);

    return 0;
}