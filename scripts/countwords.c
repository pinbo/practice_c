#include <stdio.h>
#include <string.h>
//#include "../htslib/khash.h"
#include "../klib/khash.h"
//KHASH_SET_INIT_STR(str)
KHASH_MAP_INIT_STR(str, int)      // instantiate structs and methods

int main(int argc, char *argv[])
{
    khash_t(str) *h;
    khint_t k;
    int i, absent;
    h = kh_init(str);
    for (i = 1; i < argc; ++i){
        k = kh_put(str, h, argv[i], &absent);
        if (!absent) {
        	kh_value(h, k) += 1; // set the value
        } else {
        	kh_value(h, k) = 1;
    	}
    }
    printf("# of distinct words: %d\n", kh_size(h));
    for (k = kh_begin(h); k != kh_end(h); ++k) { // traverse
    	if (!kh_exist(h,k)) continue;
    	char *kk = kh_key(h, k);
    	int vv = kh_val(h, k);
    	printf("%s has count %d\n", kk, vv);
    }

    kh_destroy(str, h);
    return 0;
}
