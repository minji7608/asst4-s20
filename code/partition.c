#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "crun.h"


/* This file is where you will write your partitioner */

/*
  This function should assign a zone id to every region in the graph.
  Zone IDs should range between 0 and Z-1.
  The updates should be made to the zone_id field of
  each entry in the region array

  Note that each region entry contains information about the position,
  node count, and edge count of the region.

  You should feel free to try out different partitioning schemes.
*/

int compare(const void *s1, const void *s2) {
  region_t *r1 = (region_t *)s1;
  region_t *r2 = (region_t *)s2;

  int r1_edge_count = r1->edge_count;
  int r2_edge_count = r2->edge_count;

  if (r1_edge_count < r2_edge_count) return -1;
  else if (r1_edge_count > r2_edge_count) return 1;
  else return 0;
}

void assign_zones(region_t *region_list, int nregion, int nzone) {
    // TODO:  This partitioner is very naive.  You can do better!
    // for (int rid = 0; rid < nregion; rid++) {
	  //   region_list[rid].zone_id = rid % nzone;
    // }
    qsort(region_list, nregion, sizeof(region_t), compare);


    double weights[nregion];
    int zones[nzone];
    int rid, zid;

    for (rid=0; rid<nregion; rid++) {
      weights[rid] = (double)(region_list[rid].edge_count);
    }

    find_partition(nregion, nzone, weights, zones);
    //fprintf(stderr, "%d\n", zones[0]);
    //fprintf(stderr, "%d\n", zones[1]);
    //fprintf(stderr, "%d\n", zones[2]);
    //fprintf(stderr, "%d\n", zones[3]);
    
    int curr_rid = 0;
    for (zid=0; zid<nzone; zid++) {
      int end_rid = curr_rid + zones[zid];
      // fpri ntf(stderr, "%d\n", end_rid);
      while (curr_rid < end_rid) {
        region_list[curr_rid].zone_id = zid;
        curr_rid++;
      }
    }
}



