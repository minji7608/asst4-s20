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



void assign_zones(region_t *region_list, int nregion, int nzone) {
    // TODO.  This partitioner is very naive.  You can do better!
    for (int rid = 0; rid < nregion; rid++) {
	region_list[rid].zone_id = rid % nzone;
    }
}



