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

compare(const void *s1, const void *s2) {
  region_t *r1 = (region_t *)s1;
  region_t *r2 = (region_t *)s2;

  int r1_edge_count = r1->edge_count;
  int r2_edge_count = r2->edge_count;

  if (r1_edge_count < r2_edge_count) return -1;
  else if (r1_edge_count > r2_edge_count) return 1;
  else return 0;
}


float stdev(region_t data[], int len, bool by_node) {
    int sum = 0;
    float mean;
    float SD = 0.0;
    int i;

    // stdev of node counts
    if (by_node) {
      for (i = 0; i < len; ++i) {
        sum += data[i].node_count;
      }
    
      mean = (float)(sum / len);
      
      for (i = 0; i < len; ++i)
          SD += pow(data[i].node_count - mean, 2);
          
      return sqrt(SD / len);
    }
    
    // stdev of edge counts
    else {
      for (i = 0; i < len; ++i) {
        sum += data[i].edge_count;
      }
    
      mean = (float)(sum / len);
      
      for (i = 0; i < len; ++i)
          SD += pow(data[i].edge_count - mean, 2);
          
      return sqrt(SD / len);
    }

}

void assign_zones(region_t *region_list, int nregion, int nzone) {
    // TODO:  This partitioner is very naive.  You can do better!

    qsort(region_list, nregion, sizeof(region_t), compare);


    double weights[nregion];
    int zones[nzone];
    int rid, zid;
  
    float node_stdev = stdev(region_list, nregion, true);
    float edge_stdev = stdev(region_list, nregion, false);
    
    if (node_stdev > edge_stdev) {
      for (rid=0; rid<nregion; rid++) {
        weights[rid] = (double)(region_list[rid].node_count);
      }
    }

    else {
      for (rid=0; rid<nregion; rid++) {
        weights[rid] = (double)(region_list[rid].edge_count);
      }
    }

    find_partition(nregion, nzone, weights, zones);
    
    int curr_rid = 0;
    for (zid=0; zid<nzone; zid++) {
      int end_rid = curr_rid + zones[zid];
      while (curr_rid < end_rid) {
        region_list[curr_rid].zone_id = zid;
        curr_rid++;
      }
    }
}


// int compare(const void *s1, const void *s2) {
//   region_t *r1 = (region_t *)s1;
//   region_t *r2 = (region_t *)s2;

//   int r1_edge_count = r1->edge_count;
//   int r2_edge_count = r2->edge_count;

//   if (r1_edge_count < r2_edge_count) return -1;
//   else if (r1_edge_count > r2_edge_count) return 1;
//   else return 0;
// }

// float stdev(region_t data[], int len) {
//   int sum = 0;
//   float mean;
//   float SD = 0.0;
//   int i;

//   // stdev of node counts
//   for (i = 0; i < len; ++i) {
//     sum += data[i].node_count;
//     mean = (float)(sum / len);
    
//     for (i = 0; i < len; ++i)
//       SD += pow(data[i].node_count - mean, 2);
//   } 
    
//   return sqrt(SD / len);
// }

// void assign_zones(region_t *region_list, int nregion, int nzone) {
//     // TODO:  This partitioner is very naive.  You can do better!

//     // qsort(region_list, nregion, sizeof(region_t), compare);


//     double weights[nregion];
//     double tempweights[nzone];
//     int zones[nzone];
//     int rid, zid;

//     for (rid=0; rid<nregion; rid++) {
//       int firstval = region_list[rid].edge_count;
//       int secondval = region_list[rid].node_count;
//       int finalval = firstval + secondval;
//       weights[rid] = (double)(finalval);
//     }


//     // //interval partioning
//     // int intervalnum = nregion/nzone;
//     // int rid = 0;
//     // tempweights[rid] = region_list[rid].node_count;
//     // for(zid = 0; zid < nzone; zid++){
//     //     if(rid % intervalnum != 0){
//     //       tempweights[zid] = region_list[rid].node_count;
//     //       rid++;
//     //     }
//     // }

//     find_partition(nregion, nzone, weights, zones);


//     int interval_stdev = stdev(zones, nzone);
//     int greedy_stdev = stdev(tempweights, nzone);


//     if(greedy_stdev < interval_stdev){

      
    
//     int curr_rid = 0;
//     for (zid=0; zid<nzone; zid++) {
//       int end_rid = curr_rid + zones[zid];
//       // fpri ntf(stderr, "%d\n", end_rid);
//       while (curr_rid < end_rid) {
//         region_list[curr_rid].zone_id = zid;
//         curr_rid++;
//       }
//     }
// }
