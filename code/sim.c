
#include "crun.h"

void print_array(int *array, int length)
{
    outmsg("[ ");
    for (int i = 0; i < length; i++) { 
        outmsg("%d ", array[i]);
    }
    outmsg("]\n");
    return;
}

/* Compute ideal load factor (ILF) for node */
static inline double neighbor_ilf(state_t *s, int nid) {
    graph_t *g = s->g;
    int outdegree = g->neighbor_start[nid+1] - g->neighbor_start[nid] - 1;
    int *start = &g->neighbor[g->neighbor_start[nid]+1];
    int i;
    double sum = 0.0;
    for (i = 0; i < outdegree; i++) {
	int lcount = s->rat_count[nid];
	int rcount = s->rat_count[start[i]];
	double r = imbalance(lcount, rcount);
	sum += r;
    }
    double ilf = BASE_ILF + 0.5 * (sum/outdegree);
    return ilf;
}

/* Compute weight for node nid */
static inline double compute_weight(state_t *s, int nid) {
    int count = s->rat_count[nid];
    double ilf = neighbor_ilf(s, nid);
    return mweight((double) count/s->load_factor, ilf);
}


/* Recompute all node counts according to rat population */
/*
  Function only called at start of simulation, at which point
  every zone has complete rat information.  Can therefore
  have every zone update every node count.
*/
static inline void take_census(state_t *s) {
    graph_t *g = s->g;
    int nnode = g->nnode;
    int *rat_position = s->rat_position;
    int *rat_count = s->rat_count;
    int nrat = s->nrat;

    memset(rat_count, 0, nnode * sizeof(int));
    int ri;
    for (ri = 0; ri < nrat; ri++) {
	rat_count[rat_position[ri]]++;
    }
}

/* Recompute all node weights */
static inline void compute_all_weights(state_t *s) {
    int ni, neighbor, nid, nid2;
    graph_t *g = s->g;
    int nzone = g->nzone;
    int this_zone = g->this_zone;
    START_ACTIVITY(ACTIVITY_WEIGHTS);
    for (int zi = 0; zi < nzone; zi++) {
        if(zi == this_zone) continue;
        int neighbor_count = g->import_node_count[zi];
        for(neighbor = 0; neighbor < neighbor_count; neighbor++){
            nid = g->import_node_list[zi][neighbor];
            s->node_weight[nid] =  compute_weight(s, nid);
        }
    }

    for (ni = 0; ni < g->local_node_count; ni++){
        nid2 = g->local_node_list[ni];   
	    s->node_weight[nid2] = compute_weight(s, nid2);
    }

    FINISH_ACTIVITY(ACTIVITY_WEIGHTS);

}

// /* Recompute all node weights */ 
// static inline void compute_all_weights(state_t *s) {
//     int nid;
//     graph_t *g = s->g;
//     double *node_weight = s->node_weight;
//     START_ACTIVITY(ACTIVITY_WEIGHTS);
//     for (nid = 0; nid < g->nnode; nid++)
// 	node_weight[nid] = compute_weight(s, nid);
//     FINISH_ACTIVITY(ACTIVITY_WEIGHTS);

// }



/* In synchronous or batch mode, can precompute sums for each region */
static inline void find_all_sums(state_t *s) {
    graph_t *g = s->g;
    
    START_ACTIVITY(ACTIVITY_SUMS);
    // outmsg("------------- executing find_all_sums -------------\n");
    // TODO: It doesn't make sense to compute the weights for nodes that are not
    // in the local zone
    int local_node_count = g->local_node_count;
    int ni, nid, eid;
    for (ni = 0; ni < local_node_count; ni++) {
        nid = g->local_node_list[ni];
	    double sum = 0.0;
	    for (eid = g->neighbor_start[nid]; eid < g->neighbor_start[nid+1]; eid++) {
            // outmsg("node [%d] weight: %f\n", g->neighbor[eid], s->node_weight[g->neighbor[eid]]);
	        sum += s->node_weight[g->neighbor[eid]];
	        s->neighbor_accum_weight[eid] = sum;
	    }
	    s->sum_weight[nid] = sum;
        // outmsg("nid: %d --- tsum: %f\n", nid, sum); 
        
    }
    // outmsg("-----------\n");
    FINISH_ACTIVITY(ACTIVITY_SUMS);
}

/*
  Given list of increasing numbers, and target number,
  find index of first one where target is less than list value
*/

/*
  Linear search
 */
static inline int locate_value_linear(double target, double *list, int len) {
    int i;
    for (i = 0; i < len; i++)
	if (target < list[i])
	    return i;
    /* Shouldn't get here */
    return -1;
}
/*
  Binary search down to threshold, and then linear
 */
static inline int locate_value(double target, double *list, int len) {
    int left = 0;
    int right = len-1;
    while (left < right) {
	if (right-left+1 < BINARY_THRESHOLD)
	    return left + locate_value_linear(target, list+left, right-left+1);
	int mid = left + (right-left)/2;
	if (target < list[mid])
	    right = mid;
	else
	    left = mid+1;
    }
    return right;
}


/*
  Version that can be used in synchronous or batch mode, where certain that node weights are already valid.
  And have already computed sum of weights for each node, and cumulative weight for each neighbor
  Given list of integer counts, generate real-valued weights
  and use these to flip random coin returning value between 0 and len-1
*/
static inline int fast_next_random_move(state_t *s, int r) {
    int nid = s->rat_position[r];
    graph_t *g = s->g;
    random_t *seedp = &s->rat_seed[r];
    /* Guaranteed that have computed sum of weights */
    double tsum = s->sum_weight[nid];   

    
    double val = next_random_float(seedp, tsum);
    // outmsg("rat: %d seedp: %u tsum: %f ---> val: %f\n", r, *seedp, tsum, val); 

    int estart = g->neighbor_start[nid];
    int elen = g->neighbor_start[nid+1] - estart;
    int offset = locate_value(val, &s->neighbor_accum_weight[estart], elen);
#if DEBUG
    if (offset < 0) {
	/* Shouldn't get here */
	outmsg("Internal error.  fast_next_random_move.  Didn't find valid move.  Target = %.2f/%.2f.\n",
	       val, tsum);
	return 0;
    }
#endif
#if DEBUG
    outmsg("Computing rat %d: node %d-->%d (%.3f/%.3f)", r, nid, g->neighbor[estart+offset], val, tsum);
#endif
    return g->neighbor[estart + offset];
}

/* Process single batch */
// TODO: Here's where things get interesting!
//    * Process rats currently in this zone
//    * Exchange rats (use function exchange_rats):
//      - Export rats that move out of this zone
//      - Import rats that move into this zone
//    * Exchange node states (use function exchange_node_states)
//      - Export counts for internal nodes adjacent to other zones
//      - Import counts for external nodes adjacent to this zone
//    * Compute weights for nodes in this zone
//    * Exchange weights (use function exchange_node_weights)
//      - Export weights for internal nodes adjacent to other zones
//      - Import weights for external nodes adjacent to this zone
static inline void do_batch(state_t *s, int batch, int bstart, int bcount) {
    int rid, ri, zi, numrats;
    find_all_sums(s);
    int *zone_id = s->g->zone_id;
    int this_zone = s->g->this_zone;
    int nzone = s->g->nzone; 

    // initialize export counts to 0 for all other zones
    for (zi = 0; zi < nzone; zi++) {
        s->export_numrats[zi] = 0;
    }

    for (ri = 0; ri < bcount; ri++) {
        rid = ri + bstart;
        char in_this_zone = s->zone_rat_bitvector[rid];
        
        // outmsg("zone: %d rid: %d in_this_zone: %d", this_zone, rid, in_this_zone);

        // only process rats currently in the zone
        if ((int)in_this_zone == 1) {
            
            int onid = s->rat_position[rid];

            // outmsg("processing rat %d in node %d\n", rid, onid);
            int nnid = fast_next_random_move(s, rid);
            int new_zone = zone_id[nnid];

           

            // if moving within the zone
            if (new_zone == this_zone) {
                s->rat_position[rid] = nnid;
                s->rat_count[onid] -= 1;
                s->rat_count[nnid] += 1;
            }

            // if moving to a new zone            
            else {

                s->rat_count[onid] -= 1;
                // clear this zone's bitvector
                s->zone_rat_bitvector[rid] = 0;
                
                numrats = s->export_numrats[new_zone];

                s->export_rat_info[new_zone][numrats * 3] = rid;
                s->export_rat_info[new_zone][numrats * 3 + 1] = nnid;
                s->export_rat_info[new_zone][numrats * 3 + 2] = (int)(s->rat_seed[rid]);

                s->export_numrats[new_zone]++;
            }


        }
    }

    // for (int i=0; i<nzone; i++) {
        
    //     // if (s->export_numrats[i] != 0){
    //     //     outmsg("EXPORT_RID (%d ---> %d)\n", this_zone, i);
    //     //     print_array(s->export_rid[i], s->export_numrats[i]);
    //     //     print_array(s->export_nid[i], s->export_numrats[i]);
    //     // }
        
    // }
    
    
    /* Update weights */
#if MPI
    exchange_rats(s);
    exchange_node_states(s);
    compute_all_weights(s);
    exchange_node_weights(s);
#endif
}

static void batch_step(state_t *s) {
    int bstart = 0;
    int bsize = s->batch_size;
    int nrat = s->nrat;
    int bcount;
    int batch = 0;
    while (bstart < nrat) {
	bcount = nrat - bstart;
	if (bcount > bsize)
	    bcount = bsize;
	do_batch(s, batch, bstart, bcount);
	batch++;
	bstart += bcount;
    }
}

double simulate(state_t *s, int count, int dinterval, bool display) {
    int i;
    /* Compute and show initial state */
    bool show_counts = true;
    double start = currentSeconds();
    take_census(s);
    compute_all_weights(s);
    if (display) {
#if MPI
	if (s->g->this_zone == 0)
	    // Process 0 has a copy of the initial counts for all nodes.
	    show(s, show_counts);
#else
	show(s, show_counts);
#endif
    }
    for (i = 0; i < count; i++) {
	batch_step(s);
	if (display) {
	    show_counts = (((i+1) % dinterval) == 0) || (i == count-1);
#if MPI
	    if (s->g->this_zone == 0) {
		// Process 0 needs to call function show on each simulation step.
		// When show_counts is true, it will need to have
		// the counts for all other zones.
		// These must be gathered from the other processes.
		if (show_counts)
		    gather_node_state(s);
		show(s, show_counts);
	    } else {
		if (show_counts)
		    // Send counts to process 0
		    send_node_state(s);
	    }
#else
	    show(s, show_counts);
#endif
	}
    }
    double delta = currentSeconds() - start;
    done(s);
    return delta;
}

