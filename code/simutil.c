#include "crun.h"
#include "assert.h"

void outmsg(char *fmt, ...) {
#if MPI
    int process_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
    if (process_id != 0)
	fprintf(stderr, "Process %.2d|", process_id);
#endif    
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    bool got_newline = fmt[strlen(fmt)-1] == '\n';
    if (!got_newline)
	fprintf(stderr, "\n");
}

/* Allocate n int's and zero them out.  Maybe you could use multiple threads ... */
int *int_alloc(size_t n) {
    return (int *) calloc(n, sizeof(int));
}

/* Allocate n doubles's and zero them out.  Maybe you could use multiple threads ... */
double *double_alloc(size_t n) {
    return (double *) calloc(n, sizeof(double));
}

/* Allocate n random number seeds and zero them out.  */
static random_t *rt_alloc(size_t n) {
    return (random_t *) calloc(n, sizeof(random_t));
}

/* Allocate simulation state */
static state_t *new_rats(graph_t *g, int nrat, random_t global_seed) {
    int nnode = g->nnode;
    
    state_t *s = malloc(sizeof(state_t));
    if (s == NULL) {
	outmsg("Couldn't allocate storage for state\n");
	return NULL;
    }

    s->g = g;
    s->nrat = nrat;
    s->global_seed = global_seed;
    s->load_factor = (double) nrat / nnode;

    /* Compute batch size as max(BATCH_FRACTION * R, sqrt(R)) */
    int rpct = (int) (BATCH_FRACTION * nrat);
    int sroot = (int) sqrt(nrat);
    if (rpct > sroot)
	s->batch_size = rpct;
    else
	s->batch_size = sroot;

    // Allocate data structures
    bool ok = true;
    s->rat_position = int_alloc(nrat);
    ok = ok && s->rat_position != NULL;
    s->rat_seed = rt_alloc(nrat);
    ok = ok && s->rat_seed != NULL;
    s->rat_count = int_alloc(nnode);
    ok = ok && s->rat_count != NULL;

    s->node_weight = double_alloc(nnode);
    ok = ok && s->node_weight != NULL;
    s->sum_weight = double_alloc(g->nnode);
    ok = ok && s->sum_weight != NULL;
    s->neighbor_accum_weight = double_alloc(g->nnode + g->nedge);
    ok = ok && s->neighbor_accum_weight != NULL;

    if (!ok) {
	    outmsg("Couldn't allocate space for %d rats", nrat);
	    return NULL;
    }
    return s;
}

/* Set seed values for the rats.  Maybe you could use multiple threads ... */
static void seed_rats(state_t *s) {
    random_t global_seed = s->global_seed;
    int nrat = s->nrat;
    int r;
    for (r = 0; r < nrat; r++) {
	random_t seeds[2];
	seeds[0] = global_seed;
	seeds[1] = r;
	reseed(&s->rat_seed[r], seeds, 2);
#if DEBUG
	if (r == TAG)
	    outmsg("Rat %d.  Setting seed to %u\n", r, (unsigned) s->rat_seed[r]);
#endif
    }
}

/* See whether line of text is a comment */
static inline bool is_comment(char *s) {
    int i;
    int n = strlen(s);
    for (i = 0; i < n; i++) {
	char c = s[i];
	if (!isspace(c))
	    return c == '#';
    }
    return false;
}

/* Read in rat file */
state_t *read_rats(graph_t *g, FILE *infile, random_t global_seed) {
    char linebuf[MAXLINE];
    int r, nnode, nid, nrat;

    // Read header information
    while (fgets(linebuf, MAXLINE, infile) != NULL) {
	if (!is_comment(linebuf))
	    break;
    }
    if (sscanf(linebuf, "%d %d", &nnode, &nrat) != 2) {
	outmsg("ERROR. Malformed rat file header (line 1)\n");
	return false;
    }
    if (nnode != g->nnode) {
	outmsg("Graph contains %d nodes, but rat file has %d\n", g->nnode, nnode);
	return NULL;
    }
    
    state_t *s = new_rats(g, nrat, global_seed);


    for (r = 0; r < nrat; r++) {
        while (fgets(linebuf, MAXLINE, infile) != NULL) {
            if (!is_comment(linebuf))
            break;
        }
        if (sscanf(linebuf, "%d", &nid) != 1) {
            outmsg("Error in rat file.  Line %d\n", r+2);
            return false;
        }
        if (nid < 0 || nid >= nnode) {
            outmsg("ERROR.  Line %d.  Invalid node number %d\n", r+2, nid);
            return false;
        }
        s->rat_position[r] = nid;
    }
    fclose(infile);

    seed_rats(s);
    outmsg("Loaded %d rats\n", nrat);
#if DEBUG
    outmsg("Load factor = %f\n", s->load_factor);
#endif
    return s;
}

/* print state of nodes */
void show(state_t *s, bool show_counts) {
    int nid;
    graph_t *g = s->g;
    printf("STEP %d %d %d\n", g->width, g->height, s->nrat);
    if (show_counts) {
	    for (nid = 0; nid < g->nnode; nid++)
		printf("%d\n", s->rat_count[nid]);
    }
    printf("END\n");
}

/* Print final output */
void done(state_t *s) {
#if MPI
    if (s == NULL || s->g->this_zone != 0)
	return;
#endif
    printf("DONE\n");
}

//TODO: Write function to initialize zone
bool init_zone(state_t *s, int zid) {

    int nzone = s->g->nzone;
    int nrat = s->nrat;
    int nnode = s->g->nnode;
    bool ok = true;

    s->import_rat_info = calloc(nzone, sizeof(int*));
    ok = ok && s->import_rat_info != NULL;
    s->export_rat_info = calloc(nzone, sizeof(int*));
    ok = ok && s->export_rat_info != NULL;

    s->import_nid = calloc(nzone, sizeof(int*));
    ok = ok && s->import_nid != NULL;
    s->export_nid = calloc(nzone, sizeof(int*));
    ok = ok && s->export_nid != NULL;

    //s->import_numrats = int_alloc(nzone);
    //ok = ok && s->import_numrats != NULL;
    s->export_numrats = int_alloc(nzone);
    ok = ok && s->export_numrats != NULL;

    s->import_rat_count = calloc(nzone, sizeof(int*));
    ok = ok && s->import_rat_count != NULL;
    s->export_rat_count = calloc(nzone, sizeof(int*));
    ok = ok && s->export_rat_count != NULL;

    s->import_node_state = calloc(nzone, sizeof(int*));
    ok = ok && s->import_node_state != NULL;
    s->export_node_state = calloc(nzone, sizeof(int*));
    ok = ok && s->export_node_state != NULL;

    s->import_node_weight = calloc(nzone, sizeof(double *));
    ok = ok && s->import_node_weight != NULL;
    s->export_node_weight = calloc(nzone, sizeof(double *));
    ok = ok && s->export_node_weight != NULL;
    
    s->zone_rat_list = int_alloc(nrat);
    ok = ok && s->zone_rat_list != NULL;

    s->zone_rat_bitvector = calloc(nrat, sizeof(unsigned char));
    ok = ok && s->zone_rat_bitvector != NULL;

    int i;
    int num = s->batch_size;
    // int local_nnode = s->g->local_node_count;
    for (i=0; i<nzone; i++) {
        s->import_rat_info[i] = int_alloc(num * 3);
        s->export_rat_info[i] = int_alloc(num * 3);

        s->import_nid[i] = int_alloc(nnode);
        s->export_nid[i] = int_alloc(nnode);

        s->import_node_state[i] = int_alloc(nnode);
        s->export_node_state[i] = int_alloc(nnode);

        s->export_rat_count[i] = int_alloc(nnode);
        s->import_rat_count[i] = int_alloc(nnode);
        s->import_node_weight[i] = double_alloc(nnode);
        s->export_node_weight[i] = double_alloc(nnode);

        ok = ok && 
             (s->import_rat_info[i] != NULL) && 
             (s->export_rat_info[i] != NULL) &&
             (s->import_rat_count[i] != NULL) &&
             (s->export_rat_count[i] != NULL) &&
             (s->import_node_state[i] != NULL) &&
             (s->export_node_state[i] != NULL) &&
             (s->import_node_weight[i] != NULL) &&
             (s->export_node_weight[i] != NULL);
    }

    if (!ok) return false;

    int ri;
    int count = 0;

    for (ri=0; ri<nrat; ri++) {
        int ni = s->rat_position[ri];
        if (s->g->zone_id[ni] == zid) {
            s->zone_rat_list[count] = ri;
            s->zone_rat_bitvector[ri] = 1;
            count++;
        }
    }
    
    return true;
}

//TODO: Implement these communication-support functions
#if MPI
/* Called by process 0 to distribute rat state to all nodes */
void send_rats(state_t *s) {
    int nrat = s->nrat;

    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);

    MPI_Bcast(&(nrat), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(s->rat_position, nrat, MPI_INT, 0, MPI_COMM_WORLD);

    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Called by other nodes to get rat state from master and set up state data structure */
state_t *get_rats(graph_t *g, random_t global_seed) {
    
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    int nrat = 0;
    /*
      Your code should go here.
      It should receive information about all rats, including how many there are,
      from process 0.
     */
    MPI_Bcast(&(nrat), 1, MPI_INT, 0 , MPI_COMM_WORLD);

    state_t *s = new_rats(g, nrat, global_seed);

    // int rat_position[nrat];
    MPI_Bcast(s->rat_position, nrat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    seed_rats(s);

    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);

    
    // memcpy(s->rat_position, rat_position, nrat * sizeof(int));

    //s->rat_position = rat_position;
     // reseeding the rats
                                                  
    return s;
}

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(state_t *s) {
    int zi, i;
    int nzone = s->g->nzone;
    int nnode = s->g->nnode;
    int len = nzone-1;
    MPI_Request request[nzone];
    int import_numrats[nzone];
    
    graph_t *g = s->g;

    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);

    /* Your code should go here */



    for (zi = 1; zi < nzone; zi++){ //goes through all the processes
        MPI_Probe(zi, zi, MPI_COMM_WORLD, &(request[zi]));
        MPI_Get_count(&(request[zi]), MPI_INT, &(import_numrats[zi]));

        MPI_Recv(s->import_nid[zi], import_numrats[zi], MPI_INT, zi, zi, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        MPI_Recv(s->import_rat_count[zi], import_numrats[zi], MPI_INT, zi, zi, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  

        for (i = 0; i < import_numrats[zi]; i++) {
            s->rat_count[s->import_nid[zi][i]] = s->import_rat_count[zi][i];
        }      
    }   

    // MPI_Waitall(len, request, MPI_STATUS_IGNORE);
    
    // for(i = 0; i < nnode; i++){
    //     int* zone_id_list = g->zone_id;
    //     int zid = zone_id_list[i]; // zone that the node is in
    //     if (zid != 0) {
    //         s->rat_count[i] = (s->import_rat_count)[zid][i];
    //     }
    // }

    // MPI_Status status[nzone];
    // int *import_rat

    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);

}

/* Called by other processes to send their node states to process 0 */
void send_node_state(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    int i, ni, zi;
    graph_t *g = s->g;
    int nnode = g->nnode;
    zi = g->this_zone;

    int count = 0;
    for (i=0; i<g->local_node_count; i++) {
        ni = g->local_node_list[i];
        s->export_nid[zi][i] = ni;
        s->export_rat_count[zi][count] = s->rat_count[ni];
        count++;
    }
    
    MPI_Request request[(g->nzone) * 2];

    MPI_Isend(s->export_nid[zi], count, MPI_INT, 0, zi, MPI_COMM_WORLD, &(request[zi*2]));
    MPI_Isend(s->export_rat_count[zi], count, MPI_INT, 0, zi, MPI_COMM_WORLD, &(request[zi*2+1]));
    
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Move rats between zones as they migrate */
void exchange_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    
    int zi, ri;
    int nzone = s->g->nzone;
    int this_zone = s->g->this_zone;
    int len = nzone;
    int export_numrats;
    MPI_Request request[len];
    MPI_Status status[len];
    int probe_ncount[len];

    // send data to all the other zones (async)
    for (zi = 0; zi < nzone; zi++) {

        if (zi == this_zone) continue;

        export_numrats = s->export_numrats[zi];
        // MPI_Isend(&(s->export_numrats[zi]), 1, MPI_INT, zi, zi*2, MPI_COMM_WORLD, &(request[zi*2]));
        // if (export_numrats != 0) {                    
        MPI_Isend(s->export_rat_info[zi], export_numrats * 3, MPI_INT, zi, zi, MPI_COMM_WORLD, &(request[zi]));
        // }
    }
    
    // receive data from all the other zones (synch)
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;

        // /MPI_Recv(&(s->import_numrats[zi]), 1, MPI_INT, zi, this_zone*2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //int probe_ncount;

        // probe size of the message
        MPI_Probe(zi, this_zone, MPI_COMM_WORLD, &(status[zi]));
        MPI_Get_count(&(status[zi]), MPI_INT, &(probe_ncount[zi]));


        // assert(s->import_numrats[zi] * 3 == probe_ncount);
        // outmsg("probed: %d\n", probe_ncount[zi]);
        // int ncount = probe_ncount[zi];

        // if (import_numrats != 0) {
        MPI_Recv(s->import_rat_info[zi], probe_ncount[zi], MPI_INT, zi, this_zone, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }
    }

    // move received data
    // int old_rat_count = s->zone_rat_count;
    int rid, nid, R;
    random_t seed;
    for (zi = 0; zi < nzone; zi++) {
        
        // only read from other zones' import buffer
        if (zi == this_zone) continue;
        R = (int)(probe_ncount[zi] / 3);

        if (R != 0) {

            for (ri = 0; ri < R; ri++) {
                rid = s->import_rat_info[zi][ri * 3];
                nid = s->import_rat_info[zi][ri * 3 + 1];
                seed = (random_t)(s->import_rat_info[zi][ri * 3 + 2]);
                                
                s->rat_position[rid] = nid;
                s->rat_count[nid]++;

                // update bitvector membership
                s->zone_rat_bitvector[rid] = 1;

                // update new rat seed
                s->rat_seed[rid] = seed;
            }
        }
    }


    for (zi = 0; zi < nzone; zi++) {
        // export_numrats = s->export_numrats[zi];
        // MPI_Wait(&(request[zi]), MPI_STATUS_IGNORE);
        // if (export_numrats != 0) {
        if (zi == this_zone) continue;
        MPI_Wait(&(request[zi]), MPI_STATUS_IGNORE);
        // }
    }


    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange node counts for boundary nodes between zones */
void exchange_node_states(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    graph_t *g = s->g;
    int nzone = g->nzone;
    int nnode = g->nnode;
    int this_zone = g->this_zone;
    MPI_Request request[nzone];

    // send to all other zones (async)
    int zi, ni, nid;
    for (zi = 0; zi < nzone; zi++) {

        int ncount = g->export_node_count[zi];

        if (ncount == 0) continue;

        for (ni = 0; ni < ncount; ni++) {
            nid = g->export_node_list[zi][ni];
            s->export_node_state[zi][ni] = s->rat_count[nid];
        }

        MPI_Isend(s->export_node_state[zi], ncount, MPI_INT, zi, zi, MPI_COMM_WORLD, &(request[zi]));
    }

    // receive from all other zones (sync)
    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->import_node_count[zi];

        if (ncount != 0) {
            MPI_Recv(s->import_node_state[zi], ncount, MPI_INT, zi, this_zone, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->export_node_count[zi];
        if (ncount != 0) {
            assert(ncount > 0);
            MPI_Wait(&(request[zi]), MPI_STATUS_IGNORE);
        }
    }


    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->import_node_count[zi];
        for (ni = 0; ni < ncount; ni++) {
            nid = g->import_node_list[zi][ni];
            s->rat_count[nid] = s->import_node_state[zi][ni];
        }
    }

    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange weights of nodes on boundaries */
void exchange_node_weights(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    graph_t *g = s->g;
    int nzone = g->nzone;
    int nnode = g->nnode;
    int this_zone = g->this_zone;
    MPI_Request request[nzone];

    // send to all other zones (async)
    int zi, ni, nid;
    for (zi = 0; zi < nzone; zi++) {

        int ncount = g->export_node_count[zi];
        if (ncount == 0) continue;

        assert(ncount > 0);

        for (ni = 0; ni < ncount; ni++) {
            nid = g->export_node_list[zi][ni];
            s->export_node_weight[zi][ni] = s->node_weight[nid];
        }

        MPI_Isend(s->export_node_weight[zi], ncount, MPI_DOUBLE, zi, zi, MPI_COMM_WORLD, &(request[zi]));
    }

    // receive from all other zones (sync)
    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->import_node_count[zi];
        if (ncount != 0) {
            MPI_Recv(s->import_node_weight[zi], ncount, MPI_DOUBLE, zi, this_zone, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    

    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->export_node_count[zi];
        if (ncount != 0) {
            assert(ncount > 0);
            MPI_Wait(&(request[zi]), MPI_STATUS_IGNORE);
        }
    }

    for (zi = 0; zi < nzone; zi++) {
        int ncount = g->import_node_count[zi];
        for (ni = 0; ni < ncount; ni++) {
            nid = g->import_node_list[zi][ni];
            s->node_weight[nid] = s->import_node_weight[zi][ni];
        }
    }

    FINISH_ACTIVITY(ACTIVITY_COMM);
}
#endif // MPI

/* Function suitable for sorting arrays of int's */
int comp_int(const void *ap, const void *bp) {
    int a = *(int *) ap;
    int b = *(int *) bp;
    int lt = a < b;
    int gt = a > b;
    return -lt + gt;
}

