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

    s->import_nid = calloc(nzone, sizeof(int*));
    ok = ok && s->import_nid != NULL;
    s->export_nid = calloc(nzone, sizeof(int*));
    ok = ok && s->export_nid != NULL;

    s->import_rid = calloc(nzone, sizeof(int*));
    ok = ok && s->import_rid != NULL;
    s->export_rid = calloc(nzone, sizeof(int*));
    ok = ok && s->export_rid != NULL;

    s->import_numrats = int_alloc(nzone);
    ok = ok && s->import_numrats != NULL;
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

    s->import_seed = calloc(nzone, sizeof(random_t *));
    ok = ok && s->import_seed != NULL;
    s->export_seed = calloc(nzone, sizeof(random_t *));
    ok = ok && s->export_seed != NULL;    
    
    s->import_node_weight = calloc(nzone, sizeof(double *));
    ok = ok && s->import_node_weight != NULL;
    s->export_node_weight = calloc(nzone, sizeof(double *));
    ok = ok && s->export_node_weight != NULL;
    
    s->zone_rat_list = int_alloc(nrat);
    ok = ok && s->zone_rat_list != NULL;

    s->zone_rat_bitvector = calloc(nrat, sizeof(unsigned char));
    ok = ok && s->zone_rat_bitvector != NULL;

    int i;

    // int local_nnode = s->g->local_node_count;
    for (i=0; i<nzone; i++) {
        s->import_nid[i] = int_alloc(nrat);
        s->export_nid[i] = int_alloc(nrat);
        s->import_rid[i] = int_alloc(nrat);
        s->export_rid[i] = int_alloc(nrat);
        s->import_seed[i] = rt_alloc(nrat);
        s->export_seed[i] = rt_alloc(nrat);
        s->import_node_state[i] = int_alloc(nnode);
        s->export_node_state[i] = int_alloc(nnode);

        s->export_rat_count[i] = int_alloc(nnode);
        s->import_rat_count[i] = int_alloc(nnode);
        s->import_node_weight[i] = double_alloc(nnode);
        s->export_node_weight[i] = double_alloc(nnode);

        ok = ok && 
             (s->import_nid[i] != NULL) &&
             (s->export_nid[i] != NULL) &&
             (s->import_rid[i] != NULL) && 
             (s->export_rid[i] != NULL) &&
             (s->import_rat_count[i] != NULL) &&
             (s->export_rat_count[i] != NULL) &&
             (s->import_node_state[i] != NULL) &&
             (s->export_node_state[i] != NULL) &&
             (s->import_seed[i] != NULL) &&
             (s->export_seed[i] != NULL) &&
             (s->import_node_weight[i] != NULL) &&
             (s->export_node_weight[i] != NULL);
    }

    if (!ok) return false;

    int ri;
    s->zone_rat_count = 0;

    for (ri=0; ri<nrat; ri++) {
        int ni = s->rat_position[ri];
        if (s->g->zone_id[ni] == zid) {
            s->zone_rat_list[s->zone_rat_count] = ri;
            assert(s->zone_rat_bitvector[ri] == 0);
            s->zone_rat_bitvector[ri] = 1;
            s->zone_rat_count++;
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

    int rat_position[nrat];
    MPI_Bcast(rat_position, nrat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);

    state_t *s = new_rats(g, nrat, global_seed);
    memcpy(s->rat_position, rat_position, nrat * sizeof(int));

    //s->rat_position = rat_position;
    seed_rats(s); // reseeding the rats
                                                  
    return s;
}

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(state_t *s) {
    int zi, i;
    int nzone = s->g->nzone;
    int nnode = s->g->nnode;
    int len = nzone-1;
    MPI_Request request[len];
    graph_t *g = s->g;

    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    int count;
    count = 0;
    for (zi = 1; zi < nzone; zi++){ //goes through all the processes
        MPI_Irecv(s->import_rat_count[zi], nnode, MPI_INT, zi, 0, MPI_COMM_WORLD, &(request[count++]));        
    }   

    MPI_Waitall(len, request, MPI_STATUS_IGNORE);

    for(i = 0; i < nnode; i++){
        int* zone_id_list = g->zone_id;
        int zid = zone_id_list[i]; // zone that the node is in
        if (zid != 0) {
            s->rat_count[i]=(s->import_rat_count)[zid][i];
        }
    }

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
    for (i=0; i<g->local_node_count; i++) {
        ni = g->local_node_list[i];
        s->export_rat_count[zi][ni] = s->rat_count[ni];
    }
    
    MPI_Request request;

    MPI_Isend(s->export_rat_count[zi], nnode, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
    
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Move rats between zones as they migrate */
void exchange_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    
    int zi, ri;
    int nzone = s->g->nzone;
    int this_zone = s->g->this_zone;
    int len = (nzone-1) * 4;
    int export_numrats, import_numrats;
    MPI_Request request[len];

    // send data to all the other zones (async)
    int count;
    count = 0;
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        export_numrats = s->export_numrats[zi];
        if (this_zone == 0) {
            // outmsg("--------------------------------------------------------\n");
            // outmsg("sending data from zone %d ---> zone %d\n", this_zone, zi);
            // outmsg("# of rats sending: %d\n", export_numrats);
        }
        
        
        
        
        MPI_Isend(&(export_numrats), 1, MPI_INT, zi, 0, MPI_COMM_WORLD, &(request[count*4]));
        MPI_Isend(s->export_nid[zi], export_numrats, MPI_INT, zi, 1, MPI_COMM_WORLD, &(request[count*4+1]));
        MPI_Isend(s->export_rid[zi], export_numrats, MPI_INT, zi, 2, MPI_COMM_WORLD, &(request[count*4+2]));
        MPI_Isend(s->export_seed[zi], export_numrats, MPI_UNSIGNED, zi, 3, MPI_COMM_WORLD, &(request[count*4+3]));
        count++;
    }
    
    // receive data from all the other zones (synch)
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        MPI_Recv(&(s->import_numrats[zi]), 1, MPI_INT, zi, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        import_numrats = s->import_numrats[zi];
        MPI_Recv(s->import_nid[zi], import_numrats, MPI_INT, zi, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(s->import_rid[zi], import_numrats, MPI_INT, zi, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(s->import_seed[zi], import_numrats, MPI_UNSIGNED, zi, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(len, request, MPI_STATUS_IGNORE);
    
    // move received data
    // int old_rat_count = s->zone_rat_count;
    int rid, nid, R;
    for (zi = 0; zi < nzone; zi++) {
        
        // only read from other zones' import buffer
        if (zi == this_zone) continue;
        R = s->import_numrats[zi];
        s->zone_rat_count += R;

        for (ri = 0; ri < R; ri++) {
            rid = s->import_rid[zi][ri];
            nid = s->import_nid[zi][ri];

            s->rat_position[rid] = nid;
            s->rat_count[nid]++;

            // update bitvector membership
            assert(s->zone_rat_bitvector[rid] == 0);
            s->zone_rat_bitvector[rid] = 1;

            // update new rat seed
            s->rat_seed[rid] = s->import_seed[zi][ri];
        }
    }

    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange node counts for boundary nodes between zones */
void exchange_node_states(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    int len, this_zone;
    graph_t *g = s->g;
    int nzone = g->nzone;
    int nnode = g->nnode;
    len = nzone-1;
    this_zone = g->this_zone;
    MPI_Request request[len];

    // send to all other zones (async)
    int count;
    int zi, ni, nid;
    count = 0;
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        int ncount;
        ncount = g->export_node_count[zi];
        for (ni = 0; ni < ncount; ni++) {
            nid = g->export_node_list[zi][ni];
            s->export_node_state[zi][nid] = s->rat_count[nid];
        }
        MPI_Isend(s->export_node_state[zi], nnode, MPI_INT, zi, 0, MPI_COMM_WORLD, &(request[count++]));
    }

    // receive from all other zones (sync)
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        int ncount;
        ncount = g->import_node_count[zi];
        MPI_Recv(s->import_node_state[zi], nnode, MPI_INT, zi, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (ni = 0; ni < ncount; ni++) {
            nid = g->import_node_list[zi][ni];
            s->rat_count[nid] = s->import_node_state[zi][nid];
        }
    }

    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange weights of nodes on boundaries */
void exchange_node_weights(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    int len, this_zone;
    graph_t *g = s->g;
    int nzone = g->nzone;
    int nnode = g->nnode;
    len = nzone-1;
    this_zone = g->this_zone;
    MPI_Request request[len];

    // send to all other zones (async)
    int count;
    int zi, ni, nid;
    count = 0;
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        int ncount;
        ncount = g->export_node_count[zi];
        for (ni = 0; ni < ncount; ni++) {
            nid = g->export_node_list[zi][ni];
            s->export_node_weight[zi][nid] = s->node_weight[nid];
        }
        MPI_Isend(s->export_node_weight[zi], nnode, MPI_INT, zi, 0, MPI_COMM_WORLD, &(request[count++]));
    }

    // receive from all other zones (sync)
    for (zi = 0; zi < nzone; zi++) {
        if (zi == this_zone) continue;
        int ncount;
        ncount = g->import_node_count[zi];
        MPI_Recv(s->import_node_weight[zi], nnode, MPI_INT, zi, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (ni = 0; ni < ncount; ni++) {
            nid = g->import_node_list[zi][ni];
            s->node_weight[nid] = s->import_node_weight[zi][nid];
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

