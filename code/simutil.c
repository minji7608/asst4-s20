#include "crun.h"

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
    return true;
}

//TODO: Implement these communication-support functions
#if MPI
/* Called by process 0 to distribute rat state to all nodes */
void send_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Called by other nodes to get rat state from master and set up state data structure */
state_t *get_rats(graph_t *g, random_t global_seed) {
    int nrat = 0;
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /*
      Your code should go here.
      It should receive information about all rats, including how many there are,
      from process 0.
     */
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    state_t *s = new_rats(g, nrat, global_seed);
    return s;
}

/* Called by process 0 to collect node states from all other processes */
void gather_node_state(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Called by other processes to send their node states to process 0 */
void send_node_state(state_t *s) {
    START_ACTIVITY(ACTIVITY_GLOBAL_COMM);
    /* Your code should go here */
    FINISH_ACTIVITY(ACTIVITY_GLOBAL_COMM);
}

/* Move rats between zones as they migrate */
void exchange_rats(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange node counts for boundary nodes between zones */
void exchange_node_states(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

/* Exchange weights of nodes on boundaries */
void exchange_node_weights(state_t *s) {
    START_ACTIVITY(ACTIVITY_COMM);
    /* Your code should go here */
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

