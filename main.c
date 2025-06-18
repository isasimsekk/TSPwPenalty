#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <stdbool.h>

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    int* items;
    size_t len;
    size_t capacity;
} Tour;

typedef struct {
    int* neighbors;
    int count;
} NeighborEntry;

typedef struct {
    NeighborEntry* entries;
    int num_cities;
} NeighborList;


void read_input(const char *path, int *penalty, Point **coords, int *n);
int** build_distance_matrix(const Point *coords, int n);
void free_distance_matrix(int **D, int n);
NeighborList build_neighbor_list(const Point *coords, int n, int k, int w);
void free_neighbor_list(NeighborList nbrs);
Tour create_tour(int capacity);
void tour_add(Tour *tour, int city);
void tour_pop(Tour *tour, int index);
void free_tour(Tour *tour);
Tour initial_tour_nn(const Point *coords, int n, int** D);
long long tour_cost(const Tour *tour, int** D, int penalty, int total_cities);
Tour two_opt_neighbor_dlb(Tour tour, int** D, const NeighborList* nbrs, int total_cities, double time_limit);
Tour skip_penalty(Tour tour, int** D, int penalty);
Tour optimize(const Point *coords, int n, int penalty, int k_neighbors, int window, double time_2opt);
void write_output(const char *path, const Tour *tour, long long cost);

typedef struct {
    int index;
    int value;
} SortItem;

int compare_sort_items(const void *a, const void *b) {
    return ((SortItem*)a)->value - ((SortItem*)b)->value;
}

int compare_doubles(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

bool has_duplicates(const Tour* tour, int n) {
    bool* seen = calloc(n, sizeof(bool));
    for (size_t i = 0; i < tour->len; ++i) {
        if (seen[tour->items[i]]) {
            free(seen);
            return true;
        }
        seen[tour->items[i]] = true;
    }
    free(seen);
    return false;
}

int main() {
    char *input_file = "test-input-2.txt";
    char *output_file = "output.txt";
    int neighbors = 120;
    int window = 400;
    double time_2opt = 2000.0;


    printf("-> Settings: Input='%s', Output='%s', Neighbors=%d, Window=%d, Time=%.1f\n",
           input_file, output_file, neighbors, window, time_2opt);

    srand(time(NULL));

    int penalty;
    Point *coords = NULL;
    int n_coords = 0;
    read_input(input_file, &penalty, &coords, &n_coords);
    printf("-> Input file read: %d cities, Penalty Value: %d\n", n_coords, penalty);

    printf("\n neighbor = %d,  window = %d", neighbors, window);

    printf("\n--- Optimization Starting ---\n");
    Tour best_tour = optimize(coords, n_coords, penalty, neighbors, window, time_2opt);
    printf("--- Optimization Complete ---\n\n");

    int** D = build_distance_matrix(coords, n_coords);
    long long cost = tour_cost(&best_tour, D, penalty, n_coords);

    printf("-> Results are being written to '%s'...\n", output_file);
    write_output(output_file, &best_tour, cost);
    printf("-> Final Cost: %lld, Visited Cities: %zu\n", cost, best_tour.len);

    free(coords);
    free_distance_matrix(D, n_coords);
    free_tour(&best_tour);

    return 0;
}


void read_input(const char *path, int *penalty, Point **coords, int *n) {
    FILE *f = fopen(path, "r");
    if (!f) {
        perror("Could not open input file");
        exit(1);
    }

    fscanf(f, "%d", penalty);

    int capacity = 100;
    *coords = malloc(capacity * sizeof(Point));
    *n = 0;
    int idx, x, y;
    while (fscanf(f, "%d %d %d", &idx, &x, &y) == 3) {
        if (idx >= capacity) {
            capacity = idx * 2;
            *coords = realloc(*coords, capacity * sizeof(Point));
        }
        (*coords)[idx] = (Point){x, y};
        if (idx + 1 > *n) {
            *n = idx + 1;
        }
    }
    fclose(f);
}

int** build_distance_matrix(const Point *coords, int n) {
    int** D = malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        D[i] = malloc(n * sizeof(int));
    }
    for (int i = 0; i < n; i++) {
        D[i][i] = 0;
        for (int j = i + 1; j < n; j++) {
            double dx = coords[i].x - coords[j].x;
            double dy = coords[i].y - coords[j].y;
            int d = round(hypot(dx, dy));
            D[i][j] = D[j][i] = d;
        }
    }
    return D;
}

void free_distance_matrix(int **D, int n) {
    if (!D) return;
    for (int i = 0; i < n; i++) {
        free(D[i]);
    }
    free(D);
}

typedef struct {
    int city;
    double dist;
} Candidate;

int compare_candidates(const void *a, const void *b) {
    double diff = ((Candidate*)a)->dist - ((Candidate*)b)->dist;
    return (diff > 0) - (diff < 0);
}

NeighborList build_neighbor_list(const Point *coords, int n, int k, int w) {
    NeighborList nbrs;
    nbrs.num_cities = n;
    nbrs.entries = malloc(n * sizeof(NeighborEntry));

    SortItem* idx_x = malloc(n * sizeof(SortItem));
    SortItem* idx_y = malloc(n * sizeof(SortItem));
    for(int i=0; i<n; ++i) {
        idx_x[i] = (SortItem){i, coords[i].x};
        idx_y[i] = (SortItem){i, coords[i].y};
    }
    qsort(idx_x, n, sizeof(SortItem), compare_sort_items);
    qsort(idx_y, n, sizeof(SortItem), compare_sort_items);

    int* pos_x = malloc(n * sizeof(int));
    int* pos_y = malloc(n * sizeof(int));
    for(int i=0; i<n; ++i) {
        pos_x[idx_x[i].index] = i;
        pos_y[idx_y[i].index] = i;
    }

    for (int city = 0; city < n; ++city) {
        bool* in_cand = calloc(n, sizeof(bool));
        Candidate* candidates = malloc(n * sizeof(Candidate));
        int cand_count = 0;

        int px = pos_x[city];
        int py = pos_y[city];

        int lo_x = (px - w > 0) ? px - w : 0;
        int hi_x = (px + w < n) ? px + w : n;
        for (int i = lo_x; i < hi_x; ++i) {
            int other = idx_x[i].index;
            if (other != city && !in_cand[other]) {
                in_cand[other] = true;
            }
        }

        int lo_y = (py - w > 0) ? py - w : 0;
        int hi_y = (py + w < n) ? py + w : n;
        for (int i = lo_y; i < hi_y; ++i) {
            int other = idx_y[i].index;
            if (other != city && !in_cand[other]) {
                in_cand[other] = true;
            }
        }

        for(int j=0; j<n; ++j){
            if(in_cand[j]){
                 double dx = coords[city].x - coords[j].x;
                 double dy = coords[city].y - coords[j].y;
                 candidates[cand_count++] = (Candidate){j, hypot(dx, dy)};
            }
        }

        qsort(candidates, cand_count, sizeof(Candidate), compare_candidates);

        int n_neighbors = (cand_count < k) ? cand_count : k;
        nbrs.entries[city].count = n_neighbors;
        nbrs.entries[city].neighbors = malloc(n_neighbors * sizeof(int));
        for (int i = 0; i < n_neighbors; ++i) {
            nbrs.entries[city].neighbors[i] = candidates[i].city;
        }
        free(in_cand);
        free(candidates);
    }

    free(idx_x);
    free(idx_y);
    free(pos_x);
    free(pos_y);
    return nbrs;
}

void free_neighbor_list(NeighborList nbrs) {
    for (int i = 0; i < nbrs.num_cities; i++) {
        free(nbrs.entries[i].neighbors);
    }
    free(nbrs.entries);
}


Tour create_tour(int capacity) {
    Tour t;
    t.len = 0;
    t.capacity = capacity;
    t.items = malloc(capacity * sizeof(int));
    return t;
}

void tour_add(Tour *tour, int city) {
    if (tour->len == tour->capacity) {
        tour->capacity *= 2;
        tour->items = realloc(tour->items, tour->capacity * sizeof(int));
    }
    tour->items[tour->len++] = city;
}

void tour_pop(Tour *tour, int index) {
    if (index < 0 || index >= tour->len) return;
    for (size_t i = index; i < tour->len - 1; ++i) {
        tour->items[i] = tour->items[i + 1];
    }
    tour->len--;
}

void free_tour(Tour *tour) {
    free(tour->items);
    tour->items = NULL;
    tour->len = tour->capacity = 0;
}

Tour initial_tour_nn(const Point *coords, int n, int** D) {
    Tour tour = create_tour(n);
    bool *used = calloc(n, sizeof(bool));

    int start = rand() % n;
    tour_add(&tour, start);
    used[start] = true;
    int cur = start;

    for (int i = 0; i < n - 1; i++) {
        int next_city = -1;
        int min_dist = INT_MAX;
        for (int j = 0; j < n; j++) {
            if (!used[j] && D[cur][j] < min_dist) {
                min_dist = D[cur][j];
                next_city = j;
            }
        }
        tour_add(&tour, next_city);
        used[next_city] = true;
        cur = next_city;
    }
    free(used);
    return tour;
}


long long tour_cost(const Tour *tour, int** D, int penalty, int total_cities) {
    long long cost = 0;
    if (tour->len == 0) {
        return (long long)penalty * total_cities;
    }
    for (size_t i = 0; i < tour->len - 1; i++) {
        cost += D[tour->items[i]][tour->items[i+1]];
    }
    cost += D[tour->items[tour->len - 1]][tour->items[0]];
    cost += (long long)(total_cities - tour->len) * penalty;
    return cost;
}

Tour two_opt_neighbor_dlb(Tour tour, int** D, const NeighborList* nbrs, int total_cities, double time_limit) {
    clock_t start_time = clock();
    double time_limit_ticks = time_limit * CLOCKS_PER_SEC;

    int *inv = malloc(total_cities * sizeof(int));
    bool *dlb = calloc(total_cities, sizeof(bool));

    bool improved = true;

    while (improved && (clock() - start_time) < time_limit_ticks) {
        improved = false;

        for (int i = 0; i < total_cities; ++i) inv[i] = -1;
        for (size_t i = 0; i < tour.len; ++i) {
            int city = tour.items[i];
            if (city >= 0 && city < total_cities) {
                inv[city] = i;
            }
        }

        for (size_t u_pos = 0; u_pos < tour.len; ++u_pos) {
            int u_city = tour.items[u_pos];
            if (dlb[u_city]) continue;

            bool found_improvement = false;
            NeighborEntry u_neighbors = nbrs->entries[u_city];

            for (int i = 0; i < u_neighbors.count; ++i) {
                int v_city = u_neighbors.neighbors[i];
                int v_pos = inv[v_city];

                if (v_pos == -1) continue; // Skip removed or invalid cities

                if (abs((int)u_pos - v_pos) <= 1 ||
                    (u_pos == 0 && v_pos == (int)tour.len - 1) ||
                    (u_pos == (int)tour.len - 1 && v_pos == 0)) {
                    continue;
                }

                int a = tour.items[u_pos];
                int b = tour.items[(u_pos + 1) % tour.len];
                int c = tour.items[v_pos];
                int d = tour.items[(v_pos + 1) % tour.len];

                int old_cost = D[a][b] + D[c][d];
                int new_cost = D[a][c] + D[b][d];

                if (new_cost < old_cost) {

                    int i = u_pos;
                    int j = v_pos;
                    if (i > j) {
                        int temp = i;
                        i = j;
                        j = temp;
                    }

                    int start_rev = i + 1;
                    int end_rev = j;


                    if (start_rev < end_rev && end_rev < (int)tour.len) {
                        while (start_rev < end_rev) {
                            int temp = tour.items[start_rev];
                            tour.items[start_rev] = tour.items[end_rev];
                            tour.items[end_rev] = temp;
                            start_rev++;
                            end_rev--;
                        }
                        improved = true;
                        found_improvement = true;
                        for (int k = 0; k < total_cities; ++k) dlb[k] = false; // Reset all
                        break;
                    }
                }
            }

            if (found_improvement) break;  // Outer loop restart
            dlb[u_city] = true;  // Mark as no-improvement
        }
    }

    free(inv);
    free(dlb);
    return tour;
}


Tour skip_penalty(Tour tour, int** D, int penalty) {
    while (true) {
        long long best_sav = penalty; // Initial best saving is just the penalty for skipping one city
        int best_i = -1;
        int m = tour.len;
        if (m <= 2) break; // Cannot remove cities if 2 or fewer remain (would break tour structure)

        for (int i = 0; i < m; i++) {
            int p_idx = (i == 0) ? m - 1 : i - 1; // Previous city index
            int c_idx = i;                      // Current city index to consider removing
            int nx_idx = (i + 1) % m;           // Next city index

            int p = tour.items[p_idx];
            int c = tour.items[c_idx];
            int nx = tour.items[nx_idx];

            // Calculate saving: (cost of edges involving c) - (cost of new edge if c is removed)
            // Plus the penalty value recovered
            long long sav = (long long)D[p][c] + D[c][nx] - D[p][nx];
            if (sav > best_sav) { // Is this saving better than just the penalty?
                best_sav = sav;
                best_i = i;
            }
        }

        if (best_i < 0) { // No city found that improves cost by skipping it
            break;
        }
        tour_pop(&tour, best_i); // Remove the city that provides the best saving
    }
    return tour;
}



Tour optimize(const Point *coords, int n, int penalty, int k_neighbors, int window, double time_2opt) {
    printf("  [Optimize] Building distance matrix and neighbor list...\n");
    int** D = build_distance_matrix(coords, n);
    NeighborList nbrs = build_neighbor_list(coords, n, k_neighbors, window);

    printf("  [Optimize] Generating initial tours...\n");
    Tour best_t = create_tour(n);
    long long best_c = __LONG_LONG_MAX__;

    for (int i = 0; i < 50; i++) {
        Tour t0 = initial_tour_nn(coords, n, D);
        if (has_duplicates(&t0, n)) {
            printf("exit i = %d", i);
            exit(0);
        }
        long long c0 = tour_cost(&t0, D, penalty, n);
        printf("    - nn tour attempt %d, cost: %lld\n", i + 1, c0);
        if (c0 < best_c) {
            free_tour(&best_t); // Free the old best tour
            best_t = t0;        // Assign the new best tour
            best_c = c0;
            printf("    -> New best initial tour found, cost: %lld\n", best_c);
        } else {
            free_tour(&t0); // Free the less optimal tour
        }
    }

        printf("  [Optimize] Starting 2-opt local search (Time: %.1fs)...\n", time_2opt);

        best_t = two_opt_neighbor_dlb(best_t, D, &nbrs, n, time_2opt);
        if (has_duplicates(&best_t, n)) {
            printf("exit 2");
            exit(0);
        }
        best_c = tour_cost(&best_t, D, penalty, n);
        printf("    - Cost after 2-opt: %lld\n", best_c);

        printf("  [Optimize] Skipping penalized cities (skip_penalty)...\n");
        best_t = skip_penalty(best_t, D, penalty);
        if (has_duplicates(&best_t, n)) {
            printf("exit 3");
            exit(0);
        }
        best_c = tour_cost(&best_t, D, penalty, n);
        printf("    - Cost after skip_penalty: %lld, Visited: %zu\n", best_c, best_t.len);

        // One more 2-opt pass after skipping cities, to re-optimize the new shorter tour
        printf("  [Optimize] Startin 2-opt local search (Time: %.1fs)...\n" ,time_2opt ); // i+2 for the second 2-opt in the iteration
        best_t = two_opt_neighbor_dlb(best_t, D, &nbrs, n, time_2opt);
        if (has_duplicates(&best_t, n)) {
            printf("exit 4");
            exit(0);
        }

        best_c = tour_cost(&best_t, D, penalty, n);
        printf("    - Cost after 2-opt: %lld\n", best_c);

    free_distance_matrix(D, n);
    free_neighbor_list(nbrs);

    return best_t;
}

void write_output(const char *path, const Tour *tour, long long cost) {
    FILE *f = fopen(path, "w");
    if (!f) {
        perror("Could not write output file");
        return;
    }
    fprintf(f, "%lld %zu\n", cost, tour->len);
    for (size_t i = 0; i < tour->len; i++) {
        fprintf(f, "%d\n", tour->items[i]);
    }
    fclose(f);
}
