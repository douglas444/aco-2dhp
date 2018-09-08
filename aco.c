#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "aco.h"

enum direction
{
    LEFT = 0,
    RIGHT = 1,
    STRAIGHT = 2
};

enum pm_type
{
    ORIGINAL = 0,
    INVERSE = 1
};

struct coord
{
    int x;
    int y;
};

struct ant
{
    int energy;
    int *energy_by_edge;
    struct coord *positions;
};

struct pm_config
{
    int amino_acid_index;
    struct coord curr;
    struct coord prev;
    struct coord next;
    struct coord f;
    struct coord c;
    enum pm_type pm_type;

};

struct candidate
{
    int heuristic;
    struct coord position;
    struct coord move;
    int readjusted;

};

void* smalloc(int mem_size)
/* ====================================
 * Allocates memory
 * ====================================
 */
{
    void *mem_pos = (void*) malloc(mem_size);

    if (mem_pos == NULL)
    {
        printf("ERROR: aco.c/smalloc(): \"Unable to allocate memory\"\n");
        exit(1);
    }
    else
    {
        return mem_pos;
    }
}

void init_solution(struct solution *solution, int seq_len)
/* ====================================================
 * Allocates memory to struct solution structure variables
 * ====================================================
 */
{
    solution->directions = (char*) smalloc(sizeof(char) * (seq_len + 1));
}

void free_solution(struct solution solution)
/* ===========================================
 * Free memory of struct solution structure variables
 * ===========================================
 */
{
    free(solution.directions);
}

void init_ant(struct ant *ant, int seq_len)

/* ===============================================
 * Allocates memory to struct ant structure variables
 * ===============================================
 */
{
    ant->energy_by_edge = (int*) smalloc(sizeof(int) * (seq_len - 1));
    ant->positions = (struct coord*) smalloc(sizeof(struct coord) * seq_len);
}

void free_ant(struct ant ant)
/* =======================================
 * Frees memory of struct ant structure variables
 * =======================================
 */
{
    free(ant.energy_by_edge);
    free(ant.positions);
}

void init_variables
(
    struct aco_config aco_config,
    struct pm_config **pm_configs,
    struct ant *pm_best_ant,
    struct ant *pm_ant,
    struct ant *best_ant,
    int seq_len,
    int **best_ant_by_edge,
    int ***lattice,
    double ***pheromone,
    struct ant **ants
)
/* =======================================
 * Allocates all aco.c exclusive variables
 * =======================================
 */
{
    int i;
    int j;

    *pm_configs = (struct pm_config*) smalloc(sizeof(struct pm_config) * 4 * (seq_len - 2));

    init_ant(pm_best_ant, seq_len);
    init_ant(pm_ant, seq_len);
    init_ant(best_ant, seq_len);

    *pheromone = (double**) smalloc(sizeof(double*) * seq_len);
    for (i = 0; i < seq_len; ++i)
    {
        (*pheromone)[i] = (double*) smalloc(sizeof(double) * 3);
        for (j = 0; j < 3; ++j)
        {
            (*pheromone)[i][j] = aco_config.ini_pheromone;
        }
    }

    *best_ant_by_edge = (int*) smalloc(sizeof(int) * (seq_len - 1));
    for (i = 0; i < seq_len - 1; ++i)
    {
        (*best_ant_by_edge)[i] = -1;
    }

    *lattice = (int**) smalloc(sizeof(int*) * (2 * seq_len + 1));
    for (i = 0; i < 2 * seq_len + 1; ++i)
    {
        (*lattice)[i] = (int*) smalloc(sizeof(int) * (2 * seq_len + 1));
        for (j = 0; j < 2 * seq_len + 1; ++j)
        {
            (*lattice)[i][j] = -1;
        }
    }

    *ants = (struct ant*) malloc(sizeof(struct ant) * aco_config.population);
    for (i = 0; i < aco_config.population; ++i)
    {
        init_ant(&((*ants)[i]), seq_len);
    }

}

void free_variables
(
    struct aco_config aco_config,
    int seq_len,
    int **lattice,
    int *best_ant_by_edge,
    struct ant best_ant,
    struct ant pm_best_ant,
    struct ant pm_ant,
    double **pheromone,
    struct pm_config *pm_configs,
    struct ant *ants
)
/* ==========================================
 * Frees all aco.c exclusive variables memory
 * ==========================================
 */
{
    int i;

    for (i = 0; i < seq_len; ++i)
    {
        free(pheromone[i]);
    }
    for (i = 0; i < 2 * seq_len + 1; ++i)
    {
        free(lattice[i]);
    }
    for (i = 0; i < aco_config.population; ++i)
    {
        free_ant(ants[i]);
    }

    free(ants);
    free(best_ant_by_edge);
    free(lattice);
    free(pheromone);
    free_ant(pm_best_ant);
    free_ant(pm_ant);
    free_ant(best_ant);
    free(pm_configs);
}

struct coord create_new_coord(int x, int y)
/* ===========================================================
 * Returns a new variables of type struct coord with the given values
 * ===========================================================
 */
{
    struct coord coord;
    coord.x = x;
    coord.y = y;
    return coord;
}

struct coord subtract_coord
(
    struct coord c1,
    struct coord c2
)
/* =======================================================
 * Calculates the difference between two given coordinates
 * =======================================================
 */
{
    struct coord c3;

    c3.x = c1.x - c2.x;
    c3.y = c1.y - c2.y;

    return c3;
}

int lateral_adj(struct coord c1, struct coord c2)
/* ===============================================
 * Check if two coordinates are laterally adjacent
 * ===============================================
 */
{
    if (abs(c1.x - c2.x) + abs(c1.y - c2.y) == 1) {
        return 1;
    } else {
        return 0;
    }

}

int abs_direction_by_move(struct coord move)
/* =====================================
 * Calculates absolute direction by move
 * =====================================
 */
{
    if (move.x == 0)
    {
        if (move.y == 1 || move.y == -1)
        {
            return STRAIGHT;
        }
        else
        {
            printf("ERROR: aco.c/abs_direction_by_move(): \"Invalid value for parameter move\"\n");
            exit(1);
        }
    }
    else if (move.y == 0)
    {
        if (move.x == 1)
        {
            return LEFT;
        }
        else if (move.x == -1)
        {
            return RIGHT;
        }
        else
        {
            printf("ERROR: aco.c/abs_direction_by_move(): \"Invalid value for parameter move\"\n");
            exit(1);
        }
    }
    else
    {
        printf("ERROR: aco.c/abs_direction_by_move(): \"Invalid value for parameter move\"\n");
        exit(1);
    }
}

int direction_by_move(struct coord prev_move, struct coord move)
/* ===========================================================
 * Calculates relative direction between two consecutive moves
 * ===========================================================
 */
{

    if (prev_move.x == move.x && prev_move.y == move.y)
    {
        return STRAIGHT;
    }
    else if (prev_move.y == 0)
    {
        if (move.y == prev_move.x)
        {
            return LEFT;
        }
        else if (abs(move.y) == abs(prev_move.x))
        {
            return RIGHT;
        }
        else
        {
            printf("ERROR: aco.c/direction_by_move(): \"Invalid values for parameters\"\n");
            exit(1);
        }
    }
    else
    {
        if (move.x == prev_move.y)
        {
            return RIGHT;
        }
        else if (abs(move.x) == abs(prev_move.y))
        {
            return LEFT;
        }
        else
        {
            printf("ERROR: aco.c/direction_by_move(): \"Invalid values for parameters\"\n");
            exit(1);
        }
    }
}

struct coord straight(struct coord prev_move)
/* =================================================================
 * Calculates straight move corresponding to the given previous move
 * =================================================================
 */
{
    struct coord move;

    if (prev_move.x == 1)
    {
        move.x = 1;
        move.y = 0;
    }
    else if (prev_move.x == -1)
    {
        move.x = -1;
        move.y = 0;
    }
    else if (prev_move.y == -1)
    {
        move.x = 0;
        move.y = -1;
    }
    else if (prev_move.y == 1)
    {
        move.x = 0;
        move.y = 1;
    }
    else
    {
        printf("ERROR: aco.c/straight(): \"Invalid value for parameter prev_move\"\n");
        exit(1);
    }

    return move;
}

struct coord left(struct coord prev_move)
/* =============================================================
 * Calculates left move corresponding to the given previous move
 * =============================================================
 */
{
    struct coord move;

    if (prev_move.x == 1)
    {
        move.x = 0;
        move.y = 1;
    }
    else if (prev_move.x == -1)
    {
        move.x = 0;
        move.y = -1;
    }
    else if (prev_move.y == -1)
    {
        move.x = 1;
        move.y = 0;
    }
    else if (prev_move.y == 1)
    {
        move.x = -1;
        move.y = 0;
    }
    else
    {
        printf("ERROR: aco.c/left(): \"Invalid value for parameter prev_move\"\n");
        exit(1);
    }

    return move;
}

struct coord right(struct coord prev_move)
/* ==============================================================
 * Calculates right move corresponding to the given previous move
 * ==============================================================
 */
{
    struct coord move;

    if (prev_move.x == 1)
    {
        move.x = 0;
        move.y = -1;
    }
    else if (prev_move.x == -1)
    {
        move.x = 0;
        move.y = 1;
    }
    else if (prev_move.y == -1)
    {
        move.x = -1;
        move.y = 0;
    }
    else if (prev_move.y == 1)
    {
        move.x = 1;
        move.y = 0;
    }
    else
    {
        printf("ERROR: aco.c/right() \"Invalid value for parameter prev_move\"\n");
        exit(1);
    }

    return move;
}

int calculate_heuristic
(
    int **lattice,
    int amino_acid_index,
    struct coord pos,
    enum polarity *seq
)
/* =================================================================================
 * Calculates the number of H-H contacts if a H amino-acid occupy the given position
 * =================================================================================
 */
{
    int heuristic_value = 0;
    int right_neighbor = lattice[pos.x + 1][pos.y];
    int left_neighbor = lattice[pos.x - 1][pos.y];
    int down_neighbor = lattice[pos.x][pos.y - 1];
    int up_neighbor = lattice[pos.x][pos.y + 1];

    if (right_neighbor >= 0 && abs(right_neighbor - amino_acid_index) > 1
        && seq[right_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (left_neighbor >= 0 && abs(left_neighbor - amino_acid_index) > 1
        && seq[left_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (up_neighbor >= 0 && abs(up_neighbor - amino_acid_index) > 1
        && seq[up_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (down_neighbor >= 0 && abs(down_neighbor - amino_acid_index) > 1
        && seq[down_neighbor] == H)
    {
        ++heuristic_value;
    }

    return heuristic_value;
}

int random_select(double *probabilities, int len)
/* =====================================================================
 * Randomly chooses a integer between 0 and len with given probabilities
 * =====================================================================
 */
{
    int i = 0;
    int result = -1;
    double cumulative_probability = 0;
    double r = ((double)rand()/RAND_MAX);

    while(result == -1)
    {
        cumulative_probability += probabilities[i];
        if(r <= cumulative_probability || i == len - 1)
        {
            result = i;
        }
        ++i;
    }
    return result;
}

void pheromone_deposit
(
    double **pheromone,
    struct ant ant,
    int seq_len,
    int best_energy
)
/* ==============================
 * Aco pheromone deposit function
 * ==============================
 */
{
    int i;
    struct coord move, prev_move;
    int direction;

    for (i = 0; i < seq_len - 1; ++i)
    {
        move = subtract_coord(ant.positions[i + 1], ant.positions[i]);
        if (i == 0)
        {
            direction = abs_direction_by_move(move);
        }
        else
        {
            direction = direction_by_move(prev_move, move);
        }
        if (best_energy != 0)
        {
            pheromone[i][direction] += (double) ant.energy / (best_energy * best_energy * best_energy);
        }
        prev_move = move;
    }
}

void pheromone_evaporation
(
    double **pheromone,
    int seq_len,
    double persistence
)
/* ==================================
 * Aco pheromone evaporation function
 * ==================================
 */
{
    int i, j;

    for (i = 0; i < seq_len - 1; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            pheromone[i][j] *= persistence;
        }
    }
}

int move_amino_acid
(
    int current_energy,
    int **lattice,
    enum polarity *seq,
    int amino_acid_index,
    struct coord src,
    struct coord dest
)
/* ==========================================================================
 * Moves a amino-acid to another position in the lattice, handling the energy
 * ==========================================================================
 */
{

    int temp;

    if (seq[amino_acid_index] == H)
    {
        current_energy += calculate_heuristic(lattice, amino_acid_index, src, seq);
        current_energy -= calculate_heuristic(lattice, amino_acid_index, dest, seq);
    }

    temp = lattice[src.x][src.y];
    lattice[src.x][src.y] = lattice[dest.x][dest.y];
    lattice[dest.x][dest.y] = temp;

    return current_energy;
}

int apply_pm
(
    struct ant ant,
    struct pm_config config,
    enum polarity *seq,
    int **lattice,
    int seq_len,
    struct coord *ant_positions
)
/* ========================================================
 * Applies the pull-move configuration to the given protein
 * ========================================================
 */
{

    int i;
    int last_modified;
    int previous_index;
    int before_previous_index;
    struct coord temp_coord;

    if (config.pm_type == ORIGINAL)
    {
        previous_index = config.amino_acid_index - 1;
        before_previous_index = config.amino_acid_index - 2;
        last_modified = 0;
    }
    else
    {
        previous_index = config.amino_acid_index + 1;
        before_previous_index = config.amino_acid_index + 2;
        last_modified = seq_len - 1;
    }

    if (config.c.x == config.prev.x && config.c.y == config.prev.y)
    {
        ant.energy = move_amino_acid(ant.energy, lattice, seq,
                                     config.amino_acid_index,
                                     config.curr, config.f);

        ant.positions[config.amino_acid_index] = config.f;
        config.curr =  ant.positions[config.amino_acid_index];
    }
    else
    {
        ant.energy = move_amino_acid(ant.energy, lattice, seq,
                                     previous_index, config.prev,
                                     config.c);

        ant.positions[previous_index] = config.c;

        ant.energy = move_amino_acid(ant.energy, lattice, seq,
                                     config.amino_acid_index,
                                     config.curr, config.f);

        ant.positions[config.amino_acid_index] = config.f;

        i = before_previous_index;

        while ((config.pm_type == ORIGINAL && i < seq_len - 1 && i >= 0 &&
                lateral_adj(ant.positions[i], ant.positions[i + 1]) == 0) ||
                (config.pm_type == INVERSE && i > 0 && i <= seq_len - 1 &&
                 lateral_adj(ant.positions[i], ant.positions[i - 1]) == 0))
        {
            ant.energy = move_amino_acid(ant.energy, lattice, seq, i,
                                         ant.positions[i], config.curr);

            temp_coord = ant.positions[i];
            ant.positions[i] = config.curr;
            config.curr = config.prev;
            config.prev = temp_coord;

            last_modified = i;

            if (config.pm_type == ORIGINAL)
            {
                --i;
            }
            else if (config.pm_type == INVERSE)
            {
                ++i;
            }

        }

    }

    //DEVOLVES LATTICE TO ORIGINAL STATE

    int start;
    int end;

    if (config.pm_type == ORIGINAL)
    {
        start = last_modified;
        end = config.amino_acid_index;
    }
    else
    {
        end = last_modified;
        start = config.amino_acid_index;
    }

    for (i = start; i <= end; ++i)
    {
        lattice[ant.positions[i].x][ant.positions[i].y] = -1;
    }
    for (i = start; i <= end; ++i)
    {
        lattice[ant_positions[i].x][ant_positions[i].y] = i;
    }

    return ant.energy;
}

int generate_pm_config
(
    struct pm_config *config,
    int **lattice,
    enum pm_type pm_type,
    struct coord curr,
    struct coord prev,
    struct coord next,
    struct coord direction
)
/* =====================================================================================
 * Generates a pull-move configuration for a specific direction of a specific amino-acid
 * =====================================================================================
 */
{
    /* F exists? */
    if (lattice[next.x + direction.x][next.y + direction.y] == -1 &&
        next.x + direction.x != curr.x && next.y + direction.y != curr.y)
    {
        config->f.x = next.x + direction.x;
        config->f.y = next.y + direction.y;
        config->c.x = config->f.x + curr.x - next.x;
        config->c.y = config->f.y + curr.y - next.y;

        /*C is empty or is equals do previous amino-acid*/
        if ((config->c.x == prev.x && config->c.y == prev.y) ||
            lattice[config->c.x][config->c.y] == -1)
        {
                config->next = next;
                config->curr = curr;
                config->prev = prev;
                config->pm_type = pm_type;

                return 1;
        }

    }

    /* Impossible generate the pull-move configuration */
    return 0;
}

void generate_pm_configs
(
    int amino_acid_index,
    int **lattice,
    struct coord *ant_positions,
    struct pm_config *configs,
    int *config_index,
    int seq_len
)
/* ============================================================
 * Generates pull-move configurations for a specific amino-acid
 * ============================================================
 */
{
    int result;
    struct pm_config config;
    config.amino_acid_index = amino_acid_index;

    struct coord right_move = create_new_coord(1, 0);
    struct coord left_move = create_new_coord(-1, 0);
    struct coord up_move = create_new_coord(0, 1);
    struct coord down_move = create_new_coord(0, -1);

    struct coord next = ant_positions[amino_acid_index + 1];
    struct coord prev = ant_positions[amino_acid_index - 1];
    struct coord curr = ant_positions[amino_acid_index];

    result = generate_pm_config(&config, lattice, ORIGINAL, curr, prev, next, right_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, ORIGINAL, curr, prev, next, left_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, ORIGINAL, curr, prev, next, up_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, ORIGINAL, curr, prev, next, down_move);
    if (result) configs[(*config_index)++] = config;

    next = ant_positions[amino_acid_index - 1];
    prev = ant_positions[amino_acid_index + 1];

    result = generate_pm_config(&config, lattice, INVERSE, curr, prev, next, right_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, INVERSE, curr, prev, next, left_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, INVERSE, curr, prev, next, up_move);
    if (result) configs[(*config_index)++] = config;
    result = generate_pm_config(&config, lattice, INVERSE, curr, prev, next, down_move);
    if (result) configs[(*config_index)++] = config;

}

void pm_search
(
    enum polarity *seq,
    int seq_len,
    struct ant *original_ant,
    int **lattice,
    struct ant best_ant,
    struct ant ant,
    struct pm_config *configs
)
/* ============================================================
 * Applies pull-moves on the protein while there is improvement
 * ============================================================
 */
{
    int i;
    int j;
    int num_configs = 0;
    int previous_energy;
    struct coord lattice_adjust;

    best_ant.energy = 0;

    //APPLY PULL-MOVE WHILE THERE IS IMPROVEMENT

    do
    {
        num_configs = 0;
        previous_energy = original_ant->energy;

        for (i = 0; i < seq_len; ++i)
        {
            lattice[original_ant->positions[i].x][original_ant->positions[i].y] = i;
        }

        //GENERATE ALL POSSIBLE PULL-MOVES CONFIGURATIONS

        for (i = 1; i < seq_len - 1; ++i)
        {
            generate_pm_configs(i, lattice, original_ant->positions,
                                configs, &num_configs, seq_len);
        }

        /*Apply all possible pull-moves and save the best*/
        for (i = 0; i < num_configs; ++i)
        {
            /*Copy original original_ant*/
            ant.energy = original_ant->energy;
            for (j = 0; j < seq_len; ++j)
            {
                ant.positions[j] = original_ant->positions[j];
            }

            /*Apply pull move on the copy*/
            ant.energy = apply_pm(ant, configs[i], seq, lattice,
                                         seq_len, original_ant->positions);

            if (ant.energy < best_ant.energy)
            {
                best_ant.energy = ant.energy;

                for (j = 0; j < seq_len; ++j)
                {
                    best_ant.positions[j] = ant.positions[j];
                }
            }
        }

        //CHECKS IF NEW CONFORMATION IS BETTER THAN ORIGINAL CONFORMATION

        if (best_ant.energy != 0 && best_ant.energy <= original_ant->energy)
        {
            lattice_adjust = create_new_coord(seq_len - best_ant.positions[0].x,
                                              seq_len - best_ant.positions[0].y);

            original_ant->energy = previous_energy;

            for (i = 0; i < seq_len; ++i)
            {
                lattice[original_ant->positions[i].x][original_ant->positions[i].y] = -1;

                /*Move protein to center*/
                best_ant.positions[i].x += lattice_adjust.x;
                best_ant.positions[i].y += lattice_adjust.y;

                original_ant->energy = best_ant.energy;

                /*Update original original_ant*/
                original_ant->positions[i] = best_ant.positions[i];

            }
        }
    }
    while(best_ant.energy < previous_energy);

    for (i = 0; i < seq_len; ++i)
    {
        lattice[original_ant->positions[i].x][original_ant->positions[i].y] = -1;
    }

}

char* ant_to_string
(
    struct ant ant,
    int num_dimensions
)
/* ===============================================
 * Extracts conformation info from struct ant to struct solution
 * ===============================================
 */
{
    struct coord move, prev_move;
    int i, direction;
    char *directions;

    directions = (char*) malloc(sizeof(char) * num_dimensions);

    for (i = 0; i < num_dimensions - 1; ++i)
    {
        move = subtract_coord(ant.positions[i + 1], ant.positions[i]);

        if (i == 0)
        {
            direction = abs_direction_by_move(move);
        }
        else
        {
            direction = direction_by_move(prev_move, move);
        }

        switch(direction)
        {
        case LEFT:
            directions[i] = 'L';
            break;
        case RIGHT:
            directions[i] = 'R';
            break;
        case STRAIGHT:
            directions[i] = 'S';
            break;
        default:
            break;
        }

        prev_move = move;
    }

    directions[num_dimensions - 1] = '\0';
    return directions;
}

int execute_look_ahead(int **lattice, struct coord position, int seq_len, int amin_index)
{
    int right;
    int left;
    int down;
    int up;

    int position_is_invalid;

    right = lattice[position.x + 1][position.y];
    left = lattice[position.x - 1][position.y];
    up = lattice[position.x][position.y + 1];
    down = lattice[position.x][position.y - 1];

    if (right > -1 && left > -1 && up > -1 && down > -1) {
        position_is_invalid = 1;
    } else {
        position_is_invalid = 0;
    }

    if (amin_index == 0 || amin_index == seq_len - 1)
    {
        position_is_invalid = 0;
    }

    return !position_is_invalid;
}

void construct_0
(
    struct aco_config aco_config,
    double **pheromone,
    int **lattice,
    enum polarity *seq,
    int seq_len,
    struct ant *ant,
    int *best_ant_by_index,
    int ant_index,
    struct ant *ants
)
/* ====================================================
 * Builds ant conformation. Based in Xiao, Li & Hu 2014
 * ====================================================
 */
{
    int i;
    int j;

    int num_candidates;
    int selected_candidate;
    double sum_probabilities;
    double probabilities[3];
    struct coord candidate_moves[3];
    struct candidate candidates[3];

    struct coord curr_position;
    struct coord foward_position;
    struct coord last_move;

    struct ant copied_ant;

    ant->energy = 0;
    ant->energy_by_edge[0] = 0;

    //Defines first edge direction

    lattice[seq_len][seq_len] = 0;
    ant->positions[0] = create_new_coord(seq_len, seq_len);

    last_move = create_new_coord(0, 1);
    curr_position = create_new_coord(ant->positions[0].x + last_move.x, ant->positions[0].y + last_move.y);
    lattice[curr_position.x][curr_position.y] = 1;
    ant->positions[1] = curr_position;

    //Constructor loop

    //For each edge except the first
    for (i = 1; i < seq_len - 1; ++i)
    {
        sum_probabilities = 0;
        num_candidates = 0;

        candidate_moves[0] = left(last_move);
        candidate_moves[1] = right(last_move);
        candidate_moves[2] = straight(last_move);

        //Defines candidate directions

        //For each direction
        for (j = 0; j < 3; ++j)
        {

            foward_position = create_new_coord(curr_position.x + candidate_moves[j].x,
                                               curr_position.y + candidate_moves[j].y);

            //If the next position in this direction is not occupied, turns current direction into a candidate
            if (lattice[foward_position.x][foward_position.y] == -1)
            {
                candidates[num_candidates].move = candidate_moves[j];
                candidates[num_candidates].position = foward_position;

                if (seq[i + 1] == H)
                {
                    candidates[num_candidates].heuristic =
                        calculate_heuristic(lattice, i + 1, foward_position, seq);
                }
                else
                {
                    candidates[num_candidates].heuristic = 0;
                }

                probabilities[num_candidates] =
                    pow(pheromone[i][j], aco_config.alpha) *
                    pow(exp((double) candidates[num_candidates].heuristic / 0.3), aco_config.beta);

                sum_probabilities += probabilities[num_candidates];
                ++num_candidates;
            }

        }

        for (j = 0; j < num_candidates; ++j)
        {
            probabilities[j] = probabilities[j]/sum_probabilities;
        }

        //Selects a candidate

        if (num_candidates == 0)
        {
            selected_candidate = -1;
        }
        else if (sum_probabilities == 0)
        {
            selected_candidate = rand() % num_candidates;
        }
        else if (num_candidates > 1)
        {
            selected_candidate = random_select(probabilities, num_candidates);
        }
        else if (num_candidates == 1)
        {
            selected_candidate = 0;
        }

        //Update conformation

        if (selected_candidate != -1)
        {
            last_move = candidates[selected_candidate].move;
            curr_position = candidates[selected_candidate].position;
            ant->energy_by_edge[i] = ant->energy;
            ant->energy -= candidates[selected_candidate].heuristic;
            ant->positions[i + 1] = curr_position;
            lattice[curr_position.x][curr_position.y] = i + 1;

        }
        else
        {
            //When is impossible continue the fold process

            //If theres no another ant to copy
            if (best_ant_by_index[i] == -1)
            {
                for (j = 0; j <= i; ++j)
                {
                    lattice[ant->positions[j].x][ant->positions[j].y] = -1;
                }

                i = 0;
                ant->energy = 0;

                //Reset first edge
                lattice[seq_len][seq_len] = 0;
                ant->positions[0] = create_new_coord(seq_len, seq_len);

                last_move = create_new_coord(0, 1);
                curr_position = create_new_coord(ant->positions[0].x + last_move.x,
                                                 ant->positions[0].y + last_move.y);
                lattice[curr_position.x][curr_position.y] = 1;
                ant->positions[1] = curr_position;
            }
            else
            {
                copied_ant = ants[best_ant_by_index[i]];

                for (j = 0; j <= i; ++j)
                {
                    lattice[ant->positions[j].x][ant->positions[j].y] = -1;
                }

                //Copy best ant for index i until i th amino-acid*/
                for (j = 0; j <= i; ++j)
                {
                    lattice[copied_ant.positions[j].x][copied_ant.positions[j].y] = j;
                    ant->positions[j] = copied_ant.positions[j];
                    if (j != i)
                    {
                        ant->energy_by_edge[j] = copied_ant.energy_by_edge[j];
                    }
                }
                curr_position = copied_ant.positions[i];
                ant->energy = copied_ant.energy_by_edge[i];

                last_move = create_new_coord(curr_position.x - copied_ant.positions[i - 1].x,
                                        curr_position.y - copied_ant.positions[i - 1].y);
                --i;
            }
        }
    }

    for (i = 0; i < seq_len - 1; ++i)
    {
        if (best_ant_by_index[i] == -1 ||
            ants[best_ant_by_index[i]].energy_by_edge[i] > ant->energy_by_edge[i])
        {
            best_ant_by_index[i] = ant_index;
        }
    }
}

void construct_1
(
    struct aco_config aco_config,
    double **pheromone,
    int **lattice,
    enum polarity *seq,
    int seq_len,
    struct ant *ant,
    int ant_index,
    struct ant *ants
)
/* ======================================================
 * Builds ant conformation. Based in  Hu, Zhang & Li 2009
 * ======================================================
 */
{
    int i;

    int num_candidates;
    int selected_candidate;
    double sum_probabilities;
    double probabilities[3];
    struct coord candidate_moves[3];
    struct candidate candidates[3];

    struct coord foward_position;
    struct coord last_move;

    int curr_amin;
    int prev_amin;
    int next_amin;

    int left_extremity;
    int right_extremity;
    int *extremity;

    int side;
    int last_unfold_side;
    int unfold_size;

    ant->energy = 0;

    //Defines start edge direction

    left_extremity = ceil((double) seq_len / 2);
    right_extremity = ceil((double) seq_len / 2) + 1;

    ant->positions[left_extremity] = create_new_coord(seq_len, seq_len);
    ant->positions[right_extremity] = create_new_coord(seq_len + 1, seq_len);

    lattice[ant->positions[left_extremity].x][ant->positions[left_extremity].y] = left_extremity;
    lattice[ant->positions[right_extremity].x][ant->positions[right_extremity].y] = right_extremity;

    side = 0;//left
    last_unfold_side = -1;


    //Constructor loop

    while (left_extremity > 0 || right_extremity < seq_len - 1)
    {
        if (side == 0 && left_extremity == 0)
        {
            side = 1;
        }
        else if (side == 1 && right_extremity == seq_len - 1)
        {
            side = 0;
        }

        if (side == 0)
        {
            curr_amin = left_extremity;
            next_amin = left_extremity - 1;
            prev_amin = left_extremity + 1;
            extremity = &left_extremity;
        }
        else
        {
            curr_amin = right_extremity;
            next_amin = right_extremity + 1;
            prev_amin = right_extremity - 1;
            extremity = &right_extremity;
        }

        last_move = subtract_coord(ant->positions[curr_amin], ant->positions[prev_amin]);

        sum_probabilities = 0;
        num_candidates = 0;

        candidate_moves[0] = left(last_move);
        candidate_moves[1] = right(last_move);
        candidate_moves[2] = straight(last_move);

        //Defines candidate directions

        //For each direction
        for (i = 0; i < 3; ++i)
        {

            foward_position = create_new_coord(ant->positions[curr_amin].x + candidate_moves[i].x,
                                               ant->positions[curr_amin].y + candidate_moves[i].y);

            //If the next position in this direction is not occupied, turns current direction into a candidate
            if (lattice[foward_position.x][foward_position.y] == -1 &&
                execute_look_ahead(lattice, foward_position, seq_len, next_amin))
            {
                candidates[num_candidates].move = candidate_moves[i];
                candidates[num_candidates].position = foward_position;

                if (seq[next_amin] == H)
                {
                    candidates[num_candidates].heuristic =
                        calculate_heuristic(lattice, next_amin, foward_position, seq);
                }
                else
                {
                    candidates[num_candidates].heuristic = 0;
                }

                probabilities[num_candidates] =
                    pow(pheromone[curr_amin][i], aco_config.alpha) *
                    pow(exp((double) candidates[num_candidates].heuristic / 0.3), aco_config.beta);

                sum_probabilities += probabilities[num_candidates];
                ++num_candidates;
            }

        }

        for (i = 0; i < num_candidates; ++i)
        {
            probabilities[i] = probabilities[i] / sum_probabilities;
        }

        //Selects a candidate

        if (num_candidates == 0)
        {
            selected_candidate = -1;
        }
        else if (sum_probabilities == 0)
        {
            selected_candidate = rand() % num_candidates;
        }
        else if (num_candidates > 1)
        {
            selected_candidate = random_select(probabilities, num_candidates);
        }
        else if (num_candidates == 1)
        {
            selected_candidate = 0;
        }

        //Update conformation

        if (selected_candidate != -1)
        {
            *extremity = next_amin;

            ant->energy -= candidates[selected_candidate].heuristic;
            ant->positions[*extremity] = candidates[selected_candidate].position;
            lattice[ant->positions[*extremity].x][ant->positions[*extremity].y] = *extremity;
            side = !side;

        }
        else
        {
            //When is impossible continue the fold process

            if (right_extremity == ceil((double) seq_len / 2) + 1 ||
                (last_unfold_side == 1 && left_extremity != ceil((double) seq_len / 2)))

            {
                unfold_size = 1 + rand() % ((int) ceil((double) seq_len / 2) - left_extremity);
                extremity = &left_extremity;
                last_unfold_side = 0;

                curr_amin = *extremity;
                next_amin = curr_amin - 1;
                prev_amin = curr_amin + 1;
            }
            else
            {

                unfold_size = 1 + rand() % (right_extremity - ((int) ceil((double) seq_len / 2) + 1));
                extremity = &right_extremity;
                last_unfold_side = 1;

                curr_amin = *extremity;
                next_amin = curr_amin + 1;
                prev_amin = curr_amin - 1;
            }


            for (i = 0; i < unfold_size; ++i)
            {
                if (seq[*extremity] == H)
                {
                    ant->energy += calculate_heuristic(lattice, *extremity, ant->positions[*extremity], seq);
                }
                lattice[ant->positions[*extremity].x][ant->positions[*extremity].y] = -1;

                next_amin = *extremity;
                curr_amin = prev_amin;
                prev_amin = curr_amin - (next_amin - curr_amin);

                *extremity = curr_amin;
            }
        }
    }
}

void construct_2
(
    struct aco_config aco_config,
    double **pheromone,
    int **lattice,
    enum polarity *seq,
    int seq_len,
    struct ant *ant,
    int ant_index,
    struct ant *ants
)
/* ========================================================
 * Builds ant conformation. Based in Shmygelska & Hoos 2003
 * ========================================================
 */
{
    int i;

    int num_candidates;
    int selected_candidate;
    double sum_probabilities;
    double probabilities[3];
    struct coord candidate_moves[3];
    struct candidate candidates[3];

    struct coord last_move;
    struct coord foward_position;
    enum direction forbidden_dir;
    double left_unfolded_percentage;

    int curr_amin;
    int prev_amin;
    int next_amin;
    int left_extremity;
    int right_extremity;
    int *extremity;
    int side;
    int unfold_size;
    int start_point;

    double readjust;

    forbidden_dir = -1;
    ant->energy = 0;

    //Defines start edge direction

    start_point = rand()%(seq_len - 1);

    left_extremity = start_point;
    right_extremity = start_point + 1;

    ant->positions[left_extremity] = create_new_coord(seq_len, seq_len);
    ant->positions[right_extremity] = create_new_coord(seq_len + 1, seq_len);

    lattice[ant->positions[left_extremity].x][ant->positions[left_extremity].y] = left_extremity;
    lattice[ant->positions[right_extremity].x][ant->positions[right_extremity].y] = right_extremity;

    //Constructor loop

    while (left_extremity > 0 || right_extremity < seq_len - 1)
    {

        left_unfolded_percentage = (double) left_extremity / (left_extremity + (seq_len - right_extremity));

        if ((double) rand() / RAND_MAX < left_unfolded_percentage) {
            side = 0;//left
        } else {
            side = 1;//right
        }


        if (side == 0 && left_extremity == 0)
        {
            side = 1;
        }
        else if (side == 1 && right_extremity == seq_len - 1)
        {
            side = 0;
        }

        if (side == 0)
        {
            curr_amin = left_extremity;
            next_amin = left_extremity - 1;
            prev_amin = left_extremity + 1;
            extremity = &left_extremity;
        }
        else
        {
            curr_amin = right_extremity;
            next_amin = right_extremity + 1;
            prev_amin = right_extremity - 1;
            extremity = &right_extremity;
        }

        last_move = subtract_coord(ant->positions[curr_amin], ant->positions[prev_amin]);

        sum_probabilities = 0;
        num_candidates = 0;

        candidate_moves[0] = left(last_move);
        candidate_moves[1] = right(last_move);
        candidate_moves[2] = straight(last_move);

        //Defines candidate direction

        //For each direction
        for (i = 0; i < 3; ++i)
        {

            foward_position = create_new_coord(ant->positions[curr_amin].x + candidate_moves[i].x,
                                               ant->positions[curr_amin].y + candidate_moves[i].y);

            //If the next position in this direction is not occupied, turns current direction into a candidate
            if (lattice[foward_position.x][foward_position.y] == -1 &&
                execute_look_ahead(lattice, foward_position, seq_len, next_amin) &&
                i != forbidden_dir)
            {
                candidates[num_candidates].move = candidate_moves[i];
                candidates[num_candidates].position = foward_position;
                candidates[num_candidates].readjusted = 0;

                if (seq[next_amin] == H)
                {
                    candidates[num_candidates].heuristic =
                        calculate_heuristic(lattice, next_amin, foward_position, seq);
                }
                else
                {
                    candidates[num_candidates].heuristic = 0;
                }

                probabilities[num_candidates] =
                    pow(pheromone[curr_amin][i], aco_config.alpha) *
                    pow(exp((double) candidates[num_candidates].heuristic / 0.3), aco_config.beta);

                sum_probabilities += probabilities[num_candidates];
                ++num_candidates;
            }

        }

        forbidden_dir = -1;

        readjust = 0;

        for (i = 0; i < num_candidates; ++i)
        {
            probabilities[i] = probabilities[i] / sum_probabilities;

            if (probabilities[i] < aco_config.min_probability)
            {
                candidates[i].readjusted = 1;
                readjust += aco_config.min_probability - probabilities[i];
                probabilities[i] = aco_config.min_probability;
            }
        }

        if (readjust > 0)
        {
            for (i = 0; i < num_candidates; ++i)
            {
                if (candidates[i].readjusted == 0)
                {
                    probabilities[i] -= (readjust * probabilities[i]);
                }
            }
        }

        //Selects a candidate

        if (num_candidates == 0)
        {
            selected_candidate = -1;
        }
        else if (sum_probabilities == 0)
        {
            selected_candidate = rand() % num_candidates;
        }
        else if (num_candidates > 1)
        {
            selected_candidate = random_select(probabilities, num_candidates);
        }
        else if (num_candidates == 1)
        {
            selected_candidate = 0;
        }

        //Update conformation

        if (selected_candidate != -1)
        {
            *extremity = next_amin;

            ant->energy -= candidates[selected_candidate].heuristic;
            ant->positions[*extremity] = candidates[selected_candidate].position;
            lattice[ant->positions[*extremity].x][ant->positions[*extremity].y] = *extremity;
            side = !side;

        }
        else
        {
            unfold_size = ceil((double)(right_extremity - left_extremity) / 2);

            if (side == 0) {


                extremity = &left_extremity;
                curr_amin = *extremity;
                next_amin = curr_amin - 1;
                prev_amin = curr_amin + 1;

            } else {

                extremity = &right_extremity;
                curr_amin = *extremity;
                next_amin = curr_amin + 1;
                prev_amin = curr_amin - 1;

            }

            for (i = 0; i < unfold_size; ++i)
            {
                if (seq[*extremity] == H)
                {
                    ant->energy += calculate_heuristic(lattice, *extremity, ant->positions[*extremity], seq);
                }
                lattice[ant->positions[*extremity].x][ant->positions[*extremity].y] = -1;

                next_amin = *extremity;
                curr_amin = prev_amin;
                prev_amin = curr_amin - (next_amin - curr_amin);

                *extremity = curr_amin;
            }

            forbidden_dir = direction_by_move(
                subtract_coord(ant->positions[curr_amin], ant->positions[prev_amin]),
                subtract_coord(ant->positions[next_amin], ant->positions[curr_amin]));
        }
    }
}

void partition(struct ant* ants, int left, int right, int *i, int *j)
{
    struct ant aux;
    double pivot;
    *i = left;
    *j = right;

    pivot = ants[(*i + *j) / 2].energy;

    do
    {
        while(ants[*i].energy < pivot)
        {
            (*i)++;
        }
        while(ants[*j].energy > pivot)
        {
            (*j)--;
        }
        if(*i <= *j)
        {
            aux = ants[*i];
            ants[*i] = ants[*j];
            ants[*j] = aux;

            (*i)++;
            (*j)--;
        }
    }
    while(*i <= *j);
}

void sort(struct ant* ants, int left, int right)
{
    int i, j;
    partition(ants, left, right, &i, &j);
    if(left < j) sort(ants, left, j);
    if(i < right) sort(ants, i, right);
}

void quick_sort(struct ant* ants, int n)
{
    sort(ants, 0, n - 1);
}

struct aco_result aco_run
(
    enum polarity *seq,
    int seq_len,
    struct aco_config aco_config,
    int *seed
)
/* ====================================
 * ACO main function
 * ====================================
 */
{
    int i;
    int j;
    int k;
    int *best_ant_by_edge;
    int **lattice;
    struct ant iteration_ant;
    struct ant best_ant;
    struct ant *ants;
    struct ant pm_best_ant;
    struct ant pm_ant;
    struct pm_config* pm_configs;
    double** pheromone;
    clock_t t0;
    struct aco_result result;

    /**///Code related to output statistics
    /**/result.energy_evolution = (int*) malloc(sizeof(int) * aco_config.iterations);
    /**/t0 = clock();

    //Sets seed
    if (*seed == -1)
    {
        *seed = (unsigned) time(NULL);
    }
    srand(*seed);

    init_variables(aco_config, &pm_configs, &pm_best_ant, &pm_ant, &best_ant,
                   seq_len, &best_ant_by_edge, &lattice, &pheromone, &ants);

    best_ant.energy = 0;

    //Iteration loop
    for (i = 0; i < aco_config.iterations; ++i)
    {
        //Population loop
        for (j = 0; j < aco_config.population; ++j)
        {

            switch (aco_config.constructor) {

            case XIAO_LI_HU_2014:
                construct_0(aco_config, pheromone, lattice, seq, seq_len, &ants[j], best_ant_by_edge, j, ants);
                break;
            case HU_ZHANG_LI_2009:
                construct_1(aco_config, pheromone, lattice, seq, seq_len, &ants[j], j, ants);
                break;
            case SHMYGELSKA_HOOS_2003:
                construct_2(aco_config, pheromone, lattice, seq, seq_len, &ants[j], j, ants);
                break;
            default:
                printf("ERROR: aco.c/aco_run(): \"Invalid constructor\"\n");
                exit(1);
                break;

            }

            //Clean lattice
            for (k = 0; k < seq_len; ++k)
            {
                lattice[ants[j].positions[k].x][ants[j].positions[k].y] = -1;
            }
        }

        //reset best ants by edge
        for (j = 0; j < seq_len - 1; ++j)
        {
            best_ant_by_edge[j] = -1;
        }

        iteration_ant = ants[0];

        //enum daemon search
        for (j = 0; j < aco_config.population; ++j)
        {

            switch (aco_config.daemon)
            {
            case PULL_MOVE:
                pm_search(seq, seq_len, &ants[j], lattice, pm_best_ant, pm_ant, pm_configs);
                break;
            default:
                break;

            }

            if (ants[j].energy < iteration_ant.energy)
            {
                iteration_ant = ants[j];
            }

        }

        /**/if (i == 0 || iteration_ant.energy < best_ant.energy) {
        /**/    result.energy_evolution[i] = iteration_ant.energy;
        /**/} else {
        /**/    result.energy_evolution[i] = 1;
        /**/}

        //Update best ant
        if (iteration_ant.energy < best_ant.energy)
        {
            /**/result.found_on_iteration = i;
            best_ant.energy = iteration_ant.energy;
            for (j = 0; j < seq_len; ++j)
            {
                best_ant.positions[j] = iteration_ant.positions[j];
            }
        }

        if (aco_config.elit_percentage < 1) {
            quick_sort(ants, aco_config.population);
        }

        //Pheromone update
        pheromone_evaporation(pheromone, seq_len, aco_config.persistence);
        for (j = 0; j < aco_config.population * aco_config.elit_percentage; ++j)
        {
            pheromone_deposit(pheromone, ants[j], seq_len, best_ant.energy);
        }
    }

    /**/result.time = (clock() - t0)/(double)CLOCKS_PER_SEC;
    /**/result.energy = best_ant.energy;
    /**/result.final_population_avg = 0;
    /**/result.final_population_solution_rate = 0;
    /**/result.final_population_stddev = 0;

    /**/for (i = 0; i < aco_config.population; ++i) {
    /**/    result.final_population_avg += ants[i].energy;
    /**/    if (ants[i].energy == best_ant.energy)
    /**/    {
    /**/        ++result.final_population_solution_rate;
    /**/    }
    /**/}

    /**/result.final_population_solution_rate /= aco_config.population;
    /**/result.final_population_avg /= aco_config.population;

    /**/for (i = 0; i < aco_config.population; ++i) {
    /**/    result.final_population_stddev +=
    /**/        (result.final_population_avg - ants[i].energy) *
    /**/        (result.final_population_avg - ants[i].energy);
    /**/}

    /**/result.final_population_stddev /= aco_config.population;
    /**/result.final_population_stddev = sqrt(result.final_population_stddev);
    /**/result.final_population_solution_rate *= 100;
    /**/result.directions = ant_to_string(best_ant, seq_len);

    free_variables(aco_config, seq_len, lattice, best_ant_by_edge, best_ant,
                   pm_best_ant, pm_ant, pheromone, pm_configs, ants);

    return result;
}
