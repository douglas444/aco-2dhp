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

enum polarity
{
    P = 0,
    H = 1
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

};

typedef enum direction Direction;
typedef enum polarity Polarity;
typedef enum pm_type Pm_type;
typedef struct coord Coord;
typedef struct ant Ant;
typedef struct pm_config Pm_config;
typedef struct candidate Candidate;
typedef struct lattice Lattice;



void* smalloc(int mem_size)
/* ====================================
 * Allocates memory
 * ====================================
 */
{
    void *mem_pos = (void*) malloc(mem_size);

    if (mem_pos == NULL)
    {
        printf("Error in function smalloc: Unable to allocate memory\n");
        exit(1);
    }
    else
    {
        return mem_pos;
    }
}



void init_ant(Ant *ant, int seq_len)

/* ====================================
 * Initializes ant
 * ====================================
 */
{
    ant->energy_by_edge = (int*) smalloc(sizeof(int) * (seq_len - 1));
    ant->positions = (Coord*) smalloc(sizeof(Coord) * seq_len);
}



void free_ant(Ant ant)
/* ====================================
 * Releases ant
 * ====================================
 */
{
    free(ant.energy_by_edge);
    free(ant.positions);
}



void variables_initialization
(
    ACO_config aco_config,
    Pm_config **pm_configs,
    Ant *pm_best_ant,
    Ant *pm_ant,
    Ant **ants,
    Ant *best_ant,
    int seq_len,
    int **best_ant_by_edge,
    int ***lattice,
    double ***pheromone,
    Solution *solution
)
/* =========================================
 * Initializes all aco.c exclusive variables
 * =========================================
 */
{
    int i;
    int j;

    *pm_configs = (Pm_config*) smalloc(sizeof(Pm_config) * 4 * (seq_len - 2));

    init_ant(pm_best_ant, seq_len);
    init_ant(pm_ant, seq_len);
    init_ant(best_ant, seq_len);

    solution->directions = (char*) smalloc(sizeof(char) * seq_len);

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

    *ants = (Ant*) smalloc(sizeof(Ant) * aco_config.population);
    for (i = 0; i < aco_config.population; ++i)
    {
        init_ant(&((*ants)[i]), seq_len);
    }
}



void free_variables
(
    ACO_config aco_config,
    int seq_len,
    int **lattice,
    int *best_ant_by_edge,
    Ant pm_best_ant,
    Ant pm_ant,
    Ant *ants,
    Ant best_ant,
    double **pheromone,
    Pm_config *pm_configs
)
/* ====================================
 * Releases all aco.c exclusive variables
 * ====================================
 */
{
    int i;

    for (i = 0; i < seq_len; ++i)
    {
        free(pheromone[i]);
    }
    for (i = 0; i < aco_config.population; ++i)
    {
        free_ant(ants[i]);
    }
    for (i = 0; i < 2 * seq_len + 1; ++i)
    {
        free(lattice[i]);
    }

    free(best_ant_by_edge);
    free(lattice);
    free(pheromone);
    free(ants);
    free_ant(pm_best_ant);
    free_ant(pm_ant);
    free_ant(best_ant);
    free(pm_configs);
}



Coord subtract_coord
(
    Coord c1,
    Coord c2
)
/* =======================================================
 * Calculates the difference between two given coordinates
 * =======================================================
 */
{
    Coord c3;

    c3.x = c1.x - c2.x;
    c3.y = c1.y - c2.y;

    return c3;
}



int lateral_adj(Coord c1, Coord c2)
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



int abs_direction_by_move(Coord move)
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
            printf("Error in function abs_direction_by_move: Invalid value for parameter move\n");
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
            printf("Error in function abs_direction_by_move: Invalid value for parameter move\n");
            exit(1);
        }
    }
    else
    {
        printf("Error in function abs_direction_by_move: Invalid value for parameter move\n");
        exit(1);
    }
}



int direction_by_move(Coord prev_move, Coord move)
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
            printf("Error in function direction_by_move: Invalid values for parameters\n");
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
            printf("Error in function direction_by_move: Invalid values for parameters\n");
            exit(1);
        }
    }
}



int calculate_heuristic
(
    int **lattice,
    int amino_acid_index,
    Coord pos,
    int *seq
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



Coord create_new_coord(int x, int y)
/* ===========================================================
 * Returns a new variables of type Coord with the given values
 * ===========================================================
 */
{
    Coord coord;
    coord.x = x;
    coord.y = y;
    return coord;
}



void pheromone_deposit
(
    double **pheromone,
    Ant ant,
    int seq_len,
    int best_energy
)
/* ==============================
 * ACO pheromone deposit function
 * ==============================
 */
{
    int i;
    Coord move, prev_move;
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
            pheromone[i][direction] += (double) ant.energy / best_energy * best_energy * best_energy;
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
 * ACO pheromone evaporation function
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



Coord straight(Coord prev_move)
/* =================================================================
 * Calculates straight move corresponding to the given previous move
 * =================================================================
 */
{
    Coord move;

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
        printf("Error in function straight: Invalid value for parameter prev_move\n");
        exit(1);
    }

    return move;
}



Coord left(Coord prev_move)
/* =============================================================
 * Calculates left move corresponding to the given previous move
 * =============================================================
 */
{
    Coord move;

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
        printf("Error in function left: Invalid value for parameter prev_move\n");
        exit(1);
    }

    return move;
}



Coord right(Coord prev_move)
/* ==============================================================
 * Calculates right move corresponding to the given previous move
 * ==============================================================
 */
{
    Coord move;

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
        printf("Error in function right: Invalid value for parameter prev_move\n");
        exit(1);
    }

    return move;
}



int generate_pm_config
(
    Pm_config *config,
    int **lattice,
    Pm_type pm_type,
    Coord curr,
    Coord prev,
    Coord next,
    Coord direction
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
    Coord *ant_positions,
    Pm_config *configs,
    int *config_index,
    int seq_len
)
/* ============================================================
 * Generates pull-move configurations for a specific amino-acid
 * ============================================================
 */
{
    int result;
    Pm_config config;
    config.amino_acid_index = amino_acid_index;

    Coord right_move = create_new_coord(1, 0);
    Coord left_move = create_new_coord(-1, 0);
    Coord up_move = create_new_coord(0, 1);
    Coord down_move = create_new_coord(0, -1);

    Coord next = ant_positions[amino_acid_index + 1];
    Coord prev = ant_positions[amino_acid_index - 1];
    Coord curr = ant_positions[amino_acid_index];

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



int move_amino_acid
(
    int current_energy,
    int **lattice,
    int *seq,
    int amino_acid_index,
    Coord src,
    Coord dest
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
    Ant ant,
    Pm_config config,
    int *seq,
    int **lattice,
    int seq_len,
    Coord *ant_positions
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
    Coord tempCoord;

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

            tempCoord = ant.positions[i];
            ant.positions[i] = config.curr;
            config.curr = config.prev;
            config.prev = tempCoord;

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



void pm_search
(
    int *seq,
    int seq_len,
    Ant *original_ant,
    int **lattice,
    Ant best_ant,
    Ant ant,
    Pm_config *configs
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
    Coord lattice_adjust;

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

    original_ant->energy = best_ant.energy;

}



void extract_solution
(
    Ant ant,
    Solution *solution,
    int seq_len
)
/* ===============================================
 * Extracts conformation info from Ant to Solution
 * ===============================================
 */
{
    Coord move, prev_move;
    int i, direction;

    solution->energy = ant.energy;

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

        switch(direction)
        {
        case LEFT:
            solution->directions[i] = 'L';
            break;
        case RIGHT:
            solution->directions[i] = 'R';
            break;
        case STRAIGHT:
            solution->directions[i] = 'S';
            break;
        default:
            break;
        }

        prev_move = move;
    }

    solution->directions[seq_len - 1] = '\0';
}



void construct_conform
(
    ACO_config aco_config,
    double **pheromone,
    int **lattice,
    int *seq,
    int seq_len,
    Ant *ant,
    int *best_ant_by_index,
    int ant_index,
    Ant *ants
)
/* ====================================
 * Builds ant conformation
 * ====================================
 */
{
    int i;
    int j;
    int num_candidates;
    int selected_candidate;
    double sum_probabilities;
    double probabilities[3];
    Coord curr_position;
    Coord move;
    Coord candidate_move[3];
    Candidate candidates[3];
    Ant copied_ant;

    ant->energy = 0;
    ant->energy_by_edge[0] = 0;

    //DEFINES FIRST EDGE

    lattice[seq_len][seq_len] = 0;
    ant->positions[0] = create_new_coord(seq_len, seq_len);

    move = create_new_coord(0, 1);
    curr_position = create_new_coord(ant->positions[0].x + move.x,
                                     ant->positions[0].y + move.y);
    lattice[curr_position.x][curr_position.y] = 1;
    ant->positions[1] = curr_position;


    //CONSTRUCTOR LOOP

    //For each edge except the first
    for (i = 1; i < seq_len - 1; ++i)
    {
        sum_probabilities = 0;
        num_candidates = 0;

        candidate_move[0] = left(move);
        candidate_move[1] = right(move);
        candidate_move[2] = straight(move);

        //DEFINES CANDIDATES DIRECTION

        //For each direction
        for (j = 0; j < 3; ++j)
        {

            //If the next position in this direction is not occupied, turns current direction into a candidate
            if (lattice[curr_position.x + candidate_move[j].x][curr_position.y + candidate_move[j].y] == -1)
            {
                candidates[num_candidates].move = candidate_move[j];
                candidates[num_candidates].position.x = curr_position.x + candidates[num_candidates].move.x;
                candidates[num_candidates].position.y = curr_position.y + candidates[num_candidates].move.y;

                if (seq[i + 1] == H)
                {
                    candidates[num_candidates].heuristic =
                        calculate_heuristic(lattice, i + 1, candidates[num_candidates].position, seq);
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

        //SELECTS A CANDIDATE

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

        //UPDATE CONFORMATION

        if (selected_candidate != -1)
        {
            move = candidates[selected_candidate].move;
            curr_position = candidates[selected_candidate].position;
            ant->energy_by_edge[i] = ant->energy;
            ant->energy -= candidates[selected_candidate].heuristic;
            ant->positions[i + 1] = curr_position;
            lattice[curr_position.x][curr_position.y] = i + 1;

        }
        else
        {
            //WHEN IS IMPOSSIBLE CONTINUE THE FOLD PROCESS

            switch (aco_config.collision_handler)
            {

            case PARTIAL_COPY:

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

                    move = create_new_coord(0, 1);
                    curr_position = create_new_coord(ant->positions[0].x + move.x,
                                                     ant->positions[0].y + move.y);
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

                    move = create_new_coord(curr_position.x - copied_ant.positions[i - 1].x,
                                            curr_position.y - copied_ant.positions[i - 1].y);
                    --i;
                }
                break;

            default:
                break;
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



Solution aco_run
(
    int *seq,
    int seq_len,
    ACO_config aco_config,
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
    Ant *ants;
    Ant iteration_ant;
    Ant best_ant;
    Ant pm_best_ant;
    Ant pm_ant;
    Pm_config* pm_configs;
    double** pheromone;
    Solution solution;


    //Sets seed
    if (*seed == -1)
    {
        *seed = (unsigned) time(NULL);
    }
    srand(*seed);

    variables_initialization(aco_config, &pm_configs, &pm_best_ant, &pm_ant,
                             &ants, &best_ant, seq_len, &best_ant_by_edge,
                             &lattice, &pheromone, &solution);

    best_ant.energy = 0;

    //Iteration loop
    for (i = 0; i < aco_config.iterations; ++i)
    {
        //Population loop
        for (j = 0; j < aco_config.population; ++j)
        {
            construct_conform(aco_config, pheromone, lattice, seq, seq_len, &ants[j],
                              best_ant_by_edge, j, ants);

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

        //Daemon search
        for (j = 0; j < aco_config.population; ++j)
        {

            switch (aco_config.daemon)
            {

            case PULL_MOVE:

                pm_search(seq, seq_len, &ants[j],
                                 lattice, pm_best_ant, pm_ant,
                                 pm_configs);
                break;

            default:
                break;

            }

            if (ants[j].energy < iteration_ant.energy)
            {
                iteration_ant = ants[j];
            }

        }

        //Update best ant
        if (iteration_ant.energy < best_ant.energy)
        {
            best_ant.energy = iteration_ant.energy;
            for (j = 0; j < seq_len; ++j)
            {
                best_ant.positions[j] = iteration_ant.positions[j];
            }
        }

        //Pheromone update
        pheromone_evaporation(pheromone, seq_len, aco_config.persistence);
        for (j = 0; j < aco_config.population; ++j)
        {
            pheromone_deposit(pheromone, ants[j], seq_len, best_ant.energy);
        }
    }

    extract_solution(best_ant, &solution, seq_len);
    return solution;
}
