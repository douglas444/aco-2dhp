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

struct lmatrix
{
    int n;
    int top;
    int order;
    int *values;
    int initial_value;
    unsigned int *to;
    unsigned int *from;
};


typedef enum direction Direction;
typedef enum polarity Polarity;
typedef enum pm_type PM_type;
typedef struct coord Coord;
typedef struct ant Ant;
typedef struct pm_config PM_config;
typedef struct candidate Candidate;
typedef struct lattice Lattice;
typedef struct lmatrix Lmatrix;



void initialize_lmatrix(Lmatrix *lmatrix, int initial_value, int n)
/* ====================================
 * Initialize lmatrix type
 * ====================================
 */
{
    int num_elements = n * n;

    lmatrix->top = 0;
    lmatrix->n = n;
    lmatrix->initial_value = initial_value;

    lmatrix->values = (int*) malloc(sizeof(int) * num_elements);
    lmatrix->to = (unsigned int*) malloc(sizeof(unsigned int) * num_elements);
    lmatrix->from = (unsigned int*) malloc(sizeof(unsigned int) * num_elements);

}


void destroy_lmatrix(Lmatrix *lmatrix)
/* ====================================
 * Free lmatrix memory
 * ====================================
 */
{
    free(lmatrix->values);
    free(lmatrix->to);
    free(lmatrix->from);
}



int lmatrix_read(Lmatrix *lmatrix, Coord coord)
/* ====================================
 * Read lmatrix
 * ====================================
 */
{
    int index = coord.x + lmatrix->n * coord.y;

    if (lmatrix->from[index] < lmatrix->top && lmatrix->to[lmatrix->from[index]] == index)
    {
        return lmatrix->values[index];
    }
    else
    {
        lmatrix->from[index] = lmatrix->top;
        lmatrix->to[lmatrix->top] = index;
        lmatrix->values[index] = lmatrix->initial_value;
        lmatrix->top++;
        return lmatrix->values[index];
    }
}



void lattice_write(Lmatrix *lmatrix, Coord coord, int value)
/* ====================================
 * Write lmatrix
 * ====================================
 */
{
    int index = coord.x + lmatrix->n * coord.y;

    if (lmatrix->from[index] < lmatrix->top && lmatrix->to[lmatrix->from[index]] == index)
    {
        lmatrix->values[index] = value;
    }
    else
    {
        lmatrix->from[index] = lmatrix->top;
        lmatrix->to[lmatrix->top] = index;
        lmatrix->values[index] = value;
        lmatrix->top++;
    }
}



void* memory_allocation(int mem_size)
/* ====================================
 * Allocates memory safely
 * ====================================
 */
{
    void *mem_pos = (void*) malloc(mem_size);

    if (mem_pos == NULL)
    {
        printf("Error in function memory_allocation: Unable to allocate memory\n");
        exit(1);
    }
    else
    {
        return mem_pos;
    }
}



void initialize_ant(Ant *ant, int sequence_len)

/* ====================================
 * Initializes ant
 * ====================================
 */
{
    ant->energy_by_edge = (int*) memory_allocation(sizeof(int) * (sequence_len - 1));
    ant->positions = (Coord*) memory_allocation(sizeof(Coord) * sequence_len);
}



void destroy_ant(Ant ant)
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
    PM_config **pm_configs,
    Ant *best_pm_ant,
    Ant *pm_ant,
    Ant **ants,
    Ant *best_ant,
    int sequence_len,
    int **best_ant_by_edge,
    Lmatrix *lattice,
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

    *pm_configs = (PM_config*) memory_allocation(sizeof(PM_config) * 4 * (sequence_len - 2));

    initialize_ant(best_pm_ant, sequence_len);
    initialize_ant(pm_ant, sequence_len);
    initialize_ant(best_ant, sequence_len);
    initialize_lmatrix(lattice, -1, 2 * sequence_len + 1);
    solution->directions = (char*) memory_allocation(sizeof(char) * (sequence_len - 1));

    *pheromone = (double**) memory_allocation(sizeof(double*) * sequence_len);
    for (i = 0; i < sequence_len; ++i)
    {
        (*pheromone)[i] = (double*) memory_allocation(sizeof(double) * 3);
        for (j = 0; j < 3; ++j)
        {
            (*pheromone)[i][j] = aco_config.ini_pheromone;
        }
    }

    *best_ant_by_edge = (int*) memory_allocation(sizeof(int) * (sequence_len - 1));
    for (i = 0; i < sequence_len - 1; ++i)
    {
        (*best_ant_by_edge)[i] = -1;
    }

    *ants = (Ant*) memory_allocation(sizeof(Ant) * aco_config.population);
    for (i = 0; i < aco_config.population; ++i)
    {
        initialize_ant(&((*ants)[i]), sequence_len);
    }
}



void free_variables
(
    ACO_config aco_config,
    int sequence_len,
    Lmatrix *lattice,
    int *best_ant_by_edge,
    Ant pm_ant,
    Ant pm_best_ant,
    Ant *ants,
    double **pheromone,
    PM_config *pm_configs
)
/* ====================================
 * Releases all aco.c exclusive variables
 * ====================================
 */
{
    int i;

    for (i = 0; i < sequence_len; ++i)
    {
        free(pheromone[i]);
    }
    for (i = 0; i < aco_config.population; ++i)
    {
        destroy_ant(ants[i]);
    }

    destroy_lmatrix(lattice);

    free(best_ant_by_edge);
    free(pheromone);
    free(ants);
    destroy_ant(pm_best_ant);
    destroy_ant(pm_ant);
    free(pm_configs);
}



Coord new_coord(int x, int y)
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



Coord subtract_coord(Coord c1, Coord c2)
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



int coords_distance(Coord c1, Coord c2)
/* =====================================================
 * Calculates the distance between two given coordinates
 * =====================================================
 */
{
    return sqrt((c1.x - c2.x) * (c1.x - c2.x) +
                (c1.y - c2.y) * (c1.y - c2.y));

}



int calculate_absolute_direction_by_move(Coord move)
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
            printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move\n");
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
            printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move\n");
            exit(1);
        }
    }
    else
    {
        printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move\n");
        exit(1);
    }
}



int calculate_direction_by_move(Coord prev_move, Coord move)
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
            printf("Error in function calculate_direction_by_move: Invalid values for parameters\n");
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
            printf("Error in function calculate_direction_by_move: Invalid values for parameters\n");
            exit(1);
        }
    }
}



int calculate_absolute_heuristic_value
(
    Lmatrix *lattice,
    int amino_acid_index,
    Coord pos,
    int *sequence
)
/* =================================================================================
 * Calculates the number of H-H contacts if a H amino-acid occupy the given position
 * =================================================================================
 */
{
    int heuristic_value = 0;
    int right_neighbor = lmatrix_read(lattice, new_coord(pos.x + 1, pos.y));
    int left_neighbor = lmatrix_read(lattice, new_coord(pos.x - 1, pos.y));
    int down_neighbor = lmatrix_read(lattice, new_coord(pos.x, pos.y - 1));
    int up_neighbor = lmatrix_read(lattice, new_coord(pos.x, pos.y + 1));

    if (right_neighbor >= 0 && abs(right_neighbor - amino_acid_index) > 1 &&
            sequence[right_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (left_neighbor >= 0 && abs(left_neighbor - amino_acid_index) > 1 &&
            sequence[left_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (up_neighbor >= 0 && abs(up_neighbor - amino_acid_index) > 1 &&
            sequence[up_neighbor] == H)
    {
        ++heuristic_value;
    }
    if (down_neighbor >= 0 && abs(down_neighbor - amino_acid_index) > 1 &&
            sequence[down_neighbor] == H)
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
    Ant ant,
    int sequence_len,
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

    for (i = 0; i < sequence_len - 1; ++i)
    {
        move = subtract_coord(ant.positions[i + 1], ant.positions[i]);
        if (i == 0)
        {
            direction = calculate_absolute_direction_by_move(move);
        }
        else
        {
            direction = calculate_direction_by_move(prev_move, move);
        }

        if (best_energy != 0)
        {
            pheromone[i][direction] += (double) ant.energy / pow(best_energy, 3);
        }
        prev_move = move;
    }
}



void pheromone_evaporation(double **pheromone, int sequence_len, double persistence)
/* ==================================
 * ACO pheromone evaporation function
 * ==================================
 */
{
    int i, j;

    for (i = 0; i < sequence_len - 1; ++i)
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
        printf("Error in function calculate_move_by_direction: Invalid value for parameter prev_move\n");
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
        printf("Error in function calculate_move_by_direction: Invalid value for parameter prev_move\n");
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
        printf("Error in function calculate_move_by_direction: Invalid value for parameter prev_move\n");
        exit(1);
    }

    return move;
}



int generate_pull_move_configs

(
    int amino_acid_index,
    Lmatrix *lattice,
    Coord *ant_positions,
    PM_config *configs,
    int config_index,
    int sequence_len
)
/* ====================================
 * Generates pull-move configurations
 * ====================================
 */
{
    PM_config config;
    int num_configs;
    int right_neighbor;
    int left_neighbor;
    int down_neighbor;
    int up_neighbor;

    num_configs = 0;

    config.pm_type = ORIGINAL;
    config.curr = ant_positions[amino_acid_index];
    config.next = ant_positions[amino_acid_index + 1];
    config.prev = ant_positions[amino_acid_index - 1];
    config.amino_acid_index = amino_acid_index;

    right_neighbor = lmatrix_read(lattice, new_coord(config.next.x + 1, config.next.y));
    left_neighbor = lmatrix_read(lattice, new_coord(config.next.x - 1, config.next.y));
    up_neighbor = lmatrix_read(lattice, new_coord(config.next.x, config.next.y + 1));
    down_neighbor = lmatrix_read(lattice, new_coord(config.next.x, config.next.y - 1));

    /*If is adjacent (right) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (right_neighbor == -1 && config.next.x + 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x + 1;
        config.f.y = config.next.y;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }

    }

    /*If is adjacent (left) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (left_neighbor == -1 && config.next.x - 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x - 1;
        config.f.y = config.next.y;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }
    }

    /*If is adjacent (up) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (up_neighbor == -1 && config.next.x != config.curr.x &&
            config.next.y + 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y + 1;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;

            ++config_index;
            ++num_configs;
        }
    }

    /*If is adjacent (down) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (down_neighbor == -1 && config.next.x != config.curr.x &&
            config.next.y - 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y - 1;

        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }
    }

    config.pm_type = INVERSE;
    config.curr = ant_positions[amino_acid_index];
    config.next = ant_positions[amino_acid_index - 1];
    config.prev = ant_positions[amino_acid_index + 1];
    config.amino_acid_index = amino_acid_index;

    right_neighbor = lmatrix_read(lattice, new_coord(config.next.x + 1, config.next.y));
    left_neighbor = lmatrix_read(lattice, new_coord(config.next.x - 1, config.next.y));
    up_neighbor = lmatrix_read(lattice, new_coord(config.next.x, config.next.y + 1));
    down_neighbor = lmatrix_read(lattice, new_coord(config.next.x, config.next.y - 1));

    /*If is adjacent (right) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (right_neighbor == -1 && config.next.x + 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x + 1;
        config.f.y = config.next.y;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }
    }

    /*If is adjacent (left) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (left_neighbor == -1 && config.next.x - 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x - 1;
        config.f.y = config.next.y;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }
    }

    /*If is adjacent (up) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (up_neighbor == -1 && config.next.x != config.curr.x &&
            config.next.y + 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y + 1;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index] = config;
            ++config_index;
            ++num_configs;
        }
    }

    /*If is adjacent (down) to next amino-acid and diagonally adjacent to current amino-acid*/
    if (down_neighbor == -1 && config.next.x != config.curr.x &&
            config.next.y - 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y - 1;
        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        /*If F exists and C is empty or is equals do previous amino-acid*/
        if ((config.c.x == config.prev.x && config.c.y == config.prev.y) ||
                lmatrix_read(lattice, config.c) == -1)
        {
            configs[config_index]= config;
            ++config_index;
            ++num_configs;
        }
    }

    return num_configs;
}



int change_amino_acid_position
(
    int current_energy,
    Lmatrix *lattice,
    int *sequence,
    int amino_acid_index,
    Coord src,
    Coord dest
)
/* ===============================================================================
 * Moves a amino-acid to another position in the lattice, handling with the energy
 * ===============================================================================
 */
{

    int temp;

    if (sequence[amino_acid_index] == H)
    {
        current_energy +=
            calculate_absolute_heuristic_value(lattice, amino_acid_index,
                                               src, sequence);

        current_energy -=
            calculate_absolute_heuristic_value(lattice, amino_acid_index,
                                               dest, sequence);
    }

    temp = lmatrix_read(lattice, src);
    lattice_write(lattice, src, lmatrix_read(lattice, dest));
    lattice_write(lattice, dest, temp);

    return current_energy;
}



int apply_pull_move
(
    Ant ant,
    PM_config config,
    int *sequence,
    Lmatrix *lattice,
    int sequence_len,
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
        last_modified = sequence_len - 1;
    }

    if (config.c.x == config.prev.x && config.c.y == config.prev.y)
    {
        ant.energy = change_amino_acid_position(ant.energy, lattice, sequence,
                                                config.amino_acid_index,
                                                config.curr, config.f);

        ant.positions[config.amino_acid_index] = config.f;
        config.curr =  ant.positions[config.amino_acid_index];
    }
    else
    {
        ant.energy = change_amino_acid_position(ant.energy, lattice, sequence,
                                                previous_index, config.prev,
                                                config.c);
        ant.positions[previous_index] = config.c;
        ant.energy = change_amino_acid_position(ant.energy, lattice, sequence,
                                                config.amino_acid_index,
                                                config.curr, config.f);
        ant.positions[config.amino_acid_index] = config.f;

        i = before_previous_index;

        while ((config.pm_type == ORIGINAL && i < sequence_len - 1 && i >= 0 &&
                coords_distance(ant.positions[i], ant.positions[i + 1]) != 1) ||

                (config.pm_type == INVERSE && i > 0 && i <= sequence_len - 1 &&
                 coords_distance(ant.positions[i], ant.positions[i - 1]) != 1))
        {
            ant.energy = change_amino_acid_position(ant.energy, lattice,
                                                    sequence, i,
                                                    ant.positions[i],
                                                    config.curr);
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
        lattice_write(lattice, ant.positions[i], -1);
    }
    for (i = start; i <= end; ++i)
    {
        lattice_write(lattice, ant_positions[i], i);
    }

    return ant.energy;
}



void pull_move_search
(
    int *sequence,
    int sequence_len,
    Ant *original_ant,
    Lmatrix *lattice,
    Ant best_ant,
    Ant ant,
    PM_config *configs
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

    best_ant.energy = 0;

    //APPLY PULL-MOVE WHILE THERE IS IMPROVEMENT

    do
    {
        num_configs = 0;
        previous_energy = original_ant->energy;

        for (i = 0; i < sequence_len; ++i)
        {
            lattice_write(lattice, original_ant->positions[i], i);
        }

        //GENERATE ALL POSSIBLE PULL-MOVES CONFIGURATIONS

        for (i = 1; i < sequence_len - 1; ++i)
        {
            num_configs += generate_pull_move_configs(
                               i, lattice,
                               original_ant->positions,
                               configs, num_configs,
                               sequence_len);
        }

        /*Apply all possible pull-moves and save the best*/
        for (i = 0; i < num_configs; ++i)
        {
            /*Copy original original_ant*/
            ant.energy = original_ant->energy;
            for (j = 0; j < sequence_len; ++j)
            {
                ant.positions[j] = original_ant->positions[j];
            }

            /*Apply pull move on the copy*/
            ant.energy = apply_pull_move(ant, configs[i], sequence, lattice,
                                         sequence_len, original_ant->positions);

            if (ant.energy < best_ant.energy)
            {
                best_ant.energy = ant.energy;

                for (j = 0; j < sequence_len; ++j)
                {
                    best_ant.positions[j] = ant.positions[j];
                }
            }
        }

        //CHECKS IF NEW CONFORMATION IS BETTER THAN ORIGINAL CONFORMATION

        if (best_ant.energy != 0 && best_ant.energy < original_ant->energy)
        {
            Coord distance_to_lattice_center;
            distance_to_lattice_center = new_coord(sequence_len - best_ant.positions[0].x,
                                                   sequence_len - best_ant.positions[0].y);

            original_ant->energy = previous_energy;

            lattice->top = 0;

            for (i = 0; i < sequence_len; ++i)
            {
                /*Move protein to center*/
                best_ant.positions[i].x += distance_to_lattice_center.x;
                best_ant.positions[i].y += distance_to_lattice_center.y;

                original_ant->energy = best_ant.energy;

                /*Update original original_ant*/
                original_ant->positions[i] = best_ant.positions[i];

            }
        }
    }
    while(best_ant.energy < previous_energy);

    lattice->top = 0;

}



void extract_solution(Ant ant, Solution *solution, int sequence_len)
/* ===============================================
 * Extracts conformation info from Ant to Solution
 * ===============================================s
 */
{
    Coord move, prev_move;
    int i, direction;

    solution->energy = ant.energy;

    for (i = 0; i < sequence_len - 1; ++i)
    {
        move = subtract_coord(ant.positions[i + 1], ant.positions[i]);

        if (i == 0)
        {
            direction = calculate_absolute_direction_by_move(move);
        }
        else
        {
            direction = calculate_direction_by_move(prev_move, move);
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
}



void construct_conform
(
    ACO_config aco_config,
    double **pheromone,
    Lmatrix *lattice,
    int *sequence,
    int sequence_len,
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
    double normalized_heuristic;
    Coord curr_position;
    Coord move;
    Coord candidate_move[3];
    Candidate candidates[3];
    Ant copied_ant;


    ant->energy = 0;
    ant->energy_by_edge[0] = 0;

    //DEFINES FIRST EDGE

    lattice_write(lattice, new_coord(sequence_len, sequence_len), 0);
    ant->positions[0] = new_coord(sequence_len, sequence_len);

    move = new_coord(0, 1);
    curr_position = new_coord(ant->positions[0].x + move.x,
                              ant->positions[0].y + move.y);
    lattice_write(lattice, curr_position, 1);
    ant->positions[1] = curr_position;


    //CONSTRUCTOR LOOP

    //For each edge except the first
    for (i = 1; i < sequence_len - 1; ++i)
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
            if (lmatrix_read(lattice, new_coord(curr_position.x + candidate_move[j].x,
                                                curr_position.y + candidate_move[j].y)) == -1)
            {
                candidates[num_candidates].move = candidate_move[j];

                candidates[num_candidates].position.x = curr_position.x +
                                                        candidates[num_candidates].move.x;

                candidates[num_candidates].position.y = curr_position.y +
                                                        candidates[num_candidates].move.y;

                if (sequence[i + 1] == H)
                {
                    candidates[num_candidates].heuristic =
                        calculate_absolute_heuristic_value(lattice,
                                                           i + 1,
                                                           candidates[num_candidates].position,
                                                           sequence);
                }
                else
                {
                    candidates[num_candidates].heuristic = 0;
                }

                normalized_heuristic = exp((double) candidates[num_candidates].heuristic / 0.3);
                probabilities[num_candidates] = pow(pheromone[i][j], aco_config.alpha) *
                                                pow(normalized_heuristic, aco_config.beta);

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
            lattice_write(lattice, curr_position, i + 1);

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
                    lattice->top = 0;

                    i = 0;
                    ant->energy = 0;

                    //Reset first edge
                    lattice_write(lattice, new_coord(sequence_len, sequence_len), 0);
                    ant->positions[0] = new_coord(sequence_len, sequence_len);

                    move = new_coord(0, 1);
                    curr_position = new_coord(ant->positions[0].x + move.x,
                                              ant->positions[0].y + move.y);
                    lattice_write(lattice, curr_position, 1);
                    ant->positions[1] = curr_position;
                }
                else
                {
                    copied_ant = ants[best_ant_by_index[i]];

                    lattice->top = 0;

                    //Copy best ant for index i until i th amino-acid*/
                    for (j = 0; j <= i; ++j)
                    {
                        lattice_write(lattice, copied_ant.positions[j], j);
                        ant->positions[j] = copied_ant.positions[j];
                        if (j != i)
                        {
                            ant->energy_by_edge[j] = copied_ant.energy_by_edge[j];
                        }
                    }

                    curr_position = copied_ant.positions[i];
                    ant->energy = copied_ant.energy_by_edge[i];

                    move.x = curr_position.x - copied_ant.positions[i - 1].x;
                    move.y = curr_position.y - copied_ant.positions[i - 1].y;

                    --i;
                }
                break;

            default:
                break;
            }
        }
    }

    for (i = 0; i < sequence_len - 1; ++i)
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
    int *sequence,
    int sequence_len,
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
    int *best_ant_by_edge;
    Lmatrix lattice;
    Ant *ants;
    Ant iteration_ant;
    Ant best_ant;
    Ant pm_best_ant;
    Ant pm_ant;
    PM_config* pm_configs;
    double** pheromone;
    Solution solution;


    //Sets seed
    if (*seed == -1)
    {
        *seed = (unsigned) time(NULL);
    }
    srand(*seed);

    variables_initialization(aco_config, &pm_configs, &pm_best_ant, &pm_ant,
                             &ants, &best_ant, sequence_len, &best_ant_by_edge,
                             &lattice, &pheromone, &solution);

    best_ant.energy = 0;

    //Iteration loop
    for (i = 0; i < aco_config.iterations; ++i)
    {
        //Population loop
        for (j = 0; j < aco_config.population; ++j)
        {
            construct_conform(aco_config, pheromone, &lattice, sequence,
                              sequence_len, &ants[j], best_ant_by_edge, j, ants);

            //Clean lattice
            lattice.top = 0;
        }

        //reset best ants by edge
        for (j = 0; j < sequence_len - 1; ++j)
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

                pull_move_search(sequence, sequence_len, &ants[j],
                                 &lattice, pm_best_ant, pm_ant,
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
            for (j = 0; j < sequence_len; ++j)
            {
                best_ant.positions[j] = iteration_ant.positions[j];
            }
        }

        //Pheromone update
        pheromone_evaporation(pheromone, sequence_len, aco_config.persistence);
        for (j = 0; j < aco_config.population; ++j)
        {
            pheromone_deposit(pheromone, ants[j], sequence_len, best_ant.energy);
        }
    }

    extract_solution(best_ant, &solution, sequence_len);

    free_variables(aco_config, sequence_len, &lattice, best_ant_by_edge, pm_ant,
                   pm_best_ant, ants, pheromone, pm_configs);

    return solution;
}
