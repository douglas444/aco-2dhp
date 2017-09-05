#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aco.h"

#define DEFAULT_BEST_ENERGY 1

typedef struct pull_move_config
{

    int nr_seq;
    Coord curr, prev, next, f, c;

} Pull_move_config;

int generate_pull_move_configs(int nr_seq, int **lattice, Coord *original_positions,
                               Pull_move_config *result_array, int index_on_result_array,
                               int seq_len);

int apply_pull_move(Coord *positions, Direction *conform, int *energy_by_link,
                    int energy, Pull_move_config config, int *seq, int **lattice, int seq_len,
                    Coord *original_positions);

int calculate_best_pull_move(Direction *original_conform, Coord *original_positions,
                             int *original_energy_by_link, int seq_len, int original_energy,
                             int **lattice, int *seq, Direction *pull_conform,
                             Coord *pull_positions, int *pull_energy_by_link,
                             Direction *pull_aux_conform, Coord *pull_aux_positions,
                             int *pull_aux_energy_by_link, Pull_move_config *possible_configs);

void* memory_allocation(int mem_size);

Coord calculate_move_by_direction(Coord prev_move, Direction direction);

int calculate_absolute_heuristic_value(int **lattice, int nr_seq, Coord pos, int *seq);

int random_select(float *probabilities, int len);

void pheromone_deposit(float **pheromone, Direction *conform, int energy, int seq_len,
                       int best_energy);
void pheromone_evaporation(float **pheromone, int seq_len,
                           float evaporation_rate);

int construct_conform(int *seq, Direction *conform, int *energy_by_link, Coord *conform_positions,
                      float alpha, float beta,int seq_len, float **pheromone, int **lattice,
                      Direction *best_confom, int *best_energy_by_link, Coord *best_positions,
                      int best_energy);

Direction calculate_direction_by_move(Coord prev_move, Coord curr_move);

Coord subtract_coord(Coord c1, Coord c2);

///Main function
int aco_run(int *seq, int seq_len, float alpha, float beta, float evaporation_rate,
            float initial_phero, int population_size, int num_iterations, int stop_criterion)
{
    //Declarations
    int i, j, **lattice,
        *ants_energy, best_energy, iteration_energy, iteration_indice,
        **ants_energy_by_link, *best_energy_by_link, *best_pull_energy_by_link, *pull_energy_by_link;

    Coord **ants_positions, *best_positions, *best_pull_positions, *pull_positions;
    Direction **ants_conform, *best_conform, *best_pull_conform, *pull_conform;
    Pull_move_config *possible_configs;

    float **pheromone;

    //Allocations
    best_energy_by_link = (int*) memory_allocation(sizeof(int) * (seq_len - 1));
    best_pull_energy_by_link = (int*) memory_allocation(sizeof(int) * (seq_len - 1));
    pull_energy_by_link = (int*) memory_allocation(sizeof(int) * (seq_len - 1));

    best_positions = (Coord*) memory_allocation(sizeof(Coord) * seq_len);
    best_pull_positions = (Coord*) memory_allocation(sizeof(Coord) * seq_len);
    pull_positions = (Coord*) memory_allocation(sizeof(Coord) * seq_len);

    best_conform = (Direction*) memory_allocation(sizeof(Direction) * (seq_len - 1));
    best_pull_conform = (Direction*) memory_allocation(sizeof(Direction) * (seq_len - 1));
    pull_conform = (Direction*) memory_allocation(sizeof(Direction) * (seq_len - 1));

    possible_configs = (Pull_move_config*) memory_allocation(sizeof(Pull_move_config) * 2 *
                       (seq_len - 2));

    ants_energy = (int*) malloc(sizeof(int) * population_size);
    ants_energy_by_link = (int**) memory_allocation(sizeof(int*) * population_size);
    ants_positions = (Coord**) memory_allocation(sizeof(Coord*) * population_size);
    ants_conform = (Direction**) memory_allocation(sizeof(Direction*) * population_size);

    for (i = 0; i < population_size; ++i)
    {
        ants_positions[i] = (Coord*) memory_allocation(sizeof(Coord) * seq_len);
        ants_energy_by_link[i] = (int*) memory_allocation(sizeof(int) * (seq_len - 1));
        ants_conform[i] = (Direction*) memory_allocation(sizeof(Direction) * (seq_len - 1));
    }

    pheromone = (float**) memory_allocation(sizeof(float*) * seq_len);
    for (i = 0; i < seq_len; ++i)
    {
        pheromone[i] = (float*) memory_allocation(sizeof(float) * 3);
        for (j = 0; j < 3; ++j)
        {
            pheromone[i][j] = initial_phero;
        }
    }

    lattice = (int**) malloc(sizeof(int*) * (2 * seq_len + 1));
    for (i = 0; i < 2 * seq_len + 1; ++i)
    {
        lattice[i] = (int*) malloc(sizeof(int) * (2 * seq_len + 1));
        for (j = 0; j < 2 * seq_len + 1; ++j)
        {
            lattice[i][j] = -1;
        }
    }

    ///ALGORITHM BEGIN

    best_energy = DEFAULT_BEST_ENERGY;

    for (i = 0; i < num_iterations; ++i)
    {
        iteration_indice = 0;
        iteration_energy = DEFAULT_BEST_ENERGY;

        for (j = 0; j < population_size; ++j)
        {

            ants_energy[j] = construct_conform(seq, ants_conform[j], ants_energy_by_link[j], ants_positions[j],
                                               alpha, beta, seq_len, pheromone, lattice,
                                               ants_conform[iteration_indice], ants_energy_by_link[iteration_indice],
                                               ants_positions[iteration_indice], iteration_energy);

            ants_energy[j] = calculate_best_pull_move(ants_conform[j], ants_positions[j], ants_energy_by_link[j],
                             seq_len, ants_energy[j], lattice, seq, best_pull_conform,
                             best_pull_positions, best_pull_energy_by_link,
                             pull_conform, pull_positions,
                             pull_energy_by_link, possible_configs);

            if (ants_energy[j] < iteration_energy)
            {
                iteration_indice = j;
                iteration_energy = ants_energy[j];
            }

            if (best_energy <= stop_criterion)
            {
                break;
            }
        }

        if (iteration_energy < best_energy)
        {
            best_energy = ants_energy[iteration_indice];
            for (j = 0; j < seq_len; ++j)
            {
                best_positions[j] = ants_positions[iteration_indice][j];
                if (j < seq_len - 1)
                {
                    best_conform[j] = ants_conform[iteration_indice][j];
                    best_energy_by_link[j] = ants_energy_by_link[iteration_indice][j];
                }

            }
        }

        pheromone_evaporation(pheromone, seq_len, evaporation_rate);

        for (j = 0; j < population_size; ++j)
        {
            pheromone_deposit(pheromone, ants_conform[j], ants_energy[j], seq_len, stop_criterion);
        }


        if (best_energy <= stop_criterion)
        {
            break;
        }
    }

    ///ALGORITHM END

    //output file
    int min_x, min_y, max_x, max_y;
    FILE *f = fopen("../protein_data", "w");
    if (f == NULL)
    {
        printf("Error in function aco_run: Unable to open output file");
    }

    min_x = best_positions[0].x;
    min_y = best_positions[0].y;
    max_x = best_positions[0].x;
    max_y = best_positions[0].y;
    for (i = 1; i < seq_len; ++i)
    {
        if (best_positions[i].x < min_x)
        {
            min_x = best_positions[i].x;
        }
        if (best_positions[i].y < min_y)
        {
            min_y = best_positions[i].y;
        }
        if (best_positions[i].x > max_x)
        {
            max_x = best_positions[i].x;
        }
        if (best_positions[i].y > max_y)
        {
            max_y = best_positions[i].y;
        }
    }

    fprintf(f, "CONFORMATION DIMENSION: %d, %d\n", max_x - min_x, max_y - min_y);
    fprintf(f, "NUMBER OF AMINO_ACIDS: %d\n", seq_len);


    for (i = 0; i < seq_len; ++i)
    {
        fprintf(f, "%d %d %d\n", seq[i], best_positions[i].x - min_x,
                best_positions[i].y - min_y);
    }
    fclose(f);

    //Free memory
    free(best_conform);
    free(best_pull_conform);
    free(pull_conform);

    free(best_energy_by_link);
    free(best_pull_energy_by_link);
    free(pull_energy_by_link);

    free(best_positions);
    free(best_pull_positions);
    free(pull_positions);

    free(possible_configs);

    for (i = 0; i < seq_len; ++i)
    {
        free(pheromone[i]);
        free(ants_positions[i]);
        if (i < seq_len - 1)
        {
            free(ants_energy_by_link[i]);
            free(ants_conform[i]);
        }
    }
    for (i = 0; i < 2 * seq_len + 1; ++i)
    {
        free(lattice[i]);
    }
    free(lattice);
    free(pheromone);
    free(ants_conform);
    free(ants_energy_by_link);
    free(ants_positions);

    return best_energy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void* memory_allocation(int mem_size)
{
    void *mem_pos = (void*) malloc(mem_size);

    if (mem_pos == NULL)
    {
        printf("Error in function memory_allocation: Unable to allocate memory");
        exit(1);
    }
    else
    {
        return mem_pos;
    }
}

Direction calculate_direction_by_move(Coord prev_move, Coord curr_move)
{

    if (prev_move.x == curr_move.x && prev_move.y == curr_move.y)
    {
        return STRAIGHT;
    }
    else if (prev_move.y == 0)
    {
        if (curr_move.y == prev_move.x)
        {
            return LEFT;
        }
        else if (abs(curr_move.y) == abs(prev_move.x))
        {
            return RIGHT;
        }
        else
        {
            printf("Error in function calculate_direction_by_move: Invalid values for parameters");
            exit(1);
        }
    }
    else
    {
        if (curr_move.x == prev_move.y)
        {
            return RIGHT;
        }
        else if (abs(curr_move.x) == abs(prev_move.y))
        {
            return LEFT;
        }
        else
        {
            printf("Error in function calculate_direction_by_move: Invalid values for parameters");
            exit(1);
        }
    }

}

Coord calculate_move_by_direction(Coord prev_move, Direction curr_direction)
{
    Coord move;

    if (prev_move.x == 1)
    {
        switch (curr_direction)
        {
        case LEFT:
            move.x = 0;
            move.y = 1;
            break;
        case RIGHT:
            move.x = 0;
            move.y = -1;
            break;
        case STRAIGHT:
            move.x = 1;
            move.y = 0;
            break;
        default:
            printf("Error in function calculate_move_by_direction: Invalid value for parameter curr_direction");
            exit(1);
            break;
        }
    }
    else if (prev_move.x == -1)
    {
        switch (curr_direction)
        {

        case LEFT:
            move.x = 0;
            move.y = -1;
            break;
        case RIGHT:
            move.x = 0;
            move.y = 1;
            break;
        case STRAIGHT:
            move.x = -1;
            move.y = 0;
            break;
        default:
            printf("Error in function calculate_move_by_direction: Invalid value for parameter curr_direction");
            exit(1);
            break;
        }
    }
    else if (prev_move.y == -1)
    {
        switch (curr_direction)
        {
        case LEFT:
            move.x = 1;
            move.y = 0;
            break;
        case RIGHT:
            move.x = -1;
            move.y = 0;
            break;
        case STRAIGHT:
            move.x = 0;
            move.y = -1;
            break;
        default:
            printf("Error in function calculate_move_by_direction: Invalid value for parameter curr_direction");
            exit(1);
            break;
        }
    }
    else if (prev_move.y == 1)
    {
        switch (curr_direction)
        {
        case LEFT:
            move.x = -1;
            move.y = 0;
            break;
        case RIGHT:
            move.x = 1;
            move.y = 0;
            break;
        case STRAIGHT:
            move.x = 0;
            move.y = 1;
            break;
        default:
            printf("Error in function calculate_move_by_direction: Invalid value for parameter curr_direction");
            exit(1);
            break;
        }
    }
    else
    {
        printf("Error in function calculate_move_by_direction: Invalid value for parameter prev_move");
        exit(1);
    }

    return move;
}

int calculate_absolute_heuristic_value(int **lattice, int nr_seq, Coord pos, int *seq)
{
    int heuristic_value = 0, right_neighbor, left_neighbor, down_neighbor, up_neighbor;

    right_neighbor = lattice[pos.x + 1][pos.y];
    left_neighbor = lattice[pos.x - 1][pos.y];
    down_neighbor = lattice[pos.x][pos.y - 1];
    up_neighbor = lattice[pos.x][pos.y + 1];

    if (right_neighbor != -1 && abs(right_neighbor - nr_seq) > 1 && seq[right_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (left_neighbor != -1 && abs(left_neighbor - nr_seq) > 1 && seq[left_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (up_neighbor != -1 && abs(up_neighbor - nr_seq) > 1 && seq[up_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (down_neighbor != -1 && abs(down_neighbor - nr_seq) > 1 && seq[down_neighbor] == 1)
    {
        ++heuristic_value;
    }

    return heuristic_value;
}

int random_select(float *probabilities, int len)
{
    int i, result;
    double cumulative_probability, r;

    r = ((double)rand()/RAND_MAX);

    i = 0;
    result = -1;
    cumulative_probability = 0;

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

void pheromone_deposit(float **pheromone, Direction *conform, int energy, int seq_len,
                       int best_energy)
{
    int i, j;

    for (i = 0; i < seq_len - 1; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            if (conform[j] == j)
            {
                pheromone[i][j] += (float) energy / pow(best_energy, 3);
            }
        }
    }
}

void pheromone_evaporation(float **pheromone, int seq_len,
                           float evaporation_rate)
{
    int i, j;

    for (i = 0; i < seq_len - 1; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            pheromone[i][j] *= evaporation_rate;
        }
    }
}

int construct_conform(int *seq, Direction *conform, int *energy_by_link, Coord *conform_positions,
                      float alpha, float beta,int seq_len, float **pheromone, int **lattice,
                      Direction *best_confom, int *best_energy_by_link, Coord *best_positions,
                      int best_energy)
{

    int i, j, num_candidates, selected_candidate, energy, candidates_heuristics[3];
    float probabilities[3], sum_probabilities;
    Coord curr_position, move, candidate_move, candidates_positions[3], candidates_moves[3];
    Direction candidates_directions[3];

    energy = 0;

    //Defines first link

    curr_position.x = seq_len;
    curr_position.y = seq_len;
    lattice[curr_position.x][curr_position.y] = 0;
    conform[0] = STRAIGHT;
    energy_by_link[0] = 0;
    conform_positions[0] = curr_position;

    move.x = 0;
    move.y = 1;
    curr_position.x += move.x;
    curr_position.y += move.y;
    lattice[curr_position.x][curr_position.y] = 1;
    conform_positions[1] = curr_position;

    //for each aminoa-acid link, except the first
    for (i = 1; i < seq_len - 1; ++i)
    {
        sum_probabilities = 0;
        num_candidates = 0;

        //for each direction
        for (j = 0; j < 3; ++j)
        {
            candidate_move = calculate_move_by_direction(move, j);

            //if the next position in this direction is not occuped
            if (lattice[curr_position.x + candidate_move.x][curr_position.y + candidate_move.y] == -1)
            {
                candidates_moves[num_candidates] = candidate_move;

                candidates_directions[num_candidates] = j;

                candidates_positions[num_candidates].x = curr_position.x +
                        candidates_moves[num_candidates].x;

                candidates_positions[num_candidates].y = curr_position.y +
                        candidates_moves[num_candidates].y;

                if (seq[i + 1] == 1)
                {
                    candidates_heuristics[num_candidates] = calculate_absolute_heuristic_value(lattice,
                                                            i + 1,
                                                            candidates_positions[num_candidates],
                                                            seq);
                }
                else
                {
                    candidates_heuristics[num_candidates] = 0;
                }


                probabilities[num_candidates] = pow(pheromone[i][j], alpha) *
                                                pow(exp((float) candidates_heuristics[num_candidates] / 0.3), beta);

                sum_probabilities += probabilities[num_candidates];

                ++num_candidates;
            }

        }

        //calculate probabilities
        for (j = 0; j < num_candidates; ++j)
        {
            probabilities[j] = probabilities[j]/sum_probabilities;
        }

        //Selects the direction
        if (num_candidates == 0)
        {
            selected_candidate = -1;
        }
        else if (sum_probabilities == 0)
        {
            selected_candidate = rand()%num_candidates;
        }
        else if (num_candidates > 1)
        {
            selected_candidate = random_select(probabilities, num_candidates);
        }
        else if (num_candidates == 1)
        {
            selected_candidate = 0;
        }

        //Updates conform
        if (selected_candidate != -1)
        {
            move = candidates_moves[selected_candidate];
            curr_position = candidates_positions[selected_candidate];
            conform[i + 1] = candidates_directions[selected_candidate];
            energy_by_link[i] = energy;
            energy -= candidates_heuristics[selected_candidate];

            conform_positions[i + 1] = curr_position;

            lattice[curr_position.x][curr_position.y] = i + 1;
        }
        else//impossible continue the fold process
        {
            //if theres no another conformation to copy
            if (best_energy == DEFAULT_BEST_ENERGY)
            {
                for (j = 0; j <= i; ++j)
                {
                    lattice[conform_positions[j].x][conform_positions[j].y] = -1;
                }

                i = 0;
                energy = 0;

                //Defines first link

                curr_position.x = seq_len;
                curr_position.y = seq_len;
                lattice[curr_position.x][curr_position.y] = 0;
                conform[0] = STRAIGHT;
                energy_by_link[0] = 0;
                conform_positions[0] = curr_position;

                move.x = 0;
                move.y = 1;
                curr_position.x += move.x;
                curr_position.y += move.y;
                lattice[curr_position.x][curr_position.y] = 1;
                conform_positions[1] = curr_position;
            }
            else
            {
                for (j = 0; j <= i; ++j)
                {
                    lattice[conform_positions[j].x][conform_positions[j].y] = -1;
                }

                for (j = 0; j <= i; ++j)
                {
                    lattice[best_positions[j].x][best_positions[j].y] = j;
                    conform_positions[j] = best_positions[j];
                    if (j != i)
                    {
                        conform[j] = best_confom[j];
                        energy_by_link[j] = best_energy_by_link[j];
                    }
                }
                curr_position = best_positions[i];
                energy = best_energy_by_link[i];
                move.x = curr_position.x - best_positions[i - 1].x;
                move.y = curr_position.y - best_positions[i - 1].y;
                --i;
            }
        }
    }

    return energy;
}

Coord subtract_coord(Coord c1, Coord c2)
{

    Coord c3;

    c3.x = c1.x - c2.x;
    c3.y = c1.y - c2.y;

    return c3;

}

int coords_distance(Coord c1, Coord c2)
{

    return sqrt(pow(c1.x - c2.x, 2) + pow(c1.y - c2.y, 2));

}

Direction calculate_absolute_direction_by_move(Coord move)
{

    if (move.x == 0)
    {
        if (move.y == 1 || move.y == -1)
        {
            return STRAIGHT;
        }
        else
        {
            printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move");
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
            printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move");
            exit(1);
        }
    }
    else
    {
        printf("Error in function calculate_absolute_direction_by_move: Invalid value for parameter move");
        exit(1);
    }
}


int calculate_relative_heuristic_value(int **lattice, int nr_seq, Coord pos, int *seq)
{
    int heuristic_value = 0, right_neighbor, left_neighbor, down_neighbor, up_neighbor;

    right_neighbor = lattice[pos.x + 1][pos.y];
    left_neighbor = lattice[pos.x - 1][pos.y];
    down_neighbor = lattice[pos.x][pos.y - 1];
    up_neighbor = lattice[pos.x][pos.y + 1];

    if (right_neighbor != -1 && right_neighbor < nr_seq - 1 && seq[right_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (left_neighbor != -1 && left_neighbor < nr_seq - 1 && seq[left_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (up_neighbor != -1 && up_neighbor < nr_seq - 1 && seq[up_neighbor] == 1)
    {
        ++heuristic_value;
    }
    if (down_neighbor != -1 && down_neighbor < nr_seq - 1 && seq[down_neighbor] == 1)
    {
        ++heuristic_value;
    }

    return heuristic_value;
}

int generate_pull_move_configs(int nr_seq, int **lattice, Coord *original_positions,
                               Pull_move_config *result_array, int index_on_result_array,
                               int seq_len)
{

    int num_configs;
    Pull_move_config config;

    num_configs = 0;

    config.curr = original_positions[nr_seq];
    config.next = original_positions[nr_seq + 1];
    config.prev = original_positions[nr_seq - 1];
    config.nr_seq = nr_seq;

    //if is adjacent (right) to next amino-acid and diagonally adjacent to current amino-acid
    if (lattice[config.next.x + 1][config.next.y] == -1 && config.next.x + 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x + 1;
        config.f.y = config.next.y;

        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        //if F exists and C is empty or is equals do previous amino-acid
        if (lattice[config.c.x][config.c.y] == -1 ||
                (config.c.x == config.prev.x && config.c.y == config.prev.y))
        {
            result_array[index_on_result_array].c = config.c;
            result_array[index_on_result_array].curr = config.curr;
            result_array[index_on_result_array].f = config.f;
            result_array[index_on_result_array].next = config.next;
            result_array[index_on_result_array].nr_seq = config.nr_seq;
            result_array[index_on_result_array].prev = config.prev;

            ++index_on_result_array;
            ++num_configs;
        }

    }

    //if is adjacent (left) to next amino-acid and diagonally adjacent to current amino-acid
    if (lattice[config.next.x - 1][config.next.y] == -1 && config.next.x - 1 != config.curr.x &&
            config.next.y != config.curr.y)
    {
        config.f.x = config.next.x - 1;
        config.f.y = config.next.y;

        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        //if F exists and C is empty or is equals do previous amino-acid
        if (lattice[config.c.x][config.c.y] == -1 ||
                (config.c.x == config.prev.x && config.c.y == config.prev.y))
        {
            result_array[index_on_result_array].c = config.c;
            result_array[index_on_result_array].curr = config.curr;
            result_array[index_on_result_array].f = config.f;
            result_array[index_on_result_array].next = config.next;
            result_array[index_on_result_array].nr_seq = config.nr_seq;
            result_array[index_on_result_array].prev = config.prev;

            ++index_on_result_array;
            ++num_configs;
        }
    }

    //if is adjacent (up) to next amino-acid and diagonally adjacent to current amino-acid
    if (lattice[config.next.x][config.next.y + 1] == -1 && config.next.x != config.curr.x &&
            config.next.y + 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y + 1;

        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        //if F exists and C is empty or is equals do previous amino-acid
        if (lattice[config.c.x][config.c.y] == -1 ||
                (config.c.x == config.prev.x && config.c.y == config.prev.y))
        {
            result_array[index_on_result_array].c = config.c;
            result_array[index_on_result_array].curr = config.curr;
            result_array[index_on_result_array].f = config.f;
            result_array[index_on_result_array].next = config.next;
            result_array[index_on_result_array].nr_seq = config.nr_seq;
            result_array[index_on_result_array].prev = config.prev;

            ++index_on_result_array;
            ++num_configs;
        }
    }

    //if is adjacent (down) to next amino-acid and diagonally adjacent to current amino-acid
    if (lattice[config.next.x][config.next.y - 1] == -1 && config.next.x != config.curr.x &&
            config.next.y - 1 != config.curr.y)
    {
        config.f.x = config.next.x;
        config.f.y = config.next.y - 1;

        config.c.x = config.f.x + config.curr.x - config.next.x;
        config.c.y = config.f.y + config.curr.y - config.next.y;

        //if F exists and C is empty or is equals do previous amino-acid
        if (lattice[config.c.x][config.c.y] == -1 ||
                (config.c.x == config.prev.x && config.c.y == config.prev.y))
        {
            result_array[index_on_result_array].c = config.c;
            result_array[index_on_result_array].curr = config.curr;
            result_array[index_on_result_array].f = config.f;
            result_array[index_on_result_array].next = config.next;
            result_array[index_on_result_array].nr_seq = config.nr_seq;
            result_array[index_on_result_array].prev = config.prev;

            ++index_on_result_array;
            ++num_configs;
        }
    }

    return num_configs;

}

int apply_pull_move(Coord *positions, Direction *conform, int *energy_by_link,
                    int energy, Pull_move_config config, int *seq, int **lattice, int seq_len,
                    Coord *original_positions)
{
    int temp, i, last_modified;

    if (config.c.x == config.prev.x && config.c.y == config.prev.y)
    {
        //config.curr <-> f
        if (seq[config.nr_seq] == 1)
        {
            energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
            energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
        }

        temp = lattice[config.curr.x][config.curr.y];
        lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
        lattice[config.f.x][config.f.y] = temp;

        positions[config.nr_seq].x = config.f.x;
        positions[config.nr_seq].y = config.f.y;

        config.curr.x =  positions[config.nr_seq].x;
        config.curr.y =  positions[config.nr_seq].y;

        if (config.nr_seq == 1)
        {
            conform[0] =
                calculate_absolute_direction_by_move(
                    subtract_coord(positions[1], positions[0]));
        }
        else
        {
            conform[config.nr_seq - 1] =
                calculate_direction_by_move(subtract_coord(config.prev,  positions[config.nr_seq - 2]),
                                            subtract_coord(config.curr, config.prev));
        }

        conform[config.nr_seq] =
            calculate_direction_by_move(subtract_coord(config.curr, config.prev),
                                        subtract_coord(config.next, config.curr));

        last_modified = config.nr_seq;
    }
    else
    {
        //config.prev <-> c
        if (seq[config.nr_seq - 1] == 1)
        {
            energy += calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, config.prev, seq);
            energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, config.c, seq);
        }

        temp = lattice[config.prev.x][config.prev.y];
        lattice[config.prev.x][config.prev.y] = lattice[config.c.x][config.c.y];
        lattice[config.c.x][config.c.y] = temp;
        positions[config.nr_seq - 1].x = config.c.x;
        positions[config.nr_seq - 1].y = config.c.y;

        //config.curr <-> f
        if (seq[config.nr_seq] == 1)
        {
            energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
            energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
        }

        temp = lattice[config.curr.x][config.curr.y];
        lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
        lattice[config.f.x][config.f.y] = temp;
        positions[config.nr_seq].x = config.f.x;
        positions[config.nr_seq].y = config.f.y;

        i = config.nr_seq - 2;
        while (i >= 0)
        {
            //i <- config.curr <- config.prev
            if (seq[i] == 1)
            {
                energy += calculate_absolute_heuristic_value(lattice, i,  positions[i], seq);
                energy -= calculate_absolute_heuristic_value(lattice, i, config.curr, seq);
            }

            temp = lattice[ positions[i].x][ positions[i].y];
            lattice[ positions[i].x][ positions[i].y] =
                lattice[config.curr.x][config.curr.y];
            lattice[config.curr.x][config.curr.y] = temp;

            temp =  positions[i].x;
            positions[i].x = config.curr.x;
            config.curr.x = config.prev.x;
            config.prev.x = temp;

            temp = positions[i].y;
            positions[i].y = config.curr.y;
            config.curr.y = config.prev.y;
            config.prev.y = temp;

            if (coords_distance( positions[i],  positions[i - 1]) == 1)
            {
                last_modified = i;
                break;
            }
            --i;
        }

        if (i == -1)
        {
            conform[0] = calculate_absolute_direction_by_move(
                             subtract_coord( positions[1],  positions[0]));
            last_modified = 0;
            i = 1;
        }

        while (i < config.nr_seq + 1)
        {
            conform[i] = calculate_direction_by_move(
                             subtract_coord( positions[i],
                                             positions[i - 1]),
                             subtract_coord( positions[i + 1],
                                             positions[i]));
            ++i;
        }
    }

    for (i = last_modified; i < seq_len - 1; ++i)
    {
        if (i == 0)
        {
            energy_by_link[i] = 0;
        }
        else
        {
            if (seq[i] == 1)
            {
                energy_by_link[i] = energy_by_link[i - 1] -
                                    calculate_relative_heuristic_value(lattice, i,  positions[i],
                                            seq);
            }
            else
            {
                energy_by_link[i] = energy_by_link[i - 1];
            }
        }
    }
    for (i = last_modified; i < seq_len - 1; ++i)
    {
        lattice[positions[i].x][positions[i].y] = -1;
        lattice[original_positions[i].x][original_positions[i].y] = i;
    }
    return energy;
}

int calculate_best_pull_move(Direction *original_conform, Coord *original_positions,
                             int *original_energy_by_link, int seq_len, int original_energy,
                             int **lattice, int *seq, Direction *best_pull_conform,
                             Coord *best_pull_positions, int *best_pull_energy_by_link,
                             Direction *pull_conform, Coord *pull_positions,
                             int *pull_energy_by_link, Pull_move_config *possible_configs)
{

    int i, j, num_possible_configs, best_pull_energy, pull_energy;

    num_possible_configs = 0;

    //Generates all possible pull-moves configs
    for (i = seq_len - 2; i > 0; --i)
    {
        num_possible_configs += generate_pull_move_configs(i, lattice, original_positions,
                                possible_configs, num_possible_configs, seq_len);
    }

    best_pull_energy = 1;//Default value

    //Apply all possible pull-moves and save the best
    for (i = 0; i < num_possible_configs; ++i)
    {
        //Copy original conformation
        pull_energy = original_energy;
        for (j = 0; j < seq_len; ++j)
        {
            pull_positions[j] = original_positions[j];
            if (j < seq_len - 1)
            {
                pull_conform[j] = original_conform[j];
                pull_energy_by_link[j] = original_energy_by_link[j];
            }
        }

        //Apply pull move on the copy
        pull_energy = apply_pull_move(pull_positions, pull_conform,
                                      pull_energy_by_link, pull_energy,
                                      possible_configs[i], seq, lattice, seq_len,
                                      original_positions);

        //If the energy of the new conformation is lower than the actual best, update the best
        if (pull_energy < best_pull_energy)
        {
            best_pull_energy = pull_energy;

            for (j = 0; j < seq_len; ++j)
            {
                best_pull_positions[j] = pull_positions[j];
                if (j < seq_len - 1)
                {
                    best_pull_conform[j] = pull_conform[j];
                    best_pull_energy_by_link[j] = pull_energy_by_link[j];
                }
            }
        }

    }

    /* If the best found pull-moved conformation is better than the original conformation,
    update the original and return the new energy */
    if (best_pull_energy <= original_energy)
    {
        Coord distance_to_lattice_center;

        distance_to_lattice_center.x = seq_len - best_pull_positions[0].x;
        distance_to_lattice_center.y = seq_len - best_pull_positions[0].y;

        for (i = 0; i < seq_len; ++i)
        {
            //Clean lattice
            lattice[original_positions[i].x][original_positions[i].y] = -1;

            //Move protein to center
            best_pull_positions[i].x += distance_to_lattice_center.x;
            best_pull_positions[i].y += distance_to_lattice_center.y;

            //Update original conformation
            original_positions[i] = best_pull_positions[i];
            if (i < seq_len - 1)
            {
                original_conform[i] = best_pull_conform[i];
                original_energy_by_link[i] = best_pull_energy_by_link[i];
            }

        }
        return best_pull_energy;
    }
    else
    {
        //Clean lattice
        for (i = 0; i < seq_len; ++i)
        {
            lattice[original_positions[i].x][original_positions[i].y] = -1;
        }
        return original_energy;
    }

}
