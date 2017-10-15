#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aco.h"


typedef struct pull_move_config
{
  int nr_seq;
  Coord curr;
  Coord prev;
  Coord next;
  Coord f;
  Coord c;

} Pull_move_config;

typedef struct end_move_config
{
  int nr_seq;
  Coord adj;
  Coord away;
  Coord curr;

} End_move_config;

int       calculate_absolute_heuristic_value   (int **lattice, int nr_seq, Coord pos, int *seq);
int       random_select                        (double *probabilities, int len);
int       construct_conform                    (int *seq, int seq_len, Conformation solution, Aco_config aco_config, double **pheromone, int **lattice, int *best_partial_indexes, int ant_index, Conformation *ants_conformation, Coord *blocked_positons, int num_blocked_positions);
int       calculate_relative_heuristic_value   (int **lattice, int nr_seq, Coord pos, int *seq);
int       coords_distance                      (Coord c1, Coord c2);
void      pheromone_deposit                    (double **pheromone, Conformation conform, int seq_len, int best_energy);
void      pheromone_evaporation                (double **pheromone, int seq_len, double evaporation_rate);
void      initialize_conform                   (Conformation *c, int seq_len);
void      destroy_conform                      (Conformation c);
void*     memory_allocation                    (int mem_size);
void      output_plot_file                     (int *seq, int seq_len, Conformation best_solution);
int       generate_pull_move_configs           (int nr_seq, int **lattice, Coord *ant_positions, Pull_move_config *configs, int config_index, int seq_len);
int       apply_pull_move                      (Conformation conformation, Pull_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions);
int       apply_pull_move_inverse              (Conformation conformation, Pull_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions);
int       calculate_best_pull_move             (int *seq, int seq_len, Conformation ant_conformation,int **lattice, Conformation best_conformation, Conformation config_conformation, Pull_move_config *configs, Pull_move_config *configs_inverse);
int       apply_tail_pull_move                 (Conformation conformation, int *seq, int **lattice, int seq_len);
int       apply_tail_pull_move_inverse         (Conformation conformation, int *seq, int **lattice, int seq_len);
Coord     calculate_move_by_direction          (Coord prev_move, Direction curr_direction);
Coord     subtract_coord                       (Coord c1, Coord c2);
Direction calculate_direction_by_move          (Coord prev_move, Coord curr_move);
Direction calculate_absolute_direction_by_move (Coord move);

/*********************************************************************************************************************/


Conformation aco_run(int *seq, int seq_len, Aco_config aco_config)
{

  /**VARIABLES DECLARATIONS*/

  int i;
  int j;
  int k;
  int num_blocked_positions;
  int iteration_best_energy;
  int iteration_best_index;
  int* best_partials_indexes;
  int** lattice;
  double** pheromone;
  Conformation* ants_conformation;
  Conformation best_conformation;
  Coord *blocked_positons;

  struct
  {
    /*Used only in calculate_best_pull_move function,
    but defined here to economize memory_allocation calls*/
    Conformation best_conformation;
    Conformation config_conformation;
    Pull_move_config* configs;
    Pull_move_config* configs_inverse;

  } pull_move_vars;


  /**MEMORY ALLOCATIONS AND VARIABLES INITIALIZATIONS*/

  best_partials_indexes = (int*) memory_allocation(sizeof(int) * (seq_len - 1));
  lattice = (int**) memory_allocation(sizeof(int*) * (2 * seq_len + 1));
  pheromone = (double**) memory_allocation(sizeof(double*) * seq_len);
  ants_conformation = (Conformation*) memory_allocation(sizeof(Conformation) * aco_config.population);
  pull_move_vars.configs = (Pull_move_config*) memory_allocation(sizeof(Pull_move_config) * 2 * (seq_len - 2));
  pull_move_vars.configs_inverse = (Pull_move_config*) memory_allocation(sizeof(Pull_move_config) * 2 * (seq_len - 2));
  blocked_positons = (Coord*) memory_allocation(sizeof(Coord) * (2 * seq_len + 1) * (2 * seq_len + 1));

  initialize_conform(&pull_move_vars.best_conformation, seq_len);
  initialize_conform(&pull_move_vars.config_conformation, seq_len);
  initialize_conform(&best_conformation, seq_len);

  for (i = 0; i < seq_len; ++i)
  {
    pheromone[i] = (double*) memory_allocation(sizeof(double) * 3);
    if (i < seq_len - 1)
    {
      best_partials_indexes[i] = -1;
    }
    for (j = 0; j < 3; ++j)
    {
      pheromone[i][j] = aco_config.ini_pheromone;
    }
  }

  for (i = 0; i < 2 * seq_len + 1; ++i)
  {
    lattice[i] = (int*) memory_allocation(sizeof(int) * (2 * seq_len + 1));
    for (j = 0; j < 2 * seq_len + 1; ++j)
    {
      lattice[i][j] = -1;
    }
  }

  for (i = 0; i < aco_config.population; ++i)
  {
    initialize_conform(&ants_conformation[i], seq_len);
  }

  best_conformation.energy = 0;


  /**ALGORITHM BEGIN*/

  for (i = 0; i < aco_config.iterations; ++i)
  {

    for (j = 0; j < aco_config.population; ++j)
    {
      num_blocked_positions = 0;

      /*Constructs ant conformation*/
      ants_conformation[j].energy = construct_conform(seq,
                                                      seq_len,
                                                      ants_conformation[j],
                                                      aco_config,
                                                      pheromone,
                                                      lattice,
                                                      best_partials_indexes,
                                                      j,
                                                      ants_conformation,
                                                      blocked_positons,
                                                      num_blocked_positions);
      /*Clean lattice*/

      for (k = 0; k < num_blocked_positions; ++k)
      {
        lattice[blocked_positons[k].x][blocked_positons[k].y] = -1;
      }

      for (k = 0; k < seq_len; ++k)
      {
        lattice[ants_conformation[j].positions[k].x][ants_conformation[j].positions[k].y] = -1;
      }

    }

    for (j = 0; j < seq_len - 1; ++j)
    {
      best_partials_indexes[j] = -1;
    }

    /**LOCAL SEARCH AND SELECTION OF BEST ANT SOLUTION*/

    iteration_best_index = 0;
    iteration_best_energy = 0;
    for (j = 0; j < aco_config.population; ++j)
    {

      switch (aco_config.local_search)
      {
        case WITHOUT_LOCAL_SEARCH:
          break;

        case PULL_MOVE:
          ants_conformation[j].energy = calculate_best_pull_move(seq, seq_len, ants_conformation[j], lattice,
                                                                  pull_move_vars.best_conformation,
                                                                  pull_move_vars.config_conformation,
                                                                  pull_move_vars.configs,
                                                                  pull_move_vars.configs_inverse);
          break;

        default:
          break;
      }

      if (ants_conformation[j].energy < iteration_best_energy)
      {
        iteration_best_index = j;
        iteration_best_energy = ants_conformation[j].energy;
      }

    }

    /*Update best iteration conformation*/
    if (iteration_best_energy < best_conformation.energy)
    {
      best_conformation.energy = ants_conformation[iteration_best_index].energy;
      for (j = 0; j < seq_len; ++j)
      {
        best_conformation.positions[j] = ants_conformation[iteration_best_index].positions[j];
        if (j < seq_len - 1)
        {
          best_conformation.directions[j] =  ants_conformation[iteration_best_index].directions[j];
        }
      }
    }

    /**PHEROMONE UPDATE*/

    pheromone_evaporation(pheromone, seq_len, aco_config.persistence);

    for (j = 0; j < aco_config.population; ++j)
    {
      pheromone_deposit(pheromone, ants_conformation[j], seq_len, aco_config.best_known_solution);
    }
  }

  /**END OF ALGORITHM*/


  output_plot_file(seq, seq_len, best_conformation);


  /**FREE MEMORY*/

  for (i = 0; i < seq_len; ++i)
  {
    free(pheromone[i]);
  }
  for (i = 0; i < aco_config.population; ++i)
  {
    destroy_conform(ants_conformation[i]);
  }
  for (i = 0; i < 2 * seq_len + 1; ++i)
  {
    free(lattice[i]);
  }

  free(blocked_positons);
  free(best_partials_indexes);
  free(lattice);
  free(pheromone);
  free(ants_conformation);
  destroy_conform(pull_move_vars.best_conformation);
  destroy_conform(pull_move_vars.config_conformation);
  free(pull_move_vars.configs);
  free(pull_move_vars.configs_inverse);


  return best_conformation;
}


/*********************************************************************************************************************/


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

void initialize_conform(Conformation *c, int seq_len)
{
  c->energy_by_link = (int*) memory_allocation(sizeof(int) * (seq_len - 1));
  c->positions = (Coord*) memory_allocation(sizeof(Coord) * seq_len);
  c->directions = (Direction*) memory_allocation(sizeof(Direction) * (seq_len - 1));
  c->length = seq_len;
}

void destroy_conform(Conformation c)
{
  free(c.energy_by_link);
  free(c.positions);
  free(c.directions);
}

void output_plot_file(int *seq, int seq_len, Conformation best_solution)
{
  int i;
  int min_x;
  int min_y;
  int max_x;
  int max_y;
  char file_name[12];

  sprintf(file_name, "%d", best_solution.energy);

  FILE *f = fopen(file_name, "w");

  if (f == NULL)
  {
    printf("Error in function output_plot_file: Unable to open output file");
    exit(1);
  }

  min_x = best_solution.positions[0].x;
  min_y = best_solution.positions[0].y;
  max_x = best_solution.positions[0].x;
  max_y = best_solution.positions[0].y;

  for (i = 1; i < seq_len; ++i)
  {
    if (best_solution.positions[i].x < min_x)
    {
      min_x = best_solution.positions[i].x;
    }
    if (best_solution.positions[i].y < min_y)
    {
      min_y = best_solution.positions[i].y;
    }
    if (best_solution.positions[i].x > max_x)
    {
      max_x = best_solution.positions[i].x;
    }
    if (best_solution.positions[i].y > max_y)
    {
      max_y = best_solution.positions[i].y;
    }
  }

  fprintf(f, "CONFORMATION DIMENSION: %d, %d\n", max_x - min_x, max_y - min_y);
  fprintf(f, "NUMBER OF AMINO_ACIDS: %d\n", seq_len);

  for (i = 0; i < seq_len; ++i)
  {
    fprintf(f, "%d %d %d\n", seq[i], best_solution.positions[i].x - min_x, best_solution.positions[i].y - min_y);
  }
  fclose(f);

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
  int heuristic_value = 0;
  int right_neighbor = lattice[pos.x + 1][pos.y];
  int left_neighbor = lattice[pos.x - 1][pos.y];
  int down_neighbor = lattice[pos.x][pos.y - 1];
  int up_neighbor = lattice[pos.x][pos.y + 1];

  if (right_neighbor >= 0 && abs(right_neighbor - nr_seq) > 1 && seq[right_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (left_neighbor >= 0 && abs(left_neighbor - nr_seq) > 1 && seq[left_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (up_neighbor >= 0 && abs(up_neighbor - nr_seq) > 1 && seq[up_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (down_neighbor >= 0 && abs(down_neighbor - nr_seq) > 1 && seq[down_neighbor] == 1)
  {
    ++heuristic_value;
  }

  return heuristic_value;
}

int random_select(double *probabilities, int len)
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

void pheromone_deposit(double **pheromone, Conformation conformation, int seq_len, int best_energy)
{
  int i;
  int j;

  for (i = 0; i < seq_len - 1; ++i)
  {
    for (j = 0; j < 3; ++j)
    {
      if (conformation.directions[i] == j && best_energy > 0)
      {
        pheromone[i][j] += (double) conformation.energy / pow(best_energy, 3);
        break;
      }
    }
  }
}

void pheromone_evaporation(double **pheromone, int seq_len, double evaporation_rate)
{
  int i;
  int j;

  for (i = 0; i < seq_len - 1; ++i)
  {
    for (j = 0; j < 3; ++j)
    {
      pheromone[i][j] *= evaporation_rate;
    }
  }
}

int construct_conform(int *seq, int seq_len, Conformation solution, Aco_config aco_config, double **pheromone, int **lattice, int *best_partial_indexes, int ant_index, Conformation *ants_conformation, Coord *blocked_positions, int num_blocked_positions)
{
  /**DECLARATIONS*/

  int i;
  int j;
  int num_candidates;
  int selected_candidate;
  double sum_probabilities;
  double probabilities[3];
  Coord curr_position;
  Coord move;
  Coord candidate_move;

  struct {
    int heuristic;
    Coord position;
    Coord move;
    Direction direction;
  } candidates[3];

  solution.energy = 0;
  solution.energy_by_link[0] = 0;


  /**DEFINES FIRST LINK*/

  /*Define first amino-acid position*/
  curr_position.x = seq_len;
  curr_position.y = seq_len;
  lattice[curr_position.x][curr_position.y] = 0;
  solution.positions[0] = curr_position;

  /*Define second amino-acid position*/
  move.x = 0;
  move.y = 1;
  solution.directions[0] = STRAIGHT;
  curr_position.x += move.x;
  curr_position.y += move.y;
  lattice[curr_position.x][curr_position.y] = 1;
  solution.positions[1] = curr_position;


  /**CONSTRUCTOR LOOP*/

  /*For each amino-acid link, except the first*/
  for (i = 1; i < seq_len - 1; ++i)
  {
    sum_probabilities = 0;
    num_candidates = 0;


    /**DEFINES CANDIDATES DIRECTION*/

    /*For each direction*/
    for (j = 0; j < 3; ++j)
    {
      candidate_move = calculate_move_by_direction(move, j);

      /*If the next position in this direction is not occupied, turns current direction into a candidate*/
      if (lattice[curr_position.x + candidate_move.x][curr_position.y + candidate_move.y] == -1)
      {
        candidates[num_candidates].move = candidate_move;
        candidates[num_candidates].direction = j;
        candidates[num_candidates].position.x = curr_position.x + candidates[num_candidates].move.x;
        candidates[num_candidates].position.y = curr_position.y + candidates[num_candidates].move.y;

        if (seq[i + 1] == 1)
        {
          candidates[num_candidates].heuristic =
            calculate_absolute_heuristic_value(lattice, i + 1, candidates[num_candidates].position, seq);
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

    /*calculate probabilities*/
    for (j = 0; j < num_candidates; ++j)
    {
      probabilities[j] = probabilities[j]/sum_probabilities;
    }


    /**SELECTS A CANDIDATE*/

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


    /**UPDATE CONFORMATION*/

    if (selected_candidate != -1)
    {
      move = candidates[selected_candidate].move;
      curr_position = candidates[selected_candidate].position;
      solution.directions[i] = candidates[selected_candidate].direction;
      solution.energy_by_link[i] = solution.energy;
      solution.energy -= candidates[selected_candidate].heuristic;
      solution.positions[i + 1] = curr_position;
      lattice[curr_position.x][curr_position.y] = i + 1;

    }
    else
    {
      /**WHEN IS IMPOSSIBLE CONTINUE THE FOLD PROCESS*/

      switch (aco_config.unfeasible_conformation_handler)
      {
        case BLOCKED_POSITIONS:

          ++num_blocked_positions;
          lattice[curr_position.x][curr_position.y] = -2;
          if (seq[i - 1] == 1)
          {
            solution.energy += calculate_absolute_heuristic_value(lattice, i - 1, curr_position, seq);
          }
          curr_position = solution.positions[i - 1];
          move = subtract_coord(curr_position, solution.positions[i - 2]);
          i -= 2;
          break;

        case PARTIAL_COPY:

          /*If theres no another conformation to copy*/
          if (best_partial_indexes[i] == -1)
          {
            /*Cleans lattice*/
            for (j = 0; j <= i; ++j)
            {
              lattice[solution.positions[j].x][solution.positions[j].y] = -1;
            }

            i = 0;
            solution.energy = 0;

            /*Defines first link*/
            curr_position.x = seq_len;
            curr_position.y = seq_len;
            lattice[curr_position.x][curr_position.y] = 0;
            solution.directions[0] = STRAIGHT;
            solution.energy_by_link[0] = 0;
            solution.positions[0] = curr_position;
            move.x = 0;
            move.y = 1;
            curr_position.x += move.x;
            curr_position.y += move.y;
            lattice[curr_position.x][curr_position.y] = 1;
            solution.positions[1] = curr_position;
          }
          else
          {
            /*Cleans lattice until i th amino-acid*/
            for (j = 0; j <= i; ++j)
            {
              lattice[solution.positions[j].x][solution.positions[j].y] = -1;
            }

            /*Copy best_ant_solution until i th amino-acid*/
            for (j = 0; j <= i; ++j)
            {
              lattice[ants_conformation[best_partial_indexes[i]].positions[j].x][ants_conformation[best_partial_indexes[i]].positions[j].y] = j;
              solution.positions[j] = ants_conformation[best_partial_indexes[i]].positions[j];
              if (j != i)
              {
                solution.directions[j] = ants_conformation[best_partial_indexes[i]].directions[j];
                solution.energy_by_link[j] = ants_conformation[best_partial_indexes[i]].energy_by_link[j];
              }
            }
            curr_position = ants_conformation[best_partial_indexes[i]].positions[i];
            solution.energy = ants_conformation[best_partial_indexes[i]].energy_by_link[i];
            move.x = curr_position.x - ants_conformation[best_partial_indexes[i]].positions[i - 1].x;
            move.y = curr_position.y - ants_conformation[best_partial_indexes[i]].positions[i - 1].y;
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
    if (best_partial_indexes[i] == -1 || ants_conformation[best_partial_indexes[i]].energy_by_link[i] > solution.energy_by_link[i])
    {
      best_partial_indexes[i] = ant_index;
    }
  }

  return solution.energy;
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

  if (right_neighbor >= 0 && right_neighbor < nr_seq - 1 && seq[right_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (left_neighbor >= 0 && left_neighbor < nr_seq - 1 && seq[left_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (up_neighbor >= 0 && up_neighbor < nr_seq - 1 && seq[up_neighbor] == 1)
  {
    ++heuristic_value;
  }
  if (down_neighbor >= 0 && down_neighbor < nr_seq - 1 && seq[down_neighbor] == 1)
  {
    ++heuristic_value;
  }

  return heuristic_value;
}

int generate_pull_move_configs(int nr_seq, int **lattice, Coord *ant_positions, Pull_move_config *configs, int config_index, int seq_len)
{
  int num_configs;
  Pull_move_config config;

  num_configs = 0;
  config.curr = ant_positions[nr_seq];
  config.next = ant_positions[nr_seq + 1];
  config.prev = ant_positions[nr_seq - 1];
  config.nr_seq = nr_seq;

  /*If is adjacent (right) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x + 1][config.next.y] == -1 && config.next.x + 1 != config.curr.x &&
      config.next.y != config.curr.y)
  {
    config.f.x = config.next.x + 1;
    config.f.y = config.next.y;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;

      ++config_index;
      ++num_configs;
    }

  }

  /*If is adjacent (left) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x - 1][config.next.y] < 0 && config.next.x - 1 != config.curr.x &&
      config.next.y != config.curr.y)
  {
    config.f.x = config.next.x - 1;
    config.f.y = config.next.y;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  /*If is adjacent (up) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x][config.next.y + 1] < 0 && config.next.x != config.curr.x &&
      config.next.y + 1 != config.curr.y)
  {
    config.f.x = config.next.x;
    config.f.y = config.next.y + 1;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;

      ++config_index;
      ++num_configs;
    }
  }

  /*If is adjacent (down) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x][config.next.y - 1] < 0 && config.next.x != config.curr.x &&
      config.next.y - 1 != config.curr.y)
  {
    config.f.x = config.next.x;
    config.f.y = config.next.y - 1;

    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  return num_configs;
}

int apply_pull_move(Conformation conformation, Pull_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions)
{

  int i;
  int temp;
  int last_modified; /*Modifications in the conformation finalizes in this amino-acid*/

  if (config.c.x == config.prev.x && config.c.y == config.prev.y)
  {
    /**WHEN C == PREV: CURRENT <-> F*/

    /*Update energy*/
    if (seq[config.nr_seq] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
    }

    /*Update lattice*/
    temp = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
    lattice[config.f.x][config.f.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq].x = config.f.x;
    conformation.positions[config.nr_seq].y = config.f.y;
    config.curr.x =  conformation.positions[config.nr_seq].x;
    config.curr.y =  conformation.positions[config.nr_seq].y;

    /*Update directions*/
    if (config.nr_seq == 1)
    {
      Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
      conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
    }
    else
    {
      Coord subtraction1_result = subtract_coord(config.prev, conformation.positions[config.nr_seq - 2]);
      Coord subtraction2_result = subtract_coord(config.curr, config.prev);
      conformation.directions[config.nr_seq - 1] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
    }

    Coord subtraction1_result = subtract_coord(config.curr, config.prev);
    Coord subtraction2_result = subtract_coord(config.next, config.curr);
    conformation.directions[config.nr_seq] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
    last_modified = config.nr_seq;
  }
  else
  {

    /**CASO EM QUE C == PREVIOUS. PREVIOUS <-> C*/

    /*Update energy*/
    if (seq[config.nr_seq - 1] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, config.prev, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, config.c, seq);
    }

    /*Update lattice*/
    temp = lattice[config.prev.x][config.prev.y];
    lattice[config.prev.x][config.prev.y] = lattice[config.c.x][config.c.y];
    lattice[config.c.x][config.c.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq - 1].x = config.c.x;
    conformation.positions[config.nr_seq - 1].y = config.c.y;

    /**CURRENT <-> F*/

    if (seq[config.nr_seq] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
    }

    /*Update lattice*/
    temp = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
    lattice[config.f.x][config.f.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq].x = config.f.x;
    conformation.positions[config.nr_seq].y = config.f.y;

    i = config.nr_seq - 2; /*Amino-acid that comes before the previous*/
    last_modified = 0;     /*If all amino-acids before the previous needs be moved, the last modified will be the first*/

    /*Adjusts the conformation until reach a valid configuration*/
    while (i >= 0)
    {

      /**I <- CURRENT <- PREVIOUS*/

      /*Update energy*/
      if (seq[i] == 1)
      {
        conformation.energy += calculate_absolute_heuristic_value(lattice, i, conformation.positions[i], seq);
        conformation.energy -= calculate_absolute_heuristic_value(lattice, i, config.curr, seq);
      }

      /*Update lattice*/
      temp = lattice[conformation.positions[i].x][conformation.positions[i].y];
      lattice[conformation.positions[i].x][conformation.positions[i].y] = lattice[config.curr.x][config.curr.y];
      lattice[config.curr.x][config.curr.y] = temp;

      /*Update positions*/
      temp = conformation.positions[i].x;
      conformation.positions[i].x = config.curr.x;
      config.curr.x = config.prev.x;
      config.prev.x = temp;
      temp = conformation.positions[i].y;
      conformation.positions[i].y = config.curr.y;
      config.curr.y = config.prev.y;
      config.prev.y = temp;

      /*If reach a valid conformation, stop. Then the last modified will be de i th amino-acid*/
      if (i != 0 && coords_distance(conformation.positions[i], conformation.positions[i - 1]) == 1)
      {
        last_modified = i;
        break;
      }

      --i;
    }

    /*Update directions*/
    if (i == -1)
    {
      Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
      conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
      i = last_modified + 1;
    }
    while (i < config.nr_seq + 1)
    {
      Coord subtraction1_result = subtract_coord(conformation.positions[i], conformation.positions[i - 1]);
      Coord subtraction2_result = subtract_coord(conformation.positions[i + 1], conformation.positions[i]);
      conformation.directions[i] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
      ++i;
    }
  }

  /**DEVOLVES LATTICE TO ORIGINAL STATE*/
  for (i = last_modified; i < config.nr_seq + 1; ++i)
  {
    lattice[conformation.positions[i].x][conformation.positions[i].y] = -1;
  }
  for (i = last_modified; i < config.nr_seq + 1; ++i)
  {
    lattice[ant_positions[i].x][ant_positions[i].y] = i;
  }

  return conformation.energy;
}

int generate_pull_move_configs_inverse(int nr_seq, int **lattice, Coord *original_positions, Pull_move_config *configs, int config_index, int seq_len)
{
  int num_configs;
  Pull_move_config config;

  num_configs = 0;
  config.curr = original_positions[nr_seq];
  config.next = original_positions[nr_seq - 1];
  config.prev = original_positions[nr_seq + 1];
  config.nr_seq = nr_seq;

  /*If is adjacent (right) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x + 1][config.next.y] < 0 && config.next.x + 1 != config.curr.x &&
      config.next.y != config.curr.y)
  {
    config.f.x = config.next.x + 1;
    config.f.y = config.next.y;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  /*If is adjacent (left) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x - 1][config.next.y] < 0 && config.next.x - 1 != config.curr.x &&
      config.next.y != config.curr.y)
  {
    config.f.x = config.next.x - 1;
    config.f.y = config.next.y;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  /*If is adjacent (up) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x][config.next.y + 1] < 0 && config.next.x != config.curr.x &&
      config.next.y + 1 != config.curr.y)
  {
    config.f.x = config.next.x;
    config.f.y = config.next.y + 1;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  /*If is adjacent (down) to next amino-acid and diagonally adjacent to current amino-acid*/
  if (lattice[config.next.x][config.next.y - 1] < 0 && config.next.x != config.curr.x &&
      config.next.y - 1 != config.curr.y)
  {
    config.f.x = config.next.x;
    config.f.y = config.next.y - 1;
    config.c.x = config.f.x + config.curr.x - config.next.x;
    config.c.y = config.f.y + config.curr.y - config.next.y;

    /*If F exists and C is empty or is equals do previous amino-acid*/
    if (lattice[config.c.x][config.c.y] < 0 || (config.c.x == config.prev.x && config.c.y == config.prev.y))
    {
      configs[config_index].c = config.c;
      configs[config_index].curr = config.curr;
      configs[config_index].f = config.f;
      configs[config_index].next = config.next;
      configs[config_index].nr_seq = config.nr_seq;
      configs[config_index].prev = config.prev;
      ++config_index;
      ++num_configs;
    }
  }

  return num_configs;
}

int apply_pull_move_inverse(Conformation conformation, Pull_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions)
{
  int i;
  int temp;

  if (config.c.x == config.prev.x && config.c.y == config.prev.y)
  {
    /**CASO EM QUE C == PREV. CURRENT <-> F*/

    /*Update energy*/
    if (seq[config.nr_seq] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
    }

    /*Update lattice*/
    temp = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
    lattice[config.f.x][config.f.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq].x = config.f.x;
    conformation.positions[config.nr_seq].y = config.f.y;
    config.curr.x =  conformation.positions[config.nr_seq].x;
    config.curr.y =  conformation.positions[config.nr_seq].y;

    /*Update directions*/
    if (config.nr_seq == 1)
    {
      Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
      conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
    }
    else
    {
      Coord subtraction1_result = subtract_coord(config.next, conformation.positions[config.nr_seq - 2]);
      Coord subtraction2_result = subtract_coord(config.curr, config.next);
      conformation.directions[config.nr_seq - 1] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
    }

    Coord subtraction1_result = subtract_coord(config.curr, config.next);
    Coord subtraction2_result = subtract_coord(config.prev, config.curr);
    conformation.directions[config.nr_seq] = calculate_direction_by_move(subtraction1_result, subtraction2_result);

  }
  else
  {
    /**CASO EM QUE C == PREVIOUS. PREVIOUS <-> C*/

    /*Update energy*/
    if (seq[config.nr_seq + 1] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq + 1, config.prev, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq + 1, config.c, seq);
    }

    /*Update lattice*/
    temp = lattice[config.prev.x][config.prev.y];
    lattice[config.prev.x][config.prev.y] = lattice[config.c.x][config.c.y];
    lattice[config.c.x][config.c.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq + 1].x = config.c.x;
    conformation.positions[config.nr_seq + 1].y = config.c.y;

    /**CURRENT <-> F*/

    if (seq[config.nr_seq] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.f, seq);
    }

    /*Update lattice*/
    temp = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = lattice[config.f.x][config.f.y];
    lattice[config.f.x][config.f.y] = temp;

    /*Update positions*/
    conformation.positions[config.nr_seq].x = config.f.x;
    conformation.positions[config.nr_seq].y = config.f.y;

    i = config.nr_seq + 2;/*Amino-acid that comes before the previous*/

    /*Adjusts the conformation until reach a valid configuration*/
    while (i < seq_len)
    {
      /**I <- CURRENT <- PREVIOUS*/

      /*Update energy*/
      if (seq[i] == 1)
      {
        conformation.energy += calculate_absolute_heuristic_value(lattice, i, conformation.positions[i], seq);
        conformation.energy -= calculate_absolute_heuristic_value(lattice, i, config.curr, seq);
      }

      /*Update lattice*/
      temp = lattice[conformation.positions[i].x][conformation.positions[i].y];
      lattice[conformation.positions[i].x][conformation.positions[i].y] = lattice[config.curr.x][config.curr.y];
      lattice[config.curr.x][config.curr.y] = temp;

      /*Update positions*/
      temp = conformation.positions[i].x;
      conformation.positions[i].x = config.curr.x;
      config.curr.x = config.prev.x;
      config.prev.x = temp;
      temp = conformation.positions[i].y;
      conformation.positions[i].y = config.curr.y;
      config.curr.y = config.prev.y;
      config.prev.y = temp;

      /*If reach a valid conformation, stop. Then the last modified will be de i th amino-acid*/
      if (i != seq_len - 1 && coords_distance(conformation.positions[i], conformation.positions[i + 1]) == 1)
      {
        break;
      }

      ++i;
    }

    /*Update directions*/
    if (config.nr_seq > 1)
    {
      i = config.nr_seq - 1;
    }
    else
    {
      Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
      conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
      i = 1;
    }

    while (i < seq_len - 1)
    {
      Coord subtraction1_result = subtract_coord(conformation.positions[i], conformation.positions[i - 1]);
      Coord subtraction2_result = subtract_coord(conformation.positions[i + 1], conformation.positions[i]);
      conformation.directions[i] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
      ++i;
    }
  }

  /**DEVOLVES LATTICE TO ORIGINAL STATE*/
  for (i = config.nr_seq; i < seq_len; ++i)
  {
    lattice[conformation.positions[i].x][conformation.positions[i].y] = -1;
  }
  for (i = config.nr_seq; i < seq_len; ++i)
  {
    lattice[ant_positions[i].x][ant_positions[i].y] = i;
  }

  return conformation.energy;
}

int generate_end_move_configs(int nr_seq, int **lattice, Coord *original_positions, End_move_config *configs, int seq_len)
{
  int num_configs;
  End_move_config config;
  Coord adj;
  int config_index = 0;

  Coord up = {original_positions[nr_seq].x, original_positions[nr_seq].y + 1};
  Coord down = {original_positions[nr_seq].x, original_positions[nr_seq].y - 1};
  Coord right = {original_positions[nr_seq].x + 1, original_positions[nr_seq].y};
  Coord left = {original_positions[nr_seq].x - 1, original_positions[nr_seq].y};
  Coord curr = {original_positions[nr_seq].x, original_positions[nr_seq].y};

  num_configs = 0;
  config.nr_seq = nr_seq;
  config.curr = curr;

  if (lattice[up.x][up.y] < 0 && (up.x == curr.x || up.y == curr.y) && (abs(up.x - curr.x) == 1 || abs(up.y - curr.y) == 1))
  {
    adj = up;

    up.y = up.y + 1;

    right.x = up.x + 1;
    right.y = up.y;

    left.x = up.x - 1;
    left.y = up.y;

    if (lattice[up.x][up.y] < 0 && (up.x == adj.x || up.y == adj.y) && (abs(up.x - adj.x) == 1 || abs(up.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = up;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[right.x][right.y] < 0 && (right.x == adj.x || right.y == adj.y) && (abs(right.x - adj.x) == 1 || abs(right.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = right;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[left.x][left.y] < 0 && (left.x == adj.x || left.y == adj.y) && (abs(left.x - adj.x) == 1 || abs(left.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = left;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }

    up.x = original_positions[nr_seq].x;
    right.x = original_positions[nr_seq].x + 1;
    left.x = original_positions[nr_seq].x - 1;

    up.y = original_positions[nr_seq].y + 1;
    right.y = original_positions[nr_seq].y;
    left.y = original_positions[nr_seq].y;
  }

  if (lattice[down.x][down.y] < 0 && (down.x == curr.x || down.y == curr.y) && (abs(down.x - curr.x) == 1 || abs(down.y - curr.y) == 1))
  {
    adj = down;

    down.y = down.y - 1;

    right.x = down.x + 1;
    right.y = down.y;

    left.x = down.x - 1;
    left.y = down.y;

    if (lattice[down.x][down.y] < 0 && (down.x == adj.x || down.y == adj.y) && (abs(down.x - adj.x) == 1 || abs(down.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = down;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[right.x][right.y] < 0 && (right.x == adj.x || right.y == adj.y) && (abs(right.x - adj.x) == 1 || abs(right.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = right;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[left.x][left.y] < 0 && (left.x == adj.x || left.y == adj.y) && (abs(left.x - adj.x) == 1 || abs(left.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = left;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }

    down.x = original_positions[nr_seq].x;
    right.x = original_positions[nr_seq].x + 1;
    left.x = original_positions[nr_seq].x - 1;

    down.y = original_positions[nr_seq].y - 1;
    right.y = original_positions[nr_seq].y;
    left.y = original_positions[nr_seq].y;
  }

  if (lattice[right.x][right.y] < 0 && (right.x == curr.x || right.y == curr.y) && (abs(right.x - curr.x) == 1 || abs(right.y - curr.y) == 1))
  {
    adj = right;

    right.x = right.x + 1;

    down.x = right.x;
    down.y = right.y - 1;

    up.x = right.x;
    up.y = right.y + 1;

    if (lattice[down.x][down.y] < 0 && (down.x == adj.x || down.y == adj.y) && (abs(down.x - adj.x) == 1 || abs(down.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = down;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[right.x][right.y] < 0 && (right.x == adj.x || right.y == adj.y) && (abs(right.x - adj.x) == 1 || abs(right.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = right;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[up.x][up.y] < 0 && (up.x == adj.x || up.y == adj.y) && (abs(up.x - adj.x) == 1 || abs(up.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = up;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }

    up.x = original_positions[nr_seq].x;
    down.x = original_positions[nr_seq].x;
    right.x = original_positions[nr_seq].x + 1;

    up.y = original_positions[nr_seq].y + 1;
    down.y = original_positions[nr_seq].y - 1;
    right.y = original_positions[nr_seq].y;
  }


  if (lattice[left.x][left.y] < 0 && (left.x == curr.x || left.y == curr.y) && (abs(left.x - curr.x) == 1 || abs(left.y - curr.y) == 1))
  {
    adj = left;

    left.x = left.x - 1;

    down.x = left.x;
    down.y = left.y - 1;

    up.x = left.x;
    up.y = left.y + 1;

    if (lattice[down.x][down.y] < 0 && (down.x == adj.x || down.y == adj.y) && (abs(down.x - adj.x) == 1 || abs(down.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = down;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[up.x][up.y] < 0 && (up.x == adj.x || up.y == adj.y) && (abs(up.x - adj.x) == 1 || abs(up.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = right;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }
    if (lattice[down.x][down.y] < 0 && (down.x == adj.x || down.y == adj.y) && (abs(down.x - adj.x) == 1 || abs(down.y - adj.y) == 1))
    {
      config.adj = adj;
      config.away = down;
      configs[config_index] = config;
      ++config_index;
      ++num_configs;
    }

    up.x = original_positions[nr_seq].x;
    down.x = original_positions[nr_seq].x;
    left.x = original_positions[nr_seq].x - 1;

    up.y = original_positions[nr_seq].y + 1;
    down.y = original_positions[nr_seq].y - 1;
    left.y = original_positions[nr_seq].y;
  }

  return num_configs;
}

int apply_end_move(Conformation conformation, End_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions)
{

  int i;
  int temp;
  int last_modified; /*Modifications in the conformation finalizes in this amino-acid*/
  Coord prev;

  /**CURR <-> AWAY*/

  /*Update energy*/
  if (seq[config.nr_seq] == 1)
  {
    conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
    conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.away, seq);
  }

  /*Update lattice*/
  temp = lattice[config.curr.x][config.curr.y];
  lattice[config.curr.x][config.curr.y] = lattice[config.away.x][config.away.y];
  lattice[config.away.x][config.away.y] = temp;

  /*Update positions*/
  conformation.positions[config.nr_seq].x = config.away.x;
  conformation.positions[config.nr_seq].y = config.away.y;

  /**PREV <-> ADJ*/

  prev = conformation.positions[config.nr_seq - 1];

  if (seq[config.nr_seq - 1] == 1)
  {
    conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, prev, seq);
    conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq - 1, config.adj, seq);
  }

  /*Update lattice*/
  temp = lattice[prev.x][prev.y];
  lattice[prev.x][prev.y] = lattice[config.adj.x][config.adj.y];
  lattice[config.adj.x][config.adj.y] = temp;

  /*Update positions*/
  conformation.positions[config.nr_seq - 1].x = config.adj.x;
  conformation.positions[config.nr_seq - 1].y = config.adj.y;

  i = config.nr_seq - 2; /*Amino-acid that comes before the previous*/
  last_modified = 0;     /*If all amino-acids before the previous needs be moved, the last modified will be the first*/

  /*Adjusts the conformation until reach a valid configuration*/
  while (i >= 0)
  {

    /**I <- CURRENT <- PREVIOUS*/

    /*Update energy*/
    if (seq[i] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, i, conformation.positions[i], seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, i, config.curr, seq);
    }

    /*Update lattice*/
    temp = lattice[conformation.positions[i].x][conformation.positions[i].y];
    lattice[conformation.positions[i].x][conformation.positions[i].y] = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = temp;

    /*Update positions*/
    temp = conformation.positions[i].x;
    conformation.positions[i].x = config.curr.x;
    config.curr.x = prev.x;
    prev.x = temp;

    temp = conformation.positions[i].y;
    conformation.positions[i].y = config.curr.y;
    config.curr.y = prev.y;
    prev.y = temp;

    /*If reach a valid conformation, stop. Then the last modified will be de i th amino-acid*/
    if (i != 0 && coords_distance(conformation.positions[i], conformation.positions[i - 1]) == 1)
    {
      last_modified = i;
      break;
    }

    --i;
  }

  /*Update directions*/
  if (i == -1)
  {
    Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
    conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
    i = last_modified + 1;
  }
  while (i < config.nr_seq)
  {
    Coord subtraction1_result = subtract_coord(conformation.positions[i], conformation.positions[i - 1]);
    if (i > seq_len - 2){
      printf("\n");
    }
    Coord subtraction2_result = subtract_coord(conformation.positions[i + 1], conformation.positions[i]);
    conformation.directions[i] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
    ++i;
  }


  /**DEVOLVES LATTICE TO ORIGINAL STATE*/
  for (i = last_modified; i < config.nr_seq + 1; ++i)
  {
    lattice[conformation.positions[i].x][conformation.positions[i].y] = -1;
  }
  for (i = last_modified; i < config.nr_seq + 1; ++i)
  {
    lattice[ant_positions[i].x][ant_positions[i].y] = i;
  }

  return conformation.energy;
}


int apply_end_move_inverse(Conformation conformation, End_move_config config, int *seq, int **lattice, int seq_len, Coord *ant_positions)
{

  int i;
  int temp;
  Coord prev;

  /**CURR <-> AWAY*/

  /*Update energy*/
  if (seq[config.nr_seq] == 1)
  {
    conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq, config.curr, seq);
    conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq, config.away, seq);
  }

  /*Update lattice*/
  temp = lattice[config.curr.x][config.curr.y];
  lattice[config.curr.x][config.curr.y] = lattice[config.away.x][config.away.y];
  lattice[config.away.x][config.away.y] = temp;

  /*Update positions*/
  conformation.positions[config.nr_seq].x = config.away.x;
  conformation.positions[config.nr_seq].y = config.away.y;

  /**PREV <-> ADJ*/

  prev = conformation.positions[config.nr_seq + 1];

  if (seq[config.nr_seq + 1] == 1)
  {
    conformation.energy += calculate_absolute_heuristic_value(lattice, config.nr_seq + 1, prev, seq);
    conformation.energy -= calculate_absolute_heuristic_value(lattice, config.nr_seq + 1, config.adj, seq);
  }

  /*Update lattice*/
  temp = lattice[prev.x][prev.y];
  lattice[prev.x][prev.y] = lattice[config.adj.x][config.adj.y];
  lattice[config.adj.x][config.adj.y] = temp;

  /*Update positions*/
  conformation.positions[config.nr_seq + 1].x = config.adj.x;
  conformation.positions[config.nr_seq + 1].y = config.adj.y;

  i = config.nr_seq + 2;/*Amino-acid that comes before the previous*/

  /*Adjusts the conformation until reach a valid configuration*/
  while (i < seq_len)
  {
    /**I <- CURRENT <- PREVIOUS*/

    /*Update energy*/
    if (seq[i] == 1)
    {
      conformation.energy += calculate_absolute_heuristic_value(lattice, i, conformation.positions[i], seq);
      conformation.energy -= calculate_absolute_heuristic_value(lattice, i, config.curr, seq);
    }

    /*Update lattice*/
    temp = lattice[conformation.positions[i].x][conformation.positions[i].y];
    lattice[conformation.positions[i].x][conformation.positions[i].y] = lattice[config.curr.x][config.curr.y];
    lattice[config.curr.x][config.curr.y] = temp;

    /*Update positions*/
    temp = conformation.positions[i].x;
    conformation.positions[i].x = config.curr.x;
    config.curr.x = prev.x;
    prev.x = temp;
    temp = conformation.positions[i].y;
    conformation.positions[i].y = config.curr.y;
    config.curr.y = prev.y;
    prev.y = temp;

    /*If reach a valid conformation, stop. Then the last modified will be de i th amino-acid*/
    if (i != seq_len - 1 && coords_distance(conformation.positions[i], conformation.positions[i + 1]) == 1)
    {
      break;
    }

    ++i;
  }

  /*Update directions*/
  if (config.nr_seq > 1)
  {
    i = config.nr_seq - 1;
  }
  else
  {
    Coord subtraction_result = subtract_coord(conformation.positions[1], conformation.positions[0]);
    conformation.directions[0] = calculate_absolute_direction_by_move(subtraction_result);
    i = 1;
  }

  while (i < seq_len - 1)
  {
    Coord subtraction1_result = subtract_coord(conformation.positions[i], conformation.positions[i - 1]);
    Coord subtraction2_result = subtract_coord(conformation.positions[i + 1], conformation.positions[i]);
    conformation.directions[i] = calculate_direction_by_move(subtraction1_result, subtraction2_result);
    ++i;
  }


  /**DEVOLVES LATTICE TO ORIGINAL STATE*/
  for (i = config.nr_seq; i < seq_len; ++i)
  {
    lattice[conformation.positions[i].x][conformation.positions[i].y] = -1;
  }
  for (i = config.nr_seq; i < seq_len; ++i)
  {
    lattice[ant_positions[i].x][ant_positions[i].y] = i;
  }

  return conformation.energy;
}

int calculate_best_pull_move(int *seq, int seq_len, Conformation ant_conformation,int **lattice, Conformation best_conformation, Conformation config_conformation, Pull_move_config *configs, Pull_move_config *configs_inverse)
{
  int i;
  int j;
  int num_configs = 0;
  int num_configs_inverse = 0;
  int previous_energy;
  End_move_config end_move_configs[9];

  best_conformation.energy = 0;

  /**APPLY PULL-MOVE WHILE THERE IS IMPORVEMENT*/

  do {

    num_configs = 0;
    num_configs_inverse = 0;
    previous_energy = ant_conformation.energy;

    for (i = 0; i < seq_len; ++i)
    {
      lattice[ant_conformation.positions[i].x][ant_conformation.positions[i].y] = i;
    }

    /**GENERATE ALL POSSIBLE PULL-MOVES CONFIGURATIONS*/

    for (i = seq_len - 2; i > 0; --i)
    {
      num_configs += generate_pull_move_configs(i, lattice, ant_conformation.positions, configs, num_configs, seq_len);
    }
    for (i = 1; i < seq_len - 1; ++i)
    {
      num_configs_inverse += generate_pull_move_configs_inverse(i, lattice, ant_conformation.positions,
                                                                configs_inverse, num_configs_inverse, seq_len);
    }

    /*Apply all possible inverse pull-moves and save the best*/
    for (i = 0; i < num_configs_inverse; ++i)
    {
      /*Copy original conformation*/
      config_conformation.energy = ant_conformation.energy;
      for (j = 0; j < seq_len; ++j)
      {
        config_conformation.positions[j] = ant_conformation.positions[j];
        if (j < seq_len - 1)
        {
          config_conformation.directions[j] = ant_conformation.directions[j];
        }
      }

      /*Apply pull move on the copy*/
      config_conformation.energy = apply_pull_move_inverse(config_conformation, configs_inverse[i], seq,
                                                           lattice, seq_len, ant_conformation.positions);

      /*If the energy of the new conformation is lower than the actual best, update the best*/
      if (config_conformation.energy < best_conformation.energy)
      {
        best_conformation.energy = config_conformation.energy;

        for (j = 0; j < seq_len; ++j)
        {
          best_conformation.positions[j] = config_conformation.positions[j];
          if (j < seq_len - 1)
          {
            best_conformation.directions[j] = config_conformation.directions[j];
          }
        }
      }
    }

    /*Apply all possible pull-moves and save the best*/
    for (i = 0; i < num_configs; ++i)
    {
      /*Copy original conformation*/
      config_conformation.energy = ant_conformation.energy;
      for (j = 0; j < seq_len; ++j)
      {
        config_conformation.positions[j] = ant_conformation.positions[j];
        if (j < seq_len - 1)
        {
          config_conformation.directions[j] = ant_conformation.directions[j];
        }
      }

      /*Apply pull move on the copy*/
      config_conformation.energy = apply_pull_move(config_conformation, configs[i], seq,
                                                   lattice, seq_len, ant_conformation.positions);

      /*If the energy of the new conformation is lower than the actual best, update the best*/
      if (config_conformation.energy < best_conformation.energy)
      {
        best_conformation.energy = config_conformation.energy;

        for (j = 0; j < seq_len; ++j)
        {
          best_conformation.positions[j] = config_conformation.positions[j];
          if (j < seq_len - 1)
          {
            best_conformation.directions[j] = config_conformation.directions[j];
          }
        }
      }
    }

    /**END-MOVES*/

    num_configs = generate_end_move_configs(seq_len - 1, lattice, ant_conformation.positions, end_move_configs, seq_len);

    /*Apply all possible end-moves and save the best*/
    for (i = 0; i < num_configs; ++i)
    {
      /*Copy original conformation*/
      config_conformation.energy = ant_conformation.energy;
      for (j = 0; j < seq_len; ++j)
      {
        config_conformation.positions[j] = ant_conformation.positions[j];
        if (j < seq_len - 1)
        {
          config_conformation.directions[j] = ant_conformation.directions[j];
        }
      }

      /*Apply pull move on the copy*/
      config_conformation.energy = apply_end_move(config_conformation, end_move_configs[i], seq,
                                                  lattice, seq_len, ant_conformation.positions);

      /*If the energy of the new conformation is lower than the actual best, update the best*/
      if (config_conformation.energy < best_conformation.energy)
      {
        best_conformation.energy = config_conformation.energy;

        for (j = 0; j < seq_len; ++j)
        {
          best_conformation.positions[j] = config_conformation.positions[j];
          if (j < seq_len - 1)
          {
            best_conformation.directions[j] = config_conformation.directions[j];
          }
        }
      }
    }

    /**INVERSE END-MOVES*/

    num_configs = generate_end_move_configs(0, lattice, ant_conformation.positions, end_move_configs, seq_len);

    /*Apply all possible end-moves and save the best*/
    for (i = 0; i < num_configs; ++i)
    {
      /*Copy original conformation*/
      config_conformation.energy = ant_conformation.energy;
      for (j = 0; j < seq_len; ++j)
      {
        config_conformation.positions[j] = ant_conformation.positions[j];
        if (j < seq_len - 1)
        {
          config_conformation.directions[j] = ant_conformation.directions[j];
        }
      }

      /*Apply pull move on the copy*/
      config_conformation.energy = apply_end_move_inverse(config_conformation, end_move_configs[i], seq,
                                                          lattice, seq_len, ant_conformation.positions);

      /*If the energy of the new conformation is lower than the actual best, update the best*/
      if (config_conformation.energy < best_conformation.energy)
      {
        best_conformation.energy = config_conformation.energy;

        for (j = 0; j < seq_len; ++j)
        {
          best_conformation.positions[j] = config_conformation.positions[j];
          if (j < seq_len - 1)
          {
            best_conformation.directions[j] = config_conformation.directions[j];
          }
        }
      }
    }

    /**CHECKS IF NEW CONFORMATION IS BETTER THAN ORIGINAL CONFORMATION*/

    /*If the best found pull-moved conformation is better than the original conformation,
    update the original and return the new energy*/
    if (best_conformation.energy != 0 && best_conformation.energy <= ant_conformation.energy)
    {
      Coord distance_to_lattice_center;

      distance_to_lattice_center.x = seq_len - best_conformation.positions[0].x;
      distance_to_lattice_center.y = seq_len - best_conformation.positions[0].y;
      ant_conformation.energy = previous_energy;

      for (i = 0; i < seq_len; ++i)
      {
        /*Clean lattice*/
        lattice[ant_conformation.positions[i].x][ant_conformation.positions[i].y] = -1;
      }

      for (i = 0; i < seq_len; ++i)
      {
        /*Move protein to center*/
        best_conformation.positions[i].x += distance_to_lattice_center.x;
        best_conformation.positions[i].y += distance_to_lattice_center.y;

        ant_conformation.energy = best_conformation.energy;

        /*Update original conformation*/
        ant_conformation.positions[i] = best_conformation.positions[i];
        if (i < seq_len - 1)
        {
          ant_conformation.directions[i] = best_conformation.directions[i];
        }

      }
    }
  }
  while(best_conformation.energy < previous_energy);

  /*Clean lattice*/
  for (i = 0; i < seq_len; ++i)
  {
    lattice[ant_conformation.positions[i].x][ant_conformation.positions[i].y] = -1;
  }

  return ant_conformation.energy;
}
