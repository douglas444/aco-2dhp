#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "aco.h"
#include "file_reader.h"

void output_plot_file(int *seq, int seq_len, Conformation best_solution, char *file_name);

int main(int argc, char **argv)
{
  int i;
  int *binary_sequence;
  Aco_config aco_config;
  char *input_file;
  char *unfeasible_conformation_handler;
  char *local_search;
  int sequence_size;
  char *char_sequence;
  Conformation conformation;
  char c;
  clock_t t0;
  double time;
  char *output_file_name;

  if (argc < 3) {
    printf("Invalid parameters\n");
    exit(1);
  }

  ///Extract file content
  input_file = load_file_content(argv[1]);

  ///Get parameters from file content
  aco_config.alpha = char_to_double(get_key_value(input_file, "alpha"));
  aco_config.beta = char_to_double(get_key_value(input_file, "beta"));
  aco_config.ini_pheromone = char_to_double(get_key_value(input_file, "initial-pheromone"));
  aco_config.persistence = char_to_double(get_key_value(input_file, "persistence"));
  aco_config.iterations = char_to_int(get_key_value(input_file, "iterations"));
  aco_config.population = char_to_int(get_key_value(input_file, "population"));
  char_sequence = get_key_value(input_file, "sequence");
  unfeasible_conformation_handler = get_key_value(input_file, "unfeasible-conformation-handler");
  local_search = get_key_value(input_file, "local-search");

  ///Set local search method
  if (strcmp(local_search, "WITHOUT_LOCAL_SEARCH") == 0) {
    aco_config.local_search = WITHOUT_LOCAL_SEARCH;
  } else if (strcmp(local_search, "PULL_MOVE") == 0) {
    aco_config.local_search = PULL_MOVE;
  } else {
    printf("Error in function main: Local search parameter not recognized\n");
    exit(1);
  }

  ///Set unfeasible conformation handler method
  if (strcmp(unfeasible_conformation_handler, "PARTIAL_COPY") == 0) {
    aco_config.unfeasible_conformation_handler = PARTIAL_COPY;
  } else if (strcmp(unfeasible_conformation_handler, "BLOCKED_POSITIONS") == 0) {
    aco_config.unfeasible_conformation_handler = BLOCKED_POSITIONS;
  } else {
    printf("Error in function main: Unfeasible conformation handler parameter not recognized\n");
    exit(1);
  }

  ///generate protein binary sequence
  sequence_size = strlen(char_sequence);
  binary_sequence = (int*) malloc(sequence_size * sizeof(int));
  if (binary_sequence == NULL)
  {
    printf("Error in function main: Unable to allocate memory");
    exit(1);
  }
  for (i = 0; i < sequence_size; ++i)
  {
    binary_sequence[i] = char_sequence[i] == 'H' ? 1 : 0;
  }

  ///Run ACO
  t0 = clock();
  conformation = aco_run(binary_sequence, sequence_size, aco_config);
  time = (clock() - t0)/(double)CLOCKS_PER_SEC;

  ///Show results
  output_file_name = (char*) malloc(sizeof(char) * strlen(argv[2]));
  if (output_file_name == NULL)
  {
    printf("Error in function main: Unable to allocate memory");
    exit(1);
  }
  strcpy(output_file_name, argv[2]);
  output_plot_file(binary_sequence, sequence_size, conformation, argv[2]);

  printf("%d %f", conformation.energy, time);


  ///Free memory
  free(input_file);
  free(char_sequence);
  free(binary_sequence);
  free(conformation.directions);
  free(conformation.energy_by_link);
  free(conformation.positions);

  return 0;
}

void output_plot_file(int *seq, int seq_len, Conformation best_solution, char *file_name)
{
  int i;
  int min_x;
  int min_y;
  int max_x;
  int max_y;

  FILE *f = fopen(file_name, "w");

  if (f == NULL)
  {
    printf("Error in function output_plot_file: Unable to open output file\n");
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
