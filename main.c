#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "aco.h"
#include "file.h"

int main(int argc, char **argv)
{

    int i;
    int *binary_sequence;
    ACO_config aco_config;
    char *input_file;
    char *collision_handler;
    char *daemon;
    int seq_len;
    char *char_sequence;
    clock_t t0;
    double time;
    char *sequence_key;

    if (argc < 3)
    {
        printf("Invalid parameters\n");
        exit(1);
    }

    ///Extract file content
    input_file = load_file_content(argv[1]);

    ///Sequence key
    sequence_key = (char*) smalloc(sizeof(char) * (strlen(argv[2]) + 1));
    if (sequence_key == NULL)
    {
        printf("Error in function main: Unable to allocate memory");
        exit(1);
    }
    strcpy(sequence_key, argv[2]);

    ///Get parameters from file content
    aco_config.alpha = char_to_double(get_key_value(input_file, "alpha"));
    aco_config.beta = char_to_double(get_key_value(input_file, "beta"));
    aco_config.ini_pheromone = char_to_double(get_key_value(input_file, "initial-pheromone"));
    aco_config.persistence = char_to_double(get_key_value(input_file, "persistence"));
    aco_config.iterations = char_to_int(get_key_value(input_file, "iterations"));
    char_sequence = get_key_value(input_file, sequence_key);
    seq_len= strlen(char_sequence);

    if (seq_len <= 25)
    {
        aco_config.population = char_to_int(get_key_value(input_file, "small-instances-population"));
    } else {
        aco_config.population = char_to_int(get_key_value(input_file, "big-instances-population"));
    }
    collision_handler = get_key_value(input_file, "collision-handler");
    daemon = get_key_value(input_file, "daemon");

    ///Set local search method
    if (strcmp(daemon, "WITHOUT_DAEMON") == 0)
    {
        aco_config.daemon = WITHOUT_DAEMON;
    }
    else if (strcmp(daemon, "PULL_MOVE") == 0)
    {
        aco_config.daemon = PULL_MOVE;
    }
    else
    {
        printf("Error in function main: Local search parameter not recognized\n");
        exit(1);
    }

    ///Set collision handler method
    if (strcmp(collision_handler, "PARTIAL_COPY") == 0)
    {
        aco_config.collision_handler = PARTIAL_COPY;
    }
    else
    {
        printf("Error in function main: Collision handler parameter not recognized\n");
        exit(1);
    }

    ///generate protein binary sequence
    binary_sequence = (int*) smalloc(seq_len * sizeof(int));
    if (binary_sequence == NULL)
    {
        printf("Error in function main: Unable to allocate memory");
        exit(1);
    }
    for (i = 0; i < seq_len; ++i)
    {
        binary_sequence[i] = char_sequence[i] == 'H' ? 1 : 0;
    }

    ///Run ACO
    int seed = -1;
    Ant *ants;
    Solution best_solution;
    Solution *final_solutions;

    init_solution(&best_solution, seq_len);
    ants = (Ant*) smalloc(sizeof(Ant) * aco_config.population);
    final_solutions = (Solution*) smalloc(sizeof(Solution) * aco_config.population);
    for (i = 0; i < aco_config.population; ++i)
    {
        init_ant(&(ants[i]), seq_len);
        init_solution(&(final_solutions[i]), seq_len);

    }

    t0 = clock();
    Ant best_ant = aco_run(binary_sequence, seq_len, aco_config, &seed, ants);
    time = (clock() - t0)/(double)CLOCKS_PER_SEC;

    extract_solution(best_ant, &best_solution, seq_len);

    ///Show best_ant results
    printf("%d %f %s ", best_solution.energy, time, best_solution.directions);

    for (i = 0; i < seq_len; ++i)
    {
        printf("%c", char_sequence[i]);
    }
    printf(" ");

    for (i = 0; i < aco_config.population; ++i)
    {
        extract_solution(ants[i], &(final_solutions[i]), seq_len);
        printf("[%d][%s]\\n", final_solutions[i].energy, final_solutions[i].directions);
    }

    ///Free memory
    for (i = 0; i < aco_config.population; ++i)
    {
        free_ant(ants[i]);
        free_solution(final_solutions[i]);
    }
    free(sequence_key);
    free(input_file);
    free(char_sequence);
    free(binary_sequence);

    return 0;
}
