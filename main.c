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
    int sequence_len;
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
    sequence_key = (char*) malloc(sizeof(char) * (strlen(argv[2]) + 1));
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
    sequence_len= strlen(char_sequence);
    if (sequence_len <= 25)
    {
        aco_config.population = char_to_int(get_key_value(input_file, "small-instances-population"));
    } else {
        aco_config.population = char_to_int(get_key_value(input_file, "big-instances-population"));
    }
    collision_handler = get_key_value(input_file, "collision-handler");
    daemon = get_key_value(input_file, "daemon");

    ///Set local search method
    if (strcmp(daemon, "WITHOUT_LOCAL_SEARCH") == 0)
    {
        aco_config.daemon = WITHOUT_LOCAL_SEARCH;
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
    binary_sequence = (int*) malloc(sequence_len * sizeof(int));
    if (binary_sequence == NULL)
    {
        printf("Error in function main: Unable to allocate memory");
        exit(1);
    }
    for (i = 0; i < sequence_len; ++i)
    {
        binary_sequence[i] = char_sequence[i] == 'H' ? 1 : 0;
    }

    ///Run ACO
    int seed = -1;
    t0 = clock();
    Solution solution = aco_run(binary_sequence, sequence_len, aco_config, &seed);
    time = (clock() - t0)/(double)CLOCKS_PER_SEC;

    ///Show results
    printf("%d %f %s", solution.energy, time, solution.directions);


    ///Free memory
    free(sequence_key);
    free(input_file);
    free(char_sequence);
    free(binary_sequence);
    free(solution.directions);

    return 0;
}
