#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "aco.h"
#include "file.h"

void set_max_priority();
Daemon string_to_daemon(char *daemon);
Daemon string_to_collision_handler(char *collision_handler);
Daemon char_to_polarity(char polarity);
void read_inputs(Aco_config *aco_config, Polarity **sequence, int *sequence_len, char **argv, char **char_sequence);

int main(int argc, char **argv)
{
    set_max_priority();

    int i;
    char *char_sequence;
    Aco_config aco_config;
    Aco_result aco_result;
    Polarity *sequence;
    int sequence_len;
    int seed = -1;

    if (argc < 3)
    {
        printf("ERROR: main.c/main(): \"Invalid application parameters\"\n");
        exit(1);
    }

    seed = -1;
    read_inputs(&aco_config, &sequence, &sequence_len, argv, &char_sequence);
    aco_result = aco_run(sequence, sequence_len, aco_config, &seed);

    //output
    printf("%d|", aco_config.iterations);
    printf("%s|", char_sequence);
    printf("%s|", aco_result.directions);
    printf("%d|", aco_result.energy);
    printf("%.2f|", aco_result.final_population_avg);
    printf("%.2f|", aco_result.final_population_stddev);
    printf("%.2f|", aco_result.final_population_solution_rate);
    printf("%d|", aco_result.found_on_iteration);
    printf("%.2f|", aco_result.time);
    for (i = 0; i < aco_config.iterations; ++i) {
        if (aco_result.energy_evolution[i] != 1) {
            printf("%d,%d/", i, aco_result.energy_evolution[i]);
        }
    }

    //free
    free(aco_result.directions);
    free(aco_result.energy_evolution);
    free(sequence);
    free(char_sequence);

    return 0;
}

void set_max_priority() {

    if (setpriority(PRIO_PROCESS, 0, -20) == -1) {
        printf("ERROR: main.c/set_max_priority(): \"Setpriority() failed\"\n");
        exit(1);
    }

}

Daemon string_to_daemon(char *daemon)
{
    if (strcmp(daemon, "WITHOUT_DAEMON") == 0)
    {
        return WITHOUT_DAEMON;
    }
    else if (strcmp(daemon, "PULL_MOVE") == 0)
    {
        return PULL_MOVE;
    }
    else
    {
        printf("ERROR: main.c/string_to_daemon(): \"Invalid value for daemon parameter\"\n");
        exit(1);
    }
}

Daemon string_to_collision_handler(char *collision_handler)
{
    if (strcmp(collision_handler, "PARTIAL_COPY") == 0)
    {
        return PARTIAL_COPY;
    }
    else
    {
        printf("ERROR: main.c/string_to_collision_handler(): \"Invalid value for collision_handler parameter\"\n");
        exit(1);
    }
}

Daemon char_to_polarity(char polarity)
{
    if (polarity == 'H')
    {
        return H;
    }
    else if (polarity == 'P')
    {
        return P;
    }
    else
    {
        printf("ERROR: main.c/char_to_polarity(): \"Invalid value for polarity parameter\"\n");
        exit(1);
    }

}

void read_inputs(Aco_config *aco_config, Polarity **sequence, int *sequence_len, char **argv, char **char_sequence)
{
    int i;
    char *input_file;
    char *collision_handler;
    char *daemon;
    char *sequence_key;

    sequence_key = (char*) malloc(sizeof(char) * (strlen(argv[2]) + 1));
    if (sequence_key == NULL)
    {
        printf("ERROR: main.c/read_inputs(): \"Unable to allocate memory\"\n");
        exit(1);
    }

    //get inputs
    strcpy(sequence_key, argv[2]);
    input_file = load_file_content(argv[1]);

    //extract from inputs
    aco_config->alpha = char_to_double(get_key_value(input_file, "alpha"));
    aco_config->beta = char_to_double(get_key_value(input_file, "beta"));
    aco_config->ini_pheromone = char_to_double(get_key_value(input_file, "initial-pheromone"));
    aco_config->persistence = char_to_double(get_key_value(input_file, "persistence"));
    aco_config->iterations = char_to_int(get_key_value(input_file, "iterations"));

    *char_sequence = get_key_value(input_file, sequence_key);
    collision_handler = get_key_value(input_file, "collision-handler");
    daemon = get_key_value(input_file, "daemon");

    *sequence_len = strlen(*char_sequence);
    aco_config->population = char_to_int(get_key_value(input_file,
        (*sequence_len <= 25) ? "small-instances-population" : "big-instances-population"));
    aco_config->daemon = string_to_daemon(daemon);
    aco_config->collision_handler = string_to_collision_handler(collision_handler);

    *sequence = (Polarity*) malloc(*sequence_len * sizeof(Polarity));
    if (*sequence == NULL)
    {
        printf("ERROR: main.c/read_inputs(): \"Unable to allocate memory\"\n");
        exit(1);
    }
    for (i = 0; i < *sequence_len; ++i)
    {
        (*sequence)[i] = char_to_polarity((*char_sequence)[i]);
    }

    //free
    free(input_file);
    free(collision_handler);
    free(daemon);
    free(sequence_key);
}
