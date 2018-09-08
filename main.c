#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "aco.h"
#include "file.h"

void set_max_priority();
enum daemon string_to_daemon(char *daemon);
enum daemon string_to_constructor(char *constructor);
enum daemon char_to_polarity(char polarity);
void read_inputs(struct aco_config *aco_config, enum polarity **sequence, int *sequence_len, char **argv, char **char_sequence);

int main(int argc, char **argv)
{
    //set_max_priority();

    int i;
    char *char_sequence;
    struct aco_config aco_config;
    struct aco_result aco_result;
    enum polarity *sequence;
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

enum daemon string_to_daemon(char *daemon)
{
    if (strcmp(daemon, "NONE") == 0)
    {
        return NONE;
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

enum daemon string_to_constructor(char *constructor)
{
    if (strcmp(constructor, "XIAO_LI_HU_2014") == 0)
    {
        return XIAO_LI_HU_2014;
    }
    else if (strcmp(constructor, "HU_ZHANG_LI_2009") == 0)
    {
        return HU_ZHANG_LI_2009;
    }
    else if (strcmp(constructor, "SHMYGELSKA_HOOS_2003") == 0)
    {
        return SHMYGELSKA_HOOS_2003;
    }
    else
    {
        printf("ERROR: main.c/string_to_constructor(): \"Invalid value for constructor parameter\"\n");
        exit(1);
    }
}

enum daemon char_to_polarity(char polarity)
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

void read_inputs(struct aco_config *aco_config, enum polarity **sequence, int *sequence_len, char **argv, char **char_sequence)
{
    int i;
    char *input_file;
    char *constructor;
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
    aco_config->min_probability = char_to_double(get_key_value(input_file, "min-probability"));
    aco_config->iterations = char_to_int(get_key_value(input_file, "iterations"));
    aco_config->elit_percentage = char_to_double(get_key_value(input_file, "elit-percentage"));

    *char_sequence = get_key_value(input_file, sequence_key);
    constructor = get_key_value(input_file, "constructor");
    daemon = get_key_value(input_file, "daemon");

    *sequence_len = strlen(*char_sequence);
    aco_config->population = char_to_int(get_key_value(input_file,
        (*sequence_len <= 25) ? "small-instances-population" : "big-instances-population"));
    aco_config->daemon = string_to_daemon(daemon);
    aco_config->constructor = string_to_constructor(constructor);

    *sequence = (enum polarity*) malloc(*sequence_len * sizeof(enum polarity));
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
    free(constructor);
    free(daemon);
    free(sequence_key);
}
