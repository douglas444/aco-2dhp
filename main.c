#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aco.h"

int main()
{

    //srand((unsigned)time(NULL));

    int s1[20] = {1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0,
                  0, 1, 0, 1},

        s2[24] = {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
                  1, 0, 0, 1, 0, 0, 1, 1},

        s3[25] = {0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                  0, 1, 1, 0, 0, 0, 0, 1, 1},

        s4[36] = {0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                  1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0,
                  0, 1, 0, 0},

        s5[48] = {0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1},

        s6[50] = {1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0,
                  0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1,
                  0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                  1, 1},

        s7[60] = {0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0,
                  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                  1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0},

        s8[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
                  0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
                  0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
                  1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    int s1_len = 20,
        s2_len = 24,
        s3_len = 25,
        s4_len = 36,
        s5_len = 48,
        s6_len = 50,
        s7_len = 60,
        s8_len = 64;

    int i, energy,
        population_size = 200,
        num_iterations = 500;

    float alpha = 1,
          beta = 2,
          evaporation_rate = 0.8,
          initial_phero = (float) 1 / 3;

    printf("s1:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s1, s1_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -9);
        printf("%d\n", energy);
    }

    printf("s2:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s2, s2_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -9);
        printf("%d\n", energy);
    }

    printf("s3:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s3, s3_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -8);
        printf("%d\n", energy);
    }

    population_size = 1500;

    printf("s4:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s4, s4_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -14);
        printf("%d\n", energy);
    }

    printf("s5:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s5, s5_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -23);
        printf("%d\n", energy);
    }

    printf("s6:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s6, s6_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -21);
        printf("%d\n", energy);
    }

    printf("s7:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s7, s7_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -36);
        printf("%d\n", energy);
    }

    printf("s8:\n");
    for (i = 0; i < 20; ++i)
    {
        energy = aco_run(s8, s8_len, alpha, beta, evaporation_rate, initial_phero, population_size, num_iterations, -42);
        printf("%d\n", energy);
    }

    return 0;
}
