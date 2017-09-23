#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aco.h"

int main()
{

  int s1[20] = {1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0,
                0, 1, 0, 1};

  int s2[24] = {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
                1, 0, 0, 1, 0, 0, 1, 1};

  int s3[25] = {0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                0, 1, 1, 0, 0, 0, 0, 1, 1};

  int s4[36] = {0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1,
                1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0,
                0, 1, 0, 0};

  int s5[48] = {0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1};

  int s6[50] = {1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0,
                0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1,
                0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0,
                1, 1};

  int s7[60] = {0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
                0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0};

  int s8[64] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
                0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
                0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
                1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  int s1_len = 20;
  int s2_len = 24;
  int s3_len = 25;
  int s4_len = 36;
  int s5_len = 48;
  int s6_len = 50;
  int s7_len = 60;
  int s8_len = 64;
  double tempo;
  clock_t t0;
  Aco_config aco_config;
  int count;

  aco_config.alpha = 1;
  aco_config.beta = 2;
  aco_config.ini_pheromone = (double) 1/3;
  aco_config.persistence = 0.8;
  aco_config.iterations = 500;
  aco_config.population = 1500;
  aco_config.best_known_solution = -23;

  Conformation conformation;
  count = 0;

  while (1)
  {
    srand((unsigned) time(NULL));
    ++count;
    t0 = clock();
    conformation = aco_run(s5, s5_len, aco_config);
    tempo = (clock() - t0)/(double)CLOCKS_PER_SEC;
    printf("Tempo(s): %f; Energia: %d; Nr. exe: %d;\n", tempo, conformation.energy, count);
  }
  return 0;
}
