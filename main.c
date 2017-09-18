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

  Aco_config aco_config;

  aco_config.alpha = 1;
  aco_config.beta = 1;
  aco_config.ini_pheromone = (float) 1 / 3;
  aco_config.persistence = 0.8;
  aco_config.iterations = 500;
  aco_config.population = 1500;
  aco_config.best_known_solution = -36;

  int i;
  for (i = 0; i < 20; ++i)
  {
    Conformation conformation = aco_run(s8, s8_len, aco_config, 1);
    printf("%d\n", conformation.energy);
  }
  return 0;
}
