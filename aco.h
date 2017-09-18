typedef struct coord
{
    int x;
    int y;

} Coord;

typedef struct aco_config
{
  int population;
  int iterations;
  int best_known_solution;
  float alpha;
  float beta;
  float persistence;
  float ini_pheromone;

} Aco_config;

typedef enum direction
{
    LEFT = 0,
    RIGHT = 1,
    STRAIGHT = 2
} Direction;

typedef struct conformation
{
  int *energy_by_link; /*Reached energy before i th link*/
  int length;
  int energy;
  Coord *positions;
  Direction *directions;

} Conformation;


Conformation aco_run(int *seq, int seq_len, Aco_config aco_config, int actived_pull_move);
