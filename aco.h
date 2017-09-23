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
  double alpha;
  double beta;
  double persistence;
  double ini_pheromone;

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


Conformation aco_run(int *seq, int seq_len, Aco_config aco_config);
