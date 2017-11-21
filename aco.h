typedef struct coord
{
  int x;
  int y;
} Coord;

typedef enum local_search
{
  WITHOUT_LOCAL_SEARCH = 0,
  PULL_MOVE = 1
} Local_search;

typedef enum unfeasible_handler
{
  PARTIAL_COPY = 1,
  BLOCKED_POSITIONS = 2
} Unfeasible_conformation_handler;

typedef struct aco_config
{
  int population;
  int iterations;

  double alpha;
  double beta;
  double persistence;
  double ini_pheromone;

  Local_search local_search;
  Unfeasible_conformation_handler unfeasible_conformation_handler;

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
