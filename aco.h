typedef struct solution
{
    int energy;
    char *directions;
} Solution;

typedef enum daemon
{
    WITHOUT_LOCAL_SEARCH = 0,
    PULL_MOVE = 1
} Daemon;

typedef enum collision_handler
{
    PARTIAL_COPY = 1
} Collision_handler;

typedef struct aco_config
{
    int population, iterations;
    double alpha, beta, persistence, ini_pheromone;
    Daemon daemon;
    Collision_handler collision_handler;

} ACO_config;

Solution aco_run(int *sequence, int sequence_len, ACO_config aco_config, int *seed);
