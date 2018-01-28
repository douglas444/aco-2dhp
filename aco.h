struct solution
{
    int energy;
    char *directions;
};

enum daemon
{
    WITHOUT_DAEMON = 0,
    PULL_MOVE = 1
};

enum collision_handler
{
    PARTIAL_COPY = 0
};

struct aco_config
{
    int population, iterations;
    double alpha, beta, persistence, ini_pheromone;
    enum daemon daemon;
    enum collision_handler collision_handler;
};

typedef enum daemon Daemon;
typedef enum collision_handler Collision_handler;
typedef struct aco_config ACO_config;
typedef struct solution Solution;

void* smalloc(int mem_size);
void free_solution(Solution solution);
Solution aco_run(int *sequence, int sequence_len, ACO_config aco_config, int *seed);
