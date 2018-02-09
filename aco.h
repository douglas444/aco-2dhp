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

struct aco_result
{
    char *directions;
    int energy;
    double time;
    int found_on_iteration;
    float final_population_avg;
    float final_population_stddev;
    float final_population_solution_rate;
};

typedef enum daemon Daemon;
typedef enum collision_handler Collision_handler;
typedef struct aco_config ACO_config;
typedef struct solution Solution;
typedef struct aco_result Aco_result;

void* smalloc(int mem_size);
void free_solution(Solution solution);
Aco_result aco_run(int *sequence, int sequence_len, ACO_config aco_config, int *seed);
