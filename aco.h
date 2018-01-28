struct solution
{
    int energy;
    char *directions;
};

struct ant
{
    int energy;
    int *energy_by_edge;
    struct coord *positions;
};

enum daemon
{
    WITHOUT_DAEMON = 0,
    PULL_MOVE = 1
};

enum collision_handler
{
    PARTIAL_COPY = 1
};

struct aco_config
{
    int population, iterations;
    double alpha, beta, persistence, ini_pheromone;
    enum daemon daemon;
    enum collision_handler collision_handler;
};

typedef struct ant Ant;
typedef enum daemon Daemon;
typedef enum collision_handler Collision_handler;
typedef struct aco_config ACO_config;
typedef struct solution Solution;

void* smalloc(int mem_size);
void free_ant(Ant ant);
void init_solution(Solution *solution, int seq_len);
void free_solution(Solution solution);
void extract_solution(Ant ant, Solution *solution, int seq_len);
Ant aco_run(int *sequence, int sequence_len, ACO_config aco_config, int *seed);
