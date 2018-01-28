struct ant
{
    int energy;
    int *energy_by_edge;
    struct coord *positions;
};

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
typedef struct solution Solution;
typedef enum daemon Daemon;
typedef enum collision_handler Collision_handler;
typedef struct aco_config ACO_config;



void* smalloc(int mem_size);
void init_solution(Solution *solution, int seq_len);
void free_solution(Solution solution);
void init_ant(Ant *ant, int seq_len);
void free_ant(Ant ant);
void extract_solution(Ant ant, Solution *solution, int seq_len);
Ant aco_run(int *sequence, int sequence_len, ACO_config aco_config, int *seed, Ant *ants, int *best_energy_evolution);
