struct solution
{
    int energy;
    char *directions;
};

enum daemon
{
    NONE = 0,
    PULL_MOVE = 1
};

enum constructor
{
    XIAO_LI_HU_2014 = 0,
    HU_ZHANG_LI_2009 = 1,
    SHMYGELSKA_HOOS_2003 = 2
};

enum polarity
{
    P = 0,
    H = 1
};

struct aco_config
{
    int population, iterations;
    double alpha, beta, persistence, ini_pheromone, elit_percentage;
    enum daemon daemon;
    enum constructor constructor;
    double min_probability;
};

struct aco_result
{
    char *directions;
    int *energy_evolution;
    int energy;
    double time;
    int found_on_iteration;
    float final_population_avg;
    float final_population_stddev;
    float final_population_solution_rate;
};

void* smalloc(int mem_size);
void free_solution(struct solution solution);
struct aco_result aco_run(enum polarity *sequence, int sequence_len, struct aco_config aco_config, int *seed);
