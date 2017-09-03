typedef struct coord
{
    int x, y;
} Coord;

typedef enum direction
{
    LEFT = 0,
    RIGHT = 1,
    STRAIGHT = 2
} Direction;

int aco_run
(
    int *seq,
    int seq_len,
    float alpha,
    float beta,
    float evaporation_rate,
    float initial_phero,
    int population_size,
    int num_iterations,
    int stop_criterion
);
