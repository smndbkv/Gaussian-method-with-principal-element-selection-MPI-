#ifndef IO_STATUS_H
#define IO_STATUS_H
#define EPS_64 (1e-64)
#define EPS (1.12e-16)
#define LEN (1234)
#define CMP(x, y) \
    ((fabs(((x) - (y))) <= fabs(((0.5) * (EPS) * ((x) + (y))))) ? 1 : 0)
#define MAX_PRINT 10
#define SWAP(x, y)        \
    {                     \
        double _ = x;     \
        x = (double &&)y; \
        y = _;            \
    }
typedef enum io_status_
{
    SUCCESS,
    ERROR_OPEN,
    ERROR_READ
} io_status;
typedef enum gauss_status_
{
    DONE,
    ZERO_MATRIX,
    NOT_APPLICABLE
} gauss_status;

#endif
