#ifndef CUBE_H
#define CUBE_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# ifndef _WIN32
#  include <unistd.h>
# else
#  include <windows.h>
void usleep(__int64 usec)
{
  HANDLE timer;
  LARGE_INTEGER ft;

  ft.QuadPart = -(10 * usec); // Convert to 100 nanosecond interval, negative value indicates relative time

  timer = CreateWaitableTimer(NULL, TRUE, NULL);
  SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0);
  WaitForSingleObject(timer, INFINITE);
  CloseHandle(timer);
}
# endif

# define PI 3.14159265358979
# define TRIG_PRECISION 50
# define TABLE_SIZE 2048

void create_sin_table(float *result);

float get_sin(float x, float *table);

float get_cos(float x, float *table);

void ft_putstr(char *str);

typedef struct linear_func
{
    float a;
    float b;
} t_linear;

typedef struct s_vec2
{
    float x;
    float y;
}   t_vec2;

typedef struct s_vec3
{
    float x;
    float y;
    float z;
}   t_vec3;

typedef struct s_triangle2
{
    t_vec2  points[3];
}   t_triangle2;

typedef struct s_triangle3
{
    t_vec3  points[3];
}   t_triangle3;

#endif