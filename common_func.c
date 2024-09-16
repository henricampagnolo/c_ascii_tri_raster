#include "tri_raster.h"

float ft_power(float x, float n)
{
    if (n < 0)
        return (0);
    if (n == 0)
        return (1);
    return (x * ft_power(x, n - 1));
}

float ft_factorial(float x)
{
    if (x < 0)
        return (0);
    if (x == 0)
        return (1);
    return (x * ft_factorial(x - 1));
}

float ft_sin(float x)
{
    float result;

    while (x >= 2 * PI)
        x -= 2 * PI;
    while (x < 0)
        x += 2 * PI;
    result = 0;
    for (int i = 0; i < TRIG_PRECISION; i++)
    {
        result += (float)(2 * (i % 2 == 0) - 1) * (ft_power(x, 1 + 2 * i) / ft_factorial(1 + 2 * i));
    }
    return (result);
}

float ft_cos(float x)
{
    float result;

    x += PI / 2;
    while (x >= 2 * PI)
        x -= 2 * PI;
    while (x < 0)
        x += 2 * PI;
    result = 0;
    for (int i = 0; i < TRIG_PRECISION; i++)
    {
        result += (float)(2 * (i % 2 == 0) - 1) * (ft_power(x, 1 + 2 * i) / ft_factorial(1 + 2 * i));
    }
    return (result);
}

/*
float get_sin(float x, float *table)
{
    int index;
    float init_x;
    float result;

    init_x = x;
    if (x >= PI && x < PI * 2)
        x -= PI;
    if (x >= (PI / 2) && x < PI)
        x = PI - x;
    index = (x / (PI / 2)) * TABLE_SIZE;
    if (index >= TABLE_SIZE)
        index = TABLE_SIZE - 1;
    result = table[index];
    if (init_x >= PI && init_x < 2 * PI)
        result *= -1;
    return (result);
}
*/

float get_sin(float x, float *table)
{
    int index;
    float init_x;
    float result;

    init_x = x;
    index = (x / (PI / 2)) * TABLE_SIZE;
    if (index >= 4 * TABLE_SIZE)
        index = TABLE_SIZE - 1;
    if (index >= TABLE_SIZE * 2 && index < TABLE_SIZE * 4)
        index -= TABLE_SIZE * 2;
    if (index >= TABLE_SIZE && index < 2 * TABLE_SIZE)
        index = 2 * TABLE_SIZE - index;
    result = table[index];
    if (init_x >= PI && init_x < 2 * PI)
        result *= -1;
    return (result);
}

float get_cos(float x, float *table)
{
    x += PI / 2;
    if (x >= PI * 2)
        x -= PI * 2;

    return (get_sin(x, table));
}

void create_sin_table(float *result)
{
    double index;
    double increment;

    index = 0;
    increment = (PI / 2) / TABLE_SIZE;
    for (int i = 0; i <= TABLE_SIZE; i++)
    {
        result[i] = ft_sin(index);
        index += increment;
    }
}

void ft_putstr(char *str)
{
    while (*str != 0)
        write(1, str++, 1);
}
