#include "tri_raster.h"

#define WIDTH 80
#define HEIGHT 44

float A, B, C;

float zBuffer[WIDTH * HEIGHT];
char buffer[WIDTH * HEIGHT];
int backgroundASCIICode = '.';
float trig_table[TABLE_SIZE + 1];
float fov = 90;
float fovtan;

// yes indeed

t_triangle2 order_points(t_triangle2 tri)
{
    t_vec2 temp;

    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 3; j++)
        {
            if (tri.points[i].x > tri.points[j].x)
            {
                //printf("i: %d and j: %d\n", i, j);
                //printf("first x: %d, second x: %d\n", tri.points[i].x, tri.points[j].x);
                temp = tri.points[i];
                tri.points[i] = tri.points[j];
                tri.points[j] = temp;
                //printf("first x: %d, second x: %d\n\n", tri.points[i].x, tri.points[j].x);
            }
        }
    }
    return (tri);
}

t_triangle2 project_tri_order(t_triangle3 tri)
{
    t_triangle2 result;
    float   tempz;

    for (int i = 0; i < 3; i++)
    {
        tempz = tri.points[i].z;
        result.points[i].x = tri.points[i].x / tempz / 2 * WIDTH * fovtan + (WIDTH / 2);
        result.points[i].y = tri.points[i].y / tempz / 2 * HEIGHT * fovtan + (HEIGHT / 2);
    }
    result = order_points(result);
    return (result);
}

void get_funcs(t_linear funcs[3], t_triangle2 proj_tri)
{
    t_vec2 p1, p2;

    for (int i = 0; i < 3; i++)
    {
        if (proj_tri.points[i].x >= proj_tri.points[(i + 2) % 3].x)
        {
            p1 = proj_tri.points[(i + 2) % 3];
            p2 = proj_tri.points[i];
        }
        else
        {
            p2 = proj_tri.points[(i + 2) % 3];
            p1 = proj_tri.points[i];
        }
        funcs[i].a = (p2.y - p1.y) / ((p2.x - p1.x) + 0.0000001);
        funcs[i].b = p1.y - p1.x * funcs[i].a;
    }
}

float min2(float a, float b)
{
    if (a < b)
        return a;
    return b;
}

float max2(float a, float b)
{
    if (a > b)
        return a;
    return b;
}

void draw_tri(t_triangle3 tri)
{
    t_triangle2 proj_tri;
    t_linear funcs[3];
    t_vec2 pstart, pmid, pend;

    proj_tri = project_tri_order(tri);

    pstart = proj_tri.points[0], pmid = proj_tri.points[1], pend = proj_tri.points[2]; 

    get_funcs (funcs, proj_tri);
    if (funcs[0].a * pmid.x + funcs[0].b < pmid.y)
    {
        for (int x = pstart.x; x < pend.x; x++)
        {
            for (int y = funcs[0].a * x + funcs[0].b; y < min2(funcs[1].a * x + funcs[1].b, funcs[2].a * x + funcs[2].b); y++)
            {
                //printf("x is: %d, y is: %d\n", x, y);
                if ((x >= 0 && x < WIDTH) && (y >= 0 && y < HEIGHT))
                    buffer[(HEIGHT - y) *  WIDTH + x] = 'Z';
            }
        }
    }
    else
    {
        for (int x = pstart.x; x < pend.x; x++)
        {
            for (int y = funcs[0].a * x + funcs[0].b; y > max2(funcs[1].a * x + funcs[1].b, funcs[2].a * x + funcs[2].b); y--)
            {
                //printf("x is: %d, y is: %d\n", x, y);
                if ((x >= 0 && x < WIDTH) && (y >= 0 && y < HEIGHT))
                    buffer[(HEIGHT - y) * WIDTH + x] = 'Z';
            }
        }
    }
}

void rotate_point_x(t_vec3 *p_p, float rad)
{
    float y = p_p->y;
    float z = p_p->z;
    float sin = get_sin(rad, trig_table);
    float cos = get_cos(rad, trig_table);
    p_p->y = cos * y + sin * z;
    p_p->z = -sin * y + cos * z; 
}

void rotate_point_y(t_vec3 *p_p, float rad)
{
    float x = p_p->x;
    float z = p_p->z;
    float sin = get_sin(rad, trig_table);
    float cos = get_cos(rad, trig_table);
    p_p->x = cos * x + -sin * z;
    p_p->z = sin * x + cos * z; 
}

void rotate_point_z(t_vec3 *p_p, float rad)
{
    float x = p_p->x;
    float y = p_p->y;
    float sin = get_sin(rad, trig_table);
    float cos = get_cos(rad, trig_table);
    p_p->x = cos * x + sin * y;
    p_p->y = -sin * x + cos * y;
}

int main(void)
{
    float z;
    create_sin_table(trig_table);
    fovtan = (get_cos(fov / 360 * PI, trig_table) / get_sin(fov / 360 * PI, trig_table));

    t_triangle3 tri;
    t_triangle3 tri_cpy;
    tri.points[0].x = -10, tri.points[0].y = 0, tri.points[0].z = 80;
    tri.points[1].x = 10, tri.points[1].y = 0, tri.points[1].z = 80;
    tri.points[2].x = 0, tri.points[2].y = 70, tri.points[2].z = 80;
    printf("\x1b[2J");

    z = 0;
    while (1)
    {
        tri_cpy = tri;
        rotate_point_z(tri_cpy.points, z);
        rotate_point_z(tri_cpy.points + 1, z);
        rotate_point_z(tri_cpy.points + 2, z);

        memset(buffer, backgroundASCIICode, WIDTH * HEIGHT);
        memset(zBuffer, 0, WIDTH * HEIGHT * 4);

        draw_tri(tri_cpy);

        printf("\x1b[H");
        for (int k = 0; k < WIDTH * HEIGHT; k++) 
        {
            putchar(k % WIDTH ? buffer[k] : 10);
        }
        if (z >= PI * 2)
            z = 0;
        z += 0.002;
        printf("z is: %f,  sin of z: %f, cos of z: %f", z, get_sin(z, trig_table), get_cos(z, trig_table));
    }
}