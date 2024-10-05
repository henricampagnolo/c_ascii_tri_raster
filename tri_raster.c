#include "tri_raster.h"

#define WIDTH 80
#define HEIGHT 40

float A, B, C;

float zBuffer[WIDTH * HEIGHT];
char buffer[WIDTH * HEIGHT];
int backgroundASCIICode = '.';
float trig_table[TABLE_SIZE + 1];
float fov = 40;
float fovtan;

// yes indeed

void order_points(t_triangle3 *tri)
{
    t_vec3 *temp;

    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 3; j++)
        {
            if (tri->points[i]->x > tri->points[j]->x)
            {
                //printf("i: %d and j: %d\n", i, j);
                //printf("first x: %d, second x: %d\n", tri.points[i].x, tri.points[j].x);
                temp = tri->points[i];
                tri->points[i] = tri->points[j];
                tri->points[j] = temp;
                //printf("first x: %d, second x: %d\n\n", tri.points[i].x, tri.points[j].x);
            }
        }
    }
}

void project_tri_order(t_triangle3 tri, t_triangle3 *result)
{
    float   tempz;

    for (int i = 0; i < 3; i++)
    {
        tempz = tri.points[i]->z;
        result->points[i]->x = tri.points[i]->x / tempz / 2 * WIDTH * fovtan + (WIDTH / 2);
        result->points[i]->y = tri.points[i]->y / tempz / 2 * HEIGHT * fovtan + (HEIGHT / 2);
        result->points[i]->z = tri.points[i]->z;
    }
    order_points(result);
}

void get_slopes(float *h_slope, float *v_slope, t_vec3 p0, t_vec3 p1, t_vec3 p2)
{
    float dx1, dx2;
    float dy1, dy2;
    float dz1, dz2;

    dx1 = p1.x - p0.x;
    dx2 = p2.x - p0.x;
    dy1 = p1.y - p0.y;
    dy2 = p2.y - p0.y;
    dz1 = p1.z - p0.z;
    dz2 = p2.z - p0.z;

    *v_slope = (dz2 * dx1 - dz1 * dx2) / (dx1 * dy2 - dx2 * dy1);
    *h_slope = (dz1 - dy1 * (*v_slope)) / dx1;  
}

void get_funcs(t_linear funcs[3], t_triangle3 proj_tri)
{
    t_vec3 p1, p2;

    for (int i = 0; i < 3; i++)
    {
        if (proj_tri.points[i][0].x >= proj_tri.points[(i + 2) % 3][0].x)
        {
            p1 = proj_tri.points[(i + 2) % 3][0];
            p2 = proj_tri.points[i][0];
        }
        else
        {
            p2 = proj_tri.points[(i + 2) % 3][0];
            p1 = proj_tri.points[i][0];
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

float map_to_range(float x, float min, float max)
{
	return min2(max2(x, min), max);
}

int round_up(float n)
{
	return (int)n + 1;
}

void draw_tri(t_triangle3 tri, char symbol)
{
    t_triangle3 proj_tri;
    t_vec3 point1, point2, point3;
    proj_tri.points[0] = &point1;
    proj_tri.points[1] = &point2;
    proj_tri.points[2] = &point3;
    t_linear funcs[3];
    t_vec3 pstart, pmid, pend;
    float start_z, h_slope, current_z, v_slope;

    project_tri_order(tri, &proj_tri);

    pstart = proj_tri.points[0][0], pmid = proj_tri.points[1][0], pend = proj_tri.points[2][0]; 

    get_slopes (&h_slope, &v_slope, pstart, pmid, pend);
    start_z = pstart.z;

    get_funcs (funcs, proj_tri);
    if (funcs[0].a * pmid.x + funcs[0].b < pmid.y)
    {
        for (int x = max2(pstart.x, 0); x < min2(pend.x, WIDTH); x++)
        {
            for (int y = round_up(map_to_range(funcs[0].a * x + funcs[0].b, 0, HEIGHT - 1));
				y < min2(min2(funcs[1].a * x + funcs[1].b, funcs[2].a * x + funcs[2].b), HEIGHT);
				y++)
            {
                current_z = start_z + (x - pstart.x) * h_slope + (y - pstart.y) * v_slope;
                //printf("x is: %d, y is: %d and current: %f and zbuf: %f\n", x, y, current_z, zBuffer[(HEIGHT - y) *  WIDTH + x]);
                if (1 / current_z > (zBuffer[(HEIGHT - y) *  WIDTH + x] + 0.000001))
                {
                    //printf("printing :: x is: %d, y is: %d", x, y);
                    zBuffer[(HEIGHT - y) *  WIDTH + x] = 1 / current_z;
                    buffer[(HEIGHT - y) *  WIDTH + x] = symbol;
                }
            }
        }
    }
    else
    {
        for (int x = max2(pstart.x, 0); x < min2(pend.x, WIDTH); x++)
        {
            for (int y = map_to_range(funcs[0].a * x + funcs[0].b, 0, HEIGHT - 1);
				y >= max2(max2(funcs[1].a * x + funcs[1].b, funcs[2].a * x + funcs[2].b), 0);
				y--)
			{
                current_z = start_z + (x - pstart.x) * h_slope + (y - pstart.y) * v_slope;
                //printf("x is: %d, y is: %d and current: %f and zbuf: %f\n", x, y, current_z, zBuffer[(HEIGHT - y) *  WIDTH + x]);
                if (1 / current_z > (zBuffer[(HEIGHT - y) *  WIDTH + x] + 0.000001))
                {
                    //printf("printing :: x is: %d, y is: %d", x, y);
                    zBuffer[(HEIGHT - y) *  WIDTH + x] = 1 / current_z;
                    buffer[(HEIGHT - y) *  WIDTH + x] = symbol;
                }
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

void define_triangle(t_triangle3 *tri, t_vec3 *p1, t_vec3 *p2, t_vec3 *p3)
{
    tri->points[0] = p1;
    tri->points[1] = p2;
    tri->points[2] = p3;
}

void draw_shape(t_shape shape)
{
    for (int i = 0; i < shape.t_size; i++)
	{
		if (i != 6)
        	draw_tri(shape.tris[i], shape.texture[i]);
	}
}

void add_vectors(t_vec3 *og, t_vec3 delta)
{
    og->x += delta.x;
    og->y += delta.y;
    og->z += delta.z;
}

t_vec3 inv_vec(t_vec3 v)
{
    t_vec3 result;

    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return (result);
}

void move_shape(t_shape *shape, t_vec3 mvmt)
{
    for (int i = 0; i < shape->p_size; i++)
    {
        add_vectors(shape->points + i, mvmt);
    }
    add_vectors(&(shape->center), mvmt);
}

void normalize_angle(float *angle)
{
    while (*angle >= PI * 2)
        *angle -= PI * 2;
    while (*angle < 0)
        *angle += PI * 2;
}

void turn_shape(t_shape shape, float x, float y, float z)
{
    normalize_angle(&x);
    normalize_angle(&y);
    normalize_angle(&z);

    t_vec3 og_center = shape.center;
    //printf(" center x: %f, y: %f, z: %f\n", shape.center.x, shape.center.y, shape.center.z);
    //printf(" start: p1 x: %f, y: %f, z: %f\n", shape.points[0].x, shape.points[0].y, shape.points[0].z);
    move_shape(&shape, inv_vec(og_center));
    //printf(" mov1 : p1 x: %f, y: %f, z: %f\n", shape.points[0].x, shape.points[0].y, shape.points[0].z);
    
    for (int i = 0; i < shape.p_size; i++)
    {
        rotate_point_x(shape.points + i, x);
        rotate_point_y(shape.points + i, y);
        rotate_point_z(shape.points + i, z);
    }
    //printf(" rot : p1 x: %f, y: %f, z: %f\n", shape.points[0].x, shape.points[0].y, shape.points[0].z);
    move_shape(&shape, og_center);
    //printf(" mov 2: p1 x: %f, y: %f, z: %f\n\n", shape.points[0].x, shape.points[0].y, shape.points[0].z);
}

void assign_vec(t_vec3 *vec, float x, float y, float z)
{
	vec->x = x;
	vec->y = y;
	vec->z = z;
}

int main(void)
{
	float alpha;
    float beta;
    create_sin_table(trig_table);
    fovtan = (get_cos(fov / 360 * PI, trig_table) / get_sin(fov / 360 * PI, trig_table));

	/*
	*/
    t_vec3 tetra_center;
    tetra_center.x = 0, tetra_center.y = 0, tetra_center.z = 0;

    t_vec3 tetra_points[4];
    tetra_points[0].x = 0, tetra_points[0].y = 40, tetra_points[0].z = 0;
    tetra_points[1].x = 32.66, tetra_points[1].y = -13.34, tetra_points[1].z = -18.86;
    tetra_points[2].x = -32.66, tetra_points[2].y = -13.34, tetra_points[2].z = -18.86;
    tetra_points[3].x = 0, tetra_points[3].y = -13.34, tetra_points[3].z = 37.71;

    t_shape tetrahedron;

    tetrahedron.t_size = 4;
    tetrahedron.p_size = 4;
    tetrahedron.points = tetra_points;
    t_triangle3 tetra_tris[4];
    define_triangle(tetra_tris, tetra_points + 0, tetra_points + 1, tetra_points + 2);
    define_triangle(tetra_tris + 1, tetra_points + 0, tetra_points + 1, tetra_points + 3);
    define_triangle(tetra_tris + 2, tetra_points + 0, tetra_points + 2, tetra_points + 3);
    define_triangle(tetra_tris + 3, tetra_points + 1, tetra_points + 2, tetra_points + 3);
    tetrahedron.tris = tetra_tris;
    tetrahedron.center = tetra_center;
    tetrahedron.texture = "#a/-";
	/*
	t_vec3 hex_center;
	assign_vec(&hex_center, 0, 0, 0);

	t_vec3 hex_points[6];
	assign_vec(hex_points, 0, -40, 0);
	assign_vec(hex_points + 1, 34.64, -20, 0);
	assign_vec(hex_points + 2, 34.64, 20, 0);
	assign_vec(hex_points + 3, -34.64, -20, 0);
	assign_vec(hex_points + 4, -34.64, 20, 0);
	assign_vec(hex_points + 5, 0, 40, 0);

	t_shape hexagon;
	
	hexagon.center = hex_center;
	hexagon.p_size = 6;
	hexagon.t_size = 4;
	t_triangle3 hex_tris[4];
	define_triangle(hex_tris, hex_points + 0, hex_points + 1, hex_points + 3);
	define_triangle(hex_tris + 1, hex_points + 1, hex_points + 2, hex_points + 3);
	define_triangle(hex_tris + 2, hex_points + 4, hex_points + 2, hex_points + 3);
	define_triangle(hex_tris + 3, hex_points + 4, hex_points + 2, hex_points + 5);
	hexagon.tris = hex_tris;
	hexagon.points = hex_points;
	hexagon.texture = "ABCD";
	*/
	/*
	t_vec3 tetra_center;
    tetra_center.x = 0, tetra_center.y = 0, tetra_center.z = 0;

    t_vec3 cube_points[8];
    cube_points[0].x = 40, cube_points[0].y = 40, cube_points[0].z = 40;
    cube_points[1].x = 40, cube_points[1].y = 40, cube_points[1].z = -40;
    cube_points[2].x = 40, cube_points[2].y = -40, cube_points[2].z = 40;
    cube_points[3].x = 40, cube_points[3].y = -40, cube_points[3].z = -40;
	cube_points[4].x = -40, cube_points[4].y = 40, cube_points[4].z = 40;
	cube_points[5].x = -40, cube_points[5].y = 40, cube_points[5].z = -40;
	cube_points[6].x = -40, cube_points[6].y = -40, cube_points[6].z = 40;
	cube_points[7].x = -40, cube_points[7].y = -40, cube_points[7].z = -40;

    t_shape cube;

    cube.size = 12;
    cube.points = cube_points;
    t_triangle3 cube_tris[12];
    define_triangle(cube_tris, cube_points + 0, cube_points + 1, cube_points + 2);
    define_triangle(cube_tris + 1, cube_points + 3, cube_points + 1, cube_points + 2);
    define_triangle(cube_tris + 2, cube_points + 0, cube_points + 2, cube_points + 4);
    define_triangle(cube_tris + 3, cube_points + 6, cube_points + 2, cube_points + 4);
    define_triangle(cube_tris + 4, cube_points + 1, cube_points + 3, cube_points + 5);
    define_triangle(cube_tris + 5, cube_points + 7, cube_points + 3, cube_points + 5);
    define_triangle(cube_tris + 6, cube_points + 4, cube_points + 5, cube_points + 6);
    define_triangle(cube_tris + 7, cube_points + 7, cube_points + 5, cube_points + 6);
    define_triangle(cube_tris + 8, cube_points + 0, cube_points + 1, cube_points + 2);
    define_triangle(cube_tris + 9, cube_points + 0, cube_points + 1, cube_points + 2);
    define_triangle(cube_tris + 10, cube_points + 0, cube_points + 1, cube_points + 2);
    define_triangle(cube_tris + 11, cube_points + 0, cube_points + 1, cube_points + 2);

    tetrahedron.tris = tetra_tris;
    tetrahedron.center = tetra_center;
    tetrahedron.texture = "#L/-";
	*/

    printf("\x1b[2J");

	alpha = 0;
    beta = 0;
    int i = 0;

    t_vec3 new_pos;
    new_pos.x = 0, new_pos.y = 0, new_pos.z = 130;

    move_shape(&tetrahedron, new_pos);

    while (i < 30000)
    {
        i += 1;

        memset(buffer, backgroundASCIICode, WIDTH * HEIGHT);
        memset(zBuffer, 1024.0, WIDTH * HEIGHT * 4);

        draw_shape(tetrahedron);

        printf("\x1b[H");
        for (int k = 0; k < WIDTH * HEIGHT; k++) 
        {
            putchar(k % WIDTH ? buffer[k] : 10);
			//putchar(k % WIDTH ? ' ' : '\t');
        }
		/*
		*/
        
        turn_shape(tetrahedron, -0.002, 0.002, 0);
        alpha+=0.0005;
        beta+=0.0005;

        //printf("\nz is: %f,  sin of z: %f, cos of z: %f\n", alpha, get_sin(alpha, trig_table), get_cos(alpha, trig_table));
    }
}
