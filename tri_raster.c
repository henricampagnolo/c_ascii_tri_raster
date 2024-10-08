#include "tri_raster.h"

#define WIDTH 240
#define HEIGHT 120 // this should always be half the width
#define SPEED 40 // higher is slower, one is the fastest, you can also change turn_shape in main loop to make faster

float A, B, C;

float zBuffer[WIDTH * HEIGHT];
char buffer[WIDTH * HEIGHT];
int backgroundASCIICode = '.';
float trig_table[TABLE_SIZE + 1];
float fov = 28;
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

    *v_slope = (dz2 * dx1 - dz1 * dx2) / (dx1 * dy2 - dx2 * dy1 + 0.00001);
    *h_slope = (dz1 - dy1 * (*v_slope)) / (dx1 + 0.00001);  
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
        funcs[i].a = (p2.y - p1.y) / ((p2.x - p1.x) + 0.00001);
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
                if (current_z < zBuffer[(HEIGHT - y) *  WIDTH + x])
                {
                    //printf("printing :: x is: %d, y is: %d", x, y);
                    zBuffer[(HEIGHT - y) *  WIDTH + x] = current_z;
                    buffer[(HEIGHT - y) *  WIDTH + x] = symbol;
					//buffer[(HEIGHT - y) *  WIDTH + x] = (char)(current_z - 135 + 30) * 10 + 34;
					//symbol += 1;
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
                if (current_z < zBuffer[(HEIGHT - y) *  WIDTH + x])
                {
                    //printf("printing :: x is: %d, y is: %d", x, y);
                    zBuffer[(HEIGHT - y) *  WIDTH + x] = current_z;
					buffer[(HEIGHT - y) *  WIDTH + x] = symbol;
                    //buffer[(HEIGHT - y) *  WIDTH + x] = (char)(current_z - 135 + 30) * 10 + 34;
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

void define_triangle(t_triangle3 *tri, t_vec3 *points, int p1, int p2, int p3)
{
    tri->points[0] = points + p1;
    tri->points[1] = points + p2;
    tri->points[2] = points + p3;
}

void define_quad(t_triangle3 *tri1, t_triangle3 *tri2, t_vec3 *points, int p1, int p2, int p3, int p4)
{
    tri1->points[0] = points + p1;
    tri1->points[1] = points + p2;
    tri1->points[2] = points + p3;
	tri2->points[0] = points + p2;
    tri2->points[1] = points + p3;
    tri2->points[2] = points + p4;
}


void draw_shape(t_shape shape)
{
    for (int i = 0; i < shape.t_size; i++)
	{
		if (shape.texture[i] != '!')
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
    t_vec3 tetra_center;
    tetra_center.x = 0, tetra_center.y = 0, tetra_center.z = 0;

    t_vec3 tetra_points[4];
	assign_vec(tetra_points + 0, 0, 40, 0);
	assign_vec(tetra_points + 1, 32.66, -13.34, -18.86);
	assign_vec(tetra_points + 2, -32.66, -13.34, -18.86);
	assign_vec(tetra_points + 3, 0, -13.34, 37.71);

    t_shape tetrahedron;

    tetrahedron.t_size = 4;
    tetrahedron.p_size = 4;
    tetrahedron.points = tetra_points;
    t_triangle3 tetra_tris[4];
    define_triangle(tetra_tris, tetra_points, 0, 1, 2);
    define_triangle(tetra_tris + 1, tetra_points, 0, 1, 3);
    define_triangle(tetra_tris + 2, tetra_points, 0, 2, 3);
    define_triangle(tetra_tris + 3, tetra_points, 1, 2, 3);
    tetrahedron.tris = tetra_tris;
    tetrahedron.center = tetra_center;
    tetrahedron.texture = "#a/-";
	*/

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

    t_shape shape4;

	t_vec3 shape4_center;
    shape4_center.x = 27.5, shape4_center.y = 20, shape4_center.z = 5;

	shape4.p_size = 18 + 24;
    t_vec3 shape4_points[shape4.p_size];
    assign_vec(shape4_points + 0, 20, 0, 0);
	assign_vec(shape4_points + 1, 30, 0, 0);
	assign_vec(shape4_points + 2, 30, 20, 0);
	assign_vec(shape4_points + 3, 20, 10, 0);
	assign_vec(shape4_points + 4, 10, 20, 0);
	assign_vec(shape4_points + 5, 0, 10, 0);
	assign_vec(shape4_points + 6, 30, 40, 0);
	assign_vec(shape4_points + 7, 20, 40, 0);
	assign_vec(shape4_points + 8, 0, 20, 0);

	assign_vec(shape4_points + 9, 20, 0, 10);
	assign_vec(shape4_points + 10, 30, 0, 10);
	assign_vec(shape4_points + 11, 30, 20, 10);
	assign_vec(shape4_points + 12, 20, 10, 10);
	assign_vec(shape4_points + 13, 10, 20, 10);
	assign_vec(shape4_points + 14, 0, 10, 10);
	assign_vec(shape4_points + 15, 30, 40, 10);
	assign_vec(shape4_points + 16, 20, 40, 10);
	assign_vec(shape4_points + 17, 0, 20, 10);


	assign_vec(shape4_points + 0 + 18, 35, 30, 0);
	assign_vec(shape4_points + 1 + 18, 35, 40, 0);
	assign_vec(shape4_points + 2 + 18, 45, 40, 0);
	assign_vec(shape4_points + 3 + 18, 55, 40, 0);
	assign_vec(shape4_points + 4 + 18, 55, 30, 0);
	assign_vec(shape4_points + 5 + 18, 45, 30, 0);
	assign_vec(shape4_points + 6 + 18, 45, 20, 0);
	assign_vec(shape4_points + 7 + 18, 35, 20, 0);
	assign_vec(shape4_points + 8 + 18, 35, 10, 0);
	assign_vec(shape4_points + 9 + 18, 45, 10, 0);
	assign_vec(shape4_points + 10 + 18, 55, 10, 0);
	assign_vec(shape4_points + 11 + 18, 55, 20, 0);

	assign_vec(shape4_points + 0 + 18 + 12, 35, 30, 10);
	assign_vec(shape4_points + 1 + 18 + 12, 35, 40, 10);
	assign_vec(shape4_points + 2 + 18 + 12, 45, 40, 10);
	assign_vec(shape4_points + 3 + 18 + 12, 55, 40, 10);
	assign_vec(shape4_points + 4 + 18 + 12, 55, 30, 10);
	assign_vec(shape4_points + 5 + 18 + 12, 45, 30, 10);
	assign_vec(shape4_points + 6 + 18 + 12, 45, 20, 10);
	assign_vec(shape4_points + 7 + 18 + 12, 35, 20, 10);
	assign_vec(shape4_points + 8 + 18 + 12, 35, 10, 10);
	assign_vec(shape4_points + 9 + 18 + 12, 45, 10, 10);
	assign_vec(shape4_points + 10 + 18 + 12, 55, 10, 10);
	assign_vec(shape4_points + 11 + 18 + 12, 55, 20, 10);

    shape4.points = shape4_points;

	shape4.t_size = 30 + 12 + 12 * 2;
    t_triangle3 shape4_tris[shape4.t_size];
    define_triangle(shape4_tris + 0, shape4_points, 0, 1, 2);
    define_triangle(shape4_tris + 1, shape4_points, 0, 2, 3);
	define_triangle(shape4_tris + 2, shape4_points, 2, 3, 4);
	define_triangle(shape4_tris + 3, shape4_points, 3, 4, 5);
	define_triangle(shape4_tris + 4, shape4_points, 5, 6, 7);
	define_triangle(shape4_tris + 5, shape4_points, 5, 7, 8);

	define_triangle(shape4_tris + 6, shape4_points, 0 + 9, 1 + 9, 2 + 9);
    define_triangle(shape4_tris + 7, shape4_points, 0 + 9, 2 + 9, 3 + 9);
	define_triangle(shape4_tris + 8, shape4_points, 2 + 9, 3 + 9, 4 + 9);
	define_triangle(shape4_tris + 9, shape4_points, 3 + 9, 4 + 9, 5 + 9);
	define_triangle(shape4_tris + 10, shape4_points, 5 + 9, 6 + 9, 7 + 9);
	define_triangle(shape4_tris + 11, shape4_points, 5 + 9, 7 + 9, 8 + 9);
	
	define_quad(shape4_tris + 12, shape4_tris + 13, shape4_points, 0, 1, 0 + 9, 1 + 9);
	define_quad(shape4_tris + 14, shape4_tris + 15, shape4_points, 1, 2, 1 + 9, 2 + 9);
	define_quad(shape4_tris + 16, shape4_tris + 17, shape4_points, 2, 4, 2 + 9, 4 + 9);
	define_quad(shape4_tris + 18, shape4_tris + 19, shape4_points, 4, 6, 4 + 9, 6 + 9);
	define_quad(shape4_tris + 20, shape4_tris + 21, shape4_points, 6, 7, 6 + 9, 7 + 9);
	define_quad(shape4_tris + 22, shape4_tris + 23, shape4_points, 7, 8, 7 + 9, 8 + 9);
	define_quad(shape4_tris + 24, shape4_tris + 25, shape4_points, 8, 5, 8 + 9, 5 + 9);
	define_quad(shape4_tris + 26, shape4_tris + 27, shape4_points, 5, 3, 5 + 9, 3 + 9);
	define_quad(shape4_tris + 28, shape4_tris + 29, shape4_points, 3, 0, 3 + 9, 0 + 9);

	define_triangle(shape4_tris + 0 + 30, shape4_points + 18, 0, 1, 2);
    define_triangle(shape4_tris + 1 + 30, shape4_points + 18, 2, 3, 5);
	define_triangle(shape4_tris + 2 + 30, shape4_points + 18, 3, 4, 7);
	define_triangle(shape4_tris + 3 + 30, shape4_points + 18, 4, 7, 8);
	define_triangle(shape4_tris + 4 + 30, shape4_points + 18, 6, 8, 9);
	define_triangle(shape4_tris + 5 + 30, shape4_points + 18, 9, 10, 11);

	define_triangle(shape4_tris + 0 + 36, shape4_points + 18 + 12, 0, 1, 2);
    define_triangle(shape4_tris + 1 + 36, shape4_points + 18 + 12, 2, 3, 5);
	define_triangle(shape4_tris + 2 + 36, shape4_points + 18 + 12, 3, 4, 7);
	define_triangle(shape4_tris + 3 + 36, shape4_points + 18 + 12, 4, 7, 8);
	define_triangle(shape4_tris + 4 + 36, shape4_points + 18 + 12, 6, 8, 9);
	define_triangle(shape4_tris + 5 + 36, shape4_points + 18 + 12, 9, 10, 11);

	define_quad(shape4_tris + 12 + 30, shape4_tris + 13 + 30, shape4_points + 18, 0 , 1 , 0  + 12, 1  + 12);
	define_quad(shape4_tris + 14 + 30, shape4_tris + 15 + 30, shape4_points + 18, 1 , 3 , 1  + 12, 3  + 12);
	define_quad(shape4_tris + 16 + 30, shape4_tris + 17 + 30, shape4_points + 18, 3 , 4 , 3  + 12, 4  + 12);
	define_quad(shape4_tris + 18 + 30, shape4_tris + 19 + 30, shape4_points + 18, 4 , 6 , 4  + 12, 6  + 12);
	define_quad(shape4_tris + 20 + 30, shape4_tris + 21 + 30, shape4_points + 18, 6 , 9 , 6  + 12, 9  + 12);
	define_quad(shape4_tris + 22 + 30, shape4_tris + 23 + 30, shape4_points + 18, 9 , 11, 9  + 12, 11 + 12);
	define_quad(shape4_tris + 24 + 30, shape4_tris + 25 + 30, shape4_points + 18, 11, 10, 11 + 12, 10 + 12);
	define_quad(shape4_tris + 26 + 30, shape4_tris + 27 + 30, shape4_points + 18, 10, 8 , 10 + 12, 8  + 12);
	define_quad(shape4_tris + 28 + 30, shape4_tris + 29 + 30, shape4_points + 18, 8 , 7 , 8  + 12, 7  + 12);
	define_quad(shape4_tris + 30 + 30, shape4_tris + 31 + 30, shape4_points + 18, 7 , 5 , 7  + 12, 5  + 12);
	define_quad(shape4_tris + 32 + 30, shape4_tris + 33 + 30, shape4_points + 18, 5 , 2 , 5  + 12, 2  + 12);
	define_quad(shape4_tris + 34 + 30, shape4_tris + 35 + 30, shape4_points + 18, 2 , 0 , 2  + 12, 0  + 12);

    shape4.tris = shape4_tris;
    shape4.center = shape4_center;
    //shape4.texture = "444444444444--||--\\\\--\\\\||--||222222222222||--||\\\\||\\\\||--||\\\\||\\\\";
	//shape4.texture = "444444DDDDDD!!!!!!!!!!!!!!!!!!222222BBBBBB!!!!!!!!!!!!!!!!!!!!!!!!";
	//shape4.texture = "!!!!!!!!!!!!--||--\\\\--\\\\||--||!!!!!!!!!!!!||--||\\\\||\\\\||--||\\\\||\\\\";
	shape4.texture = "444444444444!!||!!\\\\!!\\\\||!!||222222222222||!!||\\\\||\\\\||!!||\\\\||\\\\";
	//shape4.texture = "!!!444!!!---!!!!!!!!!!\\\\||!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	/*
	*/

    printf("\x1b[2J");

	alpha = 0;
    beta = 0;
    int i = 0;

    t_vec3 new_pos;
    new_pos.x = -27.5, new_pos.y = -20, new_pos.z = 130;

    move_shape(&shape4, new_pos);

	turn_shape(shape4, 0, 0, 0.0008);

    while (i > -1)
    {
        i += 1;

        memset(buffer, backgroundASCIICode, WIDTH * HEIGHT);
        memset(zBuffer, 'a', WIDTH * HEIGHT * 4);

        draw_shape(shape4);

		if(i % SPEED == 0)
		{
			printf("\x1b[H");
			for (int k = 0; k < WIDTH * HEIGHT; k++) 
			{
				putchar(k % WIDTH ? buffer[k] : 10);
				//putchar(k % WIDTH ? ' ' : '\t');
			}
			/*
			*/
			
			turn_shape(shape4, 0, 0.008, 0);
			alpha+=0.0005;
			beta+=0.0005;
		}

        //printf("\nz is: %f,  sin of z: %f, cos of z: %f\n", alpha, get_sin(alpha, trig_table), get_cos(alpha, trig_table));
    }
}
