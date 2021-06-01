#include <inttypes.h>
#include <math.h>
#include <notcurses/notcurses.h>

#include "vecmat.h"

// Swap the doubles a with b (for the line drawing algo)
#define SWAPDBL(a, b) do { double t_m_p = a; a = b; b = t_m_p; } while (0)

struct vertex {
	double x, y, z;
};

struct edge {
	struct vertex a, b;
};

struct obj {
	double posx, posy, posz;
	double rotx, roty, rotz;
	size_t vertcount;
	struct vertex *vertices;
	size_t edgecount;
	struct edge *edges;
};

struct g {
	struct notcurses *nc;
	struct ncplane *stdp, *drawp;
	int termw, termh;
	int wpx, hpx;
	uint8_t *emptybuf, *drawbuf;
};

static void draw_line(double x0, double y0, double x1, double y1);
static inline double fpart(double x);
static void init();
static void render();
static void render_obj(const struct obj *o);
static inline double rfpart(double x);
static inline void set_pixel(int x, int y, double intensity);
static inline void set_pixel_or_more(int x, int y, double intensity);

struct g g;

void
draw_line(double x0, double y0, double x1, double y1)
{
	// Use Wu's algorithm to draw a line
	bool steep = fabs(y1 - y0) > fabs(x1 - x0); // ... > abs(i+1 - i)
	if (steep) {
		SWAPDBL(x0, y0);
		SWAPDBL(x1, y1);
	}
	if (x0 > x1) {
		SWAPDBL(x0, x1);
		SWAPDBL(y0, y1);
	}
	double dx = x1 - x0;
	double dy = y1 - y0;
	double gradient = dy / dx;
	if (dx == 0) gradient = 1;
	// First endpoint of the line
	double xend = round(x0);
	double yend = y0 + gradient * (xend - x0);
	double xgap = rfpart(x0 + 0.5);
	int xpxl1 = xend;
	int ypxl1 = floor(yend);
	if (steep) {
		set_pixel_or_more(ypxl1, xpxl1, rfpart(yend) * xgap);
		set_pixel_or_more(ypxl1+1, xpxl1, fpart(yend) * xgap);
	} else {
		set_pixel_or_more(xpxl1, ypxl1, rfpart(yend) * xgap);
		set_pixel_or_more(xpxl1, ypxl1+1, fpart(yend) * xgap);
	}
	double intery = yend + gradient;
	// Second endpoint of the line
	xend = round(x1);
	yend = y1 + gradient * (xend - x1);
	xgap = fpart(x1 + 0.5);
	int xpxl2 = xend;
	int ypxl2 = floor(yend);
	if (steep) {
		set_pixel_or_more(ypxl2, xpxl2, rfpart(yend) * xgap);
		set_pixel_or_more(ypxl2+1, xpxl2, fpart(yend) * xgap);
	} else {
		set_pixel_or_more(xpxl2, ypxl2, rfpart(yend) * xgap);
		set_pixel_or_more(xpxl2, ypxl2+1, fpart(yend) * xgap);
	}
	// Main loop
	if (steep) {
		for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
			set_pixel_or_more(floor(intery), x, rfpart(intery));
			set_pixel_or_more(floor(intery)+1, x, fpart(intery));
			intery += gradient;
		}
	} else {
		for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
			set_pixel_or_more(x, floor(intery), rfpart(intery));
			set_pixel_or_more(x, floor(intery)+1, fpart(intery));
			intery += gradient;
		}
	}
}

double
fpart(double x)
{
	return x - floor(x);
}

void
init()
{
	g.nc = notcurses_core_init(&(struct notcurses_options) {
		.flags = NCOPTION_SUPPRESS_BANNERS
	}, stdout);
	g.stdp = notcurses_stddim_yx(g.nc, &g.termh, &g.termw);
	if (notcurses_check_pixel_support(g.nc) < 1) {
		notcurses_stop(g.nc);
		fputs("No pixel support!", stderr);
		exit(2);
	}
	// We need this because Notcurses doesn't allow pixel blitting over the standard plane.
	g.drawp = ncplane_create(g.stdp, &(struct ncplane_options) {
		.x = 0, .y = 1,
		.rows = g.termh - 1,
		.cols = g.termw
	});
	ncplane_pixelgeom(g.drawp, NULL, NULL, NULL, NULL, &g.hpx, &g.wpx);
	uint8_t *buf = malloc(g.wpx * g.hpx * 4);
	for (int i = 0; i < g.wpx * g.hpx; ++i) {
		buf[4*i + 0] = 0;
		buf[4*i + 1] = 0;
		buf[4*i + 2] = 0;
		buf[4*i + 3] = 255;
	}
	g.emptybuf = buf;
	g.drawbuf = malloc(g.wpx * g.hpx * 4);
	memcpy(g.drawbuf, g.emptybuf, g.wpx * g.hpx * 4);
}

void
render()
{
	struct ncvisual *ncv = ncvisual_from_rgba(g.drawbuf, g.hpx, g.wpx * 4, g.wpx);
	ncvisual_render(g.nc, ncv, &(struct ncvisual_options) {
		.n = g.drawp,
		.x = 0, .y = 0,
		.scaling = NCSCALE_NONE,
		.blitter = NCBLIT_PIXEL
	});
	notcurses_render(g.nc);
	ncvisual_destroy(ncv);
	memcpy(g.drawbuf, g.emptybuf, g.wpx * g.hpx * 4);
}

void
render_obj(const struct obj *o)
{
	// 1. Matrix preparation
	// Object translation matrix
	mat_t Tobj = {{
		1, 0, 0, o->posx,
		0, 1, 0, o->posy,
		0, 0, 1, o->posz,
		0, 0, 0, 1,
	}};
	double s, c;
	c = cos(o->rotx);
	s = sin(o->rotx);
	// X axe rotation matrix
	mat_t Rx = {{
		1, 0, 0, 0,
		0, c, -s, 0,
		0, s, c, 0,
		0, 0, 0, 1,
	}};
	c = cos(o->roty);
	s = sin(o->roty);
	// Y axe rotation matrix
	mat_t Ry = {{
		c, 0, s, 0,
		0, 1, 0, 0,
		-s, 0, c, 0,
		0, 0, 0, 1,
	}};
	c = cos(o->rotz);
	s = sin(o->rotz);
	// Z axe rotation matrix
	mat_t Rz = {{
		c, -s, 0, 0,
		s, c, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1,
	}};
	// Multiply all the transformation matrices all vertices share into one
	// transformation matrix. (matrix multiplication is associative)
	mat_t transformation_mat = IDENTITY_MAT;
	mul_mat_mat(&Rx, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Ry, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Rz, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Tobj, &transformation_mat, &transformation_mat);
	// 2. Actual rendering
	// 2.1. Lone vertices
	for (size_t i = 0; i < o->vertcount; ++i) {
		// Transform vertex
		vec_t v = { o->vertices[i].x, o->vertices[i].y, o->vertices[i].z, 1 };
		mul_mat_vec(&transformation_mat, &v, &v);
		// Project vertex by multiplying by the projection matrix
		vec_t final;
		mul_mat_vec(&ORTHO_MAT, &v, &final);
		// Finally fill in the pixel, converting the [-1, 1] normalized coordinates
		// to screen coordinates, aka pixels.
		double x = (final.x + 1.0)/2.0 * (g.wpx-1);
		double y = (-final.y + 1.0)/2.0 * (g.hpx-1);
		set_pixel(x, y, 1);
	}
	// 2.2. Edges
	for (size_t i = 0; i < o->edgecount; ++i) {
		// Transform edge ends
		vec_t v0 = { o->edges[i].a.x, o->edges[i].a.y, o->edges[i].a.z, 1 };
		vec_t v1 = { o->edges[i].b.x, o->edges[i].b.y, o->edges[i].b.z, 1 };
		mul_mat_vec(&transformation_mat, &v0, &v0);
		mul_mat_vec(&transformation_mat, &v1, &v1);
		// Project edge ends
		mul_mat_vec(&ORTHO_MAT, &v0, &v0);
		mul_mat_vec(&ORTHO_MAT, &v1, &v1);
		// Draw the line between the pixels where the final edges end up
		double x0 = (v0.x + 1.0)/2.0 * (g.wpx - 1);
		double y0 = (-v0.y + 1.0)/2.0 * (g.hpx - 1);
		double x1 = (v1.x + 1.0)/2.0 * (g.wpx - 1);
		double y1 = (-v1.y + 1.0)/2.0 * (g.hpx - 1);
		draw_line(x0, y0, x1, y1);
	}
}

double
rfpart(double x)
{
	return 1 - fpart(x);
}

void
set_pixel(int x, int y, double intensity)
{
	// Fill the pixel with an appropriate shade of grey
	int i = y*g.wpx*4 + x*4;
	g.drawbuf[i + 0] = intensity * 255;
	g.drawbuf[i + 1] = intensity * 255;
	g.drawbuf[i + 2] = intensity * 255;
}

void
set_pixel_or_more(int x, int y, double intensity)
{
	int i = y*g.wpx*4 + x*4;
	if (g.drawbuf[i + 0] < intensity * 255) g.drawbuf[i + 0] = intensity * 255;
	if (g.drawbuf[i + 1] < intensity * 255) g.drawbuf[i + 1] = intensity * 255;
	if (g.drawbuf[i + 2] < intensity * 255) g.drawbuf[i + 2] = intensity * 255;
}

int
main()
{
	init();
	render_obj(&(struct obj) {
		.posx = 0,
		.posy = 0,
		.posz = 0,
		.rotx = 0.1 * M_PI,
		.roty = 0.1 * M_PI,
		.rotz = 0,
		.vertcount = 1,
		.vertices = (struct vertex[]){
			{ 0.7, 0.7, 0 },
		},
		.edgecount = 12,
		.edges = (struct edge[]) {
			// Front face rim
			{{ -0.5, -0.5, 0.5 }, { -0.5,  0.5, 0.5 }},
			{{ -0.5,  0.5, 0.5 }, {  0.5,  0.5, 0.5 }},
			{{  0.5,  0.5, 0.5 }, {  0.5, -0.5, 0.5 }},
			{{  0.5, -0.5, 0.5 }, { -0.5, -0.5, 0.5 }},
			// Back face rim
			{{ -0.5, -0.5, -0.5 }, { -0.5,  0.5, -0.5 }},
			{{ -0.5,  0.5, -0.5 }, {  0.5,  0.5, -0.5 }},
			{{  0.5,  0.5, -0.5 }, {  0.5, -0.5, -0.5 }},
			{{  0.5, -0.5, -0.5 }, { -0.5, -0.5, -0.5 }},
			// Side edges
			{{  0.5,  0.5, -0.5 }, {  0.5,  0.5, 0.5 }},
			{{  0.5, -0.5, -0.5 }, {  0.5, -0.5, 0.5 }},
			{{ -0.5,  0.5, -0.5 }, { -0.5,  0.5, 0.5 }},
			{{ -0.5, -0.5, -0.5 }, { -0.5, -0.5, 0.5 }},
		}
	});
	render();
	notcurses_getc_blocking(g.nc, NULL);
	notcurses_stop(g.nc);
	free(g.emptybuf);
	free(g.drawbuf);
	return 0;
}
