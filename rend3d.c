#define REND3D_INTERNAL
#include <math.h>
#include <string.h>
#include <notcurses/notcurses.h>
#include "vecmat.h"

#include "rend3d.h"

// Swap the doubles a with b (for the line drawing algo)
#define SWAPDBL(a, b) do { double t_m_p = a; a = b; b = t_m_p; } while (0)

struct rend3d {
	struct notcurses *nc;
	struct ncplane *drawp;
	int w, h, wpx, hpx;
	uint8_t *emptybuf, *drawbuf;
	int objcount;
	struct r3d_objref *objs;
};

static void draw_line(struct rend3d *r, double x0, double y0, double x1, double y1);
static int drawp_resize_cb(struct ncplane *drawp);
static inline double fpart(double x);
static inline void render_obj(struct rend3d *r, const struct r3d_obj *o);
static inline double rfpart(double x);
static inline void set_pixel(struct rend3d *r, int x, int y, double intensity);
static inline void set_pixel_or_more(struct rend3d *r, int x, int y, double intensity);

void
draw_line(struct rend3d *r, double x0, double y0, double x1, double y1)
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
		set_pixel_or_more(r, ypxl1, xpxl1, rfpart(yend) * xgap);
		set_pixel_or_more(r, ypxl1+1, xpxl1, fpart(yend) * xgap);
	} else {
		set_pixel_or_more(r, xpxl1, ypxl1, rfpart(yend) * xgap);
		set_pixel_or_more(r, xpxl1, ypxl1+1, fpart(yend) * xgap);
	}
	double intery = yend + gradient;
	// Second endpoint of the line
	xend = round(x1);
	yend = y1 + gradient * (xend - x1);
	xgap = fpart(x1 + 0.5);
	int xpxl2 = xend;
	int ypxl2 = floor(yend);
	if (steep) {
		set_pixel_or_more(r, ypxl2, xpxl2, rfpart(yend) * xgap);
		set_pixel_or_more(r, ypxl2+1, xpxl2, fpart(yend) * xgap);
	} else {
		set_pixel_or_more(r, xpxl2, ypxl2, rfpart(yend) * xgap);
		set_pixel_or_more(r, xpxl2, ypxl2+1, fpart(yend) * xgap);
	}
	// Main loop
	if (steep) {
		for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
			set_pixel_or_more(r, floor(intery), x, rfpart(intery));
			set_pixel_or_more(r, floor(intery)+1, x, fpart(intery));
			intery += gradient;
		}
	} else {
		for (int x = xpxl1 + 1; x <= xpxl2 - 1; ++x) {
			set_pixel_or_more(r, x, floor(intery), rfpart(intery));
			set_pixel_or_more(r, x, floor(intery)+1, fpart(intery));
			intery += gradient;
		}
	}
}

int
drawp_resize_cb(struct ncplane *drawp)
{
	struct rend3d *r = ncplane_userptr(drawp);
	ncplane_resize_maximize(drawp);
	ncplane_pixelgeom(drawp, NULL, NULL, NULL, NULL, &r->hpx, &r->wpx);
	ncplane_dim_yx(drawp, &r->h, &r->w);
	r->emptybuf = realloc(r->emptybuf, r->wpx * r->hpx * 4);
	r->drawbuf = realloc(r->drawbuf, r->wpx * r->hpx * 4);
	for (int i = 0; i < r->wpx * r->hpx; ++i) {
		r->emptybuf[4*i + 0] = 0;
		r->emptybuf[4*i + 1] = 0;
		r->emptybuf[4*i + 2] = 0;
		r->emptybuf[4*i + 3] = 255;
	}
	memcpy(r->drawbuf, r->emptybuf, r->wpx * r->hpx * 4);
	return 0;
}

double
fpart(double x)
{
	return x - floor(x);
}

void
render_obj(struct rend3d *r, const struct r3d_obj *o)
{
	// 1. Matrix preparation
	// Projection matrix
	double ar = (double) r->wpx / (double) r->hpx;
	double near = 0.1;
	double far = 0.9;
	double fov = M_PI/4;
	mat_t projmat = {{
		// Orthographic projection:
		// 1/ar, 0, 0, 0,
		// 0, 1, 0, 0,
		// 0, 0, 0, 0,
		// 0, 0, 0, 0,
		// Perspective projection:
		1/(tan(fov/2)*ar), 0, 0, 0,
		0, 1/tan(fov/2), 0, 0,
		0, 0, -(far+near)/(far-near), -2*near*far/(far-near),
		0, 0, -1, 0,
	}};
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
	// Scale matrix
	mat_t S = {{
		o->sclx, 0, 0, 0,
		0, o->scly, 0, 0,
		0, 0, o->sclz, 0,
		0, 0, 0, 1,
	}};
	// Multiply all the transformation matrices all vertices share into one
	// transformation matrix. (matrix multiplication is associative)
	mat_t transformation_mat = IDENTITY_MAT;
	mul_mat_mat(&Tobj, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Rz, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Ry, &transformation_mat, &transformation_mat);
	mul_mat_mat(&Rx, &transformation_mat, &transformation_mat);
	mul_mat_mat(&S, &transformation_mat, &transformation_mat);
	// 2. Actual rendering
	// 2.1. Lone vertices
	for (size_t i = 0; i < o->vertcount; ++i) {
		// Transform vertex
		vec_t v = { o->vertices[i].x, o->vertices[i].y, o->vertices[i].z, 1 };
		mul_mat_vec(&transformation_mat, &v, &v);
		// Project vertex by multiplying by the projection matrix
		vec_t final;
		mul_mat_vec(&projmat, &v, &final);
		// An (x, y, z, w) vector really means by definition (x/w, y/w, z/w)
		final.x /= final.w;
		final.y /= final.w;
		final.z /= final.w;
		// Finally fill in the pixel, converting the [-1, 1] normalized coordinates
		// to screen coordinates, aka pixels.
		double x = (final.x + 1.0)/2.0 * (r->wpx-1);
		double y = (-final.y + 1.0)/2.0 * (r->hpx-1);
		set_pixel(r, x, y, 1);
	}
	// 2.2. Edges
	for (size_t i = 0; i < o->edgecount; ++i) {
		// Transform edge ends
		vec_t v0 = { o->edges[i].a.x, o->edges[i].a.y, o->edges[i].a.z, 1 };
		vec_t v1 = { o->edges[i].b.x, o->edges[i].b.y, o->edges[i].b.z, 1 };
		mul_mat_vec(&transformation_mat, &v0, &v0);
		mul_mat_vec(&transformation_mat, &v1, &v1);
		// Project edge ends
		mul_mat_vec(&projmat, &v0, &v0);
		mul_mat_vec(&projmat, &v1, &v1);
		// An (x, y, z, w) vector really means by definition (x/w, y/w, z/w)
		v0.x /= v0.w;
		v0.y /= v0.w;
		v0.z /= v0.w;
		v1.x /= v1.w;
		v1.y /= v1.w;
		v1.z /= v1.w;
		// Draw the line between the pixels where the final edges end up
		double x0 = (v0.x + 1.0)/2.0 * (r->wpx - 1);
		double y0 = (-v0.y + 1.0)/2.0 * (r->hpx - 1);
		double x1 = (v1.x + 1.0)/2.0 * (r->wpx - 1);
		double y1 = (-v1.y + 1.0)/2.0 * (r->hpx - 1);
		draw_line(r, x0, y0, x1, y1);
	}
}

double
rfpart(double x)
{
	return 1 - fpart(x);
}

void
set_pixel(struct rend3d *r, int x, int y, double intensity)
{
	// Fill the pixel with an appropriate shade of grey
	if (x < 0 || y < 0 || x >= r->wpx || y >= r->hpx) return;
	int i = y*r->wpx*4 + x*4;
	r->drawbuf[i + 0] = intensity * 255;
	r->drawbuf[i + 1] = intensity * 255;
	r->drawbuf[i + 2] = intensity * 255;
}

void
set_pixel_or_more(struct rend3d *r, int x, int y, double intensity)
{
	if (x < 0 || y < 0 || x >= r->wpx || y >= r->hpx) return;
	int i = y*r->wpx*4 + x*4;
	if (r->drawbuf[i + 0] < intensity * 255) r->drawbuf[i + 0] = intensity * 255;
	if (r->drawbuf[i + 1] < intensity * 255) r->drawbuf[i + 1] = intensity * 255;
	if (r->drawbuf[i + 2] < intensity * 255) r->drawbuf[i + 2] = intensity * 255;
}

struct rend3d *
rend3d_create(struct ncplane *drawp)
{
	struct rend3d *r = malloc(sizeof(struct rend3d));
	r->nc = ncplane_notcurses(drawp);
	r->drawp = drawp;
	r->objcount = 0;
	r->objs = NULL;
	ncplane_set_userptr(drawp, r);
	ncplane_dim_yx(drawp, &r->h, &r->w);
	ncplane_pixelgeom(r->drawp, NULL, NULL, NULL, NULL, &r->hpx, &r->wpx);
	ncplane_set_resizecb(r->drawp, drawp_resize_cb);
	uint8_t *buf = malloc(r->wpx * r->hpx * 4);
	for (int i = 0; i < r->wpx * r->hpx; ++i) {
		buf[4*i + 0] = 0;
		buf[4*i + 1] = 0;
		buf[4*i + 2] = 0;
		buf[4*i + 3] = 255;
	}
	r->emptybuf = buf;
	r->drawbuf = malloc(r->wpx * r->hpx * 4);
	memcpy(r->drawbuf, r->emptybuf, r->wpx * r->hpx * 4);
	return r;
}

void
rend3d_destroy(struct rend3d *r)
{
	if (!r) return;
	struct r3d_objref *tmp, *ref = r->objs;
	while (ref) {
		tmp = ref->next;
		free(ref);
		ref = tmp;
	}
	free(r->emptybuf);
	free(r->drawbuf);
	ncplane_destroy(r->drawp);
	free(r);
}

struct r3d_objref *
rend3d_add_object(struct rend3d *r, const struct r3d_obj *obj)
{
	struct r3d_objref *objref = malloc(sizeof(struct r3d_objref));
	if (r->objcount == 0) {
		r->objs = objref;
		objref->next = NULL;
	} else {
		struct r3d_objref *rf = r->objs;
		r->objs = objref;
		objref->next = rf;
	}
	++r->objcount;
	memcpy(&objref->obj, obj, sizeof(struct r3d_obj));
	objref->obj.vertices = malloc(obj->vertcount * sizeof(struct r3d_vertex));
	objref->obj.edges = malloc(obj->edgecount * sizeof(struct r3d_edge));
	memcpy(objref->obj.vertices, obj->vertices, obj->vertcount * sizeof(struct r3d_vertex));
	memcpy(objref->obj.edges, obj->edges, obj->edgecount * sizeof(struct r3d_edge));
	return objref;
}

void
rend3d_render(struct rend3d *r)
{
	struct r3d_objref *ref = r->objs;
	while (ref) {
		render_obj(r, &ref->obj);
		ref = ref->next;
	}
	struct ncvisual *ncv = ncvisual_from_rgba(r->drawbuf, r->hpx, r->wpx * 4, r->wpx);
	ncvisual_render(r->nc, ncv, &(struct ncvisual_options) {
		.n = r->drawp,
		.x = 0, .y = 0,
		.scaling = NCSCALE_NONE,
		.blitter = NCBLIT_PIXEL
	});
	ncvisual_destroy(ncv);
	memcpy(r->drawbuf, r->emptybuf, r->wpx * r->hpx * 4);
}
