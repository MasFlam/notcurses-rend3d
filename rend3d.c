/* Rend3D - Terminal 3D rendering library for use with Notcurses
 *
 * Copyright (C) 2021 Łukasz "MasFlam" Drukała
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see https://www.gnu.org/licenses/.
 */
#define REND3D_INTERNAL
#include <math.h>
#include <string.h>
#include <notcurses/notcurses.h>
#include "vecmat.h"

#include "rend3d.h"

// Swap the doubles a with b (for the line drawing algo)
#define SWAPDBL(a, b) do { double t_m_p = a; a = b; b = t_m_p; } while (0)

struct rend3d {
	struct rend3d_options opts;
	struct notcurses *nc;
	struct ncplane *drawp;
	int w, h, wpx, hpx;
	uint8_t *emptybuf, *drawbuf;
	double camx, camy, camz;
	double camrx, camry, camrz;
	int objcount;
	struct r3d_objref *objs;
};

static void draw_line(struct rend3d *r, double x0, double y0, double x1, double y1);
static int drawp_resize_cb(struct ncplane *drawp);
static inline double fpart(double x);
static inline void render_obj(struct rend3d *r, const struct r3d_obj *o, const mat_t *projmat, const mat_t *viewmat);
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
render_obj(struct rend3d *r, const struct r3d_obj *o, const mat_t *projmat, const mat_t *viewmat)
{
	// Prepare transformation matrices
	mat_t obj_translation = {{
		1, 0, 0, o->posx,
		0, 1, 0, o->posy,
		0, 0, 1, o->posz,
		0, 0, 0, 1,
	}};
	double s, c;
	c = cos(o->rotx);
	s = sin(o->rotx);
	mat_t rotation_x = {
		1, 0,  0, 0,
		0, c, -s, 0,
		0, s,  c, 0,
		0, 0,  0, 1,
	};
	c = cos(o->roty);
	s = sin(o->roty);
	mat_t rotation_y = {{
		 c, 0, s, 0,
		 0, 1, 0, 0,
		-s, 0, c, 0,
		 0, 0, 0, 1,
	}};
	c = cos(o->rotz);
	s = sin(o->rotz);
	mat_t rotation_z = {{
		c, -s, 0, 0,
		s,  c, 0, 0,
		0,  0, 1, 0,
		0,  0, 0, 1,
	}};
	mat_t scalemat = {{
		o->sclx, 0,       0,       0,
		0,       o->scly, 0,       0,
		0,       0,       o->sclz, 0,
		0,       0,       0,       1,
	}};
	
	// Combine transformations into one matrix called the model matrix
	mat_t modelmat = IDENTITY_MAT;
	mul_mat_mat(&obj_translation, &modelmat, &modelmat);
	mul_mat_mat(&rotation_z, &modelmat, &modelmat);
	mul_mat_mat(&rotation_y, &modelmat, &modelmat);
	mul_mat_mat(&rotation_x, &modelmat, &modelmat);
	mul_mat_mat(&scalemat, &modelmat, &modelmat);
	
	// Combine Model, View and Perspective into one matrix
	mat_t mvp_mat = IDENTITY_MAT;
	mul_mat_mat(projmat, &mvp_mat, &mvp_mat);
	mul_mat_mat(viewmat, &mvp_mat, &mvp_mat);
	mul_mat_mat(&modelmat, &mvp_mat, &mvp_mat);
	
	for (size_t i = 0; i < o->vertcount; ++i) {
		vec_t v = { o->vertices[i].x, o->vertices[i].y, o->vertices[i].z, .w = 1 };
		mul_mat_vec(&mvp_mat, &v, &v);
		
		if (r->opts.projtype == REND3D_PROJTYPE_PERSP) {
			// Clip on the Z axe
			if (v.w < r->opts.projopts.persp.near) {
				continue;
			}
		}
		
		v.x /= v.w;
		v.y /= v.w;
		v.z /= v.w;
		
		double x = ( v.x + 1.0) / 2.0 * (r->wpx-1);
		double y = (-v.y + 1.0) / 2.0 * (r->hpx-1);
		set_pixel(r, x, y, 1);
	}
	
	for (size_t i = 0; i < o->edgecount; ++i) {
		vec_t v0 = { o->edges[i].a.x, o->edges[i].a.y, o->edges[i].a.z, .w = 1 };
		mul_mat_vec(&mvp_mat, &v0, &v0);
		
		vec_t v1 = { o->edges[i].b.x, o->edges[i].b.y, o->edges[i].b.z, .w = 1 };
		mul_mat_vec(&mvp_mat, &v1, &v1);
		
		if (r->opts.projtype == REND3D_PROJTYPE_PERSP) {
			// Clip on the Z axe
			double near = r->opts.projopts.persp.near;
			double far = r->opts.projopts.persp.far;
			if (v0.w >= near && v1.w >= near) {
				// Render normally
			} else if (v0.w < near && v1.w < near) {
				continue;
			} else if (v0.w >= near && v1.w < near) {
				// https://stackoverflow.com/a/20180585/13694119
				double n = (v0.w - near) / (v0.w - v1.w);
				double xc = n * v0.x + (1-n) * v1.x;
				double yc = n * v0.y + (1-n) * v1.y;
				double zc = n * v0.z + (1-n) * v1.z;
				double wc = near;
				v1 = (vec_t) { xc, yc, zc, wc };
			} else if (v0.w < near && v1.w >= near) {
				double n = (v1.w - near) / (v1.w - v0.w);
				double xc = n * v1.x + (1-n) * v0.x;
				double yc = n * v1.y + (1-n) * v0.y;
				double zc = n * v1.z + (1-n) * v0.z;
				double wc = near;
				v0 = (vec_t) { xc, yc, zc, wc };
			}
		}
		
		v0.x /= v0.w;
		v0.y /= v0.w;
		v0.z /= v0.w;
		v1.x /= v1.w;
		v1.y /= v1.w;
		v1.z /= v1.w;
		
		double x0 = ( v0.x + 1.0) / 2.0 * (r->wpx-1);
		double y0 = (-v0.y + 1.0) / 2.0 * (r->hpx-1);
		double x1 = ( v1.x + 1.0) / 2.0 * (r->wpx-1);
		double y1 = (-v1.y + 1.0) / 2.0 * (r->hpx-1);
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
rend3d_create(struct ncplane *drawp, const struct rend3d_options *opts)
{
	struct rend3d_options defopts = {
		.projtype = REND3D_PROJTYPE_PERSP,
		.projopts.persp.fovx = M_PI / 4,
		.projopts.persp.fovy = M_PI / 4,
		.projopts.persp.near = 0.01,
		.projopts.persp.far = 1
	};
	if (!opts) {
		opts = &defopts;
	}
	struct rend3d *r = malloc(sizeof(struct rend3d));
	memcpy(&r->opts, opts, sizeof(struct rend3d_options));
	r->nc = ncplane_notcurses(drawp);
	r->drawp = drawp;
	r->camx = r->camy = r->camz = 0;
	r->camrx = r->camry = r->camrz = 0;
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
rend3d_cam_get_pos(struct rend3d *r, double *camx, double *camy, double *camz)
{
	*camx = r->camx;
	*camy = r->camy;
	*camz = r->camz;
}

void
rend3d_cam_set_pos(struct rend3d *r, double camx, double camy, double camz)
{
	r->camx = camx;
	r->camy = camy;
	r->camz = camz;
}

void
rend3d_cam_set_pos_x(struct rend3d *r, double camx)
{
	r->camx = camx;
}

void
rend3d_cam_set_pos_y(struct rend3d *r, double camy)
{
	r->camy = camy;
}

void
rend3d_cam_set_pos_z(struct rend3d *r, double camz)
{
	r->camz = camz;
}

void
rend3d_render(struct rend3d *r)
{
	double ar = (double) r->wpx / (double) r->hpx; // aspect ratio
	
	// Prepare projection matrix
	mat_t projmat;
	switch (r->opts.projtype) {
	case REND3D_PROJTYPE_ORTHO: {
		projmat = (mat_t) {{
			1/ar, 0, 0, 0,
			0,    1, 0, 0,
			0,    0, 1, 0,
			0,    0, 0, 1,
		}};
	} break;
	case REND3D_PROJTYPE_PERSP: {
		double near = r->opts.projopts.persp.near;
		double far = r->opts.projopts.persp.far;
		double fovx = r->opts.projopts.persp.fovx;
		double fovy = r->opts.projopts.persp.fovy;
		projmat = (mat_t) {{
			1/(tan(fovx/2)*ar), 0,              0,                      0,
			0,                  1/tan(fovy/2),  0,                      0,
			0,                  0,             -(far+near)/(far-near), -2*near*far/(far-near),
			0,                  0,             -1,                      0,
		}};
	} break;
	}
	
	// Prepare view matrix
	mat_t viewmat_translation = {{
		1, 0, 0, -r->camx,
		0, 1, 0, -r->camy,
		0, 0, 1, -r->camz,
		0, 0, 0,  1,
	}};
	double s, c;
	c = cos(-r->camrx);
	s = sin(-r->camrx);
	mat_t viewmat_rotx = {{
		1, 0,  0, 0,
		0, c, -s, 0,
		0, s,  c, 0,
		0, 0,  0, 1,
	}};
	c = cos(-r->camry);
	s = sin(-r->camry);
	mat_t viewmat_roty = {{
		 c, 0, s, 0,
		 0, 1, 0, 0,
		-s, 0, c, 0,
		 0, 0, 0, 1,
	}};
	c = cos(-r->camrz);
	s = sin(-r->camrz);
	mat_t viewmat_rotz = {{
		c, -s, 0, 0,
		s,  c, 0, 0,
		0,  0, 1, 0,
		0,  0, 0, 1,
	}};
	
	// Combine camera transformations
	mat_t viewmat = IDENTITY_MAT;
	mul_mat_mat(&viewmat_rotx, &viewmat, &viewmat);
	mul_mat_mat(&viewmat_roty, &viewmat, &viewmat);
	mul_mat_mat(&viewmat_rotz, &viewmat, &viewmat);
	mul_mat_mat(&viewmat_translation, &viewmat, &viewmat);
	
	struct r3d_objref *ref = r->objs;
	while (ref) {
		render_obj(r, &ref->obj, &projmat, &viewmat);
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
