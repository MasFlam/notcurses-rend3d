#ifndef REND3D_REND3D_H
#define REND3D_REND3D_H

#include <notcurses/notcurses.h>

struct rend3d;

struct r3d_vertex {
	double x, y, z;
};

struct r3d_edge {
	struct r3d_vertex a, b;
};

struct r3d_obj {
	double posx, posy, posz;
	double rotx, roty, rotz;
	double sclx, scly, sclz;
	size_t vertcount;
	struct r3d_vertex *vertices;
	size_t edgecount;
	struct r3d_edge *edges;
};

struct r3d_objref {
	struct r3d_obj obj;
#ifdef REND3D_INTERNAL
	struct r3d_objref *next;
#else
	void *none_of_your_god_damn_business;
#endif
};

struct rend3d *rend3d_create(struct ncplane *drawp);
void rend3d_destroy(struct rend3d *r);

struct r3d_objref *rend3d_add_object(struct rend3d *r, const struct r3d_obj *obj);

void rend3d_render(struct rend3d *r);

#endif
