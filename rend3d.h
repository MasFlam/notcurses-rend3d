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
#ifndef REND3D_REND3D_H
#define REND3D_REND3D_H

#include <notcurses/notcurses.h>

struct rend3d;

enum rend3d_projtype {
	REND3D_PROJTYPE_ORTHO,
	REND3D_PROJTYPE_PERSP,
};

struct rend3d_options {
	enum rend3d_projtype projtype;
	union {
		struct {
			double fovx, fovy;
			double near, far;
		} persp;
	} projopts;
};

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

struct rend3d *rend3d_create(struct ncplane *drawp, const struct rend3d_options *opts);
void rend3d_destroy(struct rend3d *r);

struct r3d_objref *rend3d_add_object(struct rend3d *r, const struct r3d_obj *obj);

void rend3d_render(struct rend3d *r);

#endif
