/* Rend3D Demo - Terminal 3D rendering demo program powered by Notcurses
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
#include <math.h>
#include <notcurses/notcurses.h>
#include "rend3d.h"

struct g {
	struct notcurses *nc;
	struct ncplane *stdp;
	struct rend3d *r3d;
};

static void cleanup();
static void init();

struct g g;

void
cleanup()
{
	rend3d_destroy(g.r3d);
	notcurses_stop(g.nc);
}

void
init()
{
	g.nc = notcurses_core_init(&(struct notcurses_options) {
		.flags = NCOPTION_SUPPRESS_BANNERS
	}, stdout);
	if (notcurses_check_pixel_support(g.nc) < 1) {
		notcurses_stop(g.nc);
		fputs("No pixel support!", stderr);
		exit(2);
	}
	int termw, termh;
	g.stdp = notcurses_stddim_yx(g.nc, &termh, &termw);
	struct ncplane *drawp = ncplane_create(g.stdp, &(struct ncplane_options) {
		.x = 0, .y = 0,
		.rows = termh,
		.cols = termw
	});
	g.r3d = rend3d_create(drawp, NULL);
}

int
main()
{
	init();
	struct r3d_obj obj = {
		.posx = 0,
		.posy = 0,
		.posz = 0,
		.rotx = 0,
		.roty = 0,
		.rotz = 0,
		.sclx = 0.2,
		.scly = 0.2,
		.sclz = 0.2,
		.vertcount = 1,
		.vertices = (struct r3d_vertex[]){
			{ 0.7, 0.7, 0 },
		},
		.edgecount = 12,
		.edges = (struct r3d_edge[]) {
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
	};
	
	static const int N = 3;
	for (int x = -N; x <= N; ++x) {
		obj.posx = x;
		for (int y = -N; y <= N; ++y) {
			obj.posy = y;
			for (int z = -N; z <= N; ++z) {
				obj.posz = z;
				rend3d_add_object(g.r3d, &obj);
			}
		}
	}
	
	rend3d_cam_move(g.r3d, 0.3, 0.4, 0.5);
	
	while (1) {
		rend3d_render(g.r3d);
		notcurses_render(g.nc);
		int i = 0;
		struct r3d_objref *ref = rend3d_get_first_object(g.r3d);
		while (ref) {
			ref->obj.rotx += (i%2 == 0 ? 1 : -1) * 0.02 * M_PI;
			ref->obj.roty += (i%3 == 0 ? 1 : -1) * 0.015 * M_PI;
			ref->obj.rotz += (i%5 == 0 ? 1 : -1) * 0.01 * M_PI;
			++i;
			ref = rend3d_get_next_object(ref);
		}
		rend3d_cam_rotate(g.r3d, 0.007 * M_PI, 0.01 * M_PI, 0.005 * M_PI);
	}
	cleanup();
	return 0;
}
