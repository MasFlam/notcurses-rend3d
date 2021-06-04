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
	g.r3d = rend3d_create(drawp);
}

int
main()
{
	init();
	struct r3d_obj obj = {
		.posx = 0,
		.posy = 0,
		.posz = -1,
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
	struct r3d_objref *ref = rend3d_add_object(g.r3d, &obj);
	while (1) {
		rend3d_render(g.r3d);
		notcurses_render(g.nc);
		ref->obj.rotx += 0.02 * M_PI;
		ref->obj.roty += 0.015 * M_PI;
		ref->obj.rotz += 0.01 * M_PI;
	}
	cleanup();
	return 0;
}
