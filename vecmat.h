#ifndef VECMAT_H
#define VECMAT_H

#define IDENTITY_MAT ((mat_t) {{ \
	1, 0, 0, 0, \
	0, 1, 0, 0, \
	0, 0, 1, 0, \
	0, 0, 0, 1, \
}})

#define ORTHO_MAT ((mat_t) {{ \
	1, 0, 0, 0, \
	0, 1, 0, 0, \
	0, 0, 0, 0, \
	0, 0, 0, 0, \
}})

typedef struct vec {
	double x, y, z, w;
} vec_t;

typedef struct mat {
	// 4x4 row-major (for easy initializing)
	double m[16];
} mat_t;

void mul_mat_vec(const mat_t *m, const vec_t *v, vec_t *out);
void mul_mat_mat(const mat_t *a, const mat_t *b, mat_t *out);

#endif
