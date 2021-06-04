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
#ifndef REND3D_VECMAT_H
#define REND3D_VECMAT_H

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
