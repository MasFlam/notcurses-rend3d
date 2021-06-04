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
#include <string.h>
#include "vecmat.h"

void
mul_mat_vec(const mat_t *m, const vec_t *v, vec_t *out)
{
	vec_t v2;
	v2.x = v->x * m->m[0]  + v->y * m->m[1]  + v->z * m->m[2]  + v->w * m->m[3];
	v2.y = v->x * m->m[4]  + v->y * m->m[5]  + v->z * m->m[6]  + v->w * m->m[7];
	v2.z = v->x * m->m[8]  + v->y * m->m[9]  + v->z * m->m[10] + v->w * m->m[11];
	v2.w = v->x * m->m[12] + v->y * m->m[13] + v->z * m->m[14] + v->w * m->m[15];
	memcpy(out, &v2, sizeof(vec_t));
}

void
mul_mat_mat(const mat_t *a, const mat_t *b, mat_t *out)
{
	mat_t m2;
	for (int i = 0; i < 4; ++i) {
		m2.m[0 + i]  = a->m[0 + i] * b->m[0]  + a->m[4 + i] * b->m[1]  + a->m[8 + i] * b->m[2]  + a->m[12 + i] * b->m[3];
		m2.m[4 + i]  = a->m[0 + i] * b->m[4]  + a->m[4 + i] * b->m[5]  + a->m[8 + i] * b->m[6]  + a->m[12 + i] * b->m[7];
		m2.m[8 + i]  = a->m[0 + i] * b->m[8]  + a->m[4 + i] * b->m[9]  + a->m[8 + i] * b->m[10] + a->m[12 + i] * b->m[11];
		m2.m[12 + i] = a->m[0 + i] * b->m[12] + a->m[4 + i] * b->m[13] + a->m[8 + i] * b->m[14] + a->m[12 + i] * b->m[15];
	}
	memcpy(out, &m2, sizeof(mat_t));
}
