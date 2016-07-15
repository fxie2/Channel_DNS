#pragma once

// Vertical position
extern double *x, *y, *z;

// Cell length
extern double *h, *dy2h1, *dy2h2, *dy2h3, *dyh1, *dyh2, *dyh3;
extern double *dy, *dy2dy1, *dy2dy2, *dy2dy3, *dydy1, *dydy2, *dydy3;
extern double dx, dz;

// Rank
extern int *im, *ip;
extern int *jm, *jp;
extern int *km, *kp;

// Mesh scale in y direction
extern double scale;

void init_mesh();

void del_mesh();