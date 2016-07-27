/*
	rastertest.c : Just learning.
	This code is a mess, be warned.
	(C)2016 Marisa Kirisame, UnSX Team.
	Released under the GNU GPLv3 (or later).
*/
#define _DEFAULT_SOURCE

#include <omp.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rastertest.h"

const char helptext[9][64] =
{
	"WASD : move horizontally",
	"QE   : move vertically",
	"↑←↓→ : rotate",
	"F    : toggle fps display",
	"X    : toggle fill/wireframe",
	"Z    : toggle auto-rotate",
	"C    : toggle no/back/front face culling",
	"F1   : toggle this text",
	"ESC  : quit"
};

vect_t cubevect[8] =
{
	{ 1.f,-1.f,-1.f,1.f},{ 1.f,-1.f, 1.f,1.f},
	{-1.f,-1.f, 1.f,1.f},{-1.f,-1.f,-1.f,1.f},
	{ 1.f, 1.f,-1.f,1.f},{ 1.f, 1.f, 1.f,1.f},
	{-1.f, 1.f, 1.f,1.f},{-1.f, 1.f,-1.f,1.f}
};
vect_t cubenorm[6] =
{
	{ 0.f,-1.f, 0.f,1.f},{ 0.f, 1.f, 0.f,1.f},
	{ 1.f, 0.f, 0.f,1.f},{ 0.f, 0.f, 1.f,1.f},
	{-1.f, 0.f, 0.f,1.f},{ 0.f, 0.f,-1.f,1.f}
};
coord_t cubecoord[4] =
{
	{4.f,4.f},{0.f,0.f},
	{4.f,0.f},{0.f,4.f}
};
color_t cubecolor[8] =
{
	{0.f,0.f,0.f,1.f},
	{1.f,0.f,0.f,1.f},
	{0.f,1.f,0.f,1.f},
	{1.f,1.f,0.f,1.f},
	{0.f,0.f,1.f,1.f},
	{1.f,0.f,1.f,1.f},
	{0.f,1.f,1.f,1.f},
	{1.f,1.f,1.f,1.f},
};
face_t cubefaces[12] =
{
	{{1,3,0},{0,0,0},{0,1,2},{1,1,1},0},
	{{7,5,4},{1,1,1},{3,2,0},{2,2,2},0},
	{{4,1,0},{2,2,2},{0,1,2},{3,3,3},0},
	{{5,2,1},{3,3,3},{3,2,1},{4,4,4},0},
	{{2,7,3},{4,4,4},{3,2,0},{5,5,5},0},
	{{0,7,4},{5,5,5},{3,2,1},{6,6,6},0},
	{{1,2,3},{0,0,0},{0,3,1},{1,1,1},0},
	{{7,6,5},{1,1,1},{3,1,2},{2,2,2},0},
	{{4,5,1},{2,2,2},{0,3,1},{3,3,3},0},
	{{5,6,2},{3,3,3},{3,0,2},{4,4,4},0},
	{{2,6,7},{4,4,4},{3,1,2},{5,5,5},0},
	{{0,3,7},{5,5,5},{3,0,2},{6,6,6},0}
};
model_t cube =
{
	cubevect,cubenorm,cubecoord,cubecolor,0,cubefaces,8,6,4,8,0,12
};
vect_t cubepos = {0.f,0.f,-4.f,1.f};
vect_t cuberot = {0.f,0.f,0.f,0.f};
mat_t projection, translation, rotation;
buffer_t screen = {0,0,0,640,480,0.01f,10.0f};

pixel_t cubepixels[4] =
{
	{0,0,0,255},{255,255,255,255},
	{255,255,255,255},{0,0,0,255}
};
texture_t cubetex =
{
	cubepixels,2,2,0,TXWRAP,TXWRAP,{0.f,0.f,0.f,1.f}
};

/* FPS helper code */
#define TMAX 64
int ti = 0;
float ts = 0.f, tl[TMAX] = {0.f};
float avg_fps( float nt )
{
	ts = ts-tl[ti]+nt;
	tl[ti] = nt;
	if ( ++ti == TMAX ) ti = 0;
	return ts/TMAX;
}
#define NANOS_SEC 1000000000L
long ticker( void )
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC_RAW,&ts);
	return ts.tv_nsec+ts.tv_sec*NANOS_SEC;
}

/* raster code */

void putpixel( px_t p, color_t c )
{
	if ( (p.x < 0) || (p.x >= screen.width ) || (p.y < 0)
		|| (p.y >= screen.height) )
		return;
	int coord = p.x+p.y*screen.width;
	screen.depth[coord] = p.d;
	screen.color[coord].r = c.r*255.f;
	screen.color[coord].g = c.g*255.f;
	screen.color[coord].b = c.b*255.f;
	screen.color[coord].a = c.a*255.f;
}

void drawline( int x0, int x1, int y0, int y1, color_t c )
{
	int steep = 0;
	if ( abs(x0-x1) < abs(y0-y1) )
	{
		swap(int,x0,y0);
		swap(int,x1,y1);
		steep = 1;
	}
	if ( x0 > x1 )
	{
		swap(int,x0,x1);
		swap(int,y0,y1);
	}
	int dx = x1-x0;
	int dy = y1-y0;
	int der = abs(dy)*2;
	int er = 0;
	int y = y0;
	for ( int x=x0; x<=x1; x++ )
	{
		px_t p = {x,y,0.f,0};
		if ( steep ) swap(int,p.x,p.y);
		putpixel(p,c);
		er += der;
		if ( er > dx )
		{
			y += (y1>y0)?1:-1;
			er -= dx*2;
		}
	}
}

void sampletexture( color_t *o, coord_t c, texture_t t )
{
	int cx = 0, cy = 0;
	if ( t.bu == TXCLAMP ) cx = c.s*t.width;
	else if ( t.bu == TXWRAP ) cx = (c.s-floorf(c.s))*t.width;
	else if ( t.bu == TXMIRROR )
		cx = (1.f-fabs(1.f-fmod(fabs(c.s),2.f)))*t.width;
	else if ( t.bu == TXBORDER )
	{
		cx = c.s*t.width;
		if ( c.s >= 1.f || c.s < 0.f )
		{
			*o = t.border;
			return;
		}
	}
	if ( t.bv == TXCLAMP ) cy = c.t*t.height;
	else if ( t.bv == TXWRAP ) cy = (c.t-floorf(c.t))*t.height;
	else if ( t.bv == TXMIRROR )
		cy = (1.f-fabs(1.f-fmod(fabs(c.t),2.f)))*t.height;
	else if ( t.bv == TXBORDER )
	{
		cy = c.t*t.height;
		if ( c.t >= 1.f || c.t < 0.f )
		{
			*o = t.border;
			return;
		}
	}
	cx = clamp(cx,0,t.width-1);
	cy = clamp(cy,0,t.height-1);
	o->r = t.data[cx+cy*t.width].r/255.f;
	o->g = t.data[cx+cy*t.width].g/255.f;
	o->b = t.data[cx+cy*t.width].b/255.f;
	o->a = t.data[cx+cy*t.width].a/255.f;
}

#define afn(a,b,c) ((c.x-a.x)*(b.y-a.y)-(c.y-a.y)*(b.x-a.x))

void filltriangle( px_t p[3], coord_t tc[3], color_t c[3] )
{
	px_t minb, maxb;
	maxb.x = clamp(max(p[0].x,max(p[1].x,p[2].x)),0,screen.width);
	maxb.y = clamp(max(p[0].y,max(p[1].y,p[2].y)),0,screen.height);
	minb.x = clamp(min(p[0].x,min(p[1].x,p[2].x)),0,screen.width);
	minb.y = clamp(min(p[0].y,min(p[1].y,p[2].y)),0,screen.height);
	if ( (maxb.x == minb.x) || (maxb.y == minb.y) ) return;
	tc[0].s /= p[0].d;
	tc[0].t /= p[0].d;
	tc[1].s /= p[1].d;
	tc[1].t /= p[1].d;
	tc[2].s /= p[2].d;
	tc[2].t /= p[2].d;
	c[0].r /= p[0].d;
	c[0].g /= p[0].d;
	c[0].b /= p[0].d;
	c[0].a /= p[0].d;
	c[1].r /= p[1].d;
	c[1].g /= p[1].d;
	c[1].b /= p[1].d;
	c[1].a /= p[1].d;
	c[2].r /= p[2].d;
	c[2].g /= p[2].d;
	c[2].b /= p[2].d;
	c[2].a /= p[2].d;
	p[0].d = 1.f/p[0].d;
	p[1].d = 1.f/p[1].d;
	p[2].d = 1.f/p[2].d;
	float area = afn(p[0],p[1],p[2]);
	const int this_is_invariant_i_swear = maxb.y;
	const int this_is_also_invariant_i_swear = maxb.x;
	#pragma omp parallel for collapse(2)
	for ( int y=minb.y; y<this_is_invariant_i_swear; y++ )
	{
		for ( int x=minb.x; x<this_is_also_invariant_i_swear; x++ )
		{
			px_t px = {x+0.5f,y+0.5f,0,0};
			vect_t q =
			{
				afn(p[1],p[2],px)/area,
				afn(p[2],p[0],px)/area,
				afn(p[0],p[1],px)/area,
				0.f
			};
			if ( (q.x<0.f) || (q.y<0.f) || (q.z<0.f) ) continue;
			int coord = x+y*screen.width;
			float dep = screen.depth[coord];
			float pd = 1.f/(q.x*p[0].d+q.y*p[1].d+q.z*p[2].d);
			if ( (pd > dep) ) continue;
			color_t pc =
			{
				q.x*c[0].r+q.y*c[1].r+q.z*c[2].r,
				q.x*c[0].g+q.y*c[1].g+q.z*c[2].g,
				q.x*c[0].b+q.y*c[1].b+q.z*c[2].b,
				q.x*c[0].a+q.y*c[1].a+q.z*c[2].a,
			};
			pc.r *= pd;
			pc.g *= pd;
			pc.b *= pd;
			pc.a *= pd;
			coord_t tx =
			{
				q.x*tc[0].s+q.y*tc[1].s+q.z*tc[2].s,
				q.x*tc[0].t+q.y*tc[1].t+q.z*tc[2].t
			};
			tx.s *= pd;
			tx.t *= pd;
			sampletexture(&pc,tx,cubetex);
			px.x = x;
			px.y = y;
			px.d = pd;
			putpixel(px,pc);
		}
	}
}

int culling = 1;

void drawclippedtriangle( tri_t t )
{
	px_t pts[3] =
	{
		{(1.f+t.v[0].x/t.v[0].z)*0.5f*screen.width,
			(1.f+t.v[0].y/t.v[0].z)*0.5f*screen.height,
			t.v[0].z,0.f},
		{(1.f+t.v[1].x/t.v[1].z)*0.5f*screen.width,
			(1.f+t.v[1].y/t.v[1].z)*0.5f*screen.height,
			t.v[1].z,0.f},
		{(1.f+t.v[2].x/t.v[2].z)*0.5f*screen.width,
			(1.f+t.v[2].y/t.v[2].z)*0.5f*screen.height,
			t.v[2].z,0.f}
	};
	vect_t facet, ab, ac;
	ab.x = pts[1].x-pts[0].x;
	ab.y = pts[1].y-pts[0].y;
	ab.z = pts[1].d-pts[0].d;
	ac.x = pts[2].x-pts[0].x;
	ac.y = pts[2].y-pts[0].y;
	ac.z = pts[2].d-pts[0].d;
	cross(&facet,ab,ac);
	if ( (culling == 1) && (facet.z < 0.f) ) return;
	if ( (culling == 2) && facet.z > 0.f ) return;
	color_t cols[3] = {t.c[0],t.c[1],t.c[2]};
	coord_t coords[3] = {t.t[0],t.t[1],t.t[2]};
	filltriangle(pts,coords,cols);
}

void clipanddrawwire( vect_t a, vect_t b )
{
	vect_t p;
	float t;
	a.z *= -1;
	b.z *= -1;
	if ( (a.z < screen.znear) && (b.z < screen.znear) ) return;
	if ( (a.z > screen.zfar) && (b.z > screen.zfar) ) return;
	if ( a.z < screen.znear )
	{
		vsub(&p,a,b);
		t = (screen.znear-b.z)/p.z;
		a.x = b.x+p.x*t;
		a.y = b.y+p.y*t;
		a.z = screen.znear;
	}
	else if ( b.z < screen.znear )
	{
		vsub(&p,b,a);
		t = (screen.znear-a.z)/p.z;
		b.x = a.x+p.x*t;
		b.y = a.y+p.y*t;
		b.z = screen.znear;
	}
	if ( a.z > screen.zfar )
	{
		vsub(&p,a,b);
		t = (screen.zfar-b.z)/p.z;
		a.x = b.x+p.x*t;
		a.y = b.y+p.y*t;
		a.z = screen.zfar;
	}
	else if ( b.z > screen.zfar )
	{
		vsub(&p,b,a);
		t = (screen.zfar-a.z)/p.z;
		b.x = a.x+p.x*t;
		b.y = a.y+p.y*t;
		b.z = screen.zfar;
	}
	a.x = (1.f+a.x/a.z)*0.5f*screen.width;
	a.y = (1.f+a.y/a.z)*0.5f*screen.height;
	b.x = (1.f+b.x/b.z)*0.5f*screen.width;
	b.y = (1.f+b.y/b.z)*0.5f*screen.height;
	color_t wire = {0.f,0.6f,0.f,1.f};
	drawline(a.x,b.x,a.y,b.y,wire);
}

void drawtriangle( tri_t t, unsigned char wireframe )
{
	if ( wireframe )
	{
		for ( int i=0; i<3; i++ ) clipanddrawwire(t.v[i],t.v[(i+1)%3]);
		return;
	}
	t.v[0].z *= -1;
	t.v[1].z *= -1;
	t.v[2].z *= -1;
	if ( (t.v[0].z < screen.znear) && (t.v[1].z < screen.znear)
		&& (t.v[2].z < screen.znear) ) return;
	if ( (t.v[0].z > screen.zfar) && (t.v[1].z > screen.zfar)
		&& (t.v[2].z > screen.zfar) ) return;
	int ff = 0, fb = 0, fc = 0;
	for ( int i=0; i<3; i++ )
	{
		if ( t.v[i].z < screen.znear )
		{
			fb = i;
			fc++;
		}
		else ff = i;
	}
	if ( fc == 2 )
	{
		tri_t newt;
		vect_t p1, p2;
		float t1, t2;
		int a = (ff+1)%3, b = (ff+2)%3, c = ff;
		vsub(&p1,t.v[a],t.v[c]);
		vsub(&p2,t.v[b],t.v[c]);
		t1 = (screen.znear-t.v[c].z)/p1.z;
		t2 = (screen.znear-t.v[c].z)/p2.z;
		newt.m = t.m;
		newt.v[0].x = t.v[c].x+p1.x*t1;
		newt.v[0].y = t.v[c].y+p1.y*t1;
		newt.v[0].z = screen.znear;
		newt.v[0].w = 1.f;
		vlerp(&newt.n[0],t.n[c],t.n[a],t1);
		clerp(&newt.c[0],t.c[c],t.c[a],t1);
		newt.t[0].s = lerp(t.t[c].s,t.t[a].s,t1);
		newt.t[0].t = lerp(t.t[c].t,t.t[a].t,t1);
		newt.v[1].x = t.v[c].x+p2.x*t2;
		newt.v[1].y = t.v[c].y+p2.y*t2;
		newt.v[1].z = screen.znear;
		newt.v[1].w = 1.f;
		vlerp(&newt.n[1],t.n[c],t.n[b],t2);
		clerp(&newt.c[1],t.c[c],t.c[b],t2);
		newt.t[1].s = lerp(t.t[c].s,t.t[b].s,t2);
		newt.t[1].t = lerp(t.t[c].t,t.t[b].t,t2);
		newt.v[2] = t.v[c];
		newt.n[2] = t.n[c];
		newt.c[2] = t.c[c];
		newt.t[2] = t.t[c];
		drawclippedtriangle(newt);
	}
	else if ( fc == 1 )
	{
		tri_t newt1, newt2;
		vect_t p1, p2;
		float t1, t2;
		int a = (fb+1)%3, b = (fb+2)%3, c = fb;
		vsub(&p1,t.v[c],t.v[a]);
		vsub(&p2,t.v[c],t.v[b]);
		t1 = (screen.znear-t.v[a].z)/p1.z;
		t2 = (screen.znear-t.v[b].z)/p2.z;
		newt1.m = t.m;
		newt1.v[0] = t.v[a];
		newt1.n[0] = t.n[a];
		newt1.c[0] = t.c[a];
		newt1.t[0] = t.t[a];
		newt1.v[1] = t.v[b];
		newt1.n[1] = t.n[b];
		newt1.c[1] = t.c[b];
		newt1.t[1] = t.t[b];
		newt1.v[2].x = t.v[a].x+p1.x*t1;
		newt1.v[2].y = t.v[a].y+p1.y*t1;
		newt1.v[2].z = screen.znear;
		newt1.v[2].w = 1.f;
		vlerp(&newt1.n[2],t.n[a],t.n[c],t1);
		clerp(&newt1.c[2],t.c[a],t.c[c],t1);
		newt1.t[2].s = lerp(t.t[a].s,t.t[c].s,t1);
		newt1.t[2].t = lerp(t.t[a].t,t.t[c].t,t1);
		drawclippedtriangle(newt1);
		newt2.m = t.m;
		newt2.v[0].x = t.v[a].x+p1.x*t1;
		newt2.v[0].y = t.v[a].y+p1.y*t1;
		newt2.v[0].z = screen.znear;
		newt2.v[0].w = 1.f;
		vlerp(&newt2.n[0],t.n[a],t.n[c],t1);
		clerp(&newt2.c[0],t.c[a],t.c[c],t1);
		newt2.t[0].s = lerp(t.t[a].s,t.t[c].s,t1);
		newt2.t[0].t = lerp(t.t[a].t,t.t[c].t,t1);
		newt2.v[1] = t.v[b];
		newt2.n[1] = t.n[b];
		newt2.c[1] = t.c[b];
		newt2.t[1] = t.t[b];
		newt2.v[2].x = t.v[b].x+p2.x*t2;
		newt2.v[2].y = t.v[b].y+p2.y*t2;
		newt2.v[2].z = screen.znear;
		newt2.v[2].w = 1.f;
		vlerp(&newt2.n[2],t.n[b],t.n[c],t2);
		clerp(&newt2.c[2],t.c[b],t.c[c],t2);
		newt2.t[2].s = lerp(t.t[b].s,t.t[c].s,t2);
		newt2.t[2].t = lerp(t.t[b].t,t.t[c].t,t2);
		drawclippedtriangle(newt2);
	}
	else drawclippedtriangle(t);
}

pixel_t clearcolor = {16,16,16,255};

int showfps = 1;
int showhelp = 1;
int drawwire = 0;
int autorot = 0;

void rendermodel( void )
{
	const int this_is_invariant_i_swear = screen.width*screen.height;
	#pragma omp parallel for
	for ( int i=0; i<this_is_invariant_i_swear; i++ )
		screen.color[i] = clearcolor;
	#pragma omp parallel for
	for ( int i=0; i<this_is_invariant_i_swear; i++ )
		screen.depth[i] = screen.zfar;
	mat_t fulltransform;
	vect_t basev[3], transv[3];
	vect_t basen[3], transn[3];
	fulltransform = rotation;
	mmul(&fulltransform,fulltransform,translation);
	mmul(&fulltransform,fulltransform,projection);
	for ( int i=0; i<cube.ntri; i++ )
	{
		for ( int j=0; j<3; j++ )
		{
			basev[j].x = cube.vertices[cube.triangles[i].v[j]].x;
			basev[j].y = cube.vertices[cube.triangles[i].v[j]].y;
			basev[j].z = cube.vertices[cube.triangles[i].v[j]].z;
			basev[j].w = 1.f;
			vmat(&transv[j],fulltransform,basev[j]);
			basen[j].x = cube.normals[cube.triangles[i].n[j]].x;
			basen[j].y = cube.normals[cube.triangles[i].n[j]].y;
			basen[j].z = cube.normals[cube.triangles[i].n[j]].z;
			basen[j].w = 1.f;
			vmat(&transn[j],rotation,basen[j]);
		}
		tri_t face =
		{
			{transv[0],transv[1],transv[2]},
			{transn[0],transn[1],transn[2]},
			{cube.txcoords[cube.triangles[i].t[0]],
			cube.txcoords[cube.triangles[i].t[1]],
			cube.txcoords[cube.triangles[i].t[2]]},
			{cube.colors[cube.triangles[i].c[0]],
			cube.colors[cube.triangles[i].c[1]],
			cube.colors[cube.triangles[i].c[2]]},
			&cube.materials[cube.triangles[i].mat]
		};
		if ( drawwire < 2 ) drawtriangle(face,0);
		if ( drawwire > 0 ) drawtriangle(face,1);
	}
}

int main( void )
{
	SDL_Init(SDL_INIT_VIDEO|SDL_INIT_EVENTS);
	SDL_Window *win = SDL_CreateWindow("RasterTest",SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,screen.width,screen.height,
		SDL_WINDOW_SHOWN);
	SDL_Surface *scr = SDL_GetWindowSurface(win);
	screen.color = malloc(sizeof(pixel_t)*screen.width*screen.height);
	screen.depth = malloc(sizeof(float)*screen.width*screen.height);
	SDL_Surface *fb = SDL_CreateRGBSurfaceFrom(screen.color,screen.width,
		screen.height,32,sizeof(pixel_t)*screen.width,0xFF0000,0xFF00,
		0xFF,0xFF000000);
	TTF_Init();
	TTF_Font *fon = TTF_OpenFont("/usr/share/fonts/misc/unifont.bdf",16);
	SDL_Event e;
	SDL_Color fpscol = {160,0,0,255};
	int active = 1;
	float t = 0.f;
	float frame = 0.f, fps = NAN, fpsmin = INFINITY, fpsmax = -INFINITY,
		fpsavg = 0.f;
	long tick, tock;
	float fh = tanf(90.f/360.f*PI)*screen.znear,
		fw = fh*(screen.width/(float)screen.height);
	frustum(&projection,-fw,fw,-fh,fh,screen.znear,screen.zfar);
	char fpst[16];
	float rx = 0, ry = 0;
	while ( active )
	{
		tick = ticker();
		while ( SDL_PollEvent(&e) )
		{	if ( e.type == SDL_QUIT ) active = 0;
			if ( e.type == SDL_KEYDOWN )
			{
				if ( e.key.keysym.sym == SDLK_ESCAPE )
					active = 0;
				else if ( e.key.keysym.sym == SDLK_a )
					cubepos.x -= 0.1f;
				else if ( e.key.keysym.sym == SDLK_d )
					cubepos.x += 0.1f;
				else if ( e.key.keysym.sym == SDLK_q )
					cubepos.y -= 0.1f;
				else if ( e.key.keysym.sym == SDLK_e )
					cubepos.y += 0.1f;
				else if ( e.key.keysym.sym == SDLK_w )
					cubepos.z += 0.1f;
				else if ( e.key.keysym.sym == SDLK_s )
					cubepos.z -= 0.1f;
				else if ( e.key.keysym.sym == SDLK_LEFT )
					ry += 0.025f;
				else if ( e.key.keysym.sym == SDLK_RIGHT )
					ry -= 0.025f;
				else if ( e.key.keysym.sym == SDLK_UP )
					rx -= 0.025f;
				else if ( e.key.keysym.sym == SDLK_DOWN )
					rx += 0.025f;
				else if ( e.key.keysym.sym == SDLK_z )
					autorot = !autorot;
				else if ( e.key.keysym.sym == SDLK_x )
					drawwire = (drawwire<2)?(drawwire+1):0;
				else if ( e.key.keysym.sym == SDLK_f )
					showfps = !showfps;
				else if ( e.key.keysym.sym == SDLK_F1 )
					showhelp = !showhelp;
				else if ( e.key.keysym.sym == SDLK_c )
					culling = (culling<2)?(culling+1):0;
			}
		}
		mat_t rotx, roty;
		if ( autorot )
		{
			rx += frame;
			ry += frame;
		}
		rotate(&rotx,50.f*rx,ROT_X);
		rotate(&roty,80.f*ry,ROT_Y);
		mmul(&rotation,roty,rotx);
		translate(&translation,cubepos);
		rendermodel();
		SDL_BlitSurface(fb,0,scr,0);
		SDL_Surface *txt;
		if ( showfps )
		{
			snprintf(fpst,15,"%.2f FPS",fps);
			txt = TTF_RenderText_Blended(fon,fpst,fpscol);
			SDL_BlitSurface(txt,0,scr,0);
			SDL_FreeSurface(txt);
		}
		if ( showhelp )
		{
			SDL_Rect dr =
				{0,screen.height-TTF_FontHeight(fon)*9,-1,-1};
			for ( int i=0; i<9; i++ )
			{
				txt = TTF_RenderUTF8_Blended(fon,helptext[i],
					fpscol);
				SDL_BlitSurface(txt,0,scr,&dr);
				SDL_FreeSurface(txt);
				dr.y += TTF_FontHeight(fon);
			}
		}
		SDL_UpdateWindowSurface(win);
		tock = ticker();
		frame = (float)(tock-tick)/NANOS_SEC;
		fps = 1.f/frame;
		fpsavg = avg_fps(fps);
		if ( fps > fpsmax ) fpsmax = fps;
		if ( fps < fpsmin ) fpsmin = fps;
		t += frame;
	}
	free(screen.color);
	free(screen.depth);
	TTF_CloseFont(fon);
	TTF_Quit();
	SDL_FreeSurface(fb);
	SDL_DestroyWindow(win);
	SDL_Quit();
	printf("FPS: %.2f min, %.2f max, %.2f avg\n",fpsmin,fpsmax,fpsavg);
	return 0;
}
