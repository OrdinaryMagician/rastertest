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
	{{1,3,0},{0,0,0},{0,1,2},{1,3,0},0},
	{{7,5,4},{1,1,1},{3,2,0},{7,5,4},0},
	{{4,1,0},{2,2,2},{0,1,2},{4,1,0},0},
	{{5,2,1},{3,3,3},{3,2,1},{5,2,1},0},
	{{2,7,3},{4,4,4},{3,2,0},{2,7,3},0},
	{{0,7,4},{5,5,5},{3,2,1},{0,7,4},0},
	{{1,2,3},{0,0,0},{0,3,1},{1,2,3},0},
	{{7,6,5},{1,1,1},{3,1,2},{7,6,5},0},
	{{4,5,1},{2,2,2},{0,3,1},{4,5,1},0},
	{{5,6,2},{3,3,3},{3,0,2},{5,6,2},0},
	{{2,6,7},{4,4,4},{3,1,2},{2,6,7},0},
	{{0,3,7},{5,5,5},{3,0,2},{0,3,7},0}
};
model_t cube =
{
	cubevect,cubenorm,cubecoord,cubecolor,0,cubefaces,8,6,4,8,0,12
};
vect_t cubepos = {0.f,0.f,-5.f,1.f};
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
	screen.color[coord].r = saturate(c.r)*255.f;
	screen.color[coord].g = saturate(c.g)*255.f;
	screen.color[coord].b = saturate(c.b)*255.f;
	screen.color[coord].a = saturate(c.a)*255.f;
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

void fragmentprogram( color_t *out, frag_t data )
{
	color_t res;
	sampletexture(&res,data.txcoord,cubetex);
	cmul(out,res,data.color);
	vect_t light = {0,0,16,0};
	vect_t eye = {0,0,1,0};
	vsub(&light,light,cubepos);
	vsub(&eye,eye,cubepos);
	vsub(&light,light,data.position);
	vsub(&eye,eye,data.position);
	normalize(&light);
	normalize(&eye);
	vect_t ref;
	reflect(&ref,light,data.normal);
	vscale(&ref,ref,-1.0);
	normalize(&ref);
	color_t amb = {0.03,0.03,0.03,1.0};
	color_t diff = {1.0,1.0,1.0,1.0};
	cscale(&diff,diff,max(dot(data.normal,light),0.0));
	color_t spec = {0.3,0.3,0.3,1.0};
	cscale(&spec,spec,powf(max(dot(ref,eye),0.0),4.0));
	color_t lit;
	cadd(&lit,amb,diff);
	cmul(out,*out,lit);
	cadd(out,*out,spec);
}

void filltriangle( px_t p[3], coord_t tc[3], color_t c[3], vect_t vp[3],
	vect_t vn[3] )
{
	px_t minb, maxb;
	maxb.x = clamp(max(p[0].x,max(p[1].x,p[2].x)),0,screen.width);
	maxb.y = clamp(max(p[0].y,max(p[1].y,p[2].y)),0,screen.height);
	minb.x = clamp(min(p[0].x,min(p[1].x,p[2].x)),0,screen.width);
	minb.y = clamp(min(p[0].y,min(p[1].y,p[2].y)),0,screen.height);
	if ( (maxb.x == minb.x) || (maxb.y == minb.y) ) return;
	for ( int i=0; i<3; i++ )
	{
		tc[i].s /= p[i].d;
		tc[i].t /= p[i].d;
		c[i].r /= p[i].d;
		c[i].g /= p[i].d;
		c[i].b /= p[i].d;
		c[i].a /= p[i].d;
		vp[i].x /= p[i].d;
		vp[i].y /= p[i].d;
		vp[i].z /= p[i].d;
		vp[i].w /= p[i].d;
		vn[i].x /= p[i].d;
		vn[i].y /= p[i].d;
		vn[i].z /= p[i].d;
		vn[i].w /= p[i].d;
	}
	for ( int i=0; i<3; i++ ) p[i].d = 1.f/p[i].d;
	float area = afn(p[0],p[1],p[2]);
	for ( int y=minb.y; y<maxb.y; y++ )
	{
		for ( int x=minb.x; x<maxb.x; x++ )
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
			coord_t tx =
			{
				pd*(q.x*tc[0].s+q.y*tc[1].s+q.z*tc[2].s),
				pd*(q.x*tc[0].t+q.y*tc[1].t+q.z*tc[2].t)
			};
			color_t pc =
			{
				pd*(q.x*c[0].r+q.y*c[1].r+q.z*c[2].r),
				pd*(q.x*c[0].g+q.y*c[1].g+q.z*c[2].g),
				pd*(q.x*c[0].b+q.y*c[1].b+q.z*c[2].b),
				pd*(q.x*c[0].a+q.y*c[1].a+q.z*c[2].a)
			};
			vect_t fp =
			{
				pd*(q.x*vp[0].x+q.y*vp[1].x+q.z*vp[2].x),
				pd*(q.x*vp[0].y+q.y*vp[1].y+q.z*vp[2].y),
				pd*(q.x*vp[0].z+q.y*vp[1].z+q.z*vp[2].z),
				pd*(q.x*vp[0].w+q.y*vp[1].w+q.z*vp[2].w)
			};
			vect_t fn =
			{
				pd*(q.x*vn[0].x+q.y*vn[1].x+q.z*vn[2].x),
				pd*(q.x*vn[0].y+q.y*vn[1].y+q.z*vn[2].y),
				pd*(q.x*vn[0].z+q.y*vn[1].z+q.z*vn[2].z),
				pd*(q.x*vn[0].w+q.y*vn[1].w+q.z*vn[2].w)
			};
			px.x = x;
			px.y = y;
			px.d = pd;
			frag_t frag = {px,fp,fn,tx,pc,0};
			fragmentprogram(&pc,frag);
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
	vect_t pos[3] = {t.v[0],t.v[1],t.v[2]};
	vect_t norm[3] = {t.n[0],t.n[1],t.n[2]};
	filltriangle(pts,coords,cols,pos,norm);
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
int autorot = 1;

void rendermodel( void )
{
	/* and here we have buffer clearing, where 50% of the fps go die */
	for ( int i=0; i<screen.width*screen.height; i++ )
		screen.color[i] = clearcolor;
	for ( int i=0; i<screen.width*screen.height; i++ )
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
			&cube.materials[cube.triangles[i].mat],
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
	float fh = tanf(70.f/360.f*PI)*screen.znear,
		fw = fh*(screen.width/(float)screen.height);
	frustum(&projection,-fw,fw,-fh,fh,screen.znear,screen.zfar);
	char fpst[16];
	float rx = 0, rz = 0;
	while ( active )
	{
		tick = ticker();
		while ( SDL_PollEvent(&e) )
		{	if ( e.type == SDL_QUIT ) active = 0;
			if ( e.type == SDL_KEYDOWN )
			{
				if ( e.key.keysym.sym == SDLK_a )
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
					rx += 0.025f;
				else if ( e.key.keysym.sym == SDLK_RIGHT )
					rx -= 0.025f;
				else if ( e.key.keysym.sym == SDLK_UP )
					rz -= 0.025f;
				else if ( e.key.keysym.sym == SDLK_DOWN )
					rz += 0.025f;
				/* no repeated input beyond this point due to a
				   SDL2 bug */
				if ( e.key.repeat ) break;
				if ( e.key.keysym.sym == SDLK_ESCAPE )
					active = 0;
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
		mat_t rotx, rotz;
		if ( autorot )
		{
			rx -= frame;
			rz -= frame;
		}
		rotate(&rotx,90.f*rx,ROT_Y);
		rotate(&rotz,90.f*rz,ROT_Z);
		mmul(&rotation,rotz,rotx);
		translate(&translation,cubepos);
		rendermodel();
		SDL_BlitSurface(fb,0,scr,0);
		SDL_Surface *txt;
		if ( showfps )
		{
			snprintf(fpst,15,"%.2f FPS",fpsavg);
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
