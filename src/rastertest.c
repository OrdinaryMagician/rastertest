/*
	rastertest.c : Just learning.
	This code is a mess, be warned.
	(C)2016 Marisa Kirisame, UnSX Team.
	Released under the GNU GPLv3 (or later).
*/
#define _DEFAULT_SOURCE

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rastertest.h"

vect_t cubevect[8] =
{
	{-1.f,-1.f, 1.f,1.f},{ 1.f,-1.f, 1.f,1.f},
	{-1.f, 1.f, 1.f,1.f},{ 1.f, 1.f, 1.f,1.f},
	{-1.f,-1.f,-1.f,1.f},{ 1.f,-1.f,-1.f,1.f},
	{-1.f, 1.f,-1.f,1.f},{ 1.f, 1.f,-1.f,1.f}
};
vect_t cubenorm[6] =
{
	{1.f,0.f,0.f,1.f},{-1.f,0.f,0.f,1.f},
	{0.f,1.f,0.f,1.f},{0.f,-1.f,0.f,1.f},
	{0.f,0.f,1.f,1.f},{0.f,0.f,-1.f,1.f}
};
coord_t cubecoord[4] =
{
	{0.f,0.f},{1.f,0.f},
	{0.f,1.f},{1.f,1.f}
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
	{{0,6,4},{0,0,0},{0,3,2},{0,6,4},0},
	{{0,2,6},{0,0,0},{0,1,3},{0,2,6},0},
	{{0,3,2},{0,0,0},{0,3,2},{0,3,2},0},
	{{0,1,3},{0,0,0},{0,1,3},{0,1,3},0},
	{{2,7,6},{0,0,0},{0,3,2},{2,7,6},0},
	{{2,3,7},{0,0,0},{0,1,3},{2,3,7},0},
	{{4,6,7},{0,0,0},{0,1,3},{4,6,7},0},
	{{4,7,5},{0,0,0},{0,3,2},{4,7,5},0},
	{{0,4,5},{0,0,0},{0,1,3},{0,4,5},0},
	{{0,5,1},{0,0,0},{0,3,2},{0,5,1},0},
	{{1,5,7},{0,0,0},{0,1,3},{1,5,7},0},
	{{1,7,3},{0,0,0},{0,3,2},{1,7,3},0}
};
model_t cube =
{
	cubevect,cubenorm,cubecoord,cubecolor,0,cubefaces,8,6,4,8,0,12
};
vect_t cubepos = {0.f,0.f,-10.f,1.f};
vect_t cuberot = {0.f,0.f,0.f,0.f};
mat_t projection, translation, rotation;
buffer_t screen = {0,0,0,640,480,0.1f,100.0f};

pixel_t cubepixels[256] =
{
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,  0,  0,255},{  0,  0,255,255},{  0,  0,255,255},
	{  0,  0,255,255},{  0,  0,255,255},{  0,  0,255,255},{  0,  0,255,255},
	{  0,  0,255,255},{  0,  0,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,255,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,  0,  0,255},{  0,  0,  0,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{  0,  0,  0,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,  0,  0,255},{  0,  0,  0,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{  0,  0,  0,255},{255,255,255,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{  0,  0,  0,255},{  0,  0,  0,255},
	{255,255,255,255},{  0,  0,  0,255},{  0,  0,  0,255},{  0,  0,  0,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255},
	{255,255,255,255},{255,255,255,255},{255,255,255,255},{255,255,255,255}
};
texture_t cubetex =
{
	cubepixels,16,16,0,TXCLAMP,TXCLAMP,{0.f,0.f,0.f,1.f}
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
	int cx = floorf(c.s*t.width), cy = floorf(c.t*t.height);
	cx = clamp(cx,0,t.width-1);
	cy = clamp(cy,0,t.height-1);
	o->r = t.data[cx+cy*t.width].r/255.f;
	o->g = t.data[cx+cy*t.width].g/255.f;
	o->b = t.data[cx+cy*t.width].b/255.f;
	o->a = t.data[cx+cy*t.width].a/255.f;
}

void filltriangle( px_t p[3], coord_t tc[3], color_t c[3] )
{
	px_t minb, maxb;
	maxb.x = clamp(max(p[0].x,max(p[1].x,p[2].x)),0,screen.width);
	maxb.y = clamp(max(p[0].y,max(p[1].y,p[2].y)),0,screen.height);
	minb.x = clamp(min(p[0].x,min(p[1].x,p[2].x)),0,screen.width);
	minb.y = clamp(min(p[0].y,min(p[1].y,p[2].y)),0,screen.height);
	if ( (maxb.x == minb.x) || (maxb.y == minb.y) ) return;
	vect_t d1 = {p[1].x-p[0].x,p[1].y-p[0].y,0.f,0.f},
		d2 = {p[2].x-p[0].x,p[2].y-p[0].y,0.f,0.f};
	color_t c1 = {c[1].r-c[0].r,c[1].g-c[0].g,c[1].b-c[0].b,c[1].a-c[0].a},
		c2 = {c[2].r-c[0].r,c[2].g-c[0].g,c[2].b-c[0].b,c[2].a-c[0].a};
	coord_t t1 = {tc[1].s-tc[0].s,tc[1].t-tc[0].t},
		t2 = {tc[2].s-tc[0].s,tc[2].t-tc[0].t};
	float dp1 = p[1].d-p[0].d, dp2 = p[2].d-p[0].d;
	#pragma omp parallel for
	for ( int y=minb.y; y<maxb.y; y++ )
	{
		for ( int x=minb.x; x<maxb.x; x++ )
		{
			vect_t q = {x-p[0].x,y-p[0].y,0.f,0.f};
			float u = d1.x*d2.y-d1.y*d2.x,
				s = (q.x*d2.y-q.y*d2.x)/u,
				t = (d1.x*q.y-d1.y*q.x)/u;
			if ( (s>=0.f) && (t>=0.f) && (s+t<=1.f) )
			{
				int coord = x+y*screen.width;
				float dep = screen.depth[coord];
				float pd = p[0].d+s*dp1+t*dp2;
				if ( (pd < screen.znear) || (pd > screen.zfar)
					|| (pd > dep) ) continue;
				color_t pc =
				{
					c[0].r+s*c1.r+t*c2.r,
					c[0].g+s*c1.g+t*c2.g,
					c[0].b+s*c1.b+t*c2.b,
					c[0].a+s*c1.a+t*c2.a
				};
				coord_t tx =
				{
					tc[0].s+s*t1.s+t*t2.s,
					tc[0].t+s*t1.t+t*t2.t
				};
				sampletexture(&pc,tx,cubetex);
				px_t px = {x,y,pd,0};
				putpixel(px,pc);
			}
		}
	}
}

void drawtriangle( tri_t t, unsigned char wireframe )
{
	px_t pts[3] =
	{
		{floorf((t.v[0].x*0.5f+0.5f)*screen.width),
			floorf((t.v[0].y*0.5f+0.5f)*screen.height),-t.v[0].z,
			0.f},
		{floorf((t.v[1].x*0.5f+0.5f)*screen.width),
			floorf((t.v[1].y*0.5f+0.5f)*screen.height),-t.v[1].z,
			0.f},
		{floorf((t.v[2].x*0.5f+0.5f)*screen.width),
			floorf((t.v[2].y*0.5f+0.5f)*screen.height),-t.v[2].z,
			0.f},
	};
	if ( wireframe )
	{
		color_t wire = {0.f,0.6f,0.f,1.f};
		drawline(pts[0].x,pts[1].x,pts[0].y,pts[1].y,wire);
		drawline(pts[0].x,pts[2].x,pts[0].y,pts[2].y,wire);
		drawline(pts[1].x,pts[2].x,pts[1].y,pts[2].y,wire);
		return;
	}
	vect_t facet, ab, ac;
	ab.x = t.v[1].x-t.v[0].x;
	ab.y = t.v[1].y-t.v[0].y;
	ab.z = t.v[1].z-t.v[0].z;
	ac.x = t.v[2].x-t.v[0].x;
	ac.y = t.v[2].y-t.v[0].y;
	ac.z = t.v[2].z-t.v[0].z;
	cross(&facet,ab,ac);
	if ( facet.z < 0.f ) return;
	color_t cols[3] = {t.c[0],t.c[1],t.c[2]};
	coord_t coords[3] = {t.t[0],t.t[1],t.t[2]};
	filltriangle(pts,coords,cols);
}

pixel_t clearcolor = {16,16,16,255};
float cleardepth = INFINITY;

void rendermodel( void )
{
	int i;
	for ( i=0; i<screen.width*screen.height; i++ )
	{
		screen.color[i] = clearcolor;
		screen.depth[i] = cleardepth;
	}
	mat_t fulltransform;
	vect_t basev[3], transv[3];
	fulltransform = rotation;
	mmul(&fulltransform,fulltransform,translation);
	mmul(&fulltransform,fulltransform,projection);
	for ( int i=0; i<cube.ntri; i++ )
	{
		#pragma omp parallel for
		for ( int j=0; j<3; j++ )
		{
			basev[j].x = cube.vertices[cube.triangles[i].v[j]].x;
			basev[j].y = cube.vertices[cube.triangles[i].v[j]].y;
			basev[j].z = cube.vertices[cube.triangles[i].v[j]].z;
			basev[j].w = 1.f;
			vmat(&transv[j],fulltransform,basev[j]);
			transv[j].x /= transv[j].w;
			transv[j].y /= transv[j].w;
			transv[j].z /= transv[j].w;
		}
		tri_t face =
		{
			{transv[0],transv[1],transv[2]},
			{{0,0,0,0},{0,0,0,0},{0,0,0,0}},
			{cube.txcoords[cube.triangles[i].t[0]],
			cube.txcoords[cube.triangles[i].t[1]],
			cube.txcoords[cube.triangles[i].t[2]]},
			{cube.colors[cube.triangles[i].c[0]],
			cube.colors[cube.triangles[i].c[1]],
			cube.colors[cube.triangles[i].c[2]]},
			0
		};
		drawtriangle(face,0);
		drawtriangle(face,1);
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
			}
		}
		mat_t rotx, roty;
		rotate(&rotx,50.f*t,ROT_X);
		rotate(&roty,80.f*t,ROT_Y);
		mmul(&rotation,rotx,roty);
		translate(&translation,cubepos);
		rendermodel();
		SDL_BlitSurface(fb,0,scr,0);
		snprintf(fpst,15,"%.2f FPS",fps);
		SDL_Surface *txt = TTF_RenderText_Blended(fon,fpst,fpscol);
		SDL_BlitSurface(txt,0,scr,0);
		SDL_FreeSurface(txt);
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
