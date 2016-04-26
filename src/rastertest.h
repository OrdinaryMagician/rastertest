/*
	rastertest.h : Just learning.
	This code is a mess, be warned.
	(C)2016 Marisa Kirisame, UnSX Team.
	Released under the GNU GPLv3 (or later).
*/
typedef struct
{
	float x,y,z,w;
} vect_t;
typedef struct
{
	float s,t;
} coord_t;
typedef struct
{
	float r,g,b,a;
} color_t;
typedef struct
{
	float c[4][4];
} mat_t;
typedef struct
{
	float x,y,d;
	int s;
} px_t;
typedef struct
{
	unsigned char b,g,r,a;
} __attribute__((packed)) pixel_t;
typedef struct
{
	int v[3], n[3], t[3], c[3], mat;
} face_t;
typedef struct
{
	pixel_t *color;
	float *depth;
	int *stencil;
	int width, height;
	float znear, zfar;
} buffer_t;
#define TXCLAMP  0
#define TXWRAP   1
#define TXMIRROR 2
#define TXBORDER 3
typedef struct
{
	pixel_t *data;
	int width, height, depth;
	int bu, bv;
	color_t border;
} texture_t;
/*
   material reference:
    diffuse: rgba color
    normal: rgb tangent, a parallax
    lighting: r roughness, g metalness, b glowness, a unused
    cubemap: rgba color
*/
#define BLNONE  0
#define BLALPHA 1
#define BLADD   2
#define BLMUL   3
#define BLSUB   4
#define BLDIV   5
#define BLOVER  6
#define BLMOD   7
typedef struct
{
	texture_t *diffuse, *normal, *lighting, *cubemap;
	color_t diffmul, lightmul;
	int blend;
} material_t;
typedef struct
{
	vect_t v[3], n[3];
	coord_t t[3];
	color_t c[3];
	material_t *m;
} tri_t;
typedef struct
{
	vect_t *vertices;
	vect_t *normals;
	coord_t *txcoords;
	color_t *colors;
	material_t *materials;
	face_t *triangles;
	int nvert, nnorm, ncoord, ncolor, nmat, ntri;
} model_t;
void vadd( vect_t *o, vect_t a, vect_t b )
{
	o->x = a.x+b.x;
	o->y = a.y+b.y;
	o->z = a.z+b.z;
	o->w = a.w+b.w;
}
void vmul( vect_t *o, vect_t a, vect_t b )
{
	o->x = a.x*b.x;
	o->y = a.y*b.y;
	o->z = a.z*b.z;
	o->w = a.w*b.w;
}
void vscale( vect_t *o, vect_t a, float b )
{
	o->x = a.x*b;
	o->y = a.y*b;
	o->z = a.z*b;
	o->w = a.w*b;
}
void vmat( vect_t *o, mat_t a, vect_t b )
{
	o->x = a.c[0][0]*b.x+a.c[1][0]*b.y+a.c[2][0]*b.z+a.c[3][0]*b.w;
	o->y = a.c[0][1]*b.x+a.c[1][1]*b.y+a.c[2][1]*b.z+a.c[3][1]*b.w;
	o->z = a.c[0][2]*b.x+a.c[1][2]*b.y+a.c[2][2]*b.z+a.c[3][2]*b.w;
	o->w = a.c[0][3]*b.x+a.c[1][3]*b.y+a.c[2][3]*b.z+a.c[3][3]*b.w;
}
void mident( mat_t *o )
{
	int i,j;
	for ( i=0; i<4; i++ ) for ( j=0; j<4; j++ )
		o->c[i][j] = (i==j)?1.f:0.f;
}
void mmul( mat_t *o, mat_t a, mat_t b )
{
	int i,j;
	for ( i=0; i<4; i++ ) for ( j=0; j<4; j++ )
		o->c[i][j] = a.c[i][0]*b.c[0][j]+a.c[i][1]*b.c[1][j]
			+a.c[i][2]*b.c[2][j]+a.c[i][3]*b.c[3][j];
}
void mscale( mat_t *o, mat_t a, float b )
{
	int i,j;
	for ( i=0; i<4; i++ ) for ( j=0; j<4; j++ )
		o->c[i][j] = a.c[i][j]*b;
}
void cadd( color_t *o, color_t a, color_t b )
{
	o->r = a.r+b.r;
	o->g = a.g+b.g;
	o->b = a.b+b.b;
	o->a = a.a+b.a;
}
void cmul( color_t *o, color_t a, color_t b )
{
	o->r = a.r*b.r;
	o->g = a.g*b.g;
	o->b = a.b*b.b;
	o->a = a.a*b.a;
}
void cscale( color_t *o, color_t a, float b )
{
	o->r = a.r*b;
	o->g = a.g*b;
	o->b = a.b*b;
	o->a = a.a*b;
}
#define PI 3.14159265359f
void frustum( mat_t *o, float left, float right, float bottom, float top,
	float near, float far )
{
	o->c[0][0] = (2.f*near)/(right-left);
	o->c[0][1] = 0.f;
	o->c[0][2] = (right+left)/(right-left);
	o->c[0][3] = 0.f;
	o->c[1][0] = 0.f;
	o->c[1][1] = (2.f*near)/(top-bottom);
	o->c[1][2] = (top+bottom)/(top-bottom);
	o->c[1][3] = 0.f;
	o->c[2][0] = 0.f;
	o->c[2][1] = 0.f;
	o->c[2][2] = (far+near)/(far-near);
	o->c[2][3] = -(2.f*far*near)/(far-near);
	o->c[3][0] = 0.f;
	o->c[3][1] = 0.f;
	o->c[3][2] = 1.f;
	o->c[3][3] = 0.f;
}
#define ROT_X 0
#define ROT_Y 1
#define ROT_Z 2
void rotate( mat_t *o, float angle, int axis )
{
	float theta = (angle/180.f)*PI;
	vect_t n = {0.f,0.f,0.f,0.f};
	if ( axis == 0 ) n.x = 1.f;
	else if ( axis == 1 ) n.y = 1.f;
	else n.z = 1.f;
	float s = sinf(theta);
	float c = cosf(theta);
	float oc = 1.f-c;
	o->c[0][0] = oc*n.x*n.x+c;
	o->c[0][1] = oc*n.x*n.y-n.z*s;
	o->c[0][2] = oc*n.z*n.x+n.y*s;
	o->c[0][3] = 0.f;
	o->c[1][0] = oc*n.x*n.y+n.z*s;
	o->c[1][1] = oc*n.y*n.y+c;
	o->c[1][2] = oc*n.y*n.z-n.x*s;
	o->c[1][3] = 0.f;
	o->c[2][0] = oc*n.z*n.x-n.y*s;
	o->c[2][1] = oc*n.y*n.z+n.x*s;
	o->c[2][2] = oc*n.z*n.z+c;
	o->c[2][3] = 0.f;
	o->c[3][0] = 0.f;
	o->c[3][1] = 0.f;
	o->c[3][2] = 0.f;
	o->c[3][3] = 1.f;
}
void translate( mat_t *o, vect_t offset )
{
	o->c[0][0] = 1.f;
	o->c[0][1] = 0.f;
	o->c[0][2] = 0.f;
	o->c[0][3] = 0.f;
	o->c[1][0] = 0.f;
	o->c[1][1] = 1.f;
	o->c[1][2] = 0.f;
	o->c[1][3] = 0.f;
	o->c[2][0] = 0.f;
	o->c[2][1] = 0.f;
	o->c[2][2] = 1.f;
	o->c[2][3] = 0.f;
	o->c[3][0] = offset.x;
	o->c[3][1] = offset.y;
	o->c[3][2] = offset.z;
	o->c[3][3] = 1.f;
}
#define saturate(x) (((x)>0)?((x)>1)?1:(x):0)
#define swap(t,x,y) {t tmp=x;x=y;y=tmp;}
#define abs(x) (((x)>0)?(x):-(x))
#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))
#define clamp(a,b,c) (((a)>(b))?((a)>(c))?(c):(a):(b))
#define lerp(a,b,f) ((a)*(1.f-(f))+(b)*(f))
void vlerp( vect_t *o, vect_t a, vect_t b, float f )
{
	o->x = lerp(a.x,b.x,f);
	o->y = lerp(a.y,b.y,f);
	o->z = lerp(a.z,b.z,f);
	o->w = lerp(a.w,b.w,f);
}
void clerp( color_t *o, color_t a, color_t b, float f )
{
	o->r = lerp(a.r,b.r,f);
	o->g = lerp(a.g,b.g,f);
	o->b = lerp(a.b,b.b,f);
	o->a = lerp(a.a,b.a,f);
}
float vsize( vect_t v )
{
	return sqrt(powf(v.x,2.f)+powf(v.y,2.f)+powf(v.z,2.f));
}
void cross( vect_t *a, vect_t b, vect_t c )
{
	a->x = b.y*c.z-b.z*c.y;
	a->y = b.z*c.x-b.x*c.z;
	a->z = b.x*c.y-b.y*c.x;
}
float dot( vect_t a, vect_t b )
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
void normalize( vect_t *a )
{
	float scale = sqrtf(powf(a->x,2.f)+powf(a->y,2.f)+powf(a->z,2.f));
	a->x /= scale;
	a->y /= scale;
	a->z /= scale;
}
