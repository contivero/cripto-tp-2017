#include <dirent.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "util.h"

#define BMP_MAGIC_NUMBER     0x424D
#define BMP_HEADER_SIZE      14
#define DIB_HEADER_SIZE      40
#define PALETTE_SIZE         1024
#define PIXEL_ARRAY_OFFSET   (BMP_HEADER_SIZE + DIB_HEADER_SIZE + PALETTE_SIZE)
#define UNUSED2_OFFSET       8
#define WIDTH_OFFSET         18
#define HEIGHT_OFFSET        22
#define BITS_PER_PIXEL       8
#define PRIME                257
#define DEFAULT_SEED         691
#define RIGHTMOST_BIT_ON(x)  ((x) |= 0x01)
#define RIGHTMOST_BIT_OFF(x) ((x) &= 0xFE)
#define DIR_MAX              (PATH_MAX - NAME_MAX)

typedef struct {
    uint8_t  id[2];   /* magic number to identify the BMP format */
    uint32_t size;    /* size of the BMP file in bytes */
    uint16_t unused1; /* key (seed) */
    uint16_t unused2; /* shadow number */
    uint32_t offset;  /* starting address of the pixel array (bitmap data) */
} BMPheader;

/* 40 bytes BITMAPINFOHEADER */
typedef struct {
    uint32_t size;           /* the size of this header (40 bytes) */
    uint32_t width;          /* the bitmap width in pixels */
    int32_t  height;         /* the bitmap height in pixels; can be negative */
    uint16_t nplanes;        /* number of color planes used; Must set to 1 */
    uint16_t depth;          /* bpp number. Usually: 1, 4, 8, 16, 24 or 32 */
    uint32_t compression;    /* compression method used */
    uint32_t pixelarraysize; /* size of the raw bitmap (pixel) data */
    uint32_t hres;           /* horizontal resolution (pixel per meter) */
    uint32_t vres;           /* vertical resolution (pixel per meter) */
    uint32_t ncolors;        /* colors in the palette. 0 means 2^n */
    uint32_t nimpcolors;     /* important colors used, usually ignored */
} DIBheader;

typedef struct {
    BMPheader bmpheader;             /* 14 bytes BMP starting header */
    DIBheader dibheader;             /* 40 bytes DIB header */
    uint8_t   palette[PALETTE_SIZE]; /* color palette; mandatory for depth <= 8 */
    uint8_t   *imgpixels;            /* array of bytes representing each pixel */
} Bitmap;

/* prototypes */
static void     setseed(int64_t s);
static int      nextbyte(void);
static int      countfiles(const char *dirname);
static void     usage(void);
static uint32_t get32bitsfromheader(FILE *fp, int offset);
static uint32_t bmpfilewidth(FILE *fp);
static uint32_t bmpfileheight(FILE *fp);
static uint32_t bmpimagesize(Bitmap *bp);
static void     initpalette(uint8_t palette[]);
static Bitmap   *newbitmap(uint32_t width, int32_t height, uint16_t seed);
static void     freebitmap(Bitmap *bp);
static Bitmap   *newbitmaphelper(uint32_t width, int32_t height, uint16_t seed, uint16_t shadnum, uint32_t pixelarraysize);
static void     changeheaderendianness(BMPheader *h);
static void     changedibendianness(DIBheader *h);
static void     readbmpheader(Bitmap *bp, FILE *fp);
static void     writebmpheader(Bitmap *bp, FILE *fp);
static void     readdibheader(Bitmap *bp, FILE *fp);
static void     writedibheader(Bitmap *bp, FILE *fp);
static Bitmap   *bmpfromfile(const char *filename);
static bool     isvalidbmpsize(FILE *fp, uint16_t k, uint32_t secretsize);
static int      kdivisiblesize(FILE *fp, uint16_t k);
static void     bmptofile(Bitmap *bp, const char *filename);
static void     findclosestpair(int value, uint32_t *v1, int32_t *v2);
static Bitmap   *newshadow(uint32_t width, int32_t height, uint16_t seed, uint16_t shadownumber);
static Bitmap   **formshadows(Bitmap *bp, uint16_t seed, uint16_t k, uint16_t n);
static void     findcoefficients(int **mat, uint16_t k);
static void     revealsecret(Bitmap **shadows, uint16_t k, uint32_t width, int32_t height, const char *filename);
static void     hideshadow(Bitmap *bp, Bitmap *shadow);
static Bitmap   *retrieveshadow(Bitmap *bp, uint32_t width, int32_t height, uint16_t k);
static bool     isbmp(FILE *fp);
static bool     isvalidshadow(FILE *fp, uint16_t k, uint32_t size);
static bool     isvalidbmp(FILE *fp, uint16_t k, uint32_t ignored);
static void     getvalidfilenames(char **filenames, char *dir, uint16_t k, uint16_t n, bool (*isvalid)(FILE *, uint16_t, uint32_t), uint32_t);
static void     getbmpfilenames(char **filenames, char *dir, uint16_t k, uint16_t n, uint32_t size);
static void     getshadowfilenames(char **filenames, char *dir, uint16_t k, uint32_t size);
static void     distributeimage(uint16_t k, uint16_t n, uint16_t seed, char *imgpath, char *dir);
static void     recoverimage(uint16_t k, uint32_t width, int32_t height, char *filename, char *dir);
static uint32_t calculatepixelarraysize(uint32_t width, int32_t height);
static void     decreasecoeff(uint8_t *coeff);
static uint8_t  *randomtable(uint32_t tablesize, uint16_t seed);
static void     xorbmpwithrandomtable(Bitmap *bmp, uint16_t seed);

/* globals */
static char *argv0;                /* program name for usage() */
static int64_t rseed;              /* seed to use for the random table */
static const int modinv[PRIME] = { /* modular multiplicative inverses */
    0, 1, 129, 86, 193, 103, 43, 147, 225, 200, 180, 187, 150, 178, 202, 120,
    241, 121, 100, 230, 90, 49, 222, 190, 75, 72, 89, 238, 101, 195, 60, 199,
    249, 148, 189, 235, 50, 132, 115, 145, 45, 163, 153, 6, 111, 40, 95, 175,
    166, 21, 36, 126, 173, 97, 119, 243, 179, 248, 226, 61, 30, 59, 228, 102,
    253, 87, 74, 234, 223, 149, 246, 181, 25, 169, 66, 24, 186, 247, 201, 244,
    151, 165, 210, 96, 205, 127, 3, 65, 184, 26, 20, 209, 176, 152, 216, 46, 83,
    53, 139, 135, 18, 28, 63, 5, 215, 164, 177, 245, 188, 224, 250, 44, 218,
    116, 124, 38, 113, 134, 159, 54, 15, 17, 158, 140, 114, 220, 51, 85, 255, 2,
    172, 206, 37, 143, 117, 99, 240, 242, 203, 98, 123, 144, 219, 133, 141, 39,
    213, 7, 33, 69, 12, 80, 93, 42, 252, 194, 229, 239, 122, 118, 204, 174, 211,
    41, 105, 81, 48, 237, 231, 73, 192, 254, 130, 52, 161, 47, 92, 106, 13, 56,
    10, 71, 233, 191, 88, 232, 76, 11, 108, 34, 23, 183, 170, 4, 155, 29, 198,
    227, 196, 31, 9, 78, 14, 138, 160, 84, 131, 221, 236, 91, 82, 162, 217, 146,
    251, 104, 94, 212, 112, 142, 125, 207, 22, 68, 109, 8, 58, 197, 62, 156, 19,
    168, 185, 182, 67, 35, 208, 167, 27, 157, 136, 16, 137, 55, 79, 107, 70, 77,
    57, 32, 110, 214, 154, 64, 171, 128, 256
};

int
countfiles(const char *dirname){
    struct dirent *d;
    int filecount = 0;
    DIR *dp = xopendir(dirname);

    while((d = readdir(dp)))
        if(d->d_type == DT_REG) /* If the entry is a regular file */
            filecount++;
    xclosedir(dp);

    return filecount;
}

/* Algorithm based on Java's Random, which itself was defined by D. H. Lehmer
 * and described by Knuth in The Art of Computer Programming, Volume 3:
 * Seminumerical Algorithms, section 3.2.1.
 * See: https://docs.oracle.com/javase/8/docs/api/java/util/Random.html#setSeed-long-
 */
void
setseed(int64_t s){
    rseed = (s ^ 25214903917) & 281474976710655;
}

int
nextbyte(void) {
    rseed = (rseed * 25214903917 + 11) & 281474976710655;
    int n = rseed >> (48 - 31);

    return (int64_t)256 * (int64_t)n >> 31;
}

void
usage(void){
    die("usage: %s -(d|r) --secret image -k number -w width -h height -s seed"
            "[-n number] [--dir directory]\n", argv0);
}

/* Calculates needed pixelarraysize, accounting for padding.
 * See: https://en.wikipedia.org/wiki/BMP_file_format#Pixel_storage */
inline uint32_t
calculatepixelarraysize(uint32_t width, int32_t height){
    return ((BITS_PER_PIXEL * width + 31)/32) * 4 * height;
}

uint32_t
get32bitsfromheader(FILE *fp, int offset){
    uint32_t value;
    long pos = ftell(fp);

    xfseek(fp, offset, SEEK_SET);
    xfread(&value, sizeof(value), 1, fp);
    xfseek(fp, pos, SEEK_SET);

    return value;
}

uint32_t
bmpfilewidth(FILE *fp){
    return get32bitsfromheader(fp, WIDTH_OFFSET);
}

uint32_t
bmpfileheight(FILE *fp){
    return get32bitsfromheader(fp, HEIGHT_OFFSET);
}

/* initialize palette with default 8-bit greyscale values */
void
initpalette(uint8_t palette[]){
    for(uint32_t i = 0; i < 256; i++){
        uint32_t j = i * 4;
        palette[j++] = i;
        palette[j++] = i;
        palette[j++] = i;
        palette[j] = 0;
    }
}

/* If no seed is needed, just pass 0 */
Bitmap *
newbitmap(uint32_t width, int32_t height, uint16_t seed){
    uint32_t pixelarraysize = calculatepixelarraysize(width, height);
    Bitmap *bmp = newbitmaphelper(width, height, seed, 0, pixelarraysize);

    return bmp;
}

/* Helper function to build a BMP, used by newbitmap() and newshadow() */
Bitmap*
newbitmaphelper(uint32_t width, int32_t height, uint16_t seed, uint16_t shadnum, uint32_t pixelarraysize){
    Bitmap *bmp = xmalloc(sizeof(*bmp));

    bmp->imgpixels = xmalloc(pixelarraysize);
    initpalette(bmp->palette);

    bmp->bmpheader.id[0]   = 'B';
    bmp->bmpheader.id[1]   = 'M';
    bmp->bmpheader.size    = PIXEL_ARRAY_OFFSET + pixelarraysize;
    bmp->bmpheader.unused1 = seed;
    bmp->bmpheader.unused2 = shadnum;
    bmp->bmpheader.offset  = PIXEL_ARRAY_OFFSET;

    bmp->dibheader.size           = DIB_HEADER_SIZE;
    bmp->dibheader.width          = width;
    bmp->dibheader.height         = height;
    bmp->dibheader.nplanes        = 1;
    bmp->dibheader.depth          = BITS_PER_PIXEL;
    bmp->dibheader.compression    = 0;
    bmp->dibheader.pixelarraysize = pixelarraysize;
    bmp->dibheader.hres           = 0;
    bmp->dibheader.vres           = 0;
    bmp->dibheader.ncolors        = 0;
    bmp->dibheader.nimpcolors     = 0;

    return bmp;

}

void
freebitmap(Bitmap *bp){
    free(bp->imgpixels);
    free(bp);
}

void
changeheaderendianness(BMPheader *h){
    uint32swap(&h->size);
    uint16swap(&h->unused1);
    uint16swap(&h->unused2);
    uint32swap(&h->offset);
}

void
changedibendianness(DIBheader *h){
    uint32swap(&h->size);
    uint32swap(&h->width);
    int32swap(&h->height);
    uint16swap(&h->nplanes);
    uint16swap(&h->depth);
    uint32swap(&h->compression);
    uint32swap(&h->pixelarraysize);
    uint32swap(&h->hres);
    uint32swap(&h->vres);
    uint32swap(&h->ncolors);
    uint32swap(&h->nimpcolors);
}

void
readbmpheader(Bitmap *bp, FILE *fp){
    BMPheader *h = &bp->bmpheader;

    xfread(h->id, sizeof(h->id), 1, fp);
    xfread(&h->size, sizeof(h->size), 1, fp);
    xfread(&h->unused1, sizeof(h->unused1), 1, fp);
    xfread(&h->unused2, sizeof(h->unused2), 1, fp);
    xfread(&h->offset, sizeof(h->offset), 1, fp);

    if(isbigendian())
        changeheaderendianness(&bp->bmpheader);
}

void
writebmpheader(Bitmap *bp, FILE *fp){
    BMPheader h = bp->bmpheader;

    if(isbigendian())
        changeheaderendianness(&h);

    xfwrite(h.id, sizeof(h.id), 1, fp);
    xfwrite(&(h.size), sizeof(h.size), 1, fp);
    xfwrite(&(h.unused1), sizeof(h.unused1), 1, fp);
    xfwrite(&(h.unused2), sizeof(h.unused2), 1, fp);
    xfwrite(&(h.offset), sizeof(h.offset), 1, fp);
}

void
readdibheader(Bitmap *bp, FILE *fp){
    DIBheader *h = &bp->dibheader;

    xfread(&h->size, sizeof(h->size), 1, fp);
    xfread(&h->width, sizeof(h->width), 1, fp);
    xfread(&h->height, sizeof(h->height), 1, fp);
    xfread(&h->nplanes, sizeof(h->nplanes), 1, fp);
    xfread(&h->depth, sizeof(h->depth), 1, fp);
    xfread(&h->compression, sizeof(h->compression), 1, fp);
    xfread(&h->pixelarraysize, sizeof(h->pixelarraysize), 1, fp);
    xfread(&h->hres, sizeof(h->hres), 1, fp);
    xfread(&h->vres, sizeof(h->vres), 1, fp);
    xfread(&h->ncolors, sizeof(h->ncolors), 1, fp);
    xfread(&h->nimpcolors, sizeof(h->nimpcolors), 1, fp);

    if(isbigendian())
        changedibendianness(&bp->dibheader);
}

void
writedibheader(Bitmap *bp, FILE *fp){
    DIBheader h = bp->dibheader;

    if(isbigendian())
        changedibendianness(&h);

    xfwrite(&(h.size), sizeof(h.size), 1, fp);
    xfwrite(&(h.width), sizeof(h.width), 1, fp);
    xfwrite(&(h.height), sizeof(h.height), 1, fp);
    xfwrite(&(h.nplanes), sizeof(h.nplanes), 1, fp);
    xfwrite(&(h.depth), sizeof(h.depth), 1, fp);
    xfwrite(&(h.compression), sizeof(h.compression), 1, fp);
    xfwrite(&(h.pixelarraysize), sizeof(h.pixelarraysize), 1, fp);
    xfwrite(&(h.hres), sizeof(h.hres), 1, fp);
    xfwrite(&(h.vres), sizeof(h.vres), 1, fp);
    xfwrite(&(h.ncolors), sizeof(h.ncolors), 1, fp);
    xfwrite(&(h.nimpcolors), sizeof(h.nimpcolors), 1, fp);
}

Bitmap *
bmpfromfile(const char *filename){
    FILE *fp = xfopen(filename, "r");
    Bitmap *bp = xmalloc(sizeof(*bp));

    readbmpheader(bp, fp);
    readdibheader(bp, fp);
    xfread(bp->palette, sizeof(bp->palette), 1, fp);

    /* read pixel data */
    uint32_t imagesize = bmpimagesize(bp);
    bp->imgpixels = xmalloc(imagesize);
    xfread(bp->imgpixels, sizeof(bp->imgpixels[0]), imagesize, fp);
    xfclose(fp);

    return bp;
}

bool
isvalidbmpsize(FILE *fp, uint16_t k, uint32_t secretsize){
    uint32_t shadowsize = (secretsize * 8)/k;
    uint32_t imgsize    = bmpfilewidth(fp) * bmpfileheight(fp);

    return imgsize >= shadowsize;
}

int
kdivisiblesize(FILE *fp, uint16_t k){
    int pixels = bmpfilewidth(fp) * bmpfileheight(fp);
    int aux    = pixels / k;

    return pixels == aux * k;
}

uint32_t
bmpimagesize(Bitmap *bp){
    uint32_t sz = bp->bmpheader.size;

    if(sz == 0)
        return bp->dibheader.pixelarraysize;
    else
        return bp->bmpheader.size - bp->bmpheader.offset;
}

void
bmptofile(Bitmap *bp, const char *filename){
    FILE *fp = xfopen(filename, "w");

    writebmpheader(bp, fp);
    writedibheader(bp, fp);
    xfwrite(bp->palette, PALETTE_SIZE, 1, fp);
    xfwrite(bp->imgpixels, bmpimagesize(bp), 1, fp);
    xfclose(fp);
}

/* find closest pair of values that when multiplied, give x.
 * Used to make the shadows as 'squared' as possible */
void
findclosestpair(int x, uint32_t *width, int32_t *height){
    int y = floor(sqrt(x));

    for(; y > 2; y--)
        if(x % y == 0){
            *width  = y;
            *height = x / y;
            break;
        }
}

Bitmap *
newshadow(uint32_t width, int32_t height, uint16_t seed, uint16_t shadownumber){
    Bitmap *bmp = newbitmaphelper(width, height, seed, shadownumber, width * height);
    bmp->bmpheader.unused2 = shadownumber;

    return bmp;
}

Bitmap **
formshadows(Bitmap *bp, uint16_t seed, uint16_t k, uint16_t n){
    unsigned int i, j;
    uint8_t *coeff;
    uint32_t width;
    int32_t height;
    uint32_t pixelarraysize = bmpimagesize(bp);
    Bitmap **shadows = xmalloc(sizeof(*shadows) * n);
    uint16_t *pixels = xmalloc(sizeof(*pixels) * n);

    findclosestpair(pixelarraysize/k, &width, &height);

    /* allocate shadows */
    for(i = 0; i < n; i++)
        shadows[i] = newshadow(width, height, seed, i+1);

    /* generate shadow image pixels */
    for(j = 0; j*k < pixelarraysize; j++){
        coeff = &bp->imgpixels[j*k];

        /* Paper's 4th step, mixed with the 3rd one */
        step4:
        for(i = 0; i < n; i++){
            /* uses coeff[0] to coeff[k-1] (where k-1 is the degree of the
             * polynomial) to evaluate the corresponding section polynomial and
             * generate a pixel for the i-th shadow image */
            uint64_t value = 0;

            for(int r = 0; r < k; r++)
                value += coeff[r] * pow(i+1, r);

            pixels[i] = value % PRIME;
        }

        for(i = 0; i < n; i++){
            if(pixels[i] == 256){
                decreasecoeff(coeff);
                goto step4;
            }
        }

        for(i = 0; i < n; i++)
            shadows[i]->imgpixels[j] = pixels[i];
    }
    free(pixels);

    return shadows;
}

/* decrease the first non-zero coefficient by one */
void
decreasecoeff(uint8_t *coeff){
    int i = 0;

    /* We can assume some value is 0, as it is proved in the paper */
    while(coeff[i] == 0)
        i++;

    coeff[i]--;
}

void
findcoefficients(int **mat, uint16_t k){
    int i, j, t, a, temp;

    /* take matrix to echelon form */
    for(j = 0; j < k-1; j++){
        for(i = k-1; i > j; i--){
            a = mat[i][j] * modinv[mat[i-1][j]];
            for(t = j; t < k+1; t++){
                temp = mat[i][t] - ((mat[i-1][t] * a) % PRIME);
                mat[i][t] = mod(temp, PRIME);
            }
        }
    }

    /* take matrix to reduced row echelon form */
    for(i = k-1; i > 0; i--){
        mat[i][k] = (mat[i][k] * modinv[mat[i][i]]) % PRIME;
        mat[i][i] = (mat[i][i] * modinv[mat[i][i]]) % PRIME;
        for(t = i-1; t >= 0; t--){
            temp = mat[t][k] - ((mat[i][k] * mat[t][i]) % PRIME);
            mat[t][k] = mod(temp, PRIME);
            mat[t][i] = 0;
        }
    }
}

void
revealsecret(Bitmap **shadows, uint16_t k, uint32_t width, int32_t height, const char *filename){
    int i, j, t, value;
    int pixels = (*shadows)->dibheader.pixelarraysize;
    Bitmap *sp;
    Bitmap *bmp = newbitmap(width, height, (*shadows)->bmpheader.unused1);

    int **mat = xmalloc(sizeof(*mat) * k);
    for(i = 0; i < k; i++)
        mat[i] = xmalloc(sizeof(**mat) * (k+1));

    for(i = 0; i < pixels; i++){
        for(j = 0; j < k; j++){
            sp = shadows[j];
            value = sp->bmpheader.unused2;
            mat[j][0] = 1;
            for(t = 1; t < k; t++){
                mat[j][t] = value;
                value *= sp->bmpheader.unused2;
            }
            mat[j][k] = sp->imgpixels[i];
        }
        findcoefficients(mat, k);
        for(j = i * k; j < (i+1) * k; j++){
            bmp->imgpixels[j] = mat[j % k][k];
        }
    }

    xorbmpwithrandomtable(bmp, (*shadows)->bmpheader.unused1);
    bmptofile(bmp, filename);

    for(i = 0; i < k; i++)
        free(mat[i]);
    free(mat);
    freebitmap(bmp);
}

void
hideshadow(Bitmap *bp, Bitmap *shadow){
    unsigned int i, j;
    char shadowfilename[20] = {0};
    uint32_t pixels = bmpimagesize(shadow);

    bp->bmpheader.unused1 = shadow->bmpheader.unused1;
    bp->bmpheader.unused2 = shadow->bmpheader.unused2;
    xsnprintf(shadowfilename, 20, "shadow%d.bmp", shadow->bmpheader.unused2);

    for(i = 0; i < pixels; i++){
        uint8_t byte = shadow->imgpixels[i];
        for(j = i*8; j < 8*(i+1); j++){
            if(byte & 0x80) /* 1000 0000 */
                RIGHTMOST_BIT_ON(bp->imgpixels[j]);
            else
                RIGHTMOST_BIT_OFF(bp->imgpixels[j]);
            byte <<= 1;
        }
    }
    bmptofile(bp, shadowfilename);
}

/* width and height parameters needed because the image hiding the shadow could
 * be bigger than necessary */
Bitmap *
retrieveshadow(Bitmap *bp, uint32_t width, int32_t height, uint16_t k){
    uint16_t key          = bp->bmpheader.unused1;
    uint16_t shadownumber = bp->bmpheader.unused2;

    findclosestpair(calculatepixelarraysize(width, height)/k, &width, &height);
    Bitmap *shadow = newshadow(width, height, key, shadownumber);
    uint32_t shadowpixels = shadow->dibheader.pixelarraysize;

    for(uint32_t i = 0; i < shadowpixels; i++){
        uint8_t byte = 0;
        uint8_t mask = 0x80; /* 1000 0000 */
        for(uint32_t j = i*8; j < 8*(i+1); j++){
            if(bp->imgpixels[j] & 0x01)
                byte |= mask;
            mask >>= 1;
        }
        shadow->imgpixels[i] = byte;
    }

    return shadow;
}

bool
isbmp(FILE *fp){
    char magicnumber[2];
    long pos = ftell(fp);

    xfseek(fp, 0, SEEK_SET);
    xfread(magicnumber, sizeof(magicnumber), 1, fp);
    xfseek(fp, pos, SEEK_SET);

    return magicnumber[0] == 'B' && magicnumber[1] == 'M';
}

bool
isvalidshadow(FILE *fp, uint16_t k, uint32_t secretsize){
    uint16_t shadownumber;
    long pos = ftell(fp);

    xfseek(fp, UNUSED2_OFFSET, SEEK_SET);
    xfread(&shadownumber, sizeof(shadownumber), 1, fp);
    xfseek(fp, pos, SEEK_SET);

    return shadownumber && isbmp(fp) && isvalidbmpsize(fp, k, secretsize);
}

/* the last parameter is ignored, and is only present so that the function
 * prototype can be used with getvalidfilenames() */
bool
isvalidbmp(FILE *fp, uint16_t k, uint32_t ignoredparameter){
    return isbmp(fp) && kdivisiblesize(fp, k);
}

void
getvalidfilenames(char **filenames, char *dir, uint16_t k, uint16_t n, bool (*isvalid)(FILE *, uint16_t, uint32_t), uint32_t size){
    struct dirent *d;
    FILE *fp;
    DIR *dp = xopendir(dir);
    int i = 0;
    char filepath[PATH_MAX] = {0};

    while((d = readdir(dp)) && i < n){
        if(d->d_type == DT_REG){
            int len = xsnprintf(filepath, PATH_MAX, "%.*s/%.*s", DIR_MAX, dir, NAME_MAX, d->d_name);
            fp = xfopen(filepath, "r");
            if(isvalid(fp, k, size)){
                filenames[i] = xmalloc(len + 1);
                strncpy(filenames[i], filepath, len);
                filenames[i][len] = '\0'; /* NULL terminate string */
                i++;
            }
            xfclose(fp);
        }
    }
    xclosedir(dp);

    if(i < n)
        die("not enough valid bmps for a (%d,%d) threshold scheme in dir %s\n", k, n, dir);
}

void
getbmpfilenames(char **filenames, char *dir, uint16_t k, uint16_t n, uint32_t size){
    getvalidfilenames(filenames, dir, k, n, isvalidbmp, size);
}

void
getshadowfilenames(char **filenames, char *dir, uint16_t k, uint32_t size){
    getvalidfilenames(filenames, dir, k, k, isvalidshadow, size);
}

void
distributeimage(uint16_t k, uint16_t n, uint16_t seed, char *imgpath, char *dir){
    char **filepaths = xmalloc(sizeof(*filepaths) * n);
    Bitmap *bmp, **shadows;
    int i;

    bmp = bmpfromfile(imgpath);
    getbmpfilenames(filepaths, dir, k, n, bmpimagesize(bmp));
    xorbmpwithrandomtable(bmp, seed);
    shadows = formshadows(bmp, seed, k, n);
    freebitmap(bmp);

    for(i = 0; i < n; i++){
        bmp = bmpfromfile(filepaths[i]);
        hideshadow(bmp, shadows[i]);
        freebitmap(bmp);
    }

    for(i = 0; i < n; i++){
        free(filepaths[i]);
        freebitmap(shadows[i]);
    }
    free(filepaths);
    free(shadows);
}

void
recoverimage(uint16_t k, uint32_t width, int32_t height, char *filename, char *dir){
    char **filepaths = xmalloc(sizeof(*filepaths) * k);
    Bitmap **shadows = xmalloc(sizeof(*shadows) * k);
    unsigned int i;

    getshadowfilenames(filepaths, dir, k, width * height);
    for(i = 0; i < k; i++){
        Bitmap *bp = bmpfromfile(filepaths[i]);
        shadows[i] = retrieveshadow(bp, width, height, k);
        freebitmap(bp);
    }

    revealsecret(shadows, k, width, height, filename);

    for(i = 0; i < k; i++){
        free(filepaths[i]);
        freebitmap(shadows[i]);
    }
    free(filepaths);
    free(shadows);
}

uint8_t *
randomtable(uint32_t tablesize, uint16_t seed) {
    uint8_t *table = xmalloc(tablesize * sizeof(*table));
    setseed(seed);

    if(table)
        for(unsigned int i = 0; i < tablesize; i++)
            table[i] = nextbyte();

    return table;
}

void
xorbmpwithrandomtable(Bitmap *bmp, uint16_t seed){
    uint32_t imgsize = bmpimagesize(bmp);
    uint8_t *table   = randomtable(imgsize, seed);

    for(unsigned int i = 0; i < imgsize; i++)
        bmp->imgpixels[i] ^= table[i];

    free(table);
}

int
main(int argc, char *argv[]){
    char *filename = NULL;
    char *dir = "./";
    int dflag = 0;
    int rflag = 0;
    int kflag = 0;
    int wflag = 0;
    int hflag = 0;
    int nflag = 0;
    int secretflag = 0;
    uint16_t seed = DEFAULT_SEED;
    uint16_t k = 0, n = 0;
    uint32_t width = 0;
    int32_t height = 0;
    int i;

    argv0 = argv[0]; /* save program name for usage() */

    for(i = 1; i < argc; i++){
        if(strcmp(argv[i], "-d") == 0){
            dflag = 1;
        } else if(strcmp(argv[i], "-r") == 0) {
            rflag = 1;
        } else if(strcmp(argv[i], "--secret") == 0){
            secretflag = 1;
            if(i + 1 < argc){
                filename = argv[++i];
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "-k") == 0){
            kflag = 1;
            if(i + 1 < argc){
                k = atoi(argv[++i]);
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "-w") == 0){
            wflag = 1;
            if(i + 1 < argc){
                width = atoi(argv[++i]);
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "-h") == 0){
            hflag = 1;
            if(i + 1 < argc){
                height = atoi(argv[++i]);
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "-s") == 0) {
            if(i + 1 < argc){
                seed = atoi(argv[++i]);
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "-n") == 0){
            nflag = 1;
            if(i + 1 < argc){
                n = atoi(argv[++i]);
            } else {
                usage();
            }
        } else if(strcmp(argv[i], "--dir") == 0){
            if(i + 1 < argc){
                dir = argv[++i];
            } else{
                usage();
            }
        } else {
            die("invalid %s parameter \n", argv[i]);
        }
    }

    if(!(dflag || rflag) || !secretflag || !kflag)
        usage();
    if((rflag && !(wflag && hflag)) || !width || !height)
        die("specify a positive width and height with -w -h for the revealed image\n");

    if(!nflag)
        n = countfiles(dir);

    if(k > n || k < 2 || n < 2)
        die("k and n must be: 2 <= k <= n\n");
    if(dflag && rflag)
        die("can't use -d and -r flags simultaneously\n");

    if(dflag)
        distributeimage(k, n, seed, filename, dir);
    else if(rflag)
        recoverimage(k, width, height, filename, dir);

    return EXIT_SUCCESS;
}
