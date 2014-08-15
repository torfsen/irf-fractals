/*
 * Fractals from iterated root finding methods.
 *
 * Compile and run this C program to generate fractals from a complex function
 * and a root finding method. See the comments below for more information.
 *
 * Note: Must be compiled using the -lm switch. Also requires a compiler
 * supporting C99, use -std=c99 for gcc.
 *
 * Copyright (c) 2010, Florian Brucker (mail@florianbrucker.de)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Allocate memory and check for failure
#define MALLOC(p, n) if ((p = malloc((n) * sizeof(*p))) == NULL) { \
    fputs("Could not allocate memory.\n", stderr); \
    exit(EXIT_FAILURE); \
}
#define REALLOC(p, n) if ((p = realloc(p, (n) * sizeof(*p))) == NULL) { \
    fputs("Could not allocate memory.\n", stderr); \
    exit(EXIT_FAILURE); \
}


/********************************************************************
 * METHOD WRAPPER
 *******************************************************************/

// Type of a function for which a fractal should be generated
typedef complex double (*fun)(complex double);

/*
 * Each root finding method has the same signature: A double function that takes
 * three function pointers (the function and its first two derivatives), and a
 * starting value. The return value is the next approximation.
 */
#define SIG fun f, fun df, fun ddf, complex double x

/*
 * Type of a root finding method, taking pointers to a function and its first
 * two derivatives, and a starting point. Returns the next point in the
 * iteration.
 */
typedef complex double (*method)(SIG);

/*
 * This wrapper includes all the necessary step counting and accuracy
 * thresholding, so that the method implementations themselves are
 * as simple as possible.
 *
 * Arguments:
 *   m: Root finding method
 *   f: Function to generate fractal for
 *   df: First derivative
 *   ddf: Second derivative
 *   x: Starting point for root finding iteration
 *   e: Accuracy threshold
 *   ms: Maximum step size
 */
int iterative_root_finding(method m, fun f, fun df, fun ddf, complex double *x,
                           double e, int ms) {
    int i;
    for (i = 0; i <= ms; i++) {
        if (cabs((*f)(*x)) < e) {
            // Accuracy threshold reached
            break;
        }
        *x = (*m)(f, df, ddf, *x);
    }
    return i;
}


/********************************************************************
 * METHODS
 *******************************************************************/

/*
 * This is a collection of various root finding methods. See the website given
 * for each one for a more detailed description. Note that not all of these
 * were designed for use with complex numbers, but our goal is to create nice
 * fractals, not mathematical rigor.
 */

// http://mathworld.wolfram.com/NewtonsMethod.html
complex double newton(SIG) {
    return x - (*f)(x) / (*df)(x);
}

// http://mathworld.wolfram.com/HouseholdersMethod.html 
complex double householder(SIG) {
    complex double fx = (*f)(x);
    complex double dfx = (*df)(x);
    complex double ddfx = (*ddf)(x);
    return x - fx / dfx * (1 + f(x) * ddfx / (2 * dfx * dfx));
}

// http://mathworld.wolfram.com/HalleysIrrationalFormula.html
complex double halley_irrational(SIG) {
    complex double fx = (*f)(x);
    complex double dfx = (*df)(x);
    complex double ddfx = (*ddf)(x);
    return x + (-dfx + csqrt(dfx * dfx - 2 * fx * ddfx)) / ddf(x);
}

// http://mathworld.wolfram.com/HalleysMethod.html
complex double halley(SIG) {
    complex double fx = (*f)(x);
    complex double dfx = (*df)(x);
    complex double ddfx = (*ddf)(x);
    return x - 2 * fx * dfx / (2 * dfx * dfx - fx * ddfx);
}

// http://mathworld.wolfram.com/SchroedersMethod.html
complex double schroeder(SIG) {
    complex double fx = (*f)(x);
    complex double dfx = (*df)(x);
    complex double ddfx = (*ddf)(x);
    return x - fx * dfx / (dfx * dfx - fx * ddfx);
}

// http://en.wikipedia.org/wiki/Steffensen%27s_method
complex double steffensen(SIG) {
    complex double fx = (*f)(x);
    return x - fx * fx / ((*f)(x + fx) - fx);
}


/********************************************************************
 * RASTERING
 *******************************************************************/

/*
 * This method performs the actual rastering. That is, it loops through all
 * pixels, performs the root finding, manages the list of found roots, etc.
 *
 * Arguments:
 *   m: IRF method to use
 *   f: Function to search roots of
 *   df: First derivative of f
 *   ddf: Second derivative of f
 *   e: Accuracy threshold
 *   ms: Maximal step count
 *   x: real coordinates
 *   nx: number of real coordinates
 *   y: imaginary coordinates
 *   ny: number of imaginary coordinates
 *   r: 2D array of size nx times ny, the indices of the roots
 *   maxr: The highest index of a root
 *   s: 2D array of size nx times ny, the number of steps used
 *   maxs: The highest number of steps used in a converging approximation
 *   roots: List of roots, allocated by the function
 */

void raster(method m, fun f, fun df, fun ddf, double e, int ms, double *x,
            int nx, double *y, int ny, int **r, int *maxr, int **s, int *maxs,
            complex double **roots) {
    int nroots = 10;  // Initial number of root-slots
    int iroots = 0;  // Index of first unused root-slot
    MALLOC(*roots, nroots)
    *maxs = 0;
    int last_disp = -1;
    for (int i = 0; i < nx; i++) {
        int p = (10 * i) / nx;
        if (p > last_disp) {
            printf("\t%d0%%\n", p);
            last_disp = p;
        }
        for (int j = 0; j < ny; j++) {
            complex double p = x[j] + y[i] * I;
            s[i][j] = iterative_root_finding(m, f, df, ddf, &p, e, ms);
            if (s[i][j] < ms) {
                if (s[i][j] > *maxs) {
                    *maxs = s[i][j];
                }
                // Check which root we found
                int best_index = -1;
                double best_dist = 1e-2;
                for (int k = 0; k < iroots; k++) {
                    double dist = cabs((*roots)[k] - p);
                    if (dist < best_dist) {
                        best_index = k;
                        best_dist = dist;
                    }
                }
                if (best_index > -1) {
                    r[i][j] = best_index + 1;  // 0 means no convergence
                } else {
                    // This seems to be a new root
                    if (iroots >= nroots) {
                        // No free slots available, get some more
                        nroots *= 2;
                        REALLOC(*roots, nroots)
                    }
                    (*roots)[iroots] = p;
                    r[i][j] = iroots;
                    printf("\tRoot #%d: %f + %fi\n", iroots + 1, creal(p), cimag(p));
                    iroots++;
                }
            } else {
                // No convergence
                r[i][j] = 0;
            }
        }
    }
    *maxr = iroots - 1;
}


/*
 * C version of MATLAB's linspace method. Returns an array of 
 * <steps> doubles, equally spaced between <a> and <b>.
 */
double* linspace(double a, double b, int steps) {
    double *x;
    MALLOC(x, steps)
    if (steps > 1) {
        for (int i = 0; i < steps; i++) {
            x[i] = a + i * (b - a) / (steps - 1);
        }
    } else {
        x[0] = a;
    }
    return x;
}


/********************************************************************
 * IMAGE OUTPUT
 *******************************************************************/

/**
 * Writes image data in TGA format to a file.
 *
 * Uses information from Paul Bourke's website at
 * http://local.wasp.uwa.edu.au/~pbourke/dataformats/tga/
 *
 * Arguments:
 *   red    2D image data array (Red part)
 *   green  2D image data array (Green part)
 *   blue   2D image data array (Blue part)
 *   width  Width of the image
 *   height Height of the image
 *   file   Filename to save data to
 *
 * Returns 1 if successful, 0 otherwise
 */
int write_tga(unsigned char **red, unsigned char **green, unsigned char **blue, int width, int height, char *file) {
    int i, j;
    FILE *tga_file = fopen(file, "wb");
    if (tga_file == NULL) return 0;

    // Write TGA header
    putc(0, tga_file);
    putc(0, tga_file);
    putc(2, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc(0, tga_file);
    putc((width & 0x00ff), tga_file);
    putc((width & 0xff00) / 256, tga_file);
    putc((height & 0x00ff), tga_file);
    putc((height & 0xff00) / 256, tga_file);
    putc(24, tga_file);
    putc(0, tga_file);

    // Write image data
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            putc(blue[i][j], tga_file);
            putc(green[i][j], tga_file);
            putc(red[i][j], tga_file);
        }
    }

    if (fclose(tga_file) != 0) return 0;
    return 1;
}


/*
 * Converts a color from HSV to RGB space. Input arguments are the HSV
 * coordinates, which then get overwritten with the corresponding RGB
 * coordinates.
 */
void hsv2rgb(double *x, double *y, double *z) {
    *x *= 6;
    int i = floor(*x);
    double f = *x - i;
    if (!(i & 1)) {
        f = 1 - f;
    }
    double m = *z * (1 - *y);
    double n = *z * (1 - *y * f);
    switch (i) {
        case 6:
        case 0: *x = *z; *y = n; *z = m; break;
        case 1: *x = n; *y = *z; *z = m; break;
        case 2: *x = m; *y = *z; *z = n; break;
        case 3: *x = m; *y = n; *z = *z; break;
        case 4: *x = n; *y = m; *z = *z; break;
        case 5: *x = *z; *y = m; *z = n; break;
    }
}


/********************************************************************
 * TARGET FUNCTIONS AND SETTINGS
 *******************************************************************/

/*
 * These are some nice setups for image generation. Make sure that only one
 * of them is uncommented.
 */

/*
 * Some shortcuts for defining whole setups in a compact way
 */
#define F(X) complex double f(complex double x) { return X; }
#define DF(X) complex double df(complex double x) { return X; }
#define DDF(X) complex double ddf(complex double x) { return X; }
#define G(A, B, C, X1, X2, Y1, Y2) F(A) DF(B) DDF(C) double x1 = X1; \
    double x2 = X2; double y1 = Y1; double y2 = Y2;


// 2*x^3 - 2*x + 2
//G(2 * cpow(x, 3) - 2 * x + 2, 6 * x * x - 2, 12 * x, 1.53, 2.05, -0.26, 0.26) // Newton
G(2 * cpow(x, 3) - 2 * x + 2, 6 * x * x - 2, 12 * x, -1.47, 0.23, -0.5, 1.2) // Householder
//G(2 * cpow(x, 3) - 2 * x + 2, 6 * x * x - 2, 12 * x, 0.2, 0.4, -0.1, 0.1) // Halley
//G(2 * cpow(x, 3) - 2 * x + 2, 6 * x * x - 2, 12 * x, 1, 1.7, -0.35, 0.35) // Schroeder
//G(2 * cpow(x, 3) - 2 * x + 2, 6 * x * x - 2, 12 * x, 0.29, 0.57, -0.75, -0.47) // Steffensen

// x^3 - 3^x
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), 4.4, 5.4, 7.2, 8.2) // Newton (Some black areas remain, they need more than 1000000 steps for convergence)
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), -6, 6, -6, 6) // Steffensen doesn't converge outside of this
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), 1.48, 3.58, -9.61, -7.51) // Schroeder
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), 5.14, 6.14, 5.82, 6.82) // Halley
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), 3.9, 5.8, 6.55, 8.45) // Householder
//G(cpow(x, 3) - cpow(3, x), 3 * cpow(x, 2) - clog(3) * cexp(clog(3) * x), 6 * x - cpow(clog(3), 2) * cexp(clog(3) * x), -0.23, 0.09, -0.16, 0.16) // Steffensen (black areas don't converge even fafter 1000000 steps)

// cos(x) (make sure SINGLE_COLOR is defined)
//G(ccos(x), -csin(x), -ccos(x), -0.435, 0.435, -0.435, 0.435) // Newton
//G(ccos(x), -csin(x), -ccos(x), -0.65, 0.65, -0.65, 0.65) // Householder
//G(ccos(x), -csin(x), -ccos(x), -0.4, 0.4, -0.4 + 2.925, 0.4 + 2.925) // Halley
//G(ccos(x), -csin(x), -ccos(x), -3.25, 3.25, -3.25, 3.25) // Schroeder
//G(ccos(x), -csin(x), -ccos(x), -8 + 3.14/2, 8 + 3.14/2, -8, 8) // Steffensen


/*
 * Settings
 */

// Maximum number of steps per pixel
#define MAX_STEPS 1000

// Side-lengh of the image in pixels
#define IMAGE_SIZE 1000

//Available methods: newton householder halley_irrational halley
//                   schroeder steffensen
#define METHOD schroeder

// Use gray-scale?
//#define SINGLE_COLOR


/********************************************************************
 * MAIN ROUTINE
 *******************************************************************/

/*
 * This function is a weighting function for the saturation of the
 * colors. If one just scales the saturation linearly from 0 to 1
 * according to the number of steps needed for convergence, then
 * the image will be very bright, due to a low number of points needing
 * a high number of steps. Thus, the scaling is done adaptively,
 * depending on the average number of steps needed. The function below
 * evaluates a polygon P at the point x, where P has the following
 * properties:
 *
 * P(0) = 0
 * P'(0) = 0
 * P(1) = 1
 * P'(1) = 1
 * P(x_0) = y_0
 */
double blend(double x, double x0, double y0) {
    double x02 = x0 * x0;
    double x03 = x02 * x0;
    double x04 = x03 * x0;

    double c = (y0 + 3 * x04 - 4 * x03) / (x04 - 2 * x03 + x02);
    double a = -3 + c;
    double b = 4 - 2 * c;

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;

    return a * x4 + b * x3 + c * x2;
}


/*
 * Compares two complex numbers lexicographically. Returns
 * 0 if both are equal, -1 if the first is smaller and 1
 * if the second is smaller.
 */
int compare_complex(complex double a, complex double b) {
    if (a == b) {
        return 0;
    } else if (creal(a) < creal(b)) {
        return -1;
    } else if (creal(b) < creal(a)) {
        return 1;
    } else if (cimag(a) < cimag(b)) {
        return -1;
    }
    return 1;
}


/**
 * Converts convergence data to RGB data.
 *
 * Arguments:
 *   r: Array of roots algorithm converged to for that cell
 *   roots: List of roots
 *   maxr: Number of roots
 *   s: Step count array
 *   maxs: Maximal number of steps required for convergence
 *   avgs: Average number of steps required for convergence
 *   red: Pointer to array for red component
 *   green: Pointer to array for green component
 *   blue: Pointer to array for blue component
 */
void convert_to_rgb(int **r, complex double *roots, int maxr, int **s,
                    int maxs, float avgs, unsigned char ***red,
                    unsigned char ***green, unsigned char ***blue) {

    // Sort the roots array to ensure each root gets the same color every time
    int *root_index;
    MALLOC(root_index, maxr + 1)
    for (int i = 0; i <= maxr; i++) {
        root_index[i] = 0;
        for (int j = 0; j <= maxr; j++) {
            if (i != j && compare_complex(roots[i], roots[j]) == 1) {
                root_index[i] = root_index[i] + 1;
            }
        }
    }

    double a, b, c;
    for (int i = 0; i < IMAGE_SIZE; i++) {
        for (int j = 0; j < IMAGE_SIZE; j++) {
            if (r[i][j] == 0) {
                (*red)[i][j] = 0;
                (*green)[i][j] = 0;
                (*blue)[i][j] = 0;
            } else {
                #ifdef SINGLE_COLOR
                    a = 0;
                    b = 0;
                    c = 1 - 0.75 * blend((s[i][j] - 1) / ((double) maxs),
                                         (avgs - 1) / (maxs), 0.5);
                #else
                    a = root_index[r[i][j] - 1] / ((double) maxr + 1);
                    b = 0.75 * blend((s[i][j] - 1) / ((double) maxs),
                                     (avgs - 1) / (maxs), 0.5);
                    c = 1;
                #endif
                hsv2rgb(&a, &b, &c);
                (*red)[i][j] = (unsigned char) (255 * a);
                (*green)[i][j] = (unsigned char) (255 * b);
                (*blue)[i][j] = (unsigned char) (255 * c);
            }
        }
    }

    free(root_index);
}


/*
 * TODO: Allow the specification of options at the command line
 */
int main(int argc, char *argv[]) {
    printf("Initializing memory\n");
    double *x = linspace(x1, x2, IMAGE_SIZE);
    double *y = linspace(y1, y2, IMAGE_SIZE);
    int **r;
    int **s;
    unsigned char **red;
    unsigned char **green;
    unsigned char **blue;
    int maxr, maxs;
    MALLOC(r, IMAGE_SIZE)
    MALLOC(s, IMAGE_SIZE)
    MALLOC(red, IMAGE_SIZE)
    MALLOC(green, IMAGE_SIZE)
    MALLOC(blue, IMAGE_SIZE)
    for (int k = 0; k < IMAGE_SIZE; k++) {
        MALLOC(r[k], IMAGE_SIZE)
        MALLOC(s[k], IMAGE_SIZE)
        MALLOC(red[k], IMAGE_SIZE)
        MALLOC(green[k], IMAGE_SIZE)
        MALLOC(blue[k], IMAGE_SIZE)
    }
    complex double *roots;

    printf("Rastering\n");
    time_t start_time = time(NULL);
    raster(&METHOD, &f, &df, &ddf, 1e-10, MAX_STEPS, x, IMAGE_SIZE, y,
           IMAGE_SIZE, r, &maxr, s, &maxs, &roots);
    double duration = difftime(time(NULL), start_time);

    if (maxr == -1) {
        // No convergence anywhere
        fputs("Could not find any roots.\n", stderr);
        exit(EXIT_FAILURE);
    }

    long tsc = 0;  // Total step count for converged pixels
    long ts = 0;  // Total step count
    long conv = 0;  // Number of pixels where algorithm converged
    for (int i = 0; i < IMAGE_SIZE; i++) {
        for (int j = 0; j < IMAGE_SIZE; j++) {
            ts += s[i][j];
            if (s[i][j] < MAX_STEPS) {
                tsc += s[i][j];
                conv++;
            }
        }
    }
    double avgs = tsc / ((double) conv);

    printf("Found %d roots\n", maxr + 1);
    printf("Algorithm converged on %d of %d pixels (%.2f%%)\n", (int) conv,
           (int) (IMAGE_SIZE * IMAGE_SIZE),
           100 * (conv / ((double) IMAGE_SIZE * IMAGE_SIZE)));
    printf("Total number of steps: %d\n", (int) ts);
    printf("Total number of steps until convergence: %d\n", (int) tsc);
    printf("Average number of steps until convergence: %.4f\n", avgs);
    printf("Maximal number of steps until convergence: %d\n", maxs);
    printf("Seconds for rastering: %.2f (%e per pixel, %e per step)\n",
           duration, duration / (IMAGE_SIZE * IMAGE_SIZE), duration / ts);

    if (maxs == MAX_STEPS - 1) {
        printf("WARNING: MAX_STEPS might be set too low\n");
    }

    printf("Converting to colors\n");
    convert_to_rgb(r, roots, maxr, s, maxs, avgs, &red, &green, &blue);

    printf("Writing image file\n");
    if (!write_tga(red, green, blue, IMAGE_SIZE, IMAGE_SIZE, "fractal.tga")) {
        fputs("Could not write image file.\n", stderr);
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}
