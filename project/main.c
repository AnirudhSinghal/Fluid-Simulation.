#include <math.h>
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdlib.h>
int N = 840;
int iter = 6;

#define IX( x, y) (( x) + ( y) * N )


struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;
    float *s;
    float *density;
    float *Vx;
    float *Vy;
    float *Vx0;
    float *Vy0;
};
typedef struct FluidCube FluidCube;



FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, float dt)
{
    FluidCube *cube = malloc(sizeof(*cube));
    N = size;
    
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    
    cube->s = calloc(N * N , sizeof(float));
    cube->density = calloc(N * N , sizeof(float));
    
    cube->Vx = calloc(N * N , sizeof(float));
    cube->Vy = calloc(N * N , sizeof(float));
    
    cube->Vx0 = calloc(N * N , sizeof(float));
    cube->Vy0 = calloc(N * N , sizeof(float));
    
    return cube;
}

void FluidCubeFree(FluidCube *cube)
{
    free(cube->s);
    free(cube->density);
    
    free(cube->Vx);
    free(cube->Vy);
    
    free(cube->Vx0);
    free(cube->Vy0);
    
    free(cube);
}


static void set_bnd(int b, float *x, int N){

    for ( int i = 1 ; i <= N ; i++ ) {
        x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
        x[IX(N-1,i)] = b==1 ? -x[IX(N-1,i)] : x[IX(N-1,i)];
        x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N-1)] = b==2 ? -x[IX(i,N-1)] : x[IX(i,N-1)];
    }
    x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
    x[IX(0  ,N)] = 0.5f*(x[IX(1,N)]+x[IX(0  ,N-1)]);
    x[IX(N,0  )] = 0.5f*(x[IX(N-1,0  )]+x[IX(N,1)]);
    x[IX(N,N)] = 0.5f*(x[IX(N,N-1)]+x[IX(N,N)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j)] =
                        (x0[IX(i, j)]
                            + a*(    x[IX(i+1 , j   )]
                                    +x[IX(i  , j+1  )]
                                    +x[IX(i  , j-1  )]
                                    +x[IX(i - 1  , j  )]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
}

static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j, k;
    
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j)];
                tmp2 = dty * velocY[IX(i, j)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5f) x = 0.5f; 
                if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                i0 = floorf(x); 
                i1 = i0 + 1.0f;
                if(y < 0.5f) y = 0.5f; 
                if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                j0 = floorf(y);
                j1 = j0 + 1.0f; 
               
                
                s1 = x - i0; 
                s0 = 1.0f - s1; 
                t1 = y - j0; 
                t0 = 1.0f - t1;
        
                
                int i0i = (int) i0;
                int i1i = (int) i1;
                int j0i = (int) j0;
                int j1i = (int) j1;
    
                
                d[IX(i, j)] = 
                
                    s0 * ( t0 * (d0[IX(i0i, j0i)])
                                
                        +( t1 * d0[IX(i0i, j1i)]))

                   +s1 * ( t0 * ( d0[IX(i1i, j0i)])
                                
                        +( t1 * d0[IX(i1i, j1i)]));
                               
            }
        }
    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N)
{
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j)] = -0.5f*(
                         velocX[IX(i+1, j    )]
                        -velocX[IX(i-1, j   )]
                        +velocY[IX(i  , j+1  )]
                        -velocY[IX(i  , j-1  )]
                    )/N;
                p[IX(i, j)] = 0;
            }
        }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                                -p[IX(i-1, j)]) * N;
                velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                                -p[IX(i, j-1)]) * N;
            }
        }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

void FluidCubeStep(FluidCube *cube)
{
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;
    float dt       = cube->dt;
    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *s       = cube->s;
    float *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);

    project(Vx0, Vy0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);
    
    project(Vx, Vy, Vx0, Vy0, 4, N);
    
    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, dt, N);
}

void FluidCubeAddDensity(FluidCube *cube, int x, int y, int amount )
{
    int N = cube->size;
    cube->density[IX(x, y)] += amount;
}

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, float amountX, float amountY)
{
    int N = cube->size;
    int index = IX(x, y);
    
    cube->Vx[index] += amountX;
    cube->Vy[index] += amountY;
}

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 640
#define GRID_ROWS 64
#define GRID_COLS 64
#define SQUARE_SIZE (SCREEN_WIDTH / GRID_COLS)

SDL_Window* gWindow = NULL;
SDL_Renderer* gRenderer = NULL;

void drawSquare(int row, int col, int brightness) {
    SDL_Rect squareRect = {col * SQUARE_SIZE, row * SQUARE_SIZE, SQUARE_SIZE, SQUARE_SIZE};
    SDL_RenderFillRect(gRenderer, &squareRect);
}

SDL_Color mapDensityToColor(float density) {
    SDL_Color color;
    color.r = density; 
    color.g = density;
    color.b = density;
    color.a = 255; 
    return color;
}

int main(int argc, char* args[]) {

   
    float deltaTime = 0.1;
    FluidCube *cube = FluidCubeCreate(N , 0.00001 , 0.00005 , deltaTime);
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
        return 1;
    }

    gWindow = SDL_CreateWindow("Grid Hover Effect", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
    if (gWindow == NULL) {
        printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
        return 2;
    }

    gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);
    if (gRenderer == NULL) {
        printf("Renderer could not be created! SDL Error: %s\n", SDL_GetError());
        return 3;
    }


    int quit = 0;
    SDL_Event e;


    while (!quit) {

           while (SDL_PollEvent(&e) != 0) {
                if (e.type == SDL_QUIT) {
                    quit = 1;
                }
                else if (e.type == SDL_MOUSEMOTION) {
                float PmouseX = e.motion.x / SQUARE_SIZE;
                float PmouseY = e.motion.y / SQUARE_SIZE;
                int mouseX, mouseY;
                SDL_GetMouseState(&mouseX, &mouseY);
                int PemouseX = mouseX / SQUARE_SIZE, PemouseY = mouseY / SQUARE_SIZE ;
                FluidCubeAddDensity(cube , PemouseX , PemouseY , 255 );
                FluidCubeAddVelocity(cube , PemouseX , PemouseY , PmouseX, PmouseY);
                }
         
        }

        SDL_SetRenderDrawColor(gRenderer, 0, 0, 0, 255);
        SDL_RenderClear(gRenderer);

        FluidCubeAddDensity(cube , 100 , 146 , 600);
        FluidCubeAddVelocity(cube , 100 , 150 , 0 , -10);
        FluidCubeStep(cube);


        for (int row = 0; row < GRID_ROWS; row++) {
            for (int col = 0; col < GRID_COLS; col++) {
            SDL_Color squareColor = mapDensityToColor((cube->density[IX(col, row)]>255) ? 255 : cube->density[IX(col , row)]);
            SDL_SetRenderDrawColor(gRenderer, squareColor.r, squareColor.g, squareColor.b, squareColor.a);
            drawSquare(row, col, cube->density[IX(col,row)]);
            }

            

        

        }
        SDL_RenderPresent(gRenderer);
    }
    FluidCubeFree(cube);

    SDL_DestroyRenderer(gRenderer);
    SDL_DestroyWindow(gWindow);
    gRenderer = NULL;
    gWindow = NULL;

    SDL_Quit();

    return 0;
}
