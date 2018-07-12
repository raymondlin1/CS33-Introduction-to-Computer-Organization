//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Raymond Lin
 * UCLA ID: 304937942
 * Email id: superlintendo@g.ucla.edu
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int OMP_xMax;
#define xMax OMP_xMax
int OMP_yMax;
#define yMax OMP_yMax
int OMP_zMax;
#define zMax OMP_zMax

int BLOCK_SIZE = 8;

int OMP_Index(int x, int y, int z)
{
  return ((z * yMax + y) * xMax + x);
}
#define Index(x, y, z) OMP_Index(x, y, z)

double OMP_SQR(double x)
{
  return pow(x, 2.0);
}
#define SQR(x) OMP_SQR(x)

double* OMP_conv;
double* OMP_g;

void OMP_Initialize(int xM, int yM, int zM)
{
  xMax = xM;
  yMax = yM;
  zMax = zM;
  assert(OMP_conv = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
  assert(OMP_g = (double*)malloc(sizeof(double) * xMax * yMax * zMax));
}
void OMP_Finish()
{
  free(OMP_conv);
  free(OMP_g);
}
void OMP_GaussianBlur(double *u, double Ksigma, int stepCount)
{
  double lambda = (Ksigma * Ksigma) / (double)(2 * stepCount);
  double nu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
  int x, y, z, step;
  int x1, y1, z1, xlim, ylim, zlim;
  double boundryScale = 1.0 / (1.0 - nu);
  double postScale = pow(nu / lambda, (double)(3 * stepCount));
  
  for(step = 0; step < stepCount; step++)
    {
      for(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(y = 0; y < yMax; y+=BLOCK_SIZE)
	    {
	      zlim = z + BLOCK_SIZE;
	      for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		{
		  ylim = y + BLOCK_SIZE;
		  for(y1 = y; y1 < ylim && y1 < yMax; y1++)
		    {
		       u[Index(0, y1, z1)] *= boundryScale;
		    }
		}
	    }
	}
      /*LATESTfor(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(y = 0; y < yMax; y+=BLOCK_SIZE)
	    {
	      for( x = 1; x < xMax; x+=BLOCK_SIZE)
		{
		  zlim = z + BLOCK_SIZE;
		  for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		    {
		      ylim = y + BLOCK_SIZE;
		      for(y1 = y; y1 < ylim && y1 < yMax; y1++)
			{
			  xlim = x + BLOCK_SIZE;
			  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
			    {
			      u[Index(x1, y1, z1)] += u[Index(x1 - 1, y1, z1)] * nu;
			    }
			}
		    }
		}
	    }
	}*/
      for(z = 0; z < zMax; z+=2)
	{
	  for(y = 0; y < yMax; y++)
	    {
	      u[Index(0, y, z)] *= boundryScale;
	      u[Index(0, y, z + 1)] *= boundryScale;
	    }
	}
      /*for(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(y = 0; y < yMax; y+=BLOCK_SIZE)
	    {
	      for(x = xMax - 2; x >= 0; x-=BLOCK_SIZE)
		{
		  zlim = z + BLOCK_SIZE;
		  for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		    {
		      ylim = y + BLOCK_SIZE;
		      for(y1 = y; y1 < ylim && y1 < yMax; y1++)
			{
			  xlim = x + BLOCK_SIZE;
			  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
			    {
			      u[Index(x1, y1, z1)] += u[Index(x1 + 1, y1, z1)] * nu;
			    }
			}
		    }
		}
	    }
	    }*/
      for(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(x = 0; x < xMax; x+=BLOCK_SIZE)
	    {
	      zlim = z + BLOCK_SIZE;
	      for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		{
		  xlim = x + BLOCK_SIZE;
		  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
		    {
		      u[Index(x1, 0, z1)] *= boundryScale;
		    }
	        }
	    }
	}
      /*for(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(y = 1; y < yMax; y+=BLOCK_SIZE)
	    {
	      for(x = 0; x < xMax; x+=BLOCK_SIZE)
		{
		  zlim = z + BLOCK_SIZE;
		  for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		    {
		      ylim = y + BLOCK_SIZE;
		      for(y1 = y; y1 < ylim && y1 < yMax; y1++)
			{
			  xlim = x + BLOCK_SIZE;
			  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
			    {
			      u[Index(x, y, z)] += u[Index(x, y - 1, z)] * nu;
			    }
			}
		    }
		}
	    }
	}*/
      for(z = 0; z < zMax; z+=2)
	{
	  for(x = 0; x < xMax; x++)
	    {
	      u[Index(x, yMax - 1, z)] *= boundryScale;
	      u[Index(x, yMax - 1, z + 1)] *= boundryScale;
	    }
	}
      /*for(z = 0; z < zMax; z+=BLOCK_SIZE)
	{
	  for(y = yMax - 2; y >= 0; y-=BLOCK_SIZE)
	    {
	      for(x = 0; x < xMax; x+=BLOCK_SIZE)
		{
		  zlim = z + BLOCK_SIZE;
		  for(z1 = z; z1 < zlim && z1 < zMax; z1++)
		    {
		      ylim = y - BLOCK_SIZE;
		      for(y1 = y; y1 > ylim && y1 >= 0; y1--)
			{
			  xlim = x + BLOCK_SIZE;
			  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
			    {
			      u[Index(x, y, z)] += u[Index(x, y + 1, z)] * nu;
			    }
			}
		    }
		}
	    }
	}*/
      for(y = 0; y < yMax; y++)
	{
	  for(x = 0; x < xMax; x++)
	    {
	      u[Index(x, y, 0)] *= boundryScale;
	    }
	}
      for(z = 1; z < zMax - 1; z+=BLOCK_SIZE)
	{
	  for(y = 0; y < yMax; y+=BLOCK_SIZE)
	    {
	      for(x = 0; x < xMax; x+=BLOCK_SIZE)
		{
		  zlim = z + BLOCK_SIZE;
		  for(z1 = z; z1 < zlim && z1 < zMax - 1; z1++)
		    {
		      ylim = y + BLOCK_SIZE;
		      for(y1 = y; y1 < ylim && y1 < yMax; y1++)
			{
			  xlim = x + BLOCK_SIZE;
			  for(x1 = x; x1 < xlim && x1 < xMax; x1++)
			    {
			      u[Index(x1, y1, z1)] = u[Index(x1, y1, z1 - 1)] * nu;
			    }
			}
		    }
		}
	    }
	 }
      for(y = 0; y < yMax; y++)
	{
	  for(x = 0; x < xMax; x++)
	    {
	      u[Index(x, y, zMax - 1)] *= boundryScale;
	    }
	}
      /*for(z = zMax - 2; z > 1; z-=2)
	{
	  for(y = 0; y < yMax; y++)
	    {
	      for(x = 0; x < xMax; x++)
		{
		  u[Index(x, y, z)] += u[Index(x, y, z + 1)] * nu;
		  u[Index(x, y, z - 1)] += u[Index(x, y, z)] * nu;
		}
	    }
	    }*/
    }
  for(z = 0; z < zMax; z+=2)
    {
      for(y = 0; y < yMax; y++)
	{
	  for(x = 0; x < xMax; x++)
	    {
	      u[Index(x, y, z)] *= postScale;
	      u[Index(x, y, z + 1)] *= postScale;
	    }
	}
    }
}
void OMP_Deblur(double* u, const double* f, int maxIterations, double dt, double gamma, double sigma, double Ksigma)
{
  double epsilon = 1.0e-7;
  double sigma2 = SQR(sigma);
  int x, y, z, iteration;
  int converged = 0;
  int lastConverged = 0;
  int fullyConverged = (xMax - 1) * (yMax - 1) * (zMax - 1);
  double* conv = OMP_conv;
  double* g = OMP_g;

  for(iteration = 0; iteration < maxIterations && converged != fullyConverged; iteration++)
    {
      for(z = 1; z < zMax - 1; z++)
	{
	  for(y = 1; y < yMax - 1; y++)
	    {
	      for(x = 1; x < xMax - 1; x++)
		{
		  g[Index(x, y, z)] = 1.0 / sqrt(epsilon + 
						 SQR(u[Index(x, y, z)] - u[Index(x + 1, y, z)]) + 
						 SQR(u[Index(x, y, z)] - u[Index(x - 1, y, z)]) + 
						 SQR(u[Index(x, y, z)] - u[Index(x, y + 1, z)]) + 
						 SQR(u[Index(x, y, z)] - u[Index(x, y - 1, z)]) + 
						 SQR(u[Index(x, y, z)] - u[Index(x, y, z + 1)]) + 
						 SQR(u[Index(x, y, z)] - u[Index(x, y, z - 1)]));
		}
	    }
	}
      memcpy(conv, u, sizeof(double) * xMax * yMax * zMax);
      OMP_GaussianBlur(conv, Ksigma, 1);
      for(z = 0; z < zMax; z++)
	{
	  for(y = 0; y < yMax; y++)
	    {
	      for(x = 0; x < xMax; x++)
		{
		  double r = conv[Index(x, y, z)] * f[Index(x, y, z)] / sigma2;
		  r = (r * (2.38944 + r * (0.950037 + r))) / (4.65314 + r * (2.57541 + r * (1.48937 + r)));
		  conv[Index(x, y, z)] -= f[Index(x, y, z)] * r;
		}
	    }
	}
      OMP_GaussianBlur(conv, Ksigma, 1);
      converged = 0;
      for(z = 1; z < zMax - 1; z++)
	{
	  for(y = 1; y < yMax - 1; y++)
	    {
	      for(x = 1; x < xMax - 1; x++)
		{
		  double oldVal = u[Index(x, y, z)];
		  double newVal = (u[Index(x, y, z)] + dt * ( 
							     u[Index(x - 1, y, z)] * g[Index(x - 1, y, z)] + 
							     u[Index(x + 1, y, z)] * g[Index(x + 1, y, z)] + 
							     u[Index(x, y - 1, z)] * g[Index(x, y - 1, z)] + 
							     u[Index(x, y + 1, z)] * g[Index(x, y + 1, z)] + 
							     u[Index(x, y, z - 1)] * g[Index(x, y, z - 1)] + 
							     u[Index(x, y, z + 1)] * g[Index(x, y, z + 1)] - gamma * conv[Index(x, y, z)])) /
		    (1.0 + dt * (g[Index(x + 1, y, z)] + g[Index(x - 1, y, z)] + g[Index(x, y + 1, z)] + g[Index(x, y - 1, z)] + g[Index(x, y, z + 1)] + g[Index(x, y, z - 1)]));
		  if(fabs(oldVal - newVal) < epsilon)
		    {
		      converged++;
		    }
		  u[Index(x, y, z)] = newVal;
		}
	    }
	}
      if(converged > lastConverged)
	{
	  printf("%d pixels have converged on iteration %d\n", converged, iteration);
	  lastConverged = converged;
	}
    }
}

