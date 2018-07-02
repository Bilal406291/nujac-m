#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
//#include <gsl/include/gsl_linalg.h>
using namespace  std;

////////////////////////////////////////////////////////////////////////////////
////////////////////// diff = difference of two vectors ////////////////////////
////////////////////////////////////////////////////////////////////////////////

   void diff(double *yn, double *xn, double *vn,int N)
           {

           for (int i=0;i<N;i++)
            {
            vn[i]=(yn[i]-xn[i]);
            }
           }


void gauss(double *F, double **J, double *xn, int neq)

  {
   
   // creating augmented matrix

        double **A = new double*[neq];

         for (int j=0;j<neq;j++) 
             {
             A[j] = new double[neq+1];
             }
 
      // Read input data

             for (int j=0; j<neq; j++)
                {
              for (int k=0; k<neq; k++)
                  {
                 A[j][k]=J[j][k];
                  }
                }

            for (int j=0; j<neq; j++)
                {
                A[j][neq]=F[j];
                }

      // solving AX=B 

   for (int i=0; i<neq; i++)
      {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<neq; k++)
          {
            if (abs(A[k][i]) > maxEl)
            {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
          }

        // Swap maximum row with current row (column by column)

        for (int k=i; k<neq+1;k++)
           {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
           }

        // Make all rows below this one 0 in current column

        for (int k=i+1; k<neq; k++)
        {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<neq+1; j++)
            {
                if (i==j)
                {
                    A[k][j] = 0;
                } else
                   {
                    A[k][j] += c * A[i][j];
                   }
            }
        }
     }  // i loop

    // Solve equation Ax=b for an upper triangular matrix A

   
    for (int i=neq-1; i>=0; i--)
      {
        xn[i] = A[i][neq]/A[i][i];
        for (int k=i-1;k>=0; k--)
        {
            A[k][neq] -= A[k][i] * xn[i];
        }
      }

    for (int i=0; i<neq; i++)
      {
       delete [] A[i];
      }
  delete [] A;
  
}


void Rho(double *y, double *rho_sm, int neq)

    {

      int N=neq/2;
        
      double *u = new double[N];

      double *x = new double[N];

            for (int j=0;j<N;j++) 
             {
             u[j] =y[j];
             x[j] =y[j+N];
             } 

    // double *rho = new double[N];

     //double *v = new double[N];

     vector<double> rho(N);

     vector<double> v(N);

     for (int j=1;j<N-1;j++)
         {
         v[j]=2.0/(x[j+1]-x[j-1])*( (u[j+1]-u[j])/(x[j+1]-x[j])-(u[j]-u[j-1])/(x[j]-x[j-1]) );
         }
     v[0] = 2.0*((x[1]-x[0])*(u[2]-u[0])-(x[2]-x[0])*(u[1]-u[0]))/((x[2]-x[0])*(x[1]-x[0])*(x[2]-x[1]));
            
     v[N-1] = 2.0*((x[N-2]-x[N-1])*(u[N-3]-u[N-1])-(x[N-3]-x[N-1])*(u[N-2]-u[N-1]))/((x[N-3]-x[N-1])*(x[N-2]-x[N-1])*(x[N-3]-x[N-2]));
     
     for (int j=0;j<N;j++)

       {
       rho[j]=rho[j]+pow(v[j],2);      // rho = rho + v.^2;
       }

                    
   // alpha calculation
   
   double gamma = 1.0/3;
   
   double Alpha =0.0;  // eps;

   for (int j=1;j<N;j++)

       {
       Alpha=Alpha+ (1.0/2)*(pow(rho[j],gamma) + pow(rho[j-1],gamma))*(x[j]-x[j-1]);
       }

   Alpha = pow(Alpha,3);

   
   //alpha=max(1,Alpha);

   //mesh density function

  // double *rh = new double[N];

   vector<double> rh(neq);

     for (int j=0;j<N;j++)

       {
       rh[j]=pow( (1.0+(1.0/Alpha)*rho[j]),1.0/3);  //  rho = (1+(1/alpha)*rho).^(1/3);
       }
   
  // smoothing mesh density function
   
  for (int j=1;j<N-1;j++)
       {
       rho_sm[j] = 1.0/4*(rh[j-1]+rh[j+1])+1.0/2*rh[j];  // 1/4*(rho(j-1)+rho(j+1))+1/2*rho(j);
       }
   rho_sm[0] = 1.0/2*(rh[0]+rh[1]);
   
   rho_sm[N-1] = 1.0/2*(rh[N-1]+rh[N-2]);

  delete u;
  delete x;
 // delete rho;
 // delete v;
 //  delete rh;

   }   


void mass(double *y, double **M, int neq)
   
   {
     int N=neq/2;

      double *u = new double[N];

      double *x = new double[N];

            for (int j=0;j<N;j++) 
             {
             u[j] =y[j];
             x[j] =y[j+N];
             } 

      /////// IC ///////

      double x0,u0,xNP1,uNP1;

      x0=0.0;
      u0=0.0;
      xNP1=1.0;
      uNP1=0.0;

      double **M1 = new double*[N];

         for (int j=0;j<N;j++) 
             {
             M1[j] = new double[N];
             }

         for (int i=0;i<N;i++)
            {
             M1[i][i]=1;
             }

      double **M2 = new double*[N];

         for (int j=0;j<N;j++) 
             {
             M2[j] = new double[N];
             }
         M2[0][0]= -(u[1] - u0)/(x[1] - x0);

         M2[N-1][N-1]=- (uNP1 - u[N-2])/(xNP1 - x[N-2]);
         
         for (int i=0;i<N-1;i++)
             {
             M2[i][i]= - (u[i+1] - u[i-1])/(x[i+1] - x[i-1]);
             }

      double **M3 = new double*[N];

         for (int j=0;j<N;j++) 
             {
             M3[j] = new double[N];
             }

       for (int i=0;i<N;i++) 
             {
             for (int j=0;j<N;j++)
               {
               M3[i][j]=0;
               }
             }

       double **M4 = new double*[N];

       for (int j=0;j<N;j++) 
             {
             M4[j] = new double[N];
             }

       M4[0][0]=-2;
       M4[0][1]=1;
       for (int i=1;i<N-1;i++)
          {
          M4[i][i]=-2;
          M4[i][i-1]=1;
          M4[i][i+1]=1;
          }
       M4[N-1][N-1]=-2;
       M4[N-1][N-2]=1;

        for (int i=0;i<N;i++) 
             {
             for (int j=0;j<N;j++)

               {
               M[i][j]=M1[i][j];
               M[i][j+N]=M2[i][j];
               M[i+N][j]=M3[i][j];
               M[i+N][j+N]=M4[i][j];
               }
            
             }

      delete [] u;

      delete [] x;

     for (int i=0; i<N; i++)
      {
       delete [] M1[i];
       delete [] M2[i];
       delete [] M3[i];
       delete [] M4[i];
      }
  
     delete [] M1;
     delete [] M2;
     delete [] M3;
     delete [] M4;   

     
   }


void rhs(double *y, double *f, int neq)

      {

      int N=neq/2;

      double ep, tau, x0,u0,xNP1,uNP1;

      ep=0.01;

      tau=0.01;

      double *u = new double[N];

      double *x = new double[N];

            for (int j=0;j<N;j++) 
             {
             u[j] =y[j];
             x[j] =y[j+N];
             } 

      /////// IC ///////

      x0=0.0;
      u0=0.0;
      xNP1=1.0;
      uNP1=0.0;

     ////////////////

      double *g = new double[neq];

      double dx;

      for (int i=1;i<(N-1);i++)
         {
         dx = x[i+1] - x[i-1];
         g[i] = (2.0*ep)/dx*( (u[i+1] - u[i])/(x[i+1] - x[i] ) - (u[i] - u[i-1])/(x[i] - x[i-1]) )- (1.0/2.0)*(pow(u[i+1],2)-pow(u[i-1],2))/dx;
         }

     dx = x[1] - x0;    
    
     g[0] = (2.0*ep)/dx*((u[1] - u[0])/(x[1] - x[0]) - (u[0] - u0)/(x[0] - x0)) - (1.0/2.0)*(pow(u[1],2) - pow(u0,2))/dx;

     dx = xNP1 - x[N-2];   
    
     g[N-1]=(2.0*ep)/dx*((uNP1 - u[N-1])/(xNP1 - x[N-1]) - (u[N-1] - u[N-2])/(x[N-1] - x[N-2]))/dx - (1.0/2.0)*(pow(uNP1,2) - pow(u[N-2],2))/dx;

     double *rho_sm = new double[N]; 
   
     Rho(y,rho_sm,neq); 

    //  for (int i=0;i<N;i++){
    //  printf("%16.12f\n",rho_sm[i]); }     

     double *w = new double[N];

     for (int i=1;i<(N-1);i++)
         {
         w[i] =( (rho_sm[i+1] + rho_sm[i])*(x[i+1] - x[i]) - (rho_sm[i] + rho_sm[i-1])*(x[i] - x[i-1]) );
         }

    w[0] =  ( (rho_sm[1] + rho_sm[0])*(x[1] - x[0]) - (rho_sm[0] + rho_sm[0])*(x[0] - x0) );
    
    w[N-1] =( (rho_sm[N-1] + rho_sm[N-1])*(xNP1 - x[N-1]) - (rho_sm[N-1] + rho_sm[N-2])*(x[N-1] - x[N-2]) );

    for (int i=0;i<N;i++)
        {
        g[i+N]=-1.0/(2.0*tau)*w[i];  
        }

     double **L = new double*[neq];

         for (int j=0;j<neq;j++) 
             {
             L[j] = new double[neq];
             }
          
      // matrix L(y) calling 

     mass(y,L,neq); 

     double *xn = new double[neq];

     gauss(g,L,xn,neq);

      for (int j=0;j<neq;j++) 
             {
             f[j] = xn[j];
             }

    
     delete [] u;

    delete [] x;

    delete [] g;
    
    delete [] rho_sm;

    delete [] w;

    delete [] xn;

    for (int i=0;i<neq;i++)
        {
        delete [] L[i];
        }
     delete [] L;
     

   }

/////////////////////////////////////////////////////////////////////////////////////
////////////////////  function in the form F(x)=0  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    void fun(double *p, double *u, double *Fn, int neq)

      {

        double* frh = new double[neq];

         rhs(p,frh,neq);

         double *v = new double[neq];

        // v=p-u

        diff(p,u,v,neq);

       
        for (int j=0;j<neq;j++)

               {
               Fn[j]=v[j]-dt*frh[j];
               }

      delete [] frh;  

      delete [] v;

     } 



////////////////////////////////////////////////////////////////////////////////
//////////////////// Numerical jacobian of F(x)=0 //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void jac(double *x, double **J, int neq)
  {

  double dx=0.00000001;

   
  for (int i=0;i<neq;i++)

    {

      vector<double> xx(neq);

      for (int j=0;j<neq;j++)

       {
       xx[j]=x[j];
       }

     xx[i]=x[i]+dx; 

    double *xn = new double[neq];
    
    for (int j=0;j<neq;j++)

       {
       xn[j]=xx[j];
       }
     
    double *fx = new double[neq];

    rhs(x,fx,neq); 
    
    double *fxx = new double[neq];

    rhs(xn,fxx,neq);

      for (int k=0;k<neq;k++)
      
        {
         J[k][i]=(fxx[k]-fx[k])/dx;
        }

    delete [] fx;
   
    delete [] xn;

    delete [] fxx;

    }
 

  }

int  main() 
  {
  int start_s=clock(); 
  int order, nt;
  double *sol;

  int neq = 38;
    
  int N=neq/2;

  double h = 1.0/(N+1);

  double *xint = new double[N];

  double *uint = new double[N];

   for (int i=0;i<N;i++)
       {
       xint[i]=h*(i+1);
       uint[i]=sin(2*M_PI*xint[i]) + (1.0/2.0)*sin(M_PI*xint[i]);
       }

   sol = new double[neq];

   for (int i=0;i<N;i++)
       {
       sol[i]=uint[i];
       sol[i+N]=xint[i];
       }

  // call ridc 

  double *f = new double[neq]; 

  rhs(sol,f,neq);

   

  //  Rho(sol,rh,neq);

    for (int i=0;i<neq;i++){
      printf("%16.12f\n",f[i]); }


  double **J =new double*[neq];
      for (int i=0;i<neq;i++)
          {
          J[i]=new double[neq];
          }

 // jac(sol,J,neq);

 // for (int i = 0; i < neq; i++) {
 //    for (int j=0;j< neq;j++)   {
 //      printf("%10.2f ", L[i][j]);  } 
  //    printf("\n");
  //    printf("\n");   }

for (int i=0;i<neq;i++)

    {
    delete [] J[i];
    }
  
  delete [] J;

printf("\nTime taken: %.12fs\n", (double)(clock() - start_s)/CLOCKS_PER_SEC); 

  delete [] sol;
  delete [] xint;
  delete [] uint;
  delete [] f;
  
 
 } 
