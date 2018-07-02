#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <gsl/gsl_linalg.h>
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

////////////////////////////////////////////////////////////////////////////////
////////////////////// Rho = Mesh density function /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

   void Rho(double *y, double *rho_sm)

    {
  
     int neq=38;

      int N=neq/2;
        
      vector<double> u(N);
      vector<double> x(N);

            for (int j=0;j<N;j++) 
             {
             u[j] =y[j];
             x[j] =y[j+N];
             } 

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
   
   double gamma = 1.0/3.0;
   
   double Alpha =0.0;  // eps;

   for (int j=1;j<N;j++)

       {
       Alpha=Alpha+ (1.0/2.0)*(pow(rho[j],gamma) + pow(rho[j-1],gamma))*(x[j]-x[j-1]);
       }

   Alpha = pow(Alpha,3);

   
   vector<double> rh(neq);

     for (int j=0;j<N;j++)

       {
       rh[j]=pow( (1.0+(1.0/Alpha)*rho[j]),1.0/3.0);  //  rho = (1+(1/alpha)*rho).^(1/3);
       }
   
  // smoothing mesh density function
   
  for (int j=1;j<N-1;j++)
       {
       rho_sm[j] = 1.0/4.0*(rh[j-1]+rh[j+1])+1.0/2.0*rh[j];  // 1/4*(rho(j-1)+rho(j+1))+1/2*rho(j);
       }
   rho_sm[0] = 1.0/2.0*(rh[0]+rh[1]);
   
   rho_sm[N-1] = 1.0/2.0*(rh[N-1]+rh[N-2]);
  
 }   
////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Mass = mass matrix /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void mass(double t, double *y, double **M)
   
   {
	 int neq=38;
	   
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

      vector< vector<double> > M1(N,vector<double>(N));
      vector< vector<double> > M2(N,vector<double>(N));
      vector< vector<double> > M3(N,vector<double>(N));
      vector< vector<double> > M4(N,vector<double>(N));
      
      
      // M1
         for (int i=0;i<N;i++)
             {
             M1[i][i]=1;
             }
      // M2
        M2[0][0]= -(u[1] - u0)/(x[1] - x0);

         M2[N-1][N-1]=- (uNP1 - u[N-2])/(xNP1 - x[N-2]);
         
         for (int i=0;i<N-1;i++)
             {
             M2[i][i]= - (u[i+1] - u[i-1])/(x[i+1] - x[i-1]);
             }

      // M3

       for (int i=0;i<N;i++) 
             {
             for (int j=0;j<N;j++)
               {
               M3[i][j]=0;
               }
             }

       // M4

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
    
   }

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// gauss elimination //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void gauss(double *F, double **J, double *xn)

{
neq=38;	
   
double *Jn = new double[neq*neq];
  
   
  for(int i=0;i<neq;i++) {
    for(int j=0;j<neq;j++)
      Jn[i*neq + j]=J[i][j];
  }
  
     	  
  gsl_matrix_view m = gsl_matrix_view_array (Jn, neq, neq);

  gsl_vector_view b = gsl_vector_view_array (F, neq);

  gsl_vector *x = gsl_vector_alloc (neq);
  
  int s;

  gsl_permutation * p = gsl_permutation_alloc (neq);

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

   
  for (int i=0;i<neq;i++)
      {
	xn[i]=x->data[i];  
	  }	  
  
 // printf ("x = \n");
 // gsl_vector_fprintf (stdout, x, "%g");

  gsl_permutation_free (p);
  gsl_vector_free (x);
  
  delete [] Jn;

}

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////// rhs of the ode y'=L^{-1}g(t,y) /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

void rhs(double t, double *y, double *f)

      {

     int neq=38;

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
   
     Rho(y,rho_sm); 

     double *v = new double[N];

     for (int i=1;i<(N-1);i++)
         {
         v[i] =( (rho_sm[i+1] + rho_sm[i])*(x[i+1] - x[i]) - (rho_sm[i] + rho_sm[i-1])*(x[i] - x[i-1]) );
         }

    v[0] =  ( (rho_sm[1] + rho_sm[0])*(x[1] - x[0]) - (rho_sm[0] + rho_sm[0])*(x[0] - x0) );
    
    v[N-1] =( (rho_sm[N-1] + rho_sm[N-1])*(xNP1 - x[N-1]) - (rho_sm[N-1] + rho_sm[N-2])*(x[N-1] - x[N-2]) );

    for (int i=0;i<N;i++)
        {
        g[i+N]=-1.0/(2.0*tau)*v[i];  
        }
    
      double **L = new double*[neq];

         for (int j=0;j<neq;j++) 
             {
             L[j] = new double[neq];
             }
          
       mass(t,y,L);

      double *w = new double[neq];

      gauss(g,L,w);

       for (int j=0;j<neq;j++)

         {
         f[j]=w[j];
         }

     for (int j=0; j<neq; j++)
      {
       delete [] L[j];
      }
  
    delete [] L;

    delete [] g;
    
    delete [] rho_sm;

    delete w;

    delete v; 
    
    delete u;

    delete x; 
  
   }

/////////////////////////////////////////////////////////////////////////////////////
////////////////////  function in the form F(x)=0  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

    void fun(double t, double *p, double *u, double *Fn)

      {
		double dt=0.01;
		  
		int neq=38;

        double* frh = new double[neq];

           rhs(t,p,frh);
   
        for (int j=0;j<neq;j++)

               {
               Fn[j]=p[j]-u[j]-dt*frh[j];
               }

      delete [] frh;  

   } 

////////////////////////////////////////////////////////////////////////////////
//////////////////// Numerical jacobian of F(x)=0 //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void jac(double t, double *x, double *xold, double **J)
  {
  int neq=38;	  

  double dx=0.00000001;
  
  vector<double> xx(neq);
   
  for (int i=0;i<neq;i++)

    {
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

    fun(t,x,xold,fx);

   double *fxx = new double[neq];

    fun(t,xn,xold,fxx);

       for (int k=0;k<neq;k++)
      
        {
         J[k][i]=(fxx[k]-fx[k])/dx;
        }

    delete [] fx;
   
    delete [] xn;

    delete [] fxx;
    
  }
  
  vector<double>().swap(xx);
 

  }

////////////////////////////////////////////////////////////////////////////////
////////////////////////// l2 norm of vector ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void l2_norm(double *w, int n,double *norm)
        {
         double accum = 0.0;
         for (int i = 0; i < n; ++i)
          {
            accum += w[i] * w[i];
          }
         norm[0]=sqrt(accum);

        }

////////////////////////////////////////////////////////////////////////////////
/////////////////////////// BE with Newton's Solver ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void euler_be(double *u, double *un, int neq, int nt, double dt)

   {
	   
	double ti=0.0;
	   
   //////// setting parameters for Newton's iteration  ////////// 

          int max_iter=100;

          double tol=0.00000001;
          
///////////////////////////////////////////////////////////////
////////////////////// empty vectors //////////////////////////
///////////////////////////////////////////////////////////////
    
     vector<double> z(neq);
     
     vector<double> p(neq);
     
     vector<double> u_st(neq);
     
     
///////////////////////////////////////////////////////////////
////////////////////// initial guess //////////////////////////
///////////////////////////////////////////////////////////////
     
        for (int j=0;j<neq;j++)
              {
               p[j]=u[j];
              }
              
///////////////////////////////////////////////////////////////
///////////////////// Time loop start ////////////////////////
///////////////////////////////////////////////////////////////             
	
	for (int i=0;i<nt;i++)

        {   
		  
         double t=ti+i*dt;
	 
        double *w = new double[neq];   
        
        for (int j=0;j<neq;j++)
            {
             w[j]=p[j];
            }
/////////////////////////////////////////////////////////////////
/////////////////////// Newton's iteration starts ///////////////
///////////////////////////////////////////////////////////////// 
     
        for (int l=0;l<max_iter;l++)
            {
			       
             double *v = new double[neq];

             for (int j=0;j<neq;j++)
                 {
                 v[j]=p[j];
                 }
             
             double *Fn = new double[neq]; 

             fun(t,v,w,Fn);                // calling function F(X)=0
             
             double **Jn = new double*[neq];

               for (int j=0;j<neq;j++) 
                  {
                  Jn[j] = new double[neq];
                  }
             
              jac(t,v,w,Jn);               // calling jacobian 

         double *x=new double[neq];

         gauss(Fn,Jn,x);                  // calling linear solver for x=Jn/Fn
         
         for (int j=0;j<neq;j++)

             {
              z[j]=v[j]-x[j];
             }
             
         double *dv=new double[neq];

         for (int j=0;j<neq;j++)

            {
            dv[j]=z[j]-v[j];
            }

        double *norm=new double[1];
       
        l2_norm(dv,neq,norm);
       
          if (norm[0]<tol)

            {
		  //  printf("converged\n");		
		    break;
            }
            
            p=z;

////////////////////////////////////////////////////////////////////            
//////////// deconstructor /////////////////////////////////////////
////////////////////////////////////////////////////////////////////            

         delete [] x;

         delete [] Fn;

          delete [] v;
          
          delete [] dv;
          
          delete [] norm;

          for (int i=0; i<neq; i++)
            {
            delete [] Jn[i];
            }
           delete [] Jn;

     } 

/////////////////////////////////////////////////////////////////
/////////////////////// Ends of the Newton's iteration //////////
/////////////////////////////////////////////////////////////////        
       
       delete [] w;
       
       u_st=z;
      
   } 
///////////////////////////////////////////////////////////////
///////////////////// End of the time loop ///////////////////////
/////////////////////////////////////////////////////////////// 
     
  
   for (int j=0;j<neq;j++)
      {
      un[j]=u_st[j];
      } 
 
} 

int main(int argc, char *argv[]) 

 {
  int start_s=clock();

  int neq = 38;
  int ti = 0;
  int tf = 1;
  
  double dt=0.01; 
  
  int nt=100;

   int N=neq/2;

   double h = 1.0/(double(N)+1);

   double *xint = new double[N];

   double *uintv = new double[N];
   
   double *u = new double[neq];
   double *un = new double[neq];

    for (int i=0;i<N;i++)
       {
       xint[i]=h*(i+1);
       uintv[i]=sin(2*M_PI*xint[i]) + (1.0/2.0)*sin(M_PI*xint[i]);
       }

    for (int i=0;i<N;i++)
       {
       u[i]=uintv[i];
       u[i+N]=xint[i];
       }
       
   // euler_be(u, un, neq, nt, dt);   
    
    double **Jc = new double*[neq];

               for (int j=0;j<neq;j++) 
                  {
                  Jc[j] = new double[neq];
                  }
                  
     mass(0.01,u,Jc); 
     double *fv = new double[neq];
     double *xv = new double[neq];
     rhs(0.01,u,fv);
     gauss(fv,Jc,xv);
     
                   
            
            for (int i=0; i<neq; i++)
            {
            delete [] Jc[i];
            }
           delete [] Jc;
      
    
      
     for (int i = 0; i < neq; i++)   {
     printf("%16.14f\n", xv[i]);     }
    
   delete [] xint;
   
   delete [] uintv;
   
   delete [] u;
   
   delete [] un;
   
   delete [] fv;
   
   delete [] xv;
 
  printf("Time taken: %.12fs\n", (double)(clock() - start_s)/CLOCKS_PER_SEC); 

}


