//code for RUNGA KUTTA 4-th ORDER
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define PI 3.141592653589793
#define B 2.
#define Wo 10.
#define M 0.1
#define Fo 1.
#define cte 0.4
#define TOL 1e-10

double g(double v);
double f(double x, double v, double t);
double solExacta(double t);
void rk_4(double h,int N, double t[N], double x[N],double v[N], double tn);



int main(){
    FILE *fp_x;
    FILE *fp_v;
    FILE *fp_ES;
    FILE *fp2;
    FILE *fp3;
    fp_x=fopen("RK4-points_x.txt","w"); //we create the txt file that contains the points of x forthe RK4 algorithm
    fp_v=fopen("RK4-points_v.txt","w"); //we create the txt file that contains the points of v for the RK4 algorithm
    fp_ES=fopen("RK4-points_ES.txt","w");//we create the txt file that contains the points of the analytical solution
    fp2=fopen("RK4Errors.txt","w"); //we create the txt file that contains the exact errors.
    fp3=fopen("RK4-rel-errors.txt","w"); //we create the txt file that contains the relative errors.
    printf("Introduce the value for h:\n");
    double h;
    scanf("%lf",&h);
    double tn=10.*PI/Wo;
    int N;
    double step=tn/h;
    N=(int)step+1;
    double t[N],x[N],v[N];
    x[0]=0.;
    v[0]=1.;
    t[0]=0.;
    

    rk_4(h,N,t,x,v,tn);
    int i;
    for(i=1;i<=(tn/h)+h;i++){//save points (t,x(t))
        fprintf(fp_x,"%.16G  \t %.16G \n",t[i],x[i]);
    }
    for(i=1;i<=(tn/h)+h;i++){//save points (t,v(t))
        fprintf(fp_v,"%.16G  \t %.16G \n",t[i],v[i]);
    }
    for(i=1;i<=(tn/h)+h;i++){//save points (t,solExacta(t))
        fprintf(fp_ES,"%.16G  \t %.16G \n",t[i],solExacta(t[i]));
    }
    for(i=1;i<=(tn/h)+h;i++){ //save errors
        fprintf(fp2,"%.16G  \t %.16G \n",t[i], fabs(solExacta(t[i])-x[i]));
    }
    for(i=1;i<=(tn/h)+h;i++){ //save errors
        fprintf(fp3,"%.16G  \t %.16G \n",t[i], fabs(solExacta(t[i])-x[i])/fabs(solExacta(t[i])));
    }
    fclose(fp_x);
    fclose(fp_v);
    fclose(fp_ES);
    fclose(fp2);
    fclose(fp3);
    return 0;
}


double g(double v){
    return v;
}

double f(double x, double v, double t){
    double w=cte*Wo;
    return -(B/M)*v-Wo*Wo*x+(Fo/M)*cos(w*t);
}

double solExacta(double t){
    if(B==0){ //simple harmonic oscillator
        return (1./Wo)*sin(Wo*t);
        
    }
    double alpha,theta;
    double w=cte*Wo;
    double inside_sqrt=B*B-4.*M*M*Wo*Wo;
    
    double inside_sqrt_2=M*M*(Wo*Wo-w*w)*(Wo*Wo-w*w)+B*B*w*w;
    if(inside_sqrt_2<=TOL){//should never enter here but just in case
        inside_sqrt_2=TOL;
    }
    alpha=Fo/sqrt(inside_sqrt_2);
    if(cte==1){ //case when w_0 = w
        theta=PI/2.; //because arctan(infty) = pi/2
    }
    else{
        theta=atan((B*w)/(M*(Wo*Wo-w*w)));
    }
                
        double c1,c2;
    if(fabs(inside_sqrt)<=TOL){ //case when lambda1 and lambda2 are equal, critically damped case
        double lambda;
        lambda=-B/(2.*M);
        c1=-alpha*cos(theta);
        c2=alpha*lambda*cos(theta)-alpha*w*sin(theta)+1.;
        return exp(lambda*t)*(c1+c2*t)+alpha*cos(w*t-theta);
        
        
    }
    else if (inside_sqrt>TOL){ //overdamped case
        double lambda1,lambda2;
        lambda1=(B+sqrt(inside_sqrt))/(2.*M);
        lambda2=(B-sqrt(inside_sqrt))/(2.*M);
        c1=(alpha*lambda2*cos(theta)+alpha*w*sin(theta)-1.)/(lambda1-lambda2);
        c2=(alpha*lambda1*cos(theta)+alpha*w*sin(theta)-1.)/(lambda2-lambda1);

        return c1*exp(-t*lambda1)+c2*exp(-t*lambda2)+alpha*cos(w*t-theta);
        
    }
    else{ //underdamped case
        double real, imag;
        real=-B/(2.*M);
        imag=sqrt(fabs(B*B-4.*M*M*Wo*Wo))/(2.*M);
        c1=-alpha*cos(theta);
        c2=-(alpha*B*cos(theta)+2.*alpha*M*w*sin(theta)-2.*M)/(sqrt(fabs(B*B-4.*M*M*Wo*Wo)));
        return exp(real*t)*(c1*cos(imag*t)+c2*sin(imag*t))+alpha*cos(w*t-theta);
    }
}

void rk_4(double h,int N, double t[N], double x[N], double v[N], double tn){//runge-kutta
    int i=0;
    double k1,k2,k3,k4;
    double l1,l2,l3,l4;
    l1=h*f(x[0],v[0],t[0]);
    k1=h*g(v[0]);
    l2=h*f(x[0]+k1/2.,v[0]+l1/2.,t[0]+h/2.);
    k2=h*g(v[0]+l1/2.);
    l3=h*f(x[0]+k2/2.,v[0]+l2/2.,t[0]+h/2.);
    k3=h*g(v[0]+l2/2.);
    l4=h*f(x[0]+k3,v[0]+l3,t[0]+h);
    k4=h*g(v[0]+k3);
    double aux=0.;
    while(t[i]<=tn+h){
        i++;
        aux++;
        x[i]=x[i-1]+(k1+2*k2+2*k3+k4)/6.;
        v[i]=v[i-1]+(l1+2*l2+2*l3+l4)/6.;
        t[i]=t[0]+aux*h;
        l1=h*f(x[i],v[i],t[i]);
        k1=h*g(v[i]);
        l2=h*f(x[i]+k1/2.,v[i]+l1/2.,t[i]+h/2.);
        k2=h*g(v[i]+l1/2.);
        l3=h*f(x[i]+k2/2.,v[i]+l2/2.,t[i]+h/2.);
        k3=h*g(v[i]+l2/2.);
        l4=h*f(x[i]+k3,v[i]+l3,t[i]+h);
        k4=h*g(v[i]+k3);
    }
}






