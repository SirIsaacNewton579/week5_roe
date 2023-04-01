#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;
double g=1.4;
void mul_matrix(double *A,double *B,double *ret,int m=3,int o=3,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = 0.0;
            for(int k=0;k<o;k++){
                ret[n*i+j] += A[i*o+k]*B[k*n+j];
            }
        }
    }
}
void mul_matrix(double *A,double a,int m=1,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            A[n*i+j] *= a;
        }
    }
}
void add_matrix(double *A,double *B,double *ret,int m=1,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = A[n*i+j]+B[n*i+j];
        }
    }
}
void minus_matrix(double *A,double *B,double *ret,int m=1,int n=3){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ret[n*i+j] = A[n*i+j]-B[n*i+j];
        }
    }
}
void Su(double *U,double *S){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    S[0] = 0.5*u*u - c*c/(g-1);
    S[1] = -u;
    S[2] = 1.;
    S[3] = -u -(g-1)/c*0.5*u*u;
    S[4] = 1+(g-1)/c*u;
    S[5] = -(g-1)/c;
    S[6] = -u+(g-1)/c*0.5*u*u;
    S[7] = 1-(g-1)/c*u;
    S[8] = (g-1)/c;
}
void invSu(double *U,double *invS){
    double rho=U[0],u=U[1]/U[0],p=(g-1)*(U[2]-0.5*rho*u*u);
    double c = sqrt(g*p/rho);
    double h = 0.5*u*u + g/(g-1)*p/rho;
    invS[0] = -(g-1)/(c*c);
    invS[1] = -1./(2*c);
    invS[2] = 1./(2*c);
    invS[3] = -(g-1)/(c*c)*u;
    invS[4] = -(u-c)/(2*c);
    invS[5] = (u+c)/(2*c);
    invS[6] = -(g-1)/(c*c)*0.5*u*u;
    invS[7] = -1./(2*c)*(h-u*c);
    invS[8] = 1./(2*c)*(h+u*c);
}
void fU(double *U,double *f){
    f[0] = U[1];
    f[1]=(g-1)*U[2]+0.5*(3-g)*U[1]*U[1]/U[0];
    f[2]=g*U[2]*U[1]/U[0]+0.5*(g-1)*U[1]*U[1]*U[1]/(U[0]*U[0]);
}
void Ubar(double *UL,double *UR,double *Ub){
    double rhoL=UL[0],uL=UL[1]/UL[0],pL=(g-1)*(UL[2]-0.5*rhoL*uL*uL);
    double hL = 0.5*uL*uL + g/(g-1)*pL/rhoL;

    double rhoR=UR[0],uR=UR[1]/UR[0],pR=(g-1)*(UR[2]-0.5*rhoR*uR*uR);
    double hR = 0.5*uR*uR + g/(g-1)*pR/rhoR;

    double rhob = pow((sqrt(rhoL)+sqrt(rhoR))/2,2);
    double ub = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(2*sqrt(rhob));
    double hb = (sqrt(rhoL)*hL+sqrt(rhoR)*hR)/(2*sqrt(rhob));
    double pb = (g-1)/g*rhob*(hb - 0.5*ub*ub);
    Ub[0] = rhob;Ub[1] = rhob*ub; Ub[2] = rhob*hb-pb;
}
void LBD_abs(double *U,double* LBDa){
    double rho,u,p,c;
    rho = U[0];u = U[1]/U[0];p =(g-1)*(U[2]-0.5*rho*u*u);
    c = sqrt(g*p/rho);
    LBDa[0] = abs(u);LBDa[4] = abs(u-c);LBDa[8] = abs(u+c);
}
void Flux_fds(double *U,int Nx,double *fo){
    double UL[3],UR[3],Ub[3],fUL[3],fUR[3];
    double S[9],invS[9],LBDa[9] = {0.0};
    double Ab[9];
    double tmp[9],tmp2[3],tmp3[3];
    int i,j;
    for(j=0;j<Nx-1;j++){
        for(i=0;i<3;i++) {UL[i]=U[i*Nx+j];UR[i]=U[i*Nx+j+1];}
        fU(UL,fUL);
        fU(UR,fUR);
        Ubar(UL,UR,Ub);
        Su(Ub,S);
        invSu(Ub,invS);
        LBD_abs(Ub,LBDa);
        mul_matrix(invS,LBDa,tmp); // tmp = S^-1*Lambda 
        mul_matrix(tmp,S,Ab);  //|A| = S^-1*Lambda*S
        minus_matrix(UR,UL,tmp2);  //tmp2 = UR-UL
        mul_matrix(Ab,tmp2,tmp3,3,3,1);  //tmp3 = |A|*(UR-UL)
        add_matrix(fUL,fUR,tmp2); //tmp2 = f(UR)+f(UL)
        minus_matrix(tmp2,tmp3,tmp2); //tmp2 = f(UR)+f(UL)-|A|*(UR-UL)
        mul_matrix(tmp2,0.5);
        for(i=0;i<3;i++) fo[i*(Nx-1)+j] = tmp2[i];
    }
}
void updateU(double *U,int Nx,double dx,double dt,int Nt){
    double fo[3][Nx-1];  
    double U1[3][Nx],U2[3][Nx],Unext[3][Nx];
    double cfl = dt/dx;
    int i,j,N;
    for(N=1;N<=Nt;N++){
        Flux_fds(U,Nx,fo[0]);
        for(i=0;i<3;i++) {U1[i][0] = U[i*Nx];U1[i][Nx-1] = U[i*Nx+Nx-1];}
        for(j=1;j<Nx-1;j++){
            for(i=0;i<3;i++){
                U1[i][j] = U[j+i*Nx] - cfl*(fo[i][j]-fo[i][j-1]);
            }    
        }
        Flux_fds(U1[0],Nx,fo[0]);
        for(i=0;i<3;i++) {U2[i][0] = U1[i][0];U2[i][Nx-1] = U1[i][Nx-1];}
        for(j=1;j<Nx-1;j++){
            for(i=0;i<3;i++){
                U2[i][j] = 0.75*U[j+i*Nx] + 0.25*(U1[i][j]- cfl*(fo[i][j]-fo[i][j-1]));
            }
        }
        Flux_fds(U2[0],Nx,fo[0]);
        for(i=0;i<3;i++) {Unext[i][0] = U2[i][0];Unext[i][Nx-1] = U2[i][Nx-1];}
        for(j=1;j<Nx-1;j++){
            for(i=0;i<3;i++){
                Unext[i][j] = 1.0*U[j+i*Nx]/3.0 + 2.0/3.0*(U2[i][j]- cfl*(fo[i][j]-fo[i][j-1]));
            }
        }
        for(i=0;i<3;i++){
            for(j=0;j<Nx;j++){
                U[j+i*Nx] = Unext[i][j];
            }
        }
    }
}
int main(){
    //开始计算
    int Nx = 101;
    double dx = 1./(Nx-1),dt = 0.001;
    double t_end = 0.14,Nt = round(t_end/dt);
    cout << "Nt=" << Nt << endl;
    double U[3][Nx];
    
    for(int i=0;i<Nx;i++){
        //cout << i*dx << endl;
        U[0][i] = (i>Nx/2 ? 0.125 : 1);  //rho
        U[1][i] = 0.;  // rho*u
        U[2][i] = (i>Nx/2 ? 0.1 : 1)/(g-1); //E = 1/2*rho*u^2 + p/(g-1)
    }
    updateU(U[0],Nx,dx,dt,Nt);  
    //输出
    ofstream csvfile;
    csvfile.open("FDS-Roe_t=0.14.csv", ios::out | ios::trunc);
    csvfile <<"x" << "," << "rho"<<","<< "u" <<","<<"p"<< endl;
    for(int i=0;i<Nx;i++){
        csvfile <<i*dx <<","<< U[0][i]<<","<<U[1][i]/U[0][i]<<","<<(g-1)*(U[2][i]-0.5*U[1][i]*U[1][i]/U[0][i]) << endl;
    }
    csvfile.close();
    //system("pause");
    return 0;
}