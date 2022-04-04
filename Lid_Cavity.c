#include<stdio.h>
#include<math.h>
#include<stdlib.h>

void main(){
    //Inputs
    float X = 1;          		//Size of cavity X
    float Y = 1;          		//Size of cavity Y
    int nx = 101;         	//# gridpoints of X
    int ny = 101;        		//# gridpoints of Y
    float dx = X/(nx-1);  	//grid size along X
    float dy = Y/(ny-1);  	//grid size along Y
    //General data
    float Re = 500.0;     	//Assumed Reynolds number
    float U0 = 1.0;       	//Top lid velocity
    float dt = 0.001;     	//Time step
    
//Initialising necessary functions
    float psi[ny][nx];
    float psiNext[ny][nx];
    float zeta[ny][nx];
    float zetaNext[ny][nx];
    float u[ny][nx];
    float v[ny][nx];
    float p[ny][nx];
    //Initial condition
    for(int i=0;i<ny;i++){
        for(int j=0;j<nx;j++){
            psi[i][j] = 0;
            zeta[i][j] = 0;
            u[i][j] = 0;
            v[i][j] = 0;
            p[i][j] = 0;
            zetaNext[i][j] = 0;
            psiNext[i][j] = 0;
        }
    }
    float D1,D2;        //Dummy variables
    float EZ,EP;        //Zeta and Psi Error
    int limit = 100000;
    float EZI[1][limit],EPI[1][limit],Index[1][limit];
    for(int T=0; T<limit; T++){
        EZ = 0;
        EP = 0;
        //Boundary conditions
        //Stream function at the walls = 0 {Included in the initial condition}
        //Vorticity transport equation on the walls
        for(int j=1;j<nx-1;j++){ 
            zeta[0][j] = 2*(psi[0][j]-U0*dy-psi[1][j])/(dy*dy);              	//Top lid
            zeta[ny-1][j] = -2*(psi[ny-2][j]-psi[ny-1][j])/(dy*dy);          	//Bottom wall
            u[0][j] = U0;
        }
        for(int i=1;i<ny-1;i++){
            zeta[i][0] = -2*(psi[i][1]-psi[i][0])/(dx*dx);                   		//Left wall
            zeta[i][nx-1] = -2*(psi[i][nx-2]-psi[i][nx-1])/(dx*dx);          	//Right wall
        }
        //Vorticity Calculation
        for(int i=1;i<ny-1;i++){
            for(int j=1;j<nx-1;j++){
                D1=(u[i][j+1]*zeta[i][j+1]-u[i][j-1]*zeta[i][j-1])/(2*dx)
                    +(v[i-1][j]*zeta[i-1][j]-v[i+1][j]*zeta[i+1][j])/(2*dy);
                D2=((zeta[i][j+1]-2*zeta[i][j]+zeta[i][j-1])/(dx*dx))/Re
                    +((zeta[i-1][j]-2*zeta[i][j]+zeta[i+1][j])/(dy*dy))/Re;
                zetaNext[i][j]=zeta[i][j]+dt*(D2-D1);
                EZ=EZ+abs(zetaNext[i][j]-zeta[i][j]);
            }
        }
        //Stream Function
        for(int i=1; i<ny; i++){
            for(int j=1; j<nx; j++){
                psiNext[i][j]=(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]+zetaNext[i][j]*dx*dx)/4;
                EP=EP+abs(psiNext[i][j]-psi[i][j]);
            }
            //Boundary conditions updated
            for(int j=0;j<nx;j++){
                zetaNext[0][j] = 2*(psi[0][j]-U0*dy-psi[1][j])/(dy*dy);              //Top lid
                zetaNext[ny-1][j] = -2*(psi[ny-2][j]-psi[ny-1][j])/(dy*dy);         //Bottom wall
                EZ=EZ+abs(zetaNext[0][j]-zeta[0][j]);
                EZ=EZ+abs(zetaNext[ny-1][j]-zeta[ny-1][j]);
            }
            for(int i=0;i<ny;i++){
                zetaNext[i][0] = -2*(psi[i][1]-psi[i][0])/(dx*dx);                   	//Left wall
                zetaNext[i][nx-1] = -2*(psi[i][nx-2]-psi[i][nx-1])/(dx*dx);       //Right wall
                EZ=EZ+abs(zetaNext[i][0]-zeta[i][0]);
                EZ=EZ+abs(zetaNext[i][nx-1]-zeta[i][nx-1]);
            }
        }
        //Calculating U and V
        for(int i=1;i<ny-1;i++){
            for(int j=1;j<nx-1;j++){
                u[i][j] = (psi[i][j]-psi[i+1][j])/(dy);
                v[i][j] = (psi[i][j-1]-psi[i][j])/(dx);
            }
        }
        for(int i=0;i<ny;i++){
            for(int j=0;j<nx;j++){
                zeta[i][j] = zetaNext[i][j];
            }
        }
        EPI[0][T]=EP/(ny*nx);
        EZI[0][T]=EZ/(ny*nx);
        Index[0][T]=T;
        if (T>=10){
            if (EPI[0][T]<=0.00000001 &&EZI[0][T]<=0.00000001){
                printf("Convergence of Vorticity transport occured at %d iteration\n",T);
                break;
            }
        }
    }
    //Calculating Pressure
    p[ny][0]=1;                         	//Reference pressure
    float EPr;                          	//Pressure error
    float D3,D4;                        	//Extra dummy variables
    for(int T=0; T<limit; T++){
        EPr=0;
        for(int i=ny-2;i>=0;i--){       	//Left Wall
            D1=(-3*zeta[i][0]+4*zeta[i][1]-zeta[i][2])/(2*Re);
            D2=p[i][0];
            p[i][0]=p[i+1][0]-D1;
            EPr=EPr+abs(p[i][0]-D2);
        }
        for(int i=ny-2;i>=0;i--){       	//Right Wall
            D1=(-3*zeta[i][nx-1]+4*zeta[i][nx-2]-zeta[i][nx-3])/(2*Re);
            D2=p[i][nx-1];
            p[i][nx-1]=p[i+1][nx-1]-D1;
            EPr=EPr+abs(p[i][0]-D2);
        }
        for(int j=1;j<nx;j--){       	//Bottom Wall
            D1=(-3*zeta[ny-1][j]+4*zeta[ny-2][j]-zeta[ny-3][j])/(2*Re);
            D2=p[ny-1][j];
            p[ny-1][j]=p[ny-1][j-1]-D1;
            EPr=EPr+abs(p[ny-1][j]-D2);
        }
        //Pressure Poisson
        for(int i=1;i<ny;i++){
            for(int j=1;j<nx;j++){
                D1=(psi[i][j+1]-2*psi[i][j]+psi[i][j-1])*(psi[i-1][j]-2*psi[i][j]+psi[i+1][j]);
                D2=(psi[i-1][j+1]-psi[i+1][j+1]-psi[i-1][j-1])*(psi[i-1][j+1]-psi[i+1][j+1]+psi[i+1][j-1]);
                D3=D1/(dx*dx*dy*dy)-D2/(4*dx*dy*4*dx*dy);
                D4=p[i][j];
                p[i][j]=(p[i][j+1]+p[i][j-1]+p[i+1][j]+p[i-1][j]-dx*dy*(2*D2))/4;
                EPr=EPr+abs(p[i][j]-D4);
            }
        }
        if (T>=10){
            if (EPr/(ny*nx)<=0.00000001){
                printf("Convergence of Pressure occured at %d iteration\n",T);
                break;
            }
        }
    }
}
