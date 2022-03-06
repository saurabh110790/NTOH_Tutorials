#include <iostream>
#include<fstream>
using namespace std;

int main() {
    int nx, ny, niter,i,j;
    float x,y, dx,dy;
    float b[100][100],Tn[100][100],T[100][100],e[100];
    cout<<"Prgramme for Solving Poisson's Equation Using 5 Point Formula"<<endl;
    cout<<"--------------------------------------------------------------------"<<endl;
    cout<<"                         "<<"Given Datasets"<<"              "<<endl;
    cout<<"--------------------------------------------------------------------"<<endl<<endl;
    cout<<"Case 1:"<<" Drichlet Boundary Condition with source and sink"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"All values are temperatures (Boundary Value Problem)"<<endl<<endl;
    cout<<"        40    "<<endl;
    cout<<"   ------------"<<endl;
    cout<<"  |"<<"            "<<"|"<<endl;
    cout<<"0 "<<"|""  .A    .B  "<<"|"<<" 0     "<<"Drichlet Boundary Condition"<<endl;
    cout<<"  |"<<"            "<<"|"<<endl;
    cout<<"   ------------"<<endl;
    cout<<"        40    "<<endl<<endl;
    cout<<"A(x,y) = (x/5,y/2)"<<";  "<<"Qa = -1.0"<<"   Sink"<<endl;
    cout<<"B(x,y) = (4x/5,y/2)"<<"; "<<"Qb =  1.0"<<"   Source"<<endl;
    cout<<endl<<endl;
    
    cout<<"Case 2:"<<" Neumann Boundary Condition with source and sink"<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<"All values are temperatures (Boundary Value Problem)"<<endl<<endl;
    cout<<"      T,x = 10    "<<endl;
    cout<<"   ------------"<<endl;
    cout<<"  |"<<"            "<<"|"<<endl;
    cout<<"0 "<<"|""  .A    .B  "<<"|"<<" T,y = 1   "<<"Neumann Boundary Condition"<<endl;
    cout<<"  |"<<"            "<<"|"<<endl;
    cout<<"   ------------"<<endl;
    cout<<"        40    "<<endl<<endl;
    cout<<"A(x,y) = (x/5,y/2)"<<endl;
    cout<<"B(x,y) = (4x/5,y/2)"<<endl;
    cout<<"(Qa,Qb) = (1,1) or (Qa,Qb) = (-1,-1)"<<endl<<endl;
    cout<<"--------------------------------------------------------------------"<<endl<<endl;
    
    cout<<"Assuming that the grids in x and y direction are not equally spaced"<<endl<<endl;
    cout<<"Five Point Formula for CDM is given as:"<<endl<<endl;
    cout<<"Ti,j = dy^2(Ti+1,j + Ti-1,j) + dx^2(Ti,j+i + Ti,j-1)"<<endl;
    cout<<"       ___________________________________"<<endl;
    cout<<"                    2(dx^2+dy^2)          "<<endl;
    
    cout<<"Solving Case 1 "<<endl<<endl;
    
    cout<<"Enter the distance between (0,0) and (n,0) coordinate:  ";
    cin>>x;
    cout<<"Enter the distance between (0,0) and (0,m) coordinate:  ";
    cin>>y;
    cout<<"Enter the number of steps in X direction (max 100):  ";
    cin>>nx;
    cout<<"Enter the number of steps in Y direction (max 100):  ";
    cin>>ny;
    dx = x/(nx-1);
    dy = y/(ny-1);
    cout<<"Width of space step on x is "<<dx<<endl;
    cout<<"Width of space step on y is "<<dy<<endl;
    cout<<"Number of iteration steps:  ";
    cin>>niter;
    for (i=0;i<nx;i++)
        {
            for (j=0;j<ny;j++)
                {
                    b[i][j] = 0; // Preallocating Q values
                    Tn[i][j] = 0; // Preallocating initial values of Tn(i,j) as 0
                    T[i][j] = 0; // initially allocating values of T(i,j) as 0
                }
        }
    // Boundary Conditions
    for (i=0;i<nx;i++)
        {
            for (j=0;j<ny;j++)
            {
                T[0][j] = 0;
                T[nx-1][j] = 0;
                T[i][0] = 40;
                T[i][ny-1] = 40;
            }
        }
    // Allocating Source and Sink Values
    b[nx/5][ny/2] = -1;
    b[4*nx/5][ny/2] = 1;
    
    // Loop for calculating the 5 Point Formula
    for(int k=1;k<niter;k++)
    {
        for (i=0;i<nx;i++)
        {
            for (j=0;j<ny;j++)
            {
                Tn[i][j] = T[i][j];
            }
        }
        for (i=0;i<nx;i++)
        {
            for (j=0;j<ny;j++)
            {
                T[i][j] = (((dy*dy) * (Tn[i+1][j] + Tn[i-1][j])) + ((dx*dx) * (Tn[i][j+1] + Tn[i][j-1])) - (b[i][j]*dx*dx*dy*dy))/(2*((dx*dx)+(dy*dy))) ;
                T[0][j] = 0;
                T[nx-1][j] = 0;
                T[i][0] = 40;
                T[i][ny-1] = 40;
            }
        }
        e[k] = T[nx/2][ny/2];
    }
    for (j=0;j<ny;j++)
        {
            for (i=0;i<nx;i++)
            {
                cout<<T[i][j]<<"    ";
            }
            cout<<endl;
        }
    /*
     for (j=1;j<niter;j++)
        {
            cout<<e[j]<<"   ";
        }
    cout<<endl<<endl;
    */
     // Writing the file
               
                   ofstream Myfile("/Users/saurabh/Desktop/Drichlet.csv");
                   for(i=0;i<20;i++)
                   {
                       for(j=0;j<21;j++)
                       {
                           Myfile << T[i][j]<<" ";
                       }
                       Myfile << endl;
                   }
                   Myfile.close();

    return 0;
}
