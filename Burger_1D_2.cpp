#include <iostream>
#include<fstream>
using namespace std;

int main(){
    float u0,dt,dx,a,t,dt_choice,c1,c2;
    int n0,nx;
    float C[50][10000];
    int choice,nt;
    cout<<"Solving the Burger Equation"<<endl;
    cout<<"----------------------------"<<endl;
    cout<<"C,t + u0 C,x = a C,xx"<<endl;
    cout<<"----------------------------"<<endl<<endl;
    cout<<"Following cases will be solved"<<endl<<endl;
    cout<<"Case 1: u0 = 1;   y = 0.01"<<endl;
    cout<<"Case 2: u0 = 1;   y = 0.001"<<endl;
    cout<<"Case 3: u0 = 1;   y = 0.1"<<endl<<endl;
    cout<<"Enter the total time period (eg. max 5 sec)...";
    cin>>t; // time for the full case taken as 5 secs
    cout<<"Enter the choice of case to be solved...";
    cin>>choice;
    cout<<endl;
    switch (choice)
    {
        case 1:
                cout<<"Case 1: u0 = 1;  y = 0.01"<<endl;
                u0 = 1;
                a = 0.01;
                break;
        case 2:
                cout<<"Case 2: u0 = 1;  y = 0.001"<<endl;
                u0 = 1;
                a = 0.001;
                break;
        case 3:
                cout<<"Case 3: u0 = 1;  y = 0.1"<<endl;
                u0 = 1;
                a = 0.1;
                break;
    }
    dx = 0.5;
    u0 = 1;
    dt_choice = min((2*a)/(u0*u0),(dx*dx)/(2*a));
    cout<<"dt should be less than or equal to "<<dt_choice<<" sec for stable solution"<<endl;
    cout<<"Enter dt value...";
    cin>>dt;
    nt = t/dt;
    cout<<"Number of time steps are "<<nt<<endl;
   
    c2 = (u0*dt)/(2*dx); // u0/2dx
    c1 = (a*dt)/(dx*dx); // a/dx^2
    cout<<"c2 = "<<c2;
    cout<<"c1 = "<<c1;
    nx = 42*dx;
    n0 = (nx/2);
    for (int i=0;i<=nx;i++)
        {
            for (int j=0;j<=nt;j++)
            {
                C[i][j] = 0; // initial condition
            }
        }
    /*for (int i=n0;i<=nx;i++)
        {
            for (int j=0;j<=nt;j++)
            {
                C[nx][j] = C[nx-2][j]; // right BC
                C[i][j+1] = (1-2*c2)* C[i][j] + C[i-1][j] *(c2+c1) + (c2-c1)*C[i+1][j];
                // C[0][j] = C[2][j]; // left BC
            }
        }*/
    for (int j=0;j<=nt;j++)
        {
            for (int i=1;i<=nx;i++)
            {
                C[n0][j] = 1; // centre BC
                C[2][j+1] = (1-2*c2)* C[2][j] + C[1][j] *(c2+c1) + (c2-c1)*C[3][j];
                C[0][j+1] = C[2][j+1]; // left BC
                C[nx-2][j+1] = (1-2*c2)* C[nx-2][j] + C[nx-1][j] *(c2+c1) + (c2-c1)*C[nx-3][j];
                C[nx][j+1] = C[nx-2][j+1]; // left BC
                C[i][j+1] = (1-2*c2)* C[i][j] + C[i-1][j] *(c2+c1) + (c2-c1)*C[i+1][j];
                // C[0][j] = C[2][j]; // left BC
            }
        }
    
     /*for (int i=n0;i>=0;i--)
        {
            for (int j=0;j<=nt;j++)
            {
                C[n0][j] = 1; // centre BC
                C[i][j+1] = C[i][j] - c1*(C[i+1][j]-C[i-1][j]) + c2*(C[i+1][j]-2*C[i][j]+C[i-1][j]);
                C[0][j] = C[2][j]; // left BC
                C[nx][j] = C[nx-2][j]; // right BC
            }
        }
    */
    // Writing the file
                  
    ofstream Myfile("/Users/saurabh/Desktop/A.csv");
    for(int i=0;i<nx;i++)
        {
            for(int j=0;j<=nt;j++)
            {
                Myfile << C[i][j]<<" ";
            }
        Myfile << endl;
        }
        Myfile.close();
    cout<<endl<<endl<<endl;
    return 0;
}

