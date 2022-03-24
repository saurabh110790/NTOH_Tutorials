#include <iostream>
#include<fstream>
using namespace std;

int main(){
    float u0,dt,dx,a,t,dt_choice,c1,c2;
    int n0,nx;
    float C[50][5000];
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
    dt_choice = min((2*a)/(u0*u0),(dx*dx)/(2*a));
    cout<<"dt should be less than or equal to "<<dt_choice<<" sec for stable solution"<<endl;
    cout<<"Enter dt value...";
    cin>>dt;
    nt = t/dt;
    cout<<"Number of time steps are "<<nt<<endl;
   
    c1 = (u0*dt)/(2*dx); // u0/2dx
    c2 = (a*dt)/(dx*dx); // a/dx^2
    nx = 42*dx;
    n0 = (nx/2)-1;
    for (int i=0;i<=nx;i++)
        {
            for (int j=0;j<=nt;j++)
            {
                C[i][j] = 0; // initial condition
                C[n0][0] = 1; // centre BC
                C[0][j] = C[2][j]; // left BC
                C[nx][j] = C[nx-2][j]; // right BC
            }
        }
    for (int i=0;i<=nx;i++)
        {
            for (int j=0;j<=nt;j++)
            {
                C[i][j+1] = C[i][j] - c1*(C[i+1][j]-C[i-1][j]) + c2*(C[i+1][j]-2*C[i][j]+C[i-1][j]);
                C[0][j] = C[2][j]; // left BC
                C[nx][j] = C[nx-2][j]; // right BC
            }
        }
    
    // Writing the file
                  
    ofstream Myfile("/Users/saurabh/Desktop/A.csv");
    for(int i=0;i<=nx;i++)
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

