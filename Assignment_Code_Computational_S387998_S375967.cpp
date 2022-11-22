//fn0 = f^n+1 && fn = f^n
#include <iostream>
#include <vector>
#include <cmath>
#define TEMP 79

using namespace std;

vector<double> x;
    double Min = 0.0;
    double Max = 400.0;
    double delta_x = 5.0;
    double delta_t = 0.02;
    double T = 0.5;
    double N = T / delta_t;
    double diff = Max/delta_x;
    double t;
    double c;
    double u = 250;

double calculateAnalyticalValues(vector<double> fn, double i, double t, double max){
	double k = fn.at(i);
	
	if(k >= 0 && k <= 50){
		k = 0.0;
	}
	else if (k > 50 + 250 * t && k <= 110 + 250 * t)
    {
        k = 100 * sin(3.14159265358979323846 * (k - 50 - 250 * t) / 60);
    }
    else if (k > 110 + 250 * t && k <= max)
    {
        k = 0.0;
    }
    else{
    	k = 0.0;
	}
	return k;
}

void AnalyticalSolution(){
	
	vector<double> fn0, fn1, fn;
    fn = x;
    fn.at(0) = 0;
    fn.at(diff) = 0;
    
    fn0 = fn;
    fn1 = fn0;
	
	for(auto j = 0; j<=N; j++){
		t = j * delta_t;
		
		for(auto i=0; i<x.size(); ++i)
        {
            double k = calculateAnalyticalValues(fn, i, t, Max);
            cout<<k<<" ";
        }
        
        cout<<endl;
	}
	
}

//-----------------------------------------------------------------------------------------------------

double calculateInitialBoundaryValues(vector<double> fn, double i, double max){
	double k = fn.at(i);
	
	if(k >= 0 && k <= 50){
		k = 0.0;
	}
	else if (k > 50 && k <= 110)
    {
        k = 100 * sin(3.14159265358979323846 * (k - 50) / 60);
    }
    else if (k > 110 && k <= max)
    {
        k = 0.0;
    }
    else{
    	k = 0.0;
	}
	return k;
}

//-----------------------------------------------------------------------------------------------------

double calculateExplicitUpwindFTBS(vector<double> fn, vector<double> fn0, double i, double max){
	
	return fn.at(i) - ((u * delta_t )/ delta_x ) * (fn.at(i) - fn.at(i-1));
}

void ExplicitUpwindFTBS(){
	
	vector<double> fn0, fn1, fn;
    fn = x;
    fn.at(0) = 0;
    fn.at(diff) = 0;
    
    fn0 = fn;
    fn1 = fn0;
    
	for(auto i=0; i<x.size(); ++i)
        {
            double k = calculateInitialBoundaryValues(fn, i, Max);
            fn0.at(i) = k;
            cout<<" "<<k;
        }
        
        cout<<endl;
		
		fn = fn0;
	for(auto j = 1; j<=N; j++){

		cout<<" "<<fn.at(0);
        for(auto i=1; i<fn0.size()-1; i++)
        {
        	fn.at(0) = 0.0;
        	fn.at(diff) = 0.0;

    		fn0.at(0) = 0.0;
    		fn0.at(diff) = 0.0;
            double k = calculateExplicitUpwindFTBS(fn, fn0, i, Max); // x->fn, fn->x, fn0 = initial, fn = x, fn1 = future;
			if(i==diff){
			k = 0.0;
    		fn0.at(i) = k;
			}
			else{
				fn0.at(i) = k;
			}
			cout<<" "<<k;
        }
    	cout<<" "<<fn.at(diff);
        fn = fn0;
        fn0 = fn1;
        cout<<endl;
	}

}

//-----------------------------------------------------------------------------------------------------

double calculateImplicitUpwindFTBS(vector<double> fn, vector<double> fn0, double i, double max){
	
	return u * delta_t/(delta_x + u * delta_t)* (fn0.at(i - 1)) + delta_x/(delta_x + u * delta_t) * fn.at(i);
}

void ImplicitUpwindFTBS(){
	
	vector<double> fn0, fn1, fn;
    fn = x;
    fn.at(0) = 0;
    fn.at(diff) = 0;
    
    fn0 = fn;
    fn1 = fn0;
    
    
	for(auto i=0; i<x.size(); ++i)
        {
            double k = calculateInitialBoundaryValues(fn, i, Max);
            fn0.at(i) = k;
            cout<<" "<<k;
        }
        
        cout<<endl;
		
		fn = fn0;
	for(auto j = 1; j<=N; j++){

		cout<<" "<<fn.at(0);
        for(auto i=1; i<fn0.size()-1; i++)
        {
        	fn.at(0) = 0.0;
        	fn.at(diff) = 0.0;

    		fn0.at(0) = 0.0;
    		fn0.at(diff) = 0.0;
            double k = calculateImplicitUpwindFTBS(fn, fn0, i, Max); // x->fn, fn->x, fn0 = initial, fn = x, fn1 = future;
			if(i==diff){
			k = 0.0;
    		fn0.at(i) = k;
			}
			else{
				fn0.at(i) = k;
			}
			cout<<" "<<k;
        }
    	cout<<" "<<fn.at(diff);
        fn = fn0;
        fn0 = fn1;
        cout<<endl;
	}

}

//-----------------------------------------------------------------------------------------------------

double calculateLaxWendroff(vector<double> fn, vector<double> fn0, double i, double max){
	
	return ((delta_t*delta_t)*(u*u))/(2*delta_x*delta_x) * (fn.at(i+1)-2*fn.at(i)+fn.at(i-1)) - u*delta_t/(2*delta_x) * (fn.at(i+1)-fn.at(i-1))+fn.at(i);
}

void LaxWendroff(){
	
	vector<double> fn0, fn1, fn;
    fn = x;
    fn.at(0) = 0;
    fn.at(diff) = 0;
    
    fn0 = fn;
    fn1 = fn0;
    
	for(auto i=0; i<x.size(); ++i)
        {
            double k = calculateInitialBoundaryValues(fn, i, Max);
            fn0.at(i) = k;
            cout<<" "<<k;
        }
        
        cout<<endl;
		
		fn = fn0;
		for(auto j = 1; j<=N; j++){

		cout<<" "<<fn.at(0);
        for(auto i=1; i<fn0.size()-1; i++)
        {
        	fn.at(0) = 0.0;
        	fn.at(diff) = 0.0;

    		fn0.at(0) = 0.0;
    		fn0.at(diff) = 0.0;
            double k = calculateLaxWendroff(fn, fn0, i, Max); // x->fn, fn->x, fn0 = initial, fn = x, fn1 = future;
			if(i==diff){
			k = 0.0;
    		fn0.at(i) = k;
			}
			else{
				fn0.at(i) = k;
			}
			cout<<" "<<k;
        }
    	cout<<" "<<fn.at(diff);
        fn = fn0;
        fn0 = fn1;
        cout<<endl;
	}

}

//-----------------------------------------------------------------------------------------------------

void calculateImplicitFTCS(vector<double>& a, vector<double>& b,
	vector<double>& c,vector<double>& d, vector<double>& x) {
    c.at(0) = c.at(0) / b.at(0);
    d.at(0) = d.at(0) / b.at(0);
    for (int i = 1; i < TEMP - 1; i++){
    	c.at(i) = c.at(i) / (b.at(i) - c.at(i-1) * a.at(i));
	}  
    for (int i = 1; i < TEMP; i++){
    	d.at(i) = (d.at(i) - d.at(i-1) * a.at(i)) / (b.at(i) - c.at(i-1) * a.at(i));
	}
    x[TEMP - 1] = d[TEMP - 1];
    for (int i = TEMP - 2; i >= 0; i--) {
        x.at(i) = d.at(i) - c.at(i) * x.at(i+1);
    }

}

void ImplicitFTCS(){
	N = (T + 1e-6) / delta_t;
	diff = (Max + 1e-6) / delta_x;
	
	vector<double> fn0, fn1;
    fn0 = x;
    fn1 = x;
    
    vector<double> z(diff + 1,0);
    vector<vector<double>> fn(N + 1, z);

	vector<double> a (TEMP);
	vector<double> b (TEMP);
	vector<double> c (TEMP);
	vector<double> d (TEMP);
	vector<double> x (TEMP);
	fill (a.begin(),a.begin()+TEMP,(-u * delta_t));
	fill (b.begin(),b.begin()+TEMP,(2 * delta_x));
	fill (c.begin(),c.begin()+TEMP,(u * delta_t));
    
    for(auto i=0; i<fn1.size(); ++i)
        {
            double k = calculateInitialBoundaryValues(fn1, i, Max);
            fn0.at(i) = k;
            fn[0][i] = k;
        }

    for (int i = 1; i <= N; i++)
    {
    	fill (c.begin(),c.begin()+TEMP,(u * delta_t));
    	
        for (int j = 0; j < TEMP; j++)
        {
            d.at(j) = 2 * delta_x * fn[i-1][j+1];
        }
        
        calculateImplicitFTCS(a, b, c, d, x);
        
        for (int j = 1; j <= TEMP; j++)
        {
            fn[i][j] = x.at(j-1);
        }
        
    }
	for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < diff; j++)
        {
            cout << fn[i][j] << " ";
        }
        cout << "\n";
    }
}

//-----------------------------------------------------------------------------------------------------


int main()
{
    int ch1,ch2;
    
	//adding Total no of values into the vector (0,5,10,15,20.....400)
        for(float j = 0; j<=diff; ++j){
			double temp = j;
			temp = temp * delta_x;
			x.push_back(temp);
        }
        
   //end here------------------------------------------
    
    cout<<"\n ------------------- Select Computation Time -------------------";
    cout<<"\n 1: Delta_t = 0.02sec";
    cout<<"\n 2: Delta_t = 0.01sec";
	cout<<"\n 3: Delta_t = 0.005sec";
	cout<<"\n 4: Exit";
	cout<<"\n Enter Your Choice : ";
	cin>>ch1;
    
    switch(ch1)
	{
		case 1:
				delta_t = 0.02;
				N = T / delta_t;
    			cout<<"\n   --------------------- * ---------------------   "<<endl;
    			cout<<"\n 1: Explicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 2: Implicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 3: Lax-Wendroff";
    			cout<<"\n 4: Implicit FTCS (Forward time, Central space)";
    			cout<<"\n 5: Analytical Solution";
    			cout<<"\n 6: Exit";
    			cout<<"\n Enter Your Choice : ";
    			cin>>ch2;
    			
    			switch(ch2)
				{
					case 1:
							cout<<endl;
							ExplicitUpwindFTBS();
							break;
					case 2: 
							cout<<endl;
							ImplicitUpwindFTBS();
							break;
					case 3:
							cout<<endl;
							LaxWendroff();
							break;
					case 4: 
							cout<<endl;
							ImplicitFTCS();
							break;
					case 5:
							cout<<endl;
							AnalyticalSolution();
							break;
					case 6: 
							exit(1);
					default:
							cout<<"\n Wrong Choice!! Please Try Again.. ";
				}
    			
				break;
		case 2: 
    			delta_t = 0.01;
    			N = T / delta_t;
    			cout<<"\n   --------------------- * ---------------------   "<<endl;
    			cout<<"\n 1: Explicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 2: Implicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 3: Lax-Wendroff";
    			cout<<"\n 4: Implicit FTCS (Forward time, Central space)";
    			cout<<"\n 5: Analytical Solution";
    			cout<<"\n 6: Exit";
    			cout<<"\n Enter Your Choice : ";
    			cin>>ch2;
    			
    			switch(ch2)
				{
					case 1:
							cout<<endl;
							ExplicitUpwindFTBS();
							break;
					case 2: 
							cout<<endl;
							ImplicitUpwindFTBS();
							break;
					case 3:
							cout<<endl;
							LaxWendroff();
							break;
					case 4: 
							cout<<endl;
							ImplicitFTCS();
							break;
					case 5:
							cout<<endl;
							AnalyticalSolution();
							break;
					case 6: 
							exit(1);
					default:
							cout<<"\n Wrong Choice!! Please Try Again.. ";
				}
    			
				break;
		case 3:    			
				delta_t = 0.005;
				N = T / delta_t;
    			cout<<"\n   --------------------- * ---------------------   "<<endl;
    			cout<<"\n 1: Explicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 2: Implicit Upwind FTBS (Forward time, Backward space)";
    			cout<<"\n 3: Lax-Wendroff";
    			cout<<"\n 4: Implicit FTCS (Forward time, Central space)";
    			cout<<"\n 5: Analytical Solution";
    			cout<<"\n 6: Exit";
    			cout<<"\n Enter Your Choice : ";
    			cin>>ch2;
    			
    			switch(ch2)
				{
					case 1:
							cout<<endl;
							ExplicitUpwindFTBS();
							break;
					case 2: 
							cout<<endl;
							ImplicitUpwindFTBS();
							break;
					case 3:
							cout<<endl;
							LaxWendroff();
							break;
					case 4: 
							cout<<endl;
							ImplicitFTCS();
							break;
					case 5:
							cout<<endl;
							AnalyticalSolution();
							break;
					case 6: 
							exit(1);
					default:
							cout<<"\n Wrong Choice!! Please Try Again.. ";
				}
    				break;
		case 4:
				exit(1);
		default:
				cout<<"\n Wrong Choice!! Please Try Again.. ";
	}
        
    return 0;
}

