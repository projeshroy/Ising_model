
double RBIM(int& type, int& nei_type) 
{
    //RBIM......................................................
    double J = J_matrix(type, nei_type);
/*	      
       std::normal_distribution<double> generator(J_matrix(type, nei_type), Delta_J(type, nei_type));
       dool found = false;
       while(!found){
             rand = generator(gen);
	     if( (rand > 0) && ( rand < (2*J_matrix(nei_type, type))) ){
	       J = rand;
	       found = true;
	       }
	     }
*/
         double rand = getRand(0, 1);
         if(rand <= RBIM_p(type, nei_type))
            J += Delta_J(type, nei_type);
  	 else J -= Delta_J(type, nei_type);

       return J;
} 
