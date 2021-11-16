
Mat_d RFIM()
{
    Mat_d New_H_value; 
    New_H_value.resize(total_spins, total_spins);
    New_H_value.setZero();

    //RFIM......................................
    #pragma omp parallel for schedule(static) shared (New_H_value)
      for(int spin_i = 0; spin_i < total_spins-1; spin_i++){
          for(int spin_j = spin_i+1; spin_j < total_spins; spin_j++){
              New_H_value(spin_i, spin_j) = Field;
	      New_H_value(spin_j, spin_i) = Field;

              if(Sigma_H != 0){
                std::normal_distribution<double> distribution(Field, Sigma_H);
		bool found = false;
		while(!found){
                     double H = distribution(generator);
		     if( ( std::abs(H) > 0) && ( std::abs(H) < std::abs(2*Field) ) ){
	 	       New_H_value(spin_i, spin_j) = H;
		       New_H_value(spin_j, spin_i) = H;
		       found = true;
		       }
		     }
  	         }
	     }
	 }
} 
