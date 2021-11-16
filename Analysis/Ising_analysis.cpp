#include "../includeallglobals.h"
#include "../Ising_Input.cpp"
#include "../Ising_Energy.cpp"
#include "../Ising_Thread_Initialize.cpp"
#include "../Ising_MasterThread_Initialize.cpp"


int main (int argc, char** argv)
{
   std::string directory[2];
   directory[0]= argv[1];
   directory[1]= argv[2];
   std::string outdir = directory[1];
   std::string indir = directory[0];

   std::string input_file_address = directory[0].append(std::string("ising.in"));
   directory[0] = indir;
   Read_Input(input_file_address);
   
   int max_threads = 0;
#ifdef _OPENMP
   max_threads = omp_get_max_threads();
#endif   
   MasterThread_Initialize();
   Thread_Initialize();
#pragma omp parallel num_threads(max_threads)
 {
   int thread_id = omp_get_thread_num();
   if(thread_id > 0){ 
   Read_Input(input_file_address);
   MasterThread_Initialize();
   Thread_Initialize();
 }
 }
   
   std::string energy_input_file_address = directory[0].append(std::string("ising_energy.in"));
   std::ifstream energy_input_file(energy_input_file_address.c_str());
   directory[0] = indir;
   
   if (!energy_input_file.is_open())
      std::cout << " Unable to find file " << std::endl;

   std::string energy_hist_file_address = directory[1].append(std::string("Ising_Energy_hist.dat"));
   std::ofstream energy_hist_file(energy_hist_file_address.c_str());
   energy_hist_file.precision(14);
   directory[1] = outdir;

   std::string atomistic_energy_hist_file_address = directory[1].append(std::string("Ising_Atomistic_Energy_hist.dat"));
   std::ofstream atomistic_energy_hist_file(atomistic_energy_hist_file_address.c_str());
   energy_hist_file.precision(14);
   directory[1] = outdir;

   std::cout.unsetf ( std::ios::floatfield );// floatfield not set
   std::cout.precision(10);
   std::cout.setf( std::ios::fixed, std:: ios::floatfield );
//...........................................................................................
   std::string read_string;
   std::string spin_file_address, energy_file_address;
   energy_input_file >> read_string;
   energy_input_file >> read_string >> energy_file_address;
   energy_input_file >> read_string >> spin_file_address;
   
   std::ifstream energy_file(energy_file_address.c_str());
   if (!energy_file.is_open()) 
      std::cout << " Unable to find file " << std::endl;
   std::ifstream spin_file(spin_file_address.c_str());
   if (!spin_file.is_open()) 
      std::cout << " Unable to find file " << std::endl;

   double min_E, max_E, bin_width;
   double atomistic_min_E, atomistic_max_E, atomistic_bin_width;
   energy_input_file >> read_string >> min_E;
   energy_input_file >> read_string >> max_E;
   energy_input_file >> read_string >> bin_width;
   energy_input_file >> read_string >> atomistic_min_E;
   energy_input_file >> read_string >> atomistic_max_E;
   energy_input_file >> read_string >> atomistic_bin_width;

   Vec_d Energy_bins, Energy_hist;
   int Energy_bins_size = (int)((max_E - min_E)/bin_width) + 1;
   Energy_bins.resize(Energy_bins_size);
   for(int i = 0; i < Energy_bins.size(); i++)
       Energy_bins[i] = min_E + i * bin_width;
   
   Energy_hist.resize(Energy_bins_size);
   Energy_hist.setZero();

   Vec_d atomistic_Energy_bins, atomistic_Energy_hist;
   int atomistic_Energy_bins_size = (int)((atomistic_max_E - atomistic_min_E)/atomistic_bin_width) + 1;
   atomistic_Energy_bins.resize(atomistic_Energy_bins_size);
   for(int i = 0; i < atomistic_Energy_bins.size(); i++)
       atomistic_Energy_bins[i] = atomistic_min_E + i * atomistic_bin_width;
   
   atomistic_Energy_hist.resize(atomistic_Energy_bins_size);
   atomistic_Energy_hist.setZero();

   Mat_d typewise_atomistic_Energy_hist;
   typewise_atomistic_Energy_hist.resize(total_types, atomistic_Energy_bins_size);
   typewise_atomistic_Energy_hist.setZero();
//...............................................................................................................   
std::cout << " total_steps_count " << total_steps_count << std::endl;
   for(int step = 0; step < total_steps_count; step++){
      print_progress_bar<int>(step, total_steps_count, 20);
      energy_file >> Ising_Energy;
      for(int spin = 0; spin < total_spins; spin++){
          spin_file >> Ising_Spin[spin];  
          energy_file >> Ising_Atomistic_Energy[spin];

	  for(int t = 0; t < total_types; t++){
		if(Ising_Spin[spin] == type_spins[t]){
	  	   Ising_Spin_type[spin] = t;	
		   Ising_Spin_name[spin] = type_names[t];
		   break;
		}
 	     }
	  }
      
      for(int h = 1; h < Energy_bins.size(); h++){
         if( (Ising_Energy < Energy_bins[h])&&(Ising_Energy >= Energy_bins[h-1]) ){
	     Energy_hist[h-1]++;
	     break;
	     }
	   }
      
      for(int h = 1; h < atomistic_Energy_bins.size(); h++){
	  for(int spin = 0; spin < total_spins; spin++){ 
	  int type = Ising_Spin_type[spin];
          if( (Ising_Atomistic_Energy[spin] < atomistic_Energy_bins[h])
	  &&(Ising_Atomistic_Energy[spin] >= atomistic_Energy_bins[h-1]) ){
	      atomistic_Energy_hist[h-1]++;
	      typewise_atomistic_Energy_hist(type, h-1)++;
	      typewise_atomistic_Energy_hist(flip_type[type], h-1)++;
	      break;
	      }
	    }
	  }
       }

    energy_hist_file << "# Energy    Probability" << std::endl;
    for(int i = 0; i < Energy_bins.size(); i++)
       energy_hist_file << Energy_bins[i] << "   " << Energy_hist[i]/Energy_hist.sum()<< std::endl;

    atomistic_energy_hist_file << "# Energy   Probability   typewise_probability(X" << total_types << ")" << std::endl;
    for(int i = 0; i < atomistic_Energy_bins.size(); i++){
        atomistic_energy_hist_file << atomistic_Energy_bins[i] 
	 		           << "   " << atomistic_Energy_hist[i]/atomistic_Energy_hist.sum();
	for(int j = 0; j < total_types; j++)
	    atomistic_energy_hist_file << "   " << typewise_atomistic_Energy_hist(j, i)/typewise_atomistic_Energy_hist.row(j).sum();

	atomistic_energy_hist_file << std::endl;
	}
//..................................................................................................................
}
