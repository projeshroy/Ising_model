#include "includeallglobals.h"
#include "Ising_Input.cpp"
#include "Ising_Energy.cpp"
#include "Ising_Thread_Initialize.cpp"
#include "Ising_MasterThread_Initialize.cpp"
#include "Ising_Output.cpp"

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

   std::string spin_matrix_file_address = directory[1].append(std::string("Ising_Spin_Matrix.dat"));
   std::ofstream spin_matrix_file(spin_matrix_file_address.c_str());
   directory[1] = outdir;
   
   std::string trajectory_file_address = directory[1].append(std::string("Ising_Trajectories.xyz"));
   std::ofstream trajectory_file(trajectory_file_address.c_str());
   directory[1] = outdir;

   std::string energy_file_address = directory[1].append(std::string("Ising_Energy.dat"));
   std::ofstream energy_file(energy_file_address.c_str());
   directory[1] = outdir;

//Equilibration and Run...........................................................................   
   if(RFIM_)
   H_value = RFIM();

   Ising_Atomistic_Energy = Calc_Ising_Energy(H_value, Ising_Spin, Ising_Spin_type);
   Ising_Energy = Ising_Atomistic_Energy.sum();

   for(int step = 0; step <= (equilib_steps + total_steps); step++)
      {
      print_progress_bar<int>(step, (equilib_steps + total_steps), 10);
      int random_spin = getRand(0, total_spins-1);
      Vec_d New_Ising_Spin = Ising_Spin;
      Vec_d New_Ising_Spin_type = Ising_Spin_type;	
      Vec_s New_Ising_Spin_name = Ising_Spin_name;
//....................................................................................
      bool flip = false; bool swap = false;
      double prob = getRand(0, 1);

      if(FLIP_PROBABILITY >= SWAP_PROBABILITY){
      if(prob <= FLIP_PROBABILITY)
      flip = true;
      else swap = true;
      }
      else {
      if(prob <= SWAP_PROBABILITY)
      swap = true;
      else flip = true;
      }

      if(flip){
      New_Ising_Spin_type[random_spin] = flip_type[New_Ising_Spin_type[random_spin]];
      New_Ising_Spin_name[random_spin] = type_names[New_Ising_Spin_type[random_spin]];
      New_Ising_Spin[random_spin] = type_spins[New_Ising_Spin_type[random_spin]];
      }
      else if(swap){
      int swap_spin = random_spin;

      while(swap_spin == random_spin){
      swap_spin = getRand(0, total_spins-1);
      }
      double old_spin = New_Ising_Spin[random_spin];
      double old_spin_type = New_Ising_Spin_type[random_spin];
      std::string old_spin_name = New_Ising_Spin_name[random_spin];
      New_Ising_Spin[random_spin] = New_Ising_Spin[swap_spin];
      New_Ising_Spin_type[random_spin] = New_Ising_Spin_type[swap_spin];
      New_Ising_Spin_name[random_spin] = New_Ising_Spin_name[swap_spin];
      New_Ising_Spin[swap_spin] = old_spin;
      New_Ising_Spin_type[swap_spin] = old_spin_type;
      New_Ising_Spin_name[swap_spin] = old_spin_name;
      }

//.......................................................................................
      bool accept = false;
      Vec_d New_Ising_Atomistic_Energy = Calc_Ising_Energy(H_value, New_Ising_Spin, New_Ising_Spin_type);
      double New_Ising_Energy = New_Ising_Atomistic_Energy.sum();

      if(Ising_Energy >= New_Ising_Energy)
         accept = true;
      else if (exp((Ising_Energy - New_Ising_Energy)/temperature) > getRand(0,1))
         accept = true;
      
      if(accept){
        Ising_Spin = New_Ising_Spin;
        Ising_Spin_type = New_Ising_Spin_type;
        Ising_Spin_name = New_Ising_Spin_name;
        Ising_Energy = New_Ising_Energy;
	Ising_Atomistic_Energy = New_Ising_Atomistic_Energy;
	}

      if((step == equilib_steps)||(((step-equilib_steps) > data_write_step)&&((step-equilib_steps) % data_write_step == 0)))
        {
	for(int spin = 0; spin < total_spins; spin++)
            spin_matrix_file << Ising_Spin[spin] << "  ";
	spin_matrix_file << std::endl;
	Write_Trajectory(trajectory_file, Ising_Spin_name);
        energy_file << Ising_Energy << "   ";
	for(int spin = 0; spin < total_spins; spin++)
	   energy_file << Ising_Atomistic_Energy[spin] << " ";
	energy_file << std::endl;  
	
	}
    
        }
//......................................................................................
}   



