#pragma once
#include "./include/Eigen_3.3.1/Eigen/Core"
#include "./include/Eigen_3.3.1/Eigen/Dense"

template <class T>
class EigenMultiDimArray
{
      public:
            std::vector<int> arraydims;
	    std::vector<int> arraydims_offset;
	    int total_size = 1;
	    int D = 1;
            Eigen::Matrix<T, Eigen::Dynamic, 1> DummyMD;

	    int getD(){
   	    return D;
	    }
            
	    int getMDSize()
            {
	     total_size = 1;
	     for(int i = 0; i < D; i++)
	        total_size *= arraydims[i];
             return total_size;
            }

	    Eigen::Matrix<T, Eigen::Dynamic, 1> getarraydims(){
	     Eigen::Matrix<T, Eigen::Dynamic, 1> new_arraydims; new_arraydims.resize(D);
	     for(int d = 0; d < D; d++)
		new_arraydims[d] = arraydims[d];
             return new_arraydims;
	    }            

            void resize_DummyMD(std::vector<int> idarray)
            {
	       arraydims = idarray;
	       arraydims_offset = idarray;
	       D = idarray.size();

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);
            }

	    void resize_DummyMD(std::vector<int> idarray, T value)
            {
	       arraydims = idarray;
	       arraydims_offset = idarray;
	       D = idarray.size();

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);

               #pragma omp parallel for
	       for(int i = 0; i < total_size; i++)
	           DummyMD[i] = value;
            }

            void resize_DummyMD(Eigen::Matrix<T, Eigen::Dynamic, 1>& idarray)
            {
               D = idarray.size();
	       arraydims.push_back(D);
	       arraydims_offset.push_back(D);
	       for(int i = 0; i < D; i++){
		   arraydims[i] = idarray[i];
                   arraydims_offset[i] = idarray[i];
	       }

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);
            }

	    void resize_DummyMD(Eigen::Matrix<T, Eigen::Dynamic, 1>& idarray, T value)
            {
               D = idarray.size();
	       arraydims.push_back(D);
	       arraydims_offset.push_back(D);

	       for(int i = 0; i < D; i++){
		   arraydims[i] = idarray[i];
                   arraydims_offset[i] = idarray[i];
	       }

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);

               #pragma omp parallel for
	       for(int i = 0; i < total_size; i++)
	           DummyMD[i] = value;
            }


//basics operations ............................................
	    int gid(std::vector<int> idarray)
	    {
                int id = 0;
                for(int i = 0; i < D; i++)
                   id += idarray[i] * arraydims_offset[i];
                return id;
	    }

            Eigen::Matrix<int, Eigen::Dynamic, 1> gdim(int& id){
                Eigen::Matrix<int, Eigen::Dynamic, 1> dim;
		dim.resize(D+1); dim.setZero();
		int arraydims_prod_prev = 0;
		double new_id = (double)id;

		for(int i = (D-1); i >= 0; i--){
			int arraydims_prod = 1;

			for(int j = 0; j < i; j++)
				arraydims_prod *= arraydims[j];

			new_id -= arraydims_prod_prev*dim[i+1];
			dim[i] = (int)(new_id/((double)arraydims_prod));
			arraydims_prod_prev = arraydims_prod;
		} 	                
                dim.conservativeResize(D);
		return dim;
	    }

            void s(std::vector<int> idarray, T value)
            {
              int id = gid(idarray);
	      DummyMD[id] = value;
	    }

	    void s(T value)
	    {
	      #pragma omp parallel for
	      for(int i = 0; i < total_size; i++)
	          DummyMD[i] = value;
	    }

	    T g(std::vector<int> idarray)
	    {
	      int id = gid(idarray);
              return DummyMD[id];
	    }

            T* p(std::vector<int> idarray)
	    {
              int id = gid(idarray);
              return DummyMD[id];
            }
};

