#include "util.h"



bool zeroTest(unsigned i, unsigned j)
{  
   // diagonal
   //return i == j;
   // checkerboard => 50% sparsity
   return  ((i+j) % 2 == 0);
   // dense 
   //return  true;
}



void print_header(unsigned n, char const* header)
{
   std::string s(" test_");
   s += std::to_string(n) + ": " + std::string(header);

   unsigned len(s.length());
   std::cout << std::string(len+1, '=') << std::endl;
   std::cout << s << std::endl;
   std::cout << std::string(len+1, '=') << std::endl;
   
}



std::vector<double> readFile(std::string const& filename)
{
   std::vector<double> mat;
   std::ifstream ifs(filename.c_str(), std::ios::in);
   if (!ifs.is_open()) {
      std::cerr << "Failed to open flie " << filename << std::endl;
      return mat;
   }

   double x;
   while (ifs >> x) { mat.push_back(x); }

   ifs.close();

   return mat;
}



void readMatrix(std::string const& fname, TileArray<double>& TA)
{

    std::ifstream output;
    output.open(fname, std::ios::binary);
    
    int n_csf_classes;
    output.read(reinterpret_cast<char*>(&n_csf_classes), sizeof(int));
    std::cout << n_csf_classes << std::endl;

    TA.resize(n_csf_classes,n_csf_classes);

    std::vector<int> csf_sizes(n_csf_classes);

    for (int i = 0; i < n_csf_classes; i++) {
        int size;
        output.read(reinterpret_cast<char*>(&size), sizeof(int));
        std::cout << size << std::endl;
        csf_sizes[i] = size;
    }

    for (int I = 0; I < n_csf_classes; I++) {
        for (int J = I; J < n_csf_classes; J++) {
  
            int type;
            output.read(reinterpret_cast<char*>(&type), sizeof(int));
            std::cout << "type = " << type << std::endl;

            int csf0_size = csf_sizes[I];
            int csf1_size = csf_sizes[J];
            
            // Zero block
            if (type == 0) {
               TA.set(I,J, new ZeroTile<double>(csf0_size, csf1_size));
               TA.set(J,I, new ZeroTile<double>(csf1_size, csf0_size));

            // Stripe block
            } else if (type == 2) {

               int n_stripes;
               output.read(reinterpret_cast<char*>(&n_stripes), sizeof(int));
               std::cout << "number of stripes " << n_stripes << std::endl;

               StripedTile<double>* IJ = new StripedTile<double>(csf0_size, csf1_size, n_stripes);
               StripedTile<double>* JI = new StripedTile<double>(csf1_size, csf0_size, n_stripes);
               IJ->alloc();
               JI->alloc();

               TA.set(I,J, IJ);
               TA.set(J,I, JI);

               double* IJdata = IJ->data();
               double* JIdata = JI->data();

               std::vector<int> stripes;
               int ld(std::min(csf0_size, csf1_size));

                for (int s = 0; s < n_stripes; s++) {
                    int diag_off;
                    output.read(reinterpret_cast<char*>(&diag_off), sizeof(int));
                    std::cout << diag_off << std::endl;
                    stripes.push_back(diag_off);
                    int size = std::min(csf0_size, csf1_size - diag_off) - std::max(0, - diag_off);
                    for (int i = 0; i < size; i++) {
                        double res;
                        output.read(reinterpret_cast<char*>(&res), sizeof(double));
                        std::cout << res << std::endl;

/*
               std::cout << "Setting IJ[" <<ld*s+i << "] = " << res << std::endl;
               std::cout << "Setting JI[" <<ld*(n_stripes-s-1)+i << "] = " << res << std::endl;
*/
                        IJdata[ld*s+i] = res;
                        JIdata[ld*(n_stripes-s-1)+i] = res;
                    }
                }

               IJ->setStripes(stripes);

               reverse(stripes.begin(), stripes.end());
               for (int s = 0; s < n_stripes; s++) {
                   stripes[s] = -stripes[s];
               }

               JI->setStripes(stripes);

            // Dense block
            } else if (type == 3) {
               CMTile<double> *IJ; 
               CMTile<double> *JI; 

               IJ = new CMTile<double>(csf0_size, csf1_size);
               IJ->alloc();

               if (I != J) {
                  JI = new CMTile<double>(csf1_size, csf0_size);
                  JI->alloc();
               }

               TA.set(I,J, IJ);
               if (I != J) TA.set(J,I, JI);

                for (int i = 0; i < csf_sizes[I]; i++) {
                    for (int j = 0; j < csf_sizes[J]; j++) {
                        double res;
                        output.read(reinterpret_cast<char*>(&res), sizeof(double));
                        if (I != J) {
                           IJ->set(i,j, res);
                           JI->set(j,i, res);
                        }else if (i <= j) {
                           IJ->set(i,j, res);
                           IJ->set(j,i, res);
                        }
                        std::cout << "(" << i << "," << j << ") = " << res << std::endl;
                    }
                }
            } 
        }
    }

}

