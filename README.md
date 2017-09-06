# statisticaldistances
Calculation of calpha distances of proteins from physically interacting residues

The provided C++ script (pdb_calpha_dist.cpp) can be used to identify physically interacting atoms among distinct residue pairs in a protein structure (usual PDB format required). A list of each interacting pair of atoms followed by their euclidian distances is provided. Moreover, a list of the Ca-Ca distances for this interacting residues is also provided, which can be used as an input of structural information for molecular modelling studies.

A description of files follows:
  1) pdb_calpha_dist.cpp - C++ code, compile to your OS with:  
           
           g++ -std=c++11 pdb_calpha_dist.cpp -o pdb_calpha_dist.exe
  
  2) pdb_calpha_dist.exe - binary executable file (Linux-x86_64-multicore () 
  3) pdb_distances_multiple.sh - bash script to run analysis for a group of PDB structures. Having all desired pdbs in your folder, run the multiple calculation to all models with the command above. Outputs are generated individually for each structure.
  
          bash pdb_distances_multiple.sh
  
  
  
  
