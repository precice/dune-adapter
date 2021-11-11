// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PRECICE_HH
#define DUNE_PRECICE_HH

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include <dune/istl/io.hh>
#include <dune/istl/bvector.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include <precice/SolverInterface.hpp>

namespace Dune::preCICE {

  template <int dim, typename VectorType, typename ParameterClass>
  class CouplingInterface
  {

    private:
    
      precice::SolverInterface precice;
      
      BlockVector<FieldVector<bool, dim>> isBoundary;

      const std::string mesh_name;
      const std::string read_data_name;
      const std::string write_data_name;
    
      int mesh_id;
      int read_data_id;
      int write_data_id;
      
      int n_interface_nodes;
      std::vector<double> interface_nodes_positions;
      std::vector<int>    interface_nodes_ids;
      
      std::vector<double> read_data;
      std::vector<double> write_data;
      
      std::vector<VectorType> old_state_data;
      double old_time;
      int old_iter;
      
      void format_dune_to_precice(VectorType &dune_to_precice);
      void format_precice_to_dune(VectorType &precice_to_dune);
      
      int mpiRank_;
      
    public:
          
      double precice_timestep_length;
    
      CouplingInterface(const ParameterClass &parameters, BlockVector<FieldVector<bool, dim>>& boundaryDofs, int mpi_rank, int mpi_size);

      template <typename Basis>
      void initialize(const Basis& basis);
       
      void send_blockvector_data(VectorType &dune_to_precice);
      void read_blockvector_data(VectorType &precice_to_dune);

      void advance(double computed_timestep_length);
       
      void save_current_state(std::vector<VectorType*> state_quantaties, double &current_time_value, int &current_iter_value);
      void reload_old_state(std::vector<VectorType*> state_quantaties, double &current_time_value, int &current_iter_value);
                                         
      bool is_save_required();
      bool is_load_required();
  
      void mark_save_fullfilled();
      void mark_load_fullfilled(); 
        
      bool is_coupling_ongoing();
      bool is_time_window_complete();
      void finalize();
        
  };
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  CouplingInterface<dim, VectorType, ParameterClass>::CouplingInterface(const ParameterClass &parameters,
                                                                        BlockVector<FieldVector<bool,dim>>& boundaryDofs,
                                                                        int mpi_rank, int mpi_size)
    : precice(parameters.participant_name, parameters.config_file, mpi_rank, mpi_size)
    , isBoundary(boundaryDofs)
    , mpiRank_(mpi_rank)
    , mesh_name(parameters.mesh_name)
    , read_data_name(parameters.read_data_name)
    , write_data_name(parameters.write_data_name) {}
    
  
  template <int dim, typename VectorType, typename ParameterClass>
  template <typename Basis>
  void CouplingInterface<dim, VectorType, ParameterClass>::initialize(
    const Basis& basis
  ) {
  
	mesh_id       = precice.getMeshID(mesh_name);
	read_data_id  = precice.getDataID(read_data_name, mesh_id);
	write_data_id = precice.getDataID(write_data_name, mesh_id);
            
    // interpolate boundary DOFS to coordinates
    VectorType coordinatesBoundary;
    auto coordinate = [] (auto x) { return x; };
 
    interpolate(basis, coordinatesBoundary, coordinate, isBoundary);

    // set up precice interface
    for(int i=0; i<coordinatesBoundary.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(isBoundary[i][j]) interface_nodes_positions.push_back(coordinatesBoundary[i][j]);
      }
    }
        
    FieldVector<bool, dim> indicator;
    indicator = true;
    n_interface_nodes = std::count(isBoundary.begin(), isBoundary.end(), indicator);
    
    write_data.resize(dim * n_interface_nodes);
    read_data.resize(dim * n_interface_nodes);
    interface_nodes_positions.resize(dim * n_interface_nodes);
    interface_nodes_ids.resize(n_interface_nodes);

    std::cout << "preCICE:  Mesh ID on rank " << mpiRank_ << " = " << mesh_id << std::endl;
    std::cout << "preCICE:  Read ID on rank " << mpiRank_ << " = " << read_data_id << std::endl;
    std::cout << "preCICE:  Write ID on rank " << mpiRank_ << " = " << write_data_id << std::endl;  
    std::cout << "preCICE:  Vertex size on rank " << mpiRank_ << " = " << n_interface_nodes << std::endl;
       
    precice.setMeshVertices(mesh_id, n_interface_nodes, interface_nodes_positions.data(), interface_nodes_ids.data());
    precice_timestep_length = precice.initialize();
    
  }
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::read_blockvector_data(VectorType &precice_to_dune) {
  
    // modify something here for nonoverlapping meshes ?
    precice.readBlockVectorData(read_data_id, n_interface_nodes, interface_nodes_ids.data(), read_data.data());
    format_precice_to_dune(precice_to_dune);
    
  }
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::send_blockvector_data(VectorType &dune_to_precice) {

    // modify something here for nonoverlapping meshes ?
    format_dune_to_precice(dune_to_precice);
    precice.writeBlockVectorData(write_data_id, n_interface_nodes, interface_nodes_ids.data(), write_data.data());

  }
    
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::format_precice_to_dune(VectorType &precice_to_dune) {
  
    int iteration_count = 0;
    for(int i=0; i<precice_to_dune.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(isBoundary[i][j]) { precice_to_dune[i][j] = read_data[iteration_count]; ++iteration_count; }
      }
    }
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::format_dune_to_precice(VectorType &dune_to_precice) {
   
    int iteration_count = 0;
    for(int i=0; i<dune_to_precice.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(isBoundary[i][j]) { write_data[iteration_count] = dune_to_precice[i][j]; ++iteration_count; }
      }
    }  
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::advance(double computed_timestep_length) {
    precice_timestep_length = precice.advance(computed_timestep_length);
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_coupling_ongoing() {
    return precice.isCouplingOngoing();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_time_window_complete() {
    return precice.isTimeWindowComplete();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::finalize() {
    precice.finalize();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_save_required() {
    return precice.isActionRequired(precice::constants::actionWriteIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::mark_save_fullfilled() {
    precice.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_load_required() {
    return precice.isActionRequired(precice::constants::actionReadIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::mark_load_fullfilled() {
    precice.markActionFulfilled(precice::constants::actionReadIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::save_current_state(
    std::vector<VectorType*> state_quantaties,
    double &time,
    int &iter
  ) {
  
    old_state_data.resize(state_quantaties.size());
    
    for (int i=0; i<state_quantaties.size(); i++) {
      old_state_data[i] = *(state_quantaties[i]);
    }

    old_time = time;
    old_iter = iter;
      
  }

  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::reload_old_state(
    std::vector<VectorType*> state_quantaties,
    double &time,
    int &iter
  ) {
  
    for (int i=0; i<state_quantaties.size(); i++) {
      *(state_quantaties[i]) = old_state_data[i];
    }
   
    time = old_time;
    iter = old_iter;

  }
  
} // namespace Dune

#endif // DUNE_PRECICE_HH
