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
      
      const BlockVector<FieldVector<bool, dim>> is_coupling_boundary_;

      const std::string mesh_name_;
      const std::string read_data_name_;
      const std::string write_data_name_;
    
      int mesh_id_;
      int read_data_id_;
      int write_data_id_;
      
      int n_interface_nodes_;
      std::vector<double> interface_nodes_positions_;
      std::vector<int>    interface_nodes_ids_;
      
      std::vector<double> read_data_;
      std::vector<double> write_data_;
      
      std::vector<VectorType> old_state_data_;
      double old_time_;
      int old_iter_;
     
      void format_precice_to_dune(VectorType &precice_to_dune);
      void format_dune_to_precice(const VectorType &dune_to_precice);     
 
      int mpi_rank_;  

    public:
              
      CouplingInterface(const ParameterClass &parameters, BlockVector<FieldVector<bool, dim>>& boundaryDofs, int mpi_rank, int mpi_size);

      template <typename PolynomialBasis>
      double initialize(const PolynomialBasis& polynomial_basis);
       
      void read_blockvector_data(VectorType &precice_to_dune);
      void send_blockvector_data(const VectorType &dune_to_precice);

      double advance(double computed_timestep_length);
       
      void save_current_state(std::vector<VectorType*> state_quantaties, double &current_time_value, int &current_iter_value);
      void reload_old_state(std::vector<VectorType*> state_quantaties, double &current_time_value, int &current_iter_value);
                                         
      bool is_save_required() const;
      bool is_load_required() const;
  
      void mark_save_fullfilled();
      void mark_load_fullfilled(); 
        
      bool is_coupling_ongoing() const;
      bool is_time_window_complete() const;
      void finalize();
        
  };
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  CouplingInterface<dim, VectorType, ParameterClass>::CouplingInterface(const ParameterClass &parameters,
                                                                        BlockVector<FieldVector<bool,dim>>& boundaryDofs,
                                                                        int mpi_rank, int mpi_size)
    : precice(parameters.participant_name, parameters.config_file, mpi_rank, mpi_size)
    , is_coupling_boundary_(boundaryDofs)
    , mpi_rank_(mpi_rank)
    , mesh_name_(parameters.mesh_name)
    , read_data_name_(parameters.read_data_name)
    , write_data_name_(parameters.write_data_name) {}
    
  
  template <int dim, typename VectorType, typename ParameterClass>
  template <typename PolynomialBasis>
  double CouplingInterface<dim, VectorType, ParameterClass>::initialize(
    const PolynomialBasis& polynomial_basis
  ) {
  
	mesh_id_       = precice.getMeshID(mesh_name_);
	read_data_id_  = precice.getDataID(read_data_name_, mesh_id_);
	write_data_id_ = precice.getDataID(write_data_name_, mesh_id_);
            
    VectorType coordinates_of_boundary;
    auto coordinate = [] (auto x) { return x; };
 
    interpolate(polynomial_basis, coordinates_of_boundary, coordinate, is_coupling_boundary_);

    for(int i=0; i<coordinates_of_boundary.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(is_coupling_boundary_[i][j]) interface_nodes_positions_.push_back(coordinates_of_boundary[i][j]);
      }
    }
        
    {
      FieldVector<bool, dim> indicator;
      indicator = true;
      n_interface_nodes_ = std::count(is_coupling_boundary_.begin(), is_coupling_boundary_.end(), indicator);
    }

    write_data_.resize(dim * n_interface_nodes_);
    read_data_.resize(dim * n_interface_nodes_);
    interface_nodes_positions_.resize(dim * n_interface_nodes_);
    interface_nodes_ids_.resize(n_interface_nodes_);

    std::cout << "preCICE:  Mesh ID on rank "     << mpi_rank_ << " = " << mesh_id_           << std::endl;
    std::cout << "preCICE:  Read ID on rank "     << mpi_rank_ << " = " << read_data_id_      << std::endl;
    std::cout << "preCICE:  Write ID on rank "    << mpi_rank_ << " = " << write_data_id_     << std::endl;  
    std::cout << "preCICE:  Vertex size on rank " << mpi_rank_ << " = " << n_interface_nodes_ << std::endl;
       
    precice.setMeshVertices(mesh_id_, n_interface_nodes_, interface_nodes_positions_.data(), interface_nodes_ids_.data());
    
    return precice.initialize();
  }
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::read_blockvector_data(VectorType &precice_to_dune) {
  
    precice.readBlockVectorData(read_data_id_, n_interface_nodes_, interface_nodes_ids_.data(), read_data_.data());
    format_precice_to_dune(precice_to_dune);
  }
  
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::send_blockvector_data(const VectorType &dune_to_precice) {

    format_dune_to_precice(dune_to_precice);
    precice.writeBlockVectorData(write_data_id_, n_interface_nodes_, interface_nodes_ids_.data(), write_data_.data());
  }
    
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::format_precice_to_dune(VectorType &precice_to_dune) {
  
    int iteration_count = 0;
    for(int i=0; i<precice_to_dune.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(is_coupling_boundary_[i][j]) { precice_to_dune[i][j] = read_data_[iteration_count]; ++iteration_count; }
      }
    }
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::format_dune_to_precice(const VectorType &dune_to_precice) {
   
    int iteration_count = 0;
    for(int i=0; i<dune_to_precice.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(is_coupling_boundary_[i][j]) { write_data_[iteration_count] = dune_to_precice[i][j]; ++iteration_count; }
      }
    }  
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  double CouplingInterface<dim, VectorType, ParameterClass>::advance(double computed_timestep_length) {
    return precice.advance(computed_timestep_length);
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_coupling_ongoing() const {
    return precice.isCouplingOngoing();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_time_window_complete() const {
    return precice.isTimeWindowComplete();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::finalize() {
    precice.finalize();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_save_required() const {
    return precice.isActionRequired(precice::constants::actionWriteIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::mark_save_fullfilled() {
    precice.markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_load_required() const {
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
  
    old_state_data_.resize(state_quantaties.size());
    
    for (int i=0; i<state_quantaties.size(); i++) {
      old_state_data_[i] = *(state_quantaties[i]);
    }

    old_time_ = time;
    old_iter_ = iter;
      
  }

  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::reload_old_state(
    std::vector<VectorType*> state_quantaties,
    double &time,
    int &iter
  ) {
  
    for (int i=0; i<state_quantaties.size(); i++) {
      *(state_quantaties[i]) = old_state_data_[i];
    }
   
    time = old_time_;
    iter = old_iter_;

  }
  
} // namespace Dune

#endif // DUNE_PRECICE_HH
