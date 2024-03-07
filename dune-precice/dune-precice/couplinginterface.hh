
#ifndef DUNE_PRECICE_HH
#define DUNE_PRECICE_HH

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>

#include <dune/functions/functionspacebases/interpolate.hh>

#include <precice/precice.hpp>

namespace Dune::preCICE {

  namespace Utilities {

    template <int dim, typename VectorType>
    static void copy_from_precice_to_dune(VectorType &dune_data,
                                          const std::vector<double> &precice_data,
                                          const BlockVector<FieldVector<bool, dim>> is_coupling_boundary)
    {
      int iteration_count = 0;
      for(int i=0; i<dune_data.N(); i++) {
        for(int j=0; j<dim; j++) {
          if(is_coupling_boundary[i][j])
          {
            dune_data[i][j] = precice_data[iteration_count];
            ++iteration_count;
          }
        }
      }
    }
  
    template <int dim, typename VectorType>
    static void copy_from_dune_to_precice(const VectorType &dune_data,
                                          std::vector<double> &precice_data,
                                          const BlockVector<FieldVector<bool, dim>> is_coupling_boundary)
    { 
      int iteration_count = 0;
      for(int i=0; i<dune_data.N(); i++) {
        for(int j=0; j<dim; j++) {
          if(is_coupling_boundary[i][j]) 
          { 
            precice_data[iteration_count] = dune_data[i][j];
            ++iteration_count;
          }
        }
      }
    }
  } // namespace Utilities  

  template <int dim, typename VectorType, typename ParameterClass>
  class CouplingInterface
  {

    private:
    
      precice::Participant precice;
      
      const BlockVector<FieldVector<bool, dim>> is_coupling_boundary_;

      const std::string mesh_name_;
      const std::string read_data_name_;
      const std::string write_data_name_;

      int n_interface_nodes_;
      std::vector<precice::VertexID> interface_nodes_ids_;
      precice::span<precice::VertexID> interface_nodes_ids_span;
      
      std::vector<double> read_data_;
      std::vector<double> write_data_;
      precice::span<double> read_data_span;
      precice::span<double> write_data_span;
      const int mpi_rank_;  
      
      struct {
        std::vector<VectorType> data;
        double time;
        int iter;
      } previous_solver_state_;
     
    public:
              
      CouplingInterface(const ParameterClass &parameters, 
                        const BlockVector<FieldVector<bool, dim>> &is_coupling_boundary,
                        int mpi_rank, int mpi_size);

      template <typename PolynomialBasis>
      void initialize(const PolynomialBasis &polynomial_basis);

      double get_max_time_step_size() const;
       
      void read_blockvector_data(VectorType &dune_data, double dt);
      void write_blockvector_data(const VectorType &dune_data);

      void advance(double computed_timestep_length);

      void save_current_state(const std::vector<VectorType*>& state_quantaties,
                              const double current_time_value, 
                              const int current_iter_value);

      void reload_old_state(std::vector<VectorType*>& state_quantaties,
                            double &current_time_value,
                            int &current_iter_value);
                                         
      bool is_save_required() ;
      bool is_load_required() ;
  
      bool is_coupling_ongoing() const;
      bool is_time_window_complete() const;
      void finalize();
  };
  
  template <int dim, typename VectorType, typename ParameterClass>
  CouplingInterface<dim, VectorType, ParameterClass>::CouplingInterface(
    const ParameterClass &parameters,
    const BlockVector<FieldVector<bool, dim>> &is_coupling_boundary,
    const int mpi_rank, const int mpi_size)
    : precice(parameters.participant_name, parameters.config_file, mpi_rank, mpi_size)
    , is_coupling_boundary_(is_coupling_boundary)
    , mpi_rank_(mpi_rank)
    , mesh_name_(parameters.mesh_name)
    , read_data_name_(parameters.read_data_name)
    , write_data_name_(parameters.write_data_name)
  {
  }
    
  template <int dim, typename VectorType, typename ParameterClass>
  template <typename PolynomialBasis>
  void CouplingInterface<dim, VectorType, ParameterClass>::initialize(
    const PolynomialBasis &polynomial_basis) 
  {
    VectorType coordinates_of_boundary;
    auto coordinate = [] (auto x) { return x; };
 
    interpolate(polynomial_basis, coordinates_of_boundary, coordinate, is_coupling_boundary_);

    std::vector<double> interface_nodes_positions;
    const precice::span<double> interface_nodes_positions_span = precice::span(interface_nodes_positions);
    for(int i=0; i<coordinates_of_boundary.N(); i++) {
      for(int j=0; j<dim; j++) {
        if(is_coupling_boundary_[i][j]) {
          interface_nodes_positions.push_back(coordinates_of_boundary[i][j]);
          interface_nodes_ids_.push_back(j*coordinates_of_boundary.N()+i);
          }
      }
    }
        
    {
      FieldVector<bool, dim> indicator;
      indicator = true;
      n_interface_nodes_ = std::count(is_coupling_boundary_.begin(), is_coupling_boundary_.end(), indicator);
    }

    write_data_.resize(dim * n_interface_nodes_);
    read_data_.resize(dim * n_interface_nodes_);
    interface_nodes_positions.resize(dim * n_interface_nodes_);
    interface_nodes_ids_.resize(n_interface_nodes_);

    std::cout << "preCICE:  Mesh on rank "     << mpi_rank_ << " = " << mesh_name_           << std::endl;
    std::cout << "preCICE:  Read on rank "     << mpi_rank_ << " = " << read_data_name_      << std::endl;
    std::cout << "preCICE:  Write on rank "    << mpi_rank_ << " = " << write_data_name_     << std::endl;  
    std::cout << "preCICE:  Vertex size on rank " << mpi_rank_ << " = " << n_interface_nodes_ << std::endl;

    read_data_span = precice::span(read_data_);
    write_data_span = precice::span(write_data_);
    interface_nodes_ids_span = precice::span(interface_nodes_ids_);
    

    precice.setMeshVertices(mesh_name_, interface_nodes_positions, interface_nodes_ids_span);
    
    precice.initialize();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  double CouplingInterface<dim, VectorType, ParameterClass>::get_max_time_step_size() const
  {
    return precice.getMaxTimeStepSize();
  }
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::read_blockvector_data(VectorType &dune_data, double dt)
  { 
    precice.readData(mesh_name_, read_data_name_, interface_nodes_ids_, dt, read_data_span);
    Utilities::copy_from_precice_to_dune<dim, VectorType>(dune_data, read_data_, is_coupling_boundary_);
  }

  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::write_blockvector_data(const VectorType &dune_data)
  {
    Utilities::copy_from_dune_to_precice<dim, VectorType>(dune_data, write_data_, is_coupling_boundary_);
    precice.writeData(mesh_name_, write_data_name_, interface_nodes_ids_span, write_data_span);
  }
      
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::advance(const double computed_timestep_length)
  {
    precice.advance(computed_timestep_length);
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_coupling_ongoing() const
  {
    return precice.isCouplingOngoing();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_time_window_complete() const
  {
    return precice.isTimeWindowComplete();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::finalize()
  {
    precice.finalize();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_save_required() 
  {
    return precice.requiresWritingCheckpoint();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  bool CouplingInterface<dim, VectorType, ParameterClass>::is_load_required() 
  {
    return precice.requiresReadingCheckpoint();
  }
  
  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::save_current_state(
    const std::vector<VectorType*>& state_quantaties,
    const double time,
    const int iter) 
  {  
    previous_solver_state_.data.resize(state_quantaties.size());
    
    for (int i=0; i<state_quantaties.size(); i++) {
      previous_solver_state_.data[i] = *(state_quantaties[i]);
    }

    previous_solver_state_.time = time;
    previous_solver_state_.iter = iter; 
  }

  template <int dim, typename VectorType, typename ParameterClass>
  void CouplingInterface<dim, VectorType, ParameterClass>::reload_old_state(
    std::vector<VectorType*>& state_quantaties,
    double &time,
    int &iter)
  {
    assert( state_quantaties.size() == previous_solver_state_.data.size() );

    for (int i=0; i<state_quantaties.size(); i++) {
      *(state_quantaties[i]) = previous_solver_state_.data[i];
    }
   
    time = previous_solver_state_.time;
    iter = previous_solver_state_.iter;
  }
} // namespace Dune

#endif // DUNE_PRECICE_HH
