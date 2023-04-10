#ifndef WRITE_TO_VTK_H
#define WRITE_TO_VTK_H
///////////////////////////////////////////////////////////////////////////
#include<stdexcept>
#include<string>
#include<vtkCellArray.h>
#include<vtkDoubleArray.h>
#include<vtkCellData.h>
#include<vtkPointData.h>
#include<vtkPoints.h>
#include<vtkXMLStructuredGridWriter.h>
#include<vtkStructuredGrid.h>
#include<vtkSmartPointer.h>
#include"Output_Messages.h"
#include<string>
#include<cstring>

namespace homedf {

/////////////////////////////////////////////////////////////////////////
template<typename Func_type>
void Write_Distribution_Function_to_data(Func_type& F,
					 int num_points_v,
					 std::string filename) {

  Output_Begin_Status("Starting process to Output Distribution Function to data file");
  
  Output_Detail_Status("Outputing Distribution Function to data file");

  if(filename.size() < 4){
    filename = filename + ".dat";
  } else if(filename.compare(filename.size()-4, 4, ".dat") != 0){
    filename = filename + ".dat";
  }
  
  std::ofstream outputFile;
  outputFile.open(filename);
  outputFile << std::left <<
    std::setw(15) << "v" <<
    std::setw(15) << "F" <<
    std::endl << std::endl;

  double dv = (F.v_max - F.v_min)/static_cast<double>(num_points_v - 1);

  for(int i = 0; i < num_points_v-1; i++){
    double v = F.v_min + static_cast<double>(i)*dv;
    outputFile << std::left <<
      std::setw(15) << v <<
      std::setw(15) << F(v) <<
      std::endl;
  }

  outputFile.close();
  
  Output_End_Status("The Distribution Function was Successfully Outputed to "+filename+".dat");
  
}
  
/////////////////////////////////////////////////////////////////////////
template<typename Func_type>
void Write_Distribution_Function_to_VTK(Func_type& F,
					int num_points_vx,
					int num_points_vy,
					int num_points_vz,
					std::string filename) {

  Output_Begin_Status("Starting process to Output Distribution Function to VTK");
  
  Output_Detail_Status("Outputing Distribution Function to VTK");
  write_field_3D_Cartesian(F.vx_min, F.vx_max, num_points_vx,
			   F.vy_min, F.vy_max, num_points_vy,
			   F.vz_min, F.vz_max, num_points_vz,
			   F, filename);
  Output_End_Status("The Distribution Function was Successfully Outputed to "+filename+".vts");
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Func_type>
inline void write_field_3D_Cartesian(double x_min, double x_max, int nx,
				     double y_min, double y_max, int ny,
				     double z_min, double z_max, int nz,
				     const Func_type& func,
				     std::string filename) {

  if(nx < 2 || ny < 2 || nz < 2) {
    throw std::range_error("nx, ny and nz must be larger than 2 in \"write_field_3D_Cartesian\".");
  }

  if(x_max <= x_min) {
    throw std::range_error("x_max must be larger than x_min in \"write_field_3D_Cartesian\".");
  }

  if(y_max <= y_min) {
    throw std::range_error("y_max must be larger than y_min in \"write_field_3D_Cartesian\".");
  }

  if(z_max <= z_min) {
    throw std::range_error("z_max must be larger than z_min in \"write_field_3D_Cartesian\".");
  }

  if(filename.size() < 4){
    filename = filename + ".vts";
  } else if(filename.compare(filename.size()-4, 4, ".vts") != 0){
    filename = filename + ".vts";
  }

  auto structuredGrid  = vtkSmartPointer<vtkStructuredGrid>::New();
  auto points          = vtkSmartPointer<vtkPoints>::New();
  auto function_values = vtkSmartPointer<vtkDoubleArray>::New();

  function_values->SetName("Value");

  const auto dx = (x_max-x_min)/static_cast<double>(nx-1);
  const auto dy = (y_max-y_min)/static_cast<double>(ny-1);
  const auto dz = (z_max-z_min)/static_cast<double>(nz-1);

  for(int i=0; i<nx; i++){
    for(int j=0; j<ny; j++){
      for(int k=0; k<nz; k++){

	const auto x = x_min + static_cast<double>(i)*dx;
	const auto y = y_min + static_cast<double>(j)*dy;
	const auto z = z_min + static_cast<double>(k)*dz;
    
	points->InsertNextPoint(x, y, z);

	if(i>0 && j>0 && k>0){
	  const auto a   = x - dx/2.0;
	  const auto b   = y - dy/2.0;
	  const auto c   = z - dz/2.0;
	  const auto val = func(a, b, c);
	  function_values->InsertNextValue(val);
	}

      }
    }
  }

  structuredGrid->SetDimensions(nx,ny,nz);
  structuredGrid->SetPoints(points);
  structuredGrid->GetCellData()->SetScalars(function_values);

  auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(structuredGrid);
  writer->Write();

 }


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Vec_type, typename Sol_type>
inline void Write_2D_field(Vec_type x,
			   Vec_type y,
			   const Sol_type& nf,
			   std::string filename) {

  Output_Begin_Status("Starting process to Output Distribution Function to VTK");
  Output_Detail_Status("Outputing Distribution Function to VTK");
  
  int nx = x.size();
  int ny = y.size();
  int num_cells = nx*ny;
  int num_nodes = (nx+1)*(ny+1);
  double x_min = -10.0;
  double y_min = -10.0;
  double delta_x = 20.0/nx;
  double delta_y = 20.0/ny;
  
  if(filename.size() < 4){
    filename = filename + ".vtk";
  } else if(filename.compare(filename.size()-4, 4, ".vtk") != 0){
    filename = filename + ".vtk";
  }

  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  fout.precision(16);


  fout << "# vtk DataFile Version 2.0\n"
       << "Distribution function\n"
       << "ASCII\n"
       << "DATASET STRUCTURED_GRID\n"
       << "DIMENSIONS " << nx+1 << " " << ny+1 << " 1\n"
       << "POINTS " << num_nodes << " double\n";

  for(int j = 0; j <= ny; ++j) {
    for(int i = 0; i <= nx; ++i) {
      fout << x_min+static_cast<double>(i)*delta_x << " "
	   << y_min+static_cast<double>(j)*delta_y << " 0.0\n";
    }
  }

  // temperature
  fout << "\nCELL_DATA " << num_cells
       << "\nSCALARS n double 1"
       << "\nLOOKUP_TABLE default\n";
  for(int j = 0; j < ny; ++j) {
    for(int i = 0; i < nx; ++i) {
        fout << nf[i][j] << '\n';
    }
  }
  Output_End_Status("The Distribution Function was Successfully Outputed to "+filename+".vtk");
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Vec_type, typename Sol_type>
inline void Write_2D_scatter(Vec_type x,
			     Vec_type y,
			     Vec_type s,
			     const Sol_type& ev,
			     int nx,
			     int ny,
			     std::string filename) {

  Output_Begin_Status("Starting process to Output 2D Scatter to VTK");
  Output_Detail_Status("Outputing 2D Scatter to VTK");
  
  int n = x.size();
  int num_scalars = ev[0].size() + 1;
  
  if(filename.size() < 4){
    filename = filename + ".vts";
  } else if(filename.compare(filename.size()-4, 4, ".vts") != 0){
    filename = filename + ".vts";
  }

  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  fout.precision(16);


  auto structuredGrid  = vtkSmartPointer<vtkStructuredGrid>::New();
  auto points          = vtkSmartPointer<vtkPoints>::New();
  
  std::vector<vtkSmartPointer<vtkDoubleArray>> scalar_values;
  scalar_values.resize(num_scalars);
  for(int i = 0; i < num_scalars; ++i) {
    scalar_values[i] = vtkSmartPointer<vtkDoubleArray>::New();
    if(i == 0){
      scalar_values[i]->SetName("s");
    } else {
      scalar_values[i]->SetName(("EV" + std::to_string(i)).c_str());
    }
  }
  
  
  for(int i=0; i<n; i++){
    points->InsertNextPoint(x[i], y[i], 0.0);
    scalar_values[0]->InsertNextValue(s[i]);
    for(int j = 1; j < num_scalars; ++j) {
      scalar_values[j]->InsertNextValue(ev[i][j]);
    }
  }

  structuredGrid->SetDimensions(nx,ny,1);
  structuredGrid->SetPoints(points);
 
  structuredGrid->GetPointData()->SetScalars(scalar_values[0]);
  for(int i = 1; i < num_scalars; ++i) {
    structuredGrid->GetPointData()->AddArray(scalar_values[i]);
  }

  auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(structuredGrid);
  writer->Write();

  Output_End_Status("The 2D Scatter was Successfully Outputed to "+filename+".vtk");
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template<typename Vec_type>
inline void Write_2D_scatter(Vec_type x,
			     Vec_type y,
			     Vec_type s,
			     int nx,
			     int ny,
			     std::string filename) {

  Output_Begin_Status("Starting process to Output 2D Scatter to VTK");
  Output_Detail_Status("Outputing 2D Scatter to VTK");
  
  int n = x.size();
  
  if(filename.size() < 4){
    filename = filename + ".vts";
  } else if(filename.compare(filename.size()-4, 4, ".vts") != 0){
    filename = filename + ".vts";
  }

  std::ofstream fout(filename);
  if(!fout) {
    throw std::runtime_error("Could not open file: " + filename);
  }
  fout.precision(16);


  auto structuredGrid  = vtkSmartPointer<vtkStructuredGrid>::New();
  auto points          = vtkSmartPointer<vtkPoints>::New();
  
  auto scalar_value = vtkSmartPointer<vtkDoubleArray>::New();
  scalar_value->SetName("s");
  
  for(int i=0; i<n; i++){
    points->InsertNextPoint(x[i], y[i], 0.0);
    scalar_value->InsertNextValue(s[i]);
  }

  structuredGrid->SetDimensions(nx,ny,1);
  structuredGrid->SetPoints(points);
  structuredGrid->GetPointData()->SetScalars(scalar_value);

  auto writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName(filename.c_str());
  writer->SetInputData(structuredGrid);
  writer->Write();

  Output_End_Status("The 2D Scatter was Successfully Outputed to "+filename+".vtk");
}

  
} //namespace homedf

#endif //#ifndef WRITE_TO_VTK_H

