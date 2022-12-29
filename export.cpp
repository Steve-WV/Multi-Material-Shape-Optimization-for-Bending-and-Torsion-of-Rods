#include <iomanip>
#include <fstream>
#include <unordered_map>
#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"

#ifdef HAVE_VTK
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkFloatArray.h>

 void ImplicitPDESystem::Export(std::string fname) const 
{ 
	std::lock_guard<std::mutex> lock(mut);
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(m->VertexSize());
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd(); 
	float p[3]; 
	for (int j=0; vi!=ve; ++vi) {
		p[0] = (*vi)->x();
		p[1] = (*vi)->y();
		p[2] = 0.0;
		points->SetPoint(j++, p);
	}
	
	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) Vt_list_idx[*vi] = j++;
 	
	auto Conn = vtkSmartPointer<vtkIdTypeArray>::New();
	Conn->SetNumberOfValues(4*m->TriSize());
	
	auto ti = m->TriBegin(), te = m->TriEnd(); 
	for (int j = 0; ti!=te; ++ti) {
		Conn->SetValue( j++, 3 );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->a()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->b()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->c()] );
	}

	auto Tris = vtkSmartPointer<vtkCellArray>::New();
	Tris->SetCells(m->TriSize(),Conn);
	
	auto poly = vtkSmartPointer<vtkPolyData>::New();
	poly->SetPoints(points);
	poly->SetPolys(Tris);
	
	auto u = vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(1);
	u->SetName("u");
	u->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) u->SetTuple1(j++, (*vi)->u());
	poly->GetPointData()->AddArray(u);
	
	auto phi = vtkSmartPointer<vtkFloatArray>::New();
	phi->SetNumberOfComponents(1);
	phi->SetName("phi");
	phi->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) phi->SetTuple1(j++, (*vi)->Phi() );
	poly->GetPointData()->AddArray(phi);
	
	auto phase = vtkSmartPointer<vtkFloatArray>::New();
	phase->SetNumberOfComponents(1);
	phase->SetName("phase");
	phase->SetNumberOfValues(m->VertexSize());
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) phase->SetTuple1(j++, 0.5*(1.0+std::erf( 6.0*((*vi)->u()-0.5))));
	poly->GetPointData()->AddArray(phase);
	
	// Write file
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	std::string fn = fname+".vtp";
	writer->SetFileName(fn.c_str());
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(poly);
#else
	writer->SetInputData(poly);
#endif
	writer->Write();
	  
}
#endif
#ifndef HAVE_VTK
void ImplicitPDESystem::Export(std::string fname) const {
	
	std::lock_guard<std::mutex> lock(mut);
	std::ofstream os;
	os.open(fname+".vtk");
	os << std::setprecision(16);
	os << "# vtk DataFile Version 2.0" << std::endl;
	os << "s2d output data" << std::endl;
	os << "ASCII" << std::endl;
	os << "DATASET POLYDATA" << std::endl;
	
	os << "POINTS " << m->VertexSize() << " float" << std::endl;
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi)
		os << (*vi)->x() << " " << (*vi)->y() << " " << "0.0" << std::endl;

	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); int idx = 0;
	for (; vi!=ve; ++vi) {
		Vt_list_idx[*vi] = idx;
		++idx;
	}

	os << "POLYGONS " << m->TriSize() << " " << m->TriSize()*4 << std::endl;
	Mesh::const_TriIt ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti)
		os << "3 " << Vt_list_idx[(*ti)->a()] << " " << Vt_list_idx[(*ti)->b()] << " "
			<< Vt_list_idx[(*ti)->c()] << " " << std::endl;
    
	os << "POINT_DATA " << m->VertexSize() << std::endl;
	os << "SCALARS u float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi)
		os << (*vi)->u() << std::endl;
	
	
    
	os << "CELL_DATA " << m->TriSize() << std::endl;
	os << "SCALARS Tri_area float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TriBegin();
	for (; ti!=te; ++ti)
		os << (*ti)->Area() << std::endl;
	os.close();
    
}
#endif
