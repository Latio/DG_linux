#include "MeshUnion_dim.h"


int MeshUnion_dim::K = 0;
int MeshUnion_dim::Nv = 0;
int	MeshUnion_dim::Ne_inner = 0;
int MeshUnion_dim::Ne_boundary = 0;
int MeshUnion_dim::Nfp = 0;
int MeshUnion_dim::Np = 0;
int MeshUnion_dim::cell_Nv = 0;
int MeshUnion_dim::cell_Nq = 0;
int MeshUnion_dim::two = 2;
int MeshUnion_dim::one = 1;

MeshUnion_dim::MeshUnion_dim()
{
	ncdim_read();
}

MeshUnion_dim::~MeshUnion_dim()
{
}

void MeshUnion_dim::ncdim_read()
{
	static const int NC_ERR = 2;

	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	if (!dataFile.is_valid())
	{
		std::cout << "Couldn't open file!\n";
	}

	NcVar *K_v = dataFile.get_var("K");
	NcVar *Nv_v = dataFile.get_var("Nv");
	NcVar *Ne_inner_v = dataFile.get_var("InnerEdge_Ne");
	NcVar *Ne_boundary_v = dataFile.get_var("BoundaryEdge_Ne");
	NcVar *Nfp_v = dataFile.get_var("InnerEdge_Nfp");
	NcVar *Np_v = dataFile.get_var("cell_Np");
	NcVar *cell_Nv_v = dataFile.get_var("cell_Nv");
	NcVar *cell_Nq_v = dataFile.get_var("cell_Nq");

	K_v->get(&K, 1);
	Nv_v->get(&Nv, 1);
	Ne_inner_v->get(&Ne_inner, 1);
	Ne_boundary_v->get(&Ne_boundary, 1);
	Nfp_v->get(&Nfp, 1);
	Np_v->get(&Np, 1);
	cell_Nv_v->get(&cell_Nv, 1);
	cell_Nq_v->get(&cell_Nq, 1);
}



void MeshUnion_dim::ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = new double[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}

void MeshUnion_dim::ncvar_read(double *&meshunion_data, const char* ncvarname, int &dim1)
{
	meshunion_data = new double[dim1];
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = new int[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}

void MeshUnion_dim::ncvar_read(int *&meshunion_data, const char* ncvarname, int &dim1)
{
	meshunion_data = new int[dim1];
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1);
}

void MeshUnion_dim::ncvar_read(signed char *&meshunion_data, const char* ncvarname, int &dim1, int &dim2)
{
	meshunion_data = new signed char[dim1*dim2];
	NcFile dataFile("meshUnion.nc", NcFile::ReadOnly);
	NcVar *temp_v = dataFile.get_var(ncvarname);
	temp_v->get(meshunion_data, dim1, dim2);
}