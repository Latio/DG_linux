#include"NdgPhysMat.h"
//#define _CRT_SECURE_NO_WARNINGS
MeshUnion mesh;
const MeshUnion *meshunion = &mesh;
int main()
{
	//for (size_t i = 0; i < 1080; i++)
	//{
	//	std::cout << *(meshunion->EToE+i) << std::endl;
	//}
	NdgPhysMat Solver;
	Solver.matSolver();

	system("pause");
	return 0;
}
