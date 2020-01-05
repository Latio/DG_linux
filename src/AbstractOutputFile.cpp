#include "AbstractOutputFile.h"

AbstractOutputFile::AbstractOutputFile(const char* NCfile_name, double timeInerval, int StepPerFile) :timePrevious(0), outputStep(0), resultFile(NCfile_name, NcFile::Replace)
{
	this->timeInerval = timeInerval;
	this->StepPerFile = StepPerFile;
}


AbstractOutputFile::~AbstractOutputFile()
{
	resultFile.close();
}

void AbstractOutputFile::outputIntervalResult(double &time, double *field, int Nvar, int *Np, int *K)
{
	if ((time - timePrevious) > timeInerval)
	{
		outputResult(time, field, Nvar, Np, K);
		timePrevious = time;
	}
};

void AbstractOutputFile::outputResult(double time, double *field, int Nvar, int *Np, int *K)
{
	//int start_time[1] = { outputStep };
	//int count_time[1] = { 1 };
	//std::vector<size_t> startInd_time(start_time, start_time + 1);
	////startInd_time.push_back(outputStep);
	//std::vector<size_t> countInd_time(count_time, count_time + 1);
	////countInd_time.push_back(1);
	//output_time.putVar(startInd_time, countInd_time, &time);

	//int start[4] = { outputStep,0,0,0 };
	//int count[4] = { 1,Nvar,*K,*Np };
	//std::vector<size_t> startInd_fphys(start, start + sizeof(start) / sizeof(start[0]));
	//std::vector<size_t> countInd_fphys(count, count + sizeof(count) / sizeof(count[0]));
	//output_fphys.putVar(startInd_fphys, countInd_fphys, field);

	output_time->set_cur(outputStep);
	output_time->put(&time, 1);

	output_fphys->set_cur(outputStep, 0, 0, 0);
	output_fphys->put(field, 1, Nvar, *K, *Np);

	//int start[4] = { outputStep,0,0,0 };
	//int count[4] = { 1,Nvar,*K,*Np };

	//std::vector<size_t> startInd_fphys(start, start + sizeof(start) / sizeof(start[0]));
	//std::vector<size_t> countInd_fphys(count, count + sizeof(count) / sizeof(count[0]));

	//output_fphys->put(resultFile.get_dim("Nt"), field, outputStep);

	//output_fphys = resultFile.add_var("fphys", ncDouble, dimtime
	//, dimNfield, dimK, dimNp);

	if (outputStep == StepPerFile + 1)
	{
		//resultFile.close();
	}
	else
	{
		outputStep++;
	}

};

void AbstractOutputFile::ncFile_create(int *Np, int *K, int Nvar)
{
	NcDim *dimNp = resultFile.add_dim("Np", *Np);
	NcDim *dimK = resultFile.add_dim("K", *K);
	NcDim *dimNfield = resultFile.add_dim("Nvar", Nvar);
	NcDim *dimtime = resultFile.add_dim("Nt");

	//std::vector<NcDim> dims_time;
	//dims_time.push_back(dimtime);

	//std::vector<NcDim> dims_fphys;
	//dims_fphys.push_back(dimtime);
	//dims_fphys.push_back(dimNfield);
	//dims_fphys.push_back(dimK);
	//dims_fphys.push_back(dimNp);

	output_time = resultFile.add_var("time", ncDouble, dimtime);
	output_time->add_att("units", "s");
	output_fphys = resultFile.add_var("fphys", ncDouble, dimtime, dimNfield, dimK, dimNp);

	//resultFile.enddef();
};