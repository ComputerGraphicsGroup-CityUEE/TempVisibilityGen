#include "assis_func.h"

void save_MyVertexArr(const char *spath, std::vector<MyVertex> &vbuf, int nv)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nv; i++)
	{
		fprintf(f0, "%f\t%f\t%f\t\n",vbuf[i].pos.x,vbuf[i].pos.y,vbuf[i].pos.z);
		fprintf(f0, "%f\t%f\t%f\t\n",vbuf[i].normal.x,vbuf[i].normal.y,vbuf[i].normal.z);
	}

	fclose(f0);
}

void save_intArr(char *spath, std::vector<int> &ibuf, int nbuf, int npl)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nbuf; i++)
	{
		fprintf(f0, "%d\t",ibuf[i]);
		if((i+1)/npl==0)
			fprintf(f0, "\n");
	}

	fclose(f0);
}
void save_uintArr(char *spath, std::vector<uint32_t> &ibuf, int nbuf, int npl)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<nbuf; i++)
	{
		fprintf(f0, "%d\t",ibuf[i]);
		if((i+1)/npl==0)
			fprintf(f0, "\n");
	}

	fclose(f0);
}
void save_vec3Arr(char *spath, std::vector<glm::vec3> &arr, int narr)
{
	FILE *f0 = fopen(spath, "wt");

	int i;
	for(i=0; i<narr; i++)
		fprintf(f0, "%f\t%f\t%f\t\n",arr[i].x,arr[i].y,arr[i].z);

	fclose(f0);
}
