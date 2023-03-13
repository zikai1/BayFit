#include "mex.h"
//#include "RDPM.h"
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>
//#include"opencv2/opencv.hpp"
//#include"matrix.h"
#include <iostream>
#include <cmath>
#include <utility>
//#include <vector>
//#include <stdexcept>
//#include<opencv2/opencv.hpp>

using namespace std;
//using namespace cv;


/*
double PerpendicularDistance(const Point& pt, const Point& lineStart, const Point& lineEnd)
{
	double dx = lineEnd.x - lineStart.x;
	double dy = lineEnd.y - lineStart.y;

	//Normalise
	double mag = pow(pow(dx, 2.0) + pow(dy, 2.0), 0.5);
	if (mag > 0.0)
	{
		dx /= mag; dy /= mag;
	}

	double pvx = pt.x - lineStart.x;
	double pvy = pt.y - lineStart.y;

	//Get dot product (project pv onto normalized direction)
	double pvdot = dx * pvx + dy * pvy;

	//Scale line direction vector
	double dsx = pvdot * dx;
	double dsy = pvdot * dy;

	//Subtract this from pv
	double ax = pvx - dsx;
	double ay = pvy - dsy;

	return pow(pow(ax, 2.0) + pow(ay, 2.0), 0.5);
}
*/

/*
void RamerDouglasPeucker(const vector<Point>& pointList, const float& epsilon, vector<Point>& out)// const vector<Point> &pointList
{
	if (pointList.size() < 2)
		throw invalid_argument("Not enough points to simplify");

	// Find the point with the maximum distance from line between start and end
	double dmax = 0.0;
	int index = 0;
	int end = pointList.size() - 1;
	for (int i = 1; i < end; ++i)
	{
		double d = PerpendicularDistance(pointList[i], pointList[0], pointList[end]);
		if (d > dmax)
		{
			index = i;
			dmax = d;
		}
	}

	// If max distance is greater than epsilon, recursively simplify
	if (dmax > epsilon)
	{
		// Recursive call
		vector<Point> recResults1;
		vector<Point> recResults2;
		vector<Point> firstLine(pointList.begin(), pointList.begin() + index + 1);
		vector<Point> lastLine(pointList.begin() + index, pointList.end());
		RamerDouglasPeucker(firstLine, epsilon, recResults1);
		RamerDouglasPeucker(lastLine, epsilon, recResults2);

		// Build the result list
		out.assign(recResults1.begin(), recResults1.end() - 1);
		out.insert(out.end(), recResults2.begin(), recResults2.end());
		if (out.size() < 2)
			throw runtime_error("Problem assembling output");
	}
	else
	{
		//Just return start and end points
		out.clear();
		out.push_back(pointList[0]);//pointList[0]
		out.push_back(pointList[end]);
	}
}
*/

// vector<vector<Point>> LineApproximate(vector<vector<Point>> edge_contour, const float& epsilon)
// {
// 	vector<vector<Point>> contour_seg;
//     for(int i=0;i<edge_contour.size();i++){
//         vector<Point> segTemp;
//         RamerDouglasPeucker(edge_contour[i], epsilon, segTemp);// distance threshold is 2.5
//         contour_seg[i]=segTemp;
//     }
// 	return contour_seg;
// 
// }

// vector<Point> LineApproximate(vector<Point> edge_contour, const float& epsilon)
// {
// 	//vector<vector<Point>> contour_seg;
//     //for(int i=0;i<edge_contour.size();i++){
//         vector<Point> segTemp;
//         RamerDouglasPeucker(edge_contour, epsilon, segTemp);// distance threshold is 2.5
//         //contour_seg[i]=segTemp;
//     //}
// 	//return contour_seg;
//       return segTemp;
// 
// }
//#define h 1e-1

/*计算核密度*/
double kernel_comp(double* data, vector<int> total_inx, int inx, int point_num)
{  
    
    
    int i, j, N=total_inx.size();// 点inx的总共邻近点数（包含自身）
    double sum=0,h=1e-1,razn,diff;
    double temp=1/((N+1)*pow(2*3.14159265358979*h*h,1.5));
    
    //mexPrintf("Temp value:%f\n",temp);
    /* 输出没有问题 done!
    mexPrintf("Sum of near point:%d,%f\n",N,temp);
    mexPrintf("inx: %d, point_num: %d, bandwidth: %f\n",inx,point_num,h);
    */
    for(i=0;i<N;i++)
    {   
        razn=0;
        // 两点之间距离
        for (j=0;j<3;j++)
        {
           diff=*(data+inx+j*point_num)-*(data+total_inx[i]+j*point_num);
           diff=diff*diff;
           razn+=diff;
        }
        sum=sum+exp(-0.5*razn/(h*h)); 
    }
    sum=sum+1;// 点自身计算的概率密度
    //mexPrintf("Sum of density:%f\n",sum);
    sum=temp*sum;
    //mexPrintf("Final sum of density:%f\n",sum);
    return sum;//概率密度
}


void score_comp(
        double* data,
        double* knn_map_id,
        double* score,
        int M,
        int K
        )
{
   
   int m,i,j,k,l;
   double pdf;//概率密度
   vector<double> point_pdf;//存储每个点的概率密度
   vector<vector<int>> all_total_inx;// 所有点周围的K近邻点的索引
   for(i=0;i<M;i++){// 构建所有点的K近邻
       vector<int> total_inx;// 记录某点的所有K邻近点索引
       j=i+1;
       
       //以第i个点为其K近邻的点的索引
       for(m=0;m<M;m++){//从第一行遍历到末行
           for(k=0;k<K;k++){// column：第一列到第K列
               if(*(knn_map_id+m+k*M)==j)// 判断索引是否为i（从1开始）
               {
                   //加入邻近点索引
                   total_inx.push_back(m);
               }
           }
       }
       
       //第i个点的K近邻索引
       for(k=0;k<K;k++)
       {
           l=*(knn_map_id+i+k*M)-1;// 下标从0开始 
           total_inx.push_back(l);
           // 对每个点继续它自己的K近邻索引
           int l1;
           for(k=0;k<K;k++)
           {
             l1=*(knn_map_id+l+k*M)-1; 
             total_inx.push_back(l1);  
           }
       }
       //删除vector中重复的元素
       //先排序
       sort(total_inx.begin(),total_inx.end());
       //将重复的数据移到后面
       vector<int>:: iterator ite=unique(total_inx.begin(),total_inx.end());
       total_inx.erase(ite,total_inx.end());
       all_total_inx.push_back(total_inx);//记录i的所有K近邻点的索引
       //mexPrintf("%d\n",total_inx.back());//应该没问题 done!
       //计算核密度
       pdf=kernel_comp(data,total_inx,i,M);
       //mexPrintf("Pdf of each point: %f\n",pdf);
       point_pdf.push_back(pdf);//存储每个点的概率密度
            
 }
   // 计算每个点和其K近邻点平均概率密度的比值
   
   for(i=0;i<M;i++)
   {
      
       int num=all_total_inx[i].size(),temp; 
       double ave,pdf_sum=0;
       for(int j=0;j<num;j++)
       {  
          temp=all_total_inx[i][j]-1;// 索引从0开始
          pdf_sum+=point_pdf[temp]; 
       }
       //pdf_sum=pdf_sum-point_pdf[i];//不连第i个点本身的pdf
       ave=pdf_sum/(num-1);
       //mexPrintf("%f,%f\n",ave,point_pdf[i]);
       *(score+i)=ave/point_pdf[i];
   }
  return;
}


// Input paramters
#define IN_x    prhs[0]// data
#define IN_y    prhs[1]//knn_map_id
#define IN_z    prhs[2]//K

//Output paramters
#define OUT_x    plhs[0]// score for each data point


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])        
{
    if(nrhs!=3)
    {
        mexErrMsgTxt("There are 3 input parameters (data (matrix), knn_map_id(matrix), K(int)) in KNN method.");
    }
    if(nlhs!=1)
    {
        mexErrMsgTxt("Output parameter is matrix(array) class pointer only.");
    }
    
    double *data,*knn_map_id,*score;
    int M,K;
    
    
    M=mxGetM(IN_x);// the input row number
    
    OUT_x=mxCreateDoubleMatrix(M,1,mxREAL);// the output score 
    
    // assign pointers to the input parameters
    data=mxGetPr(IN_x);
    knn_map_id=mxGetPr(IN_y);
    //K=mxGetPr(IN_z);
    K=mxGetScalar(IN_z);
    
    //assign pointers to the output parameters
    score=mxGetPr(OUT_x);
    
    /* Do the actual computations in the subroutine*/
    score_comp(data,knn_map_id,score,M,K);
    return;
 /*   
    // Get the cell data
    mwSize total_num_of_cells;
    int m,n;
    mxArray *pMatrix;
    double *pVal;
    //const mxArray* cell_array_ptr=mxGetPr(prhs[0]);
    
    total_num_of_cells=mxGetNumberOfElements(prhs[0]);// elements of the cell
    //mexPrintf("total number of cells: %d\n",total_num_of_cells);
    double epsilon = *(double *)mxGetPr(prhs[1]);// distance tolerance
    mxArray *cell_array_ptr=mxCreateCellMatrix(total_num_of_cells,1);// output cell
    //vector<vector<Point>> total_contour;
    //vector<vector<Point>> contour_result;
    
    for(int i=0;i<total_num_of_cells;i++)
    {
        pMatrix=mxGetCell(prhs[0],i);
        m=mxGetM(pMatrix);
        n=mxGetN(pMatrix);
        pVal=mxGetPr(pMatrix);// store by columns 1 array
        //mexPrintf("index of element: %d, M=%d, N=%d\n",i,m,n);
        
        vector<Point> temp_contour;
        temp_contour.resize(m);
        for(int j=0;j<m;j++){
            Point temp;
            temp.x=pVal[j];
            temp.y=pVal[j+m];// need to convert to 1 dimonsional
            //mexPrintf("00000");
            temp_contour[j]=temp;
            //mexPrintf("00000");
        }
        //mexPrintf("1111111");
        vector<Point> temp_contour_result;
        temp_contour_result=LineApproximate(temp_contour,epsilon);
        //total_contour[i]=temp_contour;
        int dominant_num=temp_contour_result.size();
        mxArray *contour_result_mat=mxCreateDoubleMatrix(dominant_num,2,mxREAL);
        double *result_ptr;
        result_ptr=mxGetPr(contour_result_mat);
        //mexPrintf("2222222");
        for(int k=0;k<dominant_num;k++){
            result_ptr[k]=temp_contour_result[k].x;
            result_ptr[k+dominant_num]=temp_contour_result[k].y;
        }
        //mexPrintf("333333");
        mxSetCell(cell_array_ptr,i,mxDuplicateArray(contour_result_mat));
        //mexPrintf("44444");
        mxDestroyArray(contour_result_mat);
    }
    plhs[0]=cell_array_ptr;
*/ 
//     contour_result=LineApproximate(total_contour,epsilon);
//     for(int k=0;k<total_num_of_cells;k++){
//         
//         
//     }
//     mxSetCell(cell_array_ptr,);
//     
//     
//     
//     vector<vector<Point>> contour;
//     mexArray *cellofmatrix;// mexArray pointer
//     ptr=mexGetCell(phrs[0],0);
//     
//     
//     
//     // Get the data;
//     double  = *(double *)mxGetPr(prhs[0]);
//     double dcols = *(double *)mxGetPr(prhs[1]);
//     
//     AAMED *_aamed = new AAMED(std::round(drows), std::round(dcols));
//     int pointer_size = sizeof(_aamed);
//     
//     plhs[0] = mxCreateNumericMatrix(1, pointer_size, mxUINT8_CLASS, mxREAL);
//     unsigned char *_pt = (unsigned char *)mxGetPr(plhs[0]);
//     unsigned char *_pt_aamed = (unsigned char *)&_aamed;
//     for(int i = 0; i < pointer_size; i++)
//     {
//         _pt[i] = _pt_aamed[i];
//     }

}