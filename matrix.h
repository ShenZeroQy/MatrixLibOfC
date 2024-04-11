#pragma once
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef _Matrix_H
#define _Matrix_H
#endif // !_Matrix_H

/****************************************************
Aurther:QYShen
Email:shen99855@outlook.com
Version:1.0.2
Description:Matrixd operation and storage management
Last_Edited_Data:14/04/2021
Attention:Everytime after defining a Matrixd, it must be inited,so you need to call Matrix_init();
		For every Matrixd  defined,if you no longer need it,you must release its storage space,so you need to call Matrix_release()
Please see more at:《matrix_readme.md》
*****************************************************/
typedef enum
{
	critial_error = -2,error=-1, do_noing = 0, success = 1, optial = 2, op_optial = 3
}Cal_state;
//#define Matrix_dim 5
#define Matrix_Seek_Element(name,i,j) (name.trans_flag?*(name.data+j*name.row+i):*(name.data+i*name.col+j))//世界上最快
#define Matrix_Pointer_Seek_Element(name,i,j) (name->trans_flag?*(name->data+j*name->row+i):*(name->data+i*name->col+j))//世界上最快
#define Matrix_Pointer_Assign_Notrans_Element(name,i,j) (*(name->data+i*(name->col)+j))
#define Matrix_Assign_Notrans_Element(name, i, j) (*(name.data + i * name.col + j))
typedef struct
{
	int row;//行
	int col;//列
	
	float* data;
	//Storage* com;//带内存分配的高级操作
	int trans_flag;
	/*
	if trans_flag=0
		mat(i,j)=*(mat->data+i*col+j)
	if rans_flag=1;
		mat(i,j)=*(mat->data+j*col+i)
		#define Mat(name,i,j) (transflag?name(i,j)=*(mat->data+i*col+j):name(i,j)=*(mat->data+j*col+i))
		
		*/
}
Matrixd;

//typedef union
//{
//	Vector* vecrow;
//	Vector* veccol;
//	float* data;
//}Storage;

typedef struct
{
	 int len;
	float* data;
}Vector;

Cal_state Matrix_init( Matrixd* mat,int row, int col, int f);
Cal_state Matrix_add( Matrixd* mat, Matrixd* adder);
Cal_state Matrix_transpose(Matrixd* mat);
Cal_state Matrix_mul_num(Matrixd* mat, double num);
Cal_state Matrix_mul_Matrtix(Matrixd* mat_left, Matrixd* mat_right,Matrixd* mat_result);
Cal_state Matrix_copy(Matrixd* mat, Matrixd* copy);
Cal_state Matrix_reshape(Matrixd* mat, int row, int col);
void Matrix_release(Matrixd* mat);
void Matrix_clr_trans_flag(Matrixd* mat);
Cal_state MatrixS_inverse(Matrixd* mat, Matrixd* mat_inv);
void Matrix_show(Matrixd* mat);
double MatrixS_det(Matrixd* mat);
Cal_state Matrix_QRc(Matrixd* mat, Matrixd* Q, Matrixd* R);
Cal_state MatrixS_LUc(Matrixd* mat, Matrixd* L, Matrixd* U);
Cal_state MatrixS_cholesky_c(Matrixd* mat, Matrixd* L);
void MatrixS_Up_Triangle_inverse(Matrixd* mat, Matrixd* inv);
