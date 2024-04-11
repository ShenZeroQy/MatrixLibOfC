//#include"matrix.h"
//#include"pch.h"
//
//Cal_state Matrix_init( Matrixd* mat,int row, int col, int f)
//{
//	//mat->com = NULL;
//	mat->trans_flag = 0;
//	mat->col = col;
//	mat->row = row;
//	mat->data = (float*)malloc(col*row * sizeof(float));
//	int i, j;
//	switch (f)
//	{
//	case '0'://ȫ0
//		for (i = 0; i < row; i++)
//			for (j = 0; j < col; j++)
//				*(mat->data + i * col + j) = 0;
//		break;
//	case '1'://ȫ1
//		for (i = 0; i < row; i++)
//			for (j = 0; j < col; j++)
//				*(mat->data + i * col + j) = 1;
//		break;
//	case'i'://��λ
//		for (i = 0; i < row; i++)
//			for (j = 0; j < col; j++)
//			{
//				if (i == j)
//					*(mat->data + i * col + j) = 1;
//				else
//					*(mat->data + i * col + j) = 0;
//			}
//		break;
//	case't'://test
//		for (i = 0; i < row; i++)
//			for (j = 0; j < col; j++)
//			{
//
//				*(mat->data + i * col + j) = i * col + j;
//
//			}
//		break;
//	default:
//		break;
//	}
//
//	return success;
//
//}
//Cal_state Matrix_transpose(Matrixd* mat)
//{
//	int tem = mat->col;
//	mat->col = mat->row;
//	mat->row = tem;
//	mat->trans_flag = !(mat->trans_flag);
//	return success;
//}
//void Matrix_add(Matrixd* mat, Matrixd* adder)//�������matԭλ
//{
//	int i, j;
//	//ֻ����adder��mat��
//	if (mat->col < adder->col || mat->row < adder->row)
//		return;
//	else
//	{
//		Matrix_clr_trans_flag(mat);
//		for (i = 0; i < mat->row; i++)
//			for (j = 0; j < mat->col; j++)
//				Matrix_Pointer_Assign_Notrans_Element(mat, i, j) = Matrix_Pointer_Assign_Notrans_Element(mat, i, j)+Matrix_Pointer_Seek_Element(adder, i, j);
//		return;
//	}
//
//}
//
//void Matrix_mul_num(Matrixd* mat, float num)
//{
//	int i, j;
//	for (i = 0; i < mat->row; i++)
//		for (j = 0; j < mat->col; j++)
//			*(mat->data + i * mat->col + j) *= num;
//	return;
//}
//void Matrix_mul_Matrtix(Matrixd* mat_left, Matrixd* mat_right, Matrixd* mat_result)//f="l/row"
//{
//	if (mat_left->col != mat_right->row)
//		return;
//	else
//	{
//		int i, j, k;
//		Matrix_reshape(mat_result, mat_left->row, mat_right->col);
//		for (i = 0; i < mat_result->row; i++)
//			for (j = 0; j < mat_result->col; j++)
//			{
//				float sum = 0;
//				for (k = 0; k < mat_left->col; k++)
//					sum += Matrix_Pointer_Seek_Element(mat_left, i, k)*Matrix_Pointer_Seek_Element(mat_right, k, j);
//				Matrix_Pointer_Assign_Notrans_Element(mat_result, i, j) = sum;
//
//			}
//		return;
//	}
//
//
//}
//void Matrix_copy(Matrixd* mat, Matrixd* copy)
//{
//	//�������ƣ�ת�ñ�־λ��0���洢�ռ��������š�
//	copy->col = mat->col;
//	copy->row = mat->row;
//	copy->trans_flag = 0;
//	copy->data = (float*)realloc(copy->data, mat->col*mat->row * sizeof(float));
//	int i, j;
//	for (i = 0; i < mat->row; i++)
//		for (j = 0; j < mat->col; j++)
//			*(copy->data + i * copy->col + j) = Matrix_Pointer_Seek_Element(mat, i, j);
//	return;
//}
//Cal_state Matrix_reshape(Matrixd* mat, int row, int col)//ֻ�в����·����ַ��ʱ�򣬲ű���ԭ����
//{
//	if (row<=0||col<=0)
//		return -1;
//	Cal_state s = 0;
//	if (mat->col*mat->row < row * col)
//	{
//		mat->data = (float*)realloc(mat->data, row*col * sizeof(float));//enlarge memory
//		s = 2;
//	}
//	else
//	{
//		Matrix_clr_trans_flag(mat);//relinearlize memory
//		s = 1;
//	}
//	mat->row = row;
//	mat->col = col;
//	return s;
//
//}
//void Matrix_release(Matrixd* mat)
//{
//	free(mat->data);
//	//free(mat);
//	return;
//}
//void Matrix_clr_trans_flag(Matrixd* mat)
//{
//	//Ϊ�����Ч��,������ִ��(αת��)�߼�,���ڷ����Դ洢��Ԫ��,���ܲ�������,��reshape,���ô˺������������ת��,�ٴν��洢�ռ����Ի�
//	if (mat->trans_flag)
//	{
//		float* tem = (float*)malloc(mat->col*mat->row * sizeof(float));
//		int i, j;
//		for (i = 0; i < mat->row; i++)
//			for (j = 0; j < mat->col; j++)
//				*(tem + i * mat->col + j) = Matrix_Pointer_Seek_Element(mat, i, j);
//		mat->trans_flag = 0;
//		for (i = 0; i < mat->row; i++)
//			for (j = 0; j < mat->col; j++)
//				*(mat->data + i * mat->col + j) = *(tem + i * mat->col + j);//����洢
//		free(tem);//Ϊ�˱������,������mat->data��ָ��.������Ч��.
//
//	}
//	return;
//}
//Cal_state MatrixS_inverse(Matrixd* mat, Matrixd* mat_inv)
//{
//	//matlab ��������ֵ�ֽ������㣬����α�棻
//	//eigen ����ļ򵥴ֱ�����ֱ����Ϊ��������졣
//	if (mat->col != mat->row)
//		//���Ƿ���
//		return do_noing;
//	Matrixd matl, matu;
//	Matrix_init(&matl, mat->row, mat->col, '0');
//	Matrix_init(&matu, mat->row, mat->col, '0');
//	Matrix_reshape(mat_inv, mat->row, mat->col);
//	MatrixS_LUc(mat, &matl, &matu);
//	//Matrix_show(&matl);
//	//Matrix_show(&matu);
//	//Matrix_mul_Matrtix(&matl, &matu, mat_inv);
//	//Matrix_show(mat_inv);
//	float det = 1;
//	Cal_state s= 0;
//	int i;
//	for (i = 0; i < mat->row;i++)	
//		det *= Matrix_Assign_Notrans_Element(matl, i, i)*Matrix_Assign_Notrans_Element(matu, i, i);
//	if (det == 0)
//
//		s= error;
//	else
//	{
//		Matrixd matl_inv, matu_inv;
//		Matrix_init(&matl_inv, mat->row, mat->col, 'd');
//		Matrix_init(&matu_inv, mat->row, mat->col, 'd');
//		Matrix_transpose(&matl);
//		MatrixS_Up_Triangle_inverse(&matu, &matu_inv);
//		MatrixS_Up_Triangle_inverse(&matl, &matl_inv);
//		//Matrix_show(&matl_inv);
//		//Matrix_show(&matu_inv);
//		Matrix_transpose(&matl_inv);
//		Matrix_mul_Matrtix(&matu_inv, &matl_inv, mat_inv);
//		Matrix_release(&matl_inv);
//		Matrix_release(&matu_inv);
//		s = success;
//	}
//	Matrix_release(&matl);
//	Matrix_release(&matu);
//	return s;
//	
//}
//float MatrixS_det(Matrixd* mat)//ʹ�ÿ������ֽⷨ����
//{	//�Ż�����������ʽ��Ϊ�ݹ����
//	//10��������ʱ��û������ģ��״ι��߻ᵼ���ڴ濪������
//	/*
//	�ڴ濪����״εĺ�����ϵk��n���������¼���
//	set k��3��=9��
//	set k��2��=4��
//	k��n��=k��n-1)*(n-2)+k(n-2);
//	k(n)<n!;
//	*/
//	//���㸴�ӶȽӽ�n�������Ƽ��ô˺����жϾ����Ƿ����ȡ�
//	float result;
//	//�����ת�ò��ı�����ʽ��ֵ
//
//	if (mat->col != mat->row)//���Ƿ���
//		return -1;
//	else
//	{
//		int d = mat->col;
//		if (d == 1)
//		{
//			return *(mat->data);
//		}
//		else
//		{
//			if (d == 2)
//			{
//				return *(mat->data)**(mat->data + 3) - *(mat->data + 1)**(mat->data + 2);
//			}
//			else
//			{
//				if (d == 3)
//				{
//					float zhu = *(mat->data)**(mat->data + 4)**(mat->data + 8) +
//						*(mat->data + 1)**(mat->data + 5)**(mat->data + 6) +
//						*(mat->data + 2)**(mat->data + 3)**(mat->data + 7);
//					float fu = *(mat->data)**(mat->data + 5)**(mat->data + 7) +
//						*(mat->data + 1)**(mat->data + 3)**(mat->data + 8) +
//						*(mat->data + 2)**(mat->data + 4)**(mat->data + 6);
//					return zhu - fu;
//				}
//				else
//				{
//					result = 0;//��ͱ�������
//					Matrixd* iter;
//					Matrixd* Rsv;
//					Rsv = (Matrixd*)malloc(d * sizeof(Matrixd));//Ϊ�����󿪱ٿռ�
//					int i, j, k;
//					for (i = 0; i < d; i++)
//					{
//						Matrix_init(d - 1, d - 1, Rsv + i, 't');//��ʼ�����а�����
//						//Matrix_show(Rsv + i);
//					}
//					for (k = 0; k < d; k++)//����d��������
//					{
//						iter = Rsv + k;//ָ���k��������
//						int jump_line = k;
//
//						for (i = 0; i < d - 1; i++)
//						{
//							int line = -1;
//
//							for (j = 0; j < d - 1; j++)//��������Ԫ��
//							{
//								line++;
//								if (line == jump_line)
//									line += 1;//����һ��
//								*(iter->data + i * iter->col + j) = *(mat->data + i * mat->col + line);
//								//*(((Rsv + k)->data) + i * ((Rsv + k)->col) + j) = 1;//debug ��
//							}
//							//���������Ժ�Ҫ�����һ����Ϊ�������ϵ��
//
//						}
//						//Matrix_show(Rsv + k);//for debug
//
//					}
//					for (k = 0; k < d; k++)//����d������ʽ,ȡ���һ�зֽ⡣
//					{
//						result += (*(mat->data + d * d - d + k))*MatrixS_det(Rsv + k)*((k % 2) ? -1 : 1);//d=mat->col
//						//Matrix_release(Rsv + k);//�������������ڴ�
//						//printf_s("bolck=%x\n", (Rsv + k)->data);
//					}
//					for (k = 0; k < d; k++)//����d������ʽ,ȡ���һ�зֽ⡣
//					{
//						//Matrix_release(Rsv + k);//�������������ڴ�
//						free((Rsv + k)->data);
//						//���ﲻ���Ե���Matrix_release()���������ڴ棬��Ϊ�����Ժ�Rsv���������Ż����ˣ���͵���Ѱַʧ�ܡ�
//						//����ֻ��main�����е���Matrix_release()��������װ����ֱ��free��Ӧָ�롣
//					}
//					//printf_s("result=%f\n,order=%d", result, d);
//					free(Rsv);
//					return result;
//				}
//			}
//
//		}
//
//	}
//
//}
//void Matrix_show(Matrixd* mat)
//{
//#ifdef GUI
//	int i, j;
//	float num;
//	for (i = 0; i < mat->row; i++)
//	{
//		for (j = 0; j < mat->col; j++)
//		{
//			num = Matrix_Pointer_Seek_Element(mat, i, j);
//			printf("%4.4f ", num);
//
//		}
//		printf("\n");
//
//	}
//#endif
//	return;
//}
//char Matrix_QRc(Matrixd* mat, Matrixd* Q, Matrixd* R)
//{
//	//��������Դ���磬�ݲ��ܱ�֤Ч�ʺ��ȶ��ԡ�
//	//see:https://blog.csdn.net/baidu_41647951/article/details/84401504
//	//https ://blog.csdn.net/baidu_41647951/article/details/84401504
//	if (mat->col > mat->row)
//		return -1;//��֤����С������
//	else
//	{
//		Matrix_reshape(Q, mat->row, mat->row);
//
//		Matrix_copy(mat, R);
//		int m = mat->row;
//		int n = mat->col;
//		int i, j, k, nn, jj;
//		nn = n;
//		if (m == n)
//		{
//			nn = m - 1;
//		}
//		float u, alpha, w, t;
//		for (k = 0; k <= nn - 1; k++)//�ڴ�ѭ��k��0~m���У�����H�������⣬���Q���Լ����A
//		{
//
//			u = 0.0;
//			for (i = k; i <= m - 1; i++)
//			{
//				w = fabs(Matrix_Pointer_Assign_Notrans_Element(R, i, k));
//				if (w > u)
//					u = w;
//			}
//			alpha = 0.0;
//			for (i = k; i <= m - 1; i++)
//			{
//				t = Matrix_Pointer_Assign_Notrans_Element(R, i, k) / u;
//				alpha = alpha + t * t;
//			}
//			if (Matrix_Pointer_Assign_Notrans_Element(R, k, k) > 0.0)
//				u = -u;
//			alpha = u * sqrt(alpha);
//			if (fabs(alpha) + 1.0 == 1.0)
//			{
//				return 0;//�ֽ�ʧ��
//				//��ȫ���е��·ֽ�ʧ��
//
//			}
//
//			u = sqrt(2.0*alpha*(alpha - Matrix_Pointer_Assign_Notrans_Element(R, k, k)));
//			if ((u + 1.0) != 1.0)
//			{
//				Matrix_Pointer_Assign_Notrans_Element(R, k, k) = (Matrix_Pointer_Assign_Notrans_Element(R, k, k) - alpha) / u;
//				for (i = k + 1; i <= m - 1; i++)
//					Matrix_Pointer_Assign_Notrans_Element(R, i, k) /= u;
//
//				//���Ͼ���H�������ã�ʵ���ϳ���û�������κ����ݽṹ���洢H��
//				//�󣬶���ֱ�ӽ�u������Ԫ�ظ�ֵ��ԭA�����ԭ��������Ӧ��λ�ã�������
//				//��������Ϊ�˼�����˾���Q��A
//				for (j = 0; j <= m - 1; j++)
//				{
//					t = 0.0;
//					for (jj = k; jj <= m - 1; jj++)
//						t = t + Matrix_Pointer_Assign_Notrans_Element(R, jj, k) * Matrix_Pointer_Assign_Notrans_Element(Q, jj, j);
//					for (i = k; i <= m - 1; i++)
//						Matrix_Pointer_Assign_Notrans_Element(Q, i, j) -= 2.0*t*Matrix_Pointer_Assign_Notrans_Element(R, i, k);
//				}
//				//��˾���Q��ѭ��������õ�һ�������ٽ��������ת��һ�¾͵õ�QR�ֽ��е�Q����
//				//Ҳ������������
//
//				for (j = k + 1; j <= n - 1; j++)
//				{
//					t = 0.0;
//					for (jj = k; jj <= m - 1; jj++)
//						t = t + Matrix_Pointer_Assign_Notrans_Element(R, jj, k) * Matrix_Pointer_Assign_Notrans_Element(R, jj, j);
//					for (i = k; i <= m - 1; i++)
//						Matrix_Pointer_Assign_Notrans_Element(R, i, j) -= 2.0*t*Matrix_Pointer_Assign_Notrans_Element(R, i, k);
//				}
//				//H�������A����ѭ�����֮���������ǲ��ֵ����ݾ��������Ǿ���R
//				Matrix_Pointer_Assign_Notrans_Element(R, k, k) = alpha;
//				for (i = k + 1; i <= m - 1; i++)
//					Matrix_Pointer_Assign_Notrans_Element(R, i, k) = 0.0;
//			}
//		}
//		Matrix_transpose(Q);
//		Matrix_clr_trans_flag(Q);
//
//		//QR�ֽ����
//		//Matrix_show(Q);
//		//Matrix_show(R);
//		return 1;
//
//	}
//
//}
//char MatrixS_LUc(Matrixd* mat, Matrixd* L, Matrixd* U)
//{
//	//L ������ U ������ mat=L*U
//	//ʹ��ֱ�����Ƿֽⷨ���㷨���Ӷ�n^3
//	//https://blog.csdn.net/weixin_44116061/article/details/105628206
////https://baike.baidu.com/item/lu%E5%88%86%E8%A7%A3/764245?fr=aladdin
//	if (mat->col != mat->row)
//		//���Ƿ���
//		return 0;
//	else
//	{
//		int N = mat->col;
//		//Matrix_clr_trans_flag(mat);
//		Matrix_reshape(L, N, N);
//		Matrix_reshape(U, N, N);
//		int i, j;
//		for (i = 0; i < N - 1; i++)
//		{
//			for (j = i + 1; j < N; j++)//��Ч��
//			{
//				Matrix_Pointer_Assign_Notrans_Element(L, i, j) = 0;
//				Matrix_Pointer_Assign_Notrans_Element(U, j, i) = 0;//LU��ʼΪ0
//			}
//		}
//		//for (int i = 0; i < N ; i++)//��Ч��
//		//{
//		//	for (int j = 0; j < N; j++)
//		//	{
//		//		Matrix_Pointer_Assign_Notrans_Element(L, j, i) = 0;
//		//		Matrix_Pointer_Assign_Notrans_Element(U, i, j) = 0;//LU��ʼΪ0
//		//	}
//		//}
//		for (i = 0; i < N; i++)
//		{
//			Matrix_Pointer_Assign_Notrans_Element(U, 0, i) = Matrix_Pointer_Seek_Element(mat, 0, i);
//			Matrix_Pointer_Assign_Notrans_Element(L, i, 0) = Matrix_Pointer_Seek_Element(mat, i, 0)/Matrix_Pointer_Seek_Element(U,0,0);
//		}
//		int r, k;
//		for (r = 1; r < N; r++)
//		{
//			for (i = r; i < N; i++)
//			{
//				float sumLrkUki_r_i;
//				float re = 0.0;
//				for (k = 0; k < r; k++)
//				{
//					re += Matrix_Pointer_Assign_Notrans_Element(L, r, k)*Matrix_Pointer_Assign_Notrans_Element(U, k, i);
//				}
//				sumLrkUki_r_i = re;
//				Matrix_Pointer_Assign_Notrans_Element(U, r, i) = Matrix_Pointer_Seek_Element(mat, r, i) - sumLrkUki_r_i;
//				if (i == r)
//					Matrix_Pointer_Assign_Notrans_Element(L, r, r) = 1;
//				else
//					if (r == N)
//						Matrix_Pointer_Assign_Notrans_Element(L, N, N) = 1;
//					else
//					{
//						float sumLikUkr_r_i;
//						re = 0.0;
//						for (k = 0; k < r; k++)
//						{
//							re += Matrix_Pointer_Assign_Notrans_Element(L, i, k)*Matrix_Pointer_Assign_Notrans_Element(U, k, r);
//						}
//						sumLikUkr_r_i = re;
//						Matrix_Pointer_Assign_Notrans_Element(L, i, r) = (Matrix_Pointer_Seek_Element(mat, i, r) - sumLikUkr_r_i) / Matrix_Pointer_Assign_Notrans_Element(U, r, r);
//
//					}
//			}
//		}
//		return 1;
//
//	}
//}
//void MatrixS_Up_Triangle_inverse(Matrixd* mat, Matrixd* inv)
//{
//	//Ҫ��mat�������Ƿ��������Խ�Ԫ����Ϊ0,
//	//����mat->trans_flag=1
//	//�������Ǿ�������inverse��
//	//��װΪMatrixS_inverse���Ӻ���
//	//�㷨���Ӷ�n^2
//	int N;
//	if (mat->col != mat->row)
//		//���Ƿ���
//		return;
//	else
//		N = mat->col;
//	Matrix_reshape(inv, N, N);
//	Vector One_col;
//	One_col.len = N;
//	One_col.data = (float*)malloc(N * sizeof(float));
//	Vector solu_col;
//	solu_col.len = N;
//	solu_col.data = (float*)malloc(N * sizeof(float));
//	int i, j;
//	for (j = N - 1; j >= 0; j--)//every col
//	{
//		for (i = 0; i < N; i++)
//			if (i == j)
//				*(One_col.data + i) = 1;
//			else
//			{
//				*(One_col.data + i) = 0;
//				if (i > j)
//					Matrix_Pointer_Assign_Notrans_Element(inv, i, j) = 0;
//			}
//		//one col done
//		for (i = j; i >= 0; i--)//�������ϵ�������mat�Ĵ�СΪj*j
//		{
//			float back;
//			back = *(solu_col.data + i) = *(One_col.data + i) / Matrix_Pointer_Seek_Element(mat, i, i);
//			//////////////////////////////////////////////////////////////////////
//			int ii;
//			for (ii = i - 1; ii >= 0; ii--)//call back without change mat
//			{
//				back = *(solu_col.data + i)*Matrix_Pointer_Seek_Element(mat, ii, i);
//				*(One_col.data + ii) -= back;
//			}
//		}
//		for (i = 0; i <= j; i++)
//		{
//			//copy solution vector to inv
//			Matrix_Pointer_Assign_Notrans_Element(inv, i, j) = *(solu_col.data + i);
//		}
//
//		//one col done
//	}
//	// all col done 
//	free(solu_col.data);
//	free(One_col.data);
//	//printf("suucceed computing tri inverse");
//	return;
//
//
//
//
//}
//void MatrixS_cholesky_c(Matrixd* mat, Matrixd* L)
//{
//	//�����������ۻ�����������Ͻǲ�����Ϊ���Ϊ0���������½��ۼ������Դﵽ�����������100�����������ߡ�
//	//��������������Ϊ��ȷ��С������λ����
//	if (mat->col != mat->row)
//		//���Ƿ���
//		return;
//	else
//	{
//		Matrix_reshape(L, mat->row, mat->col);
//		int i, j, k, p;
//		for (k = 0; k < mat->col; k++)
//		{
//			for (i = 0; i < k; i++)
//				Matrix_Pointer_Assign_Notrans_Element(L, i, k) = 0;//����������
//			float sum;
//			sum = 0;
//			for (p = 0; p < k; p++)
//				sum += Matrix_Pointer_Assign_Notrans_Element(L, k, p)*Matrix_Pointer_Assign_Notrans_Element(L, k, p);
//			sum = Matrix_Pointer_Assign_Notrans_Element(mat, k, k) - sum;
//			if (sum < 0)//��������˵��sumһ��Ϊ���������������������Ӧ����Ը�ֵ��ƽ��
//				sum = 0.000001;
//			Matrix_Pointer_Assign_Notrans_Element(L, k, k) = sqrt(sum);
//			//L(k,k) done
//
//
//			for (i = k + 1; i < mat->row; i++)
//			{
//				sum = 0;
//				for (j = 0; j < k; j++)
//				{
//					sum += Matrix_Pointer_Assign_Notrans_Element(L, i, j)*Matrix_Pointer_Assign_Notrans_Element(L, k, j);
//				}
//				sum = Matrix_Pointer_Assign_Notrans_Element(mat, i, k) - sum;
//				if (Matrix_Pointer_Assign_Notrans_Element(L, k, k) == 0)//����0��Ϊ��ĸ
//					Matrix_Pointer_Assign_Notrans_Element(L, k, k) = 0.0001;
//				Matrix_Pointer_Assign_Notrans_Element(L, i, k) = sum / Matrix_Pointer_Assign_Notrans_Element(L, k, k);
//			}
//			//L(i,k) done
//
//		}
//
//
//		return;
//	}
//
//
//}
//void Matrix_assign_by_Vec(Matrixd* mat, Vector* row_base, int row)
//{
//	//һ��һ�д���
//	Matrix_reshape(mat, row, row_base->len);
//	int i, j;
//	for (i = 0; i < row; i++)
//	{
//		for (j = 0; j < row_base->len; j++)
//		{
//			Matrix_Pointer_Assign_Notrans_Element(mat, i, j) = *((row_base + i)->data + j);
//		}
//	}
//}
//void Matrix_Splic(Matrixd* mat, Matrixd* under_mat)
//{
//	//|mat|=  |mat      |
//	//|   |   |under_mat|
//	//Ҫ��under_mat���������ڵ���mat
//
//	if (under_mat->col < mat->col)
//		return;
//	else
//	{
//		Matrixd copy;
//		Matrix_init(mat->row + under_mat->row, mat->col, &copy, 'd');
//		//Matrix_clr_trans_flag(mat);
//		//Matrix_copy(mat, &copy);
//		//Matrix
//		int i, j;
//		for (i = 0; i < mat->row; i++)
//		{
//			for (j = 0; j < mat->col; j++)
//				Matrix_Assign_Notrans_Element(copy, i, j) = Matrix_Pointer_Seek_Element(mat, i, j);
//		}
//		for (i = mat->row; i < mat->row + under_mat->row; i++)
//		{
//			int ind = i - mat->row;
//			for (j = 0; j < mat->col; j++)
//				Matrix_Assign_Notrans_Element(copy, i, j) = Matrix_Pointer_Seek_Element(under_mat, ind, j);
//		}
//		Matrix_copy(&copy, mat);
//		Matrix_release(&copy);
//		return;
//	}
//
//}
//void Matrix_get_col(Matrixd* mat, Vector* col_vec, int col)
//{
//	//Ĭ��col_vec �Ѿ������˴洢
//	col_vec->len = mat->row;
//	col_vec->data = (float*)realloc(mat->row, sizeof(float));
//	int i;
//	for (i = 0; i < mat->row; i++)
//	{
//		*(col_vec->data + i) = Matrix_Pointer_Seek_Element(mat, i, col);
//	}
//	return;
//
//}
//
//
