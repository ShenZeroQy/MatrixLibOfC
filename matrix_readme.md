��ӭʹ�þ������İ����ĵ�����
���ߣ������
��Ȩ(C)��Murphy���������У���ֹ�⴫��
Bug������shen99855@outlook.com
�汾��v1.02
	���ڣ�2021 04
	�������ͣ�
			
		Matrixd���ṹ�壬�ɱ��С����
			��Ա��
				row��������
				col��������
				data*������ָ�룩
				trans_flag(ת�ñ�־)
			ʹ�÷�����
				��Ҫ׼�򣺶���󣬱������Matrix_init��������ʼ����ʹ����ɺ󣬱������Matrix_release()�ͷſռ䡣
				�߼�:��ָͨ����ڴ������Ƶ��û����Ժ�������׼��
				Ԫ�ز��ң�
					����Matrixd����ʹ�� Matrix_Seek_Element(name,i,j) �������ҵ�name[i,j]
					����Matrixd*  ��ʹ��Matrix_Pointer_Seek_Element(name,i,j)  �������ҵ�name[i,j]
				Ԫ�ظ�ֵ��
					ע�⣬ֻ�ܶ�trans_flag=0�ľ����Ԫ�ؽ��и�ֵ��
					�����֪��trans_flag�Ƿ�Ϊ0�����Ե���һ��Matrix_clr_trans_flag����������������Խ���Ӧ�����trans_flag��0��
					��֤���󲻱�ת�õ�ǰ�������£�
						����Matrixd����ʹ��  Matrix_Assign_Notrans_Element(name,i,j) �����鶨λname[i,j]
						����Matrixd*  ��ʹ�� Matrix_Pointer_Assign_Notrans_Element(name,i,j)  �����鶨λname[i,j]
					ʾ����
						Matrixd* mat ,copy;
						...
						Matrix_Pointer_Assign_Notrans_Element(copy,i,j) = Matrix_Pointer_Seek_Element(mat, i, j);
		Vector���ṹ�壬�ɱ䳤������
			��Ա��
				len�����ȣ�
				data*������ָ�룩
			��Ϊ������ȽϾ�����ԣ��򵥵ö࣬
			����û�з�װ��Matrix_init()�������Ļ���������
			�û������Լ�����data*��ָ�򣬵�Ӧע���ڴ�й©���⡣
		Cal_state:ö�٣�����ָʾ����״̬
			�μ��䶨�壬�����������ķ���ֵ���塣
	������
		����
			Matrix�������г��ִ˴ʻ㣬��Ϊ��������
			MatrixS�������г��ִ˴ʻ㣬��Ϊ��������
		�����ຯ����
			�����ຯ��ʵ���˾�����ڴ���䣬�ָ�͹���ȣ������м����ຯ��ʵ�ֵĻ�ʯ��
			Cal_state Matrix_init( Matrixd* mat,int row, int col, int f)
				������������ʼ��matΪrow��col�У��������ڴ档
				������
					mat�������������ĵ�ַ
					row������
					col������
					f����־λ�����������¼���ѡ�
						'0':ȫ0��
						'1':ȫ1��
						'i'����λ��
						't'����0��ʼԪ���������1
						'd'��ֻ�����ڴ棬����Ԫ�ؽ��г�ʼ��ֵ�������Ĳ���

				����ֵ��
					success���ɹ�
				����˵����
					�μ�Matrixdʹ�õ���Ҫ׼��
			Cal_state Matrix_transpose(Matrixd* mat)
				���������������ת�á�
				������
					mat�������������ĵ�ַ
				����ֵ��
					success���ɹ�
				����˵����
					ֻ�ı��ڴ��ȡ˳�򣬲��ı��ڴ�ʵ������˳��
			void Matrix_clr_trans_flag(Matrixd* mat)	
				������������������trans_flag��
				������
					mat�������������ĵ�ַ
				����ֵ��
					
				����˵����
					�μ�Matrixdʹ�÷�����Ԫ�ظ�ֵһ�ڡ�
			Cal_state Matrix_reshape(Matrixd* mat, int row, int col)
				���������������С���·���Ϊrow��col�С�
				������
					mat�������������ĵ�ַ
					row������
					col������
				����ֵ��
					error��row��col�Ƿ�
					success:�ɹ�����ԭ���ݱ���
					optial���ڴ汻����ԭ���ݶ�ʧ,���Ǵ洢�ռ��С�ɹ��ı䡣
				����˵����
					����successʱ�������trans_flag=0��
			Cal_state Matrix_copy(Matrixd* mat, Matrixd* copy)
				���������������ƣ�copy=mat��������mat��
				������
					mat�������������ĵ�ַ
					copy����������ַ
				����ֵ��
					do_noing��ʹ����δ��ʼ����copy
					success:�ɹ�����copy��transflag=0��					
				����˵����
					����successʱ��copy�Ĵ洢�ռ��������С�
			void Matrix_show(Matrixd* mat)
				������������ӡ���󣬸��ݲ�ͬƽ̨���ж���GUI��ʾ��
				������
					mat�������������ĵ�ַ
				����ֵ��
					��					
				����˵����
					�ú����Ե��Էǳ��а������û������Լ�д��
					����Ƕ��ʽϵͳ���������ա�
			void Matrix_release(Matrixd* mat)
				�����������ͷž����ڴ档
				������
					mat�������������ĵ�ַ
				����ֵ��
					��					
				����˵����
					�μ�Matrixdʹ�õ���Ҫ׼��
		�����ຯ����
			�����ຯ������ʵ���Ծ���Ϊ��������ļ��㷽�����㷨�ȡ�
			Cal_state Matrix_add( Matrixd* mat, Matrixd* adder)
				��������������ӣ�mat=mat+adder
				������
					mat�������������ĵ�ַ
					adder��������
				����ֵ��
					error:adder�����Ϲ淶���������������Խ����ʡ�
					success:�ɹ�
				����˵����
					�������matԭλ��
					Ҫ��adder���л����������ڵ���mat������ǵ��ڣ�matlab�����Ҫ�󣩡�
			Cal_state Matrix_mul_num(Matrixd* mat, double num)
				�����������������ˣ�mat=mat*num
				������
					mat�������������ĵ�ַ
					num���˵�ϵ��
				����ֵ��
					success:�ɹ�
				����˵����
					�������matԭλ��
			Cal_state Matrix_mul_Matrtix(Matrixd* mat_left, Matrixd* mat_right,Matrixd* mat_result)	
				��������������˾���mat_result=mat_left*mat_right
				������
					mat_left�������ĵ�ַ
					mat_right���Ҿ���ĵ�ַ
					mat_result�����
				����ֵ��
					success:�ɹ�
					error��left->rol!=right->row
				����˵����
					result�����С����Ӧ��
			Cal_state Matrix_QRc(Matrixd* mat, Matrixd* Q, Matrixd* R);
				��������������QR�ֽ⣬mat=Q*R����֪������ٶ�
				������
					mat�������������ĵ�ַ
					Q:�������
					R:�������
					
				����ֵ��
					success:�ɹ�
					error����ȫ���е��·ֽ�ʧ��
					do_nothing:������mat������С�ڵ�������
				����˵����
					Ҫ��mat������С�ڵ�������
					û�����㷨�Ż���
			Cal_state MatrixS_LUc(Matrixd* mat, Matrixd* L, Matrixd* U);
				��������������LU�ֽ⣬mat=L*U����֪������ٶ�
				������
					mat�������������ĵ�ַ
					L:�������
					U:�������
					
				����ֵ��
					success:�ɹ�					
					do_nothing:�Ƿ���
				����˵����
					�㷨���Ӷ�1/3 *N^3
			Cal_state MatrixS_cholesky_c(Matrixd* mat, Matrixd* L);		
				��������������ʵ�Գƾ���������˹���ֽ⣬mat=L*trans(L)����֪������ٶ�
				������
					mat�������������ĵ�ַ
					L:�������					
				����ֵ��
					success:�ɹ�					
					do_nothing:�Ƿ���
				����˵����
					�����������ۻ�����������Ͻǲ�����Ϊ���Ϊ0���������½��ۼ������Դﵽ�����������100�����������ߡ�
			double MatrixS_det(Matrixd* mat)
				������������������ʽ
				������
					mat�������������ĵ�ַ��Ҫ���Ƿ���					
				����ֵ��
					��������ʽ					
					do_nothing:�Ƿ���
				����˵����
					�����������ʽ������
					//���ܸ������˵���������ѣ�ʹ��LUc��������ʽ��
			Cal_state MatrixS_inverse(Matrixd* mat, Matrixd* mat_inv)
				�����������������棬mat_inv=inv(mat)
				������
					mat�������������ĵ�ַ
					mat_inv��mat����
				����ֵ��
					do_nothing:���Ƿ���
					error:mat������
					success:�ɹ�
				����˵����
					������������
					�����㷨�ʹ洢�Ż������Է��Ĵ󵨵��ã�Ч������ߵġ�