#include<stdio.h>
#include<math.h>
//#include<stdbool.h>    ����VC++��������ͷ�ļ��ʲ�ʹ��

//n��ϵ������ʹ��ȫ��define��Ϊ�˱���malloc��δ�����̬���鵼���ڴ��ͻ
#define N 2 

int main() {
	double Simplest(double coefficientMatrix[N][N]);    //������Ϊ�����Ǿ��󲢼�������ʽ��ֵ
	double getDeterminant(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1], int i);     //���ݿ���Ĭ���򣬶�ϵ����������滻	                                                                                                  
	void showMatrix(double coefficientMatrix[N][N]);    //�������
	int checkType(double nonhomogeneousTerm[N][1]);    //ͨ��ϵ����������ʽֵ������Լ������������ж���η�����Լ�������
	int solutionsInfo(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1]);    //��Է�������Է���������������
	int augmentedMatrixSimplest(double augmentedMatrix[N][N + 1]);    //�������������ͣ���Ϊ���ǵ�����ʹ��malloc�����ڴ��ͻ
	                                                                  //ʹ��define���������С�������޷�����Simplest������ֻ�ܽ�����д 
	void getHomogeneousSolutions(double augmentedMatrix[N][N], int solutionsNum);    //����������Է�����ͨ�⣬δ��ʵ��

	double coefficientMatrix[N][N];
	printf("������ϵ������: ");
	for (int i = 0; i < N; i++)    //ϵ�����󴴽�
		for (int j = 0; j < N; j++)
			scanf("%lf", &coefficientMatrix[i][j]);
	double nonhomogeneousTerm[N][1];
	printf("������������: ");
	for (int i = 0; i < N; i++)    //��������
		scanf("%lf", &nonhomogeneousTerm[i][0]);

	if (checkType(nonhomogeneousTerm) == 0) {
		printf("������Ƿ�������Է�����\n");
		double determinant[N + 1];    //����Ĭ�������������ʽֵ����
		for (int i = -1; i < N; i++) {
			determinant[(i + 1)] = getDeterminant(coefficientMatrix, nonhomogeneousTerm, i);
		}

		if (determinant[0] != 0) {
			printf("�÷�������Է��������Ψһ��\n");

			//for (int i = 0; i < N + 1; i++)
				//printf("%-7.2f", determinant[i]);

			//printf("\n");
			for (int i = 1; i < N + 1; i++) {
				printf("��%dΪ%.2f", i, (determinant[i] / determinant[0]));    //����Ĭ��������
				printf("\n");
			}
		}
		else {
			int type = solutionsInfo(coefficientMatrix, nonhomogeneousTerm);
			if (type == -1) {
				printf("�����쳣��");
			}

			//�����ߴ��鱾�����Է�����������Ľ���������޷���Ч��Ƴ���ӦC���Դ���ʵ����δ֪����ʾ����δ֪��
			//֮��д��������ʽ�ó����������ʽ��ڴ��жϽ�����

		}
	}
	else if (checkType(nonhomogeneousTerm) == 1) {
		printf("�������������Է�����\n");
		double determinant;    //����ϵ����������ʽ
		determinant = Simplest(coefficientMatrix);

		if (determinant != 0) {
			printf("��������Է�����������\n");
		}
		else {
			printf("��������Է�������ڷ����\n");
			int r = N;
			int j = 0;
			for (int i = N - 1; i >= 0; i--) {    //ͨ����ȡϵ������������жϽ�����
				for (; j < N; j++)
					if (coefficientMatrix[i][j] != 0)
						break;
				if (j == N) 
					r--;	
				else
					break;
				j = 0;
			}
			if (r != 0)
				printf("��������������%d��\n", (N - r));
			else
				printf("ϵ�����������壡");   //ϵ������ȫΪ0

			//for (int i = 1; i < N + 1; i++) {
				//printf("��%dΪ%.2f", i, (determinant[i] / determinant[0]));
				//printf("\n");
		//	}
		}
	}
	else {
		printf("�����쳣��");
	}

	
	//system("pause");
	int pause;
	scanf("%d", &pause);    //���system("pause");
	
	return 0;
}

double Simplest(double coefficientMatrix[N][N]) {

	int point = 0;    //ָ���һ��������
	double determinant = 1;    //��������ʽ��ֵ
	int flag = 1;    //�ж��ڵ�һ����λΪ�������£��Ƿ�����б任�����Ƿ��һ�о�Ϊ��
	for (int i = 0; i < N; i++) {
		if (coefficientMatrix[i][point] == 0) {
			for (int j = i + 1; j < N; j++) {    //�ڵ�һ�е�һ��Ϊ�������£�ȥ���������Ƿ������λ������
				if (coefficientMatrix[j][point] == 0)
					continue;
				else {    //�����ڣ������б任������־λ��Ϊfalse�����к�������
					for (int n = point; n < N; n++) {
						double temp = coefficientMatrix[j][n];
						coefficientMatrix[j][n] = coefficientMatrix[i][n];
						coefficientMatrix[i][n] = temp;
					}
					flag = 0;
					determinant = determinant * (-1);
				}
			}
			if (flag == 1) {    //����һ�о�Ϊ0����point���Ƶڶ��У���i--������ѭ���ظ���������
				point++;
				i--;
				continue;
			}
		}

		//showMatrix(coefficientMatrix);

		for (int j = i + 1; j < N; j++) {
			double multiple;    //��һ������������һ�е�һ���������ı���
			if (coefficientMatrix[j][point] == 0) {    //������������ƽ�����ֵ����
				for (int l = j + 1; l < N; l++) {
					if (coefficientMatrix[l][point] == 0)
						continue;
					else {
						for (int n = point; n < N; n++) {
							double temp = coefficientMatrix[l][n];
							coefficientMatrix[l][n] = coefficientMatrix[j][n];
							coefficientMatrix[j][n] = temp;
						}
						flag = 0;
						determinant = determinant * (-1);
					}
				}
				if (flag == 1) {
					break;
				}
			}
			else {
				if (point >= N)    //����point������������Խ��Ŀ��ܣ����ڴ˽����ж�
					break;
				multiple = coefficientMatrix[j][point] / coefficientMatrix[i][point];    //��ȡ��֮�䱶��
				coefficientMatrix[j][point] = 0;    //��������λ��0
				for (int k = point + 1; k < N; k++) {
					coefficientMatrix[j][k] = coefficientMatrix[j][k] - coefficientMatrix[i][k] * multiple;    //�����б任
				}
			}			
		}
		//showMatrix(coefficientMatrix);
		flag = 1;    //��ʼ��flag
		point++;    //һ��ȫ���б����仯������point����
	}

	showMatrix(coefficientMatrix);
	//printf("\n");
	
	for (int i = 0; i < N; i++) 
		determinant = determinant * coefficientMatrix[i][i];    //�����Ǿ�������ʽ����
	return determinant;
}

double getDeterminant(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1], int i) {
	double Matrix[N][N];    //ϵ������ֵ������ֱ�Ӳ���ϵ�������³�ʼϵ���������
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Matrix[i][j] = coefficientMatrix[i][j];
	double determinant = 0;
	//showMatrix(Matrix);
	if (i == -1) {    //-1��ʾΪ��ʼϵ����������ʽ��ֵ
		printf("ϵ����������Σ�\n");
		determinant = Simplest(Matrix);
		//determinant = 1;
	}
	else {    //�ӵ�һ�п�ʼ�滻��������ʽֵ
		for (int j = 0; j < N; j++)
			Matrix[j][i] = nonhomogeneousTerm[j][0];
		//showMatrix(Matrix);
		//printf("\n");
		printf("��%d���滻���������Σ�\n", (i + 1));
		determinant = Simplest(Matrix);
	}
	return determinant;
}

void showMatrix(double coefficientMatrix[N][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			printf("%-7.2f", coefficientMatrix[i][j]);
		printf("\n");
	}
	printf("\n");
}

int checkType(double nonhomogeneousTerm[N][1]) {
	int type = -1;
	for (int i = 0; i < N; i++) {
		if (nonhomogeneousTerm[i][0] != 0) {
			type = 0;    //��������Է�����
			break;
		}
		else {
			type = 1;    //������Է�����
		}
	}
	return type;
}

int solutionsInfo(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1]) {
	double augmentedMatrix[N][N + 1];    //�����������
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			augmentedMatrix[i][j] = coefficientMatrix[i][j];
	for (int i = 0; i < N; i++)
		augmentedMatrix[i][N] = nonhomogeneousTerm[i][0];	

	return augmentedMatrixSimplest(augmentedMatrix);
}

int augmentedMatrixSimplest(double augmentedMatrix[N][N + 1]) {
	int point = 0;    //ָ���һ��������

	int flag = 1;    //�ж��ڵ�һ����λΪ�������£��Ƿ�����б任�����Ƿ��һ�о�Ϊ��
	for (int i = 0; i < N; i++) {
		if (augmentedMatrix[i][point] == 0) {
			for (int j = i + 1; j < N; j++) {    //�ڵ�һ�е�һ��Ϊ�������£�ȥ���������Ƿ������λ������
				if (augmentedMatrix[j][point] == 0)
					continue;
				else {    //�����ڣ������б任������־λ��λfalse�����к�������
					for (int n = point; n < N + 1; n++) {
						double temp = augmentedMatrix[j][n];
						augmentedMatrix[j][n] = augmentedMatrix[i][n];
						augmentedMatrix[i][n] = temp;
					}
					flag = 0;
				}
			}
			if (flag == 1) {    //����һ�о�Ϊ0����point���Ƶڶ��У���i--������ѭ���ظ���������
				point++;
				i--;
				continue;
			}
		}

		//showMatrix(coefficientMatrix);

		for (int j = i + 1; j < N; j++) {
			double multiple;    //��һ������������һ�е�һ���������ı���
			if (augmentedMatrix[j][point] == 0) {    //������������ƽ�����ֵ����
				for (int l = j + 1; l < N; l++) {
					if (augmentedMatrix[l][point] == 0)
						continue;
					else {
						for (int n = point; n < N + 1; n++) {
							double temp = augmentedMatrix[l][n];
							augmentedMatrix[l][n] = augmentedMatrix[j][n];
							augmentedMatrix[j][n] = temp;
						}
						flag = 0;
					}
				}
				if (flag == 1) {
					break;
				}
			}
			else {
				if (point >= N)
					break;
				multiple = augmentedMatrix[j][point] / augmentedMatrix[i][point];    //��ȡ��֮�䱶��
				augmentedMatrix[j][point] = 0;    //��������λ��0
				for (int k = point + 1; k < N + 1; k++) {
					augmentedMatrix[j][k] = augmentedMatrix[j][k] - augmentedMatrix[i][k] * multiple;    //�����б任
				}
			}
		}
		flag = 1;    //��ʼ��flag
		point++;    //һ��ȫ���б��ӱ仯������point����
	}

	//showMatrix(augmentedMatrix);
	printf("�������Ϊ��\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N + 1; j++)
			printf("%-7.2f", augmentedMatrix[i][j]);
		printf("\n");
	}
	printf("\n");

	int r_A = N, r_AB = N;    //Ԥ��ϵ�������Լ��������Ϊ���Ⱦ���
	for (int i = N - 1; i > 0; i--) {    //��Ϊ������Σ��������������ϱ��������������ұ���
		int j = 0;
		for (; j < N + 1; j++) {
			if (augmentedMatrix[i][j] == 0)    //��0��������������0����ѭ��
				continue;
			else {
				break;
			}
		}
		if (j == N)    //�жϷ�0����λ���������ȵĲ���
			r_A--;
		else if (j == N + 1) {
			r_A--;
			r_AB--;
		}
		else
			break;
	}


	printf("ϵ���������Ϊ%d, ����������Ϊ%d\n", r_A, r_AB);
	if (r_A == r_AB) {
		printf("�÷�������Է�������������⣬����������%d��\n", (N - r_A));
		return 0;
	}
	else if (r_A < r_AB) {
		printf("�÷�������Է������޽�\n");
		return 1;
	}
	else
	    return -1;
}

void getHomogeneousSolutions(double augmentedMatrix[N][N], int solutionsNum) {
	//���õ����Լ���ٷ�������������������������Է��������������, �ӵײ����ϵ����������⣬δ��ʵ��
	double solution[N][N];    //������,����ʹ��malloc��������һ����һЩ�Ķ�ά����
	double temp[N];    //��ʱ�洢ĳ������
	for (int i = N - solutionsNum; i < N; i++) {
		for (int j = 0; j < N; j++) {
			temp[j] = augmentedMatrix[i][j];
		}
		int xNum = 0;    //��¼δ֪������
		for (int i = 0; i < N; i++)
			if (temp[i] != 0)
				xNum++;
	}
}