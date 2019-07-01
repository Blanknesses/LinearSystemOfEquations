#include<stdio.h>
#include<math.h>
//#include<stdbool.h>    考虑VC++不包含该头文件故不使用

//n阶系数矩阵，使用全局define是为了避免malloc多次创建动态数组导致内存冲突
#define N 2 

int main() {
	double Simplest(double coefficientMatrix[N][N]);    //将矩阵化为上三角矩阵并计算行列式的值
	double getDeterminant(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1], int i);     //根据克拉默法则，对系数矩阵进行替换	                                                                                                  
	void showMatrix(double coefficientMatrix[N][N]);    //输出矩阵
	int checkType(double nonhomogeneousTerm[N][1]);    //通过系数矩阵行列式值的情况以及非齐次项进行判断齐次非齐次以及解的情况
	int solutionsInfo(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1]);    //针对非齐次线性方程组分析解的数量
	int augmentedMatrixSimplest(double augmentedMatrix[N][N + 1]);    //增广矩阵化行最简型，因为考虑到避免使用malloc导致内存冲突
	                                                                  //使用define定义数组大小，所以无法复用Simplest函数，只能进行重写 
	void getHomogeneousSolutions(double augmentedMatrix[N][N], int solutionsNum);    //计算齐次线性方程组通解，未能实现

	double coefficientMatrix[N][N];
	printf("请输入系数矩阵: ");
	for (int i = 0; i < N; i++)    //系数矩阵创建
		for (int j = 0; j < N; j++)
			scanf("%lf", &coefficientMatrix[i][j]);
	double nonhomogeneousTerm[N][1];
	printf("请输入非齐次项: ");
	for (int i = 0; i < N; i++)    //非齐次项创建
		scanf("%lf", &nonhomogeneousTerm[i][0]);

	if (checkType(nonhomogeneousTerm) == 0) {
		printf("输入的是非齐次线性方程组\n");
		double determinant[N + 1];    //克拉默法则算出的行列式值数组
		for (int i = -1; i < N; i++) {
			determinant[(i + 1)] = getDeterminant(coefficientMatrix, nonhomogeneousTerm, i);
		}

		if (determinant[0] != 0) {
			printf("该非齐次线性方程组存在唯一解\n");

			//for (int i = 0; i < N + 1; i++)
				//printf("%-7.2f", determinant[i]);

			//printf("\n");
			for (int i = 1; i < N + 1; i++) {
				printf("解%d为%.2f", i, (determinant[i] / determinant[0]));    //克拉默法则计算解
				printf("\n");
			}
		}
		else {
			int type = solutionsInfo(coefficientMatrix, nonhomogeneousTerm);
			if (type == -1) {
				printf("数据异常！");
			}

			//根据线代书本上线性方程组无穷多解的解决方法，无法有效设计出对应C语言代码实现设未知数表示各个未知量
			//之后写成向量形式得出解向量，故仅在次判断解的情况

		}
	}
	else if (checkType(nonhomogeneousTerm) == 1) {
		printf("输入的是齐次线性方程组\n");
		double determinant;    //计算系数矩阵行列式
		determinant = Simplest(coefficientMatrix);

		if (determinant != 0) {
			printf("该其次线性方程组仅有零解\n");
		}
		else {
			printf("该其次线性方程组存在非零解\n");
			int r = N;
			int j = 0;
			for (int i = N - 1; i >= 0; i--) {    //通过获取系数矩阵的秩来判断解的情况
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
				printf("非零解解向量存在%d个\n", (N - r));
			else
				printf("系数矩阵无意义！");   //系数矩阵全为0

			//for (int i = 1; i < N + 1; i++) {
				//printf("解%d为%.2f", i, (determinant[i] / determinant[0]));
				//printf("\n");
		//	}
		}
	}
	else {
		printf("数据异常！");
	}

	
	//system("pause");
	int pause;
	scanf("%d", &pause);    //替代system("pause");
	
	return 0;
}

double Simplest(double coefficientMatrix[N][N]) {

	int point = 0;    //指向第一个非零列
	double determinant = 1;    //计算行列式的值
	int flag = 1;    //判断在第一行首位为零的情况下，是否进行行变换，即是否第一列均为零
	for (int i = 0; i < N; i++) {
		if (coefficientMatrix[i][point] == 0) {
			for (int j = i + 1; j < N; j++) {    //在第一行第一个为零的情况下，去其他行找是否存在首位非零行
				if (coefficientMatrix[j][point] == 0)
					continue;
				else {    //若存在，进行行变换，将标志位置为false，进行后续操作
					for (int n = point; n < N; n++) {
						double temp = coefficientMatrix[j][n];
						coefficientMatrix[j][n] = coefficientMatrix[i][n];
						coefficientMatrix[i][n] = temp;
					}
					flag = 0;
					determinant = determinant * (-1);
				}
			}
			if (flag == 1) {    //若第一列均为0，将point后移第二列，将i--，继续循环重复上述过程
				point++;
				i--;
				continue;
			}
		}

		//showMatrix(coefficientMatrix);

		for (int j = i + 1; j < N; j++) {
			double multiple;    //第一个非零数与上一行第一个非零数的倍数
			if (coefficientMatrix[j][point] == 0) {    //与上面操作类似进行零值操作
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
				if (point >= N)    //存在point自增后导致数据越界的可能，故在此进行判断
					break;
				multiple = coefficientMatrix[j][point] / coefficientMatrix[i][point];    //获取行之间倍数
				coefficientMatrix[j][point] = 0;    //被减行首位置0
				for (int k = point + 1; k < N; k++) {
					coefficientMatrix[j][k] = coefficientMatrix[j][k] - coefficientMatrix[i][k] * multiple;    //倍减行变换
				}
			}			
		}
		//showMatrix(coefficientMatrix);
		flag = 1;    //初始化flag
		point++;    //一次全部行倍减变化结束，point后移
	}

	showMatrix(coefficientMatrix);
	//printf("\n");
	
	for (int i = 0; i < N; i++) 
		determinant = determinant * coefficientMatrix[i][i];    //上三角矩阵行列式计算
	return determinant;
}

double getDeterminant(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1], int i) {
	double Matrix[N][N];    //系数矩阵赋值，避免直接操作系数矩阵导致初始系数矩阵更改
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Matrix[i][j] = coefficientMatrix[i][j];
	double determinant = 0;
	//showMatrix(Matrix);
	if (i == -1) {    //-1表示为初始系数矩阵行列式求值
		printf("系数矩阵最简形：\n");
		determinant = Simplest(Matrix);
		//determinant = 1;
	}
	else {    //从第一列开始替换计算行列式值
		for (int j = 0; j < N; j++)
			Matrix[j][i] = nonhomogeneousTerm[j][0];
		//showMatrix(Matrix);
		//printf("\n");
		printf("第%d列替换后矩阵最简形：\n", (i + 1));
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
			type = 0;    //非齐次线性方程组
			break;
		}
		else {
			type = 1;    //齐次线性方程组
		}
	}
	return type;
}

int solutionsInfo(double coefficientMatrix[N][N], double nonhomogeneousTerm[N][1]) {
	double augmentedMatrix[N][N + 1];    //定义增广矩阵
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			augmentedMatrix[i][j] = coefficientMatrix[i][j];
	for (int i = 0; i < N; i++)
		augmentedMatrix[i][N] = nonhomogeneousTerm[i][0];	

	return augmentedMatrixSimplest(augmentedMatrix);
}

int augmentedMatrixSimplest(double augmentedMatrix[N][N + 1]) {
	int point = 0;    //指向第一个非零列

	int flag = 1;    //判断在第一行首位为零的情况下，是否进行行变换，即是否第一列均为零
	for (int i = 0; i < N; i++) {
		if (augmentedMatrix[i][point] == 0) {
			for (int j = i + 1; j < N; j++) {    //在第一行第一个为零的情况下，去其他行找是否存在首位非零行
				if (augmentedMatrix[j][point] == 0)
					continue;
				else {    //若存在，进行行变换，将标志位置位false，进行后续操作
					for (int n = point; n < N + 1; n++) {
						double temp = augmentedMatrix[j][n];
						augmentedMatrix[j][n] = augmentedMatrix[i][n];
						augmentedMatrix[i][n] = temp;
					}
					flag = 0;
				}
			}
			if (flag == 1) {    //若第一列均为0，将point后移第二列，将i--，继续循环重复上述过程
				point++;
				i--;
				continue;
			}
		}

		//showMatrix(coefficientMatrix);

		for (int j = i + 1; j < N; j++) {
			double multiple;    //第一个非零数与上一行第一个非零数的倍数
			if (augmentedMatrix[j][point] == 0) {    //与上面操作类似进行零值操作
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
				multiple = augmentedMatrix[j][point] / augmentedMatrix[i][point];    //获取行之间倍数
				augmentedMatrix[j][point] = 0;    //被减行首位置0
				for (int k = point + 1; k < N + 1; k++) {
					augmentedMatrix[j][k] = augmentedMatrix[j][k] - augmentedMatrix[i][k] * multiple;    //倍减行变换
				}
			}
		}
		flag = 1;    //初始化flag
		point++;    //一次全部行倍加变化结束，point后移
	}

	//showMatrix(augmentedMatrix);
	printf("增广矩阵为：\n");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N + 1; j++)
			printf("%-7.2f", augmentedMatrix[i][j]);
		printf("\n");
	}
	printf("\n");

	int r_A = N, r_AB = N;    //预设系数矩阵以及增广矩阵为满秩矩阵
	for (int i = N - 1; i > 0; i--) {    //因为是最简形，所以行由下至上遍历，列由左至右遍历
		int j = 0;
		for (; j < N + 1; j++) {
			if (augmentedMatrix[i][j] == 0)    //遇0继续遍历，遇非0跳出循环
				continue;
			else {
				break;
			}
		}
		if (j == N)    //判断非0出现位置来进行秩的操作
			r_A--;
		else if (j == N + 1) {
			r_A--;
			r_AB--;
		}
		else
			break;
	}


	printf("系数矩阵的秩为%d, 增广矩阵的秩为%d\n", r_A, r_AB);
	if (r_A == r_AB) {
		printf("该非齐次线性方程组存在无穷多解，解向量存在%d个\n", (N - r_A));
		return 0;
	}
	else if (r_A < r_AB) {
		printf("该非齐次线性方程组无解\n");
		return 1;
	}
	else
	    return -1;
}

void getHomogeneousSolutions(double augmentedMatrix[N][N], int solutionsNum) {
	//利用迭代以及穷举法尝试求解无穷多解情况下其次线性方程组非零解解向量, 从底层往上迭代，逐层求解，未能实现
	double solution[N][N];    //解向量,避免使用malloc，创建了一个大一些的二维数组
	double temp[N];    //暂时存储某行数据
	for (int i = N - solutionsNum; i < N; i++) {
		for (int j = 0; j < N; j++) {
			temp[j] = augmentedMatrix[i][j];
		}
		int xNum = 0;    //记录未知数个数
		for (int i = 0; i < N; i++)
			if (temp[i] != 0)
				xNum++;
	}
}