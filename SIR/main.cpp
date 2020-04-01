//
//  main.cpp
//  SIR
//
//  Created by wuhan on 2020/3/29.
//  Copyright © 2020 wuhan. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

void SIR(double&, double&, double&, double&, double&, double&); // SIR模型
void createData();      //用SIR模型生成数据
void getRes();          //求解参数
int main(int argc, const char * argv[]) {
    createData();       //生成数据
    getRes();           //计算参数和R0
}

void SIR(double& I, double& S, double& N, double& Pcon, double& Pspr, double& Prec)
{
    double R = (Prec * I);//治愈人数
    I = I + (Pcon * Pspr * (I / N) * S) - R;
    S = N - I - R;
    N = N - R;
}

void createData()
{
    double I = 10.0;        //初始感染人数
    double S = 9990.0;      //初始易感染人数
    double N = I + S;       //初始总人数
    double R = 0.0;         //治愈人数
    double Pcon = 1.0;      //接触率
    double Pspr = 0.5;      //传染率
    double Prec = 0.1;      //治愈率
    int iter = 100;         //迭代次数
    
    //创建文件
    std::string I_file_name = "./data/I.txt";
    std::string S_file_name = "./data/S.txt";
    std::string R_file_name = "./data/R.txt";
    
    std::ofstream fout_I(I_file_name.c_str());
    std::ofstream fout_S(S_file_name.c_str());
    std::ofstream fout_R(R_file_name.c_str());
    
    if(!fout_I || !fout_S || !fout_R)
    {
        std::cout << "文件创建失败！"<< std::endl;
    }
    else
    {
        //写入数据
        for(int i = 0;i < iter; i++)
        {
            SIR(I, S, N, Pcon, Pspr, Prec);
            fout_I << I <<std::endl;
            fout_S << S <<std::endl;
            fout_R << (R = (Prec * I)) <<std::endl;
        }
        //关闭文件
        fout_I.close();
        fout_S.close();
        fout_R.close();
    }
}

void getRes()
{
    //打开文件
    std::string I_file_name = "./data/I.txt";
    std::string S_file_name = "./data/S.txt";
    std::string R_file_name = "./data/R.txt";
    
    std::ifstream fin_I(I_file_name.c_str());
    std::ifstream fin_S(S_file_name.c_str());
    std::ifstream fin_R(R_file_name.c_str());
    
    if(!fin_I || !fin_S || !fin_R)
    {
        std::cout<<"打开文件失败！"<<std::endl;
    }
    else
    {
        std::string line;       //读取的字符串
        
        int data_amount = 0;    //数据行数
        while(getline(fin_I,line))  //计算数据行数
        {
            data_amount++;
        }
        fin_I.clear();      //因为读到了文件尾部clear()一下
        fin_I.seekg(0);     //回到文件头
        
        //建立感染人数矩阵
        Eigen::MatrixXd I_mat(data_amount, 1);
        
        //建立易感染人数矩阵
        Eigen::MatrixXd S_mat(data_amount, 1);
        
        //建立治愈人数矩阵
        Eigen::MatrixXd R_mat(data_amount, 1);
        
        //读入感染人数数据
        for(int i = 0;i < data_amount; i++)
        {
            if(getline(fin_I,line))
            {
                I_mat(i,0) = atof(const_cast<const char *>(line.c_str()));;
            }
        }
        
        //读入易感染人数数据
        for(int i = 0;i < data_amount; i++)
        {
            if(getline(fin_S,line))
            {
                S_mat(i,0) = atof(const_cast<const char *>(line.c_str()));
            }
        }
        
        //读入治愈数据
        for(int i = 0;i < data_amount; i++)
        {
            if(getline(fin_R,line))
            {
                R_mat(i,0) = atof(const_cast<const char *>(line.c_str()));
            }
        }
        
        //关闭文件
        fin_I.close();
        fin_S.close();
        fin_R.close();

        //建立(I / (I+S) * S)矩阵
        Eigen::MatrixXd I_I_S_S(data_amount,1);
        I_I_S_S = (1.0 / (I_mat.array() + S_mat.array())) * I_mat.array() * S_mat.array();
        
        //建立I(n+1)矩阵
        Eigen::MatrixXd In_1(data_amount, 1);
        In_1=I_mat;
        for(int i = 0;i < data_amount - 1; i++)
        {
            In_1(i,0) = In_1(i + 1,0);
        }
        
        //建立I(n+1)-I(n)矩阵
        Eigen::MatrixXd In_1_In(data_amount, 1);
        In_1_In =In_1 - I_mat;

        //合并矩阵 I/(I+S)*S 和 I_mat
        Eigen::MatrixXd Combin(data_amount, 2);
        Combin<<I_I_S_S,I_mat;
        
        //建立结果矩阵
        Eigen::MatrixXd res(2,1);
        
        //解超定方程
        res = Combin.fullPivLu().solve(In_1_In);
        
        //显示结果
        std::cout<< "P_contact * P_spread: " << res(0,0) <<std::endl;
        std::cout<< "P_recover: " << -res(1,0) <<std::endl;
        std::cout<< "R0: " << -res(0,0) / res(1,0) <<std::endl;
    }
}
