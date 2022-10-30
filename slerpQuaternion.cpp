#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2/opencv.hpp>
#include <vector>


using namespace std;
using namespace cv;
using namespace Eigen;

// 四元数球面线性插值简化方法：v'=v1*cosθ' + v⊥*sinθ'，原理见公众号推送文章
 Quaterniond slerp(double t, Quaterniond &q1, Quaterniond &q2)
 {
	 q1.normalize();
	 q2.normalize();
	 double q1_array[4]={q1.w(),q1.x(),q1.y(),q1.z()};
	 double q2_array[4]={q2.w(),q2.x(),q2.y(),q2.z()};
	 Mat q1_Mat(4,1,CV_64F, q1_array);
   	 Mat q2_Mat(4, 1, CV_64F, q2_array);

	 double q1_q2_dot=q1_Mat.dot(q2_Mat);
	 double norm2 = norm(q1_Mat, NORM_L2)*norm(q2_Mat, NORM_L2); // L2为2范数 即q1的模乘以q2的模
     double cosTheta = q1_q2_dot/norm2;    // cosθ = q1.*q2/(||q1||+||q2||)
	 Mat result;
	 Eigen::Quaterniond quaternion_result;
	 if(cosTheta>0.9995){
		result = (1.0f-t)*q1_Mat + t*q2_Mat;
	 }else{
		double theta = acosf(cosTheta);
        double thetaT = theta*t;    // t为0-1的小数 thetaT即为 θ‘
        // q1 q2都为向量，现在要求q1的垂直向量qperp
        // 把q2进行向量分解 q2=qperp*sinθ + q1*cosθ
        // 解出qperp
        Mat qperp = (q2_Mat - cosTheta*q1_Mat)/sinf(theta); // qperp即为V⊥，即q1的垂直向量
        result = q1_Mat*cosf(thetaT) + qperp*sinf(thetaT);
	 }

	 result = result / norm(result, NORM_L2);
    // Mat 转化为四元数
    quaternion_result.w() = result.at<double>(0,0);
    quaternion_result.x() = result.at<double>(1,0);
    quaternion_result.y() = result.at<double>(2,0);
    quaternion_result.z() = result.at<double>(3,0);

    return quaternion_result;

 }

int main ( int argc, char** argv )
{
	double t_img(700901880170406), t1_imu(700901879318945), t2_imu(700901884127851);
	Quaterniond q1 = Eigen::Quaterniond(0.509339,0.019188, 0.049596, 0.858921);
	Quaterniond q2 = Eigen::Quaterniond(0.509443,0.018806, 0.048944, 0.858905);
	double t = (t_img - t1_imu) / (t2_imu - t1_imu);
    Quaterniond q_slerp = slerp(t, q1, q2);
	cout<<"插值后的四元数：q_slerp =\n"<< q_slerp.coeffs() <<endl;  //coeffs的顺序是(x,y,z,w)
	
    return 0;
}
