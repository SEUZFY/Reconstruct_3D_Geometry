/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>



using namespace easy3d;



/* helpful for debugging, remove later -------------------------------------------------------------*/
namespace debugger {
    void PrintVector(const Vector& v)
    {
        for (std::size_t i = 0; i != v.size(); ++i)std::cout << v[i] << " ";
        std::cout << '\n';
    }

    void PrintPoints(const std::vector<Vector2D>& points)
    {
        int count = 0;  // display the first 10 points
        for (const auto& p : points)
        {
            std::cout << p.x() << "," << p.y() << '\n';
            ++count;
            if (count > 10)break;
        }
    }

    void PrintMatrix(const Matrix& M)
    {
        auto Nrows = M.rows();
        auto Ncols = M.cols();
        for (auto i = 0; i != Nrows; ++i)
        {
            for (auto j = 0; j != Ncols; ++j)
            {
                std::cout << M(i, j) << ",";
            }
            std::cout << '\n';
        }
    }

}



/* implementations ---------------------------------------------------------------------------------*/
namespace GEO1016_A2 {
    /*
    * if the input is valid
    * CONSTRAINTS
    * (1) correspondences >= 8
    * (2) points_0.size() == points_1.size()
    */
    bool isInputValid(const std::vector<Vector2D>& points_0, const std::vector<Vector2D>& points_1)
    {
        /* at least 8 correspondences */
        if (points_0.size() < 8 || points_1.size() < 8)
        {
            LOG(ERROR) << "insufficient correspondences\n";
            return false;
        }

        /* sizes don't match */
        if (points_0.size() != points_1.size())
        {
            LOG(ERROR) << "point size MUST match\n";
            return false;
        }

        return true;
    }



    /*
    * construct transform matrix for normalization
    * explanation:
    * (tx, ty): denote the center of the image
    * pk(k = 1,2,...N): denote the centered pixel point(NOT original pixel point)
    * dc: denote the distance of each new calculated point to the origin(0, 0)
    * where dc = sqrt(SUM pk^2)
    * dc divided by the number of tracked points(N), results in the average distance to the origin:
    * avg = dc/N
    * in order to make the average distance of a point p from the origin is equal to ¡Ì2:
    * avg * s = ¡Ì2, thus s = ¡Ì2/avg, s is the scaling factor
    *
    * Then the transform matrix for one image is:
    *     s    0   -stx
    * T = 0    s   -sty     -> T is a 3 by 3 matrix
    *     0    0     1
    *
    * @return:
    * std::pair<Matrix33, bool>
    * first element is the constructed transform matrix
    * second element is a bool variable indicating whether the matrix is susseccfully constructed
    */
    std::pair<Matrix33, bool> getNormalizeTransformMatrix(
        const std::vector<Vector2D>& points)
    {
        /* elements will be returned */
        Matrix33 T;
        bool T_valid = false;


        /* get transform matrix --------------------------------------------------------------------*/
        double sumX = 0;  /* for calculating image's center */
        double sumY = 0;
        const double N = static_cast<double>(points.size());  /* here N is guaranteed to be larger than 0(ifInputValid() gets executed first) */

        for (const auto& p : points)
        {
            sumX += p.x();
            sumY += p.y();
        }
        if (sumX < 1e-8 || sumY < 1e-8)  /* sumX or sumY is considered equal to 0 */
        {
            LOG(ERROR) << "please check the divisor for x and y coordinates\n";
            return std::make_pair(T, T_valid);  /* T_valid remains false, will not trigger further process */
        }

        /* image center : (tx, ty) */
        double tx = sumX / N;
        double ty = sumY / N;

        /* scale factor */
        double sum_dist = 0;
        for (const auto& p : points)
            sum_dist += sqrt((p.x() - tx) * (p.x() - tx) + (p.y() - ty) * (p.y() - ty));
        double avg_dc = sum_dist / N;
        if (avg_dc < 1e-8)
        {
            LOG(ERROR) << "please check the average distance to the origin\n";
            return std::make_pair(T, T_valid);  /* T_valid remains false, will not trigger further process */
        }
        double s = sqrt(2) / avg_dc;

        /* construct the transform matrix for image_0, set the T0_flag to true */
        T.set_row(0, { s, 0, -s * tx });
        T.set_row(1, { 0, s, -s * ty });
        T.set_row(2, { 0, 0, 1 });
        T_valid = true;  /* T is successfully constructed */
        /* get transform matrix --------------------------------------------------------------------*/

        return std::make_pair(T, T_valid);
    }



    /*
    * NormalizePoints
    * apply transform matrix to 2D points
    *
    * @param:
    * points need to be normalized, transform matrix M
    *
    * @return:
    * std::pair<std::vector<Vector2D>, bool>
    * first element is the normalized points set
    * second element is a bool variable indicating whether the normalized points set is valid.
    */
    std::vector<Vector2D> NormalizePoints(
        const std::vector<Vector2D>& points, const Matrix33& T)
    {
        // elements will be returned
        std::vector<Vector2D> np;
        np.reserve(points.size());

        // normalize points
        for (const auto& p : points)
        {           
            Vector3D q = T * p.homogeneous();  // normalize points: q = T*p
            np.emplace_back();
            np.back() = q.cartesian();  // add the normalized 2d coordinates to the result vector
        }

        return np;
    }



    /*
    * get initial Fundamental (matrix)
    * Wf = 0
    * construct W and solve it using SVD
    * use the last column of V to form Fundamental matrix(F needs to be further processed)
    * @return:
    * std::pair<Matrix33, bool>, first is the F matrix, second is its validation status
    */
    std::pair<Matrix33, bool> getInitialFundamental(
        const std::vector<Vector2D>& normal_points_0,
        const std::vector<Vector2D>& normal_points_1)
    {
        // elements will be returned: initial Fundamental Matrix, its status
        Matrix33 F;
        bool F_valid = false;
        // elements will be returned: initial Fundamental Matrix, its status


        // initialize W matrix: Wf = 0, W: m by n matrix
        int m = static_cast<int>(normal_points_0.size());  // number of rows
        int n = 9;  // according to notes, number of columns = 9
        Matrix W(m, n);

        // define u, v for normalpoints_0, u_, v_ for normalpoints_1
        double u{}, v{}, u_{}, v_{};
        for (int i = 0; i != normal_points_0.size(); ++i)
        {
            u  = normal_points_0[i].x(); v  = normal_points_0[i].y();
            u_ = normal_points_1[i].x(); v_ = normal_points_1[i].y();
            W.set_row(i, { u * u_, v * u_, u_, u * v_, v * v_, v_, u, v, 1 });
        }


        // Now matrix W (constructed by normalized points) is constructed ----------------------------


        // solve W using SVD
        Matrix U(m, m, 0.0);   // initialized with 0s
        Matrix S(m, n, 0.0);   // initialized with 0s
        Matrix V(n, n, 0.0);   // initialized with 0s
        svd_decompose(W, U, S, V);

        // Form F using the last column of V
        Vector vlc = V.get_column(V.cols() - 1);  // get the last column of V
        F.set_row(0, { vlc[0], vlc[1], vlc[2] });
        F.set_row(1, { vlc[3], vlc[4], vlc[5] });
        F.set_row(2, { vlc[6], vlc[7], vlc[8] });
        F_valid = true;

        return std::make_pair(F, F_valid);
    }



    /*
    * get Fundamental (Matrix)
    * process:
    * (1)decompose F using SVD
    * (2)rank-2 approximation - make S(2, 2) = 0
    * (3)denormalization - scale F such that F(2, 2) = 1.0
    * 
    * @param:
    * initial_F: initial Fundamental matrix from Wf=0, solve W using SVD
    * T:  transform matrix for points_0
    * T_: transform matrix for points_1
    * @return:
    * std::pair<Matrix33, bool>, first is the processed F matrix, second is its validation status
    */
    std::pair<Matrix33, bool> getFundamental(
        const Matrix33& initial_F, const Matrix33& T, const Matrix33& T_)
    {
        // elements will be returned ----------------------------------
        Matrix33 F;
        bool F_valid = false;
        // elements will be returned ----------------------------------

        // decompose F using SVD
        int m = initial_F.rows();
        int n = initial_F.cols();

        Matrix U(m, m, 0.0);   // initialized with 0s
        Matrix S(m, n, 0.0);   // initialized with 0s
        Matrix V(n, n, 0.0);   // initialized with 0s
        svd_decompose(initial_F, U, S, V);

        // rank-2 approximation - make S(2, 2) = 0
        S(2, 2) = 0; 
        
        // update F
        F = U * S * V.transpose();
        
        // denormalize F: F = T'_transpose() * F * T scale F such that F(2, 2) = 1.0
        F = T_.transpose() * F * T;
        
        // scale F such that F(2, 2) = 1.0 (F is up to scale)
        const double threshold = 1e-10;  // check whether F(2, 2) is 0, NB: threshold need to be small enough
        if (abs(F(2, 2) < threshold))
        {
            LOG(ERROR) << "the last element in Fundamental matrix is considered equal to 0\n"
                << "please check matrix calculation or maybe the threshold is set too large\n";
            return std::make_pair(F, F_valid);  // F_valid remains false, will not trigger further process
        }
        const double scale = 1.0 / F(2, 2);
        F = F * scale;  // scale F

        if (abs(F(2, 2) - 1.0) > threshold)
        {
            LOG(ERROR) << "after scaling, F(2, 2) is not equal to 1\n"
                << "please check the scaling process or maybe the threshold is set too small\n";
        }
        F_valid = true;  // if F is correctly constructed, set its validation status to true

        return std::make_pair(F, F_valid);
    }



    /*
    * get intrinsicK matrix
    * assume skewness is 0
    */
    Matrix33 getIntrinsicK(double fx, double fy, double cx, double cy)
    {
        return Matrix33(
            fx, 0, cx,
            0, fy, cy,
            0, 0, 1
        );
    }



    /*
    * get Essential matrix
    * E = K'.transpose() * F * K
    * in this assignment K' = K, since only one camera is used
    * and assume skewness = 0
    * therefore E = K.transpose() * F * K
    * where 
    *     fx  0  cx
    * E = 0  fy  cy
    *     0   0   1
    * 
    * @param:
    * F - Fundamental matrix
    * K - intrinsic matrix
    * 
    * @return:
    * Essential matrix
    */
    Matrix33 getEssential(const Matrix33& F, const Matrix33& K)
    {
        return K.transpose() * F * K;
    }



    /*
    * struct for storing the possible R, t combinations
    */
    struct Rt {
        std::vector<Matrix33> possibleR;  // store two possible Rotation matrices
        std::vector<Vector3D> possiblet;  // store two possible translation vectors
        Rt() { possibleR.reserve(2); possiblet.reserve(2); }      
    };  

    /*
    * getPossibleRt
    * use Essential matrix to get possible R, t combinations(4 possible combinations in total)
    * 
    * @param:
    * E - Essential matrix
    * 
    * @return:
    * Rt object - contains 4 possible R and ts.
    */
    Rt getPossibleRt(const Matrix33& E)
    {
        // result Rt object
        Rt rt;

        // decompose E using SVD
        int m = E.rows();
        int n = E.cols();
        Matrix U(m, m, 0.0);   // initialized with 0s
        Matrix S(m, n, 0.0);   // initialized with 0s
        Matrix V(n, n, 0.0);   // initialized with 0s
        svd_decompose(E, U, S, V);  // E = USV.transpose()

        // define W matrix
        Matrix33 W(
            0, -1, 0,
            1, 0, 0,
            0, 0, 1
        );

        // get possible R (1)
        rt.possibleR.emplace_back();
        rt.possibleR.back() = determinant(U * W * V.transpose()) * U * W * V.transpose();

        // get possible R (2)
        rt.possibleR.emplace_back();
        rt.possibleR.back() = determinant(U * W.transpose() * V.transpose()) * U * W.transpose() * V.transpose();

        // get possible t (1)
        rt.possiblet.emplace_back();
        rt.possiblet.back() = U.get_column(U.cols() - 1);

        // get possible t (2)
        rt.possiblet.emplace_back(); 
        rt.possiblet.back() = -U.get_column(U.cols() - 1);

        return rt;
    }



    /*
    * getProjectionMatrix
    * p = MP
    * M = K[R t] - M is 3 by 4 matrix
    * for image_0, R = I, t = [0, 0, 0]
    * for image_0, R, t - 4 different combinations
    * 
    * @param:
    * K - intrinsic matrix - 3 by 3
    * R - rotation matrix - 3 by 3
    * t - translation vector - 3d vector
    * 
    * @return:
    * projection Matrix corresponding R and t - 3 by 4 matrix
    */
    Matrix34 getProjectionMatrix(
        const Matrix33& K,
        const Matrix33& R,
        const Vector3D& t)
    {
        Matrix34 rt_matrix(
            R(0, 0), R(0, 1), R(0, 2), t[0],
            R(1, 0), R(1, 1), R(1, 2), t[1],
            R(2, 0), R(2, 1), R(2, 2), t[2]
        );
        return K * rt_matrix;
    }



    /*
    * getTriangulatedPoints3D
    * get 3D points from triangulation
    * 
    * @param:
    * M -  projection matrix for image_0 (World CRS)
    * M_ - projection matrix for image_1
    * points_0 - 2D points in image_0
    * points_1 - 2D points in image_1
    * 
    * @return:
    * std::vector<Vector3D> - containing triangulated 3D points
    */
    std::vector<Vector3D> getTriangulatedPoints3D(
        const Matrix34& M,
        const Matrix34& M_,
        const std::vector<Vector2D>& points_0,
        const std::vector<Vector2D>& points_1)
    {
        std::vector<Vector3D> points_3d;  // element will be returned
        points_3d.reserve(points_0.size());

        // define M1, M2, M3 - each row in M -----------------------------------------------       
        Vector4D M1(M(0, 0), M(0, 1), M(0, 2), M(0, 3));  // M1     
        Vector4D M2(M(1, 0), M(1, 1), M(1, 2), M(1, 3));  // M2      
        Vector4D M3(M(2, 0), M(2, 1), M(2, 2), M(2, 3));  // M3

        // define M1_, M2_, M3_ - each row in M_ ------------------------------------------- 
        Vector4D M1_(M_(0, 0), M_(0, 1), M_(0, 2), M_(0, 3));  // M1_ 
        Vector4D M2_(M_(1, 0), M_(1, 1), M_(1, 2), M_(1, 3));  // M2_
        Vector4D M3_(M_(2, 0), M_(2, 1), M_(2, 2), M_(2, 3));  // M3_

        // loop through each correspondence ------------------------------------------------
        for (std::size_t i = 0; i != points_0.size(); ++i)
        {
            // construct matrix A:
            //  xM3  - M1
            //  yM3  - M2
            // x'M3' - M1'
            // y'M3' - M2'
            double x  = points_0[i].x(), y  = points_0[i].y();
            double x_ = points_1[i].x(), y_ = points_1[i].y();

            Matrix44 A;

            A.set_row(0, x  * M3  - M1 );
            A.set_row(1, y  * M3  - M2 );
            A.set_row(2, x_ * M3_ - M1_);
            A.set_row(3, y_ * M3_ - M2_);

            // solve A using SVD
            int m = A.rows();
            int n = A.cols();
            Matrix U(m, m, 0.0);   // initialized with 0s
            Matrix S(m, n, 0.0);   // initialized with 0s
            Matrix V(n, n, 0.0);   // initialized with 0s
            svd_decompose(A, U, S, V);

            Vector4D vlc = V.get_column(V.cols() - 1);  // get the last column of V
            points_3d.emplace_back();
            points_3d.back() = vlc.cartesian();  // add the 3d point to the result vector
        }

        return points_3d;
    }
    


    /*
    * struct to store the result R, t
    * and the triangulated 3d points
    */
    struct Result {
        Matrix33 R;
        Vector3D t;
        std::vector<Vector3D> points3D;
        Result(std::size_t size = 160) { points3D.reserve(size); }
    };

    /*
    * CountPositiveZ
    * count the number of positive z values of 3d points in World/Relative CRS
    *
    * @param:
    * std::vector<Vector3D> - 3d points
    *
    * @return:
    * the number of points with positive z values
    */
    std::size_t CountPositiveZ(const std::vector<Vector3D>& pts_3d)
    {
        std::size_t count{};
        for (const auto& p : pts_3d)
            count = p.z() > 0 ? ++count : count;
        return count;
    }

    /*
    * getRelativePose - find best R ant t
    * @param:
    * E - Essential matrix
    * K - intrinsic matrix
    * points_0 - 2d points for image_0
    * points_1 - 2d points for image_1
    * 
    * @return:
    * result R and t and recovered 3d points
    */
    Result getRelativePose(
        const Matrix33& E,
        const Matrix33& K,
        const std::vector<Vector2D>& points_0,
        const std::vector<Vector2D>& points_1)
    {
        // elements will be returned
        Result res;  // if points size is not 160, use: Result res(points_size) instead;
        
        // store the positive z values for each R, t combination
        // sequence of R,t: 00 01 10 11 e.g. count_z[1] - R[0]t[1]
        // first is the number in camera 1(World_CRS), second is the number in camera 2
        std::pair<std::size_t, std::size_t> count[4]{};
        
        // find possible R, t
        auto rt = getPossibleRt(E);

        // variables help to estimate the combinations -----------------------------------------------
        Matrix33 R;
        Vector3D t;

        Matrix33 I(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        );
        Vector3D zero_t(0.0, 0.0, 0.0);

        Matrix34 M = getProjectionMatrix(K, I, zero_t);  // projection matrix for image_0
        Matrix34 M_;  // projection matrix for image_1, 4 possibilities
        // variables help to estimate the combinations -----------------------------------------------
        
        
        /* 
        * define Lambda: getCount
        * because the following code is (basically) the same for 4 different R, t combinations
        * use Lambda to avoid repetition
        * @changeable params:
        * p_r - index of rt.possibleR, p_t - index of rt.possiblet, i - index of count array -------*/      
        auto getCount = [&](
            int p_r, int p_t, int i,
            const Rt& rt, const Matrix33& K, const Matrix& M,  // read - only
            Result& res, std::pair<std::size_t, std::size_t>(&count)[4])
        {
            R = rt.possibleR[p_r]; t = rt.possiblet[p_t];
            if ((determinant(R) - 1.0) < 1e-8)  // determinant(R) = 1.0 (within a tiny threshold)
            {
                M_ = getProjectionMatrix(K, R, t);

                // (1) points in camera 1 - World CRS
                res.points3D = getTriangulatedPoints3D(M, M_, points_0, points_1);  // in World_CRS
                count[i].first = CountPositiveZ(res.points3D);

                // (2) points in camera 2 - Q = R * P + t
                for (auto& p : res.points3D)
                    p = R * p + t;
                count[i].second = CountPositiveZ(res.points3D);

                res.points3D.clear();  // for next use
            }
        };
        /* Lambda definition 1 ---------------------------------------------------------------------*/
        

        // four different R, t combinations ----------------------------------------------------------       
        getCount(0, 0, 0, rt, K, M, res, count);  // combination R:0, t:0, write result to count[0]       
        getCount(0, 1, 1, rt, K, M, res, count);  // combination R:0, t:1, write result to count[1]
        getCount(1, 0, 2, rt, K, M, res, count);  // combination R:1, t:0, write result to count[2]
        getCount(1, 1, 3, rt, K, M, res, count);  // combination R:1, t:1, write result to count[3]
        // four different R, t combinations ----------------------------------------------------------


        /* define Lambda function : verify
        * number of positive z values in two camera CRS, need to be considered equal to points size
        * @changeable param:
        * i - index of count array -------------------------------------------------------------------
        */
        auto verify = [&](
            int i, std::pair<std::size_t, std::size_t>(&count)[4]) ->bool
        {
            std::size_t n = points_0.size();
            const double eps = 1e-3;  // due to possible noisy points, eps can be set larger
            return ((count[i].first - n < eps) && (count[i].second - n < eps));
        };
        /* Lambda definition 2 ---------------------------------------------------------------------*/


        // find best R, t combination ----------------------------------------------------------------
        if (verify(0, count))  // R:0, t:0
        {
            res.R = rt.possibleR[0];
            res.t = rt.possiblet[0];
            M_ = getProjectionMatrix(K, res.R, res.t);
            res.points3D = getTriangulatedPoints3D(M, M_, points_0, points_1);
        }
        else if (verify(1, count))  // R:0, t:1
        {
            res.R = rt.possibleR[0];
            res.t = rt.possiblet[1];
            M_ = getProjectionMatrix(K, res.R, res.t);
            res.points3D = getTriangulatedPoints3D(M, M_, points_0, points_1);
        }       
        else if (verify(2, count))  // R:1, t:0
        {
            res.R = rt.possibleR[1];
            res.t = rt.possiblet[0];
            M_ = getProjectionMatrix(K, res.R, res.t);
            res.points3D = getTriangulatedPoints3D(M, M_, points_0, points_1);
        }       
        else if (verify(3, count))  // R:1, t:1
        {
            res.R = rt.possibleR[1];
            res.t = rt.possiblet[1];
            M_ = getProjectionMatrix(K, res.R, res.t);
            res.points3D = getTriangulatedPoints3D(M, M_, points_0, points_1);
        }
        // find best R, t combination ----------------------------------------------------------------


        // return results
        return res;
    }



    /*
    * Evaluate: re-project recovered 3d points back to image plane
    * calculate the variances for image_0 and image_1
    * 
    * @param:
    * K - intrinsic matrix
    * R - Rotation matrix
    * t - translation vector
    * points_3d - recovered 3d points in World CRS (camera_0 CRS)
    * points_0  - known 2d points for image_0 (camera_0: World CRS)
    * points_1  - known 2d points for image_1 (camera_1: Relative CRS)
    * from camera_0(P) to camera_1(Q):
    * Q = R * P + t
    * 
    * @return:
    * std::pair<double, double> - first is the variance of image_0, second is of image_1
    */
    std::pair<double, double> Evaluate(
        const Matrix34& M,
        const Matrix34& M_,
        const std::vector<Vector3D>& points_3d,
        const std::vector<Vector2D>& points_0,
        const std::vector<Vector2D>& points_1)
    {        
        // Lambda - calculate the difference ---------------------------------------------------------
        auto getdiff = [&](
            const Matrix& M_projection,  // projection matrix
            const std::vector<Vector3D>& points_3D,  // recovered 3D points
            const std::vector<Vector2D>& points_2D   // compare with the know 2D points
            ) -> double
        {
            double diff{};
            double N = static_cast<double>(points_3D.size());
            for (std::size_t i = 0; i != points_3D.size(); ++i)
            {
                // re-project 3D point to 2D
                Vector3D p3D = M_projection * points_3D[i].homogeneous();
                Vector2D p2D = p3D.cartesian();
                
                // compare the difference
                const Vector2D& rhs = points_2D[i];  // original points_2D
                double dx = p2D.x() - rhs.x();
                double dy = p2D.y() - rhs.y();
                diff += sqrt(dx * dx + dy * dy);
            }
            return diff / N;
        };
        // Lambda - calculate the difference ---------------------------------------------------------
        

        // get the difference for two images ---------------------------------------------------------
        double diff_0 = getdiff(M, points_3d, points_0);  // image_0
        double diff_1 = getdiff(M_,points_3d, points_1);  // image_1
        // get the difference for two images ---------------------------------------------------------


        return std::make_pair(diff_0, diff_1);
    }



    /* ----------------------------------------------------------------------------------------------*/



    /*
    * The basic implementations are above.
    * The following is for optimization using the Levenberg-Marquardt method.
    * The code is inspired by LiangLiang and Nail's tutotial -> Tutorial_NonlinearLeastSquares/main.cpp
    */



    /* ----------------------------------------------------------------------------------------------*/



    /*
    * store the re-projected 3d points' x and y coordinates 
    * in a double array - to be passed to the optimize function: lm.optimize(&obj, x);
    */
    void setOptimizeVariables(std::vector<double>& x, int size, const std::vector<Vector3D>& points_3d)
    {
        int i{}, k{};
        while (i < size && k < points_3d.size())
        {
            x[i]     = points_3d[k].x();
            x[i + 1] = points_3d[k].y();
            x[i + 2] = points_3d[k].z();
            x[i + 3] = 1.0;

            i += 4;
            k += 1;
        }
    }


    /*
    * To use the API we need to provide our own data via "data" pointer
    * define Data struct to store the data we need to provide to the API
    */
    struct MyData {
        Matrix33 M;
        std::vector<double> points_0;
        MyData(const Matrix33& m)
        {
            points_0.reserve(3);
            points_0.emplace_back(1.0);
            points_0.emplace_back(1.0);
            points_0.emplace_back(1.0);

            M = m;
        }
    };

    /*
    * To use the Levenberg-Marquardt method to solve a non-linear least squares method, 
    * we need to define our own objective function that inherits 'Objective_LM'.
    */
    class MyObjective : public Objective_LM {
    public:
        MyObjective(int num_func, int num_var, MyData* data_) : Objective_LM(num_func, num_var, data_) 
        {
            data = data_;  // user-defined data
        }

        /**
         *  Calculate the values of each function at x and return the function values as a vector in fvec.
         *  @param  x           The current values of variables.
         *  @param  fvec        Return the value vector of all the functions.
         *  @return Return a negative value to terminate.
         *
         *  NOTE:
         *  the objective is: Min SUM || MPi - pi ||^2
         *  || MPi - pi || is the length of the residual vector: MPi - pi = r
         *  then || MPi - pi || = || r || = sqrt(r.x * r.x + r.y * r.y)
         *  
         *  thus, the objective is then turned into:
         *  Min SUM || r ||^2 = Min SUM (x^2 + y^2) = Min SUM (x0^2 + x1^2)
         *  
         * because we CAN NOT change the type of x parameter
         * we store the re-projected 3d points¡¯coordinates all in x:
         * thus x should contain 2 * points.size() elements
         * store each x, y coordinate for each point.
         * 
         * e.g.:
         * x[0], x[1] - indicating the first  point's x and y coordinate
         * x[2], x[3] - indicating the second point's x and y coordinate
         * ...
         * thus we should have points.size() * 2 variables (in our case 2 * 160 = 320 variables)
         * and  we should have points.size() * 2 functions, since for each point we need two objective functions.
         * therefore the objective would be:
         * MyObjective obj(160, 320, &data);  -> in our case
         * 
         * GOOD TO KNOW
         * the optimize functions do the squared residuals internally
         * so we only need to define the residuals in the evaluation() function body like so:
         * fvec[i] = x[i] - original_x
         * SPECIAL THANKS to Nail, he explains very clearly to me
         */
        int evaluate(const double* x, double* fvec)
        {
            for (int i = 0; i < 3; ++i)
                fvec[i] = x[i] - data->points_0[i];
            return 0;
        }
    protected:
        MyData* data;
    };

}



/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
    double fx, double fy,     /// input: the focal lengths (same for both cameras)
    double cx, double cy,     /// input: the principal point (same for both cameras)
    const std::vector<Vector2D>& points_0,  /// input: 2D image points in the 1st image.
    const std::vector<Vector2D>& points_1,  /// input: 2D image points in the 2nd image.
    std::vector<Vector3D>& points_3d,       /// output: reconstructed 3D points
    Matrix33& R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
    Vector3D& t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /*
    * function implementations are all above triangulation()
    */



    /* ----------------------------------------------------------------------------------------------*/



    /* check if the input is valid ------------------------------------------------------------------*/
    auto valid = GEO1016_A2::isInputValid(points_0, points_1);
    if (!valid)return false;
    /* check if the input is valid ------------------------------------------------------------------*/



    /* get transform matrix -------------------------------------------------------------------------*/
    auto trans_0 = GEO1016_A2::getNormalizeTransformMatrix(points_0);
    auto trans_1 = GEO1016_A2::getNormalizeTransformMatrix(points_1);
    if (trans_0.second == false || trans_1.second == false)
    {
        LOG(ERROR) << "construct Normalize Transform Matrix fail, please check\n"
            << "getNormalizeTransformMatrix() function in triangulation_method.cpp\n";
        return false;
    }
    Matrix33 T = trans_0.first; Matrix33 T_ = trans_1.first;
    /* get transform matrix -------------------------------------------------------------------------*/



    /* normalize points -----------------------------------------------------------------------------*/
    auto normal_points_0 = GEO1016_A2::NormalizePoints(points_0, T);
    auto normal_points_1 = GEO1016_A2::NormalizePoints(points_1, T_);
    /* normalize points -----------------------------------------------------------------------------*/



    /* get initial Fundamental matrix from Wf = 0, use SVD to solve W -------------------------------*/
    auto Funda = GEO1016_A2::getInitialFundamental(normal_points_0, normal_points_1);
    if (Funda.second == false)
    {
        LOG(ERROR) << "construct initial Fundamental Matrix from Wf=0 (solved using SVD) fail, please check\n"
            << "getInitialFundamental() function in triangulation_method.cpp\n";
        return false;
    }
    Matrix33 initial_F = Funda.first;
    /* get initial Fundamental matrix from Wf = 0, use SVD to solve W -------------------------------*/

    

    /* get Fundamental matrix -----------------------------------------------------------------------*/
    auto FF = GEO1016_A2::getFundamental(initial_F, T, T_);
    if (FF.second == false)
    {
        LOG(ERROR) << "get Fundamental Matrix fail, please check getFundamental() function\n"
            << "in triangulation_method.cpp\n";
        return false;
    }
    Matrix33 F = FF.first;
    /* get Fundamental matrix -----------------------------------------------------------------------*/



    /* get intrinsic matrix K -----------------------------------------------------------------------*/
    Matrix33 K = GEO1016_A2::getIntrinsicK(fx, fy, cx, cy);
    /* get intrinsic matrix K -----------------------------------------------------------------------*/



    /* get Essential matrix -------------------------------------------------------------------------*/
    Matrix33 E = GEO1016_A2::getEssential(F, K);
    /* get Essential matrix -------------------------------------------------------------------------*/



    /* get relative pose ----------------------------------------------------------------------------*/
    auto res = GEO1016_A2::getRelativePose(E, K, points_0, points_1);
    R = res.R;
    t = res.t;
    points_3d = res.points3D;
    /* get relative pose ----------------------------------------------------------------------------*/
    


    /* initially verify the result: if points_3d contains the same points as points_0 or points_1 ---*/
    if (points_3d.size() != points_0.size())
    {
        LOG(ERROR) << "recover points_3d fail, please check\n"
            << "getRelativePose() function in triangulation_method.cpp\n";
        return false;
    }
    /* initially verify the result ------------------------------------------------------------------*/



    /* get projection matrices for image_0 and image_1 ----------------------------------------------*/
    Matrix33 I(  /* identity matrix */
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    );
    Vector3D zero_t(0.0, 0.0, 0.0);  /* define 0 translation vector */
    Matrix34 M  = GEO1016_A2::getProjectionMatrix(K, I, zero_t);  /* projection matrix for camera_0(image_0) */
    Matrix34 M_ = GEO1016_A2::getProjectionMatrix(K, R, t);  /* projection matrix for camera_1(image_1) */
    /* get projection matrices for image_0 and image_1 ----------------------------------------------*/



    /* evaluate -------------------------------------------------------------------------------------*/
    auto avg_diff = GEO1016_A2::Evaluate(M, M_, points_3d, points_0, points_1);
    std::cout << "average difference for image 0 is: \n"; std::cout << avg_diff.first << '\n';
    std::cout << "average difference for image 1 is: \n"; std::cout << avg_diff.second << '\n';
    /* evaluate -------------------------------------------------------------------------------------*/



    /*
    * The basic implementations are above.
    * The following is for optimization using the Levenberg-Marquardt method.
    * The code is inspired by LiangLiang and Nail's tutotial -> Tutorial_NonlinearLeastSquares/main.cpp
    */



    /* ----------------------------------------------------------------------------------------------*/



    
    /*
    * initialize the objective function
    * 1st argument is the number of functions, 2nd the number of variables
    */
    GEO1016_A2::MyData data(I);
    GEO1016_A2::MyObjective obj(3, 3, &data);

    //int num_func = 2 * static_cast<int>(points_0.size());
    //int num_var  = 4 * static_cast<int>(points_3d.size());

    /* create an instance of the Levenberg - Marquardt(LM for short) optimizer */
    Optimizer_LM lm;

    /* initialized the variables.Later x will be modified after optimization. */
    std::vector<double> x = { 1.01, 1.01, 1.01 };
    //GEO1016_A2::setOptimizeVariables(x, 4, points_3d);  // add values to x

    /* optimize(i.e., minimizing the objective function).*/
    bool status = lm.optimize(&obj, x);

    /* retrieve the result. */

    std::cout << "the solution is:     " << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
    std::cout << "the expected result: 1, 1, 1" << std::endl;

    return status;

    //delete[] x;

    /* return true if all steps are executed correctly ----------------------------------------------*/
    return true;



    /* ----------------------------------------------------------------------------------------------*/

}