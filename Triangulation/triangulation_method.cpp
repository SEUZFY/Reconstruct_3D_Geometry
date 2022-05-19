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

namespace GEO1016_debugger {
    void PrintMatrix33(const Matrix33& M)
    {
        std::cout << M(0, 0) << "," << M(0, 1) << "," << M(0, 2) << '\n';
        std::cout << M(1, 0) << "," << M(1, 1) << "," << M(1, 2) << '\n';
        std::cout << M(2, 0) << "," << M(2, 1) << "," << M(2, 2) << '\n';
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
}

namespace GEO1016_A2 {
    /*
    * if the input is valid
    * CONSTRAINTS
    * (1) correspondences >= 8
    * (2) points_0.size() == points_1.size()
    */
    bool isInputValid(const std::vector<Vector2D>& points_0, const std::vector<Vector2D>& points_1)
    {
        if (points_0.size() < 8 || points_1.size() < 8)  // at least 8 correspondences
        {
            LOG(ERROR) << "insufficient correspondences\n";
            return false;
        }

        if (points_0.size() != points_1.size())  // size doesn't match
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
        // elements will be returned
        Matrix33 T;
        bool T_valid = false;
        
        // transform matrix --------------------------------------------------------------
        double sumX = 0;  // for calculating image's center
        double sumY = 0;
        const double N = static_cast<double>(points.size());  // here N is guaranteed to be larger than 0(ifInputValid() gets executed first)

        for (const auto& p : points)
        {
            sumX += p.x();
            sumY += p.y();
        }
        if (sumX < 1e-8 || sumY < 1e-8)  // sumX or sumY is considered equal to 0
        {
            LOG(ERROR) << "please check the divisor for x and y coordinates\n";
            return std::make_pair(T, T_valid);  // T_valid remains false, will not trigger further process
        }
        
        // image center: (tx, ty)
        double tx = sumX / N;
        double ty = sumY / N; 
        
        // scale factor
        double dc_squared = 0;
        for (const auto& p : points)
            dc_squared += (p.x() - tx) * (p.x() - tx) + (p.y() - ty) * (p.y() - ty);
        double dc = sqrt(dc_squared);  // dist from the origin(here origin should be the image center)
        double avg_dc = dc / N;
        if (avg_dc < 1e-8)
        {
            LOG(ERROR) << "please check the average distance to the origin\n";
            return std::make_pair(T, T_valid);  // T_valid remains false, will not trigger further process
        }
        double s = sqrt(2) / avg_dc;  // scale factor

        // construct the transform matrix for image_0, set the T0_flag to true
        T.set_row(0, { s, 0, -s * tx });
        T.set_row(1, { 0, s, -s * ty });
        T.set_row(2, { 0, 0, 1 });
        T_valid = true;  // T is successfully constructed
        // transform matrix --------------------------------------------------------------

        return std::make_pair(T, T_valid);
    }


    /*
    * NormalizePoints
    * apply transform matrix to 2D points
    * 
    * @return:
    * std::pair<std::vector<Vector2D>, bool>
    * first element is the normalized points set
    * second element is a bool variable indicating whether the normalized points set is valid.
    */
    std::pair<std::vector<Vector2D>, bool> NormalizePoints(const std::vector<Vector2D>& points)
    {
        // elements will be returned
        std::vector<Vector2D> np;
        np.reserve(points.size());
        bool np_valid = false;

        // get the transform matrix and its status
        const auto& trans = getNormalizeTransformMatrix(points);
        const auto& T = trans.first;
        const auto& T_valid = trans.second;

        if (T_valid == false)
        {
            LOG(ERROR) << "Normalize - transform matrix construction failed, please check\n"
                << "getNormalizeTransformMatrix() function in triangulation_method.cpp\n";
            return std::make_pair(np, np_valid);  // np_valid remains false and further process will not be triggered
        }

        // T_valid is true - transform matrix successfully constructed
        
        // normalize points
        for (const auto& p : points)
        {
            Vector3D q = T * p.homogeneous();  // normalize points: q = T*p
            np.push_back(q.cartesian());  // add the normalized 2d coordinates to the result vector
        }
        np_valid = true;  // result vector np is successfully constructed

        return std::make_pair(np, np_valid);
    }


    /*
    * Wf = 0
    * construct W and solve it using SVD
    * rank-2 approximation
    * denormalization
    * Fundamental matrix
    */
    Matrix33 getFundamentalMatrix(
        const std::vector<Vector2D>& points_0,
        const std::vector<Vector2D>& points_1)
    {
        // Fundamental Matrix (will be returned)
        Matrix33 F;
        
        // initialize W matrix: Wf = 0
        int Nrows = static_cast<int>(points_0.size());  // points_0.size() = points_1.size()
        int Ncols = 9;  // according to notes
        Matrix W(Nrows, Ncols);

        // Form W matrix
        for (std::size_t i = 0; i != points_0.size(); ++i)
        {

        }
        return F;
    }
}

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any provided data structures and functions. For your convenience, the\n"
                 "\tfollowing three files implement basic linear algebra data structures and operations:\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
    Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
    Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
    Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    Matrix33 A;

    /// define and initialize a 3 by 3 matrix
    Matrix33 T(1.1, 2.2, 3.3,
               0, 2.2, 3.3,
               0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    Matrix34 M(1.1, 2.2, 3.3, 0,
               0, 2.2, 3.3, 1,
               0, 0, 1, 1);

    /// set first row by a vector
    M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
    M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
    W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    int num_rows = W.rows();

    /// get the number of columns.
    int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    double value = W(1, 2);

    /// get the last column of a matrix
    Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
    Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    auto valid = GEO1016_A2::isInputValid(points_0, points_1);
    if (!valid)return false;

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    auto trans = GEO1016_A2::getNormalizeTransformMatrix(points_0);
    std::cout << "normalize transform matrix for points_0: \n";
    GEO1016_debugger::PrintMatrix33(trans.first);

    trans = GEO1016_A2::getNormalizeTransformMatrix(points_1);
    std::cout << "normalize transform matrix for points_1: \n";
    GEO1016_debugger::PrintMatrix33(trans.first);

    GEO1016_debugger::PrintPoints(points_1);  // before normalize
    auto np = GEO1016_A2::NormalizePoints(points_1);
    std::cout << "after normalize: \n";
    GEO1016_debugger::PrintPoints(np.first);  // after normalize
    
    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}