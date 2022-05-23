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

#include <iostream>
#include <easy3d/optimizer/optimizer_lm.h>
#include <Triangulation/triangulation.h>
#include <Triangulation/matrix_algo.h>


using namespace easy3d;



/// To use the Levenberg-Marquardt method to solve a non-linear least squares method, we need to define our own
/// objective function that inherits 'Objective_LM'.

// user-defined data
struct Mydata {
    double s = 2.0;
    std::vector<double> base;
    Mydata()
    {
        base.reserve(3);
        base.emplace_back(1.0);
        base.emplace_back(1.0);
        base.emplace_back(1.0);
    }
};

void setOptimizeVariables(double* x, int size, const std::vector<Vector3D>& points_3d)
{
    int i{}, k{};
    while (i < size && k < points_3d.size())
    {
        x[i] = points_3d[k].x();
        x[i + 1] = points_3d[k].y();
        x[i + 2] = points_3d[k].z();

        int a = i, b = i + 1, c = i + 2;
        std::cout << "a: " << a << "b: " << b << "c: " << c << '\n';

        i += 3;
        k += 1;

        
    }
}

class MyObjective : public Objective_LM {
public:
    MyObjective(int num_func, int num_var, Mydata* data_) : Objective_LM(num_func, num_var, data_) 
    {
        data = data_;
    }

    /**
     *  Calculate the values of each function at x and return the function values as a vector in fvec.
     *  @param  x           The current values of variables.
     *  @param  fvec        Return the value vector of all the functions.
     *  @return Return a negative value to terminate.
     *
     *  NOTE: This function implements f = (x0 - 1.0)^2 + (x1 - 1.0)^2. A client problem must implement
     *      this function to evaluate the values of each function in the expression of x.
     */
    int evaluate(const double *x, double *fvec) {
        //fvec[0] = x[0] - 1.0;
        //fvec[1] = x[1] - 1.0;
        for (int i = 0; i < 3; ++i)
            fvec[i] = x[i] - data->base[i];
        return 0;
    }
protected:
    Mydata* data;
};


int main(int argc, char **argv) {
    /// initialize the objective function
    /// 1st argument is the number of functions, 2nd the number of variables
    /// the number of functions must > 1 ?
    /// 
    /// user-defined daya
    Mydata data;
    MyObjective obj(3, 3, &data);

    /// create an instance of the Levenberg-Marquardt (LM for short) optimizer
    Optimizer_LM lm;

    /// initialized the variables. Later x will be modified after optimization.
    /// initial values will affect the results
    std::vector<double> x = { 5.0, -16.0, 3.0};

    /// optimize (i.e., minimizing the objective function).
    bool status = lm.optimize(&obj, x);

    /// retrieve the result.
    //std::cout << "the solution is:     " << x[0] << ", " << x[1] << "," << x[2] << std::endl;
    //std::cout << "the expected result: 0, 0" << std::endl;

    //return status;

    /*int size = 480;
    double x[480] = { 0 };
    std::vector<Vector3D> points_3d(160, Vector3D(0, 0, 0));
    setOptimizeVariables(x, size, points_3d);*/

    return 0;
}
