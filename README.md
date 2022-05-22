# Reconstruct_3D_Geometry

Project for 3D Geometry Reconstruction based on the open-source project - [Easy3D](
https://github.com/LiangliangNan/Easy3D).

![linear_triangulation](https://user-images.githubusercontent.com/72781910/169713033-495cbd47-6f78-48f0-a8e2-154f89df5432.gif)

**Easy3D** is a lightweight, easy-to-use, and efficient open-source C++ library for processing and rendering 3D data. **IF YOU WANT TO USE IT** please be aware of the following information:

**COPYRIGHT**:

Copyright (C) 2015 Liangliang Nan <liangliang.nan@gmail.com>

https://3d.bk.tudelft.nl/liangliang/


**CITATION INFO**:
```
@article{easy3d2021,
  title = {Easy3{D}: a lightweight, easy-to-use, and efficient {C}++ library for processing and rendering 3{D} data},
  author = {Liangliang Nan},
  journal = {Journal of Open Source Software}，
  year = {2021},
  volume = {6},
  number = {64},
  pages = {3255},
  doi = {10.21105/joss.03255},
  url = {https://doi.org/10.21105/joss.03255}
}
```

More info about **Easy3D** please refer to: https://github.com/LiangliangNan/Easy3D.

# HOW TO USE

* Clone this project at: https://github.com/SEUZFY/Reconstruct_3D_Geometry.git

* Or download the code and open the [CMakeLists.txt](https://github.com/SEUZFY/Reconstruct_3D_Geometry/blob/master/CMakeLists.txt) file in an IDE.

Build and run this project, a viewer should pop up automatically, press `space` and the reconstructed model is shown.

**Note**: 

After 3d reconstruction, move the mouse and the model can be seen from different positions using zoom in/out

<p float="left">
  <img width="320" alt="triangulation_withwindow" src="https://user-images.githubusercontent.com/72781910/169711724-cd59f57a-bbb5-473d-a6f1-98b47df225ba.PNG">
  <img width="320" alt="triangulation_zoomout_withwindow" src="https://user-images.githubusercontent.com/72781910/169691743-fa836ec5-dee1-4947-9707-d12a012e9f61.png">
</p>

Meanwhile some helpful information should be printed to the console.

# GOOD TO KNOW

* The triangulation method is described in **details** [here](https://3d.bk.tudelft.nl/courses/geo1016/handouts/04-reconstruct_3D_geometry.pdf) and [here](https://3d.bk.tudelft.nl/courses/geo1016/handouts/03-epipolar_geometry.pdf). **IT SHOULD BE NOTED** that this explanation comes from the course notes, if you want to use it in a scientific work, you are kindly asked to mention the **ORIGINAL** author: 
  
  < Liangliang Nan <liangliang.nan@gmail.com> >
  
  < https://3d.bk.tudelft.nl/liangliang/ > 

* The `triangulation implementation` is in [Triangulation_method.cpp](https://github.com/SEUZFY/Reconstruct_3D_Geometry/blob/master/Triangulation/triangulation_method.cpp), all the other files are kindly given by [Liang Liang](https://3d.bk.tudelft.nl/liangliang/).

* The `workflow` is illustrated with some code examples [here](https://github.com/SEUZFY/Reconstruct_3D_Geometry/blob/master/Triangulation_materials/3DGeometry_from_2images_workflow.pdf).

# Evaluation

**Method**:

$\frac{\sum \sqrt{\left ( \Delta x \right )^2 + \left ( \Delta y \right )^2}}{N(points)}$

$\Delta x$: obtained `x` coordinate $-$ original `x` coordinate

$\Delta y$: obtained `y` coordinate $-$ original `y` coordinate

average difference for image 0 `(points_0)` is:
`1.04006`

average difference for image 1 `(points_1)` is:
`0.912123`

