# Damped harmonical oscillator project

This github repository contains all the material that has been used for the final NSAE project.

The file “RK4.c” cointains the code for the runge kutta, and you can compile it typing:

gcc -Wall -o RK4 RK4.c -lm

./RK4

The code asks to type in the value of h, and returns 5 txt files.
1. The file “RK4-points_x.txt” contains the points (t, x(t)) of RK4.
2. The file “RK4-points_v.txt” contains the points (t, v(t)) of RK4.
3. The file “RK4-points_ES” contains the points (t, ˆx(t)) where now ˆx(t) denotes the analytical
solution at time t.
4. The file “RK4Errors.txt” contains the points (t, ε(t)) where ε(t) = |ˆx(t) − x(t)|.
5. The file “RK4-rel-errors.txt” contains the points (t, ε(t)) where ε(t) is the relative error.

The first lines of the code define the values of the parameters. Users can change them manually for
studying different cases of the problem.

Inside the folder txt_files there are all the points necessary for plotting the results.
The pdf file cointains the report of the project.
The file "plots.ipynb" is a jupyter notebook that has been used for plotting the results.
